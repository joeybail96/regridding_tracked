import os
import re
import shutil
import rioxarray
import xesmf as xe
import numpy as np
from osgeo import gdal
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
import rasterio
import pandas as pd

class Organized_regriding:
    
    def __init__(self, folder_to_save_results:str, template_file_path:str, overwrite:bool = False, descriptive_name:str = '', re_load = False):
        """_summary_

        Args:
            folder_to_save_results (str): should be a path to a folder (last char should be a /)
            overwrite (bool, optional): _description_. Defaults to False.
            descriptive_name (str, optional): _description_. Defaults to ''.
        """
        
        #set input variables
        self._folder_to_save_results = folder_to_save_results
        self._overwrite = overwrite
        self._descriptive_name = descriptive_name
            # TODO: add a method to build this template file
        self._template_path:str = template_file_path
        
        # define file paths that can be used by other methods
        #       Note: all file paths were functions are applied will be defined within functions
        self._fulrez_4326_path:str = None
        self._lowrez_4326_path:str = None
        self._output_nc_path:str = None
        self._raw_tif_5070_path:str = None
        
        # build base folder, intermediate folder, and output folder
        
        if re_load:
            self._result_folder_full_path = self._folder_to_save_results + 'results_' + self._descriptive_name + '/'
            self._intermediate_folder = self._result_folder_full_path + 'intermediate_data' + self._descriptive_name + '/'; os.makedirs(self._intermediate_folder) 
            self._output_folder = self._result_folder_full_path + 'output_data' + self._descriptive_name + '/'; os.makedirs(self._output_folder) 
            self._reload_regrid_paths()
        else:
            self._result_folder_full_path = self._create_base_folder()
            self._intermediate_folder = self._result_folder_full_path + 'intermediate_data_' + self._descriptive_name + '/'; os.makedirs(self._intermediate_folder) 
            self._output_folder = self._result_folder_full_path + 'output_data_' + self._descriptive_name + '/'; os.makedirs(self._output_folder) 
        
        # define file paths that can be used by other methods
        #       Note: all file paths were functions are applied will be defined within functions
        self._fulrez_4326_path:str = None
        self._lowrez_4326_path:str = None
        self._output_nc_path:str = None
        self._raw_tif_5070_path:str = None
        
        # define paths for tif and netcdf mutations
        self.processed_tifs = {}
        self.processed_netcdfs = {}
        self.slice_indices = {}
        self.slice_paths = {}
        self.number_of_slices = {}
        self.y_slice_regriders = {}
        self.most_recent_y_slice_index = -1
        self.complete_slice_keys = []
        
        # build save file_name
        self._save_location = self._result_folder_full_path + self._descriptive_name +  'save.pkl' 
    
    def _create_base_folder(self) -> str:
    
        results_folder_full_path = self._folder_to_save_results + 'results_' + self._descriptive_name
        
        if ((not re.search(r'__v\d+',results_folder_full_path) and re.search(r'__v',results_folder_full_path))):
            raise ValueError("Folder name inappropriately uses '__v' in a way that does not represent the version number. '__v' should always be followed by nothing but numbers")
        
        # if the folder already exists than this while statement figures out the version number
        # to append to the folder or deletes the existing folder if overwrite = True
        while os.path.exists(results_folder_full_path):
        
            if self._overwrite:
                shutil.rmtree(results_folder_full_path)
                
            else:
                # extract version number from file path
                version = re.findall(r'__v([0-9]+)', results_folder_full_path)
                
                if len(version) == 0:
                    version = 1
                else:
                    version = int(version[-1]) + 1 # the [-1] is there to ensure that the last '__v' in the string is used
                    
                results_folder_full_path = results_folder_full_path.split('__v') [0] + '__v' + str(version) # builds folder string
                self._descriptive_name = self._descriptive_name.split('__v') [0] + '__v' + str(version) #updates name to include updated version number
                
        # builds the directory   
        os.makedirs(results_folder_full_path) 
        return results_folder_full_path + '/'
    
    def _rescale_resolution(self, input_file:str, output_file:str, scaling_factor:float):
        """_summary_

        Args:
            input_file (str): name of the file you want to rescale
            output_file (str): name of file and location you want to rescale the file too
            scaling_factor (float): a parameter that determines how much to scale down the input file by
                                    example :   0.5 will make tiffs length and width half as man pixels rounded down
                                                1 will result in the output tiff being the exact same as the input
        """
        
        #open input file
        ds = gdal.Open(input_file)

        # Calculate the new dimensions
        new_width = int(ds.RasterXSize * scaling_factor)
        new_height = int(ds.RasterYSize * scaling_factor)

        # set options
        
        options = gdal.WarpOptions(
            format="GTiff",
            outputType = gdal.GDT_Float32, #outputType= gdal.GDT_Int64, this needs to be changed based on what is going to be done with the data for salt I am using a float because the data may have been fractionalized when being warped.
            width=new_width,
            height=new_height,
            resampleAlg = gdal.GRA_Average,
        )
        
        #transform,
        output_ds = gdal.Warp(output_file, ds, options=options)
        
        output_ds = None #Close the output file? I am not sure If this is necessary, gdal.Warp might do this

        ds = None  # Close the input file
            
    def _transform_5070_to_4326(self, input_tif_5070):
        """
        takes tiff input in 5070 and transforms it to 432

        Args:
            tif_path_5070 (str): path to input 5070 tif file.
            tif_path_4326 (str): path to output 4326 tif file
        """

        input_ds = gdal.Open(input_tif_5070)
        
        # set options
        options = gdal.WarpOptions(
            format="GTiff",
            outputType = gdal.GDT_Float32,
            srcSRS = "EPSG:5070", 
            dstSRS ="EPSG:4326",
            resampleAlg = gdal.GRA_Average
        )
        
        #use gdal to regrid
        gdal.Warp(self._fulrez_4326_path, input_ds, options = options) # goes from srcSRS to dstSRS
        input_ds = None #close the data set

    def _convert_to_geochem_nc(self, tif_file_path:str, nc_output_path:str, template_ds:xr.DataArray, resampling_method):
        """creates a netcdf file from a tiff input file that is gridded in a lat/lon format

        Args:
            tif_file_path (str):path to the tif file that will be regridded and saved as a netcdf file
            nc_output_path (str): the location that the netcdf file should be saved
            template_ds (xr.array): an xarry object that represents a geo-chem compatible netcdf file
            sum_or_average (str, optional): determines wether the resulting netcdf file should be regridded based on a sum or average of the input data if input is True output will regrid based on sum if it is False the regridder will regrid based on average resampling w. Defaults to True.
        """
        if (resampling_method != "sum" and resampling_method != "average"):
            raise Exception("resampling_method must be a string representing either 'sum' or 'average' and it is ", resampling_method)
        
        # I am using this with statement to try to ensure that rioxarray files are closed
        with rioxarray.open_rasterio(tif_file_path) as input_ds:
            
            # Ensure there are no nan-values and replace all possible null values
            try:
                if '_FillValue' in input_ds.attrs:
                    input_ds = input_ds.where(input_ds != input_ds._FillValue, 0) #remove fill values
            except:
                pass
                
            input_ds = input_ds.where(input_ds.notnull(),0) #remove Nan, None, or NaT
            
            #convert to data set
            input_ds = input_ds.to_dataset(name = "output")
        
            #rename axes so data will work with regridder
            input_ds = input_ds.rename({'x':'lon','y':'lat'})
        
            #create regrid object
            regridder = xe.Regridder(input_ds, template_ds, "conservative")
            # suggested method to fix Latitude outside of [-90,90] can be found @ https://stackoverflow.com/questions/66177931/how-do-i-resample-a-high-resolution-grib-grid-to-a-coarser-resolution-using-xesm
                
            #regrid xarray object
            regrided_ds = regridder(input_ds) 
            
            if (resampling_method == 'sum'):
                # Converting to sum like this only works because the input data is in kg on a 1 km^2 grid.
                # Because of this, when the average is taken it effectively gives the average units of kg/km^2.
                # All that needs to be done to convert from kg/km^2 to kg is to multiply each cell by its area.
                regrided_ds = regrided_ds * template_ds["AREA"] / 1e6 # the 1e6 is to convert from m^2 to km^2
            
            regrided_ds.to_netcdf(nc_output_path)
    
    def _reload_regrid_paths(self):
        # build file paths
        self._fulrez_4326_path:str = self._intermediate_folder + r'fulrez_4326.tif'
        self._lowrez_4326_path:str = self._intermediate_folder + r'lowrez_4326.tif'
        self._output_nc_path:str = self._output_folder + r'output.nc'
        
        paths_to_test = [
        self._fulrez_4326_path,
        self._lowrez_4326_path,
        self._output_nc_path,
        ]
        
        for path in paths_to_test:
            if (not os.path.isfile(path)):
                raise Exception(f'You cannot reload regrid data because\n\t{path}\n doesn\'t exist.')
            
    def apply_array_function_to_tif(self, input_file:str, descriptive_name_for_function:str, array_function:callable):
        """ Applies a function that takes arrays as arguments to an input geotiff file and exports a file
            that function applied. Note: this function also replaces all none values with 0 and removes _FillValue from the tiffs attributes

        Args:
            input_file (str): file path to input geotiff 
            output_file (str): file path to output geotiff 
            array_function (function): a function that takes a np.array as an argument
        """
        output_file:str = self._intermediate_folder + descriptive_name_for_function + ".tif"
        self.processed_tifs[descriptive_name_for_function] = output_file
        
        with rioxarray.open_rasterio(input_file) as input_ds:
            # Ensure there are no nan-values and replace all possible null values
            try:
                if '_FillValue' in input_ds.attrs:
                    input_ds = input_ds.where(input_ds != input_ds._FillValue, 0) #remove fill values
            except:
                pass
            input_ds = input_ds.where(input_ds.notnull(),0) #remove Nan, None, or NaT
            output_ds = array_function(input_ds) #might need to cast to array using .values
            
            # Write the processed data to the output GeoTIFF
            output_ds.rio.to_raster(output_file)
            
    def apply_array_function_to_netcdf(self, input_file:str, descriptive_name_for_function:str, array_function:callable):
        
        output_file:str = self._output_folder + descriptive_name_for_function + ".nc"
        self.processed_netcdfs[descriptive_name_for_function] = output_file
        
        with xr.open_dataset(input_file) as input_ds:
            # Make a copy of the input dataset
            # output_ds = input_ds.copy(deep=True)
            
            # Extract the 'output' layer as a NumPy array
            output_np = input_ds['output'].values[0]
            
            # Replace NaN, None, or NaT with 0
            output_np = np.nan_to_num(output_np, nan=0)
            
            # Apply array function
            modified_output_np = array_function(output_np)
            
            # Update the 'output' variable with the modified array
            input_ds['output'][:] = modified_output_np
            
            # Write the processed data to the output NetCDF
            input_ds.to_netcdf(output_file)
            # input_ds = input_ds.where(input_ds.notnull(),0) #remove Nan, None, or NaT
            # output_ds = array_function(input_ds) 
            
            # # #convert to data set
            # # output_ds = output_ds.to_dataset(name = "output")
        
            # output_ds.attrs = input_ds.attrs
            
            # # Write the processed data to the output GeoTIFF
            # output_ds.to_netcdf(output_file)
            # output_ds = None
                
    def tif_5070_to_gc_netcdf(self, input_tif_5070 ,resampling_method = "average", scaling_factor:float = 0.1):
        """
        Converts a TIFF file to a NetCDF file using a given template. Note: this will not run unless you have a lot of ram.

        Args:
            fulres_4326_path (str): The path to where the intermediate full-resolution TIFF will be stored in EPSG:4326 projection.
            lowres_4326_path (str): The path to where the intermediate low-resolution TIFF will be stored file in EPSG:4326 projection.
            output_nc_path (str): The path to save the output NetCDF file.
            tif_5070_path (str): The path to the input TIFF file in EPSG:5070 projection is stored.
            template_path (str): The path to the template NetCDF file.
            scaling_factor (float, optional): The scaling factor to apply during the resolution rescaling. Defaults to 0.01

        Returns:
            None

        Examples:
            tif_to_nc('fulres.tif', 'lowres.tif', 'result.nc', 'tif_5070.tif', 'template.nc', scaling_factor=0.1)
        """
        # check to ensure that this hasn't been ran yet
        if(self._fulrez_4326_path):
            raise Exception("tif_5070_to_gc_netcdf can only be ran once per object. Try defining new regriding object.")

        # build file paths
        self._fulrez_4326_path:str = self._intermediate_folder + r'fulrez_4326.tif'
        self._lowrez_4326_path:str = self._intermediate_folder + r'lowrez_4326.tif'
        self._output_nc_path:str = self._output_folder + r'output.nc'
        self._tif_5070_to_regrid_path:str = input_tif_5070

        #open template
        template_ds=xr.open_dataset(self._template_path)
        
        #run methods
        self._transform_5070_to_4326(self._tif_5070_to_regrid_path)
        self._rescale_resolution(self._fulrez_4326_path, self._lowrez_4326_path, scaling_factor)
        self._convert_to_geochem_nc(self._lowrez_4326_path, self._output_nc_path, template_ds, resampling_method = resampling_method)
        
        # store scaling factor
        self._scaling_factor_used = scaling_factor

    def create_mask_values_above(self, input_file:str, output_file:str, threshold:float):
        """
        creates mask with ones representing values over a threshold and zeros for values under and equal to the threshold
        Parameters

        Args:
            input_file (str): location of file you want to make mask of
            location that you want to save the mask to
            minimum value to be set at one
        """
        
        ds = gdal.Open(input_file)

        # convert to numpy array to minipulate
        data = ds.ReadAsArray()

        # Create a binary mask based on the threshold
        mask = (data > threshold).astype(np.uint8)

        # package array back into a dataset
        mask_driver = gdal.GetDriverByName('GTiff')
        mask_ds = mask_driver.Create(output_file, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte)
        mask_ds.GetRasterBand(1).WriteArray(mask)

        # Set the geo referencing information for the mask
        mask_ds.SetProjection(ds.GetProjection())
        mask_ds.SetGeoTransform(ds.GetGeoTransform())

        ds = None  # Close the input file
        mask_ds = None  # Close the mask file

        #################################################################################################################################
    
    #################################################################################################################################
    # Conservation Testing
    #################################################################################################################################
    # create slice files
    def _y_slice_raster (self, number_of_steps = "n_a", plot = False):
        
        # define input and output files paths
        
        input_file_path = self._tif_5070_to_regrid_path
        self.y_slice_conservation_test_results = self._intermediate_folder + r'y_slice_conservation_test_results' + f'_{number_of_steps}_steps/'
        self.y_sliced_tiffs = self.y_slice_conservation_test_results + r"y_sliced_tiffs/"
        
        # Check if the folder already exists
        if os.path.exists(self.y_slice_conservation_test_results):
            raise Exception(f"y slice conservation test has already been ran at {number_of_steps} steps")

        # Create the folders if they don't exist
        os.makedirs(self.y_slice_conservation_test_results)
        os.makedirs(self.y_sliced_tiffs)
        
        # create a copy of self.y_sliced_tiffs
        output_folder = self.y_sliced_tiffs
        
        #open array snd convert none values to zero
        xr_array = rioxarray.open_rasterio(input_file_path)
                    # Ensure there are no nan-values and replace all possible null values
        # try:
        #     if '_FillValue' in input_ds.attrs:
        #         input_ds = input_ds.where(input_ds != input_ds._FillValue, 0) #remove fill values
        # except:
        #     pass
        
        # extract/calculate useful quantities
        y_min = xr_array.y.min()
        y_max = xr_array.y.max()
        y_range = y_max - y_min
    
        # define the step size
        if(number_of_steps == "n_a"):
            step_size = 500000 
        else:
            step_size = int(y_range/number_of_steps)
        
        # define slice y values
        slices = np.arange(y_min, y_max+1, step_size)
        
        slice_paths = []
        for index in range(len(slices)-1):
            
            slice_output_path = f"{output_folder}raster_fulrez_5070_{slices[index]:.0e}_to_{slices[index+1]:.0e}.tif"
            slice_paths.append(slice_output_path)
            
            xr_slice = xr_array.where(xr_array.y>slices[index]).where(xr_array.y<slices[index+1])
            xr_slice = xr_slice.where(xr_array.notnull(),0)
            xr_slice.rio.to_raster(slice_output_path)
            
            print(f"finished_{slices[index]:.0e}_to_{slices[index+1]:.0e}")
            if plot:
                plt.figure()
                xr_slice.plot()
                plt.title(f"{slices[0]:.0e}_to_{slices[index+1]:.0e}")
        
        self.slice_indices[f'{number_of_steps}_steps'] = slices
        self.number_of_slices[f'{number_of_steps}_steps'] = number_of_steps
        self.slice_paths[f'{number_of_steps}_steps'] = slice_paths
        self.complete_slice_keys.append(f'{number_of_steps}_steps')
        self.most_recent_y_slice_index += 1
        
    # convert array of slices
    def _convert_y_slices (self):
        
        # define slice path array, template_path and slice_indices
        slice_path_array = self.slice_paths[self.complete_slice_keys[self.most_recent_y_slice_index]]
        slice_indices = self.slice_indices[self.complete_slice_keys[self.most_recent_y_slice_index]]
        number_of_steps = self.number_of_slices[self.complete_slice_keys[self.most_recent_y_slice_index]]
        template_path = self._template_path
        
        # define input and output files paths
        # self.y_slice_conservation_test_results = self._intermediate_folder + r'y_slice_conservation_test_results' + f'_{number_of_steps}_steps/'
        self.y_sliced_ncs = self.y_slice_conservation_test_results + r"y_sliced_ncs/"

        # Create the folders if they don't exist
        os.makedirs(self.y_sliced_ncs)
        
        regriders = []
        
        for index, slice_path in enumerate(slice_path_array):
            
            destination_folder_path = self.y_sliced_ncs + f"{slice_indices[index]:.3e}_to_{slice_indices[index+1]:.3e}" + "/"
            
            regrider = Organized_regriding(destination_folder_path, template_path, descriptive_name=f"{slice_indices[index]:.3e}_to_{slice_indices[index+1]:.3e}")

            regrider.tif_5070_to_gc_netcdf(input_tif_5070= slice_path, resampling_method='sum', scaling_factor=self._scaling_factor_used)
            
            regriders.append(regrider)
            
            print(f"finished converting {slice_indices[index]:.3e}_to_{slice_indices[index+1]:.3e}")
            
        self.y_slice_regriders[f'{number_of_steps}_steps'] = regriders
    
    def _run_stats(self):
        slice_path_array = self.slice_paths[self.complete_slice_keys[self.most_recent_y_slice_index]]
        regriders = self.y_slice_regriders[self.complete_slice_keys[self.most_recent_y_slice_index]]
        output_nc_paths = [regrider._output_nc_path for regrider in regriders]
        slice_indices = self.slice_indices[self.complete_slice_keys[self.most_recent_y_slice_index]]
        
        #save path
        self.y_slice_stats = self.y_slice_conservation_test_results + r"y_slice_stats/"

        # Create the folder
        os.makedirs(self.y_slice_stats)
        
        path_to_stats = []
        stats_dfs = []
    
    # loops each slice
        for index, slice_path in enumerate(slice_path_array):
            
            print(f"starting to take stats for slice: {slice_indices[index]:.3e}_to_{slice_indices[index+1]:.3e}")
            
            slice_name = f"{slice_indices[index]:.3e}_to_{slice_indices[index+1]:.3e}"
            output_nc_path = output_nc_paths[index]
            # lowrez_4326_tif_path = lowrez_4326_tif_paths[index]
            # fulrez_4326_path = fulrez_4326_paths[index]
            slice_type_names = ['Original', 'Output_nc']
            
            slice_data = []
            
            for index_2, file in enumerate([slice_path, output_nc_path]):
                
                xr_array = rioxarray.open_rasterio(file)
                            # Ensure there are no nan-values and replace all possible null values
                try:
                    if '_FillValue' in xr_array.attrs:
                        xr_array = xr_array.where(xr_array != xr_array._FillValue, 0) #remove fill values
                except:
                    pass
                # xr_array = xr_array.where(xr_array != xr_array._FillValue,0) #remove none values so stats can calculate ok
                
                #calculates stats of file
                row_stats = {
                    "Classification":slice_type_names[index_2],
                    "Slice":slice_name,
                    "Total_in_slice":float(xr_array.sum()),
                    "Lowest_lon":float(xr_array.x.min()),
                    "Highest_lon":float(xr_array.x.max()),
                    "Lowest_lat":float(xr_array.y.min()),
                    "Highest_lat":float(xr_array.y.max()),
                    "x_Resolution":np.abs(float(xr_array.x[1]) - float(xr_array.x[2])),
                    "y_Resolution":np.abs(float(xr_array.y[1]) - float(xr_array.y[2]))
                }
                slice_data.append(row_stats)

            slice_df = pd.DataFrame(slice_data)
            print(slice_df)
            stats_dfs.append(slice_df)
            path_to_stats.append(self.y_slice_stats + f"stats_for_slice_{slice_indices[index]:.3e}_to_{slice_indices[index+1]:.3e}.csv")
            slice_df.to_csv(path_to_stats[index])
        
        # create one mega df
        pd.concat(stats_dfs).to_csv(self._output_folder + "conservation_stats")
            
    def run_y_slice_conservation_test(self, number_of_steps = "n_a"):
        
        print('running_y_slice')
        self._y_slice_raster(number_of_steps = number_of_steps)

        print('running_y_slice_regrid')
        self._convert_y_slices()

        print('running_stats')
        self._run_stats()

def plot_netcdf_on_us_map(file_path:str, variable_name:str, title:str, cbar_label:str):
    # Read the netCDF file
    ds = xr.open_dataset(file_path)

    # Extract variable data
    variable_data = ds[variable_name]
    
    # Create a figure and axis with a Lambert Conformal projection
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.LambertConformal()})
    ax.set_extent([-125, -67, 22, 50], ccrs.Geodetic())

    # Plot the variable data on the map
    variable_data.plot(ax=ax, transform=ccrs.PlateCarree(), cmap='viridis', vmin=variable_data.min(), vmax=variable_data.max(), cbar_kwargs={'shrink': 0.5, 'label': cbar_label})
    
    # remove title
    ax.set_title(title)

    # Add coastlines and state borders
    ax.coastlines()
    ax.add_feature(cartopy.feature.STATES, linewidth=0.5, edgecolor='black')

    # Return plot
    return fig, ax

def plot_geotiff(file_path:str, title:str, cbar_label:str):
    # Open the GeoTIFF file
    with rasterio.open(file_path) as src:
        # Read the data and metadata
        data = src.read(1)  # Assuming single-band data, adjust as needed
        extent = [-125, -67, 22, 50]#[src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top]
        
            # Set all values less than 0 to 0
    data = np.where(data < 0, 0, data) #remove fill values
    data = np.where(data == np.nan, 0, data)
    data = np.where(data == np.inf, 0, data)

    # Create a figure and axis with PlateCarree projection
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent)

    # Plot the GeoTIFF data
    img = ax.imshow(data, extent=extent, cmap='viridis')

    # Add coastlines and other features
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black')

    # Add colorbar
    cbar = fig.colorbar(img, ax=ax, shrink=0.5)
    cbar.set_label(cbar_label)
    
    # at title
    ax.set_title(title)
    
    return fig, ax

def _haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    r = 6371000  # Radius of Earth in meters
    return c * r

def _calculate_cell_areas(lon, lat):
    """
    Calculate the area of each grid cell in square meters (This does not work yet)
    """
    dlon = np.radians(np.abs(lon[1] - lon[0]))
    dlat = np.radians(np.abs(lat[1] - lat[0]))
    lat_edges = np.concatenate(([lat[0] - dlat / 2], lat + dlat / 2))
    areas = np.zeros((len(lat), len(lon)))

    for i in range(len(lon)):
        for j in range(len(lat)):
            # Calculate area of cell using Haversine formula
            # This assumes the cell is a rectangle, which is an approximation for small cells
            # The actual area may vary slightly due to the Earth's curvature
            areas[j, i] = (
                _haversine(lon[i] - dlon / 2, lat_edges[j], lon[i] + dlon / 2, lat_edges[j]) *
                _haversine(lon[i], lat_edges[j] - dlat / 2, lon[i], lat_edges[j] + dlat / 2)
            )
    return areas

def build_template_netcdf(output_file, lat_resolution,lon_resolution):
    
    lat = np.arange(-90, 90, lat_resolution)
    lon = np.arange(-180, 180, lon_resolution)
    # time = pd.date_range(start=start_date, end=end_date, freq='MS')

    # Calculate cell areas
    cell_areas = _calculate_cell_areas(lon, lat)

    # Create coordinate variables
    coords = {'lat': lat, 'lon': lon}

    # Create xarray dataset
    ds = xr.Dataset(
        {
            'random': (['time', 'lat', 'lon']),
        },
        coords=coords,
    )

    # Add global attributes
    ds.attrs['title'] = 'Climate Forecast'
    ds.attrs['resolution'] = f'{lon_resolution} + X + {lat_resolution} degrees'

    # Add cell area attribute
    ds['AREA'] = (('lat', 'lon'), cell_areas)

    # Write dataset to NetCDF file
    ds.to_netcdf(output_file)
