from regrid_tools import *
print("running")
import numpy as np

destination_folder_path = r'/uufs/chpc.utah.edu/common/home/u1285966/Haskins/tclark/regrid-help-most-recent/F_playa/'
template_path = r'regrid-help-most-recent/template_files/dust_emissions_05.20210906.nc'
input_5070_path = r"/uufs/chpc.utah.edu/common/home/u1285966/Haskins/tclark/regrid-help-most-recent/F_playa/us.tif"


regrider = Organized_regriding(destination_folder_path, template_path, descriptive_name=r'f_playa')

# post-process function
def playa_mask (average_concentration_in_grid_box:np.array):
    threshold = 10
    return (average_concentration_in_grid_box > threshold).astype(np.uint8)


regrider.apply_array_function_to_tif(input_5070_path,'f_playa', playa_mask)


regrider.tif_5070_to_gc_netcdf(input_tif_5070= regrider.processed_tifs['f_playa'],
                               resampling_method='average', scaling_factor= 0.01)



regrider.run_y_slice_conservation_test(number_of_steps = 1)