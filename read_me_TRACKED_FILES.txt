<<<<<<< HEAD
I cloned the regridding.git repo to this local directory. 
The purpose of this directory is to track the changes I made to Todd's f_playa code without editing his existing files.
I only am tracking create_f_plya.py and regrid_tools.py because the original folder of all Todd's scripts was way too big to upload to git.
Also, most of the files in Todd's regrid-help-most-recent folder contained lots of files and directories that are not used.
=======
This folder contains Todd's original files in his original organizational structure for regridding soil conductivity data (f_playa). 
I (Joey) am currently tracking create_f_plya.py, regrid_tools.py, and template files in the repo regridding.git. These are the main files involved with f_playa, though I may add more files to the repo later.
I have cloned the regridding.git repo to the regridding_tracked.git repo. I did this to avoid making any changes to Todd's original files. Instead, I am tracking all the changes in this cloned repo.
If any changes are made to the files included in the regridding.git repo, you need to follow these commands to make sure they are also included in the cloned regridding_tracked.git repo

cd /uufs/chpc.utah.edu/common/home/haskins-group1/users/jbail/regridding_tracked
git remote add old-origin https://github.com/joeybail96/regridding.git
git fetch old-origin
git merge old-origin/main -m "Your custom merge message here"
git push origin main

This adds the original repository (regridding.git) as a remote named old-origin, fetches the latest changes, and merges them into your current branch (regridding_tracked.git).
>>>>>>> old-origin/main
