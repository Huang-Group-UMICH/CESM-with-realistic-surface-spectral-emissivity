# How to run the modified CESM2

A step-by-step instruction for running the modified CESM2 with the inclusion of surface spectral emissivity

1. Log into Cheyenne

2. Copy over the entire package of modified CESM2 code to your home directory. “USERNAME” should be replaced with your own real username. 

cp /glade/p/cesm/pcwg/PWS2019_DATA/day3/morning/cesm2_1_1_emis_UM.tar /glade/u/home/USERNAME

3. Go to your directory where you want to put the source code and untar the code package, e.g.

cd /glad/u/home/USERNAME/model 
tar xvf ../cesm2_1_1_emis_UM.tar
 
After it is done, you should see a new subdirectory created under “/glade/u/home/USERNAME/model”, which is cesm2_1_1_emis. 

4. Now copy over the template script for running the model to your home directory, or any directory you are used to save CESM running scripts. For example, 

cp /glade/p/cesm/pcwg/PWS2019_DATA/day3/morning/cesm211_emis_PMWG-ETEST_f19_g17.sh /glade/u/home/USERNAME/

5. This script largely resembles the default CESM2 script. A few things to note is
(1) On the first section, make sure to change following lines to 
export USER="xianglei"   # change to your own username
export email=xianglei@umich.edu # change to your own email
export CCSMROOT="/glade/u/home/hxl/model/cesm2_1_1_emis"
Note this CCSMROOT should be the path above in Step #3 (make sure to include cesm2_1_1_emis)

(2) The emissivity file that we generated is specified as
emis_surf="/glade/p/cesm/pcwg/PWS2019_DATA/day3/morning/surface_emissivity_1.9x2.5_RRTMG_53deg.nc"

(3) The rest configurations are routine, just like standard CESM2 running script. In the current script, it was set up to run for one day with 3-hourly output, with a slab-ocean configuration (ETEST compset). 

(4) The script will compile and submit the job. After the run is finished, you can find output at
The “CASEROOT” directory you specified in the script. e.g. currently it is “/glade/scratch/$USER/archive/$CASE”

For your convenience, I also include a parallel script for run the default CESM2.1.1 under /glade/p/cesm/pcwg/PWS2019_DATA/day3/morning/ as cesm211_standard_PMWG-ETEST_f19_g17.sh
Note: both scripts are set to run for only one day with 3-hourly output. If you use it for your research, you need to adjust this to your own need accordingly.

# How to use AMWG package

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Step-by-step instruction of using AMWG package
Yi-Hsuan Chen
08/05/2019
Xianglei Huang
08/06/2019
based on
Script of running AMWG package: amwg-cesm211-emis_standard.csh or diag140804.csh (provided by AMWG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

***************************
model-to-model comparison, 
  test case VS control case
***************************

IMPORTANT: BEFORE running the script, please make sure to load following modules
	module load nco
	module load ncl

This has to be outside the script and must be done before running the script. Otherwise, the script cannot produce results.
Note: you can type "module list" on command line to see all loaded modules.

Things you need to be aware in amwg-cesm211-emis_standard_PMWG.csh
1. set working director; replace "USER = hxl" with your own username

        %% set $wrkdir that the AMWG plots will be saved there
	set USER = hxl
	set scratchdir = /glade/scratch/$USER/
	set wrkdir = $scratchdir/y01-amwg-$temp/

2. Set the identifying casename and paths for the test case run. The output files MUST BE "$test_path_history/$test_casename.cam.h0.YYYY-MM.nc"
Note: for this practice, the output file is store in $test_path_history and you do not need to change it. If you use this script later with your own output file, you need to change this to the right path.

	set test_casename  = u02-cesm211_emis-ETEST-f19_g17
	set test_filetype = monthly_history
	#set test_filetype = time_series

	set test_path_history =  /glade/scratch/yihsuan/archive/$test_casename/atm/hist/
	set test_path_climo   =  $wrkdir/amwg_climo/$test_casename/
	set test_path_diag    =  $wrkdir
	set test_path_HPSS    =  /CCSM/csm/${test_casename}/atm/hist/

2. Set the identifying casename and paths for your control case run. The output files MUST BE "$cntl_path_history/$cntl_casename.cam.h0.YYYY-MM.nc"
Note: for this practice, the output file is store in $cntl_path_history and you do not need to change i
t. If you use this script later with your own output file, you need to change this to the right path.

	set cntl_casename   =  u01-cesm211_standard-ETEST-f19_g17
	set cntl_filetype = monthly_history

	set cntl_path_history = /glade/scratch/yihsuan/archive/$cntl_casename/atm/hist/
	set cntl_path_climo   = $wrkdir/amwg_climo/$cntl_casename/
	set cntl_path_HPSS    = /CCSM/csm/${cntl_casename}/atm/hist/

3. Turn on/off the computation of climatologies. If the climatologies are already in $test_path_climo and $cntl_path_climo, set to 1(=OFF).

	set test_compute_climo = 0  # (0=ON,1=OFF) 
	set cntl_compute_climo = 0    # (0=ON,1=OFF) 

   If computing climatological means for test/cntl case, specify the first year of your data, and the number of years of data to be used.

	set test_first_yr = 6           # first year (must be >= 1)
	set test_nyrs     = 5           # number of yrs (must be >= 1)

	set cntl_first_yr = $test_first_yr        # first year (must be >= 1)
	set cntl_nyrs     = $test_nyrs        # number of yrs (must be >= 1)

Note: for this practice purpose, to save year, we only compare differences for 5 years. The output here is actually for year 1-35. 

4. Select the diagnostic sets to be done. You can do one at a time or as many as you want at one time, or all at once.

	set all_sets = 0  # (0=ON,1=OFF)  Do all the CAM sets (1-16)
	set set_1  = 1    # (0=ON,1=OFF)  tables of global,regional means
	set set_2  = 1    # (0=ON,1=OFF)  implied transport plots 
	...

5. Use custom case names for the PLOTS instead of the case names encoded in the netcdf files (default). 

	set custom_names = 0     # (0=ON,1=OFF)
	set test_name = cesm211_emis              # test case name 
	set cntl_name = cesm211_standard              # control case name

6. (optional) Compute whether the means of the test case and control case are significantly different from each other at each grid point.

	set significance = 1         # (0=ON,1=OFF)
	set sig_lvl = 0.05           # level of significance
Note: 5 years are indeed not enough to give robust signficance

7. (optional: need this if you run on other machines) Set amwg diagnostic package root location 

	%% CSIL machines (geyser, caldeira, ...)
	setenv DIAG_HOME /glade/p/cesm/amwg/amwg_diagnostics

8. Change other settings in the script as necessary

9. execute the AMWG package script
	> ./amwg-cesm211-emis_standard.csh

   Once the execution finished, a tar file "$test_casename-$cntl_casename.tar" would be in $test_path_diag.

10. View AMWG plots
	untar the tar file and then a folder "$test_casename-$cntl_casename" will be created.
	> tar -xvf $test_casename-$cntl_casename.tar

	View AMWG plots
	> cd $test_casename-$cntl_casename

	Open index.html in any web browser and you will see all AMWG plots.

