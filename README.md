# How to run the modified CESM2

A step-by-step instruction for running the modified CESM2 with the inclusion of surface spectral emissivity on NCAR machines. If you need to run it on your local cluster, you can adapt it accordingly. 

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
