#!/bin/bash
#===========================
# Description:
#
# Usage:
#
# History:
#
# Author:
#
#   Yi-Hsuan Chen, yihsuan@umich.edu
# Adpated by Xianglei Huang, xianglei@umich.edu, on August 06, 2019
#===========================

# SET MACHINE TYPE AND FILEPATHS
export USER="hxl"
export MACH="cheyenne"
export email="xianglei@umich.edu"
export CCSMROOT="/glade/u/home/hxl/model/cesm2_1_1_emis"  # source code folder
export RUNDIR="/glade/scratch/$USER"     #
date=`date +%Y%m%d-%H%M`

# SET COMPONENT SET:
# for ~100 options, see: http://www.cesm.ucar.edu/models/cesm2/config/compsets.html
export COMPSET="ETEST"

# SET HORIZONTAL RESOLUTION:
# for supported resolutions, see: http://www.cesm.ucar.edu/models/cesm2/config/grids.html
export RES="f19_g17"  # atm:1.9x2.5  lnd:1.9x2.5

# SET CASE NAME
export CASE="x01-cesm211_emis-${COMPSET}-${RES}-$date"
export CASEROOT="${RUNDIR}/${CASE}"   # case root

# SET surface emissivity file path
#emis_surf="/glade/u/home/yihsuan/data/emis_lw-data/band-type/surface_emissivity_1.9x2.5_RRTMG_53deg.nc"
emis_surf="/glade/u/home/hxl/surface_emissivity_1.9x2.5_RRTMG_53deg.nc"

# if CAM code modification is needed:
#option_code="T"
#option_code="F"
#CODEDIR="/glade/u/home/yihsuan/work/research/cesm111-emis_MC6_IceScat/src.cam/E2000_rrtmg_emis_newTs_mc6_rtr2"

#JOB_QUEUE="economy"
JOB_QUEUE="regular"

# stop options
stop_option="ndays"
#stop_option="nhours"
#stop_option="nmonths"
#stop_option="nyears"
#stop_option="nsteps"

# stop_n
stop_n=1

nhtfrq=-3  # = 0, monthly 
          # > 0, each step

export PROJECT="UCGD0005"  # project code for Cheyenne for this PMWG summer school

#dd=00   # walltime day
hh=2   # walltime hour
mm=00   # walltime minute
JOB_WALLCLOCK_TIME="${hh}:${mm}:00"

option_submit="T" # set "T" if you want to submit the job.

set -x

#################
#################
#################
#
# program start
#
#################
#################
#################

#------------
# check part
#------------

# check required directories
if [ ! -d $CCSMROOT ]; then
  echo "ERROR: CESM root directory [$CCSMROOT] does not exist !!"
  echo "program stop"
  exit 1

#elif [ ! -d $DIN_LOC_ROOT ]; then
#  echo "ERROR: CESM data directory, csmdata, [$DIN_LOC_ROOT] does not exist !!"
#  echo "program stop"
#  exit 1

elif [ ! -d $RUNDIR ]; then 
  echo "ERROR: run directory [$RUNDIR] does not exist!"
  echo "program stop"
  exit 1

elif [ -d $CASEROOT ]; then 
  echo "ERROR: Case root [$CASEROOT] does exist !!"
  echo "please set a new name"
  echo "program stop"
  exit 1
fi

# check option for modifed codes
if [ $option_code ] && [ $option_code == "T" ]; then
  if [ ! -d $CODEDIR ]; then
    echo "ERROR: user modified code directory [$CODEDIR] does not exist !!"
    echo "program stop"
    exit 1
  fi
fi

# get this script name
name0=$0
pp1=`echo $name0 | grep ^/ > /dev/null ; echo $?`

if [ $pp1 -eq 1 ]; then
  path1=`pwd`
  name1=${0##*/}
  script_name="$path1/$name1"

elif [ $pp1 -eq 0 ]; then 
  script_name=$name0

else
  echo "ERROR: fail to retrive this script name"
  echo "program stop"
  exit 1
fi

#--------------------
# create a new case
#---------------------

# CREATE CASE (highest level configuration):
cd $CCSMROOT/cime/scripts || exit 101

# create new case
# add "--run-unsupported" if the grids are not scientifically supported (see it on componenet sets)
./create_newcase --case $CASEROOT \
    --compset $COMPSET \
    --res $RES \
    --project $PROJECT \
    --run-unsupported \
    || exit 105

# copy this script to CASEROOT for future reference
cp -i $script_name "$CASEROOT/00-run_script-$date.sh" && echo "Done. copy running script to CASEROOT" || exit 905

#--------------------
# user modification
#--------------------

cd $CASEROOT || exit 200

#*** user-modified source code ***
if [ $option_code ] && [ $option_code == "T" ]; then
  #/bin/cp ~/cesm1/SourceMods/SurfaceAlbedoMod.F90 $CASEROOT/SourceMods/src.clm/ || 203
  #cp -i $CODEDIR/src.clm/* $CASEROOT/SourceMods/src.clm/ && echo "Done. copy [$CODEDIR] to src.clm" || exit 203
  cp -i $CODEDIR/* $CASEROOT/SourceMods/src.cam/ && echo "Done. copy [$CODEDIR] to src.cam" || exit 203
fi 

#-----------------
# setup the case
#-----------------

# CONFIGURE CASE:

cd $CASEROOT

# CONFIGURE MODEL WITH OPTIONS SET ABOVE:
./case.setup && echo "Done. case setup" || exit 301

# link emissivity file to run directory
ln -s $emis_surf $CASEROOT/run/surface_emissivity_RRTMG_53deg.nc || exit 205

#-------------
# build case
#-------------

#*** setup model run and end time ***
./xmlchange --file env_run.xml --id STOP_OPTION --val $stop_option || exit 303
./xmlchange --file env_run.xml --id STOP_N --val $stop_n || exit 303

# number of times to resubmit:
./xmlchange --file env_run.xml --id RESUBMIT    --val 0 || exit 303

./xmlchange --file env_batch.xml --id JOB_WALLCLOCK_TIME --val $JOB_WALLCLOCK_TIME || exit 303

./xmlchange --file env_batch.xml --id JOB_QUEUE --val $JOB_QUEUE || exit 303

  # write cam namelist
  cat > $CASEROOT/user_nl_cam << EOF1 
&cam_inparm
nhtfrq = $nhtfrq
mfilt  = 12
fincl1 = "FLDSC:A"
EOF1

# BUILD THE MODEL
# qcmd is a command on Cheyenne
qcmd -- ./case.build --skip-provenance-check && echo "Done. Build CESM" || exit 303

#----------------------
# submit job 
#----------------------

if [ $option_submit ] && [ $option_submit == "T" ]; then
  ./case.submit --mail-type all --mail-user $email || exit 910
else
  echo "This job has NOT been submitted yet. Please submit it manually"
fi

exit 0

