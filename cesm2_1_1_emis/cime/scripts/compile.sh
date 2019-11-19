#! /bin/csh -fe
### This script was created 2015-11-15 by Philip Cameron-Smith (pjc@llnl.gov) and Peter Caldwell
### and incorporates some features originally from Hui Wan, Kai Zhang, and Balwinder Singh.
### Significant improvements from Michael Deakin and Chris Golaz.
###

###===================================================================
### THINGS USERS USUALLY CHANGE (SEE END OF SECTION FOR GUIDANCE)
###===================================================================
source /etc/profile.d/z00Lmod.csh
module purge
module load pgi/16.4
module load openmpi/1.10.2/pgi/16.4
module load hdf5/1.8.16/pgi/16.4
module load netcdf/4.4.1/pgi/16.4

setenv NETCDFROOT /sw/arcts/centos7/netcdf/4.4.1/pgi-16.4-hdf5-1.8.16/
setenv PATH $NETCDFROOT/bin:$PATH
setenv LD_LIBRARY_PATH $NETCDFROOT/lib:$LD_LIBRARY_PATH
setenv NETCDF_INCLUDES $NETCDFROOT/include
setenv NETCDF_LIBS $NETCDFROOT/lib
setenv MPICH_PATH /sw/arcts/centos7/openmpi/1.10.2-pgi-16.4/
setenv PNETCDFROOT $NETCDFROOT

setenv USER_FC pgi

./create_newcase --case ../../../cases/I1850Clm50Sp --compset I1850Clm50Sp --res f19_g17 --machine flux --walltime 50:00:00 -q xianglei_flux --run-unsupported
cd ../../../cases/I1850Clm50Sp
./case.setup 
./case.build
#./xmlchange --file env_mach_pes.xml --id NTASKS_ATM --val 16
#./xmlchange --file env_mach_pes.xml --id NTASKS_LND --val 16
#./xmlchange --file env_mach_pes.xml --id NTASKS_ICE --val 16
#./xmlchange --file env_mach_pes.xml --id NTASKS_OCN --val 16
#./xmlchange --file env_mach_pes.xml --id NTASKS_CPL --val 16
#./xmlchange --file env_mach_pes.xml --id NTASKS_ROF --val 16
#./xmlchange --file env_mach_pes.xml --id NTASKS_GLC --val 16   # NB: Set to 1 if using GLC compset 

./case.submit
