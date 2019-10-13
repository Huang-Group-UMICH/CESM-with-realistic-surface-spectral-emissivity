module atm_comp_mct

  use pio              , only: file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim, &
                               pio_put_att, pio_enddef, pio_initdecomp, pio_read_darray, pio_freedecomp, &
                               pio_closefile, pio_write_darray, pio_def_var, pio_inq_varid, &
	                       pio_noerr, pio_bcast_error, pio_internal_error, pio_seterrorhandling
  use mct_mod
  use seq_cdata_mod
  use esmf

  use seq_comm_mct,      only: num_inst_atm
  use seq_flds_mod
  use seq_infodata_mod
  use seq_timemgr_mod

  use shr_kind_mod     , only: r8 => shr_kind_r8, cs=>shr_kind_cs, cl=>shr_kind_cl
  use shr_file_mod     , only: shr_file_getunit,                          &
                               shr_file_setLogUnit, shr_file_setLogLevel, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
		               shr_file_setIO
  use shr_sys_mod      , only: shr_sys_flush, shr_sys_abort

  use cam_cpl_indices
  use atm_import_export
  use cam_comp,          only: cam_init, cam_run1, cam_run2, cam_run3, cam_run4, cam_final
  use cam_instance     , only: cam_instance_init, inst_suffix, inst_index
  use cam_control_mod  , only: initial_run, dart_mode, cam_ctrl_set_orbit
  use radiation        , only: radiation_nextsw_cday
  use phys_grid        , only: get_ncols_p, ngcols, get_gcol_p, get_rlat_all_p, &
	                       get_rlon_all_p, get_area_all_p
  use ppgrid           , only: pcols, begchunk, endchunk
  use dyn_grid         , only: get_horiz_grid_dim_d
  use camsrfexch       , only: cam_out_t, cam_in_t
  use cam_initfiles    , only: cam_initfiles_get_caseid, cam_initfiles_get_restdir
  use cam_abortutils   , only: endrun
  use filenames        , only: interpret_filename_spec
  use spmd_utils       , only: spmdinit, masterproc, iam
  use time_manager     , only: get_curr_calday, advance_timestep, get_curr_date, get_nstep, &
                               get_step_size
  use ioFileMod
  use perf_mod
  use cam_logfile      , only: iulog

  implicit none
  save
  private

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: atm_setgsmap_mct
  private :: atm_domain_mct
  private :: atm_read_srfrest_mct
  private :: atm_write_srfrest_mct

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(cam_in_t) , pointer :: cam_in(:)
  type(cam_out_t), pointer :: cam_out(:)

  ! Filename specifier for restart surface file
  character(len=cl) :: rsfilename_spec_cam

  integer,                 pointer :: dof(:) ! needed for pio_init decomp for restarts
  type(seq_infodata_type), pointer :: infodata

!================================================================================
CONTAINS
!================================================================================

  subroutine atm_init_mct( EClock, cdata_a, x2a_a, a2x_a, NLFilename )

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(inout)              :: EClock
    type(seq_cdata), intent(inout)              :: cdata_a
    type(mct_aVect), intent(inout)              :: x2a_a
    type(mct_aVect), intent(inout)              :: a2x_a
    character(len=*), optional,   intent(IN)    :: NLFilename ! Namelist filename
    !
    ! Locals
    !
    type(mct_gsMap), pointer   :: gsMap_atm
    type(mct_gGrid), pointer   :: dom_a
    integer :: ATMID
    integer :: mpicom_atm
    integer :: lsize

    logical :: first_time = .true.
    logical :: exists

    integer :: shrlogunit, shrloglev ! save values, restore on return

    character(len=cs) :: starttype ! infodata start type
    character(len=cl) :: caseid    ! case ID
    character(len=cl) :: ctitle    ! case title
    character(len=cl) :: model_doi_url ! DOI for CESM model run

    logical :: aqua_planet       ! Flag to run model in "aqua planet" mode
    logical :: brnch_retain_casename ! true => branch run may use same caseid as
                                     !         the run being branched from
    logical :: single_column
    real(r8):: scmlat
    real(r8):: scmlon

    real(r8) :: eccen
    real(r8) :: obliqr
    real(r8) :: lambm0
    real(r8) :: mvelpp

    logical :: perpetual_run    ! If in perpetual mode or not
    integer :: perpetual_ymd    ! Perpetual date (YYYYMMDD)

    logical :: dart_mode_in
    character(len=cl) :: atm_resume_all_inst(num_inst_atm) ! atm resume file

    real(r8):: nextsw_cday      ! calendar of next atm shortwave
    integer :: stepno           ! time step

    integer :: dtime            ! time step increment (sec)
    integer :: atm_cpl_dt       ! driver atm coupling time step
    integer :: nstep            ! CAM nstep

    real(r8):: caldayp1         ! CAM calendar day for for next cam time step

    integer :: lbnum

    integer :: hdim1_d, hdim2_d ! dimensions of rectangular horizontal grid
                                ! data structure, If 1D data structure, then
                                ! hdim2_d == 1.
    !-----------------------------------------------------------------------
    !
    ! Determine cdata points
    !
    call seq_cdata_setptrs(cdata_a, ID=ATMID, mpicom=mpicom_atm, &
         gsMap=gsMap_atm, dom=dom_a, infodata=infodata)

    if (first_time) then

       call cam_instance_init(ATMID)

       ! Set filename specifier for restart surface file
       ! (%c=caseid, $y=year, $m=month, $d=day, $s=seconds in day)
       rsfilename_spec_cam = '%c.cam' // trim(inst_suffix) // '.rs.%y-%m-%d-%s.nc'

       ! Determine attribute vector indices

       call cam_cpl_indices_set()

       ! Initialize atm use of MPI
       call spmdinit(mpicom_atm)

#if (defined _MEMTRACE)
       if(masterproc) then
          lbnum=1
          call memmon_dump_fort('memmon.out','atm_init_mct:start::',lbnum)
       endif
#endif

       ! Redirect share output to cam log

       if (masterproc) then
          inquire(file='atm_modelio.nml'//trim(inst_suffix), exist=exists)
          if (exists) then
             iulog = shr_file_getUnit()
             call shr_file_setIO('atm_modelio.nml'//trim(inst_suffix), iulog)
          endif

          write(iulog,*) "CAM atmosphere model initialization"
       endif

       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (iulog)
       !
       ! Get data from infodata object
       !
       call seq_infodata_GetData( infodata,                                           &
            case_name=caseid, case_desc=ctitle, model_doi_url=model_doi_url,          &
            start_type=starttype,                                                     &
            atm_resume=atm_resume_all_inst,                                           &
            aqua_planet=aqua_planet,                                                  &
            brnch_retain_casename=brnch_retain_casename,                              &
            single_column=single_column, scmlat=scmlat, scmlon=scmlon,                &
            orb_eccen=eccen, orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr, &
            perpetual=perpetual_run, perpetual_ymd=perpetual_ymd)

       dart_mode_in = .false.
       if (trim(atm_resume_all_inst(MIN(num_inst_atm,inst_index))) == 'TRUE') dart_mode_in = .true.

       ! Initialize CAM, allocate cam_in and cam_out and determine
       ! atm decomposition (needed to initialize gsmap)
       ! for an initial run, cam_in and cam_out are allocated in cam_init
       ! for a restart/branch run, cam_in and cam_out are allocated in restart
       !
       call cam_init(EClock, &
          caseid, ctitle, starttype, dart_mode_in,     &
          brnch_retain_casename, aqua_planet,          &
          single_column, scmlat, scmlon,               &
          eccen, obliqr, lambm0, mvelpp,               &
          perpetual_run, perpetual_ymd, model_doi_url, &
          cam_out, cam_in)

       !
       ! Initialize MCT gsMap, domain and attribute vectors (and dof)
       !
       call atm_SetgsMap_mct( mpicom_atm, ATMID, gsMap_atm )
       lsize = mct_gsMap_lsize(gsMap_atm, mpicom_atm)

       ! Set dof (module variable, needed for pio for restarts)
       call mct_gsmap_orderedpoints(gsmap_atm, iam, dof)
       !
       ! Initialize MCT domain
       !
       call atm_domain_mct( lsize, gsMap_atm, dom_a )
       !
       ! Initialize MCT attribute vectors
       !
       call mct_aVect_init(a2x_a, rList=seq_flds_a2x_fields, lsize=lsize)
       call mct_aVect_zero(a2x_a)

       call mct_aVect_init(x2a_a, rList=seq_flds_x2a_fields, lsize=lsize)
       call mct_aVect_zero(x2a_a)
       !
       ! Create initial atm export state
       !
       call atm_export( cam_out, a2x_a%rattr )
       !
       ! Set flag to specify that an extra albedo calculation is to be done (i.e. specify active)
       !
       call seq_infodata_PutData(infodata, atm_prognostic=.true.)
       call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
       call seq_infodata_PutData(infodata, atm_nx=hdim1_d, atm_ny=hdim2_d)

       ! Set flag to indicate that CAM will provide carbon and dust deposition fluxes.
       ! This is now hardcoded to .true. since the ability of CICE to read these
       ! fluxes from a file has been removed.
       call seq_infodata_PutData(infodata, atm_aero=.true.)

       !
       ! Set time step of radiation computation as the current calday
       ! This will only be used on the first timestep of an initial run
       !
       if (initial_run) then
          nextsw_cday = get_curr_calday()
          call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday )
       end if

       ! End redirection of share output to cam log

       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)

       first_time = .false.

    else

       ! For initial run, run cam radiation/clouds and return
       ! For restart run, read restart x2a_a
       ! Note - a2x_a is computed upon the completion of the previous run - cam_run1 is called
       ! only for the purposes of finishing the flux averaged calculation to compute a2x_a
       ! Note - cam_run1 is called on restart only to have cam internal state consistent with the
       ! a2x_a state sent to the coupler

       ! Redirect share output to cam log

       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (iulog)

       call seq_timemgr_EClockGetData(EClock, StepNo=StepNo)
       if (StepNo == 0) then
          call atm_import( x2a_a%rattr, cam_in )
          call cam_run1 ( cam_in, cam_out )
          call atm_export( cam_out, a2x_a%rattr )
       else
          call atm_read_srfrest_mct( EClock, x2a_a, a2x_a )
          call atm_import( x2a_a%rattr, cam_in, restart_init=.true. )
          call cam_run1 ( cam_in, cam_out )
       end if

       ! Compute time of next radiation computation, like in run method for exact restart

       call seq_timemgr_EClockGetData(Eclock,dtime=atm_cpl_dt)
       dtime = get_step_size()
       nstep = get_nstep()
       if (nstep < 1 .or. dtime < atm_cpl_dt) then
          nextsw_cday = radiation_nextsw_cday()
       else if (dtime == atm_cpl_dt) then
          caldayp1 = get_curr_calday(offset=int(dtime))
          nextsw_cday = radiation_nextsw_cday()
          if (caldayp1 /= nextsw_cday) nextsw_cday = -1._r8
       else
          call shr_sys_abort('dtime must be less than or equal to atm_cpl_dt')
       end if
       call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday )

       ! End redirection of share output to cam log

       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)

    end if

#if (defined _MEMTRACE )
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out','atm_init_mct:end::',lbnum)
      call memmon_reset_addr()
    endif
#endif

    call shr_sys_flush(iulog)

 end subroutine atm_init_mct

 !================================================================================

 subroutine atm_run_mct( EClock, cdata_a, x2a_a, a2x_a)
    use parrrtm,        only: nbndlw  ! Xiuhong Chen added on June26, 2019


    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a
    !
    ! Local variables
    !
    integer :: shrlogunit, shrloglev ! save values, restore on exit

    real(r8) :: eccen
    real(r8) :: obliqr
    real(r8) :: lambm0
    real(r8) :: mvelpp

    logical :: dosend          ! true => send data back to driver
    integer :: dtime           ! time step increment (sec)
    integer :: atm_cpl_dt      ! driver atm coupling time step
    integer :: ymd_sync        ! Sync date (YYYYMMDD)
    integer :: yr_sync         ! Sync current year
    integer :: mon_sync        ! Sync current month
    integer :: day_sync        ! Sync current day
    integer :: tod_sync        ! Sync current time of day (sec)
    integer :: ymd             ! CAM current date (YYYYMMDD)
    integer :: yr              ! CAM current year
    integer :: mon             ! CAM current month
    integer :: day             ! CAM current day
    integer :: tod             ! CAM current time of day (sec)

    real(r8):: caldayp1        ! CAM calendar day for for next cam time step
    real(r8):: nextsw_cday     ! calendar of next atm shortwave
    logical :: rstwr           ! .true. ==> write restart file before returning
    logical :: nlend           ! Flag signaling last time-step
    logical :: rstwr_sync      ! .true. ==> write restart file before returning
    logical :: nlend_sync      ! Flag signaling last time-step
    logical :: first_time = .true.
    integer :: lbnum
    character(len=*), parameter :: subname="atm_run_mct"
    !-----------------------------------------------------------------------
   !!!!!!!!!!!!  Added on Nov.10, 2018 by Xiuhong Chen

    integer :: i,j,c,ncols,sizebuf

    real(r8) :: Ts_LW              ! surface temperature derived from longwave upward flux
    real(r8) :: v1_rrtmg_lw(16 + 1)   ! RRTMG_LW band edges
    !integer ::ilats(pcols),ilons(pcols)
    real(r8) :: lats(pcols)           ! array of chunk latitudes
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real  :: pi
    real :: emis0(pcols*(endchunk-begchunk+1),nbndlw)
    real ::  desert_emis(nbndlw), water_emis(nbndlw),ice_emis(nbndlw),grass_emis(nbndlw), snow_emis(nbndlw)
!    integer
!    ::ilats2(pcols*(endchunk-begchunk+1)),ilons2(pcols*(endchunk-begchunk+1))
    real(r8) ::lats2(pcols*(endchunk-begchunk+1)),lons2(pcols*(endchunk-begchunk+1))
    data v1_rrtmg_lw /10. , 350., 500., 630., 700., 820., 980., 1080.,1180.,1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250/
    data water_emis / 0.8524,  0.9080,  0.9098,  0.9223,  0.9451,0.9780,0.9752, 0.9705, &
            0.9662, 0.9626, 0.9624,  0.9638,  0.9623, 0.9623,  0.9623, 0.9623 /
    data desert_emis /0.9085,  0.9178,  0.8747,  0.9458,  0.9488,0.9452,0.8912, &
              0.6858, 0.8513, 0.9700, 0.9504,  0.9396,  0.9373, 0.9373,0.9373,0.9373 /
!     data snow_emis / 0.9963, 0.9875, 0.9859, 0.9809, 0.9759,  0.9871,  0.9926,     0.9892, &
!     0.9862, 0.9842, 0.9838, 0.9751, 0.9720, 0.9720, 0.9720, 0.9720/  ! fine     snow
    ! data snow_emis /0.9933, 0.9905, 0.9799,  0.9717, 0.9643, 0.9819,0.9910,&
    ! 0.9863,  0.9812,  0.9776, 0.9771,  0.9721, 0.9695,  0.9695, 0.9695,
    ! 0.9695/ !median snow
    data snow_emis / 0.9917, 0.9887, 0.9771, 0.9680, 0.9597,  0.9796,0.9898,0.9845,  &
       0.9787, 0.9747, 0.9742, 0.9682, 0.9654, 0.9654, 0.9654, 0.9654/ ! coarse snow
    data ice_emis / 0.8819, 0.9517,  0.9308,  0.9197, 0.9107, 0.9534,  0.9768,0.9690, &
       0.9636, 0.9614, 0.9639, 0.9625, 0.9609, 0.9609, 0.9609, 0.9609/ ! ice
    data grass_emis / 0.9872, 0.9872, 0.9872, 0.9872, 0.9864, 0.9825,  0.9827,&
         0.9808,  0.9886, 0.9887, 0.9888, 0.9865, 0.9863, 0.9863, 0.9863,0.9863/
!!!!!!!!!!!!

#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out',SubName //':start::',lbnum)
    endif
#endif

    ! Redirect share output to cam log

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Note that sync clock time should match cam time at end of time step/loop not beginning

    call seq_timemgr_EClockGetData(EClock,curr_ymd=ymd_sync,curr_tod=tod_sync, &
       curr_yr=yr_sync,curr_mon=mon_sync,curr_day=day_sync)

    nlend_sync = seq_timemgr_StopAlarmIsOn(EClock)
    rstwr_sync = seq_timemgr_RestartAlarmIsOn(EClock)

    ! load orbital parameters
    call seq_infodata_GetData( infodata,                                           &
       orb_eccen=eccen, orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr)
    call cam_ctrl_set_orbit(eccen, obliqr, lambm0, mvelpp)

    ! Map input from mct to cam data structure

    call t_startf ('CAM_import')
    call atm_import( x2a_a%rattr, cam_in )
    call t_stopf  ('CAM_import')
    !   Xiuhong Chen added here, June26, 2019 !!!!!!!!!!!!!!!!
    call get_curr_date(yr, mon, day, tod)
    pi = 4.*atan(1.)
     sizebuf=0
     do c=begchunk, endchunk
             ncols = get_ncols_p(c)
              call get_rlat_all_p(c, ncols, lats)
              call get_rlon_all_p(c, ncols, lons)
             do i=1,ncols
                sizebuf = sizebuf+1
                 lats2(sizebuf) = lats(i)*180/pi
                 if (lons(i).lt.0) then
                    lons2(sizebuf) = lons(i)*180/pi + 180
                 else
                    lons2(sizebuf) = lons(i)*180/pi
                 endif
             enddo
     enddo
    call read_surface_emis(sizebuf,lats2,lons2,mon,emis0(1:sizebuf,:))
         sizebuf=0
         do c=begchunk, endchunk
             ncols = get_ncols_p(c)
             do i=1,ncols
                  sizebuf = sizebuf + 1
                 if (cam_in(c)%landfrac(i).gt. 0.99 .and. cam_in(c)%icefrac(i).lt. 0.01 .and. (cam_in(c)%snowhland(i) + cam_in(c)%snowhice(i)).lt. 0.001) then
                     ! if original emissivity over band 1080-1180 cm-1, ie.
                     ! cam_in%srf_emis_spec(i,8) is
                     ! smaller than 0.8, then the original surface type
                     ! of this grid is desert, otherwise is non-desert
                      cam_in(c)%srf_emis_spec(i,:)= emis0(sizebuf,:)

                     ! if the original surface type is non-desert but LAI is
                     ! smaller than 0.001, change to desert type
                      if (cam_in(c)%tlai(i).lt. 0.001 .and.emis0(sizebuf,8)>0.8) then
                        cam_in(c)%srf_emis_spec(i,:)  = desert_emis
                      endif
                     ! if the orignal surface type is desert but LAI larger than
                     ! 2, change to grass type
                      if (cam_in(c)%tlai(i).gt. 2 .and. emis0(sizebuf,8)<0.8)then
                        cam_in(c)%srf_emis_spec(i,:)  = grass_emis
                      endif
             !  cam_in(c)%srf_emis_spec(i,:)  = 1.0

           else
   cam_in(c)%srf_emis_spec(i,:)  =  emis0(sizebuf,:) * cam_in(c)%landfrac(i) + ice_emis * cam_in(c)%icefrac(i) + water_emis * cam_in(c)%ocnfrac(i)
                 endif


                call  get_Ts_from_LW_emis(v1_rrtmg_lw,real(cam_in(c)%srf_emis_spec(i,:)),cam_in(c)%lwup(i),cam_out(c)%flwds_spec(i,:), 16 + 1, 3, Ts_LW)
                cam_in(c)%ts_atm(i) =Ts_LW

               cam_in(c)%ts(i) =  cam_in(c)%ts_atm(i)
            enddo
    enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!end of add!!!!!!!!!!!!!!!




    ! Cycle over all time steps in the atm coupling interval

    dosend = .false.
    do while (.not. dosend)

       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated

       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))

       ! Determine if time to write cam restart and stop

       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.

       ! Run CAM (run2, run3, run4)

       call t_startf ('CAM_run2')
       call cam_run2( cam_out, cam_in )
       call t_stopf  ('CAM_run2')

       call t_startf ('CAM_run3')
       call cam_run3( cam_out )
       call t_stopf  ('CAM_run3')

       call t_startf ('CAM_run4')
       call cam_run4( cam_out, cam_in, rstwr, nlend, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
       call t_stopf  ('CAM_run4')

       ! Advance cam time step

       call t_startf ('CAM_adv_timestep')
       call advance_timestep()
       call t_stopf  ('CAM_adv_timestep')

       ! Run cam radiation/clouds (run1)

       call t_startf ('CAM_run1')
       call cam_run1 ( cam_in, cam_out )
       call t_stopf  ('CAM_run1')

       ! Map output from cam to mct data structures

       call t_startf ('CAM_export')
       call atm_export( cam_out, a2x_a%rattr )
       call t_stopf ('CAM_export')

    end do

    ! Get time of next radiation calculation - albedos will need to be
    ! calculated by each surface model at this time

    call seq_timemgr_EClockGetData(Eclock,dtime=atm_cpl_dt)
    dtime = get_step_size()
    if (dtime < atm_cpl_dt) then
       nextsw_cday = radiation_nextsw_cday()
    else if (dtime == atm_cpl_dt) then
       caldayp1 = get_curr_calday(offset=int(dtime))
       nextsw_cday = radiation_nextsw_cday()
       if (caldayp1 /= nextsw_cday) nextsw_cday = -1._r8
    else
       call shr_sys_abort('dtime must be less than or equal to atm_cpl_dt')
    end if

    call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday )

    ! Write merged surface data restart file if appropriate

    if (rstwr_sync) then
       call atm_write_srfrest_mct( x2a_a, a2x_a, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
    end if

    ! Check for consistency of internal cam clock with master sync clock

    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EClockGetData(EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' cam ymd=',ymd     ,'  cam tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( subname//': CAM clock is not in sync with master Sync Clock' )
    end if

    ! End redirection of share output to cam log

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out',SubName //':end::',lbnum)
      call memmon_reset_addr()
    endif
#endif

  end subroutine atm_run_mct

  !================================================================================

  subroutine atm_final_mct( EClock, cdata_a, x2a_a, a2x_a)

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a

    call cam_final( cam_out, cam_in )

  end subroutine atm_final_mct

  !================================================================================

  subroutine atm_SetgsMap_mct( mpicom_atm, ATMID, GSMap_atm )

    !-------------------------------------------------------------------
    use phys_grid, only : get_nlcols_p
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_atm
    integer        , intent(in)  :: ATMID
    type(mct_gsMap), intent(out) :: GSMap_atm
    !
    ! Local variables
    !
    integer, allocatable :: gindex(:)
    integer :: i, n, c, ncols, sizebuf, nlcols
    integer :: ier            ! error status
    !-------------------------------------------------------------------

    ! Build the atmosphere grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP

    ! Determine global seg map

    sizebuf=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          sizebuf = sizebuf+1
       end do
    end do

    allocate(gindex(sizebuf))

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          n=n+1
          gindex(n) = get_gcol_p(c,i)
       end do
    end do

    nlcols = get_nlcols_p()
    call mct_gsMap_init( gsMap_atm, gindex, mpicom_atm, ATMID, nlcols, ngcols)

    deallocate(gindex)

  end subroutine atm_SetgsMap_mct

  !===============================================================================

  subroutine atm_domain_mct( lsize, gsMap_a, dom_a )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)   :: lsize
    type(mct_gsMap), intent(in)   :: gsMap_a
    type(mct_ggrid), intent(inout):: dom_a
    !
    ! Local Variables
    !
    integer  :: n,i,c,ncols           ! indices
    real(r8) :: lats(pcols)           ! array of chunk latitudes
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real(r8) :: area(pcols)           ! area in radians squared for each grid point
    real(r8), pointer  :: data(:)     ! temporary
    integer , pointer  :: idata(:)    ! temporary
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI
    !-------------------------------------------------------------------
    !
    ! Initialize mct atm domain
    !
    call mct_gGrid_init( GGrid=dom_a, CoordChars=trim(seq_flds_dom_coord), OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Initialize attribute vector with special value
    !
    call mct_gsMap_orderedPoints(gsMap_a, iam, idata)
    call mct_gGrid_importIAttr(dom_a,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(dom_a,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_a,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_a,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_a,"aream",data,lsize)
    data(:) = 0.0_R8
    call mct_gGrid_importRAttr(dom_a,"mask" ,data,lsize)
    data(:) = 1.0_R8
    call mct_gGrid_importRAttr(dom_a,"frac" ,data,lsize)
    !
    ! Fill in correct values for domain components
    !
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       do i=1,ncols
          n = n+1
          data(n) = lats(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"lat",data,lsize)

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          n = n+1
          data(n) = lons(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"lon",data,lsize)

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_area_all_p(c, ncols, area)
       do i=1,ncols
          n = n+1
          data(n) = area(i)
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"area",data,lsize)

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          data(n) = 1._r8 ! mask
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"mask"   ,data,lsize)
    deallocate(data)

  end subroutine atm_domain_mct

  !===========================================================================================

  subroutine atm_read_srfrest_mct( EClock, x2a_a, a2x_a)

    !-----------------------------------------------------------------------
    use cam_pio_utils, only: cam_pio_openfile, cam_pio_closefile, pio_subsystem
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(inout) :: EClock
    type(mct_aVect), intent(inout) :: x2a_a
    type(mct_aVect), intent(inout) :: a2x_a

    ! Local variables
    character(len=cl)     :: fname_srf_cam  ! surface restart filename
    character(len=cl)     :: pname_srf_cam  ! surface restart full pathname
    integer               :: rcode        ! return error code
    integer               :: yr_spec      ! Current year
    integer               :: mon_spec     ! Current month
    integer               :: day_spec     ! Current day
    integer               :: sec_spec     ! Current time of day (sec)
    integer               :: nf_x2a, nf_a2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t)     :: file
    type(io_desc_t)       :: iodesc
    type(var_desc_t)      :: varid
    character(CL)         :: itemc       ! string converted to char
    type(mct_string)      :: mstring     ! mct char type
    !-----------------------------------------------------------------------

    ! Determine and open surface restart dataset

    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec )

    if (dart_mode) then
       fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, &
          yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec=sec_spec )
       pname_srf_cam = fname_srf_cam
    else
       fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, &
          case=cam_initfiles_get_caseid(), yr_spec=yr_spec, &
          mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
       pname_srf_cam = trim(cam_initfiles_get_restdir() )//fname_srf_cam
    end if
    call getfil(pname_srf_cam, fname_srf_cam)

    call cam_pio_openfile(File, fname_srf_cam, 0)
    call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), dof, iodesc)
    allocate(tmp(size(dof)))

    nf_x2a = mct_aVect_nRattr(x2a_a)
    do k=1,nf_x2a
       call mct_aVect_getRList(mstring,k,x2a_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2a_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2a_a%rattr(k,:) = tmp(:)
       else
	  if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2a_a%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_a2x = mct_aVect_nRattr(a2x_a)
    do k=1,nf_a2x
       call mct_aVect_getRList(mstring,k,a2x_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'a2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       a2x_a%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call cam_pio_closefile(File)
    deallocate(tmp)

  end subroutine atm_read_srfrest_mct

  !===========================================================================================

  subroutine atm_write_srfrest_mct( x2a_a, a2x_a, &
       yr_spec, mon_spec, day_spec, sec_spec)

    !-----------------------------------------------------------------------
    use cam_pio_utils,       only: cam_pio_createfile, cam_pio_closefile, pio_subsystem
    use cam_history_support, only: fillvalue
    !
    ! Arguments
    !
    type(mct_aVect), intent(in) :: x2a_a
    type(mct_aVect), intent(in) :: a2x_a
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    character(len=cl)         :: fname_srf_cam  ! surface restart filename
    integer                   :: rcode        ! return error code
    integer                   :: nf_x2a, nf_a2x, dimid(1), k
    type(file_desc_t)         :: file
    type(var_desc_t), pointer :: varid_x2a(:), varid_a2x(:)
    type(io_desc_t)           :: iodesc
    character(CL)             :: itemc       ! string converted to char
    type(mct_string)          :: mstring     ! mct char type
    !-----------------------------------------------------------------------

    ! Determine and open surface restart dataset

    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )

    call cam_pio_createfile(File, fname_srf_cam, 0)
    call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), dof, iodesc)

    nf_x2a = mct_aVect_nRattr(x2a_a)
    allocate(varid_x2a(nf_x2a))

    rcode = pio_def_dim(File,'x2a_nx',ngcols,dimid(1))
    do k = 1,nf_x2a
       call mct_aVect_getRList(mstring,k,x2a_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2a_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2a(k))
       rcode = pio_put_att(File,varid_x2a(k),"_fillvalue",fillvalue)
    enddo

    nf_a2x = mct_aVect_nRattr(a2x_a)
    allocate(varid_a2x(nf_a2x))

    rcode = pio_def_dim(File,'a2x_nx',ngcols,dimid(1))
    do k = 1,nf_a2x
       call mct_aVect_getRList(mstring,k,a2x_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'a2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_a2x(k))
       rcode = PIO_put_att(File,varid_a2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already


    do k=1,nf_x2a
       call pio_write_darray(File, varid_x2a(k), iodesc, x2a_a%rattr(k,:), rcode)
    end do

    do k=1,nf_a2x
       call pio_write_darray(File, varid_a2x(k), iodesc, a2x_a%rattr(k,:), rcode)
    end do

    deallocate(varid_x2a, varid_a2x)

    call pio_freedecomp(File,iodesc)
    call cam_pio_closefile(file)


  end subroutine atm_write_srfrest_mct
!  This subroutine is made by Xiuhong Chen, June26, 2019

 subroutine read_surface_emis(ncols,lats2,lons2,mn,surface_emis)
     use netcdf

     use time_manager, only: get_curr_date
     use ppgrid

     use error_messages, only : handle_ncerr
     use parrrtm,        only: nbndlw
     use cam_logfile,     only: iulog
      
      implicit none
      integer :: ncid, status, latid,lonid,bandid,timeid
      character(256) filename

      integer :: ntime, nlat, nlon, nband,i,mn
      real, allocatable :: band_emissivity(:)
      character(len = nf90_max_name) :: RecordDimName
      integer :: lat_varID,emis_varID, lon_varID
      real, allocatable:: lat(:), lon(:)
      integer ::start(4),count(4)
      integer :: ncols,j
      integer ::ilats, ilons
      real(r8) :: lats2(ncols)           ! array of chunk latitudes
      real(r8) :: lons2(ncols)           ! array of chunk longitude
      real :: surface_emis(ncols, nbndlw), minvalue
!      filename = $inputdata/emis/surface_emissivity_1.9x2.5_RRTMG.nc"
      filename = "surface_emissivity_RRTMG_53deg.nc"
      status = nf90_open(trim(filename), nf90_nowrite, ncid)
      status = nf90_inq_dimid(ncid, "time", timeID)
      status = nf90_inq_dimid(ncid, "lat", latID)
      status = nf90_inq_dimid(ncid, "lon", lonID)
      status = nf90_inq_dimid(ncid, "band", bandID)



      status =  nf90_inquire_dimension( ncid, timeID,len=ntime )
      status =  nf90_inquire_dimension( ncid, latID,len=nlat )
      status =  nf90_inquire_dimension( ncid, lonID,len=nlon )
      status =  nf90_inquire_dimension( ncid, bandID,len=nband )


      allocate(band_emissivity(nband))
      allocate(lat(nlat))
      allocate(lon(nlon))

      status = nf90_inq_varid (ncid, 'lat', lat_varID )
      status = nf90_get_var (ncid, lat_varID, lat)

      status = nf90_inq_varid (ncid, 'lon', lon_varID )
      status = nf90_get_var (ncid, lon_varID, lon)

           count =(/nbndlw,1,1,1/)
      do i = 1, ncols
         minvalue= 10000.0
         do j=1, nlat
           if (abs(lat(j) - lats2(i)) .le. minvalue)  then
             ilats = j
             minvalue = abs(lat(j) - lats2(i))
           endif
         enddo
         minvalue= 10000.0
         do j=1, nlon
           if (abs(lon(j) - lons2(i)) .le. minvalue)  then
             ilons = j
             minvalue = abs(lon(j) - lons2(i))
           endif
         enddo
       !  ilons = minloc(lon, MASK=(lon .ge. (lons2(i)-0.7031/2)))
         start =(/1,ilons,ilats,mn/)
         status = nf90_inq_varid (ncid, 'band_emissivity', emis_varID )
         status = nf90_get_var (ncid, emis_varID, band_emissivity,start =start,count = count)

       surface_emis(i,:) = band_emissivity

      enddo
       status = NF90_CLOSE( ncid )

 end subroutine read_surface_emis

  ! This subroutine is made by Xianglei Huang in MATLAB,and is translated to
! FORTRAN by Xiuhong Chen, on June26,2019
subroutine get_Ts_from_LW_emis(v1, emis, LW, LWdown, v1_num, nguass_point, Ts)
! Update on Oct 08, 2015 by Xianglei
! rewrite the iteration part to handle the case caused by numerical
! difference betwee sigma*T^4 and integration by numerics

! This is a subroutine to obtain the Ts that can give the right upward LW
! flux with given band-averaged surface emissivity
! Input variables:
! v1, band ranges.e.g. [v11, v12, v13] will indicates 2bands, with
! band1 from [v11, v12] and band2 from [v12, v13];
! emis, the band-average surface emissivity, its size should be len(v1)-1
! LW, the upward LW flux in Wm^{-2}

   use cam_logfile,     only: iulog

   integer   i, j, Count, v1_num, nguass_point
   real(8)   v1(v1_num)
   real      emis(v1_num-1)
   real(8)   LW, LW2, LW3
   real(8)    LWdown(v1_num-1)
   real      T1, T2, T3, F1, F2, F3, A1, A2, A3
   real(8)   x(nguass_point), w(nguass_point)
   real      rad1, rad2, rad3, pi
   real(8)   Ts

   pi = 4.*atan(1.)

! T1 and T2 should be two values encompassing all possible Ts
   T1 = 150;
   T2 = 400;

   F1 = 0;
   F2 = 0;
   F3 = 0;

   LW3 = LW
   do i = 1, v1_num - 1
        if (emis(i).eq.0) then
            emis(i) = 1.0
        endif
        LW3 = LW3 - (1 - emis(i)) * LWdown(i)
   enddo

   LW2 = int(LW3*10)/10.0_r8

   Count = 0

   do while (abs(F3-LW2)>0.001)

!       for count how many iteration needed

        Count = Count + 1
        T3 = (T1 + T2)/2

        do i = 1,v1_num - 1
              call gaulegf(v1(i), v1(i+1), nguass_point, x, w)

               A1 = 0
               A2 = 0
               A3 = 0

               do j = 1,nguass_point

                  A1 = A1 + planck(x(j), T1) * w(j)
                  A2 = A2 + planck(x(j), T2) * w(j)
                  A3 = A3 + planck(x(j), T3) * w(j)

                enddo

!              don't mix up i and j subscripts

               F1 = F1 + emis(i) * A1;
               F2 = F2 + emis(i) * A2;
               F3 = F3 + emis(i) * A3;
        enddo
! covert to Wm-2
        F1 = F1 * pi * 1e-3
        F2 = F2 * pi * 1e-3
        F3 = F3 * pi * 1e-3


        if (Count .eq.1 .and. (LW2 .lt. F1 .or. LW2 .gt. F2)) then
               write(iulog,*) 'CCC',LW2,LW,LWdown(1),emis(1),LWdown(6),emis(6),Count
               write(iulog, *)'the LW2 is too low or too high so it is beyond the reasonable range'
               return
        elseif (LW2 .lt. F1) then !the case due to numerical error
                T1 = T1 - 2
        elseif (LW2 .gt. F2) then ! the case due to numberical error
                T2 = T2 + 2
        elseif (LW2 .ge. F3) then
                T1 = T3
        elseif (LW2 .lt. F3) then
                T2 = T3
        endif
   enddo
!   Ts = T3
!   To correct the bias

    if (T3.lt.270) then
         Ts = T3 + (T3-160)*0.0008 + 0.13
    elseif (T3.lt.300) then
         Ts = T3 + (T3-270)*0.0006 + 0.2171
    elseif (T3.lt.320) then
         Ts = T3 + (T3-300)*0.0002 + 0.2342
    elseif (T3.lt.340) then
         Ts = T3 - (T3-320)*0.0004 + 0.2386
    else
         Ts = T3 - (T3-340)*0.0013+0.2316
    endif

end subroutine get_Ts_from_LW_emis
! gauleg.f90     P145 Numerical Recipes in Fortran
! compute x(i) and w(i)  i=1,n  Legendre ordinates and weights
! on interval -1.0 to 1.0 (length is 2.0)
! use ordinates and weights for Gauss Legendre integration
!
subroutine gaulegf(x1, x2, n, x, w)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: x1, x2
  double precision, dimension(n), intent(out) :: x, w
  integer :: i, j, m
  double precision :: p1, p2, p3, pp, xl, xm, z, z1
  double precision, parameter :: eps=3.d-14

  m = (n+1)/2
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)

  do i=1,m
    z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
    z1 = 0.0
    do while(abs(z-z1) .gt. eps)
      p1 = 1.0d0
      p2 = 0.0d0
      do j=1,n
        p3 = p2
        p2 = p1
        p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
      end do
      pp = n*(z*p1-p2)/(z*z-1.0d0)
      z1 = z
      z = z1 - p1/pp
    end do
    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
     w(n+1-i) = w(i)
  end do
end subroutine gaulegf
function  planck(freq,temp)
! freq in wavenumber, temp in Kelvin degree and
! radiance in 1e-3 W per square meter per sr per wavenumber

      real(8):: freq
      real ::  temp, ca, cb, cof, arg, zeroind
      real ::  planck

!      ca    = 3.741832e-05 / 3.14159
      ca = 1.191043934e-05

!      cb    = 1.438837896
      cb = 1.438769911

      cof   = ca * (freq **3)

      arg   = cb * freq

      planck    = cof / ( exp( (arg/ temp)) - 1.0 )



end function planck



end module atm_comp_mct
