! ================================================================================================================================
!  MODULE       : xios_orchidee
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   This module contains the initialization and interface to the XIOS code.
!!
!!\n DESCRIPTION: This module contains the interface for the use of the XIOS code. All call to XIOS are done in this module.
!!                Revision 965 of XIOS/trunk or later is needed. This version is also called XIOS2. 
!!                Older revisions and XIOS1 can not be used.
!!                
!!                Summury of subroutines
!!                      xios_orchidee_comm_init       : First call to XIOS to get the MPI communicator 
!!                      xios_orchidee_init            : Initialize variables needed for use of XIOS 
!!                                                      Deactivation of fields not calculated due specific run options
!!                      xios_orchidee_update_calendar : Update the calandar in XIOS
!!                      xios_orchidee_finalize        : Last call to XIOS for finalization
!!                      xios_orchidee_send_field      : Interface to send fields with 1, 2 or 3 dimensions to XIOS
!!                      xios_orchidee_send_field_r1d  : Internal subroutine for 1D(array) fields
!!                      xios_orchidee_send_field_r2d  : Internal subroutine for 2D fields
!!                      xios_orchidee_send_field_r3d  : Internal subroutine for 3D fields
!!
!!                It is only possible to use XIOS2. Note that compilation must be done with the preprocessing key XIOS 
!!                and CPP_PARA. Compiling without these keys makes it impossible to activate XIOS. 
!!                To activate running using XIOS, the flag XIOS_ORCHIDEE_OK=y must be set in run.def and the file iodef.xml must exist.  
!!
!! RECENT CHANGE(S): Created by Arnaud Caubel(LSCE), Josefine Ghattas (IPSL) 2013
!!                   Removed possibility to use XIOS1, 21/10/2016
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE_2_2/ORCHIDEE/src_parallel/xios_orchidee.f90 $
!! $Date: 2023-05-22 14:15:35 +0200 (Mon, 22 May 2023) $
!! $Revision: 8011 $
!! \n
!_ ================================================================================================================================

MODULE xios_orchidee

#ifdef XIOS
  USE xios
#endif
  USE defprec
  USE pft_parameters_var, ONLY : nvm
  USE constantes_var
  USE constantes_soil_var, ONLY : nstm, nscm, diaglev, check_cwrr, ok_freeze_cwrr
  USE time, ONLY : dt_sechiba
  USE vertical_soil_var, ONLY : ngrnd, nslm
  USE IOIPSL, ONLY : ioget_calendar, ju2ymds
  USE mod_orchidee_para_var
  USE mod_orchidee_transfert_para
  USE ioipsl_para

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xios_orchidee_init, xios_orchidee_change_context, &
            xios_orchidee_update_calendar, xios_orchidee_context_finalize, xios_orchidee_finalize, &
            xios_orchidee_close_definition, &
            xios_orchidee_send_field, xios_orchidee_recv_field, &
            xios_orchidee_set_file_attr, xios_orchidee_set_field_attr, xios_orchidee_set_fieldgroup_attr, xios_orchidee_setvar, xios_orchidee_addaxis


  !
  !! Declaration of public variables
  !
  LOGICAL, PUBLIC, SAVE           :: xios_orchidee_ok=.TRUE.     !! Use XIOS for diagnostic files
  !$OMP THREADPRIVATE(xios_orchidee_ok)
  LOGICAL, PUBLIC, SAVE           :: xios_interpolation          !! Do reading and interpolations with XIOS. If false, reading will be done with IOIOSL and interpolation using aggregate_p
  !$OMP THREADPRIVATE(xios_interpolation)

  REAL(r_std), PUBLIC, SAVE       :: xios_default_val=0          !! Default value (missing value) used in XIOS. The value 0 will be overwritten with the value taken from XIOS.
  !$OMP THREADPRIVATE(xios_default_val)

  !
  !! Declaration of internal variables
  !
#ifdef XIOS
  TYPE(xios_context)              :: ctx_hdl_orchidee      !! Handel for ORCHIDEE
  !$OMP THREADPRIVATE(ctx_hdl_orchidee)
#endif



  !! ==============================================================================================================================
  !! INTERFACE   : xios_orchidee_send_field
  !!
  !>\BRIEF         Send a field to XIOS.
  !!
  !! DESCRIPTION  :\n Send a field to XIOS. The field can have 1, 2 or 3 dimensions.
  !!                  This interface should be called at each time-step for each output varaiables.
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE xios_orchidee_send_field
     MODULE PROCEDURE xios_orchidee_send_field_r1d, xios_orchidee_send_field_r2d, xios_orchidee_send_field_r3d, &
                      xios_orchidee_send_field_r4d, xios_orchidee_send_field_r5d
  END INTERFACE

  INTERFACE xios_orchidee_recv_field
     MODULE PROCEDURE xios_orchidee_recv_field_r1d, xios_orchidee_recv_field_r2d, xios_orchidee_recv_field_r3d
  END INTERFACE


CONTAINS


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_init 
  !!
  !>\BRIEF         Initialize variables needed for use of XIOS.
  !!
  !! DESCRIPTION  :\n Initialization of specific varaiables needed to use XIOS such as model domain and time step. 
  !!
  !!                  In this subroutine also a section containg deactivation of some fields is found. The variables are 
  !!                  deactivated of not according to the corresponding control flag. For exemple the variables cacluated by the 
  !!                  routing scheme will be deactivated if the routing is deactivated. This is done to be able to keep the same 
  !!                  iodef.xml input file for several options without geting empty fields in the output file. Note that a field that
  !!                  is activated in the code can always be deactivated from the iodef.xml external file. 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_init(MPI_COMM_ORCH,                   &
       date0,    year,      month,             day, julian_diff, &
       lon_mpi,  lat_mpi )

    USE grid, ONLY : grid_type, unstructured, regular_lonlat, regular_xy, nvertex, &
                     longitude, latitude, bounds_lon, bounds_lat, ind_cell_glo

    USE vertical_soil_var, ONLY : znt, znh

    IMPLICIT NONE
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                            :: MPI_COMM_ORCH    !! Orchidee MPI communicator (from module mod_orchidee_mpi_data)
    REAL(r_std), INTENT(in)                               :: date0            !! Julian day at first time step
    INTEGER(i_std), INTENT(in)                            :: year, month, day !! Current date information
    REAL(r_std), INTENT(in)                               :: julian_diff      !! Current day in the year [1,365(366)]
    REAL(r_std),DIMENSION (iim_g,jj_nb), INTENT(in)       :: lon_mpi, lat_mpi !! Longitudes and latitudes on MPI local domain 2D domain
    !
    !! 0.2 Local variables
    !
#ifdef XIOS

    TYPE(xios_duration)            :: dtime_xios
    TYPE(xios_date)                :: start_date
    TYPE(xios_date)                :: time_origin
    TYPE(xios_fieldgroup)          :: fieldgroup_handle
    TYPE(xios_field)               :: field_handle
    TYPE(xios_file)                :: file_handle
#endif
    INTEGER(i_std)                 :: i
    INTEGER(i_std)                 :: year0, month0, day0 !! Time origin date information
    REAL(r_std)                    :: sec0                !! Time origin date information
    CHARACTER(LEN=20)              :: calendar_str        !! Name of current calendar
    CHARACTER(LEN=30)              :: start_str           !! Current date as character string
    CHARACTER(LEN=30)              :: startorig_str       !! Time origin date as character string

    REAL(r_std),ALLOCATABLE        :: longitude_mpi(:), latitude_mpi(:)
    REAL(r_std),ALLOCATABLE        :: bounds_lon_mpi(:,:),bounds_lat_mpi(:,:) 
    INTEGER(i_std),ALLOCATABLE     :: ind_cell_mpi(:) 
    LOGICAL                        :: xios_remap_output
    !_ ================================================================================================================================
    
    
    IF (printlev>=3) WRITE(numout,*) 'Entering xios_orchidee_init'

    !Config Key   = XIOS_ORCHIDEE_OK
    !Config Desc  = Use XIOS for writing diagnostics file
    !Config If    = 
    !Config Def   = y 
    !Config Help  = Compiling and linking with XIOS library is necessary. 
    !Config Units = [FLAG]
    CALL getin_p('XIOS_ORCHIDEE_OK',xios_orchidee_ok)
    IF (printlev>=1) WRITE(numout,*)'In xios_orchidee_init, xios_orchidee_ok=',xios_orchidee_ok

    
    ! Coherence test between flag and preprocessing key
#ifndef XIOS
    IF (xios_orchidee_ok) THEN
       CALL ipslerr_p(3,'xios_orchidee_init', 'Preprocessing key XIOS is missing to run ORCHIDEE with XIOS',&
            'Recompile with preprocessing flag XIOS or set XIOS_ORCHIDEE_OK=n in run.def', '')
    END IF
#endif



    IF (xios_orchidee_ok) THEN
      !Config Key   = XIOS_INTERPOLATION
      !Config Desc  = Actiave reading and intrepolation using XIOS
      !Config If    = XIOS_ORCHIDEE_OK
      !Config Def   = n
      !Config Help  = This flag allows the user to decide to use xios
      !Config         interpolation or standard method for reading input files
      !Config Units = [FLAG]
      xios_interpolation = .FALSE.
      CALL getin_p('XIOS_INTERPOLATION', xios_interpolation)


      !Config Key   = XIOS_REMAP_OUTPUT
      !Config Desc  = Actiave remaping of diagnostic output files to regular grid
      !Config If    = XIOS_ORCHIDEE_OK .AND. grid_type=unstructured
      !Config Def   = True
      !Config Help  = Set this flag to false to output an unstructured grid on its natvie grid without interpolation
      !Config Units = [FLAG]
      xios_remap_output=.TRUE.
      CALL getin_p("XIOS_REMAP_OUTPUT",xios_remap_output)  

   ELSE
      ! Deactivate interpolation with XIOS not possible wihtout having
      ! xios_orchidee_ok=true
      xios_interpolation = .FALSE.
   END IF

   ! Force xios_interpolation=.TRUE. if using unstructured grid
   IF (grid_type==unstructured .AND. .NOT. xios_interpolation) THEN
      WRITE(numout,*) 'xios_interpolation must be true for unstructured grid. It is now changed to true.'
      xios_interpolation=.TRUE.
   END IF
   IF (printlev>=1) WRITE(numout,*)'In xios_orchidee_init, xios_interpolation=', xios_interpolation


    !
    !! 1. Set date and calendar information on the format needed by XIOS
    !

    ! Get the calendar from IOIPSL and modify the string to correspond to what XIOS expects
    CALL ioget_calendar(calendar_str)

    IF (calendar_str == 'gregorian') THEN
       calendar_str='gregorian'
    ELSE IF (calendar_str == 'noleap') THEN
       calendar_str='noleap'
    ELSE IF (calendar_str == '360d' .OR. calendar_str == '360_day') THEN
       calendar_str='d360'
    END IF

    ! Transform the time origin from julian days into year, month, day and seconds
    CALL ju2ymds(date0, year0, month0, day0, sec0)

    IF (grid_type==unstructured) THEN
      IF (is_omp_root) THEN
        ALLOCATE(longitude_mpi(ij_nb))
        ALLOCATE(latitude_mpi(ij_nb))
        ALLOCATE(bounds_lon_mpi(ij_nb,nvertex))
        ALLOCATE(bounds_lat_mpi(ij_nb,nvertex))
        ALLOCATE(ind_cell_mpi(ij_nb))
      ELSE
        ALLOCATE(longitude_mpi(0))
        ALLOCATE(latitude_mpi(0))
        ALLOCATE(bounds_lon_mpi(0,0))
        ALLOCATE(bounds_lat_mpi(0,0))
        ALLOCATE(ind_cell_mpi(0))
      ENDIF
      
      CALL gather_unindexed_omp(longitude,longitude_mpi)
      CALL gather_unindexed_omp(latitude,latitude_mpi)
      CALL gather_unindexed_omp(bounds_lon,bounds_lon_mpi)
      CALL gather_unindexed_omp(bounds_lat,bounds_lat_mpi)
      CALL gather_unindexed_omp(ind_cell_glo,ind_cell_mpi)
    ENDIF
    
    
    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       !
       !! 2. Context initialization
       !
       CALL xios_context_initialize("orchidee",MPI_COMM_ORCH)
       CALL xios_get_handle("orchidee",ctx_hdl_orchidee)
       CALL xios_set_current_context(ctx_hdl_orchidee)

       !
       !! 2. Calendar, timstep and date definition
       !
       dtime_xios%second=dt_sechiba

       CALL xios_define_calendar(type=calendar_str, start_date=xios_date(year,month,day,0,0,0), &
            time_origin=xios_date(year0,month0,day0,0,0,0), timestep=dtime_xios)

       !
       !! 3. Domain definition
       !
       IF (grid_type==regular_lonlat) THEN
          ! Global domain
          CALL xios_set_domain_attr("domain_landpoints", ni_glo=iim_g, nj_glo=jjm_g)
          ! Local MPI domain
          CALL xios_set_domain_attr("domain_landpoints",type="rectilinear", ibegin=0, ni=iim_g, jbegin=jj_begin-1, nj=jj_nb)
          
          ! Define how data is stored on memory : 1D array for only continental points
          CALL xios_set_domain_attr("domain_landpoints",data_dim=1, data_ibegin=0, data_ni=nbp_mpi)
          CALL xios_set_domain_attr("domain_landpoints",data_ni=nbp_mpi, data_i_index=kindex_mpi-1)     
          
          ! Define longitudes and latitudes on local MPI domain
          CALL xios_set_domain_attr("domain_landpoints",lonvalue_1d=lon_mpi(:,1),latvalue_1d=lat_mpi(1,:))
          
       ELSE IF (grid_type==regular_xy ) THEN
          ! Global domain
          CALL xios_set_domain_attr("domain_landpoints", ni_glo=iim_g, nj_glo=jjm_g)
          ! Local MPI domain
          CALL xios_set_domain_attr("domain_landpoints",type="curvilinear", ibegin=0, ni=iim_g, jbegin=jj_begin-1, nj=jj_nb)

          ! Define how data is stored on memory : 1D array for only continental points
          CALL xios_set_domain_attr("domain_landpoints",data_dim=1, data_ibegin=0, data_ni=nbp_mpi)
          CALL xios_set_domain_attr("domain_landpoints",data_ni=nbp_mpi, data_i_index=kindex_mpi-1)     

          ! Define longitudes and latitudes on local MPI domain depending on grid_type
          CALL xios_set_domain_attr("domain_landpoints",lonvalue_2d=lon_mpi,latvalue_2d=lat_mpi)

       ELSE IF (grid_type==unstructured) THEN
          
          ! Global domain
          CALL xios_set_domain_attr("domain_landpoints", ni_glo=jjm_g, type="unstructured", nvertex=nvertex)
          ! Local MPI domain
          CALL xios_set_domain_attr("domain_landpoints", ibegin=ij_begin-1, ni=ij_nb)
          
          ! Define how data is stored on memory : 1D array for only continental points
          CALL xios_set_domain_attr("domain_landpoints",data_dim=1, data_ni=nbp_mpi, data_i_index=kindex_mpi-1) 
          
          ! Define longitudes and latitudes on local MPI domain
          CALL xios_set_domain_attr("domain_landpoints",lonvalue_1d=longitude_mpi,latvalue_1d=latitude_mpi)
          CALL xios_set_domain_attr("domain_landpoints",bounds_lon_1d=RESHAPE(bounds_lon_mpi,(/nvertex,ij_nb/),order=(/2,1/)))
          CALL xios_set_domain_attr("domain_landpoints",bounds_lat_1d=RESHAPE(bounds_lat_mpi,(/nvertex,ij_nb/),order=(/2,1/)))


          IF (xios_remap_output) THEN
             
             ! Define output grid as domain_landpoints_regular (grid specified in xml files)
             CALL xios_set_domain_attr("domain_landpoints_out",domain_ref="domain_landpoints_regular")
             
             CALL xios_set_fieldgroup_attr("remap_expr",expr="@this_ref")
             CALL xios_set_fieldgroup_attr("remap_1ts",   freq_op=xios_duration_convert_from_string("1ts"))
             CALL xios_set_fieldgroup_attr("remap_1800s", freq_op=xios_duration_convert_from_string("1800s"))
             CALL xios_set_fieldgroup_attr("remap_1h",    freq_op=xios_duration_convert_from_string("1h"))
             CALL xios_set_fieldgroup_attr("remap_3h",    freq_op=xios_duration_convert_from_string("3h"))
             CALL xios_set_fieldgroup_attr("remap_6h",    freq_op=xios_duration_convert_from_string("6h"))
             CALL xios_set_fieldgroup_attr("remap_1d",    freq_op=xios_duration_convert_from_string("1d"))
             CALL xios_set_fieldgroup_attr("remap_1mo",   freq_op=xios_duration_convert_from_string("1mo"))
             CALL xios_set_fieldgroup_attr("remap_1y",    freq_op=xios_duration_convert_from_string("1y"))
          ENDIF

       END IF

       !
       !! 4. Axis definition
       !
       CALL xios_set_axis_attr("nvm",n_glo=nvm ,VALUE=(/(REAL(i,r_std),i=1,nvm)/))
       CALL xios_set_axis_attr("nlut",n_glo=nlut ,VALUE=(/(REAL(i,r_std),i=1,nlut)/))
       CALL xios_set_axis_attr("ncarb",n_glo=ncarb ,VALUE=(/(REAL(i,r_std),i=1,ncarb)/))
       CALL xios_set_axis_attr("nparts",n_glo=nparts,VALUE=(/(REAL(i,r_std),i=1,nparts)/))
       CALL xios_set_axis_attr("nlaip1", n_glo=nlai+1,VALUE=(/(REAL(i,r_std),i=1,nlai+1)/))
       CALL xios_set_axis_attr("ngrnd",n_glo=ngrnd ,VALUE=znt(:))
       CALL xios_set_axis_attr("nstm", n_glo=nstm,VALUE=(/(REAL(i,r_std),i=1,nstm)/))
       CALL xios_set_axis_attr("nscm", n_glo=nscm,VALUE=(/(REAL(i,r_std),i=1,nscm)/))
       CALL xios_set_axis_attr("nnobio", n_glo=nnobio,VALUE=(/(REAL(i,r_std),i=1,nnobio)/))
       CALL xios_set_axis_attr("albtyp", n_glo=2,VALUE=(/(REAL(i,r_std),i=1,2)/))
       CALL xios_set_axis_attr("nslm", n_glo=nslm,VALUE=znh(:))
       CALL xios_set_axis_attr("P10", n_glo=10,VALUE=(/(REAL(i,r_std), i=1,10)/))
       CALL xios_set_axis_attr("P100", n_glo=100,VALUE=(/(REAL(i,r_std), i=1,100)/))
       CALL xios_set_axis_attr("P11", n_glo=11,VALUE=(/(REAL(i,r_std), i=1,11)/))
       CALL xios_set_axis_attr("P101", n_glo=101,VALUE=(/(REAL(i,r_std), i=1,101)/))
       CALL xios_set_axis_attr("nsnow", n_glo=nsnow,VALUE=(/(REAL(i,r_std),i=1,nsnow)/))
              
       !
       !! 5. Get the default value (missing value) used by XIOS. This value is set in field_def_orchidee.xml
       !
       CALL xios_get_fieldgroup_attr("field_definition", default_value=xios_default_val)
       IF (printlev>=2) WRITE(numout,*) 'Default value read from XIOS, xios_default_val=',xios_default_val

       !
       !! 5. Deactivation of some fields if they are not calculated
       !
       IF ( OFF_LINE_MODE ) THEN
          CALL xios_set_field_attr("riverflow_cpl",enabled=.FALSE.)
          CALL xios_set_field_attr("coastalflow_cpl",enabled=.FALSE.)
       END IF

       IF ( .NOT. river_routing ) THEN
          CALL xios_set_field_attr("basinmap",enabled=.FALSE.)
          CALL xios_set_field_attr("nbrivers",enabled=.FALSE.)
          CALL xios_set_field_attr("riversret",enabled=.FALSE.)
          CALL xios_set_field_attr("hydrographs",enabled=.FALSE.)
          CALL xios_set_field_attr("htuhgmon",enabled=.FALSE.)
          CALL xios_set_field_attr("fastr",enabled=.FALSE.)
          CALL xios_set_field_attr("slowr",enabled=.FALSE.)
          CALL xios_set_field_attr("streamr",enabled=.FALSE.)
          CALL xios_set_field_attr("laker",enabled=.FALSE.)
          CALL xios_set_field_attr("lake_overflow",enabled=.FALSE.)
          CALL xios_set_field_attr("mask_coast",enabled=.FALSE.)
          CALL xios_set_field_attr("pondr",enabled=.FALSE.)
          CALL xios_set_field_attr("floodr",enabled=.FALSE.)
          CALL xios_set_field_attr("slowflow",enabled=.FALSE.)
          CALL xios_set_field_attr("delfastr",enabled=.FALSE.)
          CALL xios_set_field_attr("delslowr",enabled=.FALSE.)
          CALL xios_set_field_attr("delstreamr",enabled=.FALSE.)
          CALL xios_set_field_attr("dellaker",enabled=.FALSE.)
          CALL xios_set_field_attr("delpondr",enabled=.FALSE.)
          CALL xios_set_field_attr("delfloodr",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigmap",enabled=.FALSE.)
          CALL xios_set_field_attr("swampmap",enabled=.FALSE.)
          CALL xios_set_field_attr("wbr_stream",enabled=.FALSE.)
          CALL xios_set_field_attr("wbr_fast",enabled=.FALSE.)
          CALL xios_set_field_attr("wbr_slow",enabled=.FALSE.)
          CALL xios_set_field_attr("wbr_lake",enabled=.FALSE.)
          CALL xios_set_field_attr("reinfiltration",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigation",enabled=.FALSE.)
          CALL xios_set_field_attr("netirrig",enabled=.FALSE.)
          CALL xios_set_field_attr("SurfStor",enabled=.FALSE.)
          CALL xios_set_field_attr("htutempmon",enabled=.FALSE.)
          CALL xios_set_field_attr("streamlimit",enabled=.FALSE.)
          CALL xios_set_field_attr("StreamT_TotTend",enabled=.FALSE.)
          CALL xios_set_field_attr("StreamT_AdvTend",enabled=.FALSE.)
          CALL xios_set_field_attr("StreamT_RelTend",enabled=.FALSE.)
       END IF

       IF (.NOT. ok_freeze_cwrr) THEN
          CALL xios_set_field_attr("profil_froz_hydro",enabled=.FALSE.)
       END IF

       
       IF (.NOT. check_cwrr) THEN
          CALL xios_set_field_attr("check_infilt",enabled=.FALSE.)
          CALL xios_set_field_attr("check_tr",enabled=.FALSE.)
          CALL xios_set_field_attr("check_over",enabled=.FALSE.)
          CALL xios_set_field_attr("check_under",enabled=.FALSE.)
          CALL xios_set_field_attr("check_top",enabled=.FALSE.)
          CALL xios_set_field_attr("qflux",enabled=.FALSE.)
       END IF

       IF ( .NOT. do_floodplains ) THEN
          CALL xios_set_field_attr("floodmap",enabled=.FALSE.)
          CALL xios_set_field_attr("floodh",enabled=.FALSE.)       
          CALL xios_set_field_attr("floodr",enabled=.FALSE.)
          CALL xios_set_field_attr("floodout",enabled=.FALSE.)
          CALL xios_set_field_attr("flood_frac",enabled=.FALSE.)       
       END IF

       ! Deactivate some stomate fields. 
       ! These fields were traditionally added in sechiba_history.nc output file.
       IF ( .NOT. ok_stomate ) THEN
          CALL xios_set_field_attr("nee",enabled=.FALSE.)
          CALL xios_set_field_attr("maint_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("hetero_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("growth_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("npp",enabled=.FALSE.)
       END IF

       IF ( .NOT. do_irrigation ) THEN
          CALL xios_set_field_attr("irrigation",enabled=.FALSE.)
          CALL xios_set_field_attr("netirrig",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigmap",enabled=.FALSE.)
          CALL xios_set_field_attr("irrig_deficit",enabled=.FALSE.)
          CALL xios_set_field_attr("irrig_adduct",enabled=.FALSE.)
          CALL xios_set_field_attr("irrig_gw_source",enabled=.FALSE.)
          CALL xios_set_field_attr("irrig_fast_source",enabled=.FALSE.)
          CALL xios_set_field_attr("irrig_str_source",enabled=.FALSE.)
          CALL xios_set_field_attr("Count_failure_slow",enabled=.FALSE.)
          CALL xios_set_field_attr("Count_failure_fast",enabled=.FALSE.)
          CALL xios_set_field_attr("Count_failure_stre",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigmap_dyn",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc)THEN
          CALL xios_set_field_attr("PAR",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_fertil_no",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_iso",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_mono",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_ORVOC",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_MBO",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_methanol",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetone",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetal",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_formal",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetic",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_formic",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_no_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_no",enabled=.FALSE.)
          CALL xios_set_field_attr('flx_apinen'   ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_bpinen'   ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_limonen'  ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_myrcen'   ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_sabinen'  ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_camphen'  ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_3caren'   ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_tbocimen' ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_othermono',enabled=.FALSE.)
          CALL xios_set_field_attr('flx_sesquiter',enabled=.FALSE.)
          CALL xios_set_field_attr("CRF",enabled=.FALSE.)
          CALL xios_set_field_attr("fco2",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_radcanopy ) THEN
          CALL xios_set_field_attr("PARdf",enabled=.FALSE.)
          CALL xios_set_field_attr("PARdr",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_radcanopy .OR. .NOT. ok_multilayer ) THEN
          CALL xios_set_field_attr( 'PARsuntab',enabled=.FALSE.)
          CALL xios_set_field_attr( 'PARshtab' ,enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_radcanopy .OR. ok_multilayer ) THEN
          CALL xios_set_field_attr("PARsun",enabled=.FALSE.)
          CALL xios_set_field_attr("PARsh",enabled=.FALSE.)
          CALL xios_set_field_attr("laisun",enabled=.FALSE.)
          CALL xios_set_field_attr("laish",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_bbgfertil_Nox) THEN
          CALL xios_set_field_attr("flx_co2_bbg_year",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_cropsfertil_Nox) THEN
          CALL xios_set_field_attr("N_qt_WRICE_year",enabled=.FALSE.)
          CALL xios_set_field_attr("N_qt_OTHER_year",enabled=.FALSE.)
       END IF

       ! Set record_offset for enable start in the middle of the year.
       ! julian_diff is the day of the year where the current run start
       IF (printlev>=3) WRITE(numout,*) 'In xios_orchidee_init, julian_diff, INT(julian_diff) =', &
            julian_diff, INT(julian_diff)

       IF (ok_nudge_mc .AND. nudge_interpol_with_xios) THEN
          ! Activate the input file with id="nudge_moistc" specified in file_def_orchidee.xml. 
          ! The nudging file should be called nudge_moistc.nc (see name in the xml file) and is 
          ! supposed to contain daily values for the full year for the variable moistc.
          CALL xios_set_file_attr("nudge_moistc",enabled=.TRUE.)
          ! Set record_offset to start read at correct day in the nudging file. 
          CALL xios_set_file_attr("nudge_moistc",record_offset=INT(julian_diff))
       ELSE
          ! Deactivate input file for nudging of soil moisture
          CALL xios_set_file_attr("nudge_moistc",enabled=.FALSE.)
          ! Deactivate variables related to soil moisture nudgnig
          CALL xios_set_field_attr("mask_moistc_interp",enabled=.FALSE.)
          CALL xios_set_field_attr("moistc_interp",enabled=.FALSE.)

          ! Deactivate output variables related to soil moisture nudging
          CALL xios_set_field_attr("mc_read_current",enabled=.FALSE.)
          CALL xios_set_field_attr("mc_read_prev",enabled=.FALSE.)
          CALL xios_set_field_attr("mc_read_next",enabled=.FALSE.)
          CALL xios_set_field_attr("mask_mc_interp_out",enabled=.FALSE.)
       END IF
       IF (.NOT. ok_nudge_mc ) CALL xios_set_field_attr("nudgincsm",enabled=.FALSE.)

       IF (ok_nudge_snow .AND. nudge_interpol_with_xios) THEN
          ! Activate the input file with id="nudge_snow" specified in file_def_orchidee.xml. 
          ! The nudging file should be called nudge_snow.nc (see name in the xml file) and is 
          ! supposed to contain daily values for the full year for the variables snowdz, snowtemp and snowrho.
          CALL xios_set_file_attr("nudge_snow",enabled=.TRUE.)
          ! Set record_offset to start read at correct day in the nudging file. 
          CALL xios_set_file_attr("nudge_snow",record_offset=INT(julian_diff))
       ELSE
          ! Deactivate input file for nudging of snow variables
          CALL xios_set_file_attr("nudge_snow",enabled=.FALSE.)

          ! Deactivate input variables related to snow nudging
          CALL xios_set_field_attr("mask_snow_interp",enabled=.FALSE.)
          CALL xios_set_field_attr("snowdz_interp",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho_interp",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_interp",enabled=.FALSE.)

          ! Deactivate output variables related to snow nudging
          CALL xios_set_field_attr("snowdz_read_current",enabled=.FALSE.)
          CALL xios_set_field_attr("snowdz_read_prev",enabled=.FALSE.)
          CALL xios_set_field_attr("snowdz_read_next",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho_read_current",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho_read_prev",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho_read_next",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_read_current",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_read_prev",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_read_next",enabled=.FALSE.)
          CALL xios_set_field_attr("mask_snow_interp_out",enabled=.FALSE.)
       END IF
       IF (.NOT. ok_nudge_snow) CALL xios_set_field_attr("nudgincswe",enabled=.FALSE.)

       IF (impaze) THEN
          CALL xios_set_field_attr("soilalb_vis",enabled=.FALSE.)
          CALL xios_set_field_attr("soilalb_nir",enabled=.FALSE.)
          CALL xios_set_field_attr("vegalb_vis",enabled=.FALSE.)
          CALL xios_set_field_attr("vegalb_nir",enabled=.FALSE.)
       END IF

       IF (.NOT. do_wood_harvest) THEN
          CALL xios_set_field_attr("PROD10_HARVEST",enabled=.FALSE.)
          CALL xios_set_field_attr("FLUX10_HARVEST",enabled=.FALSE.)
          CALL xios_set_field_attr("PROD100_HARVEST",enabled=.FALSE.)
          CALL xios_set_field_attr("FLUX100_HARVEST",enabled=.FALSE.)
          CALL xios_set_field_attr("CONVFLUX_HARVEST",enabled=.FALSE.)
          CALL xios_set_field_attr("CFLUX_PROD10_HARVEST",enabled=.FALSE.)
          CALL xios_set_field_attr("CFLUX_PROD100_HARVEST",enabled=.FALSE.)
          CALL xios_set_field_attr("WOOD_HARVEST",enabled=.FALSE.)
          CALL xios_set_field_attr("WOOD_HARVEST_PFT",enabled=.FALSE.)
       END IF


       IF (grid_type==unstructured) THEN
          CALL xios_set_field_attr('RESOLUTION_X',enabled=.FALSE.)
          CALL xios_set_field_attr('RESOLUTION_Y',enabled=.FALSE.)
       END IF
#endif
    END IF

    IF (xios_orchidee_ok) THEN
       ! Send variables to all OMP thredds
       CALL bcast(xios_default_val)
       CALL bcast(almaoutput)
    END IF

    IF (printlev>=3) WRITE(numout,*) 'End xios_orchidee_init'
  END SUBROUTINE xios_orchidee_init


  SUBROUTINE xios_orchidee_close_definition

    IF (printlev >=4) WRITE(numout,*) 'Start xios_orchidee_close_definition'
    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS

       !
       !! 6. Close context
       !
       CALL xios_close_context_definition()
       WRITE(numout,*) 'Done xios_orchidee_close_context'      

       !
       !! 7. Activate almaoutput if needed 
       !! Some extra calculations have to be done for the variables  
       !! delsoilmoist, delintercept, delswe and soilwet.
       !! Set almaoutput=true if at least one of these variables are defined in an output file. 
       !! If not, keep the initial value of almaoutput. 
       IF ( xios_field_is_active("delsoilmoist") .OR. xios_field_is_active("delintercept") .OR. &
            xios_field_is_active("delswe")       .OR. xios_field_is_active("soilwet")      .OR. &
            xios_field_is_active("twbr")) THEN

          almaoutput=.TRUE.
          IF (printlev >=3) WRITE(numout,*) 'The flag almaoutput has been activated in xios_orchidee_init'
       END IF
#endif
    END IF

    IF (xios_orchidee_ok) THEN
       ! Send variables to all OMP thredds
       CALL bcast(xios_default_val)
       CALL bcast(almaoutput)
    END IF
    IF (printlev >=4) WRITE(numout,*) 'End xios_orchidee_close_definition'
  END SUBROUTINE xios_orchidee_close_definition
  
  
  
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_change_context
  !!
  !>\BRIEF         Use this subroutine to switch between different context.
  !!               This subroutine must be called when running in coupled mode at each time ORCHIDEE is called, in the
  !!               begining and end of intersurf_gathered. First call is done after xios_orchidee_init is done. 
  !!
  !! DESCRIPTION  :\n 
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_change_context(new_context)
    !
    !! 0. Variable and parameter declaration
    !
    !!    Input variable
    CHARACTER(LEN=*),INTENT(IN)              :: new_context

    !! Local variables
#ifdef XIOS
    TYPE(xios_context) :: ctx_hdl
#endif
    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       CALL xios_get_handle(new_context,ctx_hdl)
       CALL xios_set_current_context(ctx_hdl)
#endif
    END IF
    
  END SUBROUTINE xios_orchidee_change_context

  !!
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_addaxis
  !!
  !>\BRIEF         Use this subroutine to add axes, needed for nbasmon and nbasmax in routing_highres.f90
  !!
  !!
  !! DESCRIPTION  :\n
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_addaxis(axname, axlen, axval)
    !
    ! INPUT variables
    CHARACTER(LEN=*), INTENT(IN)             :: axname
    INTEGER(i_std), INTENT(IN)                :: axlen
    REAL(r_std), DIMENSION(axlen), INTENT(IN) :: axval
    !
    IF (xios_orchidee_ok .AND. is_omp_root) THEN
       CALL xios_set_axis_attr(axname, n_glo=axlen, VALUE=axval)
    ENDIF
    !
  END SUBROUTINE xios_orchidee_addaxis
  !!


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_update_calendar
  !!
  !>\BRIEF          Update the calandar in XIOS.
  !!
  !! DESCRIPTION  :\n Update the calendar in XIOS : let XIOS know that ORCHIDEE avanced one time-step.
  !!                  This subroutine should be called in the beginning of each time-step. The first 
  !!                  time-step in a new execution should always start at 1. Therefore, first calculate
  !!                  an offset that is substracted to the current time step in sechiba. 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_update_calendar(itau_sechiba)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(IN) :: itau_sechiba    !! Current time step of the model
    !
    !! 0.2 Local variables
    !
    LOGICAL, SAVE         :: first=.TRUE.         !! Flag for first entering in subroutine
    INTEGER(i_std), SAVE  :: offset               !! Offset to substract from itau_sechiba
    INTEGER(i_std)        :: itau_xios            !! Current time step for XIOS

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       ! Calculate the offset
       IF (first) THEN
          offset=itau_sechiba-1
          first=.FALSE.
       END IF

       ! Substract the offset to the current time step in sechiba
       itau_xios=itau_sechiba-offset

       ! Send the new time step to XIOS
       IF (printlev>=3) WRITE(numout,*) 'xios_orchidee_update_calendar: itau_sechiba, itau_xios=',itau_sechiba,itau_xios
       CALL xios_update_calendar(itau_xios)
#endif
    END IF
  END SUBROUTINE xios_orchidee_update_calendar
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_context_finalize
  !!
  !>\BRIEF         Finalize orchidee context.
  !!
  !! DESCRIPTION  :\n This subroutine finalizes the orchidee context without finalizing XIOS. In coupled mode, the atmospheric
  !!                  modele must finalize XIOS. This subroutine is called in the end of the execution of ORCHIDEE only in 
  !!                  coupeld mode.
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_context_finalize

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
       IF (printlev>=3) WRITE(numout,*) 'Entering xios_orchidee_context_finalize'
#ifdef XIOS
       CALL xios_context_finalize()
#endif
    END IF
  END SUBROUTINE xios_orchidee_context_finalize


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_finalize
  !!
  !>\BRIEF         Last call to XIOS for finalization.
  !!
  !! DESCRIPTION  :\n Last call to XIOS for finalization of the orchidee context and XIOS.
  !!                  This subroutine is called only when ORCHIDEE is run in offline mode. In coupled mode it is the atmospheric
  !!                  model that finalizes XIOS. In that case, the context orchidee must be finalized using the 
  !!                  subroutine xios_orchidee_context_finalize
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_finalize

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
       IF (printlev>=3) WRITE(numout,*) 'Entering xios_orchidee_finalize'
#ifdef XIOS
       CALL xios_context_finalize()
       CALL xios_finalize()
#endif
    END IF
  END SUBROUTINE xios_orchidee_finalize


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r1d
  !!
  !>\BRIEF          Subroutine for sending 1D (array) fields to XIOS.
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 1D fields (array).
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r1d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)          :: field_id
    REAL(r_std), DIMENSION(:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r1d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r1d


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r2d
  !!
  !>\BRIEF          Subroutine for sending 2D fields to XIOS.
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 2D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r2d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)            :: field_id
    REAL(r_std), DIMENSION(:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r2d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r2d


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r3d
  !!
  !>\BRIEF          Subroutine for sending 3D fields to XIOS. 
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 3D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r3d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r3d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r3d

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r4d
  !!
  !>\BRIEF          Subroutine for sending 4D fields to XIOS. 
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 4D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r4d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    INTEGER :: jv
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3),size(field,4)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r4d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r4d

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r5d
  !!
  !>\BRIEF          Subroutine for sending 5D fields to XIOS. 
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 5D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r5d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    INTEGER :: jv
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3),size(field,4),size(field,5)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r5d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r5d
 
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_recv_field_r2d
  !!
  !>\BRIEF          Subroutine for receiving 1D (kjpindex) fields to XIOS. 
  !!
  !! DESCRIPTION  :\n 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_recv_field_r1d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:), INTENT(OUT)    :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi)           :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_recv_field_r1d, field_id=',field_id

       ! All master threads receive the field from XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_recv_field(field_id,field_mpi)
          IF (printlev>=5) WRITE(numout,*) 'Recieve done with xios_orchidee_recv_field_r1d, field_id=',field_id
#endif
       END IF

       ! Scatter the mpi domains on local omp domains
       CALL scatter_omp(field_mpi, field)

    END IF
  END SUBROUTINE xios_orchidee_recv_field_r1d

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_recv_field_r2d
  !!
  !>\BRIEF          Subroutine for receiving 2D(kjpindex and 1 vertical axe) fields to XIOS. 
  !!
  !! DESCRIPTION  :\n 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_recv_field_r2d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:), INTENT(OUT)  :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_recv_field_r2d, field_id=',field_id

       ! All master threads recieve the field from XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_recv_field(field_id,field_mpi)
          IF (printlev>=5) WRITE(numout,*) 'Recieve done with xios_orchidee_recv_field_r2d, field_id=',field_id
#endif
       END IF

       ! Scatter the mpi domains on local omp domains
       CALL scatter_omp(field_mpi, field)

    END IF
  END SUBROUTINE xios_orchidee_recv_field_r2d

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_recv_field_r3d
  !!
  !>\BRIEF          Subroutine for receiving 3D(kjpindex and 2 vertical axes) fields to XIOS. 
  !!
  !! DESCRIPTION  :\n 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_recv_field_r3d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(OUT) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_recv_field_r3d, field_id=',field_id

       ! All master threads receive the field from XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_recv_field(field_id,field_mpi)
          IF (printlev>=5) WRITE(numout,*) 'Recieve done with xios_orchidee_recv_field_r3d, field_id=',field_id
#endif
       END IF

       ! Scatter the mpi domains on local omp domains
       CALL scatter_omp(field_mpi, field)

    END IF
  END SUBROUTINE xios_orchidee_recv_field_r3d



  SUBROUTINE xios_orchidee_set_file_attr(attr, name, enabled)
    CHARACTER(LEN=*), INTENT(IN)            :: attr     ! Name of the attribut
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: name     ! New name
    LOGICAL, INTENT(IN), OPTIONAL             :: enabled ! Flag

    IF (xios_orchidee_ok .AND. is_omp_root) THEN

#ifdef XIOS
       IF (PRESENT(name) .AND. PRESENT(enabled)) THEN
         CALL xios_set_file_attr(attr, name=name, enabled=enabled)
       ELSE IF (PRESENT(name)) THEN
         CALL xios_set_file_attr(attr, name=name)
       ELSE IF (PRESENT(enabled)) THEN
         CALL xios_set_file_attr(attr, enabled=enabled)
       ELSE
         CALL xios_set_file_attr(attr)
       END IF
#endif

    END IF

  END SUBROUTINE xios_orchidee_set_file_attr
  
  SUBROUTINE xios_orchidee_set_field_attr(attr,name, enabled)
    CHARACTER(LEN=*), INTENT(IN)            :: attr     ! Name of the attribut
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: name     ! New name
    LOGICAL, INTENT(IN), OPTIONAL             :: enabled ! Flag

    IF (xios_orchidee_ok .AND. is_omp_root) THEN

#ifdef XIOS
       IF (PRESENT(name) .AND. PRESENT(enabled)) THEN
         CALL xios_set_field_attr(attr, name=name, enabled=enabled)
       ELSE IF (PRESENT(name)) THEN
         CALL xios_set_field_attr(attr, name=name)
       ELSE IF (PRESENT(enabled)) THEN
         CALL xios_set_field_attr(attr, enabled=enabled)
       ELSE
         CALL xios_set_field_attr(attr)
       END IF
#endif

    END IF


  END SUBROUTINE xios_orchidee_set_field_attr
  
  SUBROUTINE xios_orchidee_set_fieldgroup_attr(attr,name, enabled)
    CHARACTER(LEN=*), INTENT(IN)            :: attr     ! Name of the attribut
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL  :: name     ! New name
    LOGICAL, INTENT(IN), OPTIONAL             :: enabled ! Flag

    IF (xios_orchidee_ok .AND. is_omp_root) THEN

#ifdef XIOS
       IF (PRESENT(name) .AND. PRESENT(enabled)) THEN
         CALL xios_set_fieldgroup_attr(attr, name=name, enabled=enabled)
       ELSE IF (PRESENT(name)) THEN
         CALL xios_set_fieldgroup_attr(attr, name=name)
       ELSE IF (PRESENT(enabled)) THEN
         CALL xios_set_fieldgroup_attr(attr, enabled=enabled)
       ELSE
         CALL xios_set_fieldgroup_attr(attr)
       END IF
#endif

    END IF


  END SUBROUTINE xios_orchidee_set_fieldgroup_attr
  
  FUNCTION xios_orchidee_setvar(varname,varvalue) RESULT (out)
    CHARACTER(LEN=*), INTENT(IN) :: varname  ! Name of the variable
    REAL, INTENT(IN)               :: varvalue ! Value of the variable
    LOGICAL :: out

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
      out=xios_setvar(varname, varvalue)
#endif
    END IF

  END FUNCTION xios_orchidee_setvar

END MODULE xios_orchidee

