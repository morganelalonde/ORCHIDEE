
! =================================================================================================================================
! MODULE       : slowproc
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         Groups the subroutines that: (1) initialize all variables used in 
!! slowproc_main, (2) prepare the restart file for the next simulation, (3) Update the 
!! vegetation cover if needed, and (4) handle all slow processes if the carbon
!! cycle is activated (call STOMATE) or update the vegetation properties (LAI and 
!! fractional cover) in the case of a run with only SECHIBA.
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): Allowed reading of USDA map, Nov 2014, ADucharne
!!                   November 2020: It is possible to define soil hydraulic parameters from maps,
!!                   as needed for the SP-MIP project (Tafasca Salma and Ducharne Agnes).
!!                   Changes in slowproc_xios_initialize, slowproc_soilt and slowproc_finalize
!!                   July 2022: New irrigation scheme. Here interpolation of new maps for
!!                   the irrigation scheme
!!
!! REFERENCE(S)	:
!!- Tafasca S. (2020). Evaluation de l’impact des propriétés du sol sur l’hydrologie simulee dans le
!! modèle ORCHIDEE, PhD thesis, Sorbonne Universite. \n
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE_2_2/ORCHIDEE/src_sechiba/slowproc.f90 $
!! $Date: 2022-07-20 11:30:43 +0200 (Wed, 20 Jul 2022) $
!! $Revision: 7709 $
!! \n
!_ ================================================================================================================================

MODULE slowproc

  USE defprec
  USE constantes 
  USE constantes_soil
  USE pft_parameters
  USE ioipsl
  USE xios_orchidee
  USE ioipsl_para
  USE sechiba_io_p
  USE interpol_help
  USE stomate
  USE stomate_data
  USE grid
  USE time, ONLY : dt_sechiba, dt_stomate, one_day, FirstTsYear, LastTsDay
  USE time, ONLY : year_start, month_start, day_start, sec_start
  USE time, ONLY : month_end, day_end
  USE mod_orchidee_para

  IMPLICIT NONE

  ! Private & public routines

  PRIVATE
  PUBLIC slowproc_main, slowproc_clear, slowproc_initialize, slowproc_finalize, slowproc_change_frac, slowproc_xios_initialize

  !
  ! variables used inside slowproc module : declaration and initialisation
  !
  REAL(r_std), SAVE                                  :: slope_default = 0.1
!$OMP THREADPRIVATE(slope_default)
  INTEGER, SAVE                                      :: printlev_loc        !! Local printlev in slowproc module
!$OMP THREADPRIVATE(printlev_loc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: clayfraction        !! Clayfraction (0-1, unitless)
!$OMP THREADPRIVATE(clayfraction)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: sandfraction        !! Sandfraction (0-1, unitless)
!$OMP THREADPRIVATE(sandfraction)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: siltfraction        !! Siltfraction (0-1, unitless)
!$OMP THREADPRIVATE(siltfraction)  
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:)  :: laimap              !! LAI map when the LAI is prescribed and not calculated by STOMATE
!$OMP THREADPRIVATE(laimap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: veget_max_new       !! New year fraction of vegetation type (0-1, unitless)
!$OMP THREADPRIVATE(veget_max_new)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: irrigated_new       !! New year area equipped for irrigation  (m^{-2})
!$OMP THREADPRIVATE(irrigated_new)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)      :: woodharvest         !! New year wood harvest
!$OMP THREADPRIVATE(woodharvest)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)    :: frac_nobio_new      !! New year fraction of ice+lakes+cities+... (0-1, unitless)
!$OMP THREADPRIVATE(frac_nobio_new)
  INTEGER(i_std), SAVE                               :: lcanop              !! canopy levels used for LAI
!$OMP THREADPRIVATE(lcanop)
  INTEGER(i_std) , SAVE                              :: veget_year          !! year for vegetation update
!$OMP THREADPRIVATE(veget_year)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)       :: frac_imperv         !! Imperviousness fraction of each gridcell
!$OMP THREADPRIVATE(imperviousness)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)       :: coeff_imperv      !! Imperviousness coefficient to modify Ks
!$OMP THREADPRIVATE(imperviousness)

CONTAINS




!!  =============================================================================================================================
!! SUBROUTINE:    slowproc_xios_initialize
!!
!>\BRIEF	  Initialize xios dependant defintion before closing context defintion
!!
!! DESCRIPTION:	  Initialize xios dependant defintion before closing context defintion
!!  
!! RECENT CHANGE(S): Initialization of XIOS to read soil hydraulic parameters from maps,
!!                   as needed for the SP-MIP project (Tafasca Salma and Ducharne Agnes).
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE slowproc_xios_initialize

    CHARACTER(LEN=255) :: filename, name
    LOGICAL :: lerr
    REAL(r_std) :: slope_noreinf
    LOGICAL :: get_slope
    CHARACTER(LEN=30) :: veget_str         !! update frequency for landuse    
    INTEGER :: l

    IF (printlev>=3) WRITE(numout,*) 'In slowproc_xios_initialize'
    !! 1. Prepare for reading of soils_param file
    ! Get the file name from run.def file and set file attributes accordingly
    filename = 'soils_param.nc'
    CALL getin_p('SOILCLASS_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("soils_param_file",name=name)

    ! Determine if soils_param_file will be read. If not, deactivate the file.    
    IF (xios_interpolation .AND. restname_in=='NONE' .AND. .NOT. impsoilt) THEN
       ! Reading will be done with XIOS later
       IF (printlev>=2) WRITE(numout,*) 'Reading of soils_param file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate soils_param_file
       IF (printlev>=2) WRITE(numout,*) 'Reading of soils_param file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("soils_param_file",enabled=.FALSE.)
       CALL xios_orchidee_set_fieldgroup_attr("soil_text",enabled=.FALSE.)
    END IF

    !! See commented part below for the reading of params_sp_mip.nc if spmipexp='maps'
    !! (with a bug, but helpful)
    
    !! 2. Prepare for reading of PFTmap file
    filename = 'PFTmap.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("PFTmap_file",name=name)

    ! Get veget_update from run.def needed to know if the file needs to be read
    veget_update=0
    WRITE(veget_str,'(a)') '0Y'
    CALL getin_p('VEGET_UPDATE', veget_str)
    l=INDEX(TRIM(veget_str),'Y')
    READ(veget_str(1:(l-1)),"(I2.2)") veget_update


    ! Check if PFTmap file will be read by XIOS in this execution
    IF ( xios_interpolation .AND. .NOT. impveg .AND. &
         (veget_update>0 .OR. restname_in=='NONE')) THEN
       ! PFTmap will not be read if impveg=TRUE
       ! PFTmap file will be read each year if veget_update>0 
       ! PFTmap is read if the restart file do not exist and if impveg=F

       ! Reading will be done
       IF (printlev>=2) WRITE(numout,*) 'Reading of PFTmap file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate PFTmap file
       IF (printlev>=2) WRITE(numout,*) 'Reading of PFTmap file will not be done with XIOS.'
       
       CALL xios_orchidee_set_file_attr("PFTmap_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("frac_veget",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("frac_veget_frac",enabled=.FALSE.)
    ENDIF
    

    !! 3. Prepare for reading of topography file
    filename = 'cartepente2d_15min.nc'
    CALL getin_p('TOPOGRAPHY_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("topography_file",name=name)
    
    ! Set default values used by XIOS for the interpolation
    slope_noreinf = 0.5 ! slope in percent
    CALL getin_p('SLOPE_NOREINF',slope_noreinf)
    lerr=xios_orchidee_setvar('slope_noreinf',slope_noreinf)
    lerr=xios_orchidee_setvar('slope_default',slope_default)
    
    get_slope = .FALSE.
    CALL getin_p('GET_SLOPE',get_slope)
    IF (xios_interpolation .AND. (restname_in=='NONE' .OR. get_slope)) THEN
       ! The slope file will be read using XIOS
       IF (printlev>=2) WRITE(numout,*) 'Reading of albedo file will be done later using XIOS. The filename is ', filename
    ELSE
       ! Deactivate slope reading
       IF (printlev>=2) WRITE(numout,*) 'The slope file will not be read by XIOS'
       CALL xios_orchidee_set_file_attr("topography_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("frac_slope_interp",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("reinf_slope_interp",enabled=.FALSE.)
    END IF
     
    !! 4. Prepare for reading of lai file
    filename = 'lai2D.nc'
    CALL getin_p('LAI_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("lai_file",name=name)
    ! Determine if lai file will be read. If not, deactivate the file.    
    IF (xios_interpolation .AND. restname_in=='NONE' .AND. read_lai) THEN
       ! Reading will be done
       IF (printlev>=2) WRITE(numout,*) 'Reading of lai file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate lai file
       IF (printlev>=2) WRITE(numout,*) 'Reading of lai file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("lai_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("frac_lai_interp",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("lai_interp",enabled=.FALSE.)
    END IF
    
    !! 5. Prepare for reading of woodharvest file
    filename = 'woodharvest.nc'
    CALL getin_p('WOODHARVEST_FILE',filename)
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("woodharvest_file",name=name)
    
    IF (xios_interpolation .AND. do_wood_harvest .AND. &
         (veget_update>0 .OR. restname_in=='NONE' )) THEN
       ! Woodharvest file will be read each year if veget_update>0 or if no restart file exists

       ! Reading will be done
       IF (printlev>=2) WRITE(numout,*) 'Reading of woodharvest file will be done later using XIOS. The filename is ', filename
    ELSE
       ! No reading, deactivate woodharvest file
       IF (printlev>=2) WRITE(numout,*) 'Reading of woodharvest file will not be done with XIOS.'
       CALL xios_orchidee_set_file_attr("woodharvest_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("woodharvest_interp",enabled=.FALSE.)
    ENDIF

    !! This part was introduced to prepare the reading of params_sp_mip.nc if spmipexp='maps'
    !! but there are mistakes in the IF ELSE ENDIF and we go through ELSE
    !! each time xios_interpolation = T, even if we don't need to read this file
    !! and it is not provided by sechiba.card
    !! The corresponding part in context_input_orchidee.xml is also commented
    
!!$    !! 6. Prepare for reading of soil parameter files
!!$
!!$    ! Get the file name from run.def file and set file attributes accordingly
!!$    filename = 'params_sp_mip.nc'
!!$    CALL getin_p('PARAM_SPMIP_FILE',filename)
!!$    name = filename(1:LEN_TRIM(FILENAME)-3)
!!$    CALL xios_orchidee_set_file_attr("soilparam_file",name=name)
!!$    ! Determine if the file will be read by XIOS. If not, deactivate reading of the file.
!!$    IF (xios_interpolation .AND. restname_in=='NONE' .AND. .NOT. impsoilt) THEN
!!$       ! Reading will be done with XIOS later
!!$       IF (printlev>=2) WRITE(numout,*) 'Reading of soil hydraulic parameters file will be done later using XIOS. The filename is ', filename
!!$    ELSE
!!$       ! No reading by XIOS, deactivate soilparam_file and related variables declared in context_input_orchidee.xml.
!!$       ! If this is not done, the model will crash if the file is not available in the run directory.
!!$       IF (printlev>=2) WRITE(numout,*) 'Reading of soil parameter file will not be done with XIOS.'
!!$       CALL xios_orchidee_set_file_attr("soilparam_file",enabled=.FALSE.)
!!$       CALL xios_orchidee_set_field_attr("soilks",enabled=.FALSE.)
!!$       CALL xios_orchidee_set_field_attr("soilnvan",enabled=.FALSE.)
!!$       CALL xios_orchidee_set_field_attr("soilavan",enabled=.FALSE.)
!!$       CALL xios_orchidee_set_field_attr("soilmcr",enabled=.FALSE.)
!!$       CALL xios_orchidee_set_field_attr("soilmcs",enabled=.FALSE.)
!!$       CALL xios_orchidee_set_field_attr("soilmcfc",enabled=.FALSE.)
!!$       CALL xios_orchidee_set_field_attr("soilmcw",enabled=.FALSE.)
!!$    ENDIF

    IF (printlev_loc>=3) WRITE(numout,*) 'End slowproc_xios_intialize'
   
  END SUBROUTINE slowproc_xios_initialize


!! ================================================================================================================================
!! SUBROUTINE 	: slowproc_initialize
!!
!>\BRIEF         Initialize slowproc module and call initialization of stomate module
!!
!! DESCRIPTION : Allocate module variables, read from restart file or initialize with default values
!!               Call initialization of stomate module.
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_initialize (kjit,          kjpij,        kjpindex,                          &
                                  rest_id,       rest_id_stom, hist_id_stom,   hist_id_stom_IPCC, &
                                  IndexLand,     indexveg,     lalo,           neighbours,        &
                                  resolution,    contfrac,     temp_air,                          &
                                  soiltile,      reinf_slope,  ks,             nvan,              &
                                  avan,          mcr,          mcs,            mcfc,              &
                                  mcw,           deadleaf_cover,               assim_param,       &
                                  lai,           frac_age,     height,         veget,             &
                                  frac_nobio,    njsc,         veget_max,      fraclut,           &
                                  nwdfraclut,    tot_bare_soil,totfrac_nobio,  qsintmax,          &
                                  temp_growth,   irrigated_next, irrig_frac_next, fraction_aeirrig_sw, &
                                  reinf_slope_soil)

!! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                          :: kjit                !! Time step number
    INTEGER(i_std), INTENT(in)                          :: kjpij               !! Total size of the un-compressed grid
    INTEGER(i_std),INTENT(in)                           :: kjpindex            !! Domain size - terrestrial pixels only
    INTEGER(i_std),INTENT (in)                          :: rest_id             !! Restart file identifier
    INTEGER(i_std),INTENT (in)                          :: rest_id_stom        !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist_id_stom        !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                           :: hist_id_stom_IPCC   !! STOMATE's IPCC _history_ file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: IndexLand           !! Indices of the points on the land map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg            !! Indices of the points on the vegetation (3D map ???) 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo                !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours     !! neighbouring grid points if land.
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution          !! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: contfrac            !! Fraction of continent in the grid (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: temp_air            !! Air temperature at first atmospheric model layer (K)
    
!! 0.2 Output variables 
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: temp_growth    !! Growth temperature (°C) - Is equal to t2m_month 
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(out)       :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: lai            !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: height         !! height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(out):: frac_age   !! Age efficacity from STOMATE for isoprene
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: veget          !! Fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out)  :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: veget_max      !! Maximum fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: tot_bare_soil  !! Total evaporating bare soil fraction in the mesh  (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh  (unitless)
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)    :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)    :: fraclut        !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)    :: nwdFraclut     !! Fraction of non-woody vegetation in each landuse tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: reinf_slope    !! slope coef for reinfiltration
   REAL(r_std),DIMENSION (kjpindex, nstm), INTENT(out)     :: reinf_slope_soil  !! slope coef for reinfiltration per soil tile
    !Salma: adding soil params
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: ks             !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: nvan           !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: avan           !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcr            !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcs            !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcfc           !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcw            !! Volumetric water content at wilting point (m^{3} m^{-3})

    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (out):: assim_param    !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: deadleaf_cover !! Fraction of soil covered by dead leaves (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)     :: qsintmax       !! Maximum water storage on vegetation from interception (mm)
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: irrigated_next !! Dynamic irrig. area, calculated in slowproc and passed to routing
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: irrig_frac_next!! Dynamic irrig. fraction, calculated in slowproc and passed to routing
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)          :: fraction_aeirrig_sw  !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                                   !! 1.0 here corresponds to fraction of irrigated area, not grid cell
!! 0.3 Local variables
    INTEGER(i_std)                                     :: ji, jsl
    REAL(r_std),DIMENSION (kjpindex,nslm)              :: land_frac         !! To ouput the clay/sand/silt fractions with a vertical dim

!_ ================================================================================================================================

    !! 1. Perform the allocation of all variables, define some files and some flags. 
    !     Restart file read for Sechiba.
    CALL slowproc_init (kjit, kjpindex, IndexLand, lalo, neighbours, resolution, contfrac, &
         rest_id, lai, frac_age, veget, frac_nobio, totfrac_nobio, soiltile, fraclut, nwdfraclut, reinf_slope, &
         ks,  nvan, avan, mcr, mcs, mcfc, mcw, &
         veget_max, tot_bare_soil, njsc, &
         height, lcanop, veget_update, veget_year, fraction_aeirrig_sw)
    

    !! 2. Define Time step in days for stomate
    dt_days = dt_stomate / one_day
    

    !! 3. check time step coherence between slow processes and fast processes
    IF ( dt_stomate .LT. dt_sechiba ) THEN
       WRITE(numout,*) 'slow_processes: time step smaller than forcing time step, dt_sechiba=',dt_sechiba,' dt_stomate=',dt_stomate
       CALL ipslerr_p(3,'slowproc_initialize','Coherence problem between dt_stomate and dt_sechiba',&
            'Time step smaller than forcing time step','')
    ENDIF
    
    !! 4. Call stomate to initialize all variables manadged in stomate,
    IF ( ok_stomate ) THEN

       CALL stomate_initialize &
            (kjit,           kjpij,                  kjpindex,                        &
             rest_id_stom,   hist_id_stom,           hist_id_stom_IPCC,               &
             indexLand,      lalo,                   neighbours,   resolution,        &
             contfrac,       totfrac_nobio,          clayfraction, temp_air,          &
             lai,            veget,                  veget_max,                       &
             deadleaf_cover,         assim_param,  temp_growth )
    ENDIF
    
    !! 5. Specific run without the carbon cycle (STOMATE not called): 
    !!     Need to initialize some variables that will be used in SECHIBA:
    !!     height, deadleaf_cover, assim_param, qsintmax.
    IF (.NOT. ok_stomate ) THEN
       CALL slowproc_derivvar (kjpindex, veget, lai, &
            qsintmax, deadleaf_cover, assim_param, height, temp_growth)
    ELSE
       qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
       qsintmax(:,1) = zero
    ENDIF

    !! 5.1  Dynamic irrigation maps as output to sechiba_end
    irrigated_next(:) = zero
    irrig_frac_next(:) = zero
    IF (do_irrigation ) THEN
      irrigated_next(:) = irrigated_new(:)
      ! irrig_frac calculation
      DO ji=1,kjpindex
          IF( (resolution(ji,1)*resolution(ji,2)*contfrac(ji) ) > min_sechiba) THEN
            !SUM(routing_area(ig,:)) is totarea(ig) = m2
            irrig_frac_next(ji) = MIN( soiltile(ji,irrig_st) * SUM(veget_max(ji,:)) , &
                irrigated_new(ji) / (resolution(ji,1)*resolution(ji,2)*contfrac(ji) ) )

            ! Soiltile(fraction of vegtot) * SUM(veget_max) [SUM(veget_max) is vegtot in hydrol, calculated
            ! differently in routing] = fraction of grid cell
            ! irrigated(m2)/ SUM(routing_area(ig,:) = Fraction of grid cell
            ! irrig_frac is always fraction of grid cell
          ENDIF
      ENDDO
    END IF

    !! 5.2 Calculation of reinf_slope_soil to pass to hydrol
    reinf_slope_soil(:,:) = zero
    DO ji=1,kjpindex
        reinf_slope_soil(ji,:) = reinf_slope(ji)
        IF( Reinfiltr_IrrigField .AND.  irrig_frac_next(ji) > min_sechiba ) THEN
          reinf_slope_soil(ji, irrig_st) = MAX(reinf_slope(ji), reinf_slope_cropParam)
        ENDIF
    ENDDO

    !! 6. Output with XIOS for variables done only once per run

    DO jsl=1,nslm
       land_frac(:,jsl) = clayfraction(:)
    ENDDO
    CALL xios_orchidee_send_field("clayfraction",land_frac) ! mean fraction of clay in grid-cell
    DO jsl=1,nslm
       land_frac(:,jsl) = sandfraction(:)
    ENDDO
    CALL xios_orchidee_send_field("sandfraction",land_frac) ! mean fraction of sand in grid-cell
    DO jsl=1,nslm
       land_frac(:,jsl) = siltfraction(:)
    ENDDO
    CALL xios_orchidee_send_field("siltfraction",land_frac) ! mean fraction of silt in grid-cell
    
  END SUBROUTINE slowproc_initialize


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_main
!!
!>\BRIEF         Main routine that manage variable initialisation (slowproc_init), 
!! prepare the restart file with the slowproc variables, update the time variables 
!! for slow processes, and possibly update the vegetation cover, before calling 
!! STOMATE in the case of the carbon cycle activated or just update LAI (and possibly
!! the vegetation cover) for simulation with only SECHIBA   
!!
!!
!! DESCRIPTION  : (definitions, functional, design, flags): The subroutine manages 
!! diverses tasks:
!! (1) Initializing all variables of slowproc (first call)
!! (2) Preparation of the restart file for the next simulation with all prognostic variables
!! (3) Compute and update time variable for slow processes
!! (4) Update the vegetation cover if there is some land use change (only every years)
!! (5) Call STOMATE for the runs with the carbone cycle activated (ok_stomate) and compute the respiration
!!     and the net primary production
!! (6) Compute the LAI and possibly update the vegetation cover for run without STOMATE 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S):  ::co2_flux, ::fco2_lu,::fco2_wh, ::fco2_ha, ::lai, ::height, ::veget, ::frac_nobio,  
!! ::veget_max, ::woodharvest, ::totfrac_nobio, ::soiltype, ::assim_param, ::deadleaf_cover, ::qsintmax,
!! and resp_maint, resp_hetero, resp_growth, npp that are calculated and stored
!! in stomate is activated.  
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : 
! \latexonly 
! \includegraphics(scale=0.5){SlowprocMainFlow.eps} !PP to be finalize!!)
! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_main (kjit, kjpij, kjpindex, &
       IndexLand, indexveg, lalo, neighbours, resolution, contfrac, soiltile, fraclut, nwdFraclut, &
       temp_air, temp_sol, stempdiag, &
       humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
       deadleaf_cover, &
       assim_param, &
       lai, frac_age, height, veget, frac_nobio, veget_max, totfrac_nobio, qsintmax, &
       rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
       co2_flux, fco2_lu, fco2_wh, fco2_ha, &
       temp_growth, tot_bare_soil, &
       irrigated_next, irrig_frac_next, reinf_slope, reinf_slope_soil)
  
!! INTERFACE DESCRIPTION

!! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                          :: kjit                !! Time step number
    INTEGER(i_std), INTENT(in)                          :: kjpij               !! Total size of the un-compressed grid
    INTEGER(i_std),INTENT(in)                           :: kjpindex            !! Domain size - terrestrial pixels only
    INTEGER(i_std),INTENT (in)                          :: rest_id,hist_id     !! _Restart_ file and _history_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist2_id            !! _history_ file 2 identifier
    INTEGER(i_std),INTENT (in)                          :: rest_id_stom        !! STOMATE's _Restart_ file identifier
    INTEGER(i_std),INTENT (in)                          :: hist_id_stom        !! STOMATE's _history_ file identifier
    INTEGER(i_std),INTENT(in)                           :: hist_id_stom_IPCC   !! STOMATE's IPCC _history_ file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: IndexLand           !! Indices of the points on the land map
    INTEGER(i_std),DIMENSION (kjpindex*nvm), INTENT (in):: indexveg            !! Indices of the points on the vegetation (3D map ???) 
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo                !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in)  :: neighbours   !! neighbouring grid points if land
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution          !! size in x an y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: contfrac            !! Fraction of continent in the grid (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT (in)  :: humrel              !! Relative humidity ("moisture stress") (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex), INTENT(in)        :: temp_air            !! Temperature of first model layer (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: temp_sol            !! Surface temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)  :: stempdiag           !! Soil temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)  :: shumdiag            !! Relative soil moisture (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: litterhumdiag       !! Litter humidity  (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: precip_rain         !! Rain precipitation (mm dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: precip_snow         !! Snow precipitation (mm dt_stomate^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)    :: gpp                 !! GPP of total ground area (gC m^{-2} time step^{-1}). 
                                                                               !! Calculated in sechiba, account for vegetation cover and 
                                                                               !! effective time step to obtain gpp_d   
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)        :: reinf_slope         !! slope coef for reinfiltration

!! 0.2 Output variables 
    REAL(r_std), DIMENSION (kjpindex,nvm), INTENT(out)  :: co2_flux            !! CO2 flux per average ground area (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: fco2_lu             !! CO2 flux from land-use (without forest management) (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: fco2_wh             !! CO2 Flux to Atmosphere from Wood Harvesting (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: fco2_ha             !! CO2 Flux to Atmosphere from Crop Harvesting (gC m^{-2} dt_stomate^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)      :: temp_growth         !! Growth temperature (°C) - Is equal to t2m_month 
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)      :: tot_bare_soil       !! Total evaporating bare soil fraction in the mesh
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)      :: irrigated_next      !! Dynamic irrig. area, calculated in slowproc and passed to routing

!! 0.3 Modified variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: lai            !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: height         !! height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(inout):: frac_age   !! Age efficacity from STOMATE for isoprene
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: veget          !! Fraction of vegetation type including none biological fractionin the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (inout)  :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: veget_max      !! Maximum fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(inout)    :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(inout)    :: fraclut        !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(inout)    :: nwdFraclut     !! Fraction of non-woody vegetation in each landuse tile (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (inout):: assim_param    !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: deadleaf_cover !! Fraction of soil covered by dead leaves (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (inout)     :: qsintmax       !! Maximum water storage on vegetation from interception (mm)
    REAL(r_std), DIMENSION (kjpindex), INTENT (inout)        :: irrig_frac_next   !! Dynamic irrig. fraction, calculated in slowproc and passed to routing
    REAL(r_std),DIMENSION (kjpindex, nstm), INTENT(inout)    :: reinf_slope_soil  !!  slope coef for reinfiltration per soil tile

!! 0.4 Local variables
    INTEGER(i_std)                                     :: j, jv, ji            !! indices 
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_maint           !! Maitanance component of autotrophic respiration in (gC m^{-2} dt_stomate^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_hetero          !! heterotrophic resp. (gC/(m**2 of total ground)/time step)
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: resp_growth          !! Growth component of autotrophic respiration in gC m^{-2} dt_stomate^{-1})
    REAL(r_std), DIMENSION(kjpindex,nvm)               :: npp                  !! Net Ecosystem Exchange (gC/(m**2 of total ground)/time step)
    REAL(r_std),DIMENSION (kjpindex)                   :: totfrac_nobio_new    !! Total fraction for the next year 
    REAL(r_std),DIMENSION (kjpindex)                   :: histvar              !! Temporary variable for output

!_ ================================================================================================================================

    !! 1. Compute and update all variables linked to the date and time
    IF (printlev_loc>=5) WRITE(numout,*) 'Entering slowproc_main, year_start, month_start, day_start, sec_start=',&
         year_start, month_start,day_start,sec_start   
 
    !! 2. Activate slow processes if it is the end of the day
    IF ( LastTsDay ) THEN
       ! 3.2.2 Activate slow processes in the end of the day
       do_slow = .TRUE.
       
       ! 3.2.3 Count the number of days 
       days_since_beg = days_since_beg + 1
       IF (printlev_loc>=4) WRITE(numout,*) "New days_since_beg : ",days_since_beg
    ELSE
       do_slow = .FALSE.
    ENDIF

    !! 3. Update the vegetation if it is time to do so.
    !!    This is done at the first sechiba time step on a new year and only every "veget_update" years. 
    !!    veget_update correspond to a number of years between each vegetation updates.
    !!    Nothing is done if veget_update=0.
    !!    Update will never be done if impveg=true because veget_update=0.
    IF ( FirstTsYear ) THEN
       IF (veget_update > 0) THEN
          veget_year = veget_year + 1
   
          ! Update of the vegetation cover with Land Use only if 
          ! the current year match the requested condition (a multiple of "veget_update")
          IF ( MOD(veget_year - veget_year_orig, veget_update) == 0 ) THEN
             IF (printlev_loc>=1) WRITE(numout,*)  'We are updating the vegetation map for year =' , veget_year
             
             ! Read the new the vegetation from file. Output is veget_max_new and frac_nobio_new
             CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
                  veget_max, veget_max_new, frac_nobio_new, veget_year, .FALSE.)
             
             IF (do_wood_harvest) THEN
                ! Read the new the wood harvest map from file. Output is wood harvest
                CALL slowproc_woodharvest(kjpindex, lalo, neighbours, resolution, contfrac, woodharvest)
             ENDIF
   
             ! Set the flag do_now_stomate_lcchange to activate stomate_lcchange.
             ! This flag will be kept to true until stomate_lcchange has been done. 
             ! The variable totfrac_nobio_new will only be used in stomate when this flag is activated
             do_now_stomate_lcchange=.TRUE.
             IF ( .NOT. ok_stomate ) THEN
                ! Special case if stomate is not activated : set the variable done_stomate_lcchange=true 
                ! so that the subroutine slowproc_change_frac will be called in the end of sechiba_main.
                done_stomate_lcchange=.TRUE.
             END IF
          ENDIF
       ENDIF
       IF ( do_wood_harvest) THEN
          ! Set the flag do_now_stomate_woodharvest to activate stomate_woodharvest.
          ! This flag will be kept to true until stomate_woodharvest has been done. 
          do_now_stomate_woodharvest=.TRUE.
       ENDIF
    ENDIF

    !! 4. Main call to STOMATE
    IF ( ok_stomate ) THEN

       ! Calculate totfrac_nobio_new only for the case when the land use map has been read previously
       IF (do_now_stomate_lcchange) THEN 
          totfrac_nobio_new(:) = zero 
          DO jv = 1, nnobio 
             totfrac_nobio_new(:) = totfrac_nobio_new(:) + frac_nobio_new(:,jv)
          ENDDO 
       ELSE 
          totfrac_nobio_new(:) = zero 
       END IF 

       !! 4.1 Call stomate main routine that will call all c-cycle routines       !
       CALL stomate_main (kjit, kjpij, kjpindex, &
            IndexLand, lalo, neighbours, resolution, contfrac, totfrac_nobio, clayfraction, &
            temp_air, temp_sol, stempdiag, &
            humrel, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
            deadleaf_cover, &
            assim_param, &
            lai, frac_age, height, veget, veget_max, &
            veget_max_new, woodharvest, totfrac_nobio_new, fraclut, &
            rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
            co2_flux, fco2_lu, fco2_wh, fco2_ha, &
            resp_maint, resp_hetero, resp_growth, temp_growth)


       !! 4.2 Output the respiration terms and the net primary
       !!     production (NPP) that are calculated in STOMATE

       ! 4.2.1 Output the 3 respiration terms
       ! These variables could be output from stomate.
       ! Variables per pft
       CALL xios_orchidee_send_field("maint_resp",resp_maint/dt_sechiba)
       CALL xios_orchidee_send_field("hetero_resp",resp_hetero/dt_sechiba)
       CALL xios_orchidee_send_field("growth_resp",resp_growth/dt_sechiba)

       ! Variables on grid-cell
       CALL xios_orchidee_send_field("rh_ipcc2",SUM(resp_hetero,dim=2)/dt_sechiba)
       histvar(:)=zero
       DO jv = 2, nvm
          IF ( .NOT. is_tree(jv) .AND. natural(jv) ) THEN
             histvar(:) = histvar(:) + resp_hetero(:,jv)
          ENDIF
       ENDDO
       CALL xios_orchidee_send_field("rhGrass",histvar/dt_sechiba)

       histvar(:)=zero
       DO jv = 2, nvm
          IF ( (.NOT. is_tree(jv)) .AND. (.NOT. natural(jv)) ) THEN
             histvar(:) = histvar(:) + resp_hetero(:,jv)
          ENDIF
       ENDDO
       CALL xios_orchidee_send_field("rhCrop",histvar/dt_sechiba)

       histvar(:)=zero
       DO jv = 2, nvm
          IF ( is_tree(jv) ) THEN
             histvar(:) = histvar(:) + resp_hetero(:,jv)
          ENDIF
       ENDDO
       CALL xios_orchidee_send_field("rhTree",histvar/dt_sechiba)

       ! Output with IOIPSL
       CALL histwrite_p(hist_id, 'maint_resp', kjit, resp_maint, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'hetero_resp', kjit, resp_hetero, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'growth_resp', kjit, resp_growth, kjpindex*nvm, indexveg)
       
       ! 4.2.2 Compute the net primary production as the diff from
       ! Gross primary productin and the growth and maintenance respirations
       npp(:,1)=zero
       DO j = 2,nvm
          npp(:,j) = gpp(:,j) - resp_growth(:,j) - resp_maint(:,j)
       ENDDO
       
       CALL xios_orchidee_send_field("npp",npp/dt_sechiba)
       
       CALL histwrite_p(hist_id, 'npp', kjit, npp, kjpindex*nvm, indexveg)
       
       IF ( hist2_id > 0 ) THEN
          CALL histwrite_p(hist2_id, 'maint_resp', kjit, resp_maint, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'hetero_resp', kjit, resp_hetero, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'growth_resp', kjit, resp_growth, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'npp', kjit, npp, kjpindex*nvm, indexveg)
       ENDIF
     
    ELSE
       !! ok_stomate is not activated
       !! Define the CO2 flux from the grid point to zero (no carbone cycle)
       co2_flux(:,:) = zero
       fco2_lu(:) = zero
       fco2_wh(:) = zero
       fco2_ha(:) = zero
    ENDIF

 
    !! 5. Do daily processes if necessary
    !!
    IF ( do_slow ) THEN

       !!  5.1 Calculate the LAI if STOMATE is not activated
       IF ( .NOT. ok_stomate ) THEN
          CALL slowproc_lai (kjpindex, lcanop,stempdiag, &
               lalo,resolution,lai,laimap)
          
          frac_age(:,:,1) = un
          frac_age(:,:,2) = zero
          frac_age(:,:,3) = zero
          frac_age(:,:,4) = zero
       ENDIF

       !! 5.2 Update veget
       CALL slowproc_veget (kjpindex, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile, fraclut, nwdFraclut)

       !! 5.3 updates qsintmax and other derived variables
       IF ( .NOT. ok_stomate ) THEN
          CALL slowproc_derivvar (kjpindex, veget, lai, &
               qsintmax, deadleaf_cover, assim_param, height, temp_growth)
       ELSE
          qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)
          qsintmax(:,1) = zero
       ENDIF
    END IF

    !! 6. Calculate tot_bare_soil needed in hydrol, diffuco and condveg (fraction in the mesh)
    tot_bare_soil(:) = veget_max(:,1)
    DO jv = 2, nvm
       DO ji =1, kjpindex
          tot_bare_soil(ji) = tot_bare_soil(ji) + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
    END DO
    


    !! 6.1. Call to interpolation of dynamic irrigation map, if time to do so (as for vegetmax interpolation)
    !! Important difference: veget map is updated every veget_update years. Irrig map_pft
    !! is updated every year for now
    !! Put here to use updates values of soiltile and veget_max
    irrigated_next(:) = zero
    IF (do_irrigation .AND. irrig_map_dynamic_flag) THEN
      ! Attention: veget_year already updated, but veget_update must be >0, I.E. it must read veget maps
      ! Seems logic to update irrigation maps if vegetation maps are updated too
      IF ( FirstTsYear) THEN
        CALL slowproc_readirrigmap_dyn(kjpindex, lalo, neighbours,  resolution, contfrac,         &
             irrigated_new)
        DO ji=1,kjpindex
           IF( (resolution(ji,1)*resolution(ji,2)*contfrac(ji) ) > min_sechiba) THEN !Multiplication is total area(ig) = m2

             irrig_frac_next(ji) = MIN( soiltile(ji,irrig_st) * SUM(veget_max(ji,:)) , &
                 irrigated_new(ji) / (resolution(ji,1)*resolution(ji,2)*contfrac(ji) ) )
             ! soiltile(fraction of vegtot) * SUM(veget_max)  = fraction of grid cell
             ! irrigated(m2)/ grid_cell_area = Fraction of grid cell
             ! irrig_frac is always fraction of grid cell
           ENDIF
        ENDDO
      ENDIF
     !! Here irrigated_next from sechiba = irrigated_new from slowproc!!
     !! Dynamic irrigation maps as output to sechiba_end
    ENDIF
    irrigated_next(:) = irrigated_new(:)

    !!
    !! 6.2 Calculation of reinf_slope_soil to pass to hydrol
    IF (FirstTsYear) THEN
      reinf_slope_soil(:,:) = zero
      DO ji=1,kjpindex
          reinf_slope_soil(ji,:) = reinf_slope(ji)
          IF( Reinfiltr_IrrigField .AND.  irrig_frac_next(ji) > min_sechiba ) THEN
            reinf_slope_soil(ji, irrig_st) = MAX(reinf_slope(ji), reinf_slope_cropParam)
          ENDIF
      ENDDO
    ENDIF

    !! 7. Do some basic tests on the surface fractions updated above, only if
    !!    slowproc_veget has been done (do_slow). No change of the variables. 
    IF (do_slow) THEN
        CALL slowproc_checkveget(kjpindex, frac_nobio, veget_max, veget, tot_bare_soil, soiltile)
    END IF  

    !! 8. Write output fields
    CALL xios_orchidee_send_field("tot_bare_soil",tot_bare_soil)
    
    IF ( .NOT. almaoutput) THEN
       CALL histwrite_p(hist_id, 'tot_bare_soil', kjit, tot_bare_soil, kjpindex, IndexLand)
    END IF


    IF (printlev_loc>=3) WRITE (numout,*) ' slowproc_main done '

  END SUBROUTINE slowproc_main


!! ================================================================================================================================
!! SUBROUTINE 	: slowproc_finalize
!!
!>\BRIEF         Write to restart file variables for slowproc module and call finalization of stomate module
!!
!! DESCRIPTION :
!!
!! RECENT CHANGE(S): Add arrays of soil hydraulic parameters to the restart file.
!!                   Linked to SP-MIP project (Tafasca Salma and Ducharne Agnes).
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_finalize (kjit,       kjpindex,  rest_id,  IndexLand,  &
                                njsc,       lai,       height,   veget,      &
                                frac_nobio, veget_max, reinf_slope,          &
                                ks,  nvan, avan, mcr, mcs, mcfc, mcw,        &
                                assim_param, frac_age, fraction_aeirrig_sw)

!! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjit           !! Time step number
    INTEGER(i_std),INTENT(in)                            :: kjpindex       !! Domain size - terrestrial pixels only
    INTEGER(i_std),INTENT (in)                           :: rest_id        !! Restart file identifier
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)     :: IndexLand      !! Indices of the points on the land map
    INTEGER(i_std), DIMENSION(kjpindex), INTENT(in)      :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: lai            !! Leaf area index (m^2 m^{-2})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: height         !! height of vegetation (m)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget          !! Fraction of vegetation type including none biological fraction (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (in) :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)    :: veget_max      !! Maximum fraction of vegetation type including none biological fraction (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: reinf_slope    !! slope coef for reinfiltration
    !salma: added following soil params
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: ks             !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: nvan           !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: avan           !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: mcr            !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: mcs            !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: mcfc           !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: mcw            !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)         :: fraction_aeirrig_sw    !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                                   !! 1.0 here corresponds to fraction of irrigated area, not grid cell
    REAL(r_std),DIMENSION (kjpindex,nvm,npco2),INTENT (in):: assim_param  !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT(in):: frac_age  !! Age efficacity from STOMATE for isoprene
!! 0.4 Local variables
    REAL(r_std)                                          :: tmp_day(1)     !! temporary variable for I/O
    INTEGER                                              :: jf             !! Indice
    CHARACTER(LEN=4)                                     :: laistring      !! Temporary character string
    CHARACTER(LEN=80)                                    :: var_name       !! To store variables names for I/O
!_ ================================================================================================================================

    IF (printlev_loc>=3) WRITE (numout,*) 'Write restart file with SLOWPROC variables '

    ! 2.1 Write a series of variables controled by slowproc: day
    ! counter, vegetation fraction, max vegetation fraction, LAI
    ! variable from stomate, fraction of bare soil, soiltype
    ! fraction, clay fraction, height of vegetation, map of LAI
    
    CALL restput_p (rest_id, 'veget', nbp_glo, nvm, 1, kjit, veget, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'veget_max', nbp_glo, nvm, 1, kjit, veget_max, 'scatter',  nbp_glo, index_g)

    IF (do_wood_harvest) THEN
       CALL restput_p (rest_id, 'woodharvest', nbp_glo, 1, 1, kjit, woodharvest, 'scatter',  nbp_glo, index_g)
    END IF

    CALL restput_p (rest_id, 'lai', nbp_glo, nvm, 1, kjit, lai, 'scatter',  nbp_glo, index_g)

    CALL restput_p (rest_id, 'frac_nobio', nbp_glo, nnobio, 1, kjit, frac_nobio, 'scatter',  nbp_glo, index_g)

    IF ( do_irrigation ) THEN
          CALL restput_p (rest_id, 'irrigmap_dyn', nbp_glo, 1, 1, kjit, irrigated_new, 'scatter',  nbp_glo, index_g)
          IF ( select_source_irrig ) THEN
                CALL restput_p (rest_id, 'fraction_aeirrig_sw', nbp_glo, 1, 1, kjit, fraction_aeirrig_sw, 'scatter',  nbp_glo, index_g)
          ENDIF
    ENDIF

    IF ( do_imperviousness ) THEN
          CALL restput_p (rest_id, 'frac_imperv', nbp_glo, 1, 1, kjit, frac_imperv, 'scatter',  nbp_glo, index_g)
    ENDIF

    DO jf = 1, nleafages
       ! variable name is somewhat complicated as ioipsl does not allow 3d variables for the moment...
       WRITE(laistring,'(i4)') jf
       laistring=ADJUSTL(laistring)
       var_name='frac_age_'//laistring(1:LEN_TRIM(laistring))
       CALL restput_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, frac_age(:,:,jf), 'scatter',  nbp_glo, index_g)
    ENDDO

    ! Add the soil_classif as suffix for the variable name of njsc when it is stored in the restart file. 
    IF (soil_classif == 'zobler') THEN
       var_name= 'njsc_zobler'
    ELSE IF (soil_classif == 'usda') THEN
       var_name= 'njsc_usda'
    END IF
    CALL restput_p (rest_id, var_name, nbp_glo, 1, 1, kjit, REAL(njsc, r_std), 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'reinf_slope', nbp_glo, 1, 1, kjit, reinf_slope, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'clay_frac', nbp_glo, 1, 1, kjit, clayfraction, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'sand_frac', nbp_glo, 1, 1, kjit, sandfraction, 'scatter',  nbp_glo, index_g)
    !salma: added the following lines for restput of the soil parameters
    CALL restput_p (rest_id, 'ks', nbp_glo, 1, 1, kjit, ks, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'mcs', nbp_glo, 1, 1, kjit, mcs, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'mcr', nbp_glo, 1, 1, kjit, mcr, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'mcw', nbp_glo, 1, 1, kjit, mcw, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'mcfc', nbp_glo, 1, 1, kjit, mcfc, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'nvan', nbp_glo, 1, 1, kjit, nvan, 'scatter',  nbp_glo, index_g)
    CALL restput_p (rest_id, 'avan', nbp_glo, 1, 1, kjit, avan, 'scatter',  nbp_glo, index_g)

    ! The height of the vegetation could in principle be recalculated at the beginning of the run.
    ! However, this is very tedious, as many special cases have to be taken into account. This variable
    ! is therefore saved in the restart file.
    CALL restput_p (rest_id, 'height', nbp_glo, nvm, 1, kjit, height, 'scatter',  nbp_glo, index_g)
    !
    ! Specific case where the LAI is read and not calculated by STOMATE: need to be saved
    IF (read_lai) THEN     
       CALL restput_p (rest_id, 'laimap', nbp_glo, nvm, 12, kjit, laimap)
    ENDIF
    !
    ! If there is some land use change, write the year for the land use ??? 
    tmp_day(1) = REAL(veget_year,r_std)
    IF (is_root_prc) CALL restput (rest_id, 'veget_year', 1 , 1  , 1, kjit, tmp_day)
        
    ! 2.2 Write restart variables managed by STOMATE
    IF ( ok_stomate ) THEN
       CALL stomate_finalize (kjit, kjpindex, indexLand, clayfraction, assim_param) 
    ENDIF
    
  END SUBROUTINE slowproc_finalize


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_init
!!
!>\BRIEF         Initialisation of all variables linked to SLOWPROC
!!
!! DESCRIPTION  : (definitions, functional, design, flags): The subroutine manages 
!! diverses tasks:
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::lcanop, ::veget_update, ::veget_year,
!! ::lai, ::veget, ::frac_nobio, ::totfrac_nobio, ::veget_max, ::height, ::soiltype
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_init (kjit, kjpindex, IndexLand, lalo, neighbours, resolution, contfrac, &
       rest_id, lai, frac_age, veget, frac_nobio, totfrac_nobio, soiltile, fraclut, nwdfraclut, reinf_slope, &
       ks,  nvan, avan, mcr, mcs, mcfc, mcw, &
       veget_max, tot_bare_soil, njsc, &
       height, lcanop, veget_update, veget_year, fraction_aeirrig_sw)
    
    !! INTERFACE DESCRIPTION

    !! 0.1 Input variables
    INTEGER(i_std), INTENT (in)                           :: kjit           !! Time step number
    INTEGER(i_std), INTENT (in)                           :: kjpindex       !! Domain size - Terrestrial pixels only 
    INTEGER(i_std), INTENT (in)                           :: rest_id        !! Restart file identifier
    
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: IndexLand      !! Indices of the land points on the map
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)       :: lalo           !! Geogr. coordinates (latitude,longitude) (degrees)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours  !! Vector of neighbours for each grid point
                                                                            !! (1=North and then clockwise)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)       :: resolution     !! size in x and y of the grid (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: contfrac       !! Fraction of continent in the grid (unitless)
    
    !! 0.2 Output variables
    INTEGER(i_std), INTENT(out)                           :: lcanop         !! Number of Canopy level used to compute LAI
    INTEGER(i_std), INTENT(out)                           :: veget_update   !! update frequency in timesteps (years) for landuse
    INTEGER(i_std), INTENT(out)                           :: veget_year     !! first year for landuse   (year or index ???)
    
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: lai            !! Leaf Area index (m^2 / m^2)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: veget          !! Fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT (out) :: frac_nobio     !! Fraction of ice,lakes,cities, ... in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: totfrac_nobio  !! Total fraction of ice+lakes+cities+... in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: veget_max      !! Max fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: tot_bare_soil  !! Total evaporating bare soil fraction in the mesh
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)    :: height         !! Height of vegetation or surface in genral ??? (m)
    REAL(r_std),DIMENSION (kjpindex,nvm,nleafages), INTENT (out):: frac_age !! Age efficacity from STOMATE for isoprene
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)   :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)   :: fraclut        !! Fraction of each landuse tile
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)   :: nwdfraclut     !! Fraction of non woody vegetation in each landuse tile
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)        :: reinf_slope    !! slope coef for reinfiltration 
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: ks             !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: nvan           !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: avan           !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcr            !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcs            !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcfc           !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: mcw            !! Volumetric water content at wilting point (m^{3} m^{-3})
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: fraction_aeirrig_sw    !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                                     !! 1.0 here corresponds to fraction of irrig. area, not grid cell

    INTEGER(i_std), DIMENSION(kjpindex), INTENT(out)      :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    
    !! 0.3 Local variables
    INTEGER(i_std)                                        :: ji, jv, ier,jst   !! Indices  
    INTEGER(i_std)                                        :: l, jf             !! Indices  
    INTEGER(i_std)                                        :: vtmp(1)           !! temporary variable
    INTEGER(i_std)                                        :: njsc_imp          !! njsc to impose nvan, ks etc. if impsoil
    REAL(r_std)                                           :: tmp_veget_year(1) !! temporary variable
    REAL(r_std)                                           :: zcanop            !! ???? soil depth taken for canopy
    REAL(r_std), DIMENSION(nslm)                          :: zsoil             !! soil depths at diagnostic levels
    REAL(r_std)                                           :: frac_nobio1       !! temporary variable for frac_nobio(see above)
    REAL(r_std), DIMENSION(kjpindex)                      :: tmp_real
    REAL(r_std), DIMENSION(kjpindex,nslm)                 :: stempdiag2_bid    !! matrix to store stempdiag_bid
    REAL(r_std), DIMENSION (kjpindex,nscm)                :: soilclass         !! Fractions of each soil textural class in the grid cell (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex)                      :: frac_crop_tot     !! Total fraction occupied by crops (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex)                      :: mvan, psi_fc, psi_w  !! To calculate default wilting point and 
                                                                                  !! field capacity if impsoilt
    REAL(r_std), DIMENSION(kjpindex)                      :: mcfc_default      !! Default field capacity if impsoilt
    REAL(r_std), DIMENSION(kjpindex)                      :: mcw_default       !! Default wilting point if impsoilt
    REAL(r_std)                                           :: nvan_default      !! Default  if impsoilt
    REAL(r_std)                                           :: avan_default      !! Default  if impsoilt
    REAL(r_std)                                           :: mcr_default       !! Default  if impsoilt
    REAL(r_std)                                           :: mcs_default       !! Default  if impsoilt
    REAL(r_std)                                           :: ks_default        !! Default  if impsoilt
    REAL(r_std)                                           :: clayfraction_default  !! Default  if impsoilt
    REAL(r_std)                                           :: sandfraction_default  !! Default  if impsoilt
    REAL(r_std), DIMENSION(nbpt)                          :: height_buidings       !! fraction of imperviousness from map
    CHARACTER(LEN=4)                                      :: laistring         !! Temporary character string
    CHARACTER(LEN=80)                                     :: var_name          !! To store variables names for I/O
    CHARACTER(LEN=30), SAVE                               :: veget_str         !! update frequency for landuse
!$OMP THREADPRIVATE(veget_str)
    LOGICAL                                               :: get_slope
    LOGICAL                                               :: found_restart     !! found_restart=true if all 3 variables veget_max,  
                                                                               !! veget and frac_nobio are read from restart file
    LOGICAL                                               :: call_slowproc_soilt !! This variables will be true if subroutine
                                                                                 !! slowproc_soilt needs to be called
    !_ ================================================================================================================================
  
    !! 0. Initialize local printlev
    printlev_loc=get_printlev('slowproc')
    IF (printlev_loc>=3) WRITE (numout,*) "In slowproc_init"
    
    
    !! 1. Allocation 

    ALLOCATE (clayfraction(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable clayfraction','','')
    clayfraction(:)=undef_sechiba
    
    ALLOCATE (sandfraction(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable sandfraction','','')
    sandfraction(:)=undef_sechiba
    
    ALLOCATE (siltfraction(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable siltfraction','','')
    siltfraction(:)=undef_sechiba

    ! Allocation of last year vegetation fraction in case of land use change
    ALLOCATE(veget_max_new(kjpindex, nvm), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable veget_max_new','','')

    ! Allocation of last irrig map area in case of change
    ALLOCATE(irrigated_new(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable irrigated_new','','')

    ! Allocation of wood harvest
    ALLOCATE(woodharvest(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable woodharvest','','')

    ! Allocation of the fraction of non biospheric areas 
    ALLOCATE(frac_nobio_new(kjpindex, nnobio), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable frac_nobio_new','','')
    
    ! Allocate laimap
    IF (read_lai)THEN
       ALLOCATE (laimap(kjpindex,nvm,12),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable laimap','','')
    ELSE
       ALLOCATE (laimap(1,1,1), stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable laimap(1,1,1)','','')
    ENDIF

    ! Allocation of impervious fraction
    ALLOCATE(frac_imperv(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable frac_imperv','','')
    ! Allocation of the coefficient of imperviousness to modify ks(texture )
    ALLOCATE(coeff_imperv(kjpindex), STAT=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'slowproc_init','Problem in allocation of variable coeff_imperv','','')
    !! 2. Read soil related variables
    ! Following the trunk, we remove the dependance of impsoilt to impveg; impsoilt overrules the restart

    !! 2.a Looking first in the restart files

    ! The variables are read from restart file. If at least one of the related variables are not found, the 
    ! variable call_slowproc_soilt will be true and the variables will be read and interpolated from file in 
    ! subroutine slowproc_soilt.
    call_slowproc_soilt=.FALSE.    
        
    ! First we define the soil texture of the grid-cells

    ! Add the soil_classif as suffix for the variable name of njsc when it is stored in the restart file.
    ! AD: is this still usefull now that we use the same indexing for both kinds of texture maps??
    ! AD - Warning: if we read an "old" restart from a Zobler map with njsc in 1,2,3 instead of 3,6,9,
    !      but not so old with all the restart variables incl. ks, nvan, etc., we will badly interpret the
    !      Coarse, Medium, Fine textures of Zobler as Sandy, Loamy Sand, Sandy Loam.
    ! Fortunately, it's not likely, it only concerns runs using brach2.2 between revisions 7325 and 7337.
    ! We keep this to easily identfy if a run made with Zobler or Reynolds soil map based on restart files
    
    IF (soil_classif == 'zobler') THEN
       var_name= 'njsc_zobler'
    ELSE IF (soil_classif == 'usda') THEN
       var_name= 'njsc_usda'
    ELSE
       CALL ipslerr_p(3,'slowproc_init','Non supported soil type classification','','')
    END IF

    ! Index of the dominant soil type in the grid cell
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Index of soil type')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., tmp_real, "gather", nbp_glo, index_g)
    IF ( ALL( tmp_real(:) .EQ. val_exp) ) THEN
       njsc (:) = undef_int
       call_slowproc_soilt=.TRUE.
    ELSEIF ( ANY( ABS(tmp_real(:)) .GT. 20 ) ) THEN
       ! This is a test if some of the values for njsc are out of rang. This can be the case if
       ! there is a mismatch in the land-sea mask. njsc is the first
       ! restart variable read in current version of ORCHIDEE.
       CALL ipslerr_p(3, 'slowproc_init', 'Some values for njsc from restart file are out of range.',&
                         'There is probably a mismatch in the land/sea mask in the model and in the restart file.',&
                         'Start the model again without restart files for ORCHIDEE.')
    ELSE
       njsc = NINT(tmp_real)
    END IF

    ! Other soil parameters
    
    var_name= 'ks'
    CALL ioconf_setatt_p('UNITS', 'mm/d')
    CALL ioconf_setatt_p('LONG_NAME','Soil saturated water content')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., ks, "gather", nbp_glo, index_g)
    IF ( ALL(ks(:) .EQ. val_exp ) ) THEN
       ! ks is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    var_name= 'mcs'
    CALL ioconf_setatt_p('UNITS', 'm3/m3')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., mcs, "gather", nbp_glo, index_g)
    IF ( ALL(mcs(:) .EQ. val_exp ) ) THEN
       ! mcs is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    var_name= 'mcr'
    CALL ioconf_setatt_p('UNITS', 'm3/m3')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., mcr, "gather", nbp_glo, index_g)
    IF ( ALL(mcr(:) .EQ. val_exp ) ) THEN
       ! mcr is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    var_name= 'mcfc'
    CALL ioconf_setatt_p('UNITS', 'm3/m3')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., mcfc, "gather", nbp_glo, index_g)
    IF ( ALL(mcfc(:) .EQ. val_exp ) ) THEN
       ! mcfc is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    var_name= 'mcw'
    CALL ioconf_setatt_p('UNITS', 'm3/m3')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., mcw, "gather", nbp_glo, index_g)
    IF ( ALL(mcw(:) .EQ. val_exp ) ) THEN
       ! mcw is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    var_name= 'nvan'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., nvan, "gather", nbp_glo, index_g)
    IF ( ALL(nvan(:) .EQ. val_exp ) ) THEN
       ! nvan is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    var_name= 'avan'
    CALL ioconf_setatt_p('UNITS', 'm-1')
    CALL ioconf_setatt_p('LONG_NAME','')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., avan, "gather", nbp_glo, index_g)
    IF ( ALL(avan(:) .EQ. val_exp ) ) THEN
       ! avan is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    var_name= 'clay_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of clay in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., clayfraction, "gather", nbp_glo, index_g)
    IF ( ALL(clayfraction(:) .EQ. val_exp ) ) THEN
       ! clayfraction is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    var_name= 'sand_frac'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of sand in each mesh')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., sandfraction, "gather", nbp_glo, index_g)	
    IF ( ALL(sandfraction(:) .EQ. val_exp ) ) THEN
       ! sandfraction is not in restart file
       call_slowproc_soilt=.TRUE.
    END IF

    ! Calculate siltfraction not needed to be in restart file
    IF ( ALL( sandfraction(:) .EQ. val_exp) ) THEN
       siltfraction(:) = val_exp
    ELSE
       siltfraction(:) = 1. - clayfraction(:) - sandfraction(:)
    END IF
         
    !! 2.b If we want to prescribe the soil parameters values with a unique value (e.g. site experiments)
    ! AD: improve consistency with case spmipexp="unif"

    IF (impsoilt) THEN
       
       !Config Key   = SOIL_FRACTIONS
       !Config Desc  = Areal fraction of the 13 soil USDA textures; the dominant one is selected, Loam by default
       !Config Def   = 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
       !Config If    = IMPOSE_SOILT
       !Config Help  = Determines the fraction of the 13 USDA classes with same order as in constantes_soil_var
       !Config Units = [-]
       
       soilclass(:,:)=val_exp       
       CALL setvar_p (soilclass, val_exp, 'SOIL_FRACTIONS', soilclass_default)

       ! Simplify a heterogeneous grid-cell into an homogeneous one 
       ! with the dominant texture
       njsc(:) = 0
       DO ji = 1, kjpindex
          ! here we reduce to the dominant texture class
          njsc(ji) = MAXLOC(soilclass(ji,:),1)
       ENDDO
       njsc_imp = njsc(1) ! to prescribe the VG parameters consistely with the imposed texture    

       !Config Key   = CLAY_FRACTION
       !Config Desc  = Fraction of the clay fraction (0-dim mode)
       !Config Def   = 0.2 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the fraction of clay in the grid box.
       !Config Units = [-]
       clayfraction_default = clayfrac_usda(njsc_imp)
       CALL setvar_p (clayfraction, val_exp, 'CLAY_FRACTION', clayfraction_default)

       !Config Key   = SAND_FRACTION
       !Config Desc  = Fraction of the clay fraction (0-dim mode)
       !Config Def   = 0.4
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the fraction of clay in the grid box.
       !Config Units = [-] 
       sandfraction_default = sandfrac_usda(njsc_imp)
       CALL setvar_p (sandfraction, val_exp, 'SAND_FRACTION', sandfraction_default)

       ! Calculate silt fraction
       siltfraction(:) = 1. - clayfraction(:) - sandfraction(:)

       !Config Key   = nvan
       !Config Desc  = nvan parameter from Van genutchen equations
       !Config Def   = 1.56 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the nvan in the grid box.
       !Config Units = [-]
       nvan_default = nvan_usda(njsc_imp)
       CALL setvar_p (nvan, val_exp, 'NVAN_IMP', nvan_default)

       !Config Key   = avan
       !Config Desc  = avan parameter from Van genutchen equations
       !Config Def   = 0.0036 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the avan in the grid box.
       !Config Units = [-]
       avan_default = avan_usda(njsc_imp)
       CALL setvar_p (avan, val_exp, 'AVAN_IMP', avan_default)

       !Config Key   = mcr
       !Config Desc  = residual soil moisture
       !Config Def   = 0.078 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the mcr in the grid box.
       !Config Units = [-]
       mcr_default = mcr_usda(njsc_imp)
       CALL setvar_p (mcr, val_exp, 'MCR_IMP', mcr_default)

       !Config Key   = mcs
       !Config Desc  = saturation soil moisture
       !Config Def   = 0.43 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the mcs in the grid box.
       !Config Units = [-]
       mcs_default = mcs_usda(njsc_imp)
       CALL setvar_p (mcs, val_exp, 'MCS_IMP', mcs_default)

       !Config Key   = ks
       !Config Desc  = saturation conductivity
       !Config Def   = 249.6 if Loam
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the ks in the grid box.
       !Config Units = [mm/d]
       ks_default = ks_usda(njsc_imp)
       CALL setvar_p (ks, val_exp, 'KS_IMP', ks_default)

       ! By default, we calculate mcf and mcw from the above values, as in slowproc_soilt,
       ! but they can be overruled by values from run.def
                    
       mvan(:) = un - (un / nvan(:))
       ! Define matrix potential in mm for wilting point and field capacity (with sand vs clay-silt variation)
       psi_w(:) = 150000.
       DO ji=1, kjpindex
          IF ( ks(ji) .GE. 560 ) THEN ! Sandy soils (560 is equivalent of 2.75 at log scale of Ks, mm/d)
             psi_fc(ji) = 1000.
          ELSE ! Finer soils
             psi_fc(ji) = 3300. 
          ENDIF
       ENDDO
       mcfc_default(:) = mcr(:) + (( mcs(:) - mcr(:)) / (un + ( avan(:) * psi_fc(:))** nvan(:))** mvan(:))
       mcw_default(:)  = mcr(:) + (( mcs(:) - mcr(:)) / (un + ( avan(:) *  psi_w(:))** nvan(:))** mvan(:))
       
       !Config Key   = mcfc
       !Config Desc  = field capacity soil moisture
       !Config Def   = 0.1654 if caclulated from default 5 parameters above
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the mcfc in the grid box.
       !Config Units = [-]
       mcfc(:) = mcfc_default(:)
       CALL setvar_p (mcfc, val_exp, 'MCFC_IMP', mcfc_default)

       !Config Key   = mcw
       !Config Desc  = wilting point soil moisture
       !Config Def   = 0.0884 if caclulated from default 5 parameters above
       !Config If    = IMPOSE_SOIL
       !Config Help  = Determines the mcw in the grid box.
       !Config Units = [-]
       mcw(:) = mcw_default(:)
       CALL setvar_p (mcw, val_exp, 'MCW_IMP', mcw_default)
          
    ELSE ! impsoilt = F
       
    !! 2.c If some soil parameters are missing in the restart, we read a soil map
    ! Note that slowproc_soilt calculates mcf and lcw as a function of avan, nvan, etc.

       IF (call_slowproc_soilt) THEN
          
          ! If njsc is not in restart file, then initialize to default value
          ! tbdone - reading from run.def file 
          IF ( ALL(njsc(:) .EQ. undef_int )) THEN
             njsc(:) = usda_default ! 6 = Loam
          END IF
        
          CALL slowproc_soilt(njsc, ks, nvan, avan, mcr, mcs, mcfc, mcw, kjpindex, &
               lalo, neighbours, resolution, contfrac, soilclass, &
               clayfraction, sandfraction, siltfraction)
          
          call_slowproc_soilt=.FALSE.
       
       ENDIF

    ENDIF



















        frac_imperv(:) = zero
        coeff_imperv(:) = un

    IF ( do_imperviousness ) THEN
        var_name = 'frac_imperv'
        CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., frac_imperv, "gather", nbp_glo, index_g)
        !IF ( ALL( frac_imperv(:) .EQ. val_exp ) ) THEN
           CALL slowproc_imperviousness(kjpindex, lalo, neighbours,  resolution, contfrac)
           CALL xios_orchidee_send_field("frac_imperv",frac_imperv)
        !ENDIF
        
        coeff_imperv(:) = -0.999999 * frac_imperv(:) + un
        
        CALL xios_orchidee_send_field("coeff_imperv",coeff_imperv)
        ks(:) = ks(:)*coeff_imperv(:)
     ENDIF
        
























    ! XIOS export of Ks before changing the vertical profile
    CALL xios_orchidee_send_field("ksref",ks) ! mm/d (for CMIP6, once)

    !! 3. Read the infiltration related variables
    ! This variable helps reducing surface runuff in flat areas
    
    !! 3.a Looking first in the restart files
    
    var_name= 'reinf_slope'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Slope coef for reinfiltration')
    CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., reinf_slope, "gather", nbp_glo, index_g)

    !! 3.b We can also read from a map or prescribe, depending on IMPOSE_SLOPE

    IF (impslope) THEN

       !! Impose a constant value from the run.def (not anymore controlled by impveg)
       
       !Config Key   = REINF_SLOPE
       !Config Desc  = Fraction of reinfiltrated surface runoff 
       !Config Def   = 0.1
       !Config If    = IMPOSE_SLOPE
       !Config Help  = Determines the reinfiltration ratio in the grid box due to flat areas
       !Config Units = [-]
       slope_default=0.1
       CALL setvar_p (reinf_slope, val_exp, 'SLOPE', slope_default)

    ELSE

       !Config Key   = GET_SLOPE
       !Config Desc  = Read slopes from file and do the interpolation
       !Config Def   = n
       !Config If    =
       !Config Help  = Needed for reading the slope file and doing the interpolation. This will be
       !               used by the re-infiltration parametrization
       !Config Units = [FLAG]
       get_slope = .FALSE.
       CALL getin_p('GET_SLOPE',get_slope)

       !! Else, if not found from restart or GET_SLOPE = T, we read from a map 
       IF ( MINVAL(reinf_slope) .EQ. MAXVAL(reinf_slope) .AND. MAXVAL(reinf_slope) .EQ. val_exp .OR. get_slope) THEN
          IF (printlev_loc>=4) WRITE (numout,*) 'reinf_slope was not in restart file. Now call slowproc_slope'
          
          CALL slowproc_slope(kjpindex, lalo, neighbours, resolution, contfrac, reinf_slope)
          IF (printlev_loc>=4) WRITE (numout,*) 'After slowproc_slope'
          
       ENDIF

    ENDIF
    
    !! 4. Read the vegetation related variables

    !! 4.a General parameters 
    
    !Config Key   = SECHIBA_QSINT 
    !Config Desc  = Interception reservoir coefficient
    !Config If    = OK_SECHIBA 
    !Config Def   = 0.02
    !Config Help  = Transforms leaf area index into size of interception reservoir
    !Config         for slowproc_derivvar or stomate
    !Config Units = [m]
    CALL getin_p('SECHIBA_QSINT', qsintcst)
    IF (printlev >= 2) WRITE(numout, *)' SECHIBA_QSINT, qsintcst = ', qsintcst

    !Config Key   = SECHIBA_ZCANOP
    !Config Desc  = Soil level used for canopy development (if STOMATE disactivated)
    !Config If    = OK_SECHIBA and .NOT. OK_STOMATE  
    !Config Def   = 0.5
    !Config Help  = The temperature at this soil depth is used to determine the LAI when
    !Config         STOMATE is not activated.
    !Config Units = [m]
    zcanop = 0.5_r_std
    CALL setvar_p (zcanop, val_exp, 'SECHIBA_ZCANOP', 0.5_r_std)

    ! depth at center of the levels
    zsoil(1) = diaglev(1) / 2.
    DO l = 2, nslm
       zsoil(l) = ( diaglev(l) + diaglev(l-1) ) / 2.
    ENDDO

    ! index of this level
    vtmp = MINLOC ( ABS ( zcanop - zsoil(:) ) )
    lcanop = vtmp(1)

    !! 4.b From restart first

    found_restart=.TRUE.
    var_name= 'veget'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Vegetation fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., veget, "gather", nbp_glo, index_g)
    IF ( ALL( veget(:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    var_name= 'veget_max'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Maximum vegetation fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., veget_max, "gather", nbp_glo, index_g)
    IF ( ALL( veget_max(:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    IF (do_wood_harvest) THEN
       var_name= 'woodharvest'
       CALL ioconf_setatt_p('UNITS', 'gC m-2 yr-1')
       CALL ioconf_setatt_p('LONG_NAME','Harvest wood biomass')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., woodharvest, "gather", nbp_glo, index_g)
       IF ( ALL( woodharvest(:) .EQ. val_exp ) ) woodharvest(:)=zero
    END IF

    var_name= 'frac_nobio'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Special soil type fraction')
    CALL restget_p (rest_id, var_name, nbp_glo, nnobio, 1, kjit, .TRUE., frac_nobio, "gather", nbp_glo, index_g)
    IF ( ALL( frac_nobio(:,:) .EQ. val_exp ) ) found_restart=.FALSE.

    IF (.NOT. impveg) THEN
       var_name= 'veget_year'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Last year get in Land Use file.')
       IF (is_root_prc) THEN
          CALL restget (rest_id, var_name, 1       , 1  , 1, kjit, .TRUE., tmp_veget_year)
          IF (veget_reinit) THEN
             ! Do not take the value read from restart file
             veget_year=veget_year_orig
          ELSE IF (tmp_veget_year(1) == val_exp) THEN
             ! veget_year was not found in restart file
             veget_year=veget_year_orig
          ELSE
             ! veget_year was found in restart file, transform to integer
             veget_year=INT(tmp_veget_year(1))
          ENDIF
       ENDIF
       CALL bcast(veget_year)

       !
       !Config Key   = VEGET_UPDATE
       !Config Desc  = Update vegetation frequency: 0Y or 1Y
       !Config If    = 
       !Config Def   = 0Y
       !Config Help  = The veget datas will be update each this time step. Must be 0Y if IMPOSE_VEG=y.
       !Config Units = [years]
       !
       veget_update=0
       WRITE(veget_str,'(a)') '0Y'
       CALL getin_p('VEGET_UPDATE', veget_str)
       l=INDEX(TRIM(veget_str),'Y')
       READ(veget_str(1:(l-1)),"(I2.2)") veget_update
       IF (printlev_loc >= 2) WRITE(numout,*) "Update frequency for land use in years :",veget_update

       ! Coherence test
       IF (veget_update > 0 .AND. ok_dgvm .AND. .NOT. agriculture) THEN
          CALL ipslerr_p(3,'slowproc_init',&
               'The combination DGVM=TRUE, AGRICULTURE=FALSE and VEGET_UPDATE>0 is not possible', &
               'Set VEGET_UPDATE=0Y in run.def','')
       END IF
       
    ELSE
       ! impveg=TRUE: there can not be any land use change, veget_update must be =0
       ! Read VEGET_UPDATE from run.def and exit if it is different from 0Y
       veget_update=0
       WRITE(veget_str,'(a)') '0Y'
       CALL getin_p('VEGET_UPDATE', veget_str)
       l=INDEX(TRIM(veget_str),'Y')
       READ(veget_str(1:(l-1)),"(I2.2)") veget_update
       IF (veget_update /= 0) THEN
          WRITE(numout,*) 'veget_update=',veget_update,' is not coeherent with impveg=',impveg
          CALL ipslerr_p(3,'slowproc_init','Incoherent values between impveg and veget_update', &
               'VEGET_UPDATE must be equal to 0Y if IMPOSE_VEG=y (impveg=true)','')
       END IF

    ENDIF

    var_name= 'lai'
    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Leaf area index')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., lai, "gather", nbp_glo, index_g)

    ! The height of the vegetation could in principle be recalculated at the beginning of the run.
    ! However, this is very tedious, as many special cases have to be taken into account. This variable
    ! is therefore saved in the restart file.
    var_name= 'height'
    CALL ioconf_setatt_p('UNITS', 'm')
    CALL ioconf_setatt_p('LONG_NAME','Height of vegetation')
    CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE., height, "gather", nbp_glo, index_g)
 
    IF (read_lai)THEN
       var_name= 'laimap'
       CALL ioconf_setatt_p('UNITS', '-')
       CALL ioconf_setatt_p('LONG_NAME','Leaf area index read')
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 12, kjit, .TRUE., laimap)
    ENDIF

    CALL ioconf_setatt_p('UNITS', '-')
    CALL ioconf_setatt_p('LONG_NAME','Fraction of leaves in leaf age class ')
    DO jf = 1, nleafages
       ! variable name is somewhat complicated as ioipsl does not allow 3d variables for the moment...
       WRITE(laistring,'(i4)') jf
       laistring=ADJUSTL(laistring)
       var_name='frac_age_'//laistring(1:LEN_TRIM(laistring))
       CALL restget_p (rest_id, var_name, nbp_glo, nvm, 1, kjit, .TRUE.,frac_age(:,:,jf), "gather", nbp_glo, index_g)
    ENDDO

    !! 4.c Initialization of variables not found in restart file

    IF ( impveg ) THEN

       !! 4.1.a Case impveg=true: Initialization of variables by reading run.def
       !!       The routine setvar_p will only initialize the variable if it was not found in restart file. 
       !!       We are on a point and thus we can read the information from the run.def
       
       !Config Key   = SECHIBA_VEGMAX
       !Config Desc  = Maximum vegetation distribution within the mesh (0-dim mode)
       !Config If    = IMPOSE_VEG
       !Config Def   = 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0
       !Config Help  = The fraction of vegetation is read from the restart file. If
       !Config         it is not found there we will use the values provided here.
       !Config Units = [-]
       CALL setvar_p (veget_max, val_exp, 'SECHIBA_VEGMAX', veget_ori_fixed_test_1)

       !Config Key   = SECHIBA_FRAC_NOBIO
       !Config Desc  = Fraction of other surface types within the mesh (0-dim mode)
       !Config If    = IMPOSE_VEG
       !Config Def   = 0.0
       !Config Help  = The fraction of ice, lakes, etc. is read from the restart file. If
       !Config         it is not found there we will use the values provided here.
       !Config         For the moment, there is only ice.
       !Config Units = [-]
       frac_nobio1 = frac_nobio(1,1)
       CALL setvar_p (frac_nobio1, val_exp, 'SECHIBA_FRAC_NOBIO', frac_nobio_fixed_test_1)
       frac_nobio(:,:) = frac_nobio1

       IF (.NOT. found_restart) THEN
          ! Call slowproc_veget to correct veget_max and to calculate veget and soiltiles
          CALL slowproc_veget (kjpindex, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile, fraclut, nwdFraclut)
       END IF
       
       !Config Key   = SECHIBA_LAI
       !Config Desc  = LAI for all vegetation types (0-dim mode)
       !Config Def   = 0., 8., 8., 4., 4.5, 4.5, 4., 4.5, 4., 2., 2., 2., 2.
       !Config If    = IMPOSE_VEG
       !Config Help  = The maximum LAI used in the 0dim mode. The values should be found
       !Config         in the restart file. The new values of LAI will be computed anyway
       !Config         at the end of the current day. The need for this variable is caused
       !Config         by the fact that the model may stop during a day and thus we have not
       !Config         yet been through the routines which compute the new surface conditions.
       !Config Units = [-]
       CALL setvar_p (lai, val_exp, 'SECHIBA_LAI', llaimax)

       !Config Key   = SLOWPROC_HEIGHT
       !Config Desc  = Height for all vegetation types 
       !Config Def   = 0., 30., 30., 20., 20., 20., 15., 15., 15., .5, .6, 1.0, 1.0
       !Config If    = OK_SECHIBA
       !Config Help  = The height used in the 0dim mode. The values should be found
       !Config         in the restart file. The new values of height will be computed anyway
       !Config         at the end of the current day. The need for this variable is caused
       !Config         by the fact that the model may stop during a day and thus we have not
       !Config         yet been through the routines which compute the new surface conditions.
       !Config Units = [m]
       CALL setvar_p (height, val_exp, 'SLOWPROC_HEIGHT', height_presc)


    ELSE IF ( .NOT. found_restart .OR. vegetmap_reset ) THEN
 
       !! 4.1.b Case impveg=false and no restart files: Initialization by reading vegetation map
       
       ! Initialize veget_max and frac_nobio
       ! Case without restart file
       IF (printlev_loc>=3) WRITE(numout,*) 'Before call slowproc_readvegetmax in initialization phase without restart files'
       IF (printlev_loc>=3) WRITE(numout,*) 'veget_year=', veget_year
          
       ! Call the routine to read the vegetation from file (output is veget_max_new)
       CALL slowproc_readvegetmax(kjpindex, lalo, neighbours, resolution, contfrac, &
            veget_max, veget_max_new, frac_nobio_new, veget_year, .TRUE.)
       IF (printlev_loc>=4) WRITE (numout,*) 'After slowproc_readvegetmax in initialization phase'
       
       ! Update vegetation with values read from the file
       veget_max           = veget_max_new
       frac_nobio          = frac_nobio_new         
       
       IF (do_wood_harvest) THEN
          ! Read the new the wood harvest map from file. Output is wood harvest
          CALL slowproc_woodharvest(kjpindex, lalo, neighbours, resolution, contfrac, woodharvest)
       ENDIF
              
       !! Reset totaly or partialy veget_max if using DGVM
       IF ( ok_dgvm  ) THEN
          ! If we are dealing with dynamic vegetation then all natural PFTs should be set to veget_max = 0
          ! In case no agriculture is desired, agriculture PFTS should be set to 0 as well
          IF (agriculture) THEN
             DO jv = 2, nvm
                IF (natural(jv)) THEN 
                   veget_max(:,jv)=zero
                ENDIF
             ENDDO
             
             ! Calculate the fraction of crop for each point.
             ! Sum only on the indexes corresponding to the non_natural pfts
             frac_crop_tot(:) = zero
             DO jv = 2, nvm
                IF(.NOT. natural(jv)) THEN
                   DO ji = 1, kjpindex
                      frac_crop_tot(ji) = frac_crop_tot(ji) + veget_max(ji,jv)
                   ENDDO
                ENDIF
             END DO
            
             ! Calculate the fraction of bare soil
             DO ji = 1, kjpindex
                veget_max(ji,1) = un - frac_crop_tot(ji) - SUM(frac_nobio(ji,:))   
             ENDDO
          ELSE
             veget_max(:,:) = zero
             DO ji = 1, kjpindex
                veget_max(ji,1) = un  - SUM(frac_nobio(ji,:))
             ENDDO
          END IF
       END IF   ! end ok_dgvm
       

       ! Call slowproc_veget to correct veget_max and to calculate veget and soiltiles
       CALL slowproc_veget (kjpindex, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile, fraclut, nwdFraclut)
       
    END IF ! end impveg

    !! 4.d Continue initializing variables not found in restart file. Case for both impveg=true and false.

    ! Initialize laimap for the case read_lai if not found in restart file
    IF (read_lai) THEN
       IF ( ALL( laimap(:,:,:) .EQ. val_exp) ) THEN
          ! Interpolation of LAI
          CALL slowproc_interlai (kjpindex, lalo, resolution,  neighbours, contfrac, laimap)
       ENDIF
    ENDIF
    
    ! Initialize lai if not found in restart file and not already initialized using impveg
    IF ( MINVAL(lai) .EQ. MAXVAL(lai) .AND. MAXVAL(lai) .EQ. val_exp) THEN
       IF (read_lai) THEN
          stempdiag2_bid(1:kjpindex,1:nslm) = stempdiag_bid
          CALL slowproc_lai (kjpindex, lcanop, stempdiag2_bid, &
               lalo,resolution,lai,laimap)
       ELSE
          ! If we start from scratch, we set lai to zero for consistency with stomate
          lai(:,:) = zero
       ENDIF
       
       frac_age(:,:,1) = un
       frac_age(:,:,2) = zero
       frac_age(:,:,3) = zero
       frac_age(:,:,4) = zero
    ENDIF
    
    ! Initialize heigth if not found in restart file and not already initialized using impveg
    IF ( MINVAL(height) .EQ. MAXVAL(height) .AND. MAXVAL(height) .EQ. val_exp) THEN
       ! Impose height
       DO jv = 1, nvm
          height(:,jv) = height_presc(jv)
       ENDDO
       ! If urban config with wudapt info, update values for PFT16
       ! IF (DO_HEIGHT_BUILDING = y) THEN
       CALL slowproc_height_buidings(kjpindex, lalo, neighbours,  resolution, contfrac)
       height(:,16) = height_buidings(:)
       ! ENDIF

    ENDIF
    

    !! 4.3 Dynamic irrigation map
    !  If do_irrigation, it will look to the dynamical irrig. map in restart
    !  If not dynamic irrig. map, it will be set to zero.
    !  If not found in restart, it will try to interpolate the map

    irrigated_new(:) = zero !
    IF ( do_irrigation ) THEN
       ! It will look into restart file
       var_name = 'irrigmap_dyn'
       CALL ioconf_setatt_p('UNITS', 'm2')
       CALL ioconf_setatt_p('LONG_NAME','Dynamical area equipped for irrigation')
       CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., irrigated_new, "gather", nbp_glo, index_g)
       
       ! Now, if not found in the restart, read and interpolate from file
       IF ( ALL( irrigated_new(:) .EQ. val_exp ) ) THEN
          CALL slowproc_readirrigmap_dyn(kjpindex, lalo, neighbours,  resolution, contfrac,         &
               irrigated_new)
       ENDIF
       ! irrigated_next from sechiba (int output) = irrigated_new  in slowproc_initialize.
    ENDIF

    ! If new priorization scheme is used, it also seek into restart, or interpolate

    fraction_aeirrig_sw(:) = un
    IF ( do_irrigation .AND. select_source_irrig) THEN
      ! It will look into restart file
      var_name = 'fraction_aeirrig_sw'
      CALL ioconf_setatt_p('UNITS', '%')
      CALL ioconf_setatt_p('LONG_NAME','Fraction of area equipped for irrigation with surface water')
      CALL restget_p (rest_id, var_name, nbp_glo, 1, 1, kjit, .TRUE., fraction_aeirrig_sw, "gather", nbp_glo, index_g)

      ! Now, if not found in the restart, read and interpolate from file
      IF ( ALL( fraction_aeirrig_sw(:) .EQ. val_exp ) ) THEN
         CALL slowproc_read_aeisw_map(kjpindex, lalo, neighbours,  resolution, contfrac,         &
              fraction_aeirrig_sw)
      ENDIF
      ! irrigated_next from sechiba (int output) = irrigated_new in slowproc_initialize.
   ENDIF


    !! 5. Some calculations always done, with and without restart files
       
    ! The variables veget, veget_max and frac_nobio were all read from restart file or initialized above.
    ! Calculate now totfrac_nobio and soiltiles using these variables.
    
    ! Calculate totfrac_nobio
    totfrac_nobio(:) = zero
    DO jv = 1, nnobio
       totfrac_nobio(:) = totfrac_nobio(:) + frac_nobio(:,jv)
    ENDDO
    
    ! Calculate soiltile. This variable do not need to be in the restart file.
    ! The sum of all soiltiles makes one, and corresponds to the bio fraction
    ! of the grid cell (called vegtot in hydrol)
    soiltile(:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji = 1, kjpindex
          soiltile(ji,jst) = soiltile(ji,jst) + veget_max(ji,jv)
       ENDDO
    ENDDO
    DO ji = 1, kjpindex 
       IF (totfrac_nobio(ji) .LT. (1-min_sechiba)) THEN
          soiltile(ji,:)=soiltile(ji,:)/(1-totfrac_nobio(ji))
       ENDIF
    ENDDO
    
    ! Always calculate tot_bare_soil
    ! Fraction of bare soil in the mesh (bio+nobio) 
    tot_bare_soil(:) = veget_max(:,1)
    DO jv = 2, nvm
       DO ji =1, kjpindex
          tot_bare_soil(ji) = tot_bare_soil(ji) + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
    END DO
    

    !! Calculate fraction of landuse tiles to be used only for diagnostic variables
    fraclut(:,:)=0
    nwdFraclut(:,id_psl)=0
    nwdFraclut(:,id_crp)=1.
    nwdFraclut(:,id_urb)=xios_default_val
    nwdFraclut(:,id_pst)=xios_default_val
    DO jv=1,nvm
       IF (natural(jv)) THEN
          fraclut(:,id_psl) = fraclut(:,id_psl) + veget_max(:,jv)
          IF(.NOT. is_tree(jv)) THEN
             nwdFraclut(:,id_psl) = nwdFraclut(:,id_psl) + veget_max(:,jv) 
          ENDIF
       ELSE
          fraclut(:,id_crp) = fraclut(:,id_crp) + veget_max(:,jv)
       ENDIF
    END DO
    
    WHERE (fraclut(:,id_psl) > min_sechiba)
       nwdFraclut(:,id_psl) = nwdFraclut(:,id_psl)/fraclut(:,id_psl)
    ELSEWHERE
       nwdFraclut(:,id_psl) = xios_default_val
    END WHERE    


    IF (printlev_loc>=3) WRITE (numout,*) ' slowproc_init done '
    
  END SUBROUTINE slowproc_init

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_clear
!!
!>\BRIEF          Clear all variables related to slowproc and stomate modules  
!!
!_ ================================================================================================================================

  SUBROUTINE slowproc_clear 

  ! 1 clear all the variables defined as common for the routines in slowproc 

    IF (ALLOCATED (clayfraction)) DEALLOCATE (clayfraction)
    IF (ALLOCATED (sandfraction)) DEALLOCATE (sandfraction)
    IF (ALLOCATED (siltfraction)) DEALLOCATE (siltfraction)
    IF (ALLOCATED (laimap)) DEALLOCATE (laimap)
    IF (ALLOCATED (veget_max_new)) DEALLOCATE (veget_max_new)
    IF (ALLOCATED (irrigated_new)) DEALLOCATE (irrigated_new)
    IF (ALLOCATED (woodharvest)) DEALLOCATE (woodharvest)
    IF (ALLOCATED (frac_nobio_new)) DEALLOCATE (frac_nobio_new)
    IF (ALLOCATED (frac_imperv)) DEALLOCATE (frac_imperv)
    IF (ALLOCATED (coeff_imperv)) DEALLOCATE (coeff_imperv)
 ! 2. Clear all the variables in stomate 

    CALL stomate_clear 
    !
  END SUBROUTINE slowproc_clear

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_derivvar
!!
!>\BRIEF         Initializes variables related to the
!! parameters to be assimilated, the maximum water on vegetation, the vegetation height, 
!! and the fraction of soil covered by dead leaves and the vegetation height 
!!
!! DESCRIPTION  : (definitions, functional, design, flags):
!! (1) Initialization of the variables relevant for the assimilation parameters  
!! (2) Intialization of the fraction of soil covered by dead leaves
!! (3) Initialization of the Vegetation height per PFT
!! (3) Initialization the maximum water on vegetation for interception with a particular treatement of the PFT no.1
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::qsintmax, ::deadleaf_cover, ::assim_param, ::height  
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_derivvar (kjpindex, veget, lai, &
       qsintmax, deadleaf_cover, assim_param, height, temp_growth)

    !! INTERFACE DESCRIPTION

    !! 0.1 Input scalar and fields 
    INTEGER(i_std), INTENT (in)                                :: kjpindex       !! Domain size - terrestrial pixels only
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)          :: veget          !! Fraction of pixel covered by PFT. Fraction accounts for none-biological land covers (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (in)          :: lai            !! PFT leaf area index (m^{2} m^{-2})

    !! 0.2. Output scalar and fields 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)          :: qsintmax       !! Maximum water on vegetation for interception(mm)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)              :: deadleaf_cover !! fraction of soil covered by dead leaves (unitless)
    REAL(r_std), DIMENSION (kjpindex,nvm,npco2), INTENT (out)   :: assim_param    !! min+max+opt temperatures & vmax for photosynthesis (K, \mumol m^{-2} s^{-1})
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)          :: height         !! height of the vegetation or surface in general ??? (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)              :: temp_growth    !! growth temperature (°C)  
    !
    !! 0.3 Local declaration
    INTEGER(i_std)                                              :: jv             !! Local indices
!_ ================================================================================================================================

    !
    ! 1. Initialize (why here ??) the variables revelant for the assimilation parameters
    !
    DO jv = 1, nvm
       assim_param(:,jv,ivcmax) = vcmax_fix(jv)
    ENDDO

    !
    ! 2. Intialize the fraction of soil covered by dead leaves 
    !
    deadleaf_cover(:) = zero

    !
    ! 3. Initialize the Vegetation height per PFT
    !
    DO jv = 1, nvm
       height(:,jv) = height_presc(jv)
    ENDDO
    !
    ! 4. Initialize the maximum water on vegetation for interception
    !
    qsintmax(:,:) = qsintcst * veget(:,:) * lai(:,:)

    ! Added by Nathalie - July 2006
    !  Initialize the case of the PFT no.1 to zero 
    qsintmax(:,1) = zero

    temp_growth(:)=25.

  END SUBROUTINE slowproc_derivvar


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_mean
!!
!>\BRIEF          Accumulates field_in over a period of dt_tot.
!! Has to be called at every time step (dt). 
!! Mean value is calculated if ldmean=.TRUE.
!! field_mean must be initialized outside of this routine! 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) AcumAcuumlm 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::field_main
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_mean (npts, n_dim2, dt_tot, dt, ldmean, field_in, field_mean)

    !
    !! 0 declarations

    !! 0.1 input scalar and variables 
    INTEGER(i_std), INTENT(in)                           :: npts     !! Domain size- terrestrial pixels only 
    INTEGER(i_std), INTENT(in)                           :: n_dim2   !! Number of PFTs 
    REAL(r_std), INTENT(in)                              :: dt_tot   !! Time step of stomate (in days). The period over which the accumulation or the mean is computed 
    REAL(r_std), INTENT(in)                              :: dt       !! Time step in days 
    LOGICAL, INTENT(in)                                  :: ldmean   !! Flag to calculate the mean after the accumulation ???
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(in)      :: field_in !! Daily field 

    !! 0.3 Modified field; The computed sum or mean field over dt_tot time period depending on the flag ldmean 
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(inout)   :: field_mean !! Accumulated field at dt_tot time period or mean field over dt_tot 
 

!_ ================================================================================================================================

    !
    ! 1. Accumulation the field over dt_tot period 
    !
    field_mean(:,:) = field_mean(:,:) + field_in(:,:) * dt

    !
    ! 2. If the flag ldmean set, the mean field is computed over dt_tot period  
    !
    IF (ldmean) THEN
       field_mean(:,:) = field_mean(:,:) / dt_tot
    ENDIF

  END SUBROUTINE slowproc_mean


  
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_long
!!
!>\BRIEF        Calculates a temporally smoothed field (field_long) from
!! instantaneous input fields.Time constant tau determines the strength of the smoothing.
!! For tau -> infinity??, field_long becomes the true mean value of field_inst
!! (but  the spinup becomes infinietly long, too).
!! field_long must be initialized outside of this routine! 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) Testing the time coherence betwen the time step dt and the time tau over which
!! the rescaled of the mean is performed   
!!  (2) Computing the rescaled mean over tau period 
!! MAIN OUTPUT VARIABLE(S): field_long  
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::field_long
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_long (npts, n_dim2, dt, tau, field_inst, field_long)

    !
    ! 0 declarations
    !

    ! 0.1 input scalar and fields 

    INTEGER(i_std), INTENT(in)                                 :: npts        !! Domain size- terrestrial pixels only
    INTEGER(i_std), INTENT(in)                                 :: n_dim2      !! Second dimension of the fields, which represents the number of PFTs
    REAL(r_std), INTENT(in)                                    :: dt          !! Time step in days   
    REAL(r_std), INTENT(in)                                    :: tau         !! Integration time constant (has to have same unit as dt!)  
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(in)            :: field_inst  !! Instantaneous field 


    ! 0.2 modified field

    ! Long-term field
    REAL(r_std), DIMENSION(npts,n_dim2), INTENT(inout)         :: field_long  !! Mean value of the instantaneous field rescaled at tau time period 

!_ ================================================================================================================================

    !
    ! 1 test coherence of the time 

    IF ( ( tau .LT. dt ) .OR. ( dt .LE. zero ) .OR. ( tau .LE. zero ) ) THEN
       WRITE(numout,*) 'slowproc_long: Problem with time steps'
       WRITE(numout,*) 'dt=',dt
       WRITE(numout,*) 'tau=',tau
    ENDIF

    !
    ! 2 integration of the field over tau 

    field_long(:,:) = ( field_inst(:,:)*dt + field_long(:,:)*(tau-dt) ) / tau

  END SUBROUTINE slowproc_long


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_veget_max_limit
!!
!>\BRIEF        Set small fractions of veget_max to zero and normalize to keep the sum equal 1
!!
!! DESCRIPTION  : Set small fractions of veget_max to zero and normalize to keep the sum equal 1
!!
!! RECENT CHANGE(S): The subroutine was previously a part of slowproc_veget,
!!    but was separated to be called also from slowproc_readvegetmax in order
!!    to have limited/normalized vegetation fractions right after its reading
!!    from the file (added by V.Bastrikov, 15/06/2019)
!!
!! MAIN OUTPUT VARIABLE(S): :: frac_nobio, veget_max
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_veget_max_limit (kjpindex, frac_nobio, veget_max)
    !
    ! 0. Declarations
    !
    ! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                             :: kjpindex    !! Domain size - terrestrial pixels only

    ! 0.2 Modified variables 
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(inout) :: frac_nobio  !! Fraction of the mesh which is covered by ice, lakes, ...
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)    :: veget_max   !! Maximum fraction of vegetation type including none biological fraction (unitless)

    ! 0.4 Local scalar and varaiables 
    INTEGER(i_std)                                         :: ji, jv      !! indices 
    REAL(r_std)                                            :: SUMveg      !! Total vegetation summed across PFTs

!_ ================================================================================================================================
    IF (printlev_loc >= 3) WRITE(numout,*) 'Entering slowproc_veget_max_limit'

    !! Set to zero fractions of frac_nobio and veget_max smaller than min_vegfrac
    DO ji = 1, kjpindex
       IF ( SUM(frac_nobio(ji,:)) .LT. min_vegfrac ) THEN
          frac_nobio(ji,:) = zero
       ENDIF
    
       IF (.NOT. ok_dgvm) THEN
          DO jv = 1, nvm
             IF ( veget_max(ji,jv) .LT. min_vegfrac ) THEN
                veget_max(ji,jv) = zero
             ENDIF
          ENDDO
       END IF
 
       !! Normalize to keep the sum equal 1.
       SUMveg = SUM(frac_nobio(ji,:))+SUM(veget_max(ji,:))
       frac_nobio(ji,:) = frac_nobio(ji,:)/SUMveg
       veget_max(ji,:) = veget_max(ji,:)/SUMveg
    ENDDO

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_veget_max_limit ended'

  END SUBROUTINE slowproc_veget_max_limit


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_veget
!!
!>\BRIEF        Set small fractions to zero and normalize to keep the sum equal 1. Calucate veget and soiltile.
!!
!! DESCRIPTION  : Set small fractions to zero and normalize to keep the sum equal 1. Calucate veget and soiltile.
!! (1) Set veget_max and frac_nobio for fraction smaller than min_vegfrac.
!! (2) Calculate veget
!! (3) Calculate totfrac_nobio
!! (4) Calculate soiltile
!! (5) Calculate fraclut
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: frac_nobio, totfrac_nobio, veget_max, veget, soiltile, fraclut
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_veget (kjpindex, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile, fraclut, nwdFraclut)
    !
    ! 0. Declarations
    !
    ! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                             :: kjpindex    !! Domain size - terrestrial pixels only
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(in)       :: lai         !! PFT leaf area index (m^{2} m^{-2})

    ! 0.2 Modified variables 
    REAL(r_std), DIMENSION(kjpindex,nnobio), INTENT(inout) :: frac_nobio  !! Fraction of the mesh which is covered by ice, lakes, ...
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(inout)    :: veget_max   !! Maximum fraction of vegetation type including none biological fraction (unitless)

    ! 0.3 Output variables 
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)      :: veget       !! Fraction of pixel covered by PFT. Fraction accounts for none-biological land covers (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)         :: totfrac_nobio
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)    :: soiltile     !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)    :: fraclut      !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)    :: nwdFraclut   !! Fraction of non-woody vegetation in each landuse tile (0-1, unitless)

    ! 0.4 Local scalar and varaiables 
    INTEGER(i_std)                                         :: ji, jv, jst !! indices 

!_ ================================================================================================================================
    IF (printlev_loc > 8) WRITE(numout,*) 'Entering slowproc_veget'

    !! 1. Set to zero fractions of frac_nobio and veget_max smaller than min_vegfrac
    !!    Normalize to have the sum equal 1.
    CALL slowproc_veget_max_limit(kjpindex, frac_nobio, veget_max)

    !! 2. Calculate veget
    !!    If lai of a vegetation type (jv > 1) is small, increase soil part
    !!    stomate-like calculation
    DO ji = 1, kjpindex
       veget(ji,1)=veget_max(ji,1)
       DO jv = 2, nvm
          veget(ji,jv) = veget_max(ji,jv) * ( un - exp( - lai(ji,jv) * ext_coeff_vegetfrac(jv) ) )
       ENDDO
    ENDDO


    !! 3. Calculate totfrac_nobio
    totfrac_nobio(:) = zero
    DO jv = 1, nnobio
       totfrac_nobio(:) = totfrac_nobio(:) + frac_nobio(:,jv)
    ENDDO
    

    !! 4. Calculate soiltiles
    !! Soiltiles are only used in hydrol, but we fix them in here because some time it might depend
    !! on a changing vegetation (but then some adaptation should be made to hydrol) and be also used
    !! in the other modules to perform separated energy balances
    ! The sum of all soiltiles makes one, and corresponds to the bio fraction
    ! of the grid cell (called vegtot in hydrol)   
    soiltile(:,:) = zero
    DO jv = 1, nvm
       jst = pref_soil_veg(jv)
       DO ji = 1, kjpindex
          soiltile(ji,jst) = soiltile(ji,jst) + veget_max(ji,jv)
       ENDDO
    ENDDO
    DO ji = 1, kjpindex 
       IF (totfrac_nobio(ji) .LT. (1-min_sechiba)) THEN
          soiltile(ji,:)=soiltile(ji,:)/(1.-totfrac_nobio(ji))
       ENDIF
    ENDDO   

    !! 5. Calculate fraction of landuse tiles to be used only for diagnostic variables
    fraclut(:,:)=0
    nwdFraclut(:,id_psl)=0
    nwdFraclut(:,id_crp)=1.
    nwdFraclut(:,id_urb)=xios_default_val
    nwdFraclut(:,id_pst)=xios_default_val
    DO jv=1,nvm
       IF (natural(jv)) THEN
          fraclut(:,id_psl) = fraclut(:,id_psl) + veget_max(:,jv)
          IF(.NOT. is_tree(jv)) THEN
             nwdFraclut(:,id_psl) = nwdFraclut(:,id_psl) + veget_max(:,jv) 
          ENDIF
       ELSE
          fraclut(:,id_crp) = fraclut(:,id_crp) + veget_max(:,jv)
       ENDIF
    END DO
    
    WHERE (fraclut(:,id_psl) > min_sechiba)
       nwdFraclut(:,id_psl) = nwdFraclut(:,id_psl)/fraclut(:,id_psl)
    ELSEWHERE
       nwdFraclut(:,id_psl) = xios_default_val
    END WHERE    

  END SUBROUTINE slowproc_veget
 
 
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_lai
!!
!>\BRIEF        Do the interpolation of lai for the PFTs in case the laimap is not read   
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!! (1) Interplation by using the mean value of laimin and laimax for the PFTs    
!! (2) Interpolation between laimax and laimin values by using the temporal
!!  variations 
!! (3) If problem occurs during the interpolation, the routine stops 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::lai
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_lai (kjpindex,lcanop,stempdiag,lalo,resolution,lai,laimap)
    !
    ! 0. Declarations
    !
    !! 0.1 Input variables 
    INTEGER(i_std), INTENT(in)                          :: kjpindex   !! Domain size - terrestrial pixels only
    INTEGER(i_std), INTENT(in)                          :: lcanop     !! soil level used for LAI
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)  :: stempdiag  !! Soil temperature (K) ???
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo       !! Geogr. coordinates (latitude,longitude) (degrees)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)     :: resolution !! Size in x an y of the grid (m) - surface area of the gridbox
    REAL(r_std), DIMENSION(:,:,:), INTENT(in)           :: laimap     !! map of lai read 

    !! 0.2 Output
    REAL(r_std), DIMENSION(kjpindex,nvm), INTENT(out)   :: lai        !! PFT leaf area index (m^{2} m^{-2})LAI

    !! 0.4 Local
    INTEGER(i_std)                                      :: ji,jv      !! Local indices 
!_ ================================================================================================================================

    !
    IF  ( .NOT. read_lai ) THEN
    
       lai(: ,1) = zero
       ! On boucle sur 2,nvm au lieu de 1,nvm
       DO jv = 2,nvm
          SELECT CASE (type_of_lai(jv))
             
          CASE ("mean ")
             !
             ! 1. do the interpolation between laimax and laimin
             !
             lai(:,jv) = undemi * (llaimax(jv) + llaimin(jv))
             !
          CASE ("inter")
             !
             ! 2. do the interpolation between laimax and laimin
             !
             DO ji = 1,kjpindex
                lai(ji,jv) = llaimin(jv) + tempfunc(stempdiag(ji,lcanop)) * (llaimax(jv) - llaimin(jv))
             ENDDO
             !
          CASE default
             !
             ! 3. Problem
             !
             WRITE (numout,*) 'This kind of lai choice is not possible. '// &
                  ' We stop with type_of_lai ',jv,' = ', type_of_lai(jv) 
             CALL ipslerr_p(3,'slowproc_lai','Bad value for type_of_lai','read_lai=false','')
          END SELECT
          
       ENDDO
       !
    ELSE
       lai(: ,1) = zero
       ! On boucle sur 2,nvm au lieu de 1,nvm
       DO jv = 2,nvm

          SELECT CASE (type_of_lai(jv))
             
          CASE ("mean ")
             !
             ! 1. force MAXVAL of laimap on lai on this PFT
             !
             DO ji = 1,kjpindex
                lai(ji,jv) = MAXVAL(laimap(ji,jv,:))
             ENDDO
             !
          CASE ("inter")
             !
             ! 2. do the interpolation between laimax and laimin
             !
             !
             ! If January
             !
             IF (month_end .EQ. 1 ) THEN
                IF (day_end .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,12)*(1-(day_end+15)/30.) + laimap(:,jv,1)*((day_end+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,1)*(1-(day_end-15)/30.) + laimap(:,jv,2)*((day_end-15)/30.)
                ENDIF
                !
                ! If December
                !
             ELSE IF (month_end .EQ. 12) THEN
                IF (day_end .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,11)*(1-(day_end+15)/30.) + laimap(:,jv,12)*((day_end+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,12)*(1-(day_end-15)/30.) + laimap(:,jv,1)*((day_end-15)/30.)
                ENDIF
          !
          ! ELSE
          !
             ELSE
                IF (day_end .LE. 15) THEN
                   lai(:,jv) = laimap(:,jv,month_end-1)*(1-(day_end+15)/30.) + laimap(:,jv,month_end)*((day_end+15)/30.)
                ELSE
                   lai(:,jv) = laimap(:,jv,month_end)*(1-(day_end-15)/30.) + laimap(:,jv,month_end+1)*((day_end-15)/30.)
                ENDIF
             ENDIF
             !
          CASE default
             !
             ! 3. Problem
             !
             WRITE (numout,*) 'This kind of lai choice is not possible. '// &
                  ' We stop with type_of_lai ',jv,' = ', type_of_lai(jv) 
             CALL ipslerr_p(3,'slowproc_lai','Bad value for type_of_lai','read_lai=true','')
          END SELECT
          
       ENDDO
    ENDIF

  END SUBROUTINE slowproc_lai

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_interlai
!!
!>\BRIEF         Interpolate the LAI map to the grid of the model 
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::laimap
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_interlai(nbpt, lalo, resolution, neighbours, contfrac, laimap)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)             :: lalo(nbpt,2)          !! Vector of latitude and longitudes 
                                                                 !! (beware of the order = 1 : latitude, 2 : longitude)
    REAL(r_std), INTENT(in)             :: resolution(nbpt,2)    !! The size in km of each grid-box in X and Y
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point
                                                                 !! (1=North and then clockwise)
    REAL(r_std), INTENT(in)             :: contfrac(nbpt)        !! Fraction of land in each grid box.
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  laimap(nbpt,nvm,12)          !! lai read variable and re-dimensioned
    !
    !  0.3 LOCAL
    !
    CHARACTER(LEN=80) :: filename                               !! name of the LAI map read
    INTEGER(i_std) :: ib, ip, jp, it, jv
    REAL(r_std) :: lmax, lmin, ldelta
    LOGICAL ::           renormelize_lai  ! flag to force LAI renormelization
    INTEGER                  :: ier

    REAL(r_std), DIMENSION(nbpt)                         :: alaimap          !! availability of the lai interpolation 
    INTEGER, DIMENSION(4)                                :: invardims
    REAL(r_std), DIMENSION(nbpt,nvm,12)                  :: lairefrac        !! lai fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt,nvm,12)                  :: fraclaiinterp    !! lai fractions re-dimensioned
    REAL(r_std), DIMENSION(:), ALLOCATABLE               :: vmin, vmax       !! min/max values to use for the 
                                                                             !!   renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
!_ ================================================================================================================================

    !
    !Config Key   = LAI_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = LAI_MAP
    !Config Def   = lai2D.nc
    !Config Help  = The name of the file to be opened to read the LAI
    !Config         map is to be given here. Usualy SECHIBA runs with a 5kmx5km
    !Config         map which is derived from a Nicolas VIOVY one. 
    !Config Units = [FILE]
    !
    filename = 'lai2D.nc'
    CALL getin_p('LAI_FILE',filename)
    variablename = 'LAI'

    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_interlai: Use XIOS to read and interpolate " &
            // TRIM(filename) //" for variable " //TRIM(variablename)
    
       CALL xios_orchidee_recv_field('lai_interp',lairefrac)
       CALL xios_orchidee_recv_field('frac_lai_interp',fraclaiinterp)      
       alaimap(:) = fraclaiinterp(:,1,1)
    ELSE

      IF (printlev_loc >= 2) WRITE(numout,*) "slowproc_interlai: Start interpolate " &
           // TRIM(filename) //" for variable " //TRIM(variablename)

      ! invardims: shape of variable in input file to interpolate
      invardims = interpweight_get_var4dims_file(filename, variablename)
      ! Check coherence of dimensions read from the file
      IF (invardims(4) /= 12)  CALL ipslerr_p(3,'slowproc_interlai','Wrong dimension of time dimension in input file for lai','','')
      IF (invardims(3) /= nvm) CALL ipslerr_p(3,'slowproc_interlai','Wrong dimension of PFT dimension in input file for lai','','')

      ALLOCATE(vmin(nvm),stat=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'slowproc_interlai','Problem in allocation of variable vmin','','')

      ALLOCATE(vmax(nvm), STAT=ier)
      IF (ier /= 0) CALL ipslerr_p(3,'slowproc_interlai','Problem in allocation of variable vmax','','')


! Assigning values to vmin, vmax
      vmin = un
      vmax = nvm*un

      variabletypevals = -un

      !! Variables for interpweight
      ! Type of calculation of cell fractions
      fractype = 'default'
      ! Name of the longitude and latitude in the input file
      lonname = 'longitude'
      latname = 'latitude'
      ! Should negative values be set to zero from input file?
      nonegative = .TRUE.
      ! Type of mask to apply to the input data (see header for more details)
      maskingtype = 'mbelow'
      ! Values to use for the masking
      maskvals = (/ 20., undef_sechiba, undef_sechiba /)
      ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
      namemaskvar = ''

      CALL interpweight_4D(nbpt, nvm, variabletypevals, lalo, resolution, neighbours,        &
        contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
        maskvals, namemaskvar, nvm, invardims(4), -1, fractype,                            &
        -1., -1., lairefrac, alaimap)

      IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_interlai after interpweight_4D'

    ENDIF



    !
    !
    !Config Key   = RENORM_LAI
    !Config Desc  = flag to force LAI renormelization
    !Config If    = LAI_MAP
    !Config Def   = n
    !Config Help  = If true, the laimap will be renormalize between llaimin and llaimax parameters.
    !Config Units = [FLAG]
    !
    renormelize_lai = .FALSE.
    CALL getin_p('RENORM_LAI',renormelize_lai)

    !
    laimap(:,:,:) = zero
    !
    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_interlai before starting loop nbpt:', nbpt
    END IF 

    ! Assigning the right values and giving a value where information was not found
    DO ib=1,nbpt
      IF (alaimap(ib) < min_sechiba) THEN
        DO jv=1,nvm
          laimap(ib,jv,:) = (llaimax(jv)+llaimin(jv))/deux
        ENDDO
      ELSE
        DO jv=1, nvm
          DO it=1, 12
            laimap(ib,jv,it) = lairefrac(ib,jv,it)
          ENDDO
        ENDDO
      END IF
    ENDDO
    !
    ! Normelize the read LAI by the values SECHIBA is used to
    !
    IF ( renormelize_lai ) THEN
       DO ib=1,nbpt
          DO jv=1, nvm
             lmax = MAXVAL(laimap(ib,jv,:))
             lmin = MINVAL(laimap(ib,jv,:))
             ldelta = lmax-lmin
             IF ( ldelta < min_sechiba) THEN
                ! LAI constante ... keep it constant
                laimap(ib,jv,:) = (laimap(ib,jv,:)-lmin)+(llaimax(jv)+llaimin(jv))/deux
             ELSE
                laimap(ib,jv,:) = (laimap(ib,jv,:)-lmin)/(lmax-lmin)*(llaimax(jv)-llaimin(jv))+llaimin(jv)
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_alaimap",alaimap)
    CALL xios_orchidee_send_field("interp_diag_lai",laimap)
   
    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_interlai ended'

  END SUBROUTINE slowproc_interlai

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_readvegetmax
!!
!>\BRIEF          Read and interpolate a vegetation map (by pft)
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): The subroutine was previously called slowproc_update.
!!
!! MAIN OUTPUT VARIABLE(S): 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_readvegetmax(nbpt, lalo, neighbours,  resolution, contfrac, veget_last,         & 
       veget_next, frac_nobio_next, veget_year, init)

    USE interpweight
    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs 
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    REAL(r_std), DIMENSION(nbpt,nvm), INTENT(in)           :: veget_last      !! old max vegetfrac
    INTEGER(i_std), INTENT(in)         :: veget_year            !! first year for landuse (0 == NO TIME AXIS)
    LOGICAL, INTENT(in)                :: init                  !! initialisation : in case of dgvm, it forces update of all PFTs
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt,nvm), INTENT(out)          :: veget_next       !! new max vegetfrac
    REAL(r_std), DIMENSION(nbpt,nnobio), INTENT(out)       :: frac_nobio_next  !! new fraction of the mesh which is 
                                                                               !! covered by ice, lakes, ...
    
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib, inobio, jv
    REAL(r_std) :: sumf, err, norm
    !
    ! for DGVM case :
    REAL(r_std)                 :: sum_veg                     ! sum of vegets
    REAL(r_std)                 :: sum_nobio                   ! sum of nobios
    REAL(r_std)                 :: sumvAnthro_old, sumvAnthro  ! last an new sum of antrhopic vegets
    REAL(r_std)                 :: rapport                     ! (S-B) / (S-A)
    LOGICAL                     :: partial_update              ! if TRUE, partialy update PFT (only anthropic ones) 
                                                               ! e.g. in case of DGVM and not init (optional parameter)
    REAL(r_std), DIMENSION(nbpt,nvm)                     :: vegetrefrac      !! veget fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: aveget           !! Availability of the soilcol interpolation
    REAL(r_std), DIMENSION(nbpt,nvm)                     :: aveget_nvm       !! Availability of the soilcol interpolation
    REAL(r_std), DIMENSION(nvm)                          :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    CHARACTER(LEN=250)                                   :: msg

!_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_readvegetmax'

    !
    !Config Key   = VEGETATION_FILE
    !Config Desc  = Name of file from which the vegetation map is to be read
    !Config If    = 
    !Config Def   = PFTmap.nc
    !Config Help  = The name of the file to be opened to read a vegetation
    !Config         map (in pft) is to be given here. 
    !Config Units = [FILE]
    !
    filename = 'PFTmap.nc'
    CALL getin_p('VEGETATION_FILE',filename)
    variablename = 'maxvegetfrac'


    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readvegetmax: Use XIOS to read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       CALL xios_orchidee_recv_field('frac_veget',vegetrefrac)
       CALL xios_orchidee_recv_field('frac_veget_frac',aveget_nvm)
       aveget(:)=aveget_nvm(:,1)
       
       DO ib = 1, nbpt
          IF (aveget(ib) > min_sechiba) THEN
             vegetrefrac(ib,:) = vegetrefrac(ib,:)/aveget(ib) ! intersected area normalization
             vegetrefrac(ib,:) = vegetrefrac(ib,:)/SUM(vegetrefrac(ib,:))
          ENDIF
       ENDDO
       
    ELSE

      IF (printlev_loc >= 2) WRITE(numout,*) "slowproc_readvegetmax: Start interpolate " &
           // TRIM(filename) // " for variable " // TRIM(variablename)

      ! Assigning values to vmin, vmax
      vmin = 1
      vmax = nvm*1._r_std

      variabletypevals = -un

      !! Variables for interpweight
      ! Type of calculation of cell fractions
      fractype = 'default'
      ! Name of the longitude and latitude in the input file
      lonname = 'lon'
      latname = 'lat'
      ! Should negative values be set to zero from input file?
      nonegative = .FALSE.
      ! Type of mask to apply to the input data (see header for more details)
      maskingtype = 'msumrange'
      ! Values to use for the masking
      maskvals = (/ 1.-1.e-7, 0., 2. /)
      ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
      namemaskvar = ''

      CALL interpweight_3D(nbpt, nvm, variabletypevals, lalo, resolution, neighbours,        &
        contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
        maskvals, namemaskvar, nvm, 0, veget_year, fractype,                                 &
        -1., -1., vegetrefrac, aveget)
      IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_readvegetmax after interpeeight_3D'
    ENDIF 
    !
    ! Compute the logical for partial (only anthropic) PTFs update
    IF (ok_dgvm .AND. .NOT. init) THEN
       partial_update= .TRUE.
    ELSE
       partial_update=.FALSE.
    END IF

    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_readvegetmax before updating loop nbpt:', nbpt
    END IF

    IF ( .NOT. partial_update ) THEN
       ! Case for not DGVM or (DGVM and init)
       veget_next(:,:)=zero
       
       IF (printlev_loc >=3 .AND. ANY(aveget < min_sechiba)) THEN
          WRITE(numout,*) 'Some grid cells on the model grid did not have any points on the source grid.'
          IF (init) THEN
             WRITE(numout,*) 'Initialization with full fraction of bare soil are done for the below grid cells.'
          ELSE
             WRITE(numout,*) 'Old values are kept for the below grid cells.'
          ENDIF
          WRITE(numout,*) 'List of grid cells (ib, lat, lon):'
       END IF
 
      DO ib = 1, nbpt
          ! vegetrefrac is already normalized to sum equal one for each grid cell
          veget_next(ib,:) = vegetrefrac(ib,:)

          IF (aveget(ib) < min_sechiba) THEN
             IF (printlev_loc >=3) WRITE(numout,*) ib,lalo(ib,1),lalo(ib,2)
             IF (init) THEN
                veget_next(ib,1) = un
                veget_next(ib,2:nvm) = zero
             ELSE
                veget_next(ib,:) = veget_last(ib,:)
             ENDIF
          ENDIF
       ENDDO
    ELSE
       ! Partial update
       DO ib = 1, nbpt
          IF (aveget(ib) > min_sechiba) THEN
             ! For the case with properly interpolated grid cells (aveget>0)

             ! last veget for this point
             sum_veg=SUM(veget_last(ib,:))
             !
             ! If the DGVM is activated, only anthropic PFTs are utpdated, the others are copied from previous time-step 
             veget_next(ib,:) = veget_last(ib,:)
             
             DO jv = 2, nvm
                IF ( .NOT. natural(jv) ) THEN       
                   veget_next(ib,jv) = vegetrefrac(ib,jv)
                ENDIF
             ENDDO

             sumvAnthro_old = zero
             sumvAnthro     = zero
             DO jv = 2, nvm
                IF ( .NOT. natural(jv) ) THEN
                   sumvAnthro = sumvAnthro + veget_next(ib,jv)
                   sumvAnthro_old = sumvAnthro_old + veget_last(ib,jv)
                ENDIF
             ENDDO

             IF ( sumvAnthro_old < sumvAnthro ) THEN
                ! Increase of non natural vegetations (increase of agriculture)
                ! The proportion of natural PFT's must be preserved
                ! ie the sum of vegets is preserved
                !    and natural PFT / (sum of veget - sum of antropic veget)
                !    is preserved. 
                rapport = ( sum_veg - sumvAnthro ) / ( sum_veg - sumvAnthro_old )
                DO jv = 1, nvm
                   IF ( natural(jv) ) THEN
                      veget_next(ib,jv) = veget_last(ib,jv) * rapport
                   ENDIF
                ENDDO
             ELSE
                ! Increase of natural vegetations (decrease of agriculture)
                ! The decrease of agriculture is replaced by bare soil. The DGVM will
                ! re-introduce natural PFT's.
                DO jv = 1, nvm
                   IF ( natural(jv) ) THEN
                      veget_next(ib,jv) = veget_last(ib,jv)
                   ENDIF
                ENDDO
                veget_next(ib,1) = veget_next(ib,1) + sumvAnthro_old - sumvAnthro
             ENDIF

             ! test
             IF ( ABS( SUM(veget_next(ib,:)) - sum_veg ) > 10*EPSILON(un) ) THEN
                WRITE(numout,*) 'slowproc_readvegetmax _______'
                msg = "  No conservation of sum of veget for point "
                WRITE(numout,*) TRIM(msg), ib, ",(", lalo(ib,1),",", lalo(ib,2), ")" 
                WRITE(numout,*) "  last sum of veget ", sum_veg, " new sum of veget ",                &
                  SUM(veget_next(ib,:)), " error : ", SUM(veget_next(ib,:))-sum_veg
                WRITE(numout,*) "  Anthropic modifications : last ",sumvAnthro_old," new ",sumvAnthro     
                CALL ipslerr_p (3,'slowproc_readvegetmax',                                            &
                     &          'No conservation of sum of veget_next',                               &
                     &          "The sum of veget_next is different after reading Land Use map.",     &
                     &          '(verify the dgvm case model.)')
             ENDIF
          ELSE
             ! For the case when there was a propblem with the interpolation, aveget < min_sechiba
             WRITE(numout,*) 'slowproc_readvegetmax _______'
             WRITE(numout,*) "  No land point in the map for point ", ib, ",(", lalo(ib,1), ",",      &
               lalo(ib,2),")" 
             CALL ipslerr_p (2,'slowproc_readvegetmax',                                               &
                  &          'Problem with vegetation file for Land Use.',                            &
                  &          "No land point in the map for point",                                    & 
                  &          '(verify your land use file.)')
             veget_next(ib,:) = veget_last(ib,:)
          ENDIF
          
       ENDDO
    ENDIF

    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_readvegetmax after updating'
    !
    frac_nobio_next (:,:) = un
    !
!MM
    ! Work only for one nnobio !! (ie ice)
    DO inobio=1,nnobio
       DO jv=1,nvm
          DO ib = 1, nbpt
             frac_nobio_next(ib,inobio) = frac_nobio_next(ib,inobio) - veget_next(ib,jv)
          ENDDO
       ENDDO
    ENDDO

    DO ib = 1, nbpt
       sum_veg = SUM(veget_next(ib,:))
       sum_nobio = SUM(frac_nobio_next(ib,:))
       IF (sum_nobio < 0.) THEN
          frac_nobio_next(ib,:) = zero
          veget_next(ib,1) = veget_next(ib,1) + sum_nobio
          sum_veg = SUM(veget_next(ib,:))
       ENDIF
       sumf = sum_veg + sum_nobio
       IF (sumf > min_sechiba) THEN
          veget_next(ib,:) = veget_next(ib,:) / sumf
          frac_nobio_next(ib,:) = frac_nobio_next(ib,:) / sumf
          norm=SUM(veget_next(ib,:))+SUM(frac_nobio_next(ib,:))
          err=norm-un
          IF (printlev_loc >=5) WRITE(numout,*) "  slowproc_readvegetmax: ib ",ib,                    &
            " SUM(veget_next(ib,:)+frac_nobio_next(ib,:))-un, sumf",err,sumf
          IF (abs(err) > -EPSILON(un)) THEN
             IF ( SUM(frac_nobio_next(ib,:)) > min_sechiba ) THEN
                frac_nobio_next(ib,1) = frac_nobio_next(ib,1) - err
             ELSE
                veget_next(ib,1) = veget_next(ib,1) - err
             ENDIF
             norm=SUM(veget_next(ib,:))+SUM(frac_nobio_next(ib,:))
             err=norm-un
             IF (printlev_loc >=5) WRITE(numout,*) "  slowproc_readvegetmax: ib ", ib,                &
               " SUM(veget_next(ib,:)+frac_nobio_next(ib,:))-un",err
             IF (abs(err) > EPSILON(un)) THEN
                WRITE(numout,*) '  slowproc_readvegetmax _______'
                WRITE(numout,*) "update : Problem with point ",ib,",(",lalo(ib,1),",",lalo(ib,2),")" 
                WRITE(numout,*) "         err(sum-1.) = ",abs(err)
                CALL ipslerr_p (2,'slowproc_readvegetmax', &
                     &          'Problem with sum vegetation + sum fracnobio for Land Use.',          &
                     &          "sum not equal to 1.", &
                     &          '(verify your land use file.)')
                aveget(ib) = -0.6
             ENDIF
          ENDIF
       ELSE
          ! sumf < min_sechiba
          WRITE(numout,*) '  slowproc_readvegetmax _______'
          WRITE(numout,*)"    No vegetation nor frac_nobio for point ", ib, ",(", lalo(ib,1), ",",    &
            lalo(ib,2),")" 
          WRITE(numout,*)"    Replaced by bare_soil !! "
          veget_next(ib,1) = un
          veget_next(ib,2:nvm) = zero
          frac_nobio_next(ib,:) = zero
!!!$          CALL ipslerr_p (3,'slowproc_readvegetmax', &
!!!$               &          'Problem with vegetation file for Land Use.', &
!!!$               &          "No vegetation nor frac_nobio for point ", &
!!!$               &          '(verify your land use file.)')
       ENDIF
    ENDDO

    !! Set to zero fractions of frac_nobio and veget_max smaller than min_vegfrac
    !! Normalize to have the sum equal 1.
    CALL slowproc_veget_max_limit(nbpt, frac_nobio_next, veget_next)

    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_aveget",aveget)
    CALL xios_orchidee_send_field("interp_diag_vegetrefrac",vegetrefrac)
    CALL xios_orchidee_send_field("interp_diag_veget_next",veget_next)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_readvegetmax ended'
    
  END SUBROUTINE slowproc_readvegetmax


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_nearest
!!
!>\BRIEF         looks for nearest grid point on the fine map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::inear
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_nearest(iml, lon5, lat5, lonmod, latmod, inear)

    !! INTERFACE DESCRIPTION
    
    !! 0.1 input variables

    INTEGER(i_std), INTENT(in)                   :: iml             !! size of the vector
    REAL(r_std), DIMENSION(iml), INTENT(in)      :: lon5, lat5      !! longitude and latitude vector, for the 5km vegmap
    REAL(r_std), INTENT(in)                      :: lonmod, latmod  !! longitude  and latitude modelled

    !! 0.2 output variables
    
    INTEGER(i_std), INTENT(out)                  :: inear           !! location of the grid point from the 5km vegmap grid
                                                                    !! closest from the modelled grid point

    !! 0.4 Local variables

    REAL(r_std)                                  :: pa, p
    REAL(r_std)                                  :: coscolat, sincolat
    REAL(r_std)                                  :: cospa, sinpa
    REAL(r_std), ALLOCATABLE, DIMENSION(:)       :: cosang
    INTEGER(i_std)                               :: i
    INTEGER(i_std), DIMENSION(1)                 :: ineartab
    INTEGER                                      :: ALLOC_ERR

!_ ================================================================================================================================

    ALLOCATE(cosang(iml), STAT=ALLOC_ERR)
    IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_nearest','Error in allocation for cosang','','')

    pa = pi/2.0 - latmod*pi/180.0 ! dist. between north pole and the point a 
                                                      !! COLATITUDE, in radian
    cospa = COS(pa)
    sinpa = SIN(pa)

    DO i = 1, iml

       sincolat = SIN( pi/2.0 - lat5(i)*pi/180.0 ) !! sinus of the colatitude
       coscolat = COS( pi/2.0 - lat5(i)*pi/180.0 ) !! cosinus of the colatitude

       p = (lonmod-lon5(i))*pi/180.0 !! angle between a & b (between their meridian)in radians

       !! dist(i) = ACOS( cospa*coscolat + sinpa*sincolat*COS(p))
       cosang(i) = cospa*coscolat + sinpa*sincolat*COS(p) !! TL : cosang is maximum when angle is at minimal value  
!! orthodromic distance between 2 points : cosang = cosinus (arc(AB)/R), with
!R = Earth radius, then max(cosang) = max(cos(arc(AB)/R)), reached when arc(AB)/R is minimal, when
! arc(AB) is minimal, thus when point B (corresponding grid point from LAI MAP) is the nearest from
! modelled A point
    ENDDO

    ineartab = MAXLOC( cosang(:) )
    inear = ineartab(1)

    DEALLOCATE(cosang)
  END SUBROUTINE slowproc_nearest

!! ================================================================================================================================
!! SUBROUTINE   : slowproc_soilt
!!
!>\BRIEF         Interpolate the Zobler or Reynolds/USDA soil type map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): Nov 2014, ADucharne
!!                   Nov 2020, Salma Tafasca and Agnes Ducharne: adding a choice for spmipexp/SPMIPEXP,
!!                             and everything needed to read all maps and assign parameter values.  
!!
!! MAIN OUTPUT VARIABLE(S): ::soiltype, ::clayfraction, sandfraction, siltfraction
!!
!! REFERENCE(S) : Reynold, Jackson, and Rawls (2000). Estimating soil water-holding capacities 
!! by linking the Food and Agriculture Organization soil map of the world with global pedon
!! databases and continuous pedotransfer functions, WRR, 36, 3653-3662
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
  SUBROUTINE slowproc_soilt(njsc,  ks,  nvan, avan, mcr, mcs, mcfc, mcw, nbpt, lalo, neighbours, resolution, contfrac, soilclass, clayfraction, sandfraction, siltfraction)

    USE interpweight

    IMPLICIT NONE
    !
    !
    !   This subroutine should read the Zobler/Reynolds map and interpolate to the model grid. 
    !   The method is to get fraction of the three/12 main soiltypes for each grid box.
    !   For the Zobler case, also called FAO in the code, the soil fraction are going to be put 
    !   into the array soiltype in the following order : coarse, medium and fine.
    !   For the Reynolds/USDA case, the soiltype array follows the order defined in constantes_soil_var.f90
    !
    !
    !!  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)    :: nbpt                   !! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)       :: lalo(nbpt,2)           !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)    :: neighbours(nbpt,NbNeighb)!! Vector of neighbours for each grid point
                                                              !! (1=North and then clockwise)
    REAL(r_std), INTENT(in)       :: resolution(nbpt,2)     !! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT(in)       :: contfrac(nbpt)         !! Fraction of land in each grid box.
    !
    !  0.2 OUTPUT
    !
    !salma: added soil params and njsc because needed in the calculation of the soil params
    INTEGER(i_std),DIMENSION (nbpt), INTENT (out)      :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: ks             !! Hydraulic conductivity at saturation (mm {-1})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: nvan           !! Van Genuchten coeficients n (unitless)
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: avan           !! Van Genuchten coeficients a (mm-1})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: mcr            !! Residual volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: mcs            !! Saturated volumetric water content (m^{3} m^{-3})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: mcfc           !! Volumetric water content at field capacity (m^{3} m^{-3})
    REAL(r_std),DIMENSION (nbpt), INTENT (out)         :: mcw            !! Volumetric water content at wilting point (m^{3} m^{-3})

    REAL(r_std), INTENT(out)      :: soilclass(nbpt, nscm)  !! Soil type map to be created from the Zobler map
                                                            !! or a map defining the 12 USDA classes (e.g. Reynolds)
                                                            !! Holds the area of each texture class in the ORCHIDEE grid cells
                                                            !! Final unit = fraction of ORCHIDEE grid-cell (unitless)
    REAL(r_std), INTENT(out)      :: clayfraction(nbpt)     !! The fraction of clay as used by STOMATE
    REAL(r_std), INTENT(out)      :: sandfraction(nbpt)     !! The fraction of sand (for SP-MIP)
    REAL(r_std), INTENT(out)      :: siltfraction(nbpt)     !! The fraction of silt (for SP-MIP)
    !
    !
    !  0.3 LOCAL
    !
    !salma: added the following local variable to be used for all the soil hydraulic parameters
    REAL(r_std), DIMENSION(nbpt)        :: param            !! to be introduced in function: interpweight

    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib, ilf, nbexp, i
    INTEGER(i_std) :: fopt                                  !! Nb of pts from the texture map within one ORCHIDEE grid-cell
    INTEGER(i_std), ALLOCATABLE, DIMENSION(:) :: solt       !! Texture the different points from the input texture map 
                                                            !! in one ORCHIDEE grid cell (unitless)
    !
    ! Number of texture classes in Zobler
    !
    INTEGER(i_std), PARAMETER :: nzobler = 7                !! Nb of texture classes according in the Zobler map
    REAL(r_std),ALLOCATABLE   :: textfrac_table(:,:)        !! conversion table between the texture index
                                                            !! and the granulometric composition
    !   
    INTEGER                  :: ALLOC_ERR
    INTEGER                                              :: ntextinfile      !! number of soil textures in the in the file
    REAL(r_std), DIMENSION(:,:), ALLOCATABLE             :: textrefrac       !! text fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: atext            !! Availability of the texture interpolation
    !salma added the following 3 variables to control the SP-MIP simulations
    REAL(r_std), DIMENSION(nbpt)                         :: aparam            !! Availability of the parameter interpolation
    CHARACTER(LEN=80)                                    :: spmipexp          !! designing the number of sp-mip experiment
    CHARACTER(LEN=80)                                    :: unif_case               !! designing the model of experiment 4 (sp_mip)

    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the 

    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat name in input file
    REAL(r_std), DIMENSION(:), ALLOCATABLE               :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    INTEGER(i_std), DIMENSION(:), ALLOCATABLE            :: vecpos
    REAL(r_std)                                          :: sgn              !! sum of fractions excluding glaciers and ocean

    ! For the calculation of field capacity and wilting point
    REAL(r_std),DIMENSION (nbpt)                         :: mvan             !! Van Genuchten parameter m
    REAL(r_std),DIMENSION (nbpt)                         :: psi_w            !! Matrix potential characterizing the wilting point (mm)
    REAL(r_std),DIMENSION (nbpt)                         :: psi_fc           !! Matrix potential characterizing the field capacity (mm)
    
!_ ================================================================================================================================

    IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_soilt'

    ! The soil parameters are defined by several keywords in run.def:
    ! (a) soil_classif tells which kind of soil texture map you will read (mandatory):
    !    - usda for 12 USDA texture classes (Reynolds, SoilGrids, SPMIP, etc) updated to 13 classes
    !      for clay oxisols by Salma Tafasca
    !    - zobler to read teh Zobler map and reduce it to 3 classes (fine, medium, coarse)
    ! (b) spmipexp was introduced by Salma Tafasca for the SPMIP project  
    !   maps: Reading the soil parameter maps of SPMIP
    !   unif: Imposing uniform soil texture over the globe (4 texture options, with parameter values imposed by SP-MIP)
    ! Even with maps, some parameters (thermics) are defined based on texture.
    ! So we read a soil texture map in all experiments but unif, where soil texture is imposed by njsc(:).
    ! (c) unif_case to choose the soil texture assigned if spmipexp=maps (4 hard_coded possibilities)
    
    ! IMPORTANT: if no spmipexp is defined in run.def, the model works as before, by deriving the soil parameters 
    ! from a soil texture map, itself defined by the SOILTYPE_CLASSIF keyword, and soil_classif variable
    ! But to get a uniform texture (exp 4), you need to select a soil texture map using soil_classif, even if it's not read

    !Config Key   = SPMIPEXP
    !Config Desc  = Types of alternative hydraulic parameters
    !Config Def   = 'texture'
    !Config If    = 
    !Config Help  = possible values: maps, unif
    !Config Units = [-]
    spmipexp='texture' ! default is to define parameters from soil texture, with soil_classif = 'zobler' or 'usda'
    CALL getin_p("SPMIPEXP",spmipexp)

    IF (spmipexp == 'unif') THEN
       ! case where unif=exp4 is selected: uniform soil parameters
       ! the values of the hydraulic parameters below come from SP-MIP,
       ! and correspond to the Rosetta PTF (Schaap et al., 2001)

       ! sp_mip_experiment_4: select another level of experiment: a, b, c or d in run.def

       !Config Key   = UNIF_CASE
       !Config Desc  = Types of uniform soil textures in SPMIP
       !Config Def   = 'b'
       !Config If    = SPMIPEXP='unif'
       !Config Help  = possible values: a, b, c and d
       !Config Units = [-]
       unif_case='b' ! default = loamy soil
       CALL getin_p("UNIF_CASE",unif_case)

       SELECTCASE (unif_case)

       CASE ('a') ! loamy sand
          clayfraction=0.06
          sandfraction=0.81
          siltfraction=0.13
          DO ib=1 , nbpt
             njsc(ib) = 2
             mcr(ib) = 0.049
             mcs(ib) = 0.39
             ks(ib) = (1.41e-5)*1000*24*3600
             avan(ib) = 3.475*(1e-3)
             nvan(ib) = 1.746
             mcfc(ib) = 0.1039
             mcw(ib) = 0.05221
          ENDDO

       CASE ('b') !loam
          clayfraction=0.2
          sandfraction=0.4
          siltfraction=0.4
          DO ib=1, nbpt
             njsc(ib) = 6
             mcr(ib) = 0.061
             mcs(ib) = 0.399
             ks(ib) = (3.38e-6)*1000*24*3600
             avan(ib) = 1.112*(1e-3)
             nvan(ib) = 1.472
             mcfc(ib) = 0.236
             mcw(ib) = 0.09115
          ENDDO

       CASE ('c') !silt
          clayfraction=0.1
          sandfraction=0.06
          siltfraction=0.84
          DO ib=1, nbpt
             njsc(ib)=5
             mcr(ib) = 0.05
             mcs(ib) = 0.489
             ks(ib) = (2.81e-6)*1000*24*3600
             avan(ib) = 0.6577*(1e-3)
             nvan(ib) = 1.679
             mcfc(ib) = 0.2854
             mcw(ib) = 0.06944
          ENDDO

       CASE ('d')!clay
          clayfraction=0.55
          sandfraction=0.15
          siltfraction=0.3
          DO ib=1, nbpt
             njsc(ib)=12
             mcr(ib) = 0.098
             mcs(ib) = 0.459
             ks(ib) = (9.74e-7)*1000*24*3600
             avan(ib) = 1.496*(1e-3)
             nvan(ib) = 1.253
             mcfc(ib) = 0.3329
             mcw(ib) = 0.1897
          ENDDO

       CASE DEFAULT

          WRITE (numout,*) 'Unsupported experiment number. Choose between a, b, c or d according to sp_mip_experiment_4 number'
          CALL ipslerr_p(3,'hydrol_init','Unsupported experiment number. ',&
               'Choose between a,b,c or d','')
       ENDSELECT

    ELSE ! spmipexp is either exp1=maps, or texture for exp2 or exp3 (or typing error!)
        
       ! In these cases (maps or texture), we need to read the soil texture map
       
       !Config Key   = SOILCLASS_FILE
       !Config Desc  = Name of file from which soil types are read
       !Config Def   = soils_param.nc
       !Config If    = NOT(IMPOSE_VEG)
       !Config Help  = The name of the file to be opened to read the soil types.
       !Config         The data from this file is then interpolated to the grid of
       !Config         of the model. The aim is to get fractions for sand loam and
       !Config         clay in each grid box. This information is used for soil hydrology
       !Config         and respiration.
       !Config Units = [FILE]
       !
       ! soils_param.nc file is 1deg soil texture file (Zobler)
       ! The USDA map from Reynolds is soils_param_usda.nc (1/12deg resolution)

       filename = 'soils_param.nc'
       CALL getin_p('SOILCLASS_FILE',filename)

       variablename = 'soiltext'

       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'

       IF (xios_interpolation) THEN
          IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_soilt: Use XIOS to read and interpolate " &
               // TRIM(filename) // " for variable " // TRIM(variablename)

          SELECT CASE(soil_classif)

          CASE('none')
             ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
             IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
             DO ib=1, nbpt
                njsc(ib) = usda_default ! 6 = Loam
                clayfraction(ib) = clayfrac_usda(usda_default)
                sandfraction(ib) = sandfrac_usda(usda_default)
                siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)
             ENDDO

          CASE('zobler')
             !             !
             IF (printlev_loc>=2) WRITE(numout,*) "Using a soilclass map with Zobler classification, to be read using XIOS"
             !
             ALLOCATE(textrefrac(nbpt,nzobler))
             ALLOCATE(textfrac_table(nzobler,ntext), STAT=ALLOC_ERR)
             IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
             CALL get_soilcorr_zobler (nzobler, textfrac_table)
             CALL xios_orchidee_recv_field('soiltext1',textrefrac(:,1))
             CALL xios_orchidee_recv_field('soiltext2',textrefrac(:,2))
             CALL xios_orchidee_recv_field('soiltext3',textrefrac(:,3))
             CALL xios_orchidee_recv_field('soiltext4',textrefrac(:,4))
             CALL xios_orchidee_recv_field('soiltext5',textrefrac(:,5))
             CALL xios_orchidee_recv_field('soiltext6',textrefrac(:,6))
             CALL xios_orchidee_recv_field('soiltext7',textrefrac(:,7))

             CALL get_soilcorr_zobler (nzobler, textfrac_table)
             !             !
             DO ib =1, nbpt
                soilclass(ib,:)=0.
                soilclass(ib,fao2usda(1))=textrefrac(ib,1)
                soilclass(ib,fao2usda(2))=textrefrac(ib,2)+textrefrac(ib,3)+textrefrac(ib,4)+textrefrac(ib,7)
                soilclass(ib,fao2usda(3))=textrefrac(ib,5)

                ! clayfraction is the sum of the % of clay (as a mineral of small granulometry, and not as a texture)
                ! over the zobler pixels composing the ORCHIDEE grid-cell
                clayfraction(ib) = textfrac_table(1,3) * textrefrac(ib,1)+textfrac_table(2,3) * textrefrac(ib,2) + &
                     textfrac_table(3,3) * textrefrac(ib,3)+textfrac_table(4,3) * textrefrac(ib,4) + &
                     textfrac_table(5,3) * textrefrac(ib,5)+textfrac_table(7,3) * textrefrac(ib,7)

                sandfraction(ib) = textfrac_table(1,2) * textrefrac(ib,1)+textfrac_table(2,2) * textrefrac(ib,2) + &
                     textfrac_table(3,2) * textrefrac(ib,3)+textfrac_table(4,2) * textrefrac(ib,4) + &
                     textfrac_table(5,2) * textrefrac(ib,5)+textfrac_table(7,2) * textrefrac(ib,7)

                siltfraction(ib) = textfrac_table(1,1) * textrefrac(ib,1)+textfrac_table(2,1) * textrefrac(ib,2) + &
                     textfrac_table(3,1) * textrefrac(ib,3)+textfrac_table(4,1) * textrefrac(ib,4) + &
                     textfrac_table(5,1) * textrefrac(ib,5)+textfrac_table(7,1) * textrefrac(ib,7)

                sgn=SUM(soilclass(ib,:)) ! grid-cell fraction with texture info

                IF (sgn < min_sechiba) THEN ! if no texture info in this grid-point, we assume that texture = Loam
                   njsc(ib) = usda_default ! 6 = Loam
                   clayfraction(ib) = clayfrac_usda(usda_default)
                   sandfraction(ib) = sandfrac_usda(usda_default)
                   siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)
                   atext(ib)=0.
                ELSE
                   atext(ib)=sgn
                   clayfraction(ib) = clayfraction(ib) / sgn
                   sandfraction(ib) = sandfraction(ib) / sgn
                   siltfraction(ib) = siltfraction(ib) / sgn
                   soilclass(ib,:) = soilclass(ib,:) / sgn
                   njsc(ib) = MAXLOC(soilclass(ib,:),1) ! Dominant texture class
                ENDIF  

             ENDDO

          CASE('usda')

             IF (printlev_loc>=4) WRITE (numout,*) 'slowproc_soilt: start case usda'
             !
             WRITE(numout,*) "Using a soilclass map with usda classification, to be read using XIOS"
             !
             ALLOCATE(textrefrac(nbpt,nscm))
             ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
             IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')

             CALL get_soilcorr_usda (nscm, textfrac_table)

             IF (printlev_loc>=4) WRITE (numout,*) 'slowproc_soilt: After get_soilcorr_usda'

             CALL xios_orchidee_recv_field('soiltext1',textrefrac(:,1))
             CALL xios_orchidee_recv_field('soiltext2',textrefrac(:,2))
             CALL xios_orchidee_recv_field('soiltext3',textrefrac(:,3))
             CALL xios_orchidee_recv_field('soiltext4',textrefrac(:,4))
             CALL xios_orchidee_recv_field('soiltext5',textrefrac(:,5))
             CALL xios_orchidee_recv_field('soiltext6',textrefrac(:,6))
             CALL xios_orchidee_recv_field('soiltext7',textrefrac(:,7))
             CALL xios_orchidee_recv_field('soiltext8',textrefrac(:,8))
             CALL xios_orchidee_recv_field('soiltext9',textrefrac(:,9))
             CALL xios_orchidee_recv_field('soiltext10',textrefrac(:,10))
             CALL xios_orchidee_recv_field('soiltext11',textrefrac(:,11))
             CALL xios_orchidee_recv_field('soiltext12',textrefrac(:,12))        
             CALL xios_orchidee_recv_field('soiltext13',textrefrac(:,13))

             CALL get_soilcorr_usda (nscm, textfrac_table)
             IF (printlev_loc>=4) WRITE (numout,*) 'slowproc_soilt: After get_soilcorr_usda'

             DO ib =1, nbpt
                clayfraction(ib) = 0.0
                DO ilf = 1,nscm
                   soilclass(ib,ilf)=textrefrac(ib,ilf)
                   clayfraction(ib) = clayfraction(ib) + textfrac_table(ilf,3)*textrefrac(ib,ilf)
                   sandfraction(ib) = sandfraction(ib) + textfrac_table(ilf,2)*textrefrac(ib,ilf)
                   siltfraction(ib) = siltfraction(ib) + textfrac_table(ilf,1)*textrefrac(ib,ilf)
                   ! textfrac_table holds the %silt,%sand,%clay
                ENDDO

                sgn=SUM(soilclass(ib,:)) ! grid-cell fraction with texture info

                IF (sgn < min_sechiba) THEN ! if no texture info in this grid-point, we assume that texture = Loam
                   njsc(ib) = usda_default ! 6 = Loam
                   clayfraction(ib) = clayfrac_usda(usda_default)
                   sandfraction(ib) = sandfrac_usda(usda_default)
                   siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)
                   atext(ib)=0
                ELSE
                   soilclass(ib,:) = soilclass(ib,:) / sgn
                   clayfraction(ib) = clayfraction(ib) / sgn
                   sandfraction(ib) = sandfraction(ib) / sgn
                   siltfraction(ib) = siltfraction(ib) / sgn
                   atext(ib)=sgn
                   njsc(ib) = MAXLOC(soilclass(ib,:),1) ! Dominant texture class
                ENDIF               
                    
             ENDDO

          CASE DEFAULT
             WRITE(numout,*) 'slowproc_soilt:'
             WRITE(numout,*) '  A non supported soil type classification has been chosen'
             CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
          END SELECT

       ELSE              !    xios_interpolation
          ! Read and interpolate using stardard method with IOIPSL and aggregate

          IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_soilt: Read and interpolate " &
               // TRIM(filename) // " for variable " // TRIM(variablename)

          ! Name of the longitude and latitude in the input file
          lonname = 'nav_lon'
          latname = 'nav_lat'

          IF (printlev_loc >= 2) WRITE(numout,*) "slowproc_soilt: Start interpolate " &
               // TRIM(filename) // " for variable " // TRIM(variablename)

          IF ( TRIM(soil_classif) /= 'none' ) THEN

             ! Define a variable for the number of soil textures in the input file
             SELECTCASE(soil_classif)
             CASE('zobler')
                ntextinfile=nzobler
             CASE('usda')
                ntextinfile=nscm
             CASE DEFAULT
                WRITE(numout,*) 'slowproc_soilt:'
                WRITE(numout,*) '  A non supported soil type classification has been chosen'
                CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
             ENDSELECT

             ALLOCATE(textrefrac(nbpt,ntextinfile), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable textrefrac',&
                  '','')

             ! Assigning values to vmin, vmax
             vmin = un
             vmax = ntextinfile*un

             ALLOCATE(variabletypevals(ntextinfile), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variabletypevals','','')
             variabletypevals = -un

             !! Variables for interpweight
             ! Should negative values be set to zero from input file?
             nonegative = .FALSE.
             ! Type of mask to apply to the input data (see header for more details)
             maskingtype = 'mabove'
             ! Values to use for the masking
             maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
             ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') ( not used)

             namemaskvar = ''

             CALL interpweight_2D(nbpt, ntextinfile, variabletypevals, lalo, resolution, neighbours,        &
                  contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,    &
                  maskvals, namemaskvar, 0, 0, -1, fractype, -1., -1., textrefrac, atext)

             ALLOCATE(vecpos(ntextinfile), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable vecpos','','')
             ALLOCATE(solt(ntextinfile), STAT=ALLOC_ERR)
             IF (ALLOC_ERR /= 0) CALL ipslerr_p(3,'slowproc_soilt','Problem in allocation of variable solt','','')

             IF (printlev_loc >= 5) THEN
                WRITE(numout,*)'  slowproc_soilt after interpweight_2D'
                WRITE(numout,*)'  slowproc_soilt before starting loop nbpt:', nbpt
                WRITE(numout,*)"  slowproc_soilt starting classification '" // TRIM(soil_classif) // "'..."
             END IF
          ELSE
             IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_soilt using default values all points are propertly ' // &
                  'interpolated atext = 1. everywhere!'
             atext = 1.
          END IF

          nbexp = 0
          SELECTCASE(soil_classif)
          CASE('none')
             ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
             IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
             DO ib=1, nbpt
                njsc(ib) = usda_default ! 6 = Loam
                clayfraction(ib) = clayfrac_usda(usda_default)
                sandfraction(ib) = sandfrac_usda(usda_default)
                siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)
             ENDDO
          CASE('zobler')
             !             !
             IF (printlev_loc>=2) WRITE(numout,*) "Using a soilclass map with Zobler classification"
             !
             ALLOCATE(textfrac_table(nzobler,ntext), STAT=ALLOC_ERR)
             IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')
             CALL get_soilcorr_zobler (nzobler, textfrac_table)
                       
             IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_soilt after getting table of textures'
             DO ib =1, nbpt
                soilclass(ib,:) = zero
                clayfraction(ib) = zero
                sandfraction(ib) = zero
                siltfraction(ib) = zero
                !
                ! vecpos: List of positions where textures were not zero
                ! vecpos(1): number of not null textures found
                vecpos = interpweight_ValVecR(textrefrac(ib,:),nzobler,zero,'neq')
                fopt = vecpos(1)

                IF ( fopt .EQ. 0 ) THEN
                   ! No points were found for current grid box, use default values
                   nbexp = nbexp + 1
                   njsc(ib) = usda_default ! 6=Loam
                   clayfraction(ib) = clayfrac_usda(usda_default)
                   sandfraction(ib) = sandfrac_usda(usda_default)
                   siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)

                ELSE
                   IF (fopt == nzobler) THEN
                      ! All textures are not zero
                      solt=(/(i,i=1,nzobler)/)
                   ELSE
                      DO ilf = 1,fopt
                         solt(ilf) = vecpos(ilf+1)
                      END DO
                   END IF
                   !
                   !   Compute the fraction of each textural class
                   !
                   sgn = 0.
                   DO ilf = 1,fopt
                      !
                      ! Here we make the correspondance between the 7 zobler textures and the 3 textures in ORCHIDEE
                      ! and soilclass correspond to surfaces covered by the 3 textures of ORCHIDEE (coase,medium,fine)
                      ! For type 6 = glacier, default values are set and it is also taken into account during the normalization
                      ! of the fractions (done in interpweight_2D)
                      ! Note that type 0 corresponds to ocean but it is already removed using the mask above.
                     !
                      IF ( (solt(ilf) .LE. nzobler) .AND. (solt(ilf) .GT. 0) .AND. &
                           (solt(ilf) .NE. 6) ) THEN
                         SELECT CASE(solt(ilf))
                         CASE(1)
                            soilclass(ib,fao2usda(1)) = soilclass(ib,fao2usda(1)) + textrefrac(ib,solt(ilf))
                         CASE(2)
                            soilclass(ib,fao2usda(2)) = soilclass(ib,fao2usda(2)) + textrefrac(ib,solt(ilf))
                         CASE(3)
                            soilclass(ib,fao2usda(2)) = soilclass(ib,fao2usda(2)) + textrefrac(ib,solt(ilf))
                         CASE(4)
                            soilclass(ib,fao2usda(2)) = soilclass(ib,fao2usda(2)) + textrefrac(ib,solt(ilf))
                         CASE(5)
                            soilclass(ib,fao2usda(3)) = soilclass(ib,fao2usda(3)) + textrefrac(ib,solt(ilf))
                         CASE(7)
                            soilclass(ib,fao2usda(2)) = soilclass(ib,fao2usda(2)) + textrefrac(ib,solt(ilf))
                         CASE DEFAULT
                            WRITE(numout,*) 'We should not be here, an impossible case appeared'
                            CALL ipslerr_p(3,'slowproc_soilt','Bad value for solt','','')
                         END SELECT
                         ! clayfraction is the sum of the % of clay (as a mineral of small granulometry, and not as a texture)
                         ! over the zobler pixels composing the ORCHIDEE grid-cell
                         clayfraction(ib) = clayfraction(ib) + &
                              & textfrac_table(solt(ilf),3) * textrefrac(ib,solt(ilf))
                         sandfraction(ib) = sandfraction(ib) + &
                              & textfrac_table(solt(ilf),2) * textrefrac(ib,solt(ilf))
                         siltfraction(ib) = siltfraction(ib) + &
                              & textfrac_table(solt(ilf),1) * textrefrac(ib,solt(ilf))
                         ! Sum the fractions which are not glaciers nor ocean
                         sgn = sgn + textrefrac(ib,solt(ilf))
                      ELSE
                         IF (solt(ilf) .GT. nzobler) THEN
                            WRITE(numout,*) 'The file contains a soil color class which is incompatible with this program'
                            CALL ipslerr_p(3,'slowproc_soilt','Problem soil color class incompatible','','')
                         ENDIF
                      END IF
                   ENDDO

                   IF ( (sgn .LT. min_sechiba) .OR. (atext(ib) .LT. min_sechiba) ) THEN
                      ! Set default values if grid cells were only covered by glaciers or ocean
                      ! or if now information on the source grid was found (atext(ib)=-1).
                      nbexp = nbexp + 1
                      njsc(ib) = usda_default ! 6 = Loam
                      clayfraction(ib) = clayfrac_usda(usda_default)
                      sandfraction(ib) = sandfrac_usda(usda_default)
                      siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)
                   ELSE
                      ! Normalize using the fraction of surface not including glaciers and ocean
                      soilclass(ib,:) = soilclass(ib,:)/sgn
                      clayfraction(ib) = clayfraction(ib)/sgn
                      sandfraction(ib) = sandfraction(ib)/sgn
                      siltfraction(ib) = siltfraction(ib)/sgn              
                      njsc(ib) = MAXLOC(soilclass(ib,:),1) ! Dominant texture class
                   ENDIF
                ENDIF
             ENDDO

             ! The "USDA" case reads a map of the 12 USDA texture classes,
             ! such as to assign the corresponding soil properties
          CASE("usda")
             IF (printlev_loc>=2) WRITE(numout,*) "Using a soilclass map with usda classification"

             ALLOCATE(textfrac_table(nscm,ntext), STAT=ALLOC_ERR)
             IF (ALLOC_ERR/=0) CALL ipslerr_p(3,'slowproc_soilt','Error in allocation for textfrac_table','','')

             CALL get_soilcorr_usda (nscm, textfrac_table)

             IF (printlev_loc>=4) WRITE (numout,*) 'slowproc_soilt: After get_soilcorr_usda'
             !
             DO ib =1, nbpt
                ! GO through the point we have found
                !
                ! Provide which textures were found
                ! vecpos: List of positions where textures were not zero
                !   vecpos(1): number of not null textures found
                vecpos = interpweight_ValVecR(textrefrac(ib,:),ntextinfile,zero,'neq')
                fopt = vecpos(1)
                !
                !    Check that we found some points
                !
                soilclass(ib,:) = 0.0
                clayfraction(ib) = 0.0
                sandfraction(ib) = 0.0
                siltfraction(ib) = 0.0

                IF ( fopt .EQ. 0) THEN
                   ! No points were found for current grid box, use default values
                   IF (printlev_loc>=3) WRITE(numout,*)'slowproc_soilt: no soil class in input file found for point=', ib
                   nbexp = nbexp + 1
                   njsc(ib) = usda_default ! 6 = Loam
                   clayfraction(ib) = clayfrac_usda(usda_default)
                   sandfraction(ib) = sandfrac_usda(usda_default)
                   siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)
                ELSE
                   IF (fopt == nscm) THEN
                      ! All textures are not zero
                      solt(:) = (/(i,i=1,nscm)/)
                   ELSE
                      DO ilf = 1,fopt
                         solt(ilf) = vecpos(ilf+1)
                      END DO
                   END IF

                   !   Compute the fraction of each textural class
                   DO ilf = 1,fopt
                      IF ( (solt(ilf) .LE. nscm) .AND. (solt(ilf) .GT. 0) ) THEN
                         soilclass(ib,solt(ilf)) = textrefrac(ib,solt(ilf))
                         clayfraction(ib) = clayfraction(ib) + textfrac_table(solt(ilf),3) *                &
                              textrefrac(ib,solt(ilf))
                         sandfraction(ib) = sandfraction(ib) + textfrac_table(solt(ilf),2) * &
                              textrefrac(ib,solt(ilf))
                         siltfraction(ib) = siltfraction(ib) + textfrac_table(solt(ilf),1) * &
                              textrefrac(ib,solt(ilf))
                      ELSE
                         IF (solt(ilf) .GT. nscm) THEN
                            WRITE(*,*) 'The file contains a soil color class which is incompatible with this program'
                            CALL ipslerr_p(3,'slowproc_soilt','Problem soil color class incompatible 2','','')
                         ENDIF
                      ENDIF
                      !
                   ENDDO
                   njsc(ib) = MAXLOC(soilclass(ib,:),1) ! Dominant texture class

                   ! Set default values if the surface in source file is too small
                   ! Warning - This test is donne differently for Zobler (based on sgn, related to class 6=ice)
                   IF ( atext(ib) .LT. min_sechiba) THEN
                      nbexp = nbexp + 1
                      njsc(ib) = usda_default ! 6 = Loam
                      clayfraction(ib) = clayfrac_usda(usda_default)
                      sandfraction(ib) = sandfrac_usda(usda_default)
                      siltfraction(ib) = 1.-clayfrac_usda(usda_default)-sandfrac_usda(usda_default)
                   ENDIF
                ENDIF

             ENDDO

          IF (printlev_loc>=4) WRITE (numout,*) '  slowproc_soilt: End case usda'

          CASE DEFAULT
             WRITE(numout,*) 'slowproc_soilt _______'
             WRITE(numout,*) '  A non supported soil type classification has been chosen'
             CALL ipslerr_p(3,'slowproc_soilt','non supported soil type classification','','')
          ENDSELECT
          IF (printlev_loc >= 5 ) WRITE(numout,*)'  slowproc_soilt end of type classification'

          IF ( nbexp .GT. 0 ) THEN
             WRITE(numout,*) 'slowproc_soilt:'
             WRITE(numout,*) '  The interpolation of variable soiltext had ', nbexp
             WRITE(numout,*) '  points without data. This are either coastal points or ice covered land.'
             WRITE(numout,*) '  The problem was solved by using the default soil types.'
          ENDIF

          IF (ALLOCATED(variabletypevals)) DEALLOCATE (variabletypevals)
          IF (ALLOCATED(textrefrac)) DEALLOCATE (textrefrac)
          IF (ALLOCATED(solt)) DEALLOCATE (solt)
          IF (ALLOCATED(textfrac_table)) DEALLOCATE (textfrac_table)

       ENDIF        !      xios_interpolation

       ! End of soil texture reading, for 'maps' and classical behavior
       
       IF (spmipexp == 'maps') THEN
              IF (printlev_loc>=3) WRITE (numout,*) 'slowproc_soilt: Read soil hydraulic parameters with IOIPSL'

              ! Read using IOIPSL and interpolate using aggregate tool in ORCHIDEE

              !Config Key   = PARAM_SPMIP_FILE
              !Config Desc  = Name of file from which soil parameter  values are read
              !Config Def   = params_sp_mip.nc
              !Config If    = smipexp='maps'
              !Config Help  = The name of the file to be opened to read values of parameters.
              !Config         The data from this file is then interpolated to the grid of
              !Config         of the model.
              !Config Units = [FILE]
              !
              ! params_sp_mip.nc file is 0.5 deg soil hydraulic parameters file provided by sp_mip

              filename = 'params_sp_mip.nc'
              CALL getin_p('PARAM_SPMIP_FILE',filename)

              !! Variables for interpweight
              ! Type of calculation of cell fractions
              fractype = 'default'
              ! Name of the longitude and latitude in the input file
              lonname = 'nav_lon'
              latname = 'nav_lat'
              ! Assigning values to vmin, vmax (there are not types/categories
              vmin =0.
              vmax = 99999.
              !! Variables for interpweight
              ! Should negative values be set to zero from input file?
              nonegative = .FALSE.
              ! Type of mask to apply to the input data (see header for more details)
              maskingtype = 'mabove'
              ! Values to use for the masking
              maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
              ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') ( not used)
              namemaskvar = ''

              variablename = 'ks'
              IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_soilt: Read and interpolate " &
                   // TRIM(filename) // " for variable " // TRIM(variablename)
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   ks, aparam)
              WRITE(numout,*) 'ks map is read _______'

              variablename = 'alpha'
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   avan, aparam)
              WRITE(numout,*) 'avan map read _______'

              variablename = 'thetar'
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   mcr, aparam)
              WRITE(numout,*) 'thetar map read _______'

              variablename = 'thetas'
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   mcs, aparam)
              WRITE(numout,*) 'thetas map read _______'

              variablename = 'thetapwpvg' ! mcw
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   mcw, aparam)
              WRITE(numout,*) 'thetapwpvg map read _______'

              variablename = 'thetafcvg' !mcfc
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   mcfc, aparam)
              WRITE(numout,*) 'thetafcvg map read _______'

              variablename = 'nvg'
              CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
                   contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,     &
                   maskvals, namemaskvar, -1, fractype, 0., 0.,                              &
                   nvan, aparam)
              WRITE(numout,*) 'nvan map read _______'

       ELSE ! spmipexp is not maps nor unif, then it must be texture
          IF (spmipexp == 'texture') THEN
             ! Whichever the soil texture map, we can use the USDA parameter vectors with 13 values 
             nvan(:) = nvan_usda(njsc(:))
             avan(:) = avan_usda(njsc(:))
             mcr(:) = mcr_usda(njsc(:))
             mcs(:) = mcs_usda(njsc(:))
             ks(:) = ks_usda(njsc(:))
!!$             mcfc(:) = mcf_usda(njsc(:))
!!$             mcw(:) = mcw_usda(njsc(:))
             
             !! Calculation of FC and WP based on above 5 parameters
             mvan(:) = un - (un / nvan(:))
             ! Define matrix potential in mm for wilting point and field capacity (with sand vs clay-silt variation)
             psi_w(:) = 150000.
             DO ib=1, nbpt
                IF ( ks(ib) .GE. 560 ) THEN ! Sandy soils (560 is equivalent of 2.75 at log scale of Ks, mm/d)
                   psi_fc(ib) = 1000.
                ELSE ! Finer soils
                   psi_fc(ib) = 3300. 
                ENDIF
             ENDDO
             mcfc(:) = mcr(:) + (( mcs(:) - mcr(:)) / (un + ( avan(:) * psi_fc(:))** nvan(:))** mvan(:))
             mcw(:)  = mcr(:) + (( mcs(:) - mcr(:)) / (un + ( avan(:) *  psi_w(:))** nvan(:))** mvan(:))
             
         ELSE ! if spmipexp is not among texture or maps or unif
            WRITE(numout,*) "Unsupported spmipexp=",spmipexp
            WRITE(numout,*) "Choose between texture, maps, and unif"
            CALL ipslerr_p(3,'soilproc_soilt','Bad choice of spmipexp','Choose between texture, maps, and unif','')
         ENDIF
       ENDIF
   ENDIF ! SPMIPEXP
   
    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_atext",atext)
    CALL xios_orchidee_send_field("interp_diag_soilclass",soilclass)
    CALL xios_orchidee_send_field("interp_diag_njsc",REAL(njsc, r_std))
    CALL xios_orchidee_send_field("interp_diag_clayfraction",clayfraction)
    CALL xios_orchidee_send_field("interp_diag_sandfraction",sandfraction)
    CALL xios_orchidee_send_field("interp_diag_siltfraction",siltfraction)
    
    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_soilt ended'

  END SUBROUTINE slowproc_soilt
 
!! ================================================================================================================================
!! SUBROUTINE   : slowproc_slope
!!
!>\BRIEF         Calculate mean slope coef in each  model grid box from the slope map
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::reinf_slope
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_slope(nbpt, lalo, neighbours, resolution, contfrac, reinf_slope)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)          :: nbpt                  ! Number of points for which the data needs to be interpolated
    REAL(r_std), INTENT(in)              :: lalo(nbpt,2)          ! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), INTENT(in)          :: neighbours(nbpt,NbNeighb)! Vector of neighbours for each grid point
                                                                    ! (1=North and then clockwise)
    REAL(r_std), INTENT(in)              :: resolution(nbpt,2)    ! The size in km of each grid-box in X and Y
    REAL(r_std), INTENT (in)             :: contfrac(nbpt)         !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), INTENT(out)    ::  reinf_slope(nbpt)                   ! slope coef 
    !
    !  0.3 LOCAL
    !
    !
    REAL(r_std)  :: slope_noreinf                 ! Slope above which runoff is maximum
    CHARACTER(LEN=80) :: filename
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the 
                                                                             !!   renormalization
    REAL(r_std), DIMENSION(nbpt)                         :: aslope           !! slope availability 

    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat name in the input file
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file  (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 

!_ ================================================================================================================================
    
    !
    !Config Key   = SLOPE_NOREINF
    !Config Desc  = Slope over which surface runoff does not reinfiltrate
    !Config If    = 
    !Config Def   = 0.5
    !Config Help  = The slope above which there is no reinfiltration
    !Config Units = [%]
    !
    slope_noreinf = 0.5 ! slope in percent
    !
    CALL getin_p('SLOPE_NOREINF',slope_noreinf)
    !
    !Config Key   = TOPOGRAPHY_FILE
    !Config Desc  = Name of file from which the topography map is to be read
    !Config If    = 
    !Config Def   = cartepente2d_15min.nc
    !Config Help  = The name of the file to be opened to read the orography
    !Config         map is to be given here. Usualy SECHIBA runs with a 2'
    !Config         map which is derived from the NGDC one. 
    !Config Units = [FILE]
    !
    filename = 'cartepente2d_15min.nc'
    CALL getin_p('TOPOGRAPHY_FILE',filename)

    IF (xios_interpolation) THEN
    
      CALL xios_orchidee_recv_field('reinf_slope_interp',reinf_slope)
      CALL xios_orchidee_recv_field('frac_slope_interp',aslope)


    ELSE
    
      variablename = 'pente'
      IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_slope: Read and interpolate " &
           // TRIM(filename) // " for variable " // TRIM(variablename)

      ! For this case there are not types/categories. We have 'only' a continuos field
      ! Assigning values to vmin, vmax
      vmin = 0.
      vmax = 9999.

      !! Variables for interpweight
      ! Type of calculation of cell fractions
      fractype = 'slopecalc'
      ! Name of the longitude and latitude in the input file
      lonname = 'longitude'
      latname = 'latitude'
      ! Should negative values be set to zero from input file?
      nonegative = .FALSE.
      ! Type of mask to apply to the input data (see header for more details)
      maskingtype = 'mabove'
      ! Values to use for the masking
      maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
      ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
      namemaskvar = ''

      CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
        contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
        maskvals, namemaskvar, -1, fractype, slope_default, slope_noreinf,                              &
        reinf_slope, aslope)
      IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_slope after interpweight_2Dcont'

    ENDIF
    
      ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_aslope",aslope)
    CALL xios_orchidee_send_field("interp_diag_reinf_slope",reinf_slope)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_slope ended'

  END SUBROUTINE slowproc_slope


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_woodharvest
!!
!>\BRIEF         
!!
!! DESCRIPTION  : 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): ::
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE slowproc_woodharvest(nbpt, lalo, neighbours, resolution, contfrac, woodharvest)

    USE interpweight

    IMPLICIT NONE

    !
    !
    !
    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                           :: nbpt         !! Number of points for which the data needs to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: lalo         !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in) :: neighbours   !! Vector of neighbours for each grid point
                                                                         !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)           :: resolution   !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)             :: contfrac     !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt), INTENT(out)            ::  woodharvest !! Wood harvest
    !
    !  0.3 LOCAL
    !
    CHARACTER(LEN=80)                                    :: filename
    REAL(r_std)                                          :: vmin, vmax  
    REAL(r_std), DIMENSION(nbpt)                         :: aoutvar          !! availability of input data to
                                                                             !!   interpolate output variable 
                                                                             !!   (on the nbpt space)
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat name in the input file
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds 
                                                                             !!     of something with a temporal 
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file  (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    REAL(r_std), DIMENSION(1)                            :: variabletypevals !! 
!    REAL(r_std), DIMENSION(nbp_mpi)                      :: woodharvest_mpi  !! Wood harvest where all thredds OMP are gatherd
!_ ================================================================================================================================
    
    
    !Config Key   = WOODHARVEST_FILE
    !Config Desc  = Name of file from which the wood harvest will be read
    !Config If    = DO_WOOD_HARVEST
    !Config Def   = woodharvest.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'woodharvest.nc'
    CALL getin_p('WOODHARVEST_FILE',filename)
    variablename = 'woodharvest'


    IF (xios_interpolation) THEN
       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readwoodharvest: Use XIOS to read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       CALL xios_orchidee_recv_field('woodharvest_interp',woodharvest)

       aoutvar = 1.0
    ELSE

       IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_readwoodharvest: Read and interpolate " &
            // TRIM(filename) // " for variable " // TRIM(variablename)

       ! For this case there are not types/categories. We have 'only' a continuos field
       ! Assigning values to vmin, vmax
       vmin = 0.
       vmax = 9999.
       
       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'longitude'
       latname = 'latitude'
       ! Should negative values be set to zero from input file?
       nonegative = .TRUE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'nomask'
       ! Values to use for the masking
       maskvals = (/ min_sechiba, undef_sechiba, undef_sechiba /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''
       
       variabletypevals=-un
       CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,                                &
            contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
            maskvals, namemaskvar, -1, fractype, 0., 0., woodharvest, aoutvar)
       IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_wodharvest after interpweight_2Dcont'
       
    END IF

    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_diag_woodharvest",woodharvest)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_woodharvest ended'
  END SUBROUTINE slowproc_woodharvest


!! ================================================================================================================================
!! SUBROUTINE 	: get_soilcorr_zobler
!!
!>\BRIEF         The "get_soilcorr" routine defines the table of correspondence
!!               between the Zobler types and the three texture types known by SECHIBA and STOMATE :
!!               silt, sand and clay. 
!!
!! DESCRIPTION : get_soilcorr is needed if you use soils_param.nc .\n
!!               The data from this file is then interpolated to the grid of the model. \n
!!               The aim is to get fractions for sand loam and clay in each grid box.\n
!!               This information is used for soil hydrology and respiration.
!!
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S) : ::texfrac_table
!!
!! REFERENCE(S)	: 
!! - Zobler L., 1986, A World Soil File for global climate modelling. NASA Technical memorandum 87802. NASA 
!!   Goddard Institute for Space Studies, New York, U.S.A.
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_soilcorr_zobler (nzobler,textfrac_table)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    INTEGER(i_std),PARAMETER :: nbtypes_zobler = 7                    !! Number of Zobler types (unitless)

    !! 0.1  Input variables
    
    INTEGER(i_std),INTENT(in) :: nzobler                              !! Size of the array (unitless)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION(nzobler,ntext),INTENT(out) :: textfrac_table !! Table of correspondence between soil texture class
                                                                       !! and granulometric composition (0-1, unitless)
    
    !! 0.4 Local variables
    
    INTEGER(i_std) :: ib                                              !! Indice (unitless)
    
!_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-  
    IF (nzobler /= nbtypes_zobler) THEN 
       CALL ipslerr_p(3,'get_soilcorr', 'nzobler /= nbtypes_zobler',&
          &   'We do not have the correct number of classes', &
          &                 ' in the code for the file.')  ! Fatal error
    ENDIF

    !-
    ! 1. Textural fraction for : silt        sand         clay
    !-
    textfrac_table(1,:) = (/ 0.12, 0.82, 0.06 /)
    textfrac_table(2,:) = (/ 0.32, 0.58, 0.10 /)
    textfrac_table(3,:) = (/ 0.39, 0.43, 0.18 /)
    textfrac_table(4,:) = (/ 0.15, 0.58, 0.27 /)
    textfrac_table(5,:) = (/ 0.34, 0.32, 0.34 /)
    textfrac_table(6,:) = (/ 0.00, 1.00, 0.00 /)
    textfrac_table(7,:) = (/ 0.39, 0.43, 0.18 /)


    !-
    ! 2. Check the mapping for the Zobler types which are going into the ORCHIDEE textures classes 
    !-
    DO ib=1,nzobler ! Loop over # classes soil
       
       IF (ABS(SUM(textfrac_table(ib,:))-1.0) > EPSILON(1.0)) THEN ! The sum of the textural fractions should not exceed 1 !
          WRITE(numout,*) &
               &     'Error in the correspondence table', &
               &     ' sum is not equal to 1 in', ib
          WRITE(numout,*) textfrac_table(ib,:)
          CALL ipslerr_p(3,'get_soilcorr', 'SUM(textfrac_table(ib,:)) /= 1.0',&
               &                 '', 'Error in the correspondence table') ! Fatal error
       ENDIF
       
    ENDDO ! Loop over # classes soil

    
  END SUBROUTINE get_soilcorr_zobler

!! ================================================================================================================================
!! SUBROUTINE 	: get_soilcorr_usda
!!
!>\BRIEF         The "get_soilcorr_usda" routine defines the table of correspondence
!!               between the 12 USDA textural classes and their granulometric composition, 
!!               as % of silt, sand and clay. This is used to further defien clayfraction.
!!
!! DESCRIPTION : get_soilcorr is needed if you use soils_param.nc .\n
!!               The data from this file is then interpolated to the grid of the model. \n
!!               The aim is to get fractions for sand loam and clay in each grid box.\n
!!               This information is used for soil hydrology and respiration.
!!               The default map in this case is derived from Reynolds et al 2000, \n
!!               at the 1/12deg resolution, with indices that are consistent with the \n
!!               textures tabulated below
!!
!! RECENT CHANGE(S): Created by A. Ducharne on July 02, 2014
!!
!! MAIN OUTPUT VARIABLE(S) : ::texfrac_table
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART	: None
!! \n
!_ ================================================================================================================================

  SUBROUTINE get_soilcorr_usda (nusda,textfrac_table)

    IMPLICIT NONE

    !! 0. Variables and parameters declaration
    
    !! 0.1  Input variables
    
    INTEGER(i_std),INTENT(in) :: nusda                               !! Size of the array (unitless)
    
    !! 0.2 Output variables 
    
    REAL(r_std),DIMENSION(nusda,ntext),INTENT(out) :: textfrac_table !! Table of correspondence between soil texture class
                                                                     !! and granulometric composition (0-1, unitless)
    
    !! 0.4 Local variables

    INTEGER(i_std),PARAMETER :: nbtypes_usda = 13                    !! Number of USDA texture classes (unitless)
    INTEGER(i_std) :: n                                              !! Index (unitless)
    
!_ ================================================================================================================================

    !-
    ! 0. Check consistency
    !-  
    IF (nusda /= nbtypes_usda) THEN 
       CALL ipslerr_p(3,'get_soilcorr', 'nusda /= nbtypes_usda',&
          &   'We do not have the correct number of classes', &
          &                 ' in the code for the file.')  ! Fatal error
    ENDIF

    !! Parameters for soil type distribution :
    !! Sand, Loamy Sand, Sandy Loam, Silt Loam, Silt, Loam, Sandy Clay Loam, Silty Clay Loam, Clay Loam, Sandy Clay, Silty Clay, Clay
    ! The order comes from constantes_soil.f90
    ! The corresponding granulometric composition comes from Carsel & Parrish, 1988

    !-
    ! 1. Textural fractions for : sand, clay
    !-
    textfrac_table(1,2:3)  = (/ 0.93, 0.03 /) ! Sand
    textfrac_table(2,2:3)  = (/ 0.81, 0.06 /) ! Loamy Sand
    textfrac_table(3,2:3)  = (/ 0.63, 0.11 /) ! Sandy Loam
    textfrac_table(4,2:3)  = (/ 0.17, 0.19 /) ! Silt Loam
    textfrac_table(5,2:3)  = (/ 0.06, 0.10 /) ! Silt
    textfrac_table(6,2:3)  = (/ 0.40, 0.20 /) ! Loam
    textfrac_table(7,2:3)  = (/ 0.54, 0.27 /) ! Sandy Clay Loam
    textfrac_table(8,2:3)  = (/ 0.08, 0.33 /) ! Silty Clay Loam
    textfrac_table(9,2:3)  = (/ 0.30, 0.33 /) ! Clay Loam
    textfrac_table(10,2:3) = (/ 0.48, 0.41 /) ! Sandy Clay
    textfrac_table(11,2:3) = (/ 0.06, 0.46 /) ! Silty Clay
    textfrac_table(12,2:3) = (/ 0.15, 0.55 /) ! Clay
    textfrac_table(13,2:3) = (/ 0.15, 0.55 /) ! Clay

    ! Fraction of silt

    DO n=1,nusda
       textfrac_table(n,1) = 1. - textfrac_table(n,2) - textfrac_table(n,3)
    END DO
       
  END SUBROUTINE get_soilcorr_usda

!! ================================================================================================================================
!! FUNCTION 	: tempfunc
!!
!>\BRIEF        ! This function interpolates value between ztempmin and ztempmax
!! used for lai detection. 
!!
!! DESCRIPTION   : This subroutine calculates a scalar between 0 and 1 with the following equation :\n
!!                 \latexonly
!!                 \input{constantes_veg_tempfunc.tex}
!!                 \endlatexonly
!!
!! RECENT CHANGE(S): None
!!
!! RETURN VALUE : tempfunc_result
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  FUNCTION tempfunc (temp_in) RESULT (tempfunc_result)


    !! 0. Variables and parameters declaration

    REAL(r_std),PARAMETER    :: ztempmin=273._r_std   !! Temperature for laimin (K)
    REAL(r_std),PARAMETER    :: ztempmax=293._r_std   !! Temperature for laimax (K)
    REAL(r_std)              :: zfacteur              !! Interpolation factor   (K^{-2})

    !! 0.1 Input variables

    REAL(r_std),INTENT(in)   :: temp_in               !! Temperature (K)

    !! 0.2 Result

    REAL(r_std)              :: tempfunc_result       !! (unitless)
    
!_ ================================================================================================================================

    !! 1. Define a coefficient
    zfacteur = un/(ztempmax-ztempmin)**2
    
    !! 2. Computes tempfunc
    IF     (temp_in > ztempmax) THEN
       tempfunc_result = un
    ELSEIF (temp_in < ztempmin) THEN
       tempfunc_result = zero
    ELSE
       tempfunc_result = un-zfacteur*(ztempmax-temp_in)**2
    ENDIF !(temp_in > ztempmax)


  END FUNCTION tempfunc


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_checkveget
!!
!>\BRIEF         To verify the consistency of the various fractions defined within the grid box after having been
!!               been updated by STOMATE or the standard procedures.
!!
!! DESCRIPTION  : (definitions, functional, design, flags): 
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: none
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
!
  SUBROUTINE slowproc_checkveget(nbpt, frac_nobio, veget_max, veget, tot_bare_soil, soiltile)

    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                      :: nbpt       ! Number of points for which the data needs to be interpolated
    REAL(r_std),DIMENSION (nbpt,nnobio), INTENT(in) :: frac_nobio ! Fraction of ice,lakes,cities, ... (unitless)
    REAL(r_std),DIMENSION (nbpt,nvm), INTENT(in)    :: veget_max  ! Maximum fraction of vegetation type including none biological fraction (unitless) 
    REAL(r_std),DIMENSION (nbpt,nvm), INTENT(in)    :: veget      ! Vegetation fractions
    REAL(r_std),DIMENSION (nbpt), INTENT(in)        :: tot_bare_soil ! Total evaporating bare soil fraction within the mesh
    REAL(r_std),DIMENSION (nbpt,nstm), INTENT(in)   :: soiltile   ! Fraction of soil tiles in the gridbox (unitless)

    !  0.3 LOCAL
    !
    INTEGER(i_std) :: ji, jn, jv
    REAL(r_std)  :: epsilocal  !! A very small value
    REAL(r_std)  :: totfrac
    CHARACTER(len=80) :: str1, str2
    
!_ ================================================================================================================================
    
    !
    ! There is some margin added as the computing errors might bring us above EPSILON(un)
    !
    epsilocal = EPSILON(un)*1000.
    
    !! 1.0 Verify that none of the fractions are smaller than min_vegfrac, without beeing zero.
    !!
    DO ji=1,nbpt
       DO jn=1,nnobio
          IF ( frac_nobio(ji,jn) > epsilocal .AND. frac_nobio(ji,jn) < min_vegfrac ) THEN
             WRITE(str1,'("Occurs on grid box", I8," and nobio type ",I3 )') ji, jn
             WRITE(str2,'("The small value obtained is ", E14.4)') frac_nobio(ji,jn)
             CALL ipslerr_p (3,'slowproc_checkveget', &
                  "frac_nobio is larger than zero but smaller than min_vegfrac.", str1, str2)
          ENDIF
       ENDDO
    END DO
    
    IF (.NOT. ok_dgvm) THEN       
       DO ji=1,nbpt
          DO jv=1,nvm
             IF ( veget_max(ji,jv) > epsilocal .AND. veget_max(ji,jv) < min_vegfrac ) THEN
                WRITE(str1,'("Occurs on grid box", I8," and nobio type ",I3 )') ji, jn
                WRITE(str2,'("The small value obtained is ", E14.4)') veget_max(ji,jv)
                CALL ipslerr_p (3,'slowproc_checkveget', &
                     "veget_max is larger than zero but smaller than min_vegfrac.", str1, str2)
             ENDIF
          ENDDO
       ENDDO
    END IF
    
    !! 2.0 verify that with all the fractions we cover the entire grid box  
    !!
    DO ji=1,nbpt
       totfrac = zero
       DO jn=1,nnobio
          totfrac = totfrac + frac_nobio(ji,jn)
       ENDDO
       DO jv=1,nvm
          totfrac = totfrac + veget_max(ji,jv)
       ENDDO
       IF ( ABS(totfrac - un) > epsilocal) THEN
             WRITE(str1,'("This occurs on grid box", I8)') ji
             WRITE(str2,'("The sum over all fraction and error are ", E14.4, E14.4)') totfrac, ABS(totfrac - un)
             CALL ipslerr_p (3,'slowproc_checkveget', &
                   "veget_max + frac_nobio is not equal to 1.", str1, str2)
             WRITE(*,*) "EPSILON =", epsilocal 
       ENDIF
    ENDDO
    
    !! 3.0 Verify that veget is smaller or equal to veget_max
    !!
    DO ji=1,nbpt
       DO jv=1,nvm
          IF ( jv == ibare_sechiba ) THEN
             IF ( ABS(veget(ji,jv) - veget_max(ji,jv)) > epsilocal ) THEN
                WRITE(str1,'("This occurs on grid box", I8)') ji
                WRITE(str2,'("The difference is ", E14.4)') veget(ji,jv) - veget_max(ji,jv)
                CALL ipslerr_p (3,'slowproc_checkveget', &
                     "veget is not equal to veget_max on bare soil.", str1, str2)
             ENDIF
          ELSE
             IF ( veget(ji,jv) > veget_max(ji,jv) ) THEN
                WRITE(str1,'("This occurs on grid box", I8)') ji
                WRITE(str2,'("The values for veget and veget_max :", F8.4, F8.4)') veget(ji,jv), veget_max(ji,jv)
                CALL ipslerr_p (3,'slowproc_checkveget', &
                     "veget is greater than veget_max.", str1, str2)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    
    !! 4.0 Test tot_bare_soil in relation to the other variables
    !!
    DO ji=1,nbpt
       totfrac = zero
       DO jv=1,nvm
          totfrac = totfrac + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
       ! add the bare soil fraction to totfrac
       totfrac = totfrac + veget(ji,ibare_sechiba)
       ! do the test
       IF ( ABS(totfrac - tot_bare_soil(ji)) > epsilocal ) THEN
          WRITE(str1,'("This occurs on grid box", I8)') ji
          WRITE(str2,'("The values for tot_bare_soil, tot frac and error :", F8.4, F8.4, E14.4)') &
               &  tot_bare_soil(ji), totfrac, ABS(totfrac - tot_bare_soil(ji))
          CALL ipslerr_p (3,'slowproc_checkveget', &
               "tot_bare_soil does not correspond to the total bare soil fraction.", str1, str2)
       ENDIF
    ENDDO
    
    !! 5.0 Test that soiltile has the right sum
    !!
    DO ji=1,nbpt
       totfrac = SUM(soiltile(ji,:))
       IF ( ABS(totfrac - un) > epsilocal ) THEN
          WRITE(numout,*) "soiltile does not sum-up to one. This occurs on grid box", ji
          WRITE(numout,*) "The soiltile for ji are :", soiltile(ji,:)
          CALL ipslerr_p (2,'slowproc_checkveget', &
               "soiltile does not sum-up to one.", "", "")
       ENDIF
    ENDDO
    
  END SUBROUTINE slowproc_checkveget


!! ================================================================================================================================
!! SUBROUTINE   : slowproc_change_frac
!!
!>\BRIEF        Update the vegetation fractions
!!
!! DESCRIPTION  : Update the vegetation fractions. This subroutine is called in the same time step as lcchange in stomatelpj has
!!                has been done. This subroutine is called after the diagnostics have been written in sechiba_main.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): :: veget_max, veget, frac_nobio, totfrac_nobio, tot_bare_soil, soiltile
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================
   
  SUBROUTINE slowproc_change_frac(kjpindex, lai, &
                                  veget_max, veget, frac_nobio, totfrac_nobio, tot_bare_soil, soiltile, fraclut, nwdFraclut)
    !
    ! 0. Declarations
    !
    ! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                           :: kjpindex       !! Domain size - terrestrial pixels only
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(in)     :: lai            !! Leaf area index (m^2 m^{-2})
    
    ! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)    :: veget_max      !! Maximum fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT(out)    :: veget          !! Fraction of vegetation type in the mesh (unitless)
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(out) :: frac_nobio     !! Fraction of ice, lakes, cities etc. in the mesh
    REAL(r_std),DIMENSION (kjpindex), INTENT(out)        :: totfrac_nobio  !! Total fraction of ice+lakes+cities etc. in the mesh
    REAL(r_std), DIMENSION (kjpindex), INTENT(out)       :: tot_bare_soil  !! Total evaporating bare soil fraction in the mesh
    REAL(r_std), DIMENSION (kjpindex,nstm), INTENT(out)  :: soiltile       !! Fraction of each soil tile within vegtot (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)  :: fraclut        !! Fraction of each landuse tile (0-1, unitless)
    REAL(r_std), DIMENSION (kjpindex,nlut), INTENT(out)  :: nwdfraclut     !! Fraction of non woody vegetation in each landuse tile (0-1, unitless)
    
    ! 0.3 Local variables
    INTEGER(i_std)                                       :: ji, jv         !! Loop index
    
       
    !! Update vegetation fractions with the values coming from the vegetation file read in slowproc_readvegetmax.
    !! Partial update has been taken into account for the case with DGVM and AGRICULTURE in slowproc_readvegetmax.
    veget_max  = veget_max_new
    frac_nobio = frac_nobio_new
       
    !! Verification and correction on veget_max, calculation of veget and soiltile.
    CALL slowproc_veget (kjpindex, lai, frac_nobio, totfrac_nobio, veget_max, veget, soiltile, fraclut, nwdFraclut)
    
    !! Calculate tot_bare_soil needed in hydrol, diffuco and condveg (fraction of bare soil in the mesh)
    tot_bare_soil(:) = veget_max(:,1)
    DO jv = 2, nvm
       DO ji =1, kjpindex
          tot_bare_soil(ji) = tot_bare_soil(ji) + (veget_max(ji,jv) - veget(ji,jv))
       ENDDO
    END DO

    !! Do some basic tests on the surface fractions updated above 
    CALL slowproc_checkveget(kjpindex, frac_nobio, veget_max, veget, tot_bare_soil, soiltile)
     
  END SUBROUTINE slowproc_change_frac 

  !! ================================================================================================================================
  !! SUBROUTINE   : slowproc_readirrigmap_dyn
  !!
  !>\BRIEF        Function to interpolate irrigation maps
  !!
  !! DESCRIPTION  : This function interpolates the irrigation maps from original resolution to simul. resolution
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): :: irrigmap_new
  !!
  !! REFERENCE(S) : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================

  SUBROUTINE slowproc_readirrigmap_dyn(nbpt, lalo, neighbours,  resolution, contfrac,         &
       irrigmap_new)

    USE interpweight

    IMPLICIT NONE

    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt), INTENT(out)          :: irrigmap_new       !! new irrigation map in m2 per grid cell
    !
    !  0.3 LOCAL
    !
    !
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib
    !
    ! for irrigated_new case :

    REAL(r_std), DIMENSION(nbpt)                     :: irrigref_frac    !! irrigation fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: airrig           !! Availability of the soilcol interpolation
    REAL(r_std)                          :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds
                                                                             !!     of something with a temporal
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to
                                                                             !!   `maskingtype')
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask
    CHARACTER(LEN=250)                                   :: msg

  !_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_read irrigmap_new'

    !
    !Config Key   = IRRIGATION_DYN_FILE
    !Config Desc  = Name of file from which the DYNAMIC irrigation fraction map is to be read
    !Config If    = IRRIG_DYN
    !Config Def   = IRRIGmap.nc
    !Config Help  = The name of the file to be opened to read an irrigation
    !Config         map is to be given here.
    !Config Units = [FILE]
    !
    filename = 'IRRIGmap.nc'
    CALL getin_p('IRRIGATION_DYN_FILE',filename)
    variablename = 'irrig'

    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_read irrigmap_new: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    ! Assigning values to vmin, vmax
    vmin = 1.
    vmax = 1.

    variabletypevals = -un

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'lon'
    latname = 'lat'
    ! Should negative values be set to zero from input file?
    nonegative = .TRUE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'nomask'
    ! Values to use for the masking
    maskvals = (/ 0.05, 0.05 , un /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,        &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, -1, fractype, 0., 0.,                                 &
      irrigref_frac, airrig)

    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_read after interpweight_2Dcont'

    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_read irrigmap_new before updating loop nbpt:', nbpt
    END IF
    irrigmap_new(:) = zero
    irrigref_frac(:) = irrigref_frac(:)/100.

    DO ib=1,nbpt
      IF (irrigref_frac(ib) < 1. .AND. irrigref_frac(ib) > 0.005 ) THEN
        irrigmap_new(ib) = irrigref_frac(ib) * resolution(ib,1)*resolution(ib,2)*contfrac(ib)
      ELSEIF (irrigref_frac(ib) > 1.) THEN
        irrigmap_new(ib) = resolution(ib,1)*resolution(ib,2)*contfrac(ib)
      !ELSE  THEN irrigated_new= zero. ALready set, put to lisibility
      ENDIF
    ENDDO

    ! Write diagnostics
    !CALL xios_orchidee_send_field("airrig",airrig)

    IF (printlev_loc >= 3) WRITE(numout,*) '  slowproc_read irrigmap_new ended'

  END SUBROUTINE slowproc_readirrigmap_dyn

  !! ================================================================================================================================
  !! SUBROUTINE   : slowproc_read_aeisw_map
  !!
  !>\BRIEF        Function to interpolate irrigation maps
  !!
  !! DESCRIPTION  : This function interpolates the maps of areas equipped for irrigation with surface water
  !!                from original resolution to simul. resolution. This is used in the new irrigation scheme, that
  !!                restrain water availability according to environmental needs and type of equippment to irrigate.
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): :: fraction_aeirrig_sw
  !!
  !! REFERENCE(S) : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================


  SUBROUTINE slowproc_read_aeisw_map(nbpt, lalo, neighbours,  resolution, contfrac,         &
     fraction_aeirrig_sw)

    USE interpweight

    IMPLICIT NONE

    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
    REAL(r_std), DIMENSION(nbpt), INTENT(out)          :: fraction_aeirrig_sw       !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                                ! 1.0 here corresponds to fraction of irrig. area, not grid cell
    !
    !  0.3 LOCAL
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib
    !
    ! for irrigated_new case :

    REAL(r_std), DIMENSION(nbpt)                     :: irrigref_frac    !! irrigation fractions re-dimensioned
    REAL(r_std), DIMENSION(nbpt)                         :: airrig           !! Availability of the soilcol interpolation
    REAL(r_std)                          :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds
                                                                             !!     of something with a temporal
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to
                                                                             !!   `maskingtype')
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask
    CHARACTER(LEN=250)                                   :: msg

  !_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_read irrigmap_new'

    !
    !Config Key   = FRACTION_AEI_SW_FILE
    !Config Desc  = Name of file with AEI with SW
    !Config If    = SELECT_SOURCE_IRRIG
    !Config Def   = AEI_SW_pct.nc
    !Config Help  = The name of the file to be opened to read an AIE from SW
    !Config         map is to be given here.
    !Config Units = [FILE]
    !
    filename = 'AEI_SW_pct.nc'
    CALL getin_p('FRACTION_AEI_SW_FILE',filename)
    variablename = 'aeisw_pct'

    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_read_aeisw_map: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    ! Assigning values to vmin, vmax
    vmin = 1.
    vmax = 1.

    variabletypevals = -un

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'lon'
    latname = 'lat'
    ! Should negative values be set to zero from input file?
    nonegative = .TRUE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'mabove'
    ! Values to use for the masking
    maskvals = (/ 0.05, 0.05 , un /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,        &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, -1, fractype, 0., 0.,                                 &
      irrigref_frac, airrig)

    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_read after interpweight_2Dcont'

    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_read_aeisw_map before updating loop nbpt:', nbpt
    END IF
    fraction_aeirrig_sw(:) = irrigref_frac(:)/100.
    ! Write diagnostics
    !CALL xios_orchidee_send_field("airrig",airrig)
    DO ib=1,nbpt

      fraction_aeirrig_sw(ib) = MIN(fraction_aeirrig_sw(ib), 0.99 )
      fraction_aeirrig_sw(ib) = MAX(fraction_aeirrig_sw(ib), 0.01 )

    ENDDO


    IF (printlev_loc >= 3) WRITE(numout,*) 'slowproc_read_aeisw_map ended'

  END SUBROUTINE slowproc_read_aeisw_map

  !! ================================================================================================================================
  !! SUBROUTINE   : slowproc_imperviousness
  !!
  !>\BRIEF        Function to read and interpolate imperviousness maps
  !!
  !! DESCRIPTION  : Function to read and interpolate imperviousness maps
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): :: frac_imperv
  !!
  !! REFERENCE(S) : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================


  SUBROUTINE slowproc_imperviousness(nbpt, lalo, neighbours,  resolution, contfrac)

    USE interpweight

    IMPLICIT NONE

    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
                                                                                ! 
    !
    !  0.3 LOCAL
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib
    REAL(r_std), DIMENSION(nbpt)                         :: imperviousness       !! fraction of imperviousness from map
    REAL(r_std), DIMENSION(nbpt)                         :: aimperviousness           !! Availability of the imperviousness interpolation
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds
                                                                             !!     of something with a temporal
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to
                                                                             !!   `maskingtype')
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask
    CHARACTER(LEN=250)                                   :: msg

  !_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_read imperviousness'

    !
    !Config Key   = FRACTION_IMPERVIOUSNESS
    !Config Desc  = Name of file with frac_imperv
    !Config If    = DO_IMPERVIOUSNESS
    !Config Def   = imperviousness.nc
    !Config Help  = The name of the file to be opened to read an frac_imperv
    !Config         map is to be given here.
    !Config Units = [FILE]
    !
    filename = 'imperviousness.nc'
    CALL getin_p('FRACTION_IMPERVIOUSNESS',filename)
    variablename = 'Imperviousness' 

    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_imperviousness: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    ! Assigning values to vmin, vmax
    vmin = 1.
    vmax = 1.

    variabletypevals = -un

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'lon'
    latname = 'lat'
    ! Should negative values be set to zero from input file?
    nonegative = .TRUE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'mabove'
    ! Values to use for the masking
    maskvals = (/ 0.05, 0.05 , un /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,        &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, -1, fractype, 0., 0.,                                 &
      imperviousness, aimperviousness)

    IF (printlev_loc >= 5) WRITE(numout,*)'  slowproc_read after interpweight_2Dcont'

    IF (printlev_loc >= 5) THEN
      WRITE(numout,*)'  slowproc_imperviousness before updating loop nbpt:', nbpt
    END IF
    frac_imperv(:) = imperviousness(:)/100.
    
    ! Write diagnostics
    !CALL xios_orchidee_send_field("aimperviousness",aimperviousness)

    IF (printlev_loc >= 3) WRITE(numout,*) 'slowproc_imperviousness ended'
    DO ib=1,nbpt
      frac_imperv(ib) = MIN(frac_imperv(ib), 0.99 )
      frac_imperv(ib) = MAX(frac_imperv(ib), 0.01 )
    ENDDO



  END SUBROUTINE slowproc_imperviousness 

  !! ================================================================================================================================
  !! SUBROUTINE   : slowproc_height_buidings
  !!
  !>\BRIEF        Function to read and interpolate imperviousness maps
  !!
  !! DESCRIPTION  : Function to read and interpolate imperviousness maps
  !!
  !! RECENT CHANGE(S): None
  !!
  !! MAIN OUTPUT VARIABLE(S): :: frac_imperv
  !!
  !! REFERENCE(S) : None
  !!
  !! FLOWCHART    : None
  !! \n
  !_ ================================================================================================================================


  SUBROUTINE slowproc_height_buidings(nbpt, lalo, neighbours,  resolution, contfrac)

    USE interpweight

    IMPLICIT NONE

    !  0.1 INPUT
    !
    INTEGER(i_std), INTENT(in)                             :: nbpt            !! Number of points for which the data needs
                                                                              !! to be interpolated
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: lalo            !! Vector of latitude and longitudes (beware of the order !)
    INTEGER(i_std), DIMENSION(nbpt,NbNeighb), INTENT(in)   :: neighbours      !! Vector of neighbours for each grid point
                                                                              !! (1=North and then clockwise)
    REAL(r_std), DIMENSION(nbpt,2), INTENT(in)             :: resolution      !! The size in km of each grid-box in X and Y
    REAL(r_std), DIMENSION(nbpt), INTENT(in)               :: contfrac        !! Fraction of continent in the grid
    !
    !  0.2 OUTPUT
    !
                                                                                ! 
    !
    !  0.3 LOCAL
    CHARACTER(LEN=80) :: filename
    INTEGER(i_std) :: ib
    REAL(r_std), DIMENSION(nbpt), INTENT(out)            :: height_buidings       !! fraction of imperviousness from map
    REAL(r_std), DIMENSION(nbpt)                         :: aheight_buidings           !! Availability of the imperviousness interpolation
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the renormalization
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
    REAL(r_std), DIMENSION(nvm)                          :: variabletypevals !! Values for all the types of the variable
                                                                             !!   (variabletypevals(1) = -un, not used)
    CHARACTER(LEN=50)                                    :: fractype         !! method of calculation of fraction
                                                                             !!   'XYKindTime': Input values are kinds
                                                                             !!     of something with a temporal
                                                                             !!     evolution on the dx*dy matrix'
    LOGICAL                                              :: nonegative       !! whether negative values should be removed
    CHARACTER(LEN=50)                                    :: maskingtype      !! Type of masking
                                                                             !!   'nomask': no-mask is applied
                                                                             !!   'mbelow': take values below maskvals(1)
                                                                             !!   'mabove': take values above maskvals(1)
                                                                             !!   'msumrange': take values within 2 ranges;
                                                                             !!      maskvals(2) <= SUM(vals(k)) <= maskvals(1)
                                                                             !!      maskvals(1) < SUM(vals(k)) <= maskvals(3)
                                                                             !!        (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to
                                                                             !!   `maskingtype')
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask
    CHARACTER(LEN=250)                                   :: msg

  !_ ================================================================================================================================

    IF (printlev_loc >= 5) PRINT *,'  In slowproc_read HeightBuild'

    !
    !Config Key   = BUILDING_HEIGHT
    !Config Desc  = Name of file with HeightBuild
    !Config If    = DO_HEIGHT_BUILDING
    !Config Def   = Height_buildings.nc
    !Config Help  = The name of the file to be opened to read an HeightBuild
    !Config         map is to be given here.
    !Config Units = [FILE]
    !
    filename = 'Height_buildings.nc'
    CALL getin_p('BUILDING_HEIGHT',filename)
    variablename = 'HeightBuild' 

    IF (printlev_loc >= 1) WRITE(numout,*) "slowproc_height_buidings: Read and interpolate " &
         // TRIM(filename) // " for variable " // TRIM(variablename)

    ! Assigning values to vmin, vmax
    vmin = 1.
    vmax = 1.

    variabletypevals = -un

    !! Variables for interpweight
    ! Type of calculation of cell fractions
    fractype = 'default'
    ! Name of the longitude and latitude in the input file
    lonname = 'lon'
    latname = 'lat'
    ! Should negative values be set to zero from input file?
    nonegative = .TRUE.
    ! Type of mask to apply to the input data (see header for more details)
    maskingtype = 'mabove'
    ! Values to use for the masking
    maskvals = (/ 0.05, 0.05 , un /)
    ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
    namemaskvar = ''

    CALL interpweight_2Dcont(nbpt, 0, 0, lalo, resolution, neighbours,        &
      contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
      maskvals, namemaskvar, -1, fractype, 0., 0.,                                 &
      height_buidings, aheight_buidings)

    IF (printlev_loc >= 3) WRITE(numout,*) 'slowproc_height_buidings ended'
    DO ib=1,nbpt
      height_buidings(ib) = MIN(height_buidings(ib), 0.99 )
      height_buidings(ib) = MAX(height_buidings(ib), 0.01 )
    ENDDO



  END SUBROUTINE slowproc_height_buidings 
















END MODULE slowproc
