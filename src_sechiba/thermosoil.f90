! =================================================================================================================================
! MODULE       : thermosoil
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Calculates the soil temperatures by solving the heat
!! diffusion equation within the soil. This module is only used with CWRR hydrology.
!!
!!\n DESCRIPTION : General important informations about the numerical scheme and
!!                 the soil vertical discretization:\n
!!               - the soil is zmaxt deep (by default 10m) and divided into "ngrnd" layers. 
!!                 From 0-zmaxh(default 2m), the discretization is the same as for hydrology. 
!!                 From zmaxh(2m) and below, the depth increase linearly (by default) or geometrically. \n
!!               - "jg" is usually used as the index going from 1 to ngrnd to describe the
!!                  layers, from top (jg=1) to bottom (jg=ngrnd)\n
!!               - the thermal numerical scheme is implicit finite differences.\n
!!                 -- When it is resolved in thermosoil_profile at the present timestep t, the
!!                 dependancy from the previous timestep (t-1) is hidden in the
!!                 integration coefficients cgrnd and dgrnd, which are therefore
!!                 calculated at the very end of thermosoil_main (call to
!!                 thermosoil_coef) for use in the next timestep.\n
!!                 -- At timestep t, the system becomes :\n 
!!
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!
!!                 (the bottom boundary condition has been used to obtained this equation).\n
!!                 To solve it, the uppermost soil temperature T(1) is required.
!!                 It is obtained from the surface temperature Ts, which is
!!                 considered a linear extrapolation of T(1) and T(2)\n
!!
!!                           Ts=(1+lambda)*T(1) -lambda*T(2) \n 
!!                                      -- EQ2--\n
!!
!!                 -- caveat 1 : Ts is called 'temp_soil_new' in this routine,
!!                 don' t act.\n
!!                 -- caveat 2 : actually, the surface temperature at time t Ts
!!                 depends on the soil temperature at time t through the
!!                 ground heat flux. This is again implicitly solved, with Ts(t)
!!                 expressed as :\n
!!
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflx+otherfluxes(Ts(t))\n 
!!                                      -- EQ3 --\n
!!
!!                 and the dependency from the previous timestep is hidden in
!!                 soilcap and soilflx (apparent surface heat capacity and heat
!!                 flux respectively). Soilcap and soilflx are therefore
!!                 calculated at the previous timestep, at the very end of thermosoil
!!                 (final call to thermosoil_coef) and stored to be used at the next time step.
!!                 At timestep t, EQ3 is solved for Ts in enerbil, and Ts
!!                 is used in thermosoil to get T(1) and solve EQ1.\n
!!
!! - lambda is the @tex $\mu$ @endtex of F. Hourdin' s PhD thesis, equation (A28); ie the
!! coefficient of the linear extrapolation of Ts (surface temperature) from T1 and T2 (ptn(jg=1) and ptn(jg=2)), so that:\n
!! Ts= (1+lambda)*T(1)-lambda*T(2) --EQ2-- \n
!! lambda = (zlt(1))/((zlt(2)-zlt(1))) \n
!!
!! RECENT CHANGE(S) : - Change soil thermal properties to consider also soil texture, rev 2922.
!!                    - Change vertical discretization, rev 2917. Note: In the revised thermosoil, 
!!                    cstgrnd and lskin are not needed any more. The depth znt, zlt and dlt
!!                    are computed in vertical_soil and are in meter
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE_2_2/ORCHIDEE/src_sechiba/thermosoil.f90 $
!! $Date: 2022-07-20 13:09:05 +0200 (Wed, 20 Jul 2022) $
!! $Revision: 7710 $
!! \n
!_ ================================================================================================================================

MODULE thermosoil

  ! modules used :
  USE ioipsl
  USE ioipsl_para
  USE xios_orchidee
  USE constantes
  USE time, ONLY : one_day, dt_sechiba
  USE constantes_soil
  USE sechiba_io_p
  USE grid

  IMPLICIT NONE

  !private and public routines :
  PRIVATE
  PUBLIC :: thermosoil_main, thermosoil_clear,  thermosoil_initialize, thermosoil_finalize, thermosoil_xios_initialize

  REAL(r_std), SAVE                               :: lambda                   !! See Module description
!$OMP THREADPRIVATE(lambda)
  REAL(r_std), SAVE                               :: fz1, zalph               !! usefull constants for diverse use
!$OMP THREADPRIVATE(fz1, zalph)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: ptn                      !! vertically discretized 
                                                                              !! soil temperatures @tex ($K$) @endtex. 
!$OMP THREADPRIVATE(ptn)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: dz1                      !! numerical constant used in the thermal numerical
                                                                              !! scheme  @tex ($m^{-1}$) @endtex. ; it corresponds
                                                                              !! to the coefficient  @tex $d_k$ @endtex of equation
                                                                              !! (A.12) in F. Hourdin PhD thesis.
!$OMP THREADPRIVATE(dz1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: cgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
!$OMP THREADPRIVATE(cgrnd)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: dgrnd                    !! integration coefficient for the numerical scheme,
                                                                              !! see eq.1
!$OMP THREADPRIVATE(dgrnd)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa                    !! volumetric vertically discretized soil heat 
                                                                              !! capacity  @tex ($J K^{-1} m^{-3}$) @endtex. 
                                                                              !! It depends on the soil
                                                                              !! moisture content (shum_ngrnd_perma) and is calculated at 
                                                                              !! each time step in thermosoil_coef.
!$OMP THREADPRIVATE(pcapa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa                   !! vertically discretized soil thermal conductivity 
                                                                              !!  @tex ($W K^{-1} m^{-1}$) @endtex. Same as pcapa.
!$OMP THREADPRIVATE(pkappa)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa_snow               !! volumetric vertically discretized snow heat 
                                                                              !! capacity @tex ($J K^{-1} m^{-3}$) @endtex. 
!$OMP THREADPRIVATE(pcapa_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pkappa_snow              !! vertically discretized snow thermal conductivity 
                                                                              !! @tex ($W K^{-1} m^{-1}$) @endtex.
!$OMP THREADPRIVATE(pkappa_snow)

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcapa_en                 !! heat capacity used for surfheat_incr and 
                                                                              !! coldcont_incr 
!$OMP THREADPRIVATE(pcapa_en)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: ptn_beg                  !! ptn as it is after thermosoil_profile but before thermosoil_coef, 
                                                                              !! used in thermosoil_readjust
!$OMP THREADPRIVATE(ptn_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: temp_sol_beg             !! Surface temperature at previous timestep (K) 
!$OMP THREADPRIVATE(temp_sol_beg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: surfheat_incr            !! Change in soil heat content during the timestep 
                                                                              !!  @tex ($J$) @endtex.
!$OMP THREADPRIVATE(surfheat_incr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: coldcont_incr            !! Change in snow heat content  @tex ($J$) @endtex.
!$OMP THREADPRIVATE(coldcont_incr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: shum_ngrnd_perma         !! Saturation degree on the thermal axes (0-1, dimensionless)
!$OMP THREADPRIVATE(shum_ngrnd_perma)

  !  Variables related to soil freezing
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: profil_froz              !! Frozen fraction of the soil on hydrological levels (-)
!$OMP THREADPRIVATE(profil_froz)
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:) :: e_soil_lat                  !! Latent heat released or consumed in the freezing/thawing processes summed vertically
                                                                              !! for the whole soil (J/m2) and on the whole simulation to check/correct energy conservation
!$OMP THREADPRIVATE(e_soil_lat)
  REAL(r_std),ALLOCATABLE, SAVE, DIMENSION (:,:) :: pcappa_supp               !! Additional surfacic heat capacity due to soil freezing for each soil layer (J/K/m2)
!$OMP THREADPRIVATE(pcappa_supp)    
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: dz5                      !! Used for numerical calculation [-]
!$OMP THREADPRIVATE(dz5)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: QZ                       !! quartz content [-]
!$OMP THREADPRIVATE(QZ)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: so_capa_dry              !! Dry soil Heat capacity of soils,J.m^{-3}.K^{-1} 
!$OMP THREADPRIVATE(so_capa_dry)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: so_capa_ice              !! Heat capacity of saturated frozen soil (J/K/m3)
!$OMP THREADPRIVATE(so_capa_ice)   
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mc_layt                  !! Volumetric soil moisture (liquid+ice) (m3/m3) on the thermodynamical levels at interface
!$OMP THREADPRIVATE(mc_layt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: mcl_layt                 !! Volumetric soil moisture (liquid) (m3/m3) on the thermodynamical levels at interface
!$OMP THREADPRIVATE(mcl_layt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:) :: tmc_layt                 !! Total soil moisture content for each layer (liquid+ice) (mm) on the thermodynamical levels
!$OMP THREADPRIVATE(tmc_layt)
  INTEGER(i_std), SAVE                            :: brk_flag = 0             !! Flag to consider bedrock: 0.no; 1.yes
!$OMP THREADPRIVATE(brk_flag)
  INTEGER(i_std), SAVE                            :: nvm2 = 16                !! number of PFTs
!$OMP THREADPRIVATE(brk_flag)


CONTAINS


!! =============================================================================================================================
!! SUBROUTINE:    thermosoil_xios_initialize
!!
!>\BRIEF	  Initialize xios dependant defintion before closing context defintion
!!
!! DESCRIPTION:	  Initialize xios dependant defintion before closing context defintion
!!                Reading is deactivated if the sechiba restart file exists because the
!!                variable should be in the restart file already. 
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE thermosoil_xios_initialize

    CHARACTER(LEN=255) :: filename, name
        
    filename = 'reftemp.nc'
    CALL getin_p('REFTEMP_FILE',filename)
    
    name = filename(1:LEN_TRIM(FILENAME)-3)
    CALL xios_orchidee_set_file_attr("reftemp_file",name=name)

    ! Check if the reftemp file will be read by XIOS, by IOIPSL or not at all    
    IF (xios_interpolation .AND. read_reftemp .AND. restname_in=='NONE') THEN
       ! The reftemp file will be read using XIOS
       IF (printlev>=2) WRITE(numout,*) 'Reading of reftemp file will be done later using XIOS. The filename is ', filename
    ELSE
       IF (.NOT. read_reftemp) THEN
          IF (printlev>=2) WRITE (numout,*) 'No reading of reftemp will be done because read_reftemp=FALSE'
       ELSE IF (restname_in=='NONE') THEN
          IF (printlev>=2) WRITE (numout,*) 'The reftemp file will be read later by IOIPSL'
       ELSE
          IF (printlev>=2) WRITE (numout,*) 'The reftemp file will not be read because the restart file exists.'
       END IF

       ! The reftemp file will not be read by XIOS. Now deactivate albedo for XIOS.
       CALL xios_orchidee_set_file_attr("reftemp_file",enabled=.FALSE.)
       CALL xios_orchidee_set_field_attr("reftemp_interp",enabled=.FALSE.)
    ENDIF

  END SUBROUTINE thermosoil_xios_initialize

  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : thermosoil_initialize
  !!
  !>\BRIEF			            Allocate module variables, read from restart file or initialize with default values
  !!
  !! DESCRIPTION			    : Allocate module variables, read from restart file or initialize with default values.
  !!                                          Call thermosoil_var_init to calculate physical constants. 
  !!                                          Call thermosoil_coef to calculate thermal soil properties.
  !!
  !! RECENT CHANGE(S)			    : None
  !!
  !! REFERENCE(S)			    : None
  !! 
  !! FLOWCHART                              : None
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE thermosoil_initialize (kjit,          kjpindex,   rest_id,          mcs,       &
                                    temp_sol_new,  snow,       shumdiag_perma,              &
                                    soilcap,       soilflx,    stempdiag,        ftempdiag, &
                                    gtemp,                   &
                                    mc_layh,       mcl_layh,   tmc_layh,        njsc,     &
                                    frac_snow_veg,frac_snow_nobio,totfrac_nobio, &
                                    snowdz, snowrho, snowtemp, lambda_snow, cgrnd_snow, dgrnd_snow, pb, veget_max)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id          !! Restart file identifier (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: mcs              !! Saturated moisture content (m3/m3)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new     !! Surface temperature at the present time-step,
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: snow             !! Snow mass (kg)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma   !! Soil saturation degree (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mc_layh          !! Volumetric soil moisture content (liquid+ice) for hydrological layers, at node (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mcl_layh         !! Volumetric soil moisture content (liquid) for hydrological layers, at node (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: tmc_layh         !! Total soil moisture content(liquid+ice) for hydrological layers (mm)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)   :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: totfrac_nobio    !! Total fraction of continental ice+lakes+cities+...
                                                                              !! (unitless,0-1)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowdz           !! Snow depth
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)   :: snowrho          !! Snow density
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)   :: snowtemp         !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)        :: pb               !! Surface presure (hPa)
    REAL(r_std), DIMENSION (kjpindex,nvm2), INTENT (in)    :: veget_max       !! Fraction of PFT (unitless,0-1) 01

    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: soilcap          !! apparent surface heat capacity considering snow and soil surface (J m-2 K-1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: soilflx          !! apparent soil heat flux considering snow and soil surface (W m-2)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)   :: stempdiag        !! temperature profile on the levels in hydrol(K)
    REAL(r_std),DIMENSION (kjpindex,ngrnd), INTENT (out)  :: ftempdiag        !! temperature profile on full depth for stream temperature (K)
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)          :: gtemp            !! First soil layer temperature

    !! 0.3 Modified variables
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow     !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow      !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow      !! Integration coefficient for snow numerical scheme

    !! 0.4 Local variables
    INTEGER(i_std)                                        :: ier, i, jg, jsc
    LOGICAL                                               :: calculate_coef   !! Local flag to initialize variables by call to thermosoil_coef   
!_ ================================================================================================================================
    

    !
    !  !! Flag to consider bedrock at deeper layers
    !  !! It affects heat capacity and thermal conductivity (energy balance). 
    !
    !Config Key  = BEDROCK_FLAG
    !Config Desc = Flag to consider bedrock at deeper layers.
    !Config If   = 
    !Config Def  = 0
    !Config Help = 0, no, 1, yes. 
    !Config Units = [FLAG]
    brk_flag = 0
    CALL getin_p('BEDROCK_FLAG', brk_flag)

    IF (printlev >= 3) WRITE (numout,*) 'Start thermosoil_initialize '

    !! 1. Allocate soil temperatures variables
    ALLOCATE (ptn(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of ptn','','')

    ALLOCATE (dz1(ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of dz1','','')

    ALLOCATE (cgrnd(kjpindex,ngrnd-1),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of cgrnd','','')

    ALLOCATE (dgrnd(kjpindex,ngrnd-1),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of dgrnd','','')

    ALLOCATE (pcapa(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa','','')

    ALLOCATE (pkappa(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pkappa','','')

    ALLOCATE (pcapa_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa_snow','','')

    ALLOCATE (pkappa_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pkappa_snow','','')

    ! Temporary fix: Initialize following variable because they are output to xios before the first calculation
    pcapa  = 0
    pkappa = 0
    pcapa_snow  = 0
    pkappa_snow = 0

    ALLOCATE (surfheat_incr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of surfheat_incr','','')

    ALLOCATE (coldcont_incr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of coldcont_incr','','')

    ALLOCATE (pcapa_en(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of pcapa_en','','')
    ! Initialization to zero used at first time step in thermosoil_energy_diag, only for diagnostic variables coldcont_incr and surfheat_incr
    pcapa_en(:,:) = 0.

    ALLOCATE (ptn_beg(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of ptn_beg','','')

    ALLOCATE (temp_sol_beg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of temp_sol_beg','','')

    ALLOCATE (shum_ngrnd_perma(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of shum_ngrnd_perma','','')

    ALLOCATE (profil_froz(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of profil_froz','','')

    IF (ok_freeze_thermix) THEN
       ALLOCATE (pcappa_supp(kjpindex,ngrnd),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of ok_freeze_termix','','')
       ! Initialization to zero used at first time step only for diagnostic output. 
       ! This variable is only used in thermosoil_readajust and always calculated before in thermosoil_getdiff.
       pcappa_supp(:,:) = 0.
    END IF
    IF (ok_Ecorr) THEN
       ALLOCATE (e_soil_lat(kjpindex),stat=ier)
       IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of e_soil_lat','','')
    END IF

    ALLOCATE (dz5(ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of dz5','','')

    ALLOCATE (mc_layt(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of mc_layt','','')

    ALLOCATE (mcl_layt(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of mcl_layt','','')

    ALLOCATE (tmc_layt(kjpindex,ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of tmc_layt','','')

    ALLOCATE (QZ(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of QZ','','')

    ALLOCATE (so_capa_dry(nscm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of so_capa_dry','','')
    
    ALLOCATE (so_capa_ice(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'thermosoil_initialize', 'Error in allocation of so_capa_ice','','')


    !! Soil texture choose : Now useless since njsc defines the dominant texture within 13 classes whichever the soil map
    QZ(:) = QZ_usda(:)
    so_capa_dry(:) = so_capa_dry_usda(:)
    
    !Config Key   = DRY_SOIL_HEAT_CAPACITY
    !Config Desc  = Dry soil Heat capacity of soils
    !Config If    = OK_SECHIBA
    !Config Def   = (1.47, 1.41, 1.34, 1.27, 1.21, 1.21, 1.18, 1.32, 1.23, 1.18, 1.15, 1.09, 1.09)*e+6
    !Config Help  = Values taken from : Pielke [2002, 2013]
    !Config Units = [J.m^{-3}.K^{-1}] 
    CALL getin_p("DRY_SOIL_HEAT_CAPACITY",so_capa_dry)
    
    !! Check parameter value (correct range)
    IF ( MINVAL(so_capa_dry(:)) <= zero ) THEN
       CALL ipslerr_p(3, "thermosoil_initialize", &
            "Wrong parameter value for DRY_SOIL_HEAT_CAPACITY.", &
            "This parameter should be positive. ", &
            "Please, check parameter value in run.def or orchidee.def. ")
    END IF


    !! 2. Initialize variable from restart file or with default values 
    
    !! Reads restart files for soil temperatures only. If no restart file is
    !! found,  the initial soil temperature is by default set to 280K at all depths. The user
    !! can decide to initialize soil temperatures at an other value, in which case he should set the flag THERMOSOIL_TPRO
    !! to this specific value in the run.def.
    IF (printlev>=3) WRITE (numout,*) 'Read restart file for THERMOSOIL variables'

    CALL ioconf_setatt_p('UNITS', 'K')
    CALL ioconf_setatt_p('LONG_NAME','Soil Temperature profile')
    CALL restget_p (rest_id, 'ptn', nbp_glo, ngrnd, 1, kjit, .TRUE., ptn, "gather", nbp_glo, index_g)

    ! Initialize ptn if it was not found in restart file
    IF (ALL(ptn(:,:)==val_exp)) THEN 
       ! ptn was not found in restart file

       IF (read_reftemp) THEN
          ! Read variable ptn from file
          CALL thermosoil_read_reftempfile(kjpindex,lalo,ptn)
       ELSE
          ! Initialize ptn with a constant value which can be set in run.def

          !Config Key   = THERMOSOIL_TPRO
          !Config Desc  = Initial soil temperature profile if not found in restart
          !Config Def   = 280.
          !Config If    = OK_SECHIBA
          !Config Help  = The initial value of the temperature profile in the soil if 
          !Config         its value is not found in the restart file. Here
          !Config         we only require one value as we will assume a constant 
          !Config         throughout the column.
          !Config Units = Kelvin [K]
          CALL setvar_p (ptn, val_exp,'THERMOSOIL_TPRO',280._r_std)
       END IF
    END IF
    
    ! Initialize ptn_beg (variable needed in thermosoil_readadjust called from thermosoil_coef)
    ptn_beg(:,:) = ptn(:,:)
    
    ! Initialize temp_sol_beg with values from previous time-step
    temp_sol_beg(:) = temp_sol_new(:) 
    
    ! Read e_soil_lat from restart file or initialize
    IF (ok_Ecorr) THEN
       CALL restget_p (rest_id, 'e_soil_lat', nbp_glo, 1, 1, kjit, .TRUE., &
            e_soil_lat, "gather", nbp_glo, index_g)
       CALL setvar_p (e_soil_lat, val_exp,'NO_KEYWORD',zero)
    END IF

    ! Read gtemp from restart file
    CALL restget_p (rest_id, 'gtemp', nbp_glo, 1, 1, kjit, .TRUE., &
         gtemp, "gather", nbp_glo, index_g)
    CALL setvar_p (gtemp, val_exp,'NO_KEYWORD',zero)
    

    ! Read variables calculated in thermosoil_coef from restart file
    ! If the variables were not found in the restart file, the logical 
    ! calculate_coef will be true and thermosoil_coef will be called further below.
    ! These variables need to be in the restart file to avoid a time shift that
    ! would be done using thermosoil_coef at this stage.
    calculate_coef=.FALSE.
    CALL ioconf_setatt_p('UNITS', 'J m-2 K-1')
    CALL ioconf_setatt_p('LONG_NAME','Apparent surface heat capacity')
    CALL restget_p (rest_id, 'soilcap', nbp_glo, 1, 1, kjit, .TRUE., &
         soilcap, "gather", nbp_glo, index_g)
    IF (ALL(soilcap(:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', 'W m-2')
    CALL ioconf_setatt_p('LONG_NAME','Apparent soil heat flux')
    CALL restget_p (rest_id, 'soilflx', nbp_glo, 1, 1, kjit, .TRUE., &
         soilflx, "gather", nbp_glo, index_g)
    IF (ALL(soilflx(:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'cgrnd', nbp_glo, ngrnd-1, 1, kjit, .TRUE., &
         cgrnd, "gather", nbp_glo, index_g)
    IF (ALL(cgrnd(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'dgrnd', nbp_glo, ngrnd-1, 1, kjit, .TRUE., &
         dgrnd, "gather", nbp_glo, index_g)
    IF (ALL(dgrnd(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'cgrnd_snow', nbp_glo, nsnow, 1, kjit, .TRUE., &
         cgrnd_snow, "gather", nbp_glo, index_g)
    IF (ALL(cgrnd_snow(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Integration coefficient for the numerical scheme')
    CALL restget_p (rest_id, 'dgrnd_snow', nbp_glo, nsnow, 1, kjit, .TRUE., &
         dgrnd_snow, "gather", nbp_glo, index_g)
    IF (ALL(dgrnd_snow(:,:)==val_exp)) calculate_coef=.TRUE.

    CALL ioconf_setatt_p('UNITS', '')
    CALL ioconf_setatt_p('LONG_NAME','Coefficient of the linear extrapolation of surface temperature')
    CALL restget_p (rest_id, 'lambda_snow', nbp_glo, 1, 1, kjit, .TRUE., &
         lambda_snow, "gather", nbp_glo, index_g)
    IF (ALL(lambda_snow(:)==val_exp)) calculate_coef=.TRUE.

    !! 2.2 Computes some physical constants and arrays depending on the soil vertical discretization 

    ! Calculate so_capa_ice
    so_capa_ice(:) = capa_ice*rho_ice 
    IF (printlev>=2) WRITE(numout,*) 'Calculation of so_capa_ice(:)=', so_capa_ice(:),' and capa_ice=',capa_ice
    
    ! Computing some usefull constants for the numerical scheme
    ! Use znt(depth of nodes) and zlt(depth of deeper layer interface) from vertical_soil module.  
    DO jg=1,ngrnd-1
      dz1(jg)  = un / (znt(jg+1) - znt(jg))
      dz5(jg) = (zlt(jg) - znt(jg)) * dz1(jg)
    ENDDO
    dz5(ngrnd) = 0.0
    lambda = znt(1) * dz1(1)

    ! Send out the temperature profile on the first nslm levels(the levels treated in hydrol)
    stempdiag(:,:) = ptn(:,1:nslm)
    ftempdiag(:,:) = ptn(:,1:ngrnd)
    

    !! 2.3. Computes cgrnd, dgrnd, soilflx and soilcap coefficients only if they were not found in restart file.
    IF (calculate_coef) THEN
       ! Interpolate variables needed by thermosoil_coef to the thermal levels
       CALL thermosoil_humlev(kjpindex, shumdiag_perma, mc_layh, mcl_layh, tmc_layh)

       IF (printlev>=3) WRITE (numout,*) 'thermosoil_coef will be called in the intialization phase'
       CALL thermosoil_coef (&
            kjpindex,      temp_sol_new,    snow,           njsc, &
            mcs, frac_snow_veg, frac_snow_nobio, totfrac_nobio,        &
            snowdz,        snowrho,         snowtemp,       pb,   &
            ptn,                                                  &
            soilcap,       soilflx,         cgrnd,          dgrnd,&
	    lambda_snow,   cgrnd_snow,      dgrnd_snow, veget_max)
    END IF

  END SUBROUTINE thermosoil_initialize


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_main
!!
!>\BRIEF        Thermosoil_main computes the soil thermal properties and dynamics, ie solves
!! the heat diffusion equation within the soil. 
!!
!! DESCRIPTION : The resolution of the soil heat diffusion equation 
!! relies on a numerical finite-difference implicit scheme
!! fully described in the reference and in the header of the thermosoil module.
!! - The dependency of the previous timestep hidden in the 
!! integration coefficients cgrnd and dgrnd (EQ1), calculated in thermosoil_coef, and 
!! called at the end of the routine to prepare for the next timestep.
!! - The effective computation of the new soil temperatures is performed in thermosoil_profile. 
!!
!! - thermosoil_coef calculates the coefficients for the numerical scheme for the very first iteration of thermosoil;
!! after that, thermosoil_coef is called only at the end of the module to calculate the coefficients for the next timestep.
!! - thermosoil_profile solves the numerical scheme.\n
!!
!! - Flags : one unique flag : THERMOSOIL_TPRO (to be set to the desired initial soil in-depth temperature in K; by default 280K)
!!
!! RECENT CHANGE(S) : Change vertical discretization (consistent with hydrology layers) and soil thermal properties (taking into account soil texture effects).
!!
!! MAIN OUTPUT VARIABLE(S): vertically discretized soil temperatures ptn, soil
!! thermal properties (pcapa, pkappa), apparent surface heat capacity (soilcap)
!! and heat flux (soilflx) to be used in enerbil at the next timestep to solve
!! the surface energy balance.
!!
!! REFERENCE(S) : 
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!!  Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin' s PhD thesis relative to the thermal
!!  integration scheme has been scanned and is provided along with the documentation, with name : 
!!  Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : 
!! \latexonly
!! \includegraphics[scale = 1]{thermosoil_flowchart.png}
!! \endlatexonly
!! 
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_main (kjit, kjpindex, &
       index, indexgrnd, mcs, &
       temp_sol_new, snow, soilcap, soilflx, &
       shumdiag_perma, stempdiag, ftempdiag, ptnlev1, rest_id, hist_id, hist2_id, &
       snowdz,snowrho,snowtemp,gtemp,pb,&
       mc_layh, mcl_layh, tmc_layh, njsc, frac_snow_veg,frac_snow_nobio,totfrac_nobio,temp_sol_add, &
       lambda_snow, cgrnd_snow, dgrnd_snow, veget_max)

    !! 0. Variable and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id,hist_id  !! Restart_ file and history file identifier 
                                                                              !! (unitless)
    INTEGER(i_std),INTENT (in)                            :: hist2_id         !! history file 2 identifier (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: index            !! Indeces of the points on the map (unitless)
    INTEGER(i_std),DIMENSION (kjpindex*ngrnd), INTENT (in):: indexgrnd        !! Indeces of the points on the 3D map (vertical 
                                                                              !! dimension towards the ground) (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: mcs              !! Saturated moisture content (m3/m3)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: temp_sol_new     !! Surface temperature at the present time-step,
                                                                              !! Ts @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: snow             !! Snow mass @tex ($kg$) @endtex.
                                                                              !! Caveat: when there is snow on the
                                                                              !! ground, the snow is integrated into the soil for
                                                                              !! the calculation of the thermal dynamics. It means
                                                                              !! that the uppermost soil layers can completely or 
                                                                              !! partially consist in snow. In the second case, zx1
                                                                              !! and zx2 are the fraction of the soil layer 
                                                                              !! consisting in snow and 'normal' soil, respectively
                                                                              !! This is calculated in thermosoil_coef.
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma   !! Soil saturation degree (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowdz           !! Snow depth
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (in)    :: snowrho          !! Snow density
    REAL(r_std),DIMENSION (kjpindex,nsnow),INTENT (inout) :: snowtemp         !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex),INTENT (in)         :: pb               !! Surface presure (hPa)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mc_layh          !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid + ice) (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mcl_layh         !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid) (m3/m3)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: tmc_layh         !! Total soil moisture content for each layer in hydrol(liquid + ice) (mm)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)      :: njsc             !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)          :: frac_snow_veg    !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)   :: frac_snow_nobio  !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: totfrac_nobio    !! Total fraction of continental ice+lakes+cities+...
                                                                              !!(unitless,0-1)
    REAL(r_std),DIMENSION (kjpindex,nvm2),INTENT(in)           :: veget_max        !! Fraction of PFT (unitless,0-1)   2                                                                      
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: temp_sol_add     !! additional surface temperature due to the melt of first layer
                                                                              !! at the present time-step @tex ($K$) @endtex

    !! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)        :: ptnlev1          !! 1st level soil temperature   
    REAL(r_std),DIMENSION (kjpindex),INTENT(out)          :: gtemp            !! First soil layer temperature


    !! 0.3 Modified variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: soilcap          !! apparent surface heat capacity considering snow and soil surface
                                                                              !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)      :: soilflx          !! apparent soil heat flux  considering snow and soil surface
                                                                              !! @tex ($W m^{-2}$) @endtex
                                                                              !! , positive 
                                                                              !! towards the soil, writen as Qg (ground heat flux) 
                                                                              !! in the history files, and computed at the end of 
                                                                              !! thermosoil for the calculation of Ts in enerbil, 
                                                                              !! see EQ3.
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)   :: stempdiag        !! temperature profile @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex,ngrnd), INTENT (out)  :: ftempdiag        !! temperature profile @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow      !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow       !! Integration coefficient for snow numerical scheme

    !! 0.4 Local variables

    INTEGER(i_std)                                        :: jv,ji,ii
    REAL(r_std),DIMENSION (kjpindex)                      :: snowtemp_weighted!! Snow temperature weighted by snow density, only for diag (K)
    REAL(r_std),DIMENSION (kjpindex, nsnow)               :: pkappa_snow_diag !! Only for diag, containing xios_default_val
    REAL(r_std),DIMENSION (kjpindex, nsnow)               :: pcapa_snow_diag  !! Only for diag, containing xios_default_val
    REAL(r_std),DIMENSION (kjpindex, nsnow)               :: snowtemp_diag    !! Only for diag, containing xios_default_val

!_ ================================================================================================================================
    
  !! 3. Put the soil wetness diagnostic on the levels of the soil temperature

    !!?? this could logically be put just before the last call to
    !!thermosoil_coef, as the results are used there...
    CALL thermosoil_humlev(kjpindex, shumdiag_perma, mc_layh, mcl_layh, tmc_layh)

    
  !! 4. Effective computation of the soil temperatures profile.
  !!    cgrnd and dgrnd have been calculated in thermosoil_coef at the previous time step 
  !!    but they are correct for the actual time-step.
    CALL thermosoil_profile (kjpindex,      temp_sol_new,                   &
                             frac_snow_veg, frac_snow_nobio, totfrac_nobio, &
                             ptn,           stempdiag,       snowtemp,      &
                             cgrnd_snow,    dgrnd_snow)


  !! 5. Call to thermosoil_energy_diag for calculation of diagnostic variables
    CALL thermosoil_energy_diag(kjpindex, temp_sol_new, soilcap)

  !! Save ptn at current stage, to be used in thermosoil_readjust
    ptn_beg(:,:) = ptn(:,:)

  !! 6. Writing the history files according to the ALMA standards (or not..)

    ! Add XIOS default value where no snow
    DO ji=1,kjpindex 
       IF (snow(ji) .GT. zero) THEN
          pkappa_snow_diag(ji,:) = pkappa_snow(ji,:)
          pcapa_snow_diag(ji,:) = pcapa_snow(ji,:)
          snowtemp_diag(ji,:) = snowtemp(ji,:)
       ELSE
          pkappa_snow_diag(ji,:) = xios_default_val
          pcapa_snow_diag(ji,:) = xios_default_val
          snowtemp_diag(ji,:) = xios_default_val
       END IF
    END DO

    DO ji=1,kjpindex 
       ! Use min_sechiba instead of zero to avoid problem with division by zero
       IF (snow(ji) .GT. min_sechiba) THEN
          snowtemp_weighted(ji) = SUM(snowtemp(ji,:)*snowrho(ji,:))/SUM(snowrho(ji,:))
       ELSE
          snowtemp_weighted(ji) = xios_default_val
       END IF
    END DO
    CALL xios_orchidee_send_field("snowtemp_weighted",snowtemp_weighted)
     
    CALL xios_orchidee_send_field("ptn",ptn)
    CALL xios_orchidee_send_field("soilflx",soilflx)
    CALL xios_orchidee_send_field("surfheat_incr",surfheat_incr)
    CALL xios_orchidee_send_field("coldcont_incr",coldcont_incr)
    CALL xios_orchidee_send_field("pkappa",pkappa)
    CALL xios_orchidee_send_field("pkappa_snow",pkappa_snow_diag)
    CALL xios_orchidee_send_field("pcapa",pcapa)
    CALL xios_orchidee_send_field("pcapa_snow",pcapa_snow_diag)
    CALL xios_orchidee_send_field("snowtemp",snowtemp_diag)


    IF ( .NOT. almaoutput ) THEN
      CALL histwrite_p(hist_id, 'ptn', kjit, ptn, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'Qg', kjit, soilflx, kjpindex, index)
      CALL histwrite_p(hist_id, 'ptn_beg', kjit, ptn_beg, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'pkappa', kjit, pkappa, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'pcapa', kjit, pcapa, kjpindex*ngrnd, indexgrnd)
      
      IF (ok_freeze_thermix) THEN
         CALL histwrite_p(hist_id, 'profil_froz', kjit, profil_froz, kjpindex*ngrnd, indexgrnd)
         CALL histwrite_p(hist_id, 'pcappa_supp', kjit, pcappa_supp, kjpindex*ngrnd, indexgrnd)
      END IF
      CALL histwrite_p(hist_id, 'shum_ngrnd_perma', kjit, shum_ngrnd_perma(:,:), kjpindex*ngrnd, indexgrnd)
      
    ELSE
      CALL histwrite_p(hist_id, 'SoilTemp', kjit, ptn, kjpindex*ngrnd, indexgrnd)
      CALL histwrite_p(hist_id, 'Qg', kjit, soilflx, kjpindex, index)
      CALL histwrite_p(hist_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
      CALL histwrite_p(hist_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
    ENDIF
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
          CALL histwrite_p(hist2_id, 'ptn', kjit, ptn, kjpindex*ngrnd, indexgrnd)
       ELSE
          CALL histwrite_p(hist2_id, 'SoilTemp', kjit, ptn, kjpindex*ngrnd, indexgrnd)
          CALL histwrite_p(hist2_id, 'Qg', kjit, soilflx, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelSurfHeat', kjit, surfheat_incr, kjpindex, index)
          CALL histwrite_p(hist2_id, 'DelColdCont', kjit, coldcont_incr, kjpindex, index)
       ENDIF
    ENDIF
    
  !! 7. A last final call to thermosoil_coef
 
    !! A last final call to thermosoil_coef, which calculates the different
    !!coefficients (cgrnd, dgrnd, soilcap, soilflx) from this time step to be
    !!used at the next time step, either in the surface temperature calculation
    !!(soilcap, soilflx) or in the soil thermal numerical scheme.
    CALL thermosoil_coef (&
         kjpindex,      temp_sol_new,    snow,         njsc, &
         mcs, frac_snow_veg, frac_snow_nobio, totfrac_nobio, &
         snowdz,        snowrho,         snowtemp,     pb,   &
         ptn,                                                &
         soilcap,       soilflx,         cgrnd,        dgrnd,&
         lambda_snow,   cgrnd_snow,      dgrnd_snow, veget_max)
         

    ! Save variables for explicit snow model
    gtemp(:) = ptn(:,1)

    !! Initialize output arguments to be used in sechiba
    ptnlev1(:) = ptn(:,1)

    !! Write the new temperature in the full in a the diagnostic variable.
    ftempdiag(:,:) = ptn(:,1:ngrnd)

    !! Surface temperature is forced to zero celcius if its value is larger than melting point
    DO ji=1,kjpindex
       IF  (SUM(snowdz(ji,:)) .GT. 0.0) THEN
          IF (temp_sol_new(ji) .GE. tp_00) THEN
             temp_sol_new(ji) = tp_00
          ENDIF
       END IF
    END DO
    
    IF (printlev>=3) WRITE (numout,*) ' thermosoil_main done '

  END SUBROUTINE thermosoil_main

  !!  =============================================================================================================================
  !! SUBROUTINE		 		    : thermosoil_finalize
  !!
  !>\BRIEF                                    Write to restart file
  !!
  !! DESCRIPTION			    : This subroutine writes the module variables and variables calculated in thermosoil
  !!                                          to restart file
  !! \n
  !_ ==============================================================================================================================
  SUBROUTINE thermosoil_finalize (kjit,    kjpindex, rest_id,   gtemp, &
                                  soilcap, soilflx, lambda_snow, cgrnd_snow, dgrnd_snow)

    !! 0. Variable and parameter declaration
    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                            :: kjit             !! Time step number (unitless) 
    INTEGER(i_std), INTENT(in)                            :: kjpindex         !! Domain size (unitless)
    INTEGER(i_std),INTENT (in)                            :: rest_id          !! Restart file identifier(unitless)
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: gtemp            !! First soil layer temperature
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: soilcap
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)           :: soilflx
    REAL(r_std), DIMENSION (kjpindex), INTENT(in)         :: lambda_snow      !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (in)  :: cgrnd_snow       !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (in)  :: dgrnd_snow       !! Integration coefficient for snow numerical scheme

!_ ================================================================================================================================
    
    !! 1. Write variables to restart file to be used for the next simulation
    IF (printlev>=3) WRITE (numout,*) 'Write restart file with THERMOSOIL variables'
    
    CALL restput_p(rest_id, 'ptn', nbp_glo, ngrnd, 1, kjit, ptn, 'scatter', nbp_glo, index_g)
    
    IF (ok_Ecorr) THEN
       CALL restput_p(rest_id, 'e_soil_lat', nbp_glo, 1 , 1, kjit, e_soil_lat, 'scatter', nbp_glo, index_g)
    END IF
    
    CALL restput_p(rest_id, 'gtemp', nbp_glo, 1, 1, kjit, gtemp, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'soilcap', nbp_glo, 1, 1, kjit, soilcap, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'soilflx', nbp_glo, 1, 1, kjit, soilflx, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'cgrnd', nbp_glo, ngrnd-1, 1, kjit, cgrnd, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'dgrnd', nbp_glo, ngrnd-1, 1, kjit, dgrnd, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'cgrnd_snow', nbp_glo, nsnow, 1, kjit, cgrnd_snow, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'dgrnd_snow', nbp_glo, nsnow, 1, kjit, dgrnd_snow, 'scatter', nbp_glo, index_g)
    CALL restput_p(rest_id, 'lambda_snow', nbp_glo, 1, 1, kjit, lambda_snow, 'scatter', nbp_glo, index_g)
    
  END SUBROUTINE thermosoil_finalize


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_clear
!!
!>\BRIEF        Deallocates the allocated arrays.
!! The call of thermosoil_clear originates from sechiba_clear but the calling sequence and 
!! its purpose require further investigation.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S) : None 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None 
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

 SUBROUTINE thermosoil_clear()

        IF ( ALLOCATED (ptn)) DEALLOCATE (ptn)
        IF ( ALLOCATED (cgrnd)) DEALLOCATE (cgrnd) 
        IF ( ALLOCATED (dgrnd)) DEALLOCATE (dgrnd) 
        IF ( ALLOCATED (pcapa)) DEALLOCATE (pcapa)
        IF ( ALLOCATED (pkappa))  DEALLOCATE (pkappa)
        IF ( ALLOCATED (pcapa_snow)) DEALLOCATE (pcapa_snow)
        IF ( ALLOCATED (pkappa_snow))  DEALLOCATE (pkappa_snow)
        IF ( ALLOCATED (pcapa_en)) DEALLOCATE (pcapa_en)
        IF ( ALLOCATED (ptn_beg)) DEALLOCATE (ptn_beg)
        IF ( ALLOCATED (temp_sol_beg)) DEALLOCATE (temp_sol_beg)
        IF ( ALLOCATED (surfheat_incr)) DEALLOCATE (surfheat_incr)
        IF ( ALLOCATED (coldcont_incr)) DEALLOCATE (coldcont_incr)
        IF ( ALLOCATED (shum_ngrnd_perma)) DEALLOCATE (shum_ngrnd_perma)
        IF ( ALLOCATED (profil_froz)) DEALLOCATE (profil_froz)
        IF ( ALLOCATED (mc_layt)) DEALLOCATE (mc_layt)
        IF ( ALLOCATED (mcl_layt)) DEALLOCATE (mcl_layt)
        IF ( ALLOCATED (tmc_layt)) DEALLOCATE (tmc_layt)
  END SUBROUTINE thermosoil_clear


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_coef
!!
!>\BRIEF        Calculate soil thermal properties, integration coefficients, apparent heat flux,
!! surface heat capacity,  
!!
!! DESCRIPTION	: This routine computes : \n
!!		1. the soil thermal properties. \n 
!!		2. the integration coefficients of the thermal numerical scheme, cgrnd and dgrnd,
!!              which depend on the vertical grid and on soil properties, and are used at the next 
!!              timestep.\n
!!              3. the soil apparent heat flux and surface heat capacity (soilflx
!!              and soilcap), used by enerbil to compute the surface temperature at the next
!!              timestep.\n
!!             -  The soil thermal properties depend on water content (shum_ngrnd_perma, shumdiag_perma, 
!!              mc_layt, mcl_layt, tmc_layt), dominant soil texture(njsc), and on the presence 
!!              of snow : snow is integrated into the soil for the thermal calculations, ie if there 
!!              is snow on the ground, the first thermal layer(s) consist in snow, depending on the 
!!              snow-depth. If a layer consists out of snow and soil, wheighed soil properties are 
!!              calculated\n
!!             - The coefficients cgrnd and dgrnd are the integration
!!              coefficients for the thermal scheme \n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k) \n
!!                                      -- EQ1 -- \n
!!              They correspond respectively to $\beta$ and $\alpha$ from F. Hourdin\'s thesis and 
!!              their expression can be found in this document (eq A19 and A20)
!!             - soilcap and soilflx are the apparent surface heat capacity and flux
!!               used in enerbil at the next timestep to solve the surface
!!               balance for Ts (EQ3); they correspond to $C_s$ and $F_s$ in F.
!!               Hourdin\'s PhD thesis and are expressed in eq. A30 and A31. \n
!!                 soilcap*(Ts(t)-Ts(t-1))/dt=soilflx+otherfluxes(Ts(t)) \n
!!                                      -- EQ3 --\n
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cgrnd, dgrnd, pcapa, pkappa, soilcap, soilflx
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE thermosoil_coef (kjpindex,      temp_sol_new,    snow,           njsc, &
                              mcs,           frac_snow_veg,   frac_snow_nobio,totfrac_nobio, &
                              snowdz,        snowrho,         snowtemp,       pb,   &
                              ptn,                                                  &
                              soilcap,       soilflx,         cgrnd,          dgrnd,&
			      lambda_snow,   cgrnd_snow,      dgrnd_snow, veget_max)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                             :: kjpindex     !! Domain size (unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: temp_sol_new !! soil surface temperature @tex ($K$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: snow         !! snow mass @tex ($Kg$) @endtex
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)       :: njsc         !! Index of the dominant soil textural class
                                                                           !! in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)         :: mcs          !! Saturated moisture content (m3/m3)
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)           :: frac_snow_veg   !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)    :: frac_snow_nobio !! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)            :: totfrac_nobio   !! Total fraction of continental ice+lakes+cities+...
                                                                              !!(unitless,0-1)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowdz          !! Snow depth (m)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowrho         !! Snow density
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in)    :: snowtemp        !! Snow temperature (K)
    REAL(r_std), DIMENSION (kjpindex), INTENT (in)         :: pb              !! Surface presure (hPa)
    REAL(r_std), DIMENSION (kjpindex,nvm2), INTENT (in)         :: veget_max       !! Fraction of PFT (unitless,0-1)   3

    !! 0.2 Output variables

    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: soilcap      !! surface heat capacity considering snow and soil surface
                                                                           !! @tex ($J m^{-2} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex), INTENT (out)        :: soilflx      !! surface heat flux considering snow and soil surface @tex ($W m^{-2}$) @endtex,
                                                                           !! positive towards the 
                                                                           !! soil, writen as Qg (ground heat flux) in the history 
                                                                           !! files.
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: cgrnd        !! matrix coefficient for the computation of soil 
                                                                           !! temperatures (beta in F. Hourdin thesis)
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1), INTENT(out) :: dgrnd        !! matrix coefficient for the computation of soil 
                                                                           !! temperatures (alpha in F. Hourdin thesis)


    !! 0.3 Modified variable

    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT (inout):: ptn          !! vertically discretized soil temperatures. ptn is only modified if ok_Ecorr.
    REAL(r_std), DIMENSION (kjpindex), INTENT(inout)       :: lambda_snow  !! Coefficient of the linear extrapolation of surface temperature 
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (inout):: cgrnd_snow   !! Integration coefficient for snow numerical scheme
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT (inout):: dgrnd_snow   !! Integration coefficient for snow numerical scheme

    !! 0.4 Local variables

    INTEGER(i_std)                                         :: ji, jg
    REAL(r_std), DIMENSION (kjpindex,ngrnd-1)              :: zdz1         !! numerical (buffer) constant 
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                :: zdz2         !! numerical (buffer) constant  
                                                                           !! @tex ($W m^{-1} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex)                      :: z1           !! numerical constant @tex ($W m^{-1} K^{-1}$) @endtex
    REAL(r_std), DIMENSION (kjpindex)                      :: soilcap_nosnow      !! surface heat capacity
                                                                                  !! @tex ($J m^{-2} K^{-1}$)
                                                                                  !! @endtex
    REAL(r_std), DIMENSION (kjpindex)                      :: soilflx_nosnow      !! surface heat flux @tex ($W m^{-2}$) @endtex,
                                                                                  !! positive towards the soil, written as Qg
                                                                                  !!(ground heat flux in the history files).
    REAL(r_std), DIMENSION (kjpindex)                      :: snowcap             !! apparent snow heat capacity @tex ($J m^{-2} K^{-1}$)
    REAL(r_std), DIMENSION (kjpindex)                      :: snowflx             !! apparent snow-atmosphere heat flux @tex ($W m^{-2}$) @endtex
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: dz1_snow
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: ZSNOWDZM
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: dz2_snow
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: zdz1_snow
    REAL(r_std), DIMENSION (kjpindex,nsnow)                :: zdz2_snow
    REAL(r_std), DIMENSION (kjpindex)                      :: z1_snow
    REAL(r_std), DIMENSION (kjpindex)                      :: snowflxtot          !! Total snow flux (including snow on vegetated and bare soil and nobio areas)
                                                                                  !! @tex ($W m^{-2}$) @endtex
                                                                                  !! positive towards the soil

!_ ================================================================================================================================

  !! 1. Computation of the soil thermal properties
   
    ! Computation of the soil thermal properties; snow properties are also accounted for
    IF (ok_freeze_thermix) THEN
       CALL thermosoil_getdiff( kjpindex, snow, ptn, mcs, njsc, snowrho, snowtemp, pb, veget_max)
    ELSE
       ! Special case without soil freezing
       CALL thermosoil_getdiff_old_thermix_without_snow( kjpindex, mcs, njsc, snowrho, snowtemp, pb, veget_max)
    ENDIF

    ! Energy conservation : Correction to make sure that the same latent heat is released and 
    ! consumed during freezing and thawing
    IF (ok_Ecorr) THEN
       CALL thermosoil_readjust(kjpindex, ptn)
    ENDIF
    

    !! 2. Computation of the coefficients of the numerical integration scheme for the soil layers

    !! 2.1 Calculate numerical coefficients zdz1 and zdz2
    DO jg=1,ngrnd
      DO ji=1,kjpindex
        zdz2(ji,jg)=pcapa(ji,jg) * dlt(jg)/dt_sechiba
      ENDDO
    ENDDO
    
    DO jg=1,ngrnd-1
      DO ji=1,kjpindex
        zdz1(ji,jg) = dz1(jg) * pkappa(ji,jg)
      ENDDO
    ENDDO
    
    !! 2.2 Calculate coefficients cgrnd and dgrnd for soil
    DO ji = 1,kjpindex
      z1(ji) = zdz2(ji,ngrnd) + zdz1(ji,ngrnd-1)
      cgrnd(ji,ngrnd-1) = zdz2(ji,ngrnd) * ptn(ji,ngrnd) / z1(ji)
      dgrnd(ji,ngrnd-1) = zdz1(ji,ngrnd-1) / z1(ji)
    ENDDO

    DO jg = ngrnd-1,2,-1
      DO ji = 1,kjpindex
        z1(ji) = un / (zdz2(ji,jg) + zdz1(ji,jg-1) + zdz1(ji,jg) * (un - dgrnd(ji,jg)))
        cgrnd(ji,jg-1) = (ptn(ji,jg) * zdz2(ji,jg) + zdz1(ji,jg) * cgrnd(ji,jg)) * z1(ji)
        dgrnd(ji,jg-1) = zdz1(ji,jg-1) * z1(ji)
      ENDDO
    ENDDO


    !! 3. Computation of the coefficients of the numerical integration scheme for the snow layers

    !! 3.1 Calculate numerical coefficients zdz1_snow, zdz2_snow and lambda_snow
    DO ji = 1, kjpindex

       ! Calculate internal values
       DO jg = 1, nsnow
          ZSNOWDZM(ji,jg) = MAX(snowdz(ji,jg),psnowdzmin)
       ENDDO
       dz2_snow(ji,:)=ZSNOWDZM(ji,:)
       
       DO jg = 1, nsnow-1
          dz1_snow(ji,jg)  = 2.0 / (dz2_snow(ji,jg+1)+dz2_snow(ji,jg))
       ENDDO
       
       lambda_snow(ji) = dz2_snow(ji,1)/2.0 * dz1_snow(ji,1)
       
       DO jg=1,nsnow
          zdz2_snow(ji,jg)=pcapa_snow(ji,jg) * dz2_snow(ji,jg)/dt_sechiba
       ENDDO
       
       DO jg=1,nsnow-1
          zdz1_snow(ji,jg) = dz1_snow(ji,jg) * pkappa_snow(ji,jg)
       ENDDO
       
       ! the bottom snow
       zdz1_snow(ji,nsnow) = pkappa_snow(ji,nsnow) / ( zlt(1) + dz2_snow(ji,nsnow)/2 )
       
    ENDDO

    !! 3.2 Calculate coefficients cgrnd_snow and dgrnd_snow for snow
    DO ji = 1,kjpindex
       ! bottom level
       z1_snow(ji) = zdz2(ji,1)+(un-dgrnd(ji,1))*zdz1(ji,1)+zdz1_snow(ji,nsnow)
       cgrnd_snow(ji,nsnow) = (zdz2(ji,1) * ptn(ji,1) + zdz1(ji,1) * cgrnd(ji,1) ) / z1_snow(ji)
       dgrnd_snow(ji,nsnow) = zdz1_snow(ji,nsnow) / z1_snow(ji)
       
       ! next-to-bottom level
       z1_snow(ji) = zdz2_snow(ji,nsnow)+(un-dgrnd_snow(ji,nsnow))*zdz1_snow(ji,nsnow)+zdz1_snow(ji,nsnow-1)
       cgrnd_snow(ji,nsnow-1) = (zdz2_snow(ji,nsnow)*snowtemp(ji,nsnow)+&
            zdz1_snow(ji,nsnow)*cgrnd_snow(ji,nsnow))/z1_snow(ji)
       dgrnd_snow(ji,nsnow-1) = zdz1_snow(ji,nsnow-1) / z1_snow(ji)
       
       DO jg = nsnow-1,2,-1
          z1_snow(ji) = un / (zdz2_snow(ji,jg) + zdz1_snow(ji,jg-1) + zdz1_snow(ji,jg) * (un - dgrnd_snow(ji,jg)))
          cgrnd_snow(ji,jg-1) = (snowtemp(ji,jg) * zdz2_snow(ji,jg) + zdz1_snow(ji,jg) * cgrnd_snow(ji,jg)) * z1_snow(ji)
          dgrnd_snow(ji,jg-1) = zdz1_snow(ji,jg-1) * z1_snow(ji)
       ENDDO
    ENDDO



  !! 4. Computation of the apparent ground heat flux 
    !! Computation of apparent snow-atmosphere flux  
    DO ji = 1,kjpindex
       snowflx(ji) = zdz1_snow(ji,1) * (cgrnd_snow(ji,1) + (dgrnd_snow(ji,1)-1.) * snowtemp(ji,1))
       snowcap(ji) = (zdz2_snow(ji,1) * dt_sechiba + dt_sechiba * (un - dgrnd_snow(ji,1)) * zdz1_snow(ji,1))
       z1_snow(ji) = lambda_snow(ji) * (un - dgrnd_snow(ji,1)) + un 
       snowcap(ji) = snowcap(ji) / z1_snow(ji)
       snowflx(ji) = snowflx(ji) + &
            snowcap(ji) * (snowtemp(ji,1) * z1_snow(ji) - lambda_snow(ji) * cgrnd_snow(ji,1) - temp_sol_new(ji)) / dt_sechiba
    ENDDO

 
    !! Computation of the apparent ground heat flux (> towards the soil) and
    !! apparent surface heat capacity, used at the next timestep by enerbil to
    !! compute the surface temperature.
    DO ji = 1,kjpindex
      soilflx_nosnow(ji) = zdz1(ji,1) * (cgrnd(ji,1) + (dgrnd(ji,1)-1.) * ptn(ji,1))
      soilcap_nosnow(ji) = (zdz2(ji,1) * dt_sechiba + dt_sechiba * (un - dgrnd(ji,1)) * zdz1(ji,1))
      z1(ji) = lambda * (un - dgrnd(ji,1)) + un
      soilcap_nosnow(ji) = soilcap_nosnow(ji) / z1(ji)
      soilflx_nosnow(ji) = soilflx_nosnow(ji) + &
         & soilcap_nosnow(ji) * (ptn(ji,1) * z1(ji) - lambda * cgrnd(ji,1) - temp_sol_new(ji)) / dt_sechiba 
    ENDDO

    !! Add snow fraction
    ! Using an effective heat capacity and heat flux by a simple pondering of snow and soil fraction
    DO ji = 1, kjpindex
       soilcap(ji) = snowcap(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+    & ! weights related to snow cover fraction on vegetation  
            soilcap_nosnow(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ & ! weights related to SCF on nobio
            soilcap_nosnow(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow fraction
       soilflx(ji) = snowflx(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+    & ! weights related to snow cover fraction on vegetation  
            soilflx_nosnow(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ & ! weights related to SCF on nobio
            soilflx_nosnow(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow fraction
    ENDDO

    ! Total snow flux (including snow on vegetated and bare soil and nobio areas)
    DO ji = 1, kjpindex
    	snowflxtot(ji) = snowflx(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + &
          soilflx_nosnow(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)
    ENDDO
    CALL xios_orchidee_send_field("snowflxtot",snowflxtot(:))

    IF (printlev>=3) WRITE (numout,*) ' thermosoil_coef done '

  END SUBROUTINE thermosoil_coef
 
 
!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_profile
!!
!>\BRIEF        In this routine solves the numerical soil thermal scheme, ie calculates the new soil temperature profile; 
!! 
!!
!! DESCRIPTION	: The calculation of the new soil temperature profile is based on
!! the cgrnd and dgrnd values from the previous timestep and the surface temperature Ts aka temp_sol_new. (see detailed
!! explanation in the header of the thermosoil module or in the reference).\n
!!                              T(k+1)=cgrnd(k)+dgrnd(k)*T(k)\n
!!                                      -- EQ1 --\n
!!                           Ts=(1+lambda)*T(1) -lambda*T(2)\n 
!!                                      -- EQ2--\n
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): ptn (soil temperature profile on the thermal axis), 
!!                          stempdiag (soil temperature profile on the diagnostic axis)
!!
!! REFERENCE(S) :
!! - Hourdin, F. (1992). Study and numerical simulation of the general circulation of planetary atmospheres,
!! Ph.D. thesis, Paris VII University. Remark: the part of F. Hourdin's PhD thesis relative to the thermal
!! integration scheme has been scanned and is provided along with the documentation, with name : 
!! Hourdin_1992_PhD_thermal_scheme.pdf
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

 SUBROUTINE thermosoil_profile (kjpindex,      temp_sol_new,                   &
                                frac_snow_veg, frac_snow_nobio, totfrac_nobio, &
                                ptn,           stempdiag,       snowtemp,      &
                                cgrnd_snow,    dgrnd_snow)

  !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                               :: kjpindex       !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_sol_new   !! Surface temperature at the present time-step 
                                                                               !! @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT(in)             :: frac_snow_veg  !! Snow cover fraction on vegeted area
    REAL(r_std),DIMENSION (kjpindex,nnobio), INTENT(in)      :: frac_snow_nobio!! Snow cover fraction on non-vegeted area
    REAL(r_std),DIMENSION (kjpindex),INTENT(in)              :: totfrac_nobio  !! Total fraction of continental ice+lakes+cities+...
                                                                               !! (unitless,0-1)
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: snowtemp       !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: cgrnd_snow     !! Integration coefficient for snow numerical scheme
    REAL(r_std),DIMENSION (kjpindex,nsnow), INTENT(in)       :: dgrnd_snow     !! Integration coefficient for snow numerical scheme
 
    !! 0.3 Modified variables

 
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex,ngrnd), INTENT (out)     :: ptn            !! vertically discretized soil temperatures 
                                                                               !! @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (out)      :: stempdiag      !! diagnostic temperature profile 
                                                                               !! @tex ($K$) @endtex

    !! 0.4 Local variables

    INTEGER(i_std)                                           :: ji, jg
    REAL(r_std)                                              :: temp_sol_eff   !! effective surface temperature including snow and soil
     
!_ ================================================================================================================================

  !! 1. Computes the soil temperatures ptn.

    !! 1.1. ptn(jg=1) using EQ1 and EQ2
    DO ji = 1,kjpindex

       ! Using an effective surface temperature by a simple pondering 
       temp_sol_eff=snowtemp(ji,nsnow)*frac_snow_veg(ji)*(1-totfrac_nobio(ji))+ &      ! weights related to snow cover fraction on vegetation  
            temp_sol_new(ji)*SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji)+ &           ! weights related to SCF on nobio
            temp_sol_new(ji)*(1-(frac_snow_veg(ji)*(1-totfrac_nobio(ji))+SUM(frac_snow_nobio(ji,:))*totfrac_nobio(ji))) ! weights related to non snow fraction
       ! Soil temperature calculation with explicit snow if there is snow on the ground
       ptn(ji,1) = cgrnd_snow(ji,nsnow) + dgrnd_snow(ji,nsnow) * temp_sol_eff
    ENDDO

    !! 1.2. ptn(jg=2:ngrnd) using EQ1.
    DO jg = 1,ngrnd-1
      DO ji = 1,kjpindex
        ptn(ji,jg+1) = cgrnd(ji,jg) + dgrnd(ji,jg) * ptn(ji,jg)
      ENDDO
    ENDDO

    !! 2. Assigne the soil temperature to the output variable. It is already on the right axis. 
    stempdiag(:,:) = ptn(:,1:nslm)

    IF (printlev>=3) WRITE (numout,*) ' thermosoil_profile done '

  END SUBROUTINE thermosoil_profile

!================================================================================================================================
!! SUBROUTINE   : thermosoil_cond
!!
!>\BRIEF          Calculate soil thermal conductivity.  
!!
!! DESCRIPTION  : This routine computes soil thermal conductivity
!!                Code introduced from NOAH LSM. 
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): cnd
!!
!! REFERENCE(S) :
!!    Farouki, O.T.,1986: Thermal Properties of Soils. Series on Rock
!!            and Soil Mechanics, Vol. 11, Trans Tech, 136 PP.
!!    Johansen, O., 1975: Thermal Conductivity of Soils. Ph.D. Thesis,
!!            University of Trondheim,
!!    Peters-Lidard, C. D., Blackburn, E., Liang, X., & Wood, E. F.,
!!            1998: The effect of soil thermal conductivity 
!!            Parameterization on Surface Energy fluxes
!!            and Temperatures. J. of The Atmospheric Sciences,
!!            Vol. 55, pp. 1209-1224.
!! Modify histroy:
!!
!! FLOWCHART    : None
!! \n
!_
!================================================================================================================================

  SUBROUTINE thermosoil_cond (kjpindex, njsc, mcs, smc, qz, sh2o, cnd, veget_max)

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                                 :: kjpindex      !! Domain size (unitless)
    INTEGER(i_std), DIMENSION (kjpindex), INTENT (in)          :: njsc          !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std), DIMENSION (kjpindex), INTENT(IN)              :: mcs           !! Saturated moisture content (m3/m3)
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(IN)        :: smc           !! Volumetric Soil Moisture Content (m3/m3)
    REAL(r_std), DIMENSION (nscm), INTENT(IN)                  :: qz            !! Quartz Content (Soil Type Dependent) (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(IN)        :: sh2o          !! Unfrozen Soil Moisture Content; Frozen Soil Moisture = smc - sh2o
    REAL(r_std), DIMENSION (kjpindex,nvm2), INTENT (in)        :: veget_max       !! Fraction of PFT (unitless,0-1)  4
    !! 0.2 Output variables
    REAL(r_std), DIMENSION (kjpindex,ngrnd), INTENT(OUT)       :: cnd           !! Soil Thermal Conductivity (W/m/k)
    
    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: ake           !! Kersten Number (unitless)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: thksat        !! Saturated Thermal Conductivity (W/m/k)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: satratio      !! Degree of Saturation (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: xu            !! Unfrozen Volume For Saturation (0-1)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)                    :: xunfroz       !! Unfrozon Volume Fraction (0-1)
    REAL(r_std)                                                :: thko          !! Thermal Conductivity for Other Ssoil Components (W/m/k)
    REAL(r_std)                                                :: gammd         !! Dry Dendity (kg/m3)
    REAL(r_std)                                                :: thkdry        !! Dry Thermal Conductivity (W/m/k)
    REAL(r_std)                                                :: thks          !! Thermal Conductivity for the Solids Combined (Quartz + Other) (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKICE = 2.2  !! Ice Thermal Conductivity (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKQTZ = 7.7  !! Thermal Conductivity for Quartz (W/m/k)
    REAL(r_std), PARAMETER                                     :: THKW = 0.57   !! Water Thermal Conductivity (W/m/k)
    INTEGER(i_std)                                             :: ji, jg, jst
    
!_================================================================================================================================
    
    !! 1. Dry and Saturated Thermal Conductivity.
   
    DO ji = 1,kjpindex
      jst = njsc(ji)

      !! 1.1. Dry density (Kg/m3) and Dry thermal conductivity (W.M-1.K-1)
      gammd = (1. - mcs(ji))*2700.
      thkdry = (0.135* gammd+ 64.7)/ (2700. - 0.947* gammd)

      !! 1.2. thermal conductivity of "other" soil components
      IF (qz(jst) > 0.2) THEN
         thko = 2.0
      ELSEIF (qz(jst) <= 0.2) THEN
         thko = 3.0
      ENDIF

      !! 1.3. Thermal conductivity of solids
      thks = (THKQTZ ** qz(jst))* (thko ** (1. - qz(jst)))

      DO jg = 1,ngrnd      
        !! 1.4. saturation ratio
        satratio(ji,jg) = smc(ji,jg) / mcs(ji)
    
        !! 1.5. Saturated Thermal Conductivity (thksat)
        IF ( smc(ji,jg) > min_sechiba ) THEN
           xunfroz(ji,jg) = sh2o(ji,jg) / smc(ji,jg)   ! Unfrozen Fraction (From i.e., 100%Liquid, to 0. (100% Frozen))
           xu(ji,jg) = xunfroz(ji,jg) * mcs(ji)  ! Unfrozen volume for saturation (porosity*xunfroz)
           thksat(ji,jg) = thks ** (1. - mcs(ji))* THKICE ** (mcs(ji) - xu(ji,jg))* THKW ** (xu(ji,jg))
        ELSE
           ! this value will not be used since ake=0 for this case
           thksat(ji,jg)=0 
        END IF
      END DO ! DO jg = 1,ngrnd

      !! 2. Kersten Number (ake)
      DO jg = 1,ngrnd
        IF ( (sh2o(ji,jg) + 0.0005) <  smc(ji,jg) ) THEN
          ! Frozen
          ake(ji,jg) = satratio(ji,jg)
        ELSE
          ! Unfrozen
          ! Eq 11 in Peters-Lidard et al., 1998
          IF ( satratio(ji,jg) >  0.1 ) THEN
            IF (jst < 4 )  THEN
                ! Coarse 
                ake(ji,jg) = 0.7 * LOG10 (SATRATIO(ji,jg)) + 1.0
            ELSE
                ! Fine 
                ake(ji,jg) = LOG10 (satratio(ji,jg)) + 1.0
            ENDIF
          ELSEIF ( satratio(ji,jg) >  0.05 .AND. satratio(ji,jg) <=  0.1 ) THEN
            IF (jst < 4 )  THEN
                ! Coarse 
                ake(ji,jg) = 0.7 * LOG10 (satratio(ji,jg)) + 1.0
            ELSE
                ! Fine 
                ake(ji,jg) = 0.0
            ENDIF
          ELSE
            ake(ji,jg) = 0.0  ! use k = kdry
          END IF
        END IF
      END DO ! DO jg = 1,ngrnd

      !! 3. Thermal conductivity (cnd)
      DO jg = 1,ngrnd
        cnd(ji,jg) = ake(ji,jg) * (thksat(ji,jg) - thkdry) + thkdry
      END DO ! DO jg = 1,ngrnd

    END DO !DO ji = 1,kjpindex


    DO ji = 1,kjpindex
     DO jg = 1,ngrnd 
     
     IF (DO_URBAN_HEAT_CAPA_CONDUCT) THEN 
     IF (veget_max(ji,16) > 0.5) THEN 
       cnd(ji,jg) = 3.24
     ENDIF
     ENDIF
     
     END DO
    END DO

  END SUBROUTINE thermosoil_cond


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_humlev
!!
!>\BRIEF           Interpolate variables from the hydrology layers to the thermodynamic layers
!!
!! DESCRIPTION  :  Interpolate the volumetric soil moisture content from the node to the interface of the layer. 
!!                 The values for the deep layers in thermosoil where hydrology is not existing are constant. 
!!                 No interpolation is needed for the total soil moisture content and for the soil saturation degree.
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S): mc_layt, mcl_layt, tmc_layt, shum_ngrnd_perma
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
  SUBROUTINE thermosoil_humlev(kjpindex, shumdiag_perma, mc_layh, mcl_layh, tmc_layh)
  
  !! 0. Variables and parameter declaration

    !! 0.1 Input variables
 
    INTEGER(i_std), INTENT(in)                            :: kjpindex       !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: shumdiag_perma !! Soil saturation degree on the diagnostic axis (0-1, unitless)
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mc_layh        !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid+ice) [m/s]
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: mcl_layh       !! Volumetric soil moisture content for each layer in hydrol at nodes(liquid) [m/s]
    REAL(r_std),DIMENSION (kjpindex,nslm), INTENT (in)    :: tmc_layh       !! Total soil moisture content for each layer in hydrol(liquid+ice) [mm]
    
    !! 0.2 Output variables

    !! 0.3 Modified variables

    !! 0.4 Local variables
    INTEGER(i_std)                                       :: ji, jd

!_ ================================================================================================================================

    IF (printlev >= 4) WRITE(numout,*) 'Start thermosoil_humlev' 

    ! The values for the deep layers in thermosoil where hydrology is not existing are constant. 
    ! For exemple if thermosoil uses 8m, and hydrol uses 2m vertical discretization,
    ! the values between 2m and 8m are constant.
    ! The moisture computed in hydrol is at the nodes (except for the
    ! top and bottom layer which are at interfaces)
    ! A linear interpolation is applied to obtain the moisture values at
    ! the interfaces (mc_layt), from the mc_layh at the nodes

    DO ji=1,kjpindex
       DO jd = 1, nslm
          IF(jd == 1) THEN ! the moisture at the 1st interface mc_layh(1) is at the surface, no interpolation
             mc_layt(ji,jd) = mc_layh(ji,jd)
             mcl_layt(ji,jd) = mcl_layh(ji,jd)
          ELSEIF(jd == 2) THEN  !! the mc_layt at the 2nd interface is interpolated using mc_layh(1) at surface and mc_layh(2) at the node
             mc_layt(ji, jd) = mc_layh(ji,jd-1)*(znt(jd)-zlt(jd-1))/(znt(jd)-0.0) + &
                  mc_layh(ji, jd)*(zlt(jd-1)-0.0)/(znt(jd)-0.0)
             mcl_layt(ji, jd) = mcl_layh(ji,jd-1)*(znt(jd)-zlt(jd-1))/(znt(jd)-0.0) + &
                  mcl_layh(ji, jd)*(zlt(jd-1)-0.0)/(znt(jd)-0.0)
          ELSEIF(jd == nslm) THEN ! the mc_layt at the nslm interface is interpolated using mc_layh(nslm) and mc_layh(nslm-1)
             mc_layt(ji, jd) = mc_layh(ji,jd-1)*(zlt(jd)-zlt(jd-1))/(zlt(jd)-znt(jd-1))  + &
                  mc_layh(ji,jd)*(zlt(jd-1)-znt(jd-1))/(zlt(jd)-znt(jd-1))
             mcl_layt(ji, jd) = mcl_layh(ji,jd-1)*(zlt(jd)-zlt(jd-1))/(zlt(jd)-znt(jd-1))  + &
                  mcl_layh(ji,jd)*(zlt(jd-1)-znt(jd-1))/(zlt(jd)-znt(jd-1))
          ELSE ! the mc_layt at the other interfaces are interpolated using mc_layh at adjacent nodes.
             mc_layt(ji, jd) = mc_layh(ji, jd-1)*(1-dz5(jd-1)) + mc_layh(ji,jd)*dz5(jd-1)
             mcl_layt(ji, jd) = mcl_layh(ji, jd-1)*(1-dz5(jd-1)) + mcl_layh(ji,jd)*dz5(jd-1)
          ENDIF

          shum_ngrnd_perma(ji,jd) = shumdiag_perma(ji,jd)
          tmc_layt(ji,jd) = tmc_layh(ji,jd)
       ENDDO
       
       ! The deep layers in thermosoil where hydro is not existing
       DO jd = nslm+1, ngrnd
          shum_ngrnd_perma(ji,jd) = shumdiag_perma(ji,nslm)
          mc_layt(ji,jd) = mc_layh(ji,nslm)
          mcl_layt(ji,jd) = mcl_layh(ji,nslm)
          tmc_layt(ji,jd) = tmc_layh(ji,nslm)/dlt(nslm) *dlt(jd)
       ENDDO
    ENDDO

    IF (printlev >= 4) WRITE(numout,*) 'thermosoil_humlev done' 

  END SUBROUTINE thermosoil_humlev


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_energy_diag
!!
!>\BRIEF         Calculate diagnostics 
!!
!! DESCRIPTION  : Calculate diagnostic variables coldcont_incr and coldcont_incr
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : 
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_energy_diag(kjpindex, temp_sol_new, soilcap)
  
   !! 0. Variables and parameter declaration

    !! 0.1 Input variables

    INTEGER(i_std), INTENT(in)                     :: kjpindex     !! Domain size (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: temp_sol_new !! Surface temperature at the present time-step, Ts 
                                                                   !! @tex ($K$) @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: soilcap      !! Apparent surface heat capacity 
                                                                   !! @tex ($J m^{-2} K^{-1}$) @endtex, 
                                                                   !! see eq. A29 of F. Hourdin\'s PhD thesis.
    
    !! 0.2 Output variables

    !! 0.3 Modified variables
    
    !! 0.4 Local variables
    
    INTEGER(i_std)                                 :: ji, jg
!_ ================================================================================================================================
   
    !  Sum up the energy content of all layers in the soil.
    DO ji = 1, kjpindex
   
       IF (pcapa_en(ji,1) .LE. sn_capa) THEN
          
          ! Verify the energy conservation in the surface layer
          coldcont_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          surfheat_incr(ji) = zero
       ELSE
          
          ! Verify the energy conservation in the surface layer
          surfheat_incr(ji) = soilcap(ji) * (temp_sol_new(ji) - temp_sol_beg(ji))
          coldcont_incr(ji) = zero
       ENDIF
    ENDDO
    
    ! Save temp_sol_new to be used at next timestep
    temp_sol_beg(:)   = temp_sol_new(:)

  END SUBROUTINE thermosoil_energy_diag



!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_readjust
!!
!>\BRIEF        
!!
!! DESCRIPTION	: Energy conservation : Correction to make sure that the same latent heat is released and 
!!                consumed during freezing and thawing  
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): ptn (soil temperature profile on the thermal axis), 
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_readjust(kjpindex, ptn)

   !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                             :: kjpindex
    
    !! 0.2 Modified variables
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(inout)    :: ptn

    !! 0.3 Local variables
    INTEGER(i_std)  :: ji, jg
    INTEGER(i_std)  :: lev3m  !! Closest interface level to 3m
    REAL(r_std) :: ptn_tmp

    ! The energy is spread over the layers down to approximatly 3m
    ! Find the closest level to 3m. It can be below or above 3m.
    lev3m=MINLOC(ABS(zlt(:)-3.0),dim=1)
    IF (printlev >= 3) WRITE(numout,*) 'In thermosoil_adjust: lev3m=',lev3m, ' zlt(lev3m)=', zlt(lev3m)
    
    DO jg=1, ngrnd
       DO ji=1, kjpindex
          ! All soil latent energy is put into e_soil_lat(ji)
          ! because the variable soil layers make it difficult to keep track of all
          ! layers in this version
          ! NOTE : pcapa has unit J/K/m3 and pcappa_supp has J/K
          e_soil_lat(ji)=e_soil_lat(ji)+pcappa_supp(ji,jg)*(ptn(ji,jg)-ptn_beg(ji,jg))
       END DO
    END DO

   DO ji=1, kjpindex
      IF (e_soil_lat(ji).GT.min_sechiba.AND.MINVAL(ptn(ji,:)).GT.ZeroCelsius+fr_dT/2.) THEN
         ! The soil is thawed: we spread the excess of energy over the uppermost lev3m levels
         ! Here we increase the temperatures
         DO jg=1, lev3m
            ptn_tmp=ptn(ji,jg)
            
            ptn(ji,jg)=ptn(ji,jg)+MIN(e_soil_lat(ji)/pcapa(ji,jg)/zlt(lev3m), 0.5)
            e_soil_lat(ji)=e_soil_lat(ji)-(ptn(ji,jg)-ptn_tmp)*pcapa(ji,jg)*dlt(jg)
         ENDDO
      ELSE IF (e_soil_lat(ji).LT.-min_sechiba.AND.MINVAL(ptn(ji,:)).GT.ZeroCelsius+fr_dT/2.) THEN
         ! The soil is thawed
         ! Here we decrease the temperatures
         DO jg=1, lev3m
            ptn_tmp=ptn(ji,jg)
            ptn(ji,jg)=MAX(ZeroCelsius+fr_dT/2., ptn_tmp+e_soil_lat(ji)/pcapa(ji,jg)/zlt(lev3m))
            e_soil_lat(ji)=e_soil_lat(ji)+(ptn_tmp-ptn(ji,jg))*pcapa(ji,jg)*dlt(jg)
         END DO
      END IF
   END DO

  END SUBROUTINE thermosoil_readjust
   
!-------------------------------------------------------------------



!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_getdiff
!!
!>\BRIEF          Computes soil and snow heat capacity and conductivity    
!!
!! DESCRIPTION	: Computation of the soil thermal properties; snow properties are also accounted for
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): 
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

  SUBROUTINE thermosoil_getdiff( kjpindex, snow, ptn, mcs, njsc, snowrho, snowtemp, pb, veget_max)

   !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std),INTENT(in)				:: kjpindex
    REAL(r_std),DIMENSION(kjpindex),INTENT (in)	        :: snow       !! Snow mass
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc       !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: mcs        !! Saturated moisture content (m3/m3)
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in) :: snowrho    !! Snow density
    REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in) :: snowtemp   !! Snow temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: pb         !! Surface pressure (hPa)
    REAL(r_std),DIMENSION(kjpindex,ngrnd),INTENT(in)	  :: ptn        !! Soil temperature profile
    REAL(r_std), DIMENSION (kjpindex,nvm2), INTENT (in) :: veget_max  !! Fraction of PFT (unitless,0-1)  5
    !! 0.3 Local variables
    REAL						:: xx         !! Unfrozen fraction of the soil
    REAL(r_std), DIMENSION(kjpindex)             	:: snow_h
    REAL(r_std), DIMENSION(kjpindex,ngrnd) 		:: zx1, zx2
    INTEGER						:: ji,jg
    INTEGER                                             :: jst
    REAL(r_std), DIMENSION (kjpindex,ngrnd)             :: pcapa_tmp  !! soil heat capacity (J/m3/K)
    REAL(r_std), DIMENSION (kjpindex,ngrnd)             :: pcapa_spec !! SPECIFIC soil heat capacity (J/kg/K)
    REAL(r_std)                                         :: rho_tot    !! Soil density (kg/m3)

    pcapa_tmp(:,:) = 0.0

    !! Computes soil heat capacity and conductivity
    DO ji = 1,kjpindex

      ! Since explicitsnow module is implemented set zx1=0 and zx2=1
      zx1(ji,:) = 0.
      zx2(ji,:) = 1.

      DO jg = 1, ngrnd
         jst = njsc(ji)
         pcapa_tmp(ji, jg) = so_capa_dry(jst) * (1-mcs(ji)) + water_capa * tmc_layt(ji,jg)/mille/dlt(jg)
         !
         ! 2. Calculate volumetric heat capacity with allowance for permafrost
         ! 2.1. soil heat capacity depending on temperature and humidity
         ! For SP6MIP we also diagnose a specific heat capacity (pcapa_spec), 
         ! which requires to calculate the total density of the soil (rho_tot), always >> 0
	 
         IF (ptn(ji,jg) .LT. ZeroCelsius-fr_dT/2.) THEN
	    ! frozen soil
	    profil_froz(ji,jg) = 1.
 	    pcappa_supp(ji,jg)= 0.
            pcapa(ji, jg) = so_capa_dry(jst) * (1-mcs(ji)) + so_capa_ice(ji) * tmc_layt(ji,jg) / mille / dlt(jg)
             IF (DO_URBAN_HEAT_CAPA_CONDUCT) THEN
             IF (veget_max(ji,16) > 0.5) THEN 
               pcapa(ji,jg) = 1890000.
             ENDIF
             ENDIF
            rho_tot = rho_soil * (1-mcs(ji)) + rho_ice * tmc_layt(ji,jg) / mille / dlt(jg) 
            pcapa_spec(ji, jg) = pcapa(ji, jg) / rho_tot

     	 ELSEIF (ptn(ji,jg) .GT. ZeroCelsius+fr_dT/2.) THEN
	    ! unfrozen soil	 
            pcapa(ji, jg) = pcapa_tmp(ji, jg)
             IF (DO_URBAN_HEAT_CAPA_CONDUCT) THEN
             IF (veget_max(ji,16) > 0.5) THEN 
               pcapa(ji,jg) = 1890000.
             ENDIF
             ENDIF
	    profil_froz(ji,jg) = 0.
 	    pcappa_supp(ji,jg)= 0.
            rho_tot = rho_soil * (1-mcs(ji)) + rho_water * tmc_layt(ji,jg)/mille/dlt(jg)
            pcapa_spec(ji, jg) = pcapa(ji, jg) / rho_tot
     	 ELSE
     	   ! xx is the unfrozen fraction of soil water	   	   
     	   xx = (ptn(ji,jg)-(ZeroCelsius-fr_dT/2.)) / fr_dT
           profil_froz(ji,jg) = (1. - xx)

    	   IF (ok_freeze_thaw_latent_heat) THEN
              pcapa(ji, jg) = so_capa_dry(jst) * (1-mcs(ji)) + &
                water_capa * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                so_capa_ice(ji) * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx) + &
                shum_ngrnd_perma(ji,jg)*mcs(ji)*lhf*rho_water/fr_dT
                 IF (DO_URBAN_HEAT_CAPA_CONDUCT) THEN
                 IF (veget_max(ji,16) > 0.5) THEN 
                   pcapa(ji,jg) = 1890000.
                 ENDIF
                 ENDIF
	   ELSE
              pcapa(ji, jg) = so_capa_dry(jst) * (1-mcs(ji)) + &
                water_capa * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                so_capa_ice(ji) * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx)
                 IF (DO_URBAN_HEAT_CAPA_CONDUCT) THEN
                 IF (veget_max(ji,16) > 0.5) THEN 
                   pcapa(ji,jg) = 1890000.
                 ENDIF
                 ENDIF
           ENDIF

           rho_tot =  rho_soil* (1-mcs(ji)) + &
                rho_water * tmc_layt(ji,jg)/mille / dlt(jg) * xx + &
                rho_ice * tmc_layt(ji,jg) / mille/dlt(jg) * (1.-xx)
           pcapa_spec(ji, jg) = pcapa(ji, jg) / rho_tot

 	   pcappa_supp(ji,jg)= shum_ngrnd_perma(ji,jg)*mcs(ji)*lhf*rho_water/fr_dT*zx2(ji,jg)*dlt(jg)

         ENDIF

         ! 
	 ! 2.2. Take into account the snow and soil fractions in the layer
	 !
         pcapa(ji,jg) = zx1(ji,jg) * sn_capa + zx2(ji,jg) * pcapa(ji,jg)

	 !
	 ! 2.3. Calculate the heat capacity for energy conservation check 
	 IF ( zx1(ji,jg).GT.0. ) THEN
            pcapa_en(ji,jg) = sn_capa
	 ELSE
            pcapa_en(ji,jg) = pcapa(ji,jg)
	 ENDIF

      END DO            
    ENDDO 

    ! Output the specific heat capcaity for SP-MIP
    CALL xios_orchidee_send_field("pcapa_spec",pcapa_spec)

    !
    ! 3. Calculate the heat conductivity with allowance for permafrost
    !
    IF (ok_freeze_thaw_latent_heat) THEN
    	CALL thermosoil_cond (kjpindex, njsc, mcs, mc_layt, QZ, mcl_layt*(1-profil_froz), pkappa, veget_max)
    ELSE
    	CALL thermosoil_cond (kjpindex, njsc, mcs, mc_layt, QZ, mcl_layt, pkappa, veget_max)
    ENDIF

    !! Computes snow heat capacity and conductivity    
    DO ji = 1,kjpindex
       pcapa_snow(ji,:) = snowrho(ji,:) * xci
       pkappa_snow(ji,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*snowrho(ji,:)*snowrho(ji,:)) +      &
            MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
            ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
    END DO
   
   END SUBROUTINE thermosoil_getdiff


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_getdiff_old_thermix_without_snow
!!
!>\BRIEF          Computes soil and snow heat capacity and conductivity    
!!
!! DESCRIPTION	: Calculations of soil and snow thermal properties without effect of freezing.
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S):
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================

    SUBROUTINE thermosoil_getdiff_old_thermix_without_snow( kjpindex, mcs, njsc, snowrho, snowtemp, pb, veget_max)

   !! 0. Variables and parameter declaration

    !! 0.1 Input variables
      INTEGER(i_std), INTENT(in) :: kjpindex
      INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: njsc     !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
      REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: mcs      !! Saturated moisture content (m3/m3)
      REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in) :: snowrho  !! Snow density
      REAL(r_std), DIMENSION (kjpindex,nsnow), INTENT(in) :: snowtemp !! Snow temperature (K)
      REAL(r_std),DIMENSION (kjpindex), INTENT (in)       :: pb       !! Surface pressure (hPa)
      REAL(r_std), DIMENSION (kjpindex,nvm2), INTENT (in) :: veget_max  !! Fraction of PFT (unitless,0-1)  6

    !! 0.1 Local variables
      INTEGER(i_std)    				  :: ji,jg, jst     !! Index
      REAL(r_std), DIMENSION (kjpindex,ngrnd)             :: pcapa_tmp      !! Soil heat capacity (J/m3/K) 

      !! Computes soil heat capacity and conductivity
      DO jg = 1,ngrnd
         DO ji = 1,kjpindex
            jst = njsc(ji)
            pcapa_tmp(ji, jg) = so_capa_dry(jst) * (1-mcs(ji)) + water_capa * tmc_layt(ji,jg)/mille/dlt(jg)
            pcapa(ji,jg) = pcapa_tmp(ji, jg)
             IF (DO_URBAN_HEAT_CAPA_CONDUCT) THEN
             IF (veget_max(ji,16) > 0.5) THEN 
               pcapa(ji,jg) = 1890000.
             ENDIF
             ENDIF
            pcapa_en(ji,jg) = pcapa_tmp(ji, jg)
         ENDDO
      ENDDO

      CALL thermosoil_cond (kjpindex, njsc, mcs, mc_layt, QZ, mcl_layt, pkappa, veget_max)

      IF (brk_flag == 1) THEN
        ! Bedrock flag is activated
        DO jg = ngrnd-1,ngrnd
          DO ji = 1,kjpindex
             pcapa(ji,jg) = brk_capa
             pcapa_en(ji,jg) = brk_capa 
             pkappa(ji,jg) = brk_cond
          ENDDO
        ENDDO
      ENDIF

    !! Computes snow heat capacity and conductivity
    DO ji = 1,kjpindex
        pcapa_snow(ji,:) = snowrho(ji,:) * xci
        pkappa_snow(ji,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*snowrho(ji,:)*snowrho(ji,:)) +      &
              MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(snowtemp(ji,:)+ &
              ZSNOWTHRMCOND_CVAP)))*(XP00/(pb(ji)*100.)))
    END DO

  END SUBROUTINE thermosoil_getdiff_old_thermix_without_snow


!! ================================================================================================================================
!! SUBROUTINE   : thermosoil_read_reftempfile
!!
!>\BRIEF          
!!
!! DESCRIPTION	: Read file with longterm soil temperature
!!                
!!
!! RECENT CHANGE(S) : None
!! 
!! MAIN OUTPUT VARIABLE(S): reftemp : Reference temerature
!!                          
!! REFERENCE(S) :
!!
!! FLOWCHART    : None 
!! \n 
!_ ================================================================================================================================
  SUBROUTINE thermosoil_read_reftempfile(kjpindex,lalo,reftemp)
    
    USE interpweight

    IMPLICIT NONE

    !! 0. Variables and parameter declaration

    !! 0.1 Input variables
    INTEGER(i_std), INTENT(in) :: kjpindex
    REAL(r_std), DIMENSION(kjpindex,2), INTENT(in) :: lalo

    !! 0.2 Output variables
    REAL(r_std), DIMENSION(kjpindex, ngrnd), INTENT(out) :: reftemp

    !! 0.3 Local variables
    INTEGER(i_std) :: ib
    CHARACTER(LEN=80) :: filename
    REAL(r_std),DIMENSION(kjpindex) :: reftemp_file                          !! Horizontal temperature field interpolated from file [C]
    INTEGER(i_std),DIMENSION(kjpindex,8) :: neighbours
    REAL(r_std)                                          :: vmin, vmax       !! min/max values to use for the
                                                                             !!   renormalization
    REAL(r_std), DIMENSION(kjpindex)                     :: areftemp         !! Availability of data for  the interpolation
    CHARACTER(LEN=80)                                    :: variablename     !! Variable to interpolate
                                                                             !!   the file
    CHARACTER(LEN=80)                                    :: lonname, latname !! lon, lat names in input file
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
                                                                             !!       (normalized by maskvals(3))
                                                                             !!   'var': mask values are taken from a 
                                                                             !!     variable inside the file (>0)
    REAL(r_std), DIMENSION(3)                            :: maskvals         !! values to use to mask (according to 
                                                                             !!   `maskingtype') 
    CHARACTER(LEN=250)                                   :: namemaskvar      !! name of the variable to use to mask 
    REAL(r_std)                                          :: reftemp_norefinf
    REAL(r_std)                                          :: reftemp_default  !! Default value

    !Config Key   = SOIL_REFTEMP_FILE
    !Config Desc  = File with climatological soil temperature
    !Config If    = READ_REFTEMP
    !Config Def   = reftemp.nc
    !Config Help  = 
    !Config Units = [FILE]
    filename = 'reftemp.nc'
    CALL getin_p('REFTEMP_FILE',filename)

    variablename = 'temperature'

    IF (printlev >= 1) WRITE(numout,*) "thermosoil_read_reftempfile: Read and interpolate file " &
         // TRIM(filename) //" for variable " //TRIM(variablename)

    IF (xios_interpolation) THEN

       CALL xios_orchidee_recv_field('reftemp_interp',reftemp_file)

       DO ib=1, kjpindex
           reftemp(ib,:) = reftemp_file(ib) + ZeroCelsius 
       END DO
       areftemp = 1.0
    ELSE


       ! For this case there are not types/categories. We have 'only' a continuos field
       ! Assigning values to vmin, vmax
       
       vmin = 0.
       vmax = 9999.

       !   For this file we do not need neightbours!
       neighbours = 0
       
       !! Variables for interpweight
       ! Type of calculation of cell fractions
       fractype = 'default'
       ! Name of the longitude and latitude in the input file
       lonname = 'nav_lon'
       latname = 'nav_lat'
       ! Default value when no value is get from input file
       reftemp_default = 1.
       ! Reference value when no value is get from input file
       reftemp_norefinf = 1.
       ! Should negative values be set to zero from input file?
       nonegative = .FALSE.
       ! Type of mask to apply to the input data (see header for more details)
       maskingtype = 'nomask'
       ! Values to use for the masking (here not used)
       maskvals = (/ undef_sechiba, undef_sechiba, undef_sechiba /)
       ! Name of the variable with the values for the mask in the input file (only if maskkingtype='var') (here not used)
       namemaskvar = ''
       
       CALL interpweight_2Dcont(kjpindex, 0, 0, lalo, resolution, neighbours,                            &
            contfrac, filename, variablename, lonname, latname, vmin, vmax, nonegative, maskingtype,        &
            maskvals, namemaskvar, -1, fractype, reftemp_default, reftemp_norefinf,                         &
            reftemp_file, areftemp)
       IF (printlev >= 5) WRITE(numout,*)'  thermosoil_read_reftempfile after interpweight_2Dcont'
       
       ! Copy reftemp_file temperature to all ground levels and transform into Kelvin
       DO ib=1, kjpindex
          reftemp(ib, :) = reftemp_file(ib)+ZeroCelsius
       END DO

    END IF

    ! Write diagnostics
    CALL xios_orchidee_send_field("interp_avail_areftemp",areftemp)
    CALL xios_orchidee_send_field("interp_diag_reftemp",reftemp_file)
    
  END SUBROUTINE thermosoil_read_reftempfile

END MODULE thermosoil
