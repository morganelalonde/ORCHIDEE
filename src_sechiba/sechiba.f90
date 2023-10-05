!  ==============================================================================================================================\n
!  MODULE 	: sechiba
! 
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        Structures the calculation of atmospheric and hydrological 
!! variables by calling diffuco_main, enerbil_main, hydrol_main,
!! condveg_main and thermosoil_main. Note that sechiba_main
!! calls slowproc_main and thus indirectly calculates the biogeochemical
!! processes as well.
!!
!!\n DESCRIPTION  : :: shumdiag, :: litterhumdiag and :: stempdiag :: ftempdiag are not 
!! saved in the restart file because at the first time step because they 
!! are recalculated. However, they must be saved as they are in slowproc 
!! which is called before the modules which calculate them.
!! 
!! RECENT CHANGE(S): November 2020: It is possible to define soil hydraulic parameters from maps,
!!                    as needed for the SP-MIP project (Tafasca Salma and Ducharne Agnes).
!!                    Here, it leads to declare and allocate global variables. 
!! 
!! REFERENCE(S) : None
!!   
!! SVN     :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE_2_2/ORCHIDEE/src_sechiba/sechiba.f90 $ 
!! $Date: 2022-07-20 13:09:05 +0200 (Wed, 20 Jul 2022) $
!! $Revision: 7710 $
!! \n
!_ ================================================================================================================================
 
MODULE sechiba
 
  USE ioipsl
  USE xios_orchidee
  USE constantes
  USE time, ONLY : one_day, dt_sechiba
  USE constantes_soil
  USE pft_parameters
  USE grid
  USE diffuco
  USE condveg
  USE enerbil
  USE hydrol
  USE thermosoil
  USE sechiba_io_p
  USE slowproc
  USE routing_wrapper
  USE ioipsl_para
  USE chemistry

  IMPLICIT NONE

  PRIVATE
  PUBLIC sechiba_main, sechiba_initialize, sechiba_clear, sechiba_interface_orchidee_inca, sechiba_xios_initialize

  INTEGER(i_std), SAVE                             :: printlev_loc   !! local printlev for this module
!$OMP THREADPRIVATE(printlev_loc)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexveg       !! indexing array for the 3D fields of vegetation
!$OMP THREADPRIVATE(indexveg)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexlai       !! indexing array for the 3D fields of vegetation
!$OMP THREADPRIVATE(indexlai)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexnobio     !! indexing array for the 3D fields of other surfaces (ice,
                                                                     !! lakes, ...)
!$OMP THREADPRIVATE(indexnobio)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexsoil      !! indexing array for the 3D fields of soil types (kjpindex*nstm)
!$OMP THREADPRIVATE(indexsoil)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexgrnd      !! indexing array for the 3D ground heat profiles (kjpindex*ngrnd)
!$OMP THREADPRIVATE(indexgrnd)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexlayer     !! indexing array for the 3D fields of soil layers in CWRR (kjpindex*nslm)
!$OMP THREADPRIVATE(indexlayer)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexnslm      !! indexing array for the 3D fields of diagnostic soil layers (kjpindex*nslm)
!$OMP THREADPRIVATE(indexnslm)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexalb       !! indexing array for the 2 fields of albedo
!$OMP THREADPRIVATE(indexalb)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: indexsnow      !! indexing array for the 3D fields snow layers
!$OMP THREADPRIVATE(indexsnow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: veget          !! Fraction of vegetation type (unitless, 0-1)       
!$OMP THREADPRIVATE(veget)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: veget_max      !! Max. fraction of vegetation type (LAI -> infty, unitless)
!$OMP THREADPRIVATE(veget_max)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: tot_bare_soil  !! Total evaporating bare soil fraction 
!$OMP THREADPRIVATE(tot_bare_soil)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: height         !! Vegetation Height (m)
!$OMP THREADPRIVATE(height)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: totfrac_nobio  !! Total fraction of continental ice+lakes+cities+...
                                                                     !! (unitless, 0-1)
!$OMP THREADPRIVATE(totfrac_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: floodout       !! Flow out of floodplains from hydrol
!$OMP THREADPRIVATE(floodout)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: runoff         !! Surface runoff calculated by hydrol
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(runoff)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: drainage       !! Deep drainage calculatedd by hydrol
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(drainage)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: returnflow     !! Water flow from lakes and swamps which returns to 
                                                                     !! the grid box @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(returnflow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: reinfiltration !! Routed water which returns into the soil
!$OMP THREADPRIVATE(reinfiltration)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: irrigation     !! Irrigation flux taken from the routing reservoirs and 
                                                                     !! being put into the upper layers of the soil 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(irrigation)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: emis           !! Surface emissivity (unitless)
!$OMP THREADPRIVATE(emis)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: z0h           !! Surface roughness for heat (m)
!$OMP THREADPRIVATE(z0h)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: z0m           !! Surface roughness for momentum (m)
!$OMP THREADPRIVATE(z0m)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: roughheight    !! Effective height for roughness (m)
!$OMP THREADPRIVATE(roughheight)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: reinf_slope      !! slope coefficient (reinfiltration)
!$OMP THREADPRIVATE(reinf_slope)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: reinf_slope_soil !! slope coefficient (reinfiltration) per soil tile
!$OMP THREADPRIVATE(reinf_slope_soil)


!salma: adding soil hydraulic params
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: ks              !! Saturated soil conductivity (mm d^{-1})
!$OMP THREADPRIVATE(ks)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: nvan            !! Van Genushten n parameter (unitless)
!$OMP THREADPRIVATE(nvan)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: avan            !! Van Genushten alpha parameter (mm ^{-1})
!$OMP THREADPRIVATE(avan)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcr             !! Residual soil moisture (m^{3} m^{-3})
!$OMP THREADPRIVATE(mcr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcs             !! Saturated soil moisture (m^{3} m^{-3})
!$OMP THREADPRIVATE(mcs)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcfc            !! Volumetric water content at field capacity (m^{3} m^{-3})
!$OMP THREADPRIVATE(mcfc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)     :: mcw             !! Volumetric water content at wilting point (m^{3} m^{-3})
!$OMP THREADPRIVATE(mcw)
!end salma:adding soil hydraulic params


  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag       !! Mean relative soil moisture in the different levels used 
                                                                     !! by thermosoil.f90 (unitless, 0-1)
!$OMP THREADPRIVATE(shumdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: shumdiag_perma !! Saturation degree of the soil 
!$OMP THREADPRIVATE(shumdiag_perma)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: k_litt         !! litter cond.
!$OMP THREADPRIVATE(k_litt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: litterhumdiag  !! Litter dryness factor (unitless, 0-1)
!$OMP THREADPRIVATE(litterhumdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: stempdiag      !! Temperature which controls canopy evolution (K)
!$OMP THREADPRIVATE(stempdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: ftempdiag      !! Temperature over the full soil column for river temperature (K)
!$OMP THREADPRIVATE(ftempdiag)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: qsintveg       !! Water on vegetation due to interception 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(qsintveg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vbeta2         !! Interception resistance (unitless,0-1)
!$OMP THREADPRIVATE(vbeta2)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vbeta3         !! Vegetation resistance (unitless,0-1)
!$OMP THREADPRIVATE(vbeta3)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vbeta3pot      !! Potential vegetation resistance
!$OMP THREADPRIVATE(vbeta3pot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: gsmean         !! Mean stomatal conductance for CO2 (mol m-2 s-1) 
!$OMP THREADPRIVATE(gsmean) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: cimean         !! STOMATE: mean intercellular CO2 concentration (ppm)
!$OMP THREADPRIVATE(cimean)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vevapwet       !! Interception loss over each PFT 
                                                                     !! @tex $(kg m^{-2} days^{-1})$ @endtex
!$OMP THREADPRIVATE(vevapwet)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: transpir       !! Transpiration @tex $(kg m^{-2} days^{-1})$ @endtex
!$OMP THREADPRIVATE(transpir)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: transpot       !! Potential Transpiration (needed for irrigation)
!$OMP THREADPRIVATE(transpot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: qsintmax       !! Maximum amount of water in the canopy interception 
                                                                     !! reservoir @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(qsintmax)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: rveget         !! Surface resistance for the vegetation 
                                                                     !! @tex $(s m^{-1})$ @endtex
!$OMP THREADPRIVATE(rveget)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: rstruct        !! Vegetation structural resistance
!$OMP THREADPRIVATE(rstruct)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: snow_nobio     !! Snow mass of non-vegetative surfaces 
                                                                     !! @tex $(kg m^{-2})$ @endtex
!$OMP THREADPRIVATE(snow_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: snow_nobio_age !! Snow age on non-vegetative surfaces (days)
!$OMP THREADPRIVATE(snow_nobio_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: frac_nobio     !! Fraction of non-vegetative surfaces (continental ice, 
                                                                     !! lakes, ...) (unitless, 0-1)
!$OMP THREADPRIVATE(frac_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: assim_param    !! min+max+opt temps, vcmax, vjmax for photosynthesis
!$OMP THREADPRIVATE(assim_param)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: lai            !! Surface foliaire
!$OMP THREADPRIVATE(lai)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: gpp            !! STOMATE: GPP. gC/m**2 of total area
!$OMP THREADPRIVATE(gpp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)     :: temp_growth    !! Growth temperature (C) - Is equal to t2m_month 
!$OMP THREADPRIVATE(temp_growth) 
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: humrel         !! Relative humidity
!$OMP THREADPRIVATE(humrel)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: vegstress      !! Vegetation moisture stress (only for vegetation growth)
!$OMP THREADPRIVATE(vegstress)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:,:):: frac_age       !! Age efficacity from STOMATE for isoprene 
!$OMP THREADPRIVATE(frac_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: soiltile       !! Fraction of each soil tile (0-1, unitless)
!$OMP THREADPRIVATE(soiltile)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: fraclut        !! Fraction of each landuse tile (0-1, unitless)
!$OMP THREADPRIVATE(fraclut)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: nwdFraclut     !! Fraction of non-woody vegetation in each landuse tile (0-1, unitless)
!$OMP THREADPRIVATE(nwdFraclut)
  INTEGER(i_std), ALLOCATABLE, SAVE, DIMENSION (:) :: njsc           !! Index of the dominant soil textural class in the grid cell (1-nscm, unitless)
!$OMP THREADPRIVATE(njsc)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vbeta1         !! Snow resistance 
!$OMP THREADPRIVATE(vbeta1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vbeta4         !! Bare soil resistance
!$OMP THREADPRIVATE(vbeta4)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vbeta5         !! Floodplains resistance
!$OMP THREADPRIVATE(vbeta5)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: soilcap        !!
!$OMP THREADPRIVATE(soilcap)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: soilflx        !!
!$OMP THREADPRIVATE(soilflx)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: temp_sol       !! Soil temperature
!$OMP THREADPRIVATE(temp_sol)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: qsurf          !! near soil air moisture
!$OMP THREADPRIVATE(qsurf)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: flood_res      !! flood reservoir estimate
!$OMP THREADPRIVATE(flood_res)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: flood_frac     !! flooded fraction
!$OMP THREADPRIVATE(flood_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: snow           !! Snow mass [Kg/m^2]
!$OMP THREADPRIVATE(snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: snow_age       !! Snow age @tex ($d$) @endtex
!$OMP THREADPRIVATE(snow_age)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: drysoil_frac   !! Fraction of visibly (albedo) Dry soil (Between 0 and 1)
!$OMP THREADPRIVATE(drysoil_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evap_bare_lim  !! Bare soil stress
!$OMP THREADPRIVATE(evap_bare_lim)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: evap_bare_lim_ns !! Bare soil stress
!$OMP THREADPRIVATE(evap_bare_lim_ns)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: co2_flux       !! CO2 flux (gC/m**2 of average ground/one_day)
!$OMP THREADPRIVATE(co2_flux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evapot         !! Soil Potential Evaporation
!$OMP THREADPRIVATE(evapot)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: evapot_corr    !! Soil Potential Evaporation Correction (Milly 1992)
!$OMP THREADPRIVATE(evapot_corr)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapflo       !! Floodplains evaporation
!$OMP THREADPRIVATE(vevapflo)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapsno       !! Snow evaporation
!$OMP THREADPRIVATE(vevapsno)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vevapnu        !! Bare soil evaporation
!$OMP THREADPRIVATE(vevapnu)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: tot_melt       !! Total melt
!$OMP THREADPRIVATE(tot_melt)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: vbeta          !! Resistance coefficient
!$OMP THREADPRIVATE(vbeta)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: rau            !! Density
!$OMP THREADPRIVATE(rau)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: deadleaf_cover !! Fraction of soil covered by dead leaves
!$OMP THREADPRIVATE(deadleaf_cover)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)    :: ptnlev1        !! 1st level Different levels soil temperature
!$OMP THREADPRIVATE(ptnlev1)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: mc_layh        !! Volumetric soil moisture for each layer in hydrol(liquid + ice) (m3/m3)
!$OMP THREADPRIVATE(mc_layh)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: mcl_layh       !! Volumetric soil moisture for each layer in hydrol(liquid) (m3/m3)
!$OMP THREADPRIVATE(mcl_layh)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:,:)  :: soilmoist      !! Total soil moisture content for each layer in hydrol(liquid + ice) (mm)
!$OMP THREADPRIVATE(soilmoist)

  LOGICAL, SAVE                                    :: l_first_sechiba = .TRUE. !! Flag controlling the intialisation (true/false)
!$OMP THREADPRIVATE(l_first_sechiba)

  ! Variables related to snow processes calculations  

  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: frac_snow_veg   !! Snow cover fraction on vegetation (unitless)
!$OMP THREADPRIVATE(frac_snow_veg)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: frac_snow_nobio !! Snow cover fraction on continental ice, lakes, etc (unitless)
!$OMP THREADPRIVATE(frac_snow_nobio)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowrho      !! snow density for each layer
!$OMP THREADPRIVATE(snowrho)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowheat     !! snow heat content for each layer (J/m2)
!$OMP THREADPRIVATE(snowheat)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowgrain    !! snow grain size (m)
!$OMP THREADPRIVATE(snowgrain)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowtemp     !! snow temperature profile (K)
!$OMP THREADPRIVATE(snowtemp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: snowdz       !! snow layer thickness (m)
!$OMP THREADPRIVATE(snowdz)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: gtemp        !! soil surface temperature
!$OMP THREADPRIVATE(gtemp)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: pgflux       !! net energy into snow pack
!$OMP THREADPRIVATE(pgflux)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: cgrnd_snow   !! Integration coefficient for snow numerical scheme
!$OMP THREADPRIVATE(cgrnd_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:,:)  :: dgrnd_snow   !! Integration coefficient for snow numerical scheme
!$OMP THREADPRIVATE(dgrnd_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: lambda_snow  !! Coefficient of the linear extrapolation of surface temperature 
                                                                  !! from the first and second snow layers
!$OMP THREADPRIVATE(lambda_snow)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION(:)    :: temp_sol_add !! Additional energy to melt snow for snow ablation case (K)
!$OMP THREADPRIVATE(temp_sol_add)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: root_deficit !! water deficit to reach IRRIGATION SOIL MOISTURE TARGET
!$OMP THREADPRIVATE(root_deficit)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: irrig_frac   !! Irrig. fraction interpolated in routing, and saved to pass to slowproc if irrigated_soiltile = .TRUE.
!$OMP THREADPRIVATE(irrig_frac)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: irrigated_next         !! Dynamic irrig. area, calculated in slowproc and passed to routing
                                                                            !!  @tex $(m^{-2})$ @endtex
!$OMP THREADPRIVATE(irrigated_next)
  REAL(r_std), ALLOCATABLE, SAVE, DIMENSION (:)   :: fraction_aeirrig_sw    !! Fraction of area equipped for irrigation from surface water, of irrig_frac
                                                                            !! 1.0 here corresponds to irrig_frac, not grid cell
!$OMP THREADPRIVATE(fraction_aeirrig_sw)

CONTAINS


!!  =============================================================================================================================
!! SUBROUTINE:    sechiba_xios_initialize
!!
!>\BRIEF	  Initialize xios dependant defintion before closing context defintion
!!
!! DESCRIPTION:	  Initialize xios dependant defintion before closing context defintion
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE sechiba_xios_initialize
    LOGICAL :: lerr
    
    IF (xios_orchidee_ok) THEN  
      lerr=xios_orchidee_setvar('min_sechiba',min_sechiba)
      CALL slowproc_xios_initialize
      CALL condveg_xios_initialize
      CALL chemistry_xios_initialize
      CALL thermosoil_xios_initialize
      CALL routing_wrapper_xios_initialize
    END IF
    IF (printlev_loc>=3) WRITE(numout,*) 'End sechiba_xios_initialize'

  END SUBROUTINE sechiba_xios_initialize




!!  =============================================================================================================================
!! SUBROUTINE:    sechiba_initialize
!!
!>\BRIEF	  Initialize all prinicipal modules by calling their "_initialize" subroutines
!!
!! DESCRIPTION:	  Initialize all prinicipal modules by calling their "_initialize" subroutines
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE sechiba_initialize( &
       kjit,         kjpij,        kjpindex,     index,                   &
       lalo,         contfrac,     neighbours,   resolution, zlev,        &
       u,            v,            qair,         temp_air,    &
       petAcoef,     peqAcoef,     petBcoef,     peqBcoef,                &
       precip_rain,  precip_snow,  lwdown,       swnet,      swdown,      &
       pb,           rest_id,      hist_id,      hist2_id,                &
       rest_id_stom, hist_id_stom, hist_id_stom_IPCC,                     &
       coastalflow,  riverflow,    tsol_rad,     vevapp,       qsurf_out, &
       z0m_out,      z0h_out,      albedo,       fluxsens,     fluxlat,      emis_out,  &
       temp_sol_new, tq_cdrag)

!! 0.1 Input variables
    INTEGER(i_std), INTENT(in)                               :: kjit              !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpij             !! Total size of the un-compressed grid 
                                                                                  !! (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size - terrestrial pixels only 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id           !! _Restart_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id           !! _History_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist2_id          !! _History_ file 2 identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id_stom      !! STOMATE's _Restart_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id_stom      !! STOMATE's _History_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT(in)                                :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file file 
                                                                                  !! identifier (unitless)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)          :: lalo              !! Geographic coordinates (latitude,longitude)
                                                                                  !! for grid cells (degrees)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: contfrac          !! Fraction of continent in the grid 
                                                                                  !! (unitless, 0-1)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: index             !! Indices of the pixels on the map. 
                                                                                  !! Sechiba uses a reduced grid excluding oceans
                                                                                  !! ::index contains the indices of the 
                                                                                  !! terrestrial pixels only! (unitless)
    INTEGER(i_std), DIMENSION (kjpindex,NbNeighb), INTENT(in):: neighbours        !! Neighboring grid points if land!(unitless)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)          :: resolution        !! Size in x and y of the grid (m)
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u                 !! Lowest level wind speed in direction u 
                                                                                  !! @tex $(m.s^{-1})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v                 !! Lowest level wind speed in direction v 
                                                                                  !! @tex $(m.s^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: zlev              !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair              !! Lowest level specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_rain       !! Rain precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_snow       !! Snow precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown            !! Down-welling long-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet             !! Net surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown            !! Down-welling surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_air          !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petAcoef          !! Coefficients A for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqAcoef          !! Coefficients A for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petBcoef          !! Coefficients B for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqBcoef          !! Coefficients B for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb                !! Surface pressure (hPa)


!! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: coastalflow       !! Outflow on coastal points by small basins.
                                                                                  !! This is the water which flows in a disperse 
                                                                                  !! way into the ocean
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: riverflow         !! Outflow of the major rivers.
                                                                                  !! The flux will be located on the continental 
                                                                                  !! grid but this should be a coastal point  
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: tsol_rad          !! Radiative surface temperature 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vevapp            !! Total of evaporation 
                                                                                  !! @tex $(kg m^{-2} days^{-1})$ @endtex
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsurf_out         !! Surface specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0m_out           !! Surface roughness momentum (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0h_out           !! Surface roughness heat (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out)         :: albedo            !! Surface albedo for visible and near-infrared
                                                                                  !! (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsens          !! Sensible heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxlat           !! Latent heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: emis_out          !! Emissivity (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: temp_sol_new      !! New ground temperature (K)

!! 0.3 Modified
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: tq_cdrag          !! Surface drag coefficient (-)

!! 0.4 Local variables
    INTEGER(i_std)                                           :: ji, jv		  !! Index (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: zmaxh_glo         !! 2D field of constant soil depth (zmaxh) (m)
    CHARACTER(LEN=80)                                        :: var_name          !! To store variables names for I/O (unitless)

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Start sechiba_initialize'

    !! 1. Initialize variables on first call

    !! 1.1 Initialize most of sechiba's variables
    CALL sechiba_init (kjit, kjpij, kjpindex, index, rest_id, lalo)
    
    !! 1.3 Initialize stomate's variables

    CALL slowproc_initialize (kjit,          kjpij,        kjpindex,                          &
                              rest_id,       rest_id_stom, hist_id_stom,   hist_id_stom_IPCC, &
                              index,         indexveg,     lalo,           neighbours,        &
                              resolution,    contfrac,     temp_air,                          &
                              soiltile,      reinf_slope,  ks,             nvan,              &
                              avan,          mcr,          mcs,            mcfc,              &
                              mcw,           deadleaf_cover,               assim_param,       &
                              lai,           frac_age,     height,         veget,             &
                              frac_nobio,    njsc,         veget_max,      fraclut,           &
                              nwdfraclut,    tot_bare_soil,totfrac_nobio,  qsintmax,          &
                              temp_growth,   irrigated_next, irrig_frac,   fraction_aeirrig_sw, &
                              reinf_slope_soil)
    
    !! 1.4 Initialize diffusion coefficients
    CALL diffuco_initialize (kjit,    kjpindex, index,                  &
                             rest_id, lalo,     neighbours, resolution, &
                             rstruct, tq_cdrag)
    
    !! 1.5 Initialize remaining variables of energy budget
    CALL enerbil_initialize (kjit,     kjpindex,     index,    rest_id,  &
                             qair,                                       &
                             temp_sol, temp_sol_new, tsol_rad,           &
                             evapot,   evapot_corr,  qsurf,    fluxsens, &
                             fluxlat,  vevapp )
    
    
    !! 1.7 Initialize remaining hydrological variables
    CALL hydrol_initialize (ks,  nvan, avan, mcr, mcs, mcfc, mcw,    &
         kjit,           kjpindex,  index,         rest_id,          &
         njsc,           soiltile,  veget,         veget_max,        &
         humrel,         vegstress, drysoil_frac,                    &
         shumdiag_perma,    qsintveg,                        &
         evap_bare_lim, evap_bare_lim_ns,  snow,      snow_age,      snow_nobio,       &
         snow_nobio_age, snowrho,   snowtemp,      snowgrain,        &
         snowdz,         snowheat,  &
         mc_layh,    mcl_layh,  soilmoist)

    
    !! 1.9 Initialize surface parameters (emissivity, albedo and roughness)
    CALL condveg_initialize (kjit, kjpindex, index, rest_id, &
         lalo, neighbours, resolution, contfrac, veget, veget_max, frac_nobio, totfrac_nobio, &
         zlev, snow, snow_age, snow_nobio, snow_nobio_age, &
         drysoil_frac, height, snowdz,snowrho, tot_bare_soil, &
         temp_air, pb, u, v, lai, &
         emis, albedo, z0m, z0h, roughheight, &
         frac_snow_veg,frac_snow_nobio)
    

    !! 1.10 Initialization of soil thermodynamics
    CALL thermosoil_initialize (kjit, kjpindex, rest_id, mcs,  &
         temp_sol_new, snow,       shumdiag_perma,        &
         soilcap,      soilflx,    stempdiag, ftempdiag,             &
         gtemp,               &
         mc_layh,  mcl_layh,   soilmoist,       njsc ,     &
         frac_snow_veg,frac_snow_nobio,totfrac_nobio,     &
         snowdz, snowrho, snowtemp, lambda_snow, cgrnd_snow, dgrnd_snow, pb)


    !! 1.12 Initialize river routing
    IF ( river_routing .AND. nbp_glo .GT. 1) THEN

       !! 1.12.1 Initialize river routing
       CALL routing_wrapper_initialize( &
            kjit,        kjpindex,       index,                 &
            rest_id,     hist_id,        hist2_id,   lalo,      &
            neighbours,  resolution,     contfrac,   stempdiag, ftempdiag, &
            soiltile,    irrig_frac,     veget_max,  irrigated_next, &
            returnflow,  reinfiltration, irrigation, riverflow, &
            coastalflow, flood_frac,     flood_res )
    ELSE
       !! 1.12.2 No routing, set variables to zero
       riverflow(:) = zero
       coastalflow(:) = zero
       returnflow(:) = zero
       reinfiltration(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero
    ENDIF
    
    !! 1.13 Write internal variables to output fields
    z0m_out(:) = z0m(:)
    z0h_out(:) = z0h(:)
    emis_out(:) = emis(:) 
    qsurf_out(:) = qsurf(:)

    !! 2. Output variables only once
    zmaxh_glo(:) = zmaxh
    CALL xios_orchidee_send_field("zmaxh",zmaxh_glo)

    IF (printlev_loc>=3) WRITE(numout,*) 'sechiba_initialize done'

  END SUBROUTINE sechiba_initialize

!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_main
!!
!>\BRIEF        Main routine for the sechiba module performing three functions:
!! calculating temporal evolution of all variables and preparation of output and 
!! restart files (during the last call only)
!!
!!\n DESCRIPTION : Main routine for the sechiba module. 
!! One time step evolution consists of:
!! - call sechiba_var_init to do some initialization,
!! - call slowproc_main to do some daily calculations
!! - call diffuco_main for diffusion coefficient calculation,
!! - call enerbil_main for energy budget calculation,
!! - call hydrol_main for hydrologic processes calculation,
!! - call condveg_main for surface conditions such as roughness, albedo, and emmisivity,
!! - call thermosoil_main for soil thermodynamic calculation,
!! - call sechiba_end to swap previous to new fields.
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): Hydrological variables (:: coastalflow and :: riverflow),
!! components of the energy budget (:: tsol_rad, :: vevapp, :: fluxsens, 
!! :: temp_sol_new and :: fluxlat), surface characteristics (:: z0_out, :: emis_out, 
!! :: tq_cdrag and :: albedo) and land use related CO2 fluxes (:: netco2flux and 
!! :: fco2_lu, :: fco2_wh, ::fco2_ha)            
!!
!! REFERENCE(S)	: 
!!
!! FLOWCHART    : 
!! \latexonly 
!! \includegraphics[scale = 0.5]{sechibamainflow.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  SUBROUTINE sechiba_main (kjit, kjpij, kjpindex, index, &
       & ldrestart_read, ldrestart_write, &
       & lalo, contfrac, neighbours, resolution,&
       & zlev, u, v, qair, temp_air, epot_air, ccanopy, &
       & tq_cdrag, petAcoef, peqAcoef, petBcoef, peqBcoef, &
       & precip_rain, precip_snow, lwdown, swnet, swdown, coszang, pb, &
       & vevapp, fluxsens, fluxlat, coastalflow, riverflow, &
       & netco2flux, fco2_lu, fco2_wh, fco2_ha, &
       & tsol_rad, temp_sol_new, qsurf_out, albedo, emis_out, z0m_out, z0h_out,&
       & veget_out, lai_out, height_out, &
       & rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC)

!! 0.1 Input variables
    
    INTEGER(i_std), INTENT(in)                               :: kjit              !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpij             !! Total size of the un-compressed grid 
                                                                                  !! (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size - terrestrial pixels only 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id           !! _Restart_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id           !! _History_ file identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist2_id          !! _History_ file 2 identifier (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id_stom      !! STOMATE's _Restart_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: hist_id_stom      !! STOMATE's _History_ file identifier 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT(in)                                :: hist_id_stom_IPCC !! STOMATE's IPCC _history_ file file 
                                                                                  !! identifier (unitless)
    LOGICAL, INTENT(in)                                      :: ldrestart_read    !! Logical for _restart_ file to read 
                                                                                  !! (true/false)
    LOGICAL, INTENT(in)                                      :: ldrestart_write   !! Logical for _restart_ file to write 
                                                                                  !! (true/false)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)          :: lalo              !! Geographic coordinates (latitude,longitude)
                                                                                  !! for grid cells (degrees)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: contfrac          !! Fraction of continent in the grid 
                                                                                  !! (unitless, 0-1)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: index             !! Indices of the pixels on the map. 
                                                                                  !! Sechiba uses a reduced grid excluding oceans
                                                                                  !! ::index contains the indices of the 
                                                                                  !! terrestrial pixels only! (unitless)
    INTEGER(i_std), DIMENSION(kjpindex,NbNeighb), INTENT(in) :: neighbours        !! Neighboring grid points if land!(unitless)
    REAL(r_std), DIMENSION (kjpindex,2), INTENT(in)          :: resolution        !! Size in x and y of the grid (m)
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: u                 !! Lowest level wind speed in direction u 
                                                                                  !! @tex $(m.s^{-1})$ @endtex 
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: v                 !! Lowest level wind speed in direction v 
                                                                                  !! @tex $(m.s^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: zlev              !! Height of first layer (m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: qair              !! Lowest level specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_rain       !! Rain precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: precip_snow       !! Snow precipitation 
                                                                                  !! @tex $(kg m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: lwdown            !! Down-welling long-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: coszang           !! Cosine of the solar zenith angle (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swnet             !! Net surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: swdown            !! Down-welling surface short-wave flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: temp_air          !! Air temperature (K)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: epot_air          !! Air potential energy (??J)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: ccanopy           !! CO2 concentration in the canopy (ppm)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petAcoef          !! Coefficients A for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqAcoef          !! Coefficients A for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: petBcoef          !! Coefficients B for T from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: peqBcoef          !! Coefficients B for q from the Planetary 
                                                                                  !! Boundary Layer
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)            :: pb                !! Surface pressure (hPa)


!! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: coastalflow       !! Outflow on coastal points by small basins.
                                                                                  !! This is the water which flows in a disperse 
                                                                                  !! way into the ocean
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: riverflow         !! Outflow of the major rivers.
                                                                                  !! The flux will be located on the continental 
                                                                                  !! grid but this should be a coastal point  
                                                                                  !! @tex $(kg dt_routing^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: tsol_rad          !! Radiative surface temperature 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: vevapp            !! Total of evaporation 
                                                                                  !! @tex $(kg m^{-2} days^{-1})$ @endtex
    
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: qsurf_out         !! Surface specific humidity 
                                                                                  !! @tex $(kg kg^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0m_out           !! Surface roughness momentum (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: z0h_out           !! Surface roughness heat (output diagnostic, m)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (out)         :: albedo            !! Surface albedo for visible and near-infrared
                                                                                  !! (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxsens          !! Sensible heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fluxlat           !! Latent heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: emis_out          !! Emissivity (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: netco2flux        !! Sum CO2 flux over PFTs 
                                                                                  !! (gC/m2/dt_sechiba)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fco2_lu           !! Land Cover Change CO2 flux (gC/m2/one_day)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fco2_wh           !! Wood harvest CO2 flux (gC/m2/one_day)
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)           :: fco2_ha           !! Crop harvest CO2 flux (gC/m2/one_day)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: veget_out         !! Fraction of vegetation type (unitless, 0-1)
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: lai_out           !! Leaf area index (m^2 m^{-2}) 
    REAL(r_std),DIMENSION (kjpindex,nvm), INTENT (out)       :: height_out        !! Vegetation Height (m)

!! 0.3 Modified

    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: tq_cdrag          !! Surface drag coefficient  (-)
    REAL(r_std),DIMENSION (kjpindex), INTENT (inout)         :: temp_sol_new      !! New ground temperature (K)

!! 0.4 local variables

    INTEGER(i_std)                                           :: ji, jv		  !! Index (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: histvar           !! Computations for history files (unitless)
    REAL(r_std), DIMENSION(kjpindex,nlut)                    :: histvar2          !! Computations for history files (unitless)
    CHARACTER(LEN=80)                                        :: var_name          !! To store variables names for I/O (unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treefrac      !! Total fraction occupied by trees (0-1, uniless) 
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_grassfracC3   !! Total fraction occupied by C3 grasses (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_grassfracC4   !! Total fraction occupied by C4 grasses (0-1, unitless)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_cropfracC3    !! Total fraction occupied by C3 crops (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_cropfracC4    !! Total fraction occupied by C4 crops (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treeFracNdlEvg!! Total fraction occupied by treeFracNdlEvg (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treeFracBdlEvg!! Total fraction occupied by treeFracBdlEvg (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treeFracNdlDcd!! Total fraction occupied by treeFracNdlDcd (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: sum_treeFracBdlDcd!! Total fraction occupied by treeFracBdlDcd (0-1, unitess)
    REAL(r_std), DIMENSION(kjpindex)                         :: grndflux          !! Net energy into soil (W/m2)
    REAL(r_std), DIMENSION(kjpindex,nsnow)                   :: snowliq           !! Liquid water content (m)
    REAL(r_std), DIMENSION(kjpindex)                         :: snow_age_diag     !! Only for diag, contains xios_default_val
    REAL(r_std), DIMENSION(kjpindex,nnobio)                  :: snow_nobio_age_diag !! Only for diag, contains xios_default_val
    REAL(r_std), DIMENSION(kjpindex)                         :: snowage_glob      !! Snow age on total area including snow on vegetated and bare soil and nobio area @tex ($d$) @endtex 
    REAL(r_std), DIMENSION(kjpindex,nlut)                    :: gpplut            !! GPP on landuse tile, only for diagnostics


!_ ================================================================================================================================

    IF (printlev_loc>=3) WRITE(numout,*) 'Start sechiba_main kjpindex =',kjpindex

    !! 1. Initialize variables at each time step
    CALL sechiba_var_init (kjpindex, rau, pb, temp_air) 

    !! 2. Compute diffusion coefficients
    CALL diffuco_main (kjit, kjpindex, index, indexveg, indexlai, u, v, &
         & zlev, z0m, z0h, roughheight, temp_sol, temp_air, temp_growth, rau, tq_cdrag, qsurf, qair, pb ,  &
         & evap_bare_lim, evap_bare_lim_ns, evapot, evapot_corr, snow, flood_frac, flood_res, &
         & frac_nobio, snow_nobio, totfrac_nobio, &
         & swnet, swdown, coszang, ccanopy, humrel, veget, veget_max, lai, qsintveg, qsintmax, assim_param, &
         & vbeta, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, gsmean, rveget, rstruct, cimean, gpp, &
         & lalo, neighbours, resolution, ptnlev1, precip_rain, frac_age, tot_bare_soil, frac_snow_veg, frac_snow_nobio, &
         & hist_id, hist2_id)
   
    !! 3. Compute energy balance
    CALL enerbil_main (kjit, kjpindex, &
         & index, indexveg, zlev, lwdown, swnet, epot_air, temp_air, u, v, petAcoef, petBcoef, &
         & qair, peqAcoef, peqBcoef, pb, rau, vbeta, vbeta1, vbeta2, vbeta3, vbeta3pot, vbeta4, vbeta5, &
         & emis, soilflx, soilcap, tq_cdrag, humrel, fluxsens, fluxlat, &
         & vevapp, transpir, transpot, vevapnu, vevapwet, vevapsno, vevapflo, temp_sol, tsol_rad, &
         & temp_sol_new, qsurf, evapot, evapot_corr, rest_id, hist_id, hist2_id, &
         & precip_rain,  pgflux, snowdz, temp_sol_add)

 
    !! 4. Compute hydrology
    !! 4.1 Water balance from CWRR module (11 soil layers)
    CALL hydrol_main (ks,  nvan, avan, mcr, mcs, mcfc, mcw, kjit, kjpindex, &
         & index, indexveg, indexsoil, indexlayer, indexnslm, &
         & temp_sol_new, floodout, runoff, drainage, frac_nobio, totfrac_nobio, vevapwet, veget, veget_max, njsc, &
         & qsintmax, qsintveg, vevapnu, vevapsno, vevapflo, snow, snow_age, snow_nobio, snow_nobio_age,  &
         & tot_melt, transpir, precip_rain, precip_snow, returnflow, reinfiltration, irrigation, &
         & humrel, vegstress, drysoil_frac, evapot, evapot_corr, evap_bare_lim, evap_bare_lim_ns, flood_frac, flood_res, &
         & shumdiag,shumdiag_perma, k_litt, litterhumdiag, soilcap, soiltile, fraclut, reinf_slope_soil,&
         & rest_id, hist_id, hist2_id,&
         & contfrac, stempdiag, &
         & temp_air, pb, u, v, tq_cdrag, swnet, pgflux, &
         & snowrho,snowtemp,snowgrain,snowdz,snowheat,snowliq, &
         & grndflux,gtemp,tot_bare_soil, &
         & lambda_snow,cgrnd_snow,dgrnd_snow,frac_snow_veg,temp_sol_add, &
         & mc_layh, mcl_layh, soilmoist, root_deficit)

         
    !! 6. Compute surface variables (emissivity, albedo and roughness)
    CALL condveg_main (kjit, kjpindex, index, rest_id, hist_id, hist2_id, &
         lalo, neighbours, resolution, contfrac, veget, veget_max, frac_nobio, totfrac_nobio, &
         zlev, snow, snow_age, snow_nobio, snow_nobio_age, &
         drysoil_frac, height, snowdz, snowrho, tot_bare_soil, &
         temp_air, pb, u, v, lai, &
         emis, albedo, z0m, z0h, roughheight, &
         frac_snow_veg, frac_snow_nobio)


    !! 7. Compute soil thermodynamics
    CALL thermosoil_main (kjit, kjpindex, &
         index, indexgrnd, mcs, &
         temp_sol_new, snow, soilcap, soilflx, &
         shumdiag_perma, stempdiag, ftempdiag, ptnlev1, rest_id, hist_id, hist2_id, &
         snowdz,snowrho,snowtemp,gtemp,pb,&
         mc_layh, mcl_layh, soilmoist, njsc,frac_snow_veg,frac_snow_nobio,totfrac_nobio,temp_sol_add, &
         lambda_snow, cgrnd_snow, dgrnd_snow)


    !! 8. Compute river routing 
    IF ( river_routing .AND. nbp_glo .GT. 1) THEN
       !! 8.1 River routing
       CALL routing_wrapper_main (kjit, kjpindex, index, &
            & lalo, neighbours, resolution, contfrac, totfrac_nobio, veget_max, floodout, runoff, &
            & drainage, transpot, precip_rain, humrel, k_litt, flood_frac, flood_res, stempdiag, &
            & ftempdiag, reinf_slope, returnflow, reinfiltration, irrigation, riverflow, coastalflow, rest_id, hist_id, hist2_id, &
            & soiltile, root_deficit, irrigated_next, irrig_frac, fraction_aeirrig_sw)
    ELSE
       !! 8.2 No routing, set variables to zero
       riverflow(:) = zero
       coastalflow(:) = zero
       returnflow(:) = zero
       reinfiltration(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero

       CALL xios_orchidee_send_field("coastalflow",coastalflow/dt_sechiba)
       CALL xios_orchidee_send_field("riverflow",riverflow/dt_sechiba)
    ENDIF

    !! 9. Compute slow processes (i.e. 'daily' and annual time step)
    CALL slowproc_main (kjit, kjpij, kjpindex, &
         index, indexveg, lalo, neighbours, resolution, contfrac, soiltile, fraclut, nwdFraclut, &
         temp_air, temp_sol, stempdiag, &
         vegstress, shumdiag, litterhumdiag, precip_rain, precip_snow, gpp, &
         deadleaf_cover, &
         assim_param, &
         lai, frac_age, height, veget, frac_nobio, veget_max, totfrac_nobio, qsintmax, &
         rest_id, hist_id, hist2_id, rest_id_stom, hist_id_stom, hist_id_stom_IPCC, &
         co2_flux, fco2_lu, fco2_wh, fco2_ha, temp_growth, tot_bare_soil, &
         irrigated_next, irrig_frac, reinf_slope, reinf_slope_soil)


    !! 9.2 Compute global CO2 flux
    netco2flux(:) = zero
    DO jv = 2,nvm
      netco2flux(:) = netco2flux(:) + co2_flux(:,jv)*(1-totfrac_nobio)
    ENDDO
 
    !! 10. Update the temperature (temp_sol) with newly computed values
    CALL sechiba_end (kjpindex, temp_sol_new, temp_sol)

   
    !! 11. Write internal variables to output fields
    z0m_out(:) = z0m(:)
    z0h_out(:) = z0h(:)
    emis_out(:) = emis(:)
    qsurf_out(:) = qsurf(:)
    veget_out(:,:)  = veget(:,:)
    lai_out(:,:)    = lai(:,:)
    height_out(:,:) = height(:,:)
 
    !! 12. Write global variables to history files
    sum_treefrac(:) = zero
    sum_grassfracC3(:) = zero
    sum_grassfracC4(:) = zero
    sum_cropfracC3(:) = zero
    sum_cropfracC4(:) = zero
    sum_treeFracNdlEvg(:) = zero
    sum_treeFracBdlEvg(:) = zero
    sum_treeFracNdlDcd(:) = zero
    sum_treeFracBdlDcd(:) = zero

    DO jv = 2, nvm 
       IF (is_tree(jv) .AND. natural(jv)) THEN
          sum_treefrac(:) = sum_treefrac(:) + veget_max(:,jv)
       ELSE IF ((.NOT. is_tree(jv))  .AND. natural(jv)) THEN
          ! Grass
          IF (is_c4(jv)) THEN
             sum_grassfracC4(:) = sum_grassfracC4(:) + veget_max(:,jv)
          ELSE
             sum_grassfracC3(:) = sum_grassfracC3(:) + veget_max(:,jv)
          END IF
       ELSE 
          ! Crop and trees not natural
          IF (is_c4(jv)) THEN
             sum_cropfracC4(:) = sum_cropfracC4(:) + veget_max(:,jv)
          ELSE
             sum_cropfracC3(:) = sum_cropfracC3(:) + veget_max(:,jv)
          END IF
       ENDIF

       IF (is_tree(jv)) THEN
          IF (is_evergreen(jv)) THEN
             IF (is_needleleaf(jv)) THEN
                ! Fraction for needleleaf evergreen trees (treeFracNdlEvg)
                sum_treeFracNdlEvg(:) = sum_treeFracNdlEvg(:) + veget_max(:,jv)
             ELSE
                ! Fraction for broadleaf evergreen trees (treeFracBdlEvg)
                sum_treeFracBdlEvg(:) = sum_treeFracBdlEvg(:) + veget_max(:,jv)
             END IF
          ELSE IF (is_deciduous(jv)) THEN
             IF (is_needleleaf(jv)) THEN
                ! Fraction for needleleaf deciduous trees (treeFracNdlDcd)
                sum_treeFracNdlDcd(:) = sum_treeFracNdlDcd(:) + veget_max(:,jv)
             ELSE 
                ! Fraction for broadleafs deciduous trees (treeFracBdlDcd)
                sum_treeFracBdlDcd(:) = sum_treeFracBdlDcd(:) + veget_max(:,jv)
             END IF
          END IF
       END IF
    ENDDO          

    histvar(:)=zero
    DO jv = 2, nvm
       IF (is_deciduous(jv)) THEN
          histvar(:) = histvar(:) + veget_max(:,jv)*100*contfrac
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("treeFracPrimDec",histvar)

    histvar(:)=zero
    DO jv = 2, nvm
       IF (is_evergreen(jv)) THEN
          histvar(:) = histvar(:) + veget_max(:,jv)*100*contfrac
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("treeFracPrimEver",histvar)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( .NOT.(is_c4(jv)) ) THEN
          histvar(:) = histvar(:) + veget_max(:,jv)*100*contfrac
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("c3PftFrac",histvar)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( is_c4(jv) ) THEN
          histvar(:) = histvar(:) + veget_max(:,jv)*100*contfrac
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("c4PftFrac",histvar)

    CALL xios_orchidee_send_field("temp_sol_new",temp_sol_new)
    CALL xios_orchidee_send_field("fluxsens",fluxsens)
    CALL xios_orchidee_send_field("fluxlat",fluxlat)


    ! Add XIOS default value where no snow
    DO ji=1,kjpindex 
       IF (snow(ji) .GT. zero) THEN
          snow_age_diag(ji) = snow_age(ji)
          snow_nobio_age_diag(ji,:) = snow_nobio_age(ji,:)
       
          snowage_glob(ji) = snow_age(ji)*frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + &
               SUM(snow_nobio_age(ji,:)*frac_snow_nobio(ji,:)*frac_nobio(ji,:))
          IF (snowage_glob(ji) .NE. 0) snowage_glob(ji) = snowage_glob(ji) / &
               (frac_snow_veg(ji)*(1-totfrac_nobio(ji)) + SUM(frac_snow_nobio(ji,:)*frac_nobio(ji,:)))
       ELSE
          snow_age_diag(ji) = xios_default_val
          snow_nobio_age_diag(ji,:) = xios_default_val
          snowage_glob(ji) = xios_default_val
       END IF
    END DO
    
    CALL xios_orchidee_send_field("snow",snow)
    CALL xios_orchidee_send_field("snowage",snow_age_diag)
    CALL xios_orchidee_send_field("snownobio",snow_nobio)
    CALL xios_orchidee_send_field("snownobioage",snow_nobio_age_diag)
    CALL xios_orchidee_send_field("snowage_glob",snowage_glob)

    CALL xios_orchidee_send_field("frac_snow", SUM(frac_snow_nobio,2)*totfrac_nobio+frac_snow_veg*(1-totfrac_nobio))
    CALL xios_orchidee_send_field("frac_snow_veg", frac_snow_veg)
    CALL xios_orchidee_send_field("frac_snow_nobio", frac_snow_nobio)
    CALL xios_orchidee_send_field("reinf_slope",reinf_slope)
    CALL xios_orchidee_send_field("njsc",REAL(njsc, r_std))
    CALL xios_orchidee_send_field("vegetfrac",veget)
    CALL xios_orchidee_send_field("maxvegetfrac",veget_max)
    CALL xios_orchidee_send_field("irrigmap_dyn",irrigated_next)
    CALL xios_orchidee_send_field("aei_sw",fraction_aeirrig_sw)
    CALL xios_orchidee_send_field("nobiofrac",frac_nobio)
    CALL xios_orchidee_send_field("soiltile",soiltile)
    CALL xios_orchidee_send_field("rstruct",rstruct)
    CALL xios_orchidee_send_field("gpp",gpp/dt_sechiba)
    CALL xios_orchidee_send_field("gpp_ipcc2",SUM(gpp,dim=2)/dt_sechiba)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( .NOT. is_tree(jv) .AND. natural(jv) ) THEN
          histvar(:) = histvar(:) + gpp(:,jv)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("gppgrass",histvar/dt_sechiba)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( (.NOT. is_tree(jv)) .AND. (.NOT. natural(jv)) ) THEN
          histvar(:) = histvar(:) + gpp(:,jv)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("gppcrop",histvar/dt_sechiba)

    histvar(:)=zero
    DO jv = 2, nvm
       IF ( is_tree(jv) ) THEN
          histvar(:) = histvar(:) + gpp(:,jv)
       ENDIF
    ENDDO
    CALL xios_orchidee_send_field("gpptree",histvar/dt_sechiba)
    CALL xios_orchidee_send_field("nee",co2_flux/1.e3/one_day)
    CALL xios_orchidee_send_field("drysoil_frac",drysoil_frac)
    CALL xios_orchidee_send_field("vevapflo",vevapflo/dt_sechiba)
    CALL xios_orchidee_send_field("k_litt",k_litt)
    CALL xios_orchidee_send_field("beta",vbeta)
    CALL xios_orchidee_send_field("vbeta1",vbeta1)
    CALL xios_orchidee_send_field("vbeta2",vbeta2)
    CALL xios_orchidee_send_field("vbeta3",vbeta3)
    CALL xios_orchidee_send_field("vbeta4",vbeta4)
    CALL xios_orchidee_send_field("vbeta5",vbeta5)
    CALL xios_orchidee_send_field("gsmean",gsmean)
    CALL xios_orchidee_send_field("cimean",cimean)
    CALL xios_orchidee_send_field("rveget",rveget)
 
    histvar(:)=SUM(vevapwet(:,:),dim=2)
    CALL xios_orchidee_send_field("evspsblveg",histvar/dt_sechiba)
    histvar(:)= vevapnu(:)+vevapsno(:)
    CALL xios_orchidee_send_field("evspsblsoi",histvar/dt_sechiba)
    histvar(:)=SUM(transpir(:,:),dim=2)
    CALL xios_orchidee_send_field("tran",histvar/dt_sechiba)

    ! For CMIP6 data request: the following fractions are fractions of the total grid-cell,
    ! which explains the multiplication by contfrac
    CALL xios_orchidee_send_field("treeFrac",sum_treefrac(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("grassFracC3",sum_grassfracC3(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("grassFracC4",sum_grassfracC4(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("cropFracC3",sum_cropfracC3(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("cropFracC4",sum_cropfracC4(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("treeFracNdlEvg",sum_treeFracNdlEvg(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("treeFracBdlEvg",sum_treeFracBdlEvg(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("treeFracNdlDcd",sum_treeFracNdlDcd(:)*100*contfrac(:))
    CALL xios_orchidee_send_field("treeFracBdlDcd",sum_treeFracBdlDcd(:)*100*contfrac(:))

    histvar(:)=veget_max(:,1)*100*contfrac(:)
    CALL xios_orchidee_send_field("baresoilFrac",histvar)
    histvar(:)=SUM(frac_nobio(:,1:nnobio),dim=2)*100*contfrac(:)
    CALL xios_orchidee_send_field("residualFrac",histvar)

    ! For CMIP6 data request: cnc = canopy cover fraction over land area
    histvar(:)=zero
    DO jv=2,nvm
       histvar(:) = histvar(:) + veget_max(:,jv)*100
    END DO
    CALL xios_orchidee_send_field("cnc",histvar)
    
    CALL xios_orchidee_send_field("tsol_rad",tsol_rad-273.15)
    CALL xios_orchidee_send_field("qsurf",qsurf)
    CALL xios_orchidee_send_field("emis",emis)
    CALL xios_orchidee_send_field("z0m",z0m)
    CALL xios_orchidee_send_field("z0h",z0h)
    CALL xios_orchidee_send_field("roughheight",roughheight)
    CALL xios_orchidee_send_field("lai",lai)
    histvar(:)=zero   
    DO ji = 1, kjpindex
       IF (SUM(veget_max(ji,:)) > zero) THEN
         DO jv=2,nvm
            histvar(ji) = histvar(ji) + veget_max(ji,jv)*lai(ji,jv)/SUM(veget_max(ji,:))
         END DO
       END IF
    END DO
    CALL xios_orchidee_send_field("LAImean",histvar)
    CALL xios_orchidee_send_field("vevapsno",vevapsno/dt_sechiba)
    CALL xios_orchidee_send_field("vevapp",vevapp/dt_sechiba)
    CALL xios_orchidee_send_field("vevapnu",vevapnu/dt_sechiba)
    CALL xios_orchidee_send_field("transpir",transpir*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("transpot",transpot*one_day/dt_sechiba)
    CALL xios_orchidee_send_field("inter",vevapwet*one_day/dt_sechiba)
    histvar(:)=zero
    DO jv=1,nvm
      histvar(:) = histvar(:) + vevapwet(:,jv)
    ENDDO
    CALL xios_orchidee_send_field("ECanop",histvar/dt_sechiba)
    histvar(:)=zero
    DO jv=1,nvm
      histvar(:) = histvar(:) + transpir(:,jv)
    ENDDO
    CALL xios_orchidee_send_field("TVeg",histvar/dt_sechiba)


    !! Calculate diagnostic variables on Landuse tiles for LUMIP/CMIP6

    ! Calculate fraction of landuse tiles related to the whole grid cell
    DO jv=1,nlut
       histvar2(:,jv) = fraclut(:,jv) * contfrac(:)
    END DO
    CALL xios_orchidee_send_field("fraclut",histvar2)

    CALL xios_orchidee_send_field("nwdFraclut",nwdFraclut(:,:))
   
    ! Calculate GPP on landuse tiles
    ! val_exp is used as missing value where the values are not known i.e. where the tile is not represented 
    ! or for pasture (id_pst) or urban land (id_urb). 
    gpplut(:,:)=0
    DO jv=1,nvm
       IF (natural(jv)) THEN
          gpplut(:,id_psl) = gpplut(:,id_psl) + gpp(:,jv)
       ELSE
          gpplut(:,id_crp) = gpplut(:,id_crp) + gpp(:,jv)
       ENDIF
    END DO

    ! Transform from gC/m2/s into kgC/m2/s
    WHERE (fraclut(:,id_psl)>min_sechiba)
       gpplut(:,id_psl) = gpplut(:,id_psl)/fraclut(:,id_psl)/1000
    ELSEWHERE
       gpplut(:,id_psl) = xios_default_val
    END WHERE
    WHERE (fraclut(:,id_crp)>min_sechiba)
       gpplut(:,id_crp) = gpplut(:,id_crp)/fraclut(:,id_crp)/1000
    ELSEWHERE
       gpplut(:,id_crp) = xios_default_val
    END WHERE
    gpplut(:,id_pst) = xios_default_val
    gpplut(:,id_urb) = xios_default_val

    CALL xios_orchidee_send_field("gpplut",gpplut)


    IF ( .NOT. almaoutput ) THEN
       ! Write history file in IPSL-format
       CALL histwrite_p(hist_id, 'beta', kjit, vbeta, kjpindex, index)
       CALL histwrite_p(hist_id, 'z0m', kjit, z0m, kjpindex, index)
       CALL histwrite_p(hist_id, 'z0h', kjit, z0h, kjpindex, index)
       CALL histwrite_p(hist_id, 'roughheight', kjit, roughheight, kjpindex, index)
       CALL histwrite_p(hist_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
       CALL histwrite_p(hist_id, 'lai', kjit, lai, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'subli', kjit, vevapsno, kjpindex, index)
       CALL histwrite_p(hist_id, 'evapnu', kjit, vevapnu, kjpindex, index)
       CALL histwrite_p(hist_id, 'transpir', kjit, transpir, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'inter', kjit, vevapwet, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'vbeta1', kjit, vbeta1, kjpindex, index)
       CALL histwrite_p(hist_id, 'vbeta2', kjit, vbeta2, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'vbeta3', kjit, vbeta3, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'vbeta4', kjit, vbeta4, kjpindex, index)    
       CALL histwrite_p(hist_id, 'vbeta5', kjit, vbeta5, kjpindex, index)    
       CALL histwrite_p(hist_id, 'drysoil_frac', kjit, drysoil_frac, kjpindex, index)
       CALL histwrite_p(hist_id, 'rveget', kjit, rveget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'rstruct', kjit, rstruct, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'snow', kjit, snow, kjpindex, index)
       CALL histwrite_p(hist_id, 'snowage', kjit, snow_age, kjpindex, index)
       CALL histwrite_p(hist_id, 'snownobio', kjit, snow_nobio, kjpindex*nnobio, indexnobio)
       CALL histwrite_p(hist_id, 'snownobioage', kjit, snow_nobio_age, kjpindex*nnobio, indexnobio)

       CALL histwrite_p(hist_id, 'grndflux', kjit, grndflux, kjpindex,index)
       CALL histwrite_p(hist_id, 'snowtemp',kjit,snowtemp,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowliq', kjit,snowliq,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowdz', kjit,snowdz,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowrho', kjit,snowrho,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowgrain',kjit,snowgrain,kjpindex*nsnow,indexsnow)
       CALL histwrite_p(hist_id, 'snowheat',kjit,snowheat,kjpindex*nsnow,indexsnow)
       
       CALL histwrite_p(hist_id,'frac_snow_veg',kjit,frac_snow_veg,kjpindex,index)
       CALL histwrite_p(hist_id, 'frac_snow_nobio', kjit,frac_snow_nobio,kjpindex*nnobio, indexnobio)
       CALL histwrite_p(hist_id, 'pgflux',kjit,pgflux,kjpindex,index)
       CALL histwrite_p(hist_id, 'soiltile',  kjit, soiltile, kjpindex*nstm, indexsoil)

       CALL histwrite_p(hist_id, 'soilindex',  kjit, REAL(njsc, r_std), kjpindex, index)
       CALL histwrite_p(hist_id, 'reinf_slope',  kjit, reinf_slope, kjpindex, index)
       CALL histwrite_p(hist_id, 'k_litt', kjit, k_litt, kjpindex, index)
       
       IF ( do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'evapflo', kjit, vevapflo, kjpindex, index)
          CALL histwrite_p(hist_id, 'flood_frac', kjit, flood_frac, kjpindex, index)
       ENDIF
       
       CALL histwrite_p(hist_id, 'gsmean', kjit, gsmean, kjpindex*nvm, indexveg)    
       CALL histwrite_p(hist_id, 'gpp', kjit, gpp, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'cimean', kjit, cimean, kjpindex*nvm, indexveg)    
       
       IF ( ok_stomate ) THEN
          CALL histwrite_p(hist_id, 'nee', kjit, co2_flux/1.e3/one_day, kjpindex*nvm, indexveg)    
       ENDIF

       histvar(:)=SUM(vevapwet(:,:),dim=2)
       CALL histwrite_p(hist_id, 'evspsblveg', kjit, histvar, kjpindex, index)

       histvar(:)= vevapnu(:)+vevapsno(:)
       CALL histwrite_p(hist_id, 'evspsblsoi', kjit, histvar, kjpindex, index)

       histvar(:)=SUM(transpir(:,:),dim=2)
       CALL histwrite_p(hist_id, 'tran', kjit, histvar, kjpindex, index)

       histvar(:)= sum_treefrac(:)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'treeFrac', kjit, histvar, kjpindex, index) 

       histvar(:)= (sum_grassfracC3(:)+sum_grassfracC4(:))*100*contfrac(:)
       CALL histwrite_p(hist_id, 'grassFrac', kjit, histvar, kjpindex, index) 

       histvar(:)= (sum_cropfracC3(:)+sum_cropfracC4(:))*100*contfrac(:)
       CALL histwrite_p(hist_id, 'cropFrac', kjit, histvar, kjpindex, index)

       histvar(:)=veget_max(:,1)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'baresoilFrac', kjit, histvar, kjpindex, index)

       histvar(:)=SUM(frac_nobio(:,1:nnobio),dim=2)*100*contfrac(:)
       CALL histwrite_p(hist_id, 'residualFrac', kjit, histvar, kjpindex, index)
    ELSE
       ! Write history file in ALMA format 
       CALL histwrite_p(hist_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
       CALL histwrite_p(hist_id, 'lai', kjit, lai, kjpindex*nvm, indexveg)
       CALL histwrite_p(hist_id, 'ESoil', kjit, vevapnu, kjpindex, index)
       CALL histwrite_p(hist_id, 'EWater', kjit, vevapflo, kjpindex, index)
       CALL histwrite_p(hist_id, 'SWE', kjit, snow, kjpindex, index)
       histvar(:)=zero
       DO jv=1,nvm
          histvar(:) = histvar(:) + transpir(:,jv)
       ENDDO
       CALL histwrite_p(hist_id, 'TVeg', kjit, histvar, kjpindex, index)
       histvar(:)=zero
       DO jv=1,nvm
          histvar(:) = histvar(:) + vevapwet(:,jv)
       ENDDO
       CALL histwrite_p(hist_id, 'ECanop', kjit, histvar, kjpindex, index)
       CALL histwrite_p(hist_id, 'SnowFrac', kjit, vbeta1, kjpindex, index)
       !
       CALL histwrite_p(hist_id, 'Z0m', kjit, z0m, kjpindex, index)
       CALL histwrite_p(hist_id, 'Z0h', kjit, z0h, kjpindex, index)
       CALL histwrite_p(hist_id, 'EffectHeight', kjit, roughheight, kjpindex, index)
       !
       IF ( do_floodplains ) THEN
          CALL histwrite_p(hist_id, 'Qflood', kjit, vevapflo, kjpindex, index)
          CALL histwrite_p(hist_id, 'FloodFrac', kjit, flood_frac, kjpindex, index)
       ENDIF

       CALL histwrite_p(hist_id, 'gsmean', kjit, gsmean, kjpindex*nvm, indexveg)    
       CALL histwrite_p(hist_id, 'cimean', kjit, cimean, kjpindex*nvm, indexveg)    
       CALL histwrite_p(hist_id, 'GPP', kjit, gpp, kjpindex*nvm, indexveg)
       
       IF ( ok_stomate ) THEN
             CALL histwrite_p(hist_id, 'NEE', kjit, co2_flux, kjpindex*nvm, indexveg)    
       ENDIF
    ENDIF ! almaoutput
    
    !! 13. Write additional output file with higher frequency
    IF ( hist2_id > 0 ) THEN
       IF ( .NOT. almaoutput ) THEN
          ! Write history file in IPSL-format
          CALL histwrite_p(hist2_id, 'tsol_rad', kjit, tsol_rad, kjpindex, index)
          CALL histwrite_p(hist2_id, 'qsurf', kjit, qsurf, kjpindex, index)
          CALL histwrite_p(hist2_id, 'albedo', kjit, albedo, kjpindex*2, indexalb)
          CALL histwrite_p(hist2_id, 'emis', kjit, emis, kjpindex, index)
          CALL histwrite_p(hist2_id, 'beta', kjit, vbeta, kjpindex, index)
          CALL histwrite_p(hist2_id, 'z0m', kjit, z0m, kjpindex, index)
          CALL histwrite_p(hist2_id, 'z0h', kjit, z0h, kjpindex, index)
          CALL histwrite_p(hist2_id, 'roughheight', kjit, roughheight, kjpindex, index)
          CALL histwrite_p(hist2_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'lai', kjit, lai, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'subli', kjit, vevapsno, kjpindex, index)
          IF ( do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'vevapflo', kjit, vevapflo, kjpindex, index)
             CALL histwrite_p(hist2_id, 'flood_frac', kjit, flood_frac, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'vevapnu', kjit, vevapnu, kjpindex, index)
          CALL histwrite_p(hist2_id, 'transpir', kjit, transpir, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'inter', kjit, vevapwet, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'vbeta1', kjit, vbeta1, kjpindex, index)
          CALL histwrite_p(hist2_id, 'vbeta2', kjit, vbeta2, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'vbeta3', kjit, vbeta3, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'vbeta4', kjit, vbeta4, kjpindex, index)    
          CALL histwrite_p(hist2_id, 'vbeta5', kjit, vbeta5, kjpindex, index)    
          CALL histwrite_p(hist2_id, 'drysoil_frac', kjit, drysoil_frac, kjpindex, index)
          CALL histwrite_p(hist2_id, 'rveget', kjit, rveget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'rstruct', kjit, rstruct, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'snow', kjit, snow, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snowage', kjit, snow_age, kjpindex, index)
          CALL histwrite_p(hist2_id, 'snownobio', kjit, snow_nobio, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'snownobioage', kjit, snow_nobio_age, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'soilindex',  kjit, REAL(njsc, r_std), kjpindex, index)
          CALL histwrite_p(hist2_id, 'reinf_slope',  kjit, reinf_slope, kjpindex, index)
          
          CALL histwrite_p(hist2_id, 'gsmean', kjit, gsmean, kjpindex*nvm, indexveg)    
          CALL histwrite_p(hist2_id, 'gpp', kjit, gpp, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'cimean', kjit, cimean, kjpindex*nvm, indexveg)    
          
          IF ( ok_stomate ) THEN 
             CALL histwrite_p(hist2_id, 'nee', kjit, co2_flux/1.e3/one_day, kjpindex*nvm, indexveg)    
          ENDIF
       ELSE
          ! Write history file in ALMA format
          CALL histwrite_p(hist2_id, 'vegetfrac', kjit, veget, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'maxvegetfrac', kjit, veget_max, kjpindex*nvm, indexveg)
          CALL histwrite_p(hist2_id, 'nobiofrac', kjit, frac_nobio, kjpindex*nnobio, indexnobio)
          CALL histwrite_p(hist2_id, 'ESoil', kjit, vevapnu, kjpindex, index)
          IF ( do_floodplains ) THEN
             CALL histwrite_p(hist2_id, 'EWater', kjit, vevapflo, kjpindex, index)
             CALL histwrite_p(hist2_id, 'FloodFrac', kjit, flood_frac, kjpindex, index)
          ENDIF
          CALL histwrite_p(hist2_id, 'SWE', kjit, snow, kjpindex, index)
          histvar(:)=zero
          DO jv=1,nvm
             histvar(:) = histvar(:) + transpir(:,jv)
          ENDDO
          CALL histwrite_p(hist2_id, 'TVeg', kjit, histvar, kjpindex, index)
          histvar(:)=zero
          DO jv=1,nvm
             histvar(:) = histvar(:) + vevapwet(:,jv)
          ENDDO
          CALL histwrite_p(hist2_id, 'ECanop', kjit, histvar, kjpindex, index)
          CALL histwrite_p(hist2_id, 'SnowFrac', kjit, vbeta1, kjpindex, index)
          CALL histwrite_p(hist2_id, 'GPP', kjit, gpp, kjpindex*nvm, indexveg)
          
          IF ( ok_stomate ) THEN
             CALL histwrite_p(hist2_id, 'NEE', kjit, co2_flux, kjpindex*nvm, indexveg)    
          ENDIF
       ENDIF ! almaoutput
    ENDIF ! hist2_id


    !! Change the vegetation fractions if a new map was read in slowproc. This is done 
    !! after lcchange has been done in stomatelpj
    IF (done_stomate_lcchange) THEN
       CALL slowproc_change_frac(kjpindex, lai, &
                                 veget_max, veget, frac_nobio, totfrac_nobio, tot_bare_soil, soiltile, fraclut, nwdFraclut)
       done_stomate_lcchange=.FALSE.
    END IF

    !! 14. If it is the last time step, write restart files
    IF (ldrestart_write) THEN
       CALL sechiba_finalize( &
            kjit,     kjpij,  kjpindex, index,   rest_id, &
            tq_cdrag, vevapp, fluxsens, fluxlat, tsol_rad)
    END IF

  END SUBROUTINE sechiba_main


!!  =============================================================================================================================
!! SUBROUTINE:    sechiba_finalize
!!
!>\BRIEF	  Finalize all modules by calling their "_finalize" subroutines.
!!
!! DESCRIPTION:	  Finalize all modules by calling their "_finalize" subroutines. These subroutines will write variables to 
!!                restart file. 
!!
!! \n
!_ ==============================================================================================================================

  SUBROUTINE sechiba_finalize( &
       kjit,     kjpij,  kjpindex, index,   rest_id, &
       tq_cdrag, vevapp, fluxsens, fluxlat, tsol_rad)

!! 0.1 Input variables    
    INTEGER(i_std), INTENT(in)                               :: kjit              !! Time step number (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpij             !! Total size of the un-compressed grid 
                                                                                  !! (unitless)
    INTEGER(i_std), INTENT(in)                               :: kjpindex          !! Domain size - terrestrial pixels only 
                                                                                  !! (unitless)
    INTEGER(i_std),INTENT (in)                               :: rest_id           !! _Restart_ file identifier (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)         :: index             !! Indices of the pixels on the map. 
                                                                                  !! Sechiba uses a reduced grid excluding oceans
                                                                                  !! ::index contains the indices of the 
                                                                                  !! terrestrial pixels only! (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: tsol_rad           !! Radiative surface temperature 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: vevapp             !! Total of evaporation 
                                                                                  !! @tex $(kg m^{-2} days^{-1})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: fluxsens           !! Sensible heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: fluxlat            !! Latent heat flux 
                                                                                  !! @tex $(W m^{-2})$ @endtex
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)           :: tq_cdrag           !! Surface drag coefficient (-)

!! 0.2 Local variables
    INTEGER(i_std)                                          :: ji, jv		  !! Index (unitless)
    REAL(r_std), DIMENSION(kjpindex)                        :: histvar            !! Computations for history files (unitless)
    CHARACTER(LEN=80)                                       :: var_name           !! To store variables names for I/O (unitless)


    !! Write restart file for the next simulation from SECHIBA and other modules

    IF (printlev_loc>=3) WRITE (numout,*) 'Start sechiba_finalize for writing restart files'

    !! 1. Call diffuco_finalize to write restart files
    CALL diffuco_finalize (kjit, kjpindex, rest_id, rstruct )
    
    !! 2. Call energy budget to write restart files
    CALL enerbil_finalize (kjit,   kjpindex,    rest_id,            &
                           evapot, evapot_corr, temp_sol, tsol_rad, &
                           qsurf,  fluxsens,    fluxlat,  vevapp )
    
    !! 3. Call hydrology to write restart files
    CALL hydrol_finalize( kjit,           kjpindex, rest_id,  vegstress,  &
         qsintveg,       humrel,   snow,     snow_age, snow_nobio, &
         snow_nobio_age, snowrho,  snowtemp, snowdz,     &
         snowheat,       snowgrain,    &
         drysoil_frac,   evap_bare_lim, evap_bare_lim_ns)
    
    !! 4. Call condveg to write surface variables to restart files
    CALL condveg_finalize (kjit, kjpindex, rest_id, z0m, z0h, roughheight)
    
    !! 5. Call soil thermodynamic to write restart files
    CALL thermosoil_finalize (kjit,    kjpindex, rest_id,   gtemp, &
         soilcap, soilflx, lambda_snow, cgrnd_snow, dgrnd_snow)


    !! 6. Add river routing to restart files  
    IF ( river_routing .AND. nbp_glo .GT. 1) THEN
       !! 6.1 Call river routing to write restart files 
       CALL routing_wrapper_finalize( kjit, kjpindex, rest_id, flood_frac, flood_res )
    ELSE
       !! 6.2 No routing, set variables to zero
       reinfiltration(:) = zero
       returnflow(:) = zero
       irrigation(:) = zero
       flood_frac(:) = zero
       flood_res(:) = zero
    ENDIF
    
    !! 7. Call slowproc_main to add 'daily' and annual variables to restart file
    CALL slowproc_finalize (kjit,       kjpindex,  rest_id,  index,      &
                            njsc,       lai,       height,   veget,      &
                            frac_nobio, veget_max, reinf_slope,          &
                            ks,         nvan,      avan,     mcr,        &
                            mcs,        mcfc,      mcw,                  &
                            assim_param, frac_age, fraction_aeirrig_sw)
    
    IF (printlev_loc>=3) WRITE (numout,*) 'sechiba_finalize done'
    
  END SUBROUTINE sechiba_finalize

  
!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_init
!!
!>\BRIEF        Dynamic allocation of the variables, the dimensions of the 
!! variables are determined by user-specified settings. 
!! 
!! DESCRIPTION  : The domain size (:: kjpindex) is used to allocate the correct
!! dimensions to all variables in sechiba. Depending on the variable, its 
!! dimensions are also determined by the number of PFT's (::nvm), number of 
!! soil types (::nstm), number of non-vegetative surface types (::nnobio),
!! number of soil levels (::ngrnd), number of soil layers in the hydrological 
!! model (i.e. cwrr) (::nslm). Values for these variables are set in
!! constantes_soil.f90 and constantes_veg.f90.\n
!!
!! Memory is allocated for all Sechiba variables and new indexing tables
!! are build making use of both (::kjpij) and (::kjpindex). New indexing tables 
!! are needed because a single pixel can contain several PFTs, soil types, etc.
!! The new indexing tables have separate indices for the different
!! PFTs, soil types, etc.\n
!!
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): Strictly speaking the subroutine has no output 
!! variables. However, the routine allocates memory and builds new indexing 
!! variables for later use.
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================ 

  SUBROUTINE sechiba_init (kjit, kjpij, kjpindex, index, rest_id, lalo)

!! 0.1 Input variables
 
    INTEGER(i_std), INTENT (in)                         :: kjit               !! Time step number (unitless)
    INTEGER(i_std), INTENT (in)                         :: kjpij              !! Total size of the un-compressed grid (unitless)
    INTEGER(i_std), INTENT (in)                         :: kjpindex           !! Domain size - terrestrial pixels only (unitless)
    INTEGER(i_std), INTENT (in)                         :: rest_id            !! _Restart_ file identifier (unitless)
    INTEGER(i_std),DIMENSION (kjpindex), INTENT (in)    :: index              !! Indeces of the points on the map (unitless)
    REAL(r_std),DIMENSION (kjpindex,2), INTENT (in)     :: lalo               !! Geographical coordinates (latitude,longitude) 
                                                                              !! for pixels (degrees)
!! 0.2 Output variables

!! 0.3 Modified variables

!! 0.4 Local variables

    INTEGER(i_std)                                      :: ier                !! Check errors in memory allocation (unitless)
    INTEGER(i_std)                                      :: ji, jv             !! Indeces (unitless)
!_ ==============================================================================================================================

!! 1. Initialize variables 
    
    ! Dynamic allocation with user-specified dimensions on first call
    IF (l_first_sechiba) THEN 
       l_first_sechiba=.FALSE.
    ELSE 
       CALL ipslerr_p(3,'sechiba_init',' l_first_sechiba false . we stop ','','')
    ENDIF

    !! Initialize local printlev
    printlev_loc=get_printlev('sechiba')
    

    !! 1.1 Initialize 3D vegetation indexation table
    ALLOCATE (indexveg(kjpindex*nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexveg','','')

    ALLOCATE (indexlai(kjpindex*(nlai+1)),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexlai','','')

    ALLOCATE (indexsoil(kjpindex*nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexsoil','','')

    ALLOCATE (indexnobio(kjpindex*nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexnobio','','')

    ALLOCATE (indexgrnd(kjpindex*ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexgrnd','','')

    ALLOCATE (indexsnow(kjpindex*nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexsnow','','')

    ALLOCATE (indexlayer(kjpindex*nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexlayer','','')

    ALLOCATE (indexnslm(kjpindex*nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexnslm','','')

    ALLOCATE (indexalb(kjpindex*2),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for indexalb','','')

    !! 1.2  Initialize 1D array allocation with restartable value
    ALLOCATE (flood_res(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for flood_res','','')
    flood_res(:) = undef_sechiba

    ALLOCATE (flood_frac(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for kjpindex','','')
    flood_frac(:) = undef_sechiba

    ALLOCATE (snow(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snow','','')
    snow(:) = undef_sechiba

    ALLOCATE (snow_age(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snow_age','','')
    snow_age(:) = undef_sechiba

    ALLOCATE (drysoil_frac(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for drysoil_frac','','')

    ALLOCATE (evap_bare_lim(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for evap_bare_lim','','')

    ALLOCATE (evap_bare_lim_ns(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for evap_bare_lim_ns','','')

    ALLOCATE (evapot(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for evapot','','')
    evapot(:) = undef_sechiba

    ALLOCATE (evapot_corr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for evapot_corr','','')

    ALLOCATE (humrel(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for humrel','','')
    humrel(:,:) = undef_sechiba

    ALLOCATE (vegstress(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vegstress','','')
    vegstress(:,:) = undef_sechiba

    ALLOCATE (njsc(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for njsc','','')
    njsc(:)= undef_int

    ALLOCATE (soiltile(kjpindex,nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soiltile','','')

    ALLOCATE (fraclut(kjpindex,nlut),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for fraclut','','')

    ALLOCATE (nwdFraclut(kjpindex,nlut),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for nwdFraclut','','')

    ALLOCATE (reinf_slope(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for reinf_slope','','')

    ALLOCATE (reinf_slope_soil(kjpindex, nstm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for reinf_slope_soil','','') !

    !salma: adding soil params
    ALLOCATE (ks(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ks','','')

    ALLOCATE (nvan(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for nvan ','','')

    ALLOCATE (avan(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for avan','','')

    ALLOCATE (mcr(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcr','','')

    ALLOCATE (mcs(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcs','','')

    ALLOCATE (mcfc(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcfc','','')
    
    ALLOCATE (mcw(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcw','','')
    !end salma: adding soil params




    ALLOCATE (vbeta1(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vbeta1','','')

    ALLOCATE (vbeta4(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vbeta4','','')

    ALLOCATE (vbeta5(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vbeta5','','')

    ALLOCATE (soilcap(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soilcap','','')

    ALLOCATE (soilflx(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soilflx','','')

    ALLOCATE (temp_sol(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for temp_sol','','')
    temp_sol(:) = undef_sechiba

    ALLOCATE (qsurf(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for qsurf','','')
    qsurf(:) = undef_sechiba

    !! 1.3 Initialize 2D array allocation with restartable value
    ALLOCATE (qsintveg(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for qsintveg','','')
    qsintveg(:,:) = undef_sechiba

    ALLOCATE (vbeta2(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vbeta2','','')

    ALLOCATE (vbeta3(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vbeta3','','')

    ALLOCATE (vbeta3pot(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vbeta3pot','','')

    ALLOCATE (gsmean(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for gsmean','','')

    ALLOCATE (cimean(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for cimean','','')

    ALLOCATE (gpp(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for gpp','','')
    gpp(:,:) = undef_sechiba

 
    ALLOCATE (temp_growth(kjpindex),stat=ier) 
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for temp_growth','','')
    temp_growth(:) = undef_sechiba 

    ALLOCATE (veget(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for veget','','')
    veget(:,:)=undef_sechiba

    ALLOCATE (veget_max(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for veget_max','','')

    ALLOCATE (tot_bare_soil(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for tot_bare_soil','','')

    ALLOCATE (lai(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for lai','','')
    lai(:,:)=undef_sechiba

    ALLOCATE (frac_age(kjpindex,nvm,nleafages),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_age','','')
    frac_age(:,:,:)=undef_sechiba

    ALLOCATE (height(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for height','','')
    height(:,:)=undef_sechiba

    ALLOCATE (frac_nobio(kjpindex,nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_nobio','','')
    frac_nobio(:,:) = undef_sechiba

    ALLOCATE (snow_nobio(kjpindex,nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snow_nobio','','')
    snow_nobio(:,:) = undef_sechiba

    ALLOCATE (snow_nobio_age(kjpindex,nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snow_nobio_age','','')
    snow_nobio_age(:,:) = undef_sechiba

    ALLOCATE (assim_param(kjpindex,nvm,npco2),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for assim_param','','')

    !! 1.4 Initialize 1D array allocation 
    ALLOCATE (vevapflo(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapflo','','')
    vevapflo(:)=zero

    ALLOCATE (vevapsno(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapsno','','')

    ALLOCATE (vevapnu(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapnu','','')

    ALLOCATE (totfrac_nobio(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for totfrac_nobio','','')

    ALLOCATE (floodout(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for floodout','','')

    ALLOCATE (runoff(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for runoff','','')

    ALLOCATE (drainage(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for drainage','','')

    ALLOCATE (returnflow(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for returnflow','','')
    returnflow(:) = zero

    ALLOCATE (reinfiltration(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for reinfiltration','','')
    reinfiltration(:) = zero

    ALLOCATE (irrigation(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for irrigation','','')
    irrigation(:) = zero

    ALLOCATE (z0h(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for z0h','','')

    ALLOCATE (z0m(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for z0m','','')

    ALLOCATE (roughheight(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for roughheight','','')

    ALLOCATE (emis(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for emis','','')

    ALLOCATE (tot_melt(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for tot_melt','','')

    ALLOCATE (vbeta(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vbeta','','')

    ALLOCATE (rau(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for rau','','')

    ALLOCATE (deadleaf_cover(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for deadleaf_cover','','')

    ALLOCATE (stempdiag(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for stempdiag','','')

    ALLOCATE (ftempdiag(kjpindex, ngrnd),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ftempdiag','','')

    ALLOCATE (co2_flux(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for co2_flux','','')
    co2_flux(:,:)=zero

    ALLOCATE (shumdiag(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for shumdiag','','')
    
    ALLOCATE (shumdiag_perma(kjpindex,nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for shumdiag_perma','','')

    ALLOCATE (litterhumdiag(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for litterhumdiag','','')

    ALLOCATE (ptnlev1(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for ptnlev1','','')

    ALLOCATE (k_litt(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for k_litt','','')

    !! 1.5 Initialize 2D array allocation
    ALLOCATE (vevapwet(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for vevapwet','','')
    vevapwet(:,:)=undef_sechiba

    ALLOCATE (transpir(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for transpir','','')

    ALLOCATE (transpot(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for transpot','','')

    ALLOCATE (qsintmax(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for qsintmax','','')

    ALLOCATE (rveget(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for rveget','','')

    ALLOCATE (rstruct(kjpindex,nvm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for rstruct','','')

    ALLOCATE (pgflux(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for pgflux','','')
    pgflux(:)= 0.0

    ALLOCATE (cgrnd_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for cgrnd_snow','','')
    cgrnd_snow(:,:) = 0

    ALLOCATE (dgrnd_snow(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for dgrnd_snow','','')
    dgrnd_snow(:,:) = 0

    ALLOCATE (lambda_snow(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for lambda_snow','','')
    lambda_snow(:) = 0

    ALLOCATE (temp_sol_add(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for temp_sol_add','','')

    ALLOCATE (gtemp(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for gtemp','','')

    ALLOCATE (frac_snow_veg(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_snow_veg','','')

    ALLOCATE (frac_snow_nobio(kjpindex,nnobio),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for frac_snow_nobio','','')

    ALLOCATE (snowrho(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowrho','','')

    ALLOCATE (snowheat(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowheat','','')

    ALLOCATE (snowgrain(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowgrain','','')

    ALLOCATE (snowtemp(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowtemp','','')

    ALLOCATE (snowdz(kjpindex,nsnow),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for snowdz','','')

    ALLOCATE (mc_layh(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mc_layh','','')

    ALLOCATE (mcl_layh(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for mcl_layh','','')

    ALLOCATE (soilmoist(kjpindex, nslm),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for soilmoist','','')


    !1.5 Irrigation related variables
    ALLOCATE (root_deficit(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for root_deficit','','') !

    ALLOCATE (irrig_frac(kjpindex),stat=ier)
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for irrig_frac','','') !
    irrigation(:) = zero

    ALLOCATE (irrigated_next(kjpindex),stat=ier) !
    IF (ier /= 0) CALL ipslerr_p(3,'hydrol_init','Problem in allocate of variable irrigated_next','','') !

    ALLOCATE (fraction_aeirrig_sw(kjpindex),stat=ier) !
    IF (ier /= 0) CALL ipslerr_p(3,'sechiba_init','Pb in alloc for fraction_aeirrig_sw','','')

    !! 1.6 Initialize indexing table for the vegetation fields. 
    ! In SECHIBA we work on reduced grids but to store in the full 3D filed vegetation variable 
    ! we need another index table : indexveg, indexsoil, indexnobio and indexgrnd
    DO ji = 1, kjpindex
       !
       DO jv = 1, nlai+1
          indexlai((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, nvm
          indexveg((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !      
       DO jv = 1, nstm
          indexsoil((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !      
       DO jv = 1, nnobio
          indexnobio((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, ngrnd
          indexgrnd((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, nsnow
          indexsnow((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij
       ENDDO

       DO jv = 1, nslm
          indexnslm((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij
       ENDDO

       DO jv = 1, nslm
          indexlayer((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
       DO jv = 1, 2
          indexalb((jv-1)*kjpindex + ji) = INDEX(ji) + (jv-1)*kjpij + offset_omp - offset_mpi
       ENDDO
       !
    ENDDO

!! 2. Read the default value that will be put into variable which are not in the restart file
    CALL ioget_expval(val_exp)
    
    IF (printlev_loc>=3) WRITE (numout,*) ' sechiba_init done '

  END SUBROUTINE sechiba_init
  

!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_clear
!!
!>\BRIEF        Deallocate memory of sechiba's variables
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None 
!!
!! REFERENCE(S)	: None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================ 

  SUBROUTINE sechiba_clear()

!! 1. Initialize first run

    l_first_sechiba=.TRUE.

!! 2. Deallocate dynamic variables of sechiba

    IF ( ALLOCATED (indexveg)) DEALLOCATE (indexveg)
    IF ( ALLOCATED (indexlai)) DEALLOCATE (indexlai)
    IF ( ALLOCATED (indexsoil)) DEALLOCATE (indexsoil)
    IF ( ALLOCATED (indexnobio)) DEALLOCATE (indexnobio)
    IF ( ALLOCATED (indexsnow)) DEALLOCATE (indexsnow)
    IF ( ALLOCATED (indexgrnd)) DEALLOCATE (indexgrnd)
    IF ( ALLOCATED (indexlayer)) DEALLOCATE (indexlayer)
    IF ( ALLOCATED (indexnslm)) DEALLOCATE (indexnslm)
    IF ( ALLOCATED (indexalb)) DEALLOCATE (indexalb)
    IF ( ALLOCATED (flood_res)) DEALLOCATE (flood_res)
    IF ( ALLOCATED (flood_frac)) DEALLOCATE (flood_frac)
    IF ( ALLOCATED (snow)) DEALLOCATE (snow)
    IF ( ALLOCATED (snow_age)) DEALLOCATE (snow_age)
    IF ( ALLOCATED (drysoil_frac)) DEALLOCATE (drysoil_frac)
    IF ( ALLOCATED (evap_bare_lim)) DEALLOCATE (evap_bare_lim)
    IF ( ALLOCATED (evap_bare_lim_ns)) DEALLOCATE (evap_bare_lim_ns)
    IF ( ALLOCATED (evapot)) DEALLOCATE (evapot)
    IF ( ALLOCATED (evapot_corr)) DEALLOCATE (evapot_corr)
    IF ( ALLOCATED (humrel)) DEALLOCATE (humrel)
    IF ( ALLOCATED (vegstress)) DEALLOCATE (vegstress)
    IF ( ALLOCATED (soiltile)) DEALLOCATE (soiltile)
    IF ( ALLOCATED (fraclut)) DEALLOCATE (fraclut)
    IF ( ALLOCATED (nwdFraclut)) DEALLOCATE (nwdFraclut)
    IF ( ALLOCATED (njsc)) DEALLOCATE (njsc)
    IF ( ALLOCATED (reinf_slope)) DEALLOCATE (reinf_slope)
    IF ( ALLOCATED (reinf_slope_soil)) DEALLOCATE (reinf_slope_soil)

    !salma: adding soil hydraulic params
    IF ( ALLOCATED (ks)) DEALLOCATE (ks)
    IF ( ALLOCATED (nvan)) DEALLOCATE (nvan)
    IF ( ALLOCATED (avan)) DEALLOCATE (avan)
    IF ( ALLOCATED (mcr)) DEALLOCATE (mcr)
    IF ( ALLOCATED (mcs)) DEALLOCATE (mcs)
    IF ( ALLOCATED (mcfc)) DEALLOCATE (mcfc)
    IF ( ALLOCATED (mcw)) DEALLOCATE (mcw)
    !end salma: adding soil hydraulic params

    IF ( ALLOCATED (vbeta1)) DEALLOCATE (vbeta1)
    IF ( ALLOCATED (vbeta4)) DEALLOCATE (vbeta4)
    IF ( ALLOCATED (vbeta5)) DEALLOCATE (vbeta5)
    IF ( ALLOCATED (soilcap)) DEALLOCATE (soilcap)
    IF ( ALLOCATED (soilflx)) DEALLOCATE (soilflx)
    IF ( ALLOCATED (temp_sol)) DEALLOCATE (temp_sol)
    IF ( ALLOCATED (qsurf)) DEALLOCATE (qsurf)
    IF ( ALLOCATED (qsintveg)) DEALLOCATE (qsintveg)
    IF ( ALLOCATED (vbeta2))  DEALLOCATE (vbeta2)
    IF ( ALLOCATED (vbeta3)) DEALLOCATE (vbeta3)
    IF ( ALLOCATED (vbeta3pot)) DEALLOCATE (vbeta3pot)
    IF ( ALLOCATED (gsmean)) DEALLOCATE (gsmean)
    IF ( ALLOCATED (cimean)) DEALLOCATE (cimean)
    IF ( ALLOCATED (gpp)) DEALLOCATE (gpp)
    IF ( ALLOCATED (temp_growth)) DEALLOCATE (temp_growth) 
    IF ( ALLOCATED (veget)) DEALLOCATE (veget)
    IF ( ALLOCATED (veget_max)) DEALLOCATE (veget_max)
    IF ( ALLOCATED (tot_bare_soil)) DEALLOCATE (tot_bare_soil)
    IF ( ALLOCATED (lai)) DEALLOCATE (lai)
    IF ( ALLOCATED (frac_age)) DEALLOCATE (frac_age)
    IF ( ALLOCATED (height)) DEALLOCATE (height)
    IF ( ALLOCATED (roughheight)) DEALLOCATE (roughheight)
    IF ( ALLOCATED (frac_nobio)) DEALLOCATE (frac_nobio)
    IF ( ALLOCATED (snow_nobio)) DEALLOCATE (snow_nobio)
    IF ( ALLOCATED (snow_nobio_age)) DEALLOCATE (snow_nobio_age)
    IF ( ALLOCATED (assim_param)) DEALLOCATE (assim_param)
    IF ( ALLOCATED (vevapflo)) DEALLOCATE (vevapflo)
    IF ( ALLOCATED (vevapsno)) DEALLOCATE (vevapsno)
    IF ( ALLOCATED (vevapnu)) DEALLOCATE (vevapnu)
    IF ( ALLOCATED (totfrac_nobio)) DEALLOCATE (totfrac_nobio)
    IF ( ALLOCATED (floodout)) DEALLOCATE (floodout)
    IF ( ALLOCATED (runoff)) DEALLOCATE (runoff)
    IF ( ALLOCATED (drainage)) DEALLOCATE (drainage)
    IF ( ALLOCATED (reinfiltration)) DEALLOCATE (reinfiltration)
    IF ( ALLOCATED (irrigation)) DEALLOCATE (irrigation)
    IF ( ALLOCATED (tot_melt)) DEALLOCATE (tot_melt)
    IF ( ALLOCATED (vbeta)) DEALLOCATE (vbeta)
    IF ( ALLOCATED (rau)) DEALLOCATE (rau)
    IF ( ALLOCATED (deadleaf_cover)) DEALLOCATE (deadleaf_cover)
    IF ( ALLOCATED (stempdiag)) DEALLOCATE (stempdiag)
    IF ( ALLOCATED (ftempdiag)) DEALLOCATE (ftempdiag)
    IF ( ALLOCATED (co2_flux)) DEALLOCATE (co2_flux)
    IF ( ALLOCATED (shumdiag)) DEALLOCATE (shumdiag)
    IF ( ALLOCATED (shumdiag_perma)) DEALLOCATE (shumdiag_perma)
    IF ( ALLOCATED (litterhumdiag)) DEALLOCATE (litterhumdiag)
    IF ( ALLOCATED (ptnlev1)) DEALLOCATE (ptnlev1)
    IF ( ALLOCATED (k_litt)) DEALLOCATE (k_litt)
    IF ( ALLOCATED (vevapwet)) DEALLOCATE (vevapwet)
    IF ( ALLOCATED (transpir)) DEALLOCATE (transpir)
    IF ( ALLOCATED (transpot)) DEALLOCATE (transpot)
    IF ( ALLOCATED (qsintmax)) DEALLOCATE (qsintmax)
    IF ( ALLOCATED (rveget)) DEALLOCATE (rveget)
    IF ( ALLOCATED (rstruct)) DEALLOCATE (rstruct)
    IF ( ALLOCATED (frac_snow_veg)) DEALLOCATE (frac_snow_veg)
    IF ( ALLOCATED (frac_snow_nobio)) DEALLOCATE (frac_snow_nobio)
    IF ( ALLOCATED (snowrho)) DEALLOCATE (snowrho)
    IF ( ALLOCATED (snowgrain)) DEALLOCATE (snowgrain)
    IF ( ALLOCATED (snowtemp)) DEALLOCATE (snowtemp)
    IF ( ALLOCATED (snowdz)) DEALLOCATE (snowdz)
    IF ( ALLOCATED (snowheat)) DEALLOCATE (snowheat)
    IF ( ALLOCATED (cgrnd_snow)) DEALLOCATE (cgrnd_snow)
    IF ( ALLOCATED (dgrnd_snow)) DEALLOCATE (dgrnd_snow)
    IF ( ALLOCATED (lambda_snow)) DEALLOCATE(lambda_snow) 
    IF ( ALLOCATED (gtemp)) DEALLOCATE (gtemp)
    IF ( ALLOCATED (pgflux)) DEALLOCATE (pgflux)
    IF ( ALLOCATED (mc_layh)) DEALLOCATE (mc_layh)
    IF ( ALLOCATED (mcl_layh)) DEALLOCATE (mcl_layh)
    IF ( ALLOCATED (soilmoist)) DEALLOCATE (soilmoist)
    IF ( ALLOCATED (root_deficit)) DEALLOCATE (root_deficit)
    IF ( ALLOCATED (irrig_frac)) DEALLOCATE (irrig_frac)
    IF ( ALLOCATED (irrigated_next)) DEALLOCATE (irrigated_next)

!! 3. Clear all allocated memory

    CALL pft_parameters_clear
    CALL slowproc_clear 
    CALL diffuco_clear 
    CALL enerbil_clear  
    CALL hydrol_clear 
    CALL thermosoil_clear
    CALL condveg_clear 
    CALL routing_wrapper_clear

  END SUBROUTINE sechiba_clear


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_var_init
!!
!>\BRIEF        Calculate air density as a function of air temperature and 
!! pressure for each terrestrial pixel.
!! 
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): air density (::rau, kg m^{-3}).
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE sechiba_var_init (kjpindex, rau, pb, temp_air) 

!! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                    :: kjpindex        !! Domain size - terrestrial pixels only (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: pb              !! Surface pressure (hPa)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)  :: temp_air        !! Air temperature (K)
    
!! 0.2 Output variables

    REAL(r_std),DIMENSION (kjpindex), INTENT (out) :: rau             !! Air density @tex $(kg m^{-3})$ @endtex

!! 0.3 Modified variables

!! 0.4 Local variables

    INTEGER(i_std)                                 :: ji              !! Indices (unitless)
!_ ================================================================================================================================
    
!! 1. Calculate intial air density (::rau)
   
    DO ji = 1,kjpindex
       rau(ji) = pa_par_hpa * pb(ji) / (cte_molr*temp_air(ji))
    END DO

    IF (printlev_loc>=3) WRITE (numout,*) ' sechiba_var_init done '

  END SUBROUTINE sechiba_var_init


!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_end
!!
!>\BRIEF        Swap old for newly calculated soil temperature.
!! 
!! RECENT CHANGE(S): None
!! 
!! MAIN OUTPUT VARIABLE(S): soil temperature (::temp_sol; K)
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!! ================================================================================================================================ 

  SUBROUTINE sechiba_end (kjpindex, temp_sol_new, temp_sol)
                         

!! 0.1 Input variables

    INTEGER(i_std), INTENT (in)                       :: kjpindex           !! Domain size - terrestrial pixels only (unitless)
    REAL(r_std),DIMENSION (kjpindex), INTENT (in)     :: temp_sol_new       !! New soil temperature (K)
    
    !! 0.2 Output variables
    REAL(r_std),DIMENSION (kjpindex), INTENT (out)    :: temp_sol           !! Soil temperature (K)

!_ ================================================================================================================================
    
!! 1. Swap temperature

    temp_sol(:) = temp_sol_new(:)
    
    IF (printlev_loc>=3) WRITE (numout,*) ' sechiba_end done '

  END SUBROUTINE sechiba_end

!! ==============================================================================================================================\n
!! SUBROUTINE 	: sechiba_interface_orchidee_inca
!!
!>\BRIEF        make the interface between surface and atmospheric chemistry
!! 
!! DESCRIPTION  : This subroutine is called from INCA, the atmospheric chemistry model. It is used to transfer variables from ORCHIDEE to INCA. 
!!
!! RECENT CHANGE(S): move from chemistry module to be more generic (feb - 2017)
!! 
!! MAIN OUTPUT VARIABLE(S): emission COV to be transport by orchidee to inca in fields_out array 
!! 
!! REFERENCE(S)	: None
!! 
!! FLOWCHART    : None
!! \n
!! ================================================================================================================================ 
  SUBROUTINE sechiba_interface_orchidee_inca( &
       nvm_out, veget_max_out, veget_frac_out, lai_out, snow_out, &
       field_out_names, fields_out, field_in_names, fields_in)


    INTEGER, INTENT(out)                      :: nvm_out            !! Number of vegetation types
    REAL(r_std), DIMENSION (:,:), INTENT(out) :: veget_max_out      !! Max. fraction of vegetation type (LAI -> infty)
    REAL(r_std), DIMENSION (:,:), INTENT(out) :: veget_frac_out     !! Fraction of vegetation type (unitless, 0-1)  
    REAL(r_std), DIMENSION (:,:), INTENT(out) :: lai_out            !! Surface foliere
    REAL(r_std), DIMENSION (:)  , INTENT(out) :: snow_out           !! Snow mass [Kg/m^2]

    !
    ! Optional arguments
    !
    ! Names and fields for emission variables : to be transport by Orchidee to Inca
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN) :: field_out_names
    REAL(r_std),DIMENSION(:,:,:), OPTIONAL, INTENT(OUT) :: fields_out
    !
    ! Names and fields for deposit variables : to be transport from chemistry model by INCA to ORCHIDEE.
    CHARACTER(LEN=*),DIMENSION(:), OPTIONAL, INTENT(IN) :: field_in_names
    REAL(r_std),DIMENSION(:,:), OPTIONAL, INTENT(IN)    :: fields_in


    ! Variables always transmitted from sechiba to inca
    nvm_out = nvm 
    veget_max_out(:,:)  = veget_max(:,:) 
    veget_frac_out(:,:) = veget(:,:) 
    lai_out(:,:)  = lai(:,:) 
    snow_out(:)  = snow(:) 

    ! Call chemistry_flux_interface if at least one of variables field_out_names or
    ! field_in_names is present in the argument list of sechiba_interface_orchidee_inca when called from inca.
    IF (PRESENT(field_out_names) .AND. .NOT. PRESENT(field_in_names)) THEN 
       CALL chemistry_flux_interface(field_out_names=field_out_names, fields_out=fields_out)
    ELSE IF (.NOT. PRESENT(field_out_names) .AND. PRESENT(field_in_names)) THEN 
       CALL chemistry_flux_interface(field_in_names=field_in_names, fields_in=fields_in)
    ELSE IF (PRESENT(field_out_names) .AND. PRESENT(field_in_names)) THEN 
       CALL chemistry_flux_interface(field_out_names=field_out_names, fields_out=fields_out, &
            field_in_names=field_in_names, fields_in=fields_in)
    ENDIF

  END SUBROUTINE sechiba_interface_orchidee_inca


END MODULE sechiba

