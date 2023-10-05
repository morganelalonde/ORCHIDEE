! =================================================================================================================================
! MODULE 	: constantes_soil_var
!
! CONTACT       : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE       : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF         "constantes_soil_var" module contains the parameters related to soil and hydrology.
!!
!!\n DESCRIPTION : The non saturated hydraulic properties are defined from the  
!!                 formulations of van Genuchten (1980) and Mualem (1976), combined as  
!!                 explained in d'Orgeval (2006). \n
!!                 The related parameters for main soil textures (coarse, medium and fine if "fao", 
!!                 12 USDA testures if "usda") come from Carsel and Parrish (1988).
!!
!! RECENT CHANGE(S): AD: mcw and mcf depend now on soil texture, based on Van Genuchten equations 
!!                   and classical matric potential values, and pcent is adapted
!!                   November 2020 by Salma Tafasca and Agnes Ducharne : we introduce a new texture class
!!                   for clay oxisols (cf. Tafasca, 2020, PhD thesis; Tafasca et al., in prep for GRL).
!!                   It makes no change if we read a soil texture map with only 12 USDA classes.
!!                   Lookup tables for Zobler replaces by pointer to read the corresponding values in the
!!                   13-value USDA tables
!!
!! REFERENCE(S)	:
!!- Roger A.Pielke, (2002), Mesoscale meteorological modeling, Academic Press Inc. 
!!- Polcher, J., Laval, K., Dümenil, L., Lean, J., et Rowntree, P. R. (1996).
!! Comparing three land surface schemes used in general circulation models. Journal of Hydrology, 180(1-4), 373--394.
!!- Ducharne, A., Laval, K., et Polcher, J. (1998). Sensitivity of the hydrological cycle
!! to the parametrization of soil hydrology in a GCM. Climate Dynamics, 14, 307--327. 
!!- Rosnay, P. de et Polcher, J. (1999). Modelling root water uptake in a complex land surface
!! scheme coupled to a GCM. Hydrol. Earth Syst. Sci., 2(2/3), 239--255.
!!- d'Orgeval, T. et Polcher, J. (2008). Impacts of precipitation events and land-use changes
!! on West African river discharges during the years 1951--2000. Climate Dynamics, 31(2), 249--262. 
!!- Carsel, R. and Parrish, R.: Developing joint probability distributions of soil water
!! retention characteristics, Water Resour. Res.,24, 755–769, 1988.
!!- Mualem Y (1976). A new model for predicting the hydraulic conductivity  
!! of unsaturated porous media. Water Resources Research 12(3):513-522
!!- Van Genuchten M (1980). A closed-form equation for predicting the  
!! hydraulic conductivity of unsaturated soils. Soil Sci Soc Am J, 44(5):892-898
!!- Tafasca S. (2020). Evaluation de l’impact des propriétés du sol sur l’hydrologie simulee dans le
!! modèle ORCHIDEE, PhD thesis, Sorbonne Universite. 
!!- Tafasca S., Ducharne A. and Valentin C. Accounting for soil structure in pedo-transfer functions:
!!  swelling vs non swelling clays. In prep for GRL.
!!
!! SVN          :
!! $HeadURL: $
!! $Date: $
!! $Revision: $
!! \n
!_ ================================================================================================================================

MODULE constantes_soil_var

  USE defprec
  USE vertical_soil_var

  IMPLICIT NONE

  LOGICAL, SAVE             :: check_cwrr               !! Calculate diagnostics to check the water balance in hydrol (true/false)
!$OMP THREADPRIVATE(check_cwrr)

  !! Number of soil classes

  INTEGER(i_std), PARAMETER :: ntext=3                  !! Number of soil textures (Silt, Sand, Clay)
  INTEGER(i_std), PARAMETER :: nstm=3                   !! Number of soil tiles (unitless)
  CHARACTER(LEN=30)         :: soil_classif             !! Type of classification used for the map of soil types.
                                                        !! It must be consistent with soil file given by 
                                                        !! SOILCLASS_FILE parameter.
!$OMP THREADPRIVATE(soil_classif)
  INTEGER(i_std), PARAMETER :: nscm_fao=3               !! For FAO Classification (unitless)
  INTEGER(i_std), PARAMETER :: nscm_usda=13             !! For USDA Classification (unitless)
  INTEGER(i_std), SAVE      :: nscm=nscm_usda            !! Default value for nscm
!$OMP THREADPRIVATE(nscm)

  !! Parameters for soil thermodynamics
  REAL(r_std), SAVE :: sn_cond = 0.3                    !! Thermal Conductivity of snow 
                                                        !! @tex $(W.m^{-2}.K^{-1})$ @endtex  
!$OMP THREADPRIVATE(sn_cond)
  REAL(r_std), SAVE :: sn_dens = 330.0                  !! Snow density for the soil thermodynamics
                                                        !! (kg/m3)
!$OMP THREADPRIVATE(sn_dens)
  REAL(r_std), SAVE :: sn_capa                          !! Volumetric heat capacity for snow 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(sn_capa)
  REAL(r_std), PARAMETER :: capa_ice = 2.228*1.E3       !! Specific heat capacity of ice (J/kg/K)

  REAL(r_std), SAVE :: water_capa = 4.18e+6             !! Volumetric water heat capacity 
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(water_capa)
  REAL(r_std), SAVE :: brk_capa = 2.0e+6                !! Volumetric heat capacity of generic rock
                                                        !! @tex $(J.m^{-3}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(brk_capa)
  REAL(r_std), SAVE :: brk_cond = 3.0                   !! Thermal conductivity of saturated granitic rock
                                                        !! @tex $(W.m^{-1}.K^{-1})$ @endtex
!$OMP THREADPRIVATE(brk_cond)

  REAL(r_std), SAVE :: qsintcst = 0.02                  !! Transforms leaf area index into size of interception reservoir
                                                        !! (unitless)
!$OMP THREADPRIVATE(qsintcst)
  REAL(r_std), SAVE :: mx_eau_nobio = 150.              !! Volumetric available soil water capacity in nobio fractions
                                                        !! @tex $(kg.m^{-3} of soil)$ @endtex
!$OMP THREADPRIVATE(mx_eau_nobio)

  !! Parameters specific for the CWRR hydrology.

  !!  1. Parameters for FAO-Zobler Map

  INTEGER(i_std), PARAMETER,DIMENSION(nscm_fao) :: fao2usda = (/ 3,6,9 /) !! To find the values of Coarse, Medium, Fine in Zobler map
                                                                          !! from the USDA lookup tables
  
  !!  2. Parameters for USDA Classification

  !! Parameters for soil type distribution :
  !! Sand, Loamy Sand, Sandy Loam, Silt Loam, Silt, Loam, Sandy Clay Loam, Silty Clay Loam, Clay Loam, Sandy Clay, Silty Clay, Clay

  INTEGER(i_std), SAVE      :: usda_default = 6            !! Default USDA texture class if no value found from map
!$OMP THREADPRIVATE(usda_default)

  REAL(r_std), PARAMETER, DIMENSION(nscm_usda) :: soilclass_default = (/0.0, &
       0.0, 0.0, 0.0, 0.0, 1.0, 0.0, &  !! Areal fraction of the 13 soil USDA textures; 
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)   !! the dominant one will selected
  
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: nvan_usda = &            !! Van Genuchten coefficient n (unitless)
 & (/ 2.68_r_std, 2.28_r_std, 1.89_r_std, 1.41_r_std, &                   !  RK: 1/n=1-m
 &    1.37_r_std, 1.56_r_std, 1.48_r_std, 1.23_r_std, &
 &    1.31_r_std, 1.23_r_std, 1.09_r_std, 1.09_r_std, &
 &    1.552_r_std    /) ! oxisols

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: avan_usda = &            !! Van Genuchten coefficient a 
 & (/ 0.0145_r_std, 0.0124_r_std, 0.0075_r_std, 0.0020_r_std, &          !!  @tex $(mm^{-1})$ @endtex
 &    0.0016_r_std, 0.0036_r_std, 0.0059_r_std, 0.0010_r_std, &
 &    0.0019_r_std, 0.0027_r_std, 0.0005_r_std, 0.0008_r_std, &
 &    0.0132_r_std /) ! oxisols

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcr_usda = &             !! Residual volumetric water content 
 & (/ 0.045_r_std, 0.057_r_std, 0.065_r_std, 0.067_r_std, &              !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.034_r_std, 0.078_r_std, 0.100_r_std, 0.089_r_std, &
 &    0.095_r_std, 0.100_r_std, 0.070_r_std, 0.068_r_std, &
 &    0.068_r_std /) ! oxisols

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcs_usda = &             !! Saturated volumetric water content 
 & (/ 0.43_r_std, 0.41_r_std, 0.41_r_std, 0.45_r_std, &                  !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.46_r_std, 0.43_r_std, 0.39_r_std, 0.43_r_std, &
 &    0.41_r_std, 0.38_r_std, 0.36_r_std, 0.38_r_std, &
 &    0.503_r_std  /) ! oxisols

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: ks_usda = &              !! Hydraulic conductivity at saturation
 & (/ 7128.0_r_std, 3501.6_r_std, 1060.8_r_std, 108.0_r_std, &           !!  @tex $(mm d^{-1})$ @endtex
 &    60.0_r_std, 249.6_r_std, 314.4_r_std, 16.8_r_std, &
 &    62.4_r_std, 28.8_r_std, 4.8_r_std, 48.0_r_std, &
 &    6131.4_r_std  /) ! oxisols

! The max available water content is smaller when mcw and mcf depend on texture,
! so we increase pcent to a classical value of 80%  
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: pcent_usda = &           !! Fraction of saturated volumetric soil moisture
 & (/ 0.8_r_std, 0.8_r_std, 0.8_r_std, 0.8_r_std, &                      !! above which transpir is max (0-1, unitless)
 &    0.8_r_std, 0.8_r_std, 0.8_r_std, 0.8_r_std, &
 &    0.8_r_std, 0.8_r_std, 0.8_r_std, 0.8_r_std, &
 &    0.8_r_std /) ! oxisols

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: free_drain_max_usda = &  !! Max=default value of the permeability coeff 
 & (/ 1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std, &                      !! at the bottom of the soil (0-1, unitless)
 &    1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std, &
 &    1.0_r_std, 1.0_r_std, 1.0_r_std, 1.0_r_std,  &
 &    1.0_r_std /)
  
!! We use the VG relationships to derive mcw and mcf depending on soil texture
!! assuming that the matric potential for wilting point and field capacity is
!! -150m (permanent WP) and -3.3m respectively
!! (-1m for FC for the three sandy soils following Richards, L.A. and Weaver, L.R. (1944)
!! Note that mcw GE mcr
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcf_usda = &             !! Volumetric water content at field capacity
 & (/ 0.0493_r_std, 0.0710_r_std, 0.1218_r_std, 0.2402_r_std, &          !!  @tex $(m^{3} m^{-3})$ @endtex
      0.2582_r_std, 0.1654_r_std, 0.1695_r_std, 0.3383_r_std, &
      0.2697_r_std, 0.2672_r_std, 0.3370_r_std, 0.3469_r_std, &
      0.172_r_std  /) ! oxisols
  
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mcw_usda = &             !! Volumetric water content at wilting point
 & (/ 0.0450_r_std, 0.0570_r_std, 0.0657_r_std, 0.1039_r_std, &          !!  @tex $(m^{3} m^{-3})$ @endtex
      0.0901_r_std, 0.0884_r_std, 0.1112_r_std, 0.1967_r_std, &
      0.1496_r_std, 0.1704_r_std, 0.2665_r_std, 0.2707_r_std, &
      0.075_r_std  /) ! oxisols

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mc_awet_usda = &         !! Vol. wat. cont. above which albedo is cst
 & (/ 0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &                  !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &
 &    0.25_r_std, 0.25_r_std, 0.25_r_std, 0.25_r_std, &
 &    0.25_r_std /)

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: mc_adry_usda = &         !! Vol. wat. cont. below which albedo is cst
 & (/ 0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &                      !!  @tex $(m^{3} m^{-3})$ @endtex
 &    0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &
 &    0.1_r_std, 0.1_r_std, 0.1_r_std, 0.1_r_std, &
 &    0.1_r_std /) ! oxisols
 
  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: QZ_usda = &              !! QUARTZ CONTENT (SOIL TYPE DEPENDENT)
 & (/ 0.92_r_std, 0.82_r_std, 0.60_r_std, 0.25_r_std, &                  !! Peters et al [1998]
 &    0.10_r_std, 0.40_r_std, 0.60_r_std, 0.10_r_std, &
 &    0.35_r_std, 0.52_r_std, 0.10_r_std, 0.25_r_std, &
&     0.25_r_std /)  ! oxisols                 

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: so_capa_dry_usda = &     !! Dry soil Volumetric Heat capacity of soils,J.m^{-3}.K^{-1}
 & (/ 1.47e+6_r_std, 1.41e+6_r_std, 1.34e+6_r_std, 1.27e+6_r_std, &      !! Pielke [2002, 2013]
 &    1.21e+6_r_std, 1.21e+6_r_std, 1.18e+6_r_std, 1.32e+6_r_std, &
 &    1.23e+6_r_std, 1.18e+6_r_std, 1.15e+6_r_std, 1.09e+6_r_std, &
 &    1.09e+6_r_std /) ! oxisols

  REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: clayfrac_usda = &   !! % clay particles in the 13 USDA texture classes
       (/ 0.03_r_std, 0.06_r_std, 0.11_r_std, 0.19_r_std , &        !! values taken from get_soilcorr_usda in slowproc
          0.10_r_std, 0.20_r_std, 0.27_r_std, 0.33_r_std, &
          0.33_r_std, 0.41_r_std, 0.46_r_std, 0.55_r_std, &
          0.55_r_std /) ! oxisols                                                                  

   REAL(r_std),PARAMETER,DIMENSION(nscm_usda) :: sandfrac_usda = &  !! % sand particles in the 13 USDA texture classes 
        (/ 0.93_r_std, 0.81_r_std, 0.63_r_std, 0.17_r_std, &        !! values taken from get_soilcorr_usda in slowproc   
           0.06_r_std, 0.40_r_std, 0.54_r_std, 0.08_r_std, &
           0.30_r_std, 0.48_r_std, 0.06_r_std, 0.15_r_std, &
           0.15_r_std /) ! oxisols
   
  !! Parameters for the numerical scheme used by CWRR

  INTEGER(i_std), PARAMETER :: imin = 1                                 !! Start for CWRR linearisation (unitless)
  INTEGER(i_std), PARAMETER :: nbint = 50                               !! Number of interval for CWRR linearisation (unitless)
  INTEGER(i_std), PARAMETER :: imax = nbint+1                           !! Number of points for CWRR linearisation (unitless)
  REAL(r_std), PARAMETER    :: w_time = 1.0_r_std                       !! Time weighting for CWRR numerical integration (unitless)


  !! Variables related to soil freezing, in thermosoil : 
  LOGICAL, SAVE        :: ok_Ecorr                    !! Flag for energy conservation correction
!$OMP THREADPRIVATE(ok_Ecorr)
  LOGICAL, SAVE        :: ok_freeze_thermix           !! Flag to activate thermal part of the soil freezing scheme
!$OMP THREADPRIVATE(ok_freeze_thermix)
  LOGICAL, SAVE        :: ok_freeze_thaw_latent_heat  !! Flag to activate latent heat part of the soil freezing scheme
!$OMP THREADPRIVATE(ok_freeze_thaw_latent_heat)
  LOGICAL, SAVE        :: read_reftemp                !! Flag to initialize soil temperature using climatological temperature
!$OMP THREADPRIVATE(read_reftemp)
  REAL(r_std), SAVE    :: fr_dT                       !! Freezing window (K)
!$OMP THREADPRIVATE(fr_dT)

  !! Variables related to soil freezing, in hydrol : 
  LOGICAL, SAVE        :: ok_freeze_cwrr              !! CWRR freezing scheme by I. Gouttevin
!$OMP THREADPRIVATE(ok_freeze_cwrr)
  LOGICAL, SAVE        :: ok_thermodynamical_freezing !! Calculate frozen fraction thermodynamically
!$OMP THREADPRIVATE(ok_thermodynamical_freezing)

  
END MODULE constantes_soil_var
