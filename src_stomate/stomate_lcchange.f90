! =================================================================================================================================
! MODULE       : stomate_lcchange
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2006)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF       Impact of land cover change on carbon stocks
!!
!!\n DESCRIPTION: None
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: svn://forge.ipsl.jussieu.fr/orchidee/branches/ORCHIDEE_2_2/ORCHIDEE/src_stomate/stomate_lcchange.f90 $
!! $Date: 2021-10-20 18:39:22 +0200 (Wed, 20 Oct 2021) $
!! $Revision: 7326 $
!! \n
!_ ================================================================================================================================


MODULE stomate_lcchange

  ! modules used:
  
  USE ioipsl_para
  USE stomate_data
  USE pft_parameters
  USE constantes
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC lcchange_main
  
CONTAINS


!! ================================================================================================================================
!! SUBROUTINE   : lcchange_main
!!
!>\BRIEF        Impact of land cover change on carbon stocks
!!
!! DESCRIPTION  : This subroutine is always activate if VEGET_UPDATE>0Y in the configuration file, which means that the 
!! vegetation map is updated regulary. lcchange_main is called from stomateLpj the first time step after the vegetation 
!! map has been changed. 
!! The impact of land cover change on carbon stocks is computed in this subroutine. The land cover change is written
!! by the difference of current and previous "maximal" coverage fraction of a PFT. 
!! On the basis of this difference, the amount of 'new establishment'/'biomass export',
!! and increase/decrease of each component, are estimated.\n
!!
!! Main structure of lpj_establish.f90 is:
!! 1. Initialization
!! 2. Calculation of changes in carbon stocks and biomass by land cover change
!! 3. Update 10 year- and 100 year-turnover pool contents
!! 4. History
!!
!! RECENT CHANGE(S) : None
!!
!! MAIN OUTPUT VARIABLE(S) : ::prod10, ::prod100, ::flux10, ::flux100,
!!   :: cflux_prod10 and :: cflux_prod100 
!!
!! REFERENCES   : None
!!
!! FLOWCHART    : 
!! \latexonly 
!!     \includegraphics[scale=0.5]{lcchange.png}
!! \endlatexonly
!! \n
!_ ================================================================================================================================

  
  SUBROUTINE lcchange_main ( npts, dt_days, veget_cov_max_old, veget_cov_max_new, &
       biomass, ind, age, PFTpresent, senescence, when_growthinit, everywhere, &        
       co2_to_bm, bm_to_litter, turnover_daily, bm_sapl, cn_ind,flux10,flux100, &
       prod10,prod100,convflux,cflux_prod10,cflux_prod100,leaf_frac,&
       npp_longterm, lm_lastyearmax, litter, carbon,&
       convfluxpft, fDeforestToProduct, fLulccResidue)

    
    IMPLICIT NONE
    
  !! 0. Variable and parameter declaration 
    
    !! 0.1 Input variables
    
    INTEGER, INTENT(in)                                       :: npts             !! Domain size - number of pixels (unitless)
    REAL(r_std), INTENT(in)                                   :: dt_days          !! Time step of vegetation dynamics for stomate
                                                                                  !! (days)
    REAL(r_std), DIMENSION(nvm, nparts,nelements), INTENT(in) :: bm_sapl          !! biomass of sapling 
                                                                                  !! @tex ($gC individual^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_old!! Current "maximal" coverage fraction of a PFT (LAI
                                                                                  !! -> infinity) on ground
    REAL(r_std), DIMENSION(npts,nvm), INTENT(in)              :: veget_cov_max_new!! New "maximal" coverage fraction of a PFT (LAI ->
                                                                                  !! infinity) on ground (unitless) 
 
    !! 0.2 Output variables

    REAL(r_std), DIMENSION(npts), INTENT(out)                 :: convflux         !! release during first year following land cover
                                                                                  !! change
    REAL(r_std), DIMENSION(npts), INTENT(out)                 :: cflux_prod10     !! total annual release from the 10 year-turnover
                                                                                  !! pool @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts), INTENT(out)                 :: cflux_prod100    !! total annual release from the 100 year-
                                                                                  !! turnover pool @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: turnover_daily   !! Turnover rates 
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: convfluxpft      !! release during first year following land cover                                       
                                                                                  !! change   
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: fDeforestToProduct !!  Deforested biomass into product pool due to anthorpogenic 
                                                                                    !! land use change
    REAL(r_std), DIMENSION(npts,nvm), INTENT(out)             :: fLulccResidue      !!  carbon mass flux into soil and litter due to anthropogenic land use or land cover change

    !! 0.3 Modified variables   
    
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout):: biomass    !! biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: ind              !! Number of individuals @tex ($m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: age              !! mean age (years)
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: senescence       !! plant senescent (only for deciduous trees) Set
                                                                                  !! to .FALSE. if PFT is introduced or killed
    LOGICAL, DIMENSION(npts,nvm), INTENT(inout)               :: PFTpresent       !! Is pft there (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: everywhere       !! is the PFT everywhere in the grid box or very 
                                                                                  !! localized (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: when_growthinit  !! how many days ago was the beginning of the 
                                                                                  !! growing season (days)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: co2_to_bm        !! biomass uptaken 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm,nparts,nelements), INTENT(inout) :: bm_to_litter !! conversion of biomass to litter 
                                                                                  !! @tex ($gC m^{-2} day^{-1}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: cn_ind           !! crown area of individuals 
                                                                                  !! @tex ($m^{2}$) @endtex
    REAL(r_std), DIMENSION(npts,0:10), INTENT(inout)          :: prod10           !! products remaining in the 10 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (10 + 1 : input from year of land
                                                                                  !! cover change)
    REAL(r_std), DIMENSION(npts,0:100), INTENT(inout)         :: prod100          !! products remaining in the 100 year-turnover
                                                                                  !! pool after the annual release for each 
                                                                                  !! compartment (100 + 1 : input from year of land
                                                                                  !! cover change)
    REAL(r_std), DIMENSION(npts,10), INTENT(inout)            :: flux10           !! annual release from the 10/100 year-turnover 
                                                                                  !! pool compartments
    REAL(r_std), DIMENSION(npts,100), INTENT(inout)           :: flux100          !! annual release from the 10/100 year-turnover
                                                                                  !! pool compartments
    REAL(r_std), DIMENSION(npts,nvm,nleafages), INTENT(inout) :: leaf_frac        !! fraction of leaves in leaf age class 
                                                                                  !! (unitless)
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: lm_lastyearmax   !! last year's maximum leaf mass for each PFT 
                                                                                  !! @tex ($gC m^{-2}$) @endtex
    REAL(r_std), DIMENSION(npts,nvm), INTENT(inout)           :: npp_longterm     !! "long term" net primary productivity 
                                                                                  !! @tex ($gC m^{-2} year^{-1}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nvm,nlevs,nelements), INTENT(inout):: litter !! metabolic and structural litter, above and 
                                                                                  !! below ground @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ncarb,nvm), INTENT(inout)      :: carbon           !! carbon pool: active, slow, or passive 

                                                                                 !! @tex ($gC m^{-2}$) @endtex

    !! 0.4 Local variables

    INTEGER(i_std)                                            :: i, j, k, l, m    !! indices (unitless)
    REAL(r_std),DIMENSION(npts,nelements)                     :: bm_new           !! biomass increase @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nparts,nelements)              :: biomass_loss     !! biomass loss @tex ($gC m^{-2}$) @endtex
    REAL(r_std)                                               :: above            !! aboveground biomass @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,nlitt,nlevs,nelements)         :: dilu_lit         !! Litter dilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(npts,ncarb)                         :: dilu_soil_carbon !! Soil Carbondilution @tex ($gC m^{-2}$) @endtex
    REAL(r_std),DIMENSION(nvm)                                :: delta_veg        !! changes in "maximal" coverage fraction of PFT 
    REAL(r_std)                                               :: delta_veg_sum    !! sum of delta_veg
    REAL(r_std),DIMENSION(npts,nvm)                           :: delta_ind        !! change in number of individuals  

!_ ================================================================================================================================

    IF (printlev>=3) WRITE(numout,*) 'Entering lcchange_main'
    
  !! 1. initialization
    
    prod10(:,0)         = zero
    prod100(:,0)        = zero   
    above               = zero
    convflux(:)         = zero
    convfluxpft(:,:)    = zero
    cflux_prod10(:)     = zero
    cflux_prod100(:)    = zero
    delta_ind(:,:)      = zero
    delta_veg(:)        = zero
    fDeforestToProduct(:,:)  = zero
    fLulccResidue(:,:)       = zero


    
  !! 3. calculation of changes in carbon stocks and biomass by land cover change\n
    
    DO i = 1, npts ! Loop over # pixels - domain size
       
       !! 3.1 initialization of carbon stocks\n
       delta_veg(:) = veget_cov_max_new(i,:)-veget_cov_max_old(i,:)
       delta_veg_sum = SUM(delta_veg,MASK=delta_veg.LT.0.)
       
       dilu_lit(i,:,:,:) = zero
       dilu_soil_carbon(i,:) = zero
       biomass_loss(i,:,:) = zero
       
       !! 3.2 if vegetation coverage decreases, compute dilution of litter, soil carbon, and biomass.\n
       DO j=1, nvm
          IF ( delta_veg(j) < -min_stomate ) THEN 
             dilu_lit(i,:,:,:) = dilu_lit(i,:,:,:) + delta_veg(j)*litter(i,:,j,:,:) / delta_veg_sum
             dilu_soil_carbon(i,:) =  dilu_soil_carbon(i,:) + delta_veg(j) * carbon(i,:,j) / delta_veg_sum
             biomass_loss(i,:,:) = biomass_loss(i,:,:) + biomass(i,j,:,:)*delta_veg(j) / delta_veg_sum
          ENDIF
       ENDDO
       
       !! 3.3 
       DO j=1, nvm ! Loop over # PFTs

          !! 3.3.1 The case that vegetation coverage of PFTj increases
          IF ( delta_veg(j) > min_stomate) THEN

             !! 3.3.1.1 Initial setting of new establishment
             IF (veget_cov_max_old(i,j) .LT. min_stomate) THEN 
                IF (is_tree(j)) THEN

                   ! cn_sapl(j)=0.5; stomate_data.f90
                   cn_ind(i,j) = cn_sapl(j) 
                ELSE
                   cn_ind(i,j) = un
                ENDIF
                ind(i,j)= delta_veg(j) / cn_ind(i,j)
                PFTpresent(i,j) = .TRUE.
                everywhere(i,j) = 1.
                senescence(i,j) = .FALSE.
                age(i,j) = zero

                ! large_value = 1.E33_r_std
                when_growthinit(i,j) = large_value 
                leaf_frac(i,j,1) = 1.0
                npp_longterm(i,j) = npp_longterm_init
                lm_lastyearmax(i,j) = bm_sapl(j,ileaf,icarbon) * ind(i,j)
             ENDIF
             IF ( cn_ind(i,j) > min_stomate ) THEN
                delta_ind(i,j) = delta_veg(j) / cn_ind(i,j) 
             ENDIF
             
             !! 3.3.1.2 Update of biomass in each each carbon stock component 
             !!         Update of biomass in each each carbon stock component (leaf, sapabove, sapbelow,
             !>         heartabove, heartbelow, root, fruit, and carbres)\n
             DO k = 1, nparts ! loop over # carbon stock components, nparts = 8; stomate_constant.f90 
                DO l = 1,nelements ! loop over # elements

                   bm_new(i,l) = delta_ind(i,j) * bm_sapl(j,k,l) 
                   IF (veget_cov_max_old(i,j) .GT. min_stomate) THEN

                      ! in the case that bm_new is overestimated compared with biomass?
                      IF ((bm_new(i,l)/delta_veg(j)) > biomass(i,j,k,l)) THEN
                         bm_new(i,l) = biomass(i,j,k,l)*delta_veg(j)
                      ENDIF
                   ENDIF
                   biomass(i,j,k,l) = ( biomass(i,j,k,l) * veget_cov_max_old(i,j) + bm_new(i,l) ) / veget_cov_max_new(i,j)
                   co2_to_bm(i,j) = co2_to_bm(i,j) + (bm_new(i,icarbon)* dt_days) / (one_year * veget_cov_max_new(i,j))
                END DO ! loop over # elements
             ENDDO ! loop over # carbon stock components

             !! 3.3.1.3 Calculation of dilution in litter, soil carbon, and  input of litter
             !!        In this 'IF statement', dilu_* is zero. Formulas for litter and soil carbon
             !!         could be shortend?? Are the following formulas correct?

             ! Litter
             litter(i,:,j,:,:)=(litter(i,:,j,:,:) * veget_cov_max_old(i,j) + &
                  dilu_lit(i,:,:,:) * delta_veg(j)) / veget_cov_max_new(i,j)
            
             ! Soil carbon
             carbon(i,:,j)=(carbon(i,:,j) * veget_cov_max_old(i,j) + dilu_soil_carbon(i,:) * delta_veg(j)) / veget_cov_max_new(i,j)

             DO l = 1,nelements

                ! Litter input
                bm_to_litter(i,j,isapbelow,l) = (bm_to_litter(i,j,isapbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,isapbelow,l)*delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,iheartbelow,l) = (bm_to_litter(i,j,iheartbelow,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iheartbelow,l) *delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,iroot,l) = (bm_to_litter(i,j,iroot,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,iroot,l)*delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ifruit,l) = (bm_to_litter(i,j,ifruit,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ifruit,l)*delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,icarbres,l) = (bm_to_litter(i,j,icarbres,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,icarbres,l)   *delta_veg(j)) / veget_cov_max_new(i,j)
                bm_to_litter(i,j,ileaf,l) = (bm_to_litter(i,j,ileaf,l) * veget_cov_max_old(i,j) + &
                     biomass_loss(i,ileaf,l)*delta_veg(j)) / veget_cov_max_new(i,j)

             END DO

             age(i,j)=age(i,j)*veget_cov_max_old(i,j)/veget_cov_max_new(i,j)
             
          !! 3.3.2 The case that vegetation coverage of PFTj is no change or decreases
          ELSE 
 
             !! 3.3.2.1 Biomass export
             ! coeff_lcchange_*:  Coeff of biomass export for the year, decade, and century
             above = biomass(i,j,isapabove,icarbon) + biomass(i,j,iheartabove,icarbon)
             convflux(i)  = convflux(i)  - ( coeff_lcchange_1(j) * above * delta_veg(j) ) 
             convfluxpft(i,j)= convfluxpft(i,j) - (coeff_lcchange_1(j) * above * delta_veg(j) )
             prod10(i,0)  = prod10(i,0)  - ( coeff_lcchange_10(j) * above * delta_veg(j) )
             prod100(i,0) = prod100(i,0) - ( coeff_lcchange_100(j) * above * delta_veg(j) )

             fDeforestToProduct(i,j) = - above * delta_veg(j)
             fLulccResidue(i,j) = - ( biomass(i,j,isapbelow,icarbon) &
                  + biomass(i,j,iheartbelow,icarbon) &
                  + biomass(i,j,iroot,icarbon) &
                  + biomass(i,j,ifruit,icarbon) &
                  + biomass(i,j,icarbres,icarbon) &
                  + biomass(i,j,ileaf,icarbon) ) * delta_veg(j)
             !! 3.3.2.2 Total reduction
             !! If the vegetation is to small, it has been set to 0.
             IF ( veget_cov_max_new(i,j) .LT. min_stomate ) THEN 
                
                ind(i,j) = zero
                biomass(i,j,:,:) = zero
                PFTpresent(i,j) = .FALSE.
                senescence(i,j) = .FALSE.
                age(i,j) = zero
                when_growthinit(i,j) = undef
                everywhere(i,j) = zero
                carbon(i,:,j) = zero
                litter(i,:,j,:,:) = zero
                bm_to_litter(i,j,:,:) = zero
                turnover_daily(i,j,:,:) = zero
                
             ENDIF
 
          ENDIF ! End if PFT's coverage reduction
          
       ENDDO ! Loop over # PFTs
       
       !! 3.4 update 10 year-turnover pool content following flux emission
       !!     (linear decay (10%) of the initial carbon input)
       DO  l = 0, 8
          m = 10 - l
          cflux_prod10(i) =  cflux_prod10(i) + flux10(i,m)
          prod10(i,m)     =  prod10(i,m-1)   - flux10(i,m-1)
          flux10(i,m)     =  flux10(i,m-1)
          
          IF (prod10(i,m) .LT. 1.0) prod10(i,m) = zero
       ENDDO
       
       cflux_prod10(i) = cflux_prod10(i) + flux10(i,1) 
       flux10(i,1)     = 0.1 * prod10(i,0)
       prod10(i,1)     = prod10(i,0)
       
       !! 3.5 update 100 year-turnover pool content following flux emission\n
       DO   l = 0, 98
          m = 100 - l
          cflux_prod100(i)  =  cflux_prod100(i) + flux100(i,m)
          prod100(i,m)      =  prod100(i,m-1)   - flux100(i,m-1)
          flux100(i,m)      =  flux100(i,m-1)
          
          IF (prod100(i,m).LT.1.0) prod100(i,m) = zero
       ENDDO
       
       cflux_prod100(i)  = cflux_prod100(i) + flux100(i,1) 
       flux100(i,1)      = 0.01 * prod100(i,0)
       prod100(i,1)      = prod100(i,0)
       prod10(i,0)        = zero
       prod100(i,0)       = zero 

       
    ENDDO ! Loop over # pixels - domain size
    
  !! 4. history
    convflux        = convflux/one_year*dt_days
    convfluxpft     = convfluxpft/one_year*dt_days
    fDeforestToProduct= fDeforestToProduct/one_year*dt_days
    fLulccResidue   = fLulccResidue/one_year*dt_days
    cflux_prod10    = cflux_prod10/one_year*dt_days
    cflux_prod100   = cflux_prod100/one_year*dt_days

   
    IF (printlev>=4) WRITE(numout,*) 'Leaving lcchange_main'
    
  END SUBROUTINE lcchange_main
  
END MODULE stomate_lcchange
