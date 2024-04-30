! =================================================================================================================================
! MODULE       : pft_parameters
!
! CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      : IPSL (2011)
! This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This module initializes all the pft parameters in function of the
!!              number of vegetation types and of the values chosen by the user.
!!
!!\n DESCRIPTION:  This module allocates and initializes the pft parameters in function of the number of pfts
!!                 and the values of the parameters. \n
!!                 The number of PFTs is read in control.f90 (subroutine control_initialize). \n
!!                 Then we can initialize the parameters. \n
!!                 This module is the result of the merge of constantes_co2, constantes_veg, stomate_constants.\n
!!
!! RECENT CHANGE(S): Josefine Ghattas 2013 : The declaration part has been extracted and moved to module pft_parameters_var
!!
!! REFERENCE(S)	: None
!!
!! SVN          :
!! $HeadURL: $
!! $Date: 2019-12-16 12:11:50 +0100 (Mon, 16 Dec 2019) $
!! $Revision: 6393 $
!! \n
!_ ================================================================================================================================

MODULE pft_parameters

  USE pft_parameters_var
  USE vertical_soil_var
  USE constantes_mtc
  USE constantes
  USE ioipsl
  USE ioipsl_para 
  USE defprec

  IMPLICIT NONE

CONTAINS
  

!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_main
!!
!>\BRIEF          This subroutine initializes all the pft parameters in function of the
!! number of vegetation types chosen by the user.
!!
!! DESCRIPTION  : This subroutine is called after the reading of the number of PFTS and the options 
!!                activated by the user in the configuration files. \n
!!                The allocation is done just before reading the correspondence table  between PFTs and MTCs
!!                defined by the user in the configuration file.\n
!!                With the correspondence table, the subroutine can initialize the pft parameters in function
!!                of the flags activated (ok_sechiba, ok_stomate, routing,...) in order to
!!                optimize the memory allocation. \n
!!                If the number of PFTs and pft_to_mtc are not found, the standard configuration will be used
!!                (13 PFTs, PFT = MTC). \n 
!!                Some restrictions : the pft 1 can only be the bare soil and it is unique. \n
!!                Algorithm : Build new PFT from 13 generic-PFT or meta-classes.
!!                1. Read the number of PFTs in "run.def". If nothing is found, it is assumed that the user intend to use 
!!                   the standard of PFTs (13).
!!                2. Read the index vector in "run.def". The index vector associates one PFT to one meta-classe (or generic PFT).
!!                   When the association is done, the PFT defined by the user inherited the default values from the meta classe.
!!                   If nothing is found, it is assumed to use the standard index vector (PFT = MTC).
!!                3. Check consistency
!!                4. Memory allocation and initialization.
!!                5. The parameters are read in the configuration file in config_initialize (control module).
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_main()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variables  

    INTEGER(i_std) :: j                             !! Index (unitless)

    !_ ================================================================================================================================ 

    !
    ! PFT global
    !

    IF(l_first_pft_parameters) THEN

       !! 1. First time step
       IF(printlev>=3) THEN
          WRITE(numout,*) 'l_first_pft_parameters :we read the parameters from the def files'
       ENDIF

       !! 2. Memory allocation for the pfts-parameters
       CALL pft_parameters_alloc()

       !! 3. Correspondance table 

       !! 3.1 Initialisation of the correspondance table
       !! Initialisation of the correspondance table
       IF (nvm == nvmc) THEN
          pft_to_mtc = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 /)
       ELSE
          pft_to_mtc(:) = undef_int
       ENDIF !(nvm  == nvmc)

       !! 3.2 Reading of the conrrespondance table in the .def file
       !
       !Config Key   = PFT_TO_MTC
       !Config Desc  = correspondance array linking a PFT to MTC
       !Config if    = OK_SECHIBA or OK_STOMATE
       !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('PFT_TO_MTC',pft_to_mtc)

       !! 3.3 If the user want to use the standard configuration, he needn't to fill the correspondance array
       !!     If the configuration is wrong, send a error message to the user.
       IF(nvm /= nvmc ) THEN
          !
          IF(pft_to_mtc(1) == undef_int) THEN
             STOP ' The array PFT_TO_MTC is empty : we stop'
          ENDIF !(pft_to_mtc(1) == undef_int)
          !
       ENDIF !(nvm /= nvmc )

       !! 3.4 Some error messages

       !! 3.4.1 What happened if pft_to_mtc(j) > nvmc or pft_to_mtc(j) <=0 (if the mtc doesn't exist)?
       DO j = 1, nvm ! Loop over # PFTs  
          !
          IF( (pft_to_mtc(j) > nvmc) .OR. (pft_to_mtc(j) <= 0) ) THEN
             WRITE(numout,*) 'the metaclass chosen does not exist'
             STOP 'we stop reading pft_to_mtc'
          ENDIF !( (pft_to_mtc(j) > nvmc) .OR. (pft_to_mtc(j) <= 0) )
          !
       ENDDO  ! Loop over # PFTs  


       !! 3.4.2 Check if pft_to_mtc(1) = 1 
       IF(pft_to_mtc(1) /= 1) THEN
          !
          WRITE(numout,*) 'the first pft has to be the bare soil'
          STOP 'we stop reading next values of pft_to_mtc'
          !
       ELSE
          !
          DO j = 2,nvm ! Loop over # PFTs different from bare soil
             !
             IF(pft_to_mtc(j) == 1) THEN
                WRITE(numout,*) 'only pft_to_mtc(1) has to be the bare soil'
                STOP 'we stop reading pft_to_mtc'
             ENDIF ! (pft_to_mtc(j) == 1)
             !
          ENDDO ! Loop over # PFTs different from bare soil
          !
       ENDIF !(pft_to_mtc(1) /= 1)


       !! 4.Initialisation of the pfts-parameters
       CALL pft_parameters_init()

       !! 5. Useful data

       !! 5.1 Read the name of the PFTs given by the user
       !
       !Config Key   = PFT_NAME
       !Config Desc  = Name of a PFT
       !Config if    = OK_SECHIBA or OK_STOMATE
       !Config Def   = bare ground, tropical broad-leaved evergreen, tropical broad-leaved raingreen, 
       !Config         temperate needleleaf evergreen, temperate broad-leaved evergreen temperate broad-leaved summergreen,
       !Config         boreal needleleaf evergreen, boreal broad-leaved summergreen, boreal needleleaf summergreen,
       !Config         C3 grass, C4 grass, C3 agriculture, C4 agriculture, urban    
       !Config Help  = the user can name the new PFTs he/she introducing for new species
       !Config Units = [-]
       CALL getin_p('PFT_NAME',pft_name)

       !! 5.2 A useful message to the user: correspondance between the number of the pft
       !! and the name of the associated mtc 
       IF (printlev >=1 ) THEN
          WRITE(numout,*) ''
          DO j = 2,nvm ! Loop over # PFTs
             WRITE(numout,*) 'The PFT',j, 'called ', trim(PFT_name(j)),' corresponds to the MTC : ',trim(MTC_name(pft_to_mtc(j)))
          END DO
          WRITE(numout,*) ''
       END IF


       !! 6. End message
       IF (printlev>=3) WRITE(numout,*) 'pft_parameters_done'

       !! 8. Reset flag
       l_first_pft_parameters = .FALSE.

    ELSE 

       RETURN

    ENDIF !(l_first_pft_parameters)

  END SUBROUTINE pft_parameters_main


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_init 
!!
!>\BRIEF          This subroutine initializes all the pft parameters by the default values
!! of the corresponding metaclasse. 
!!
!! DESCRIPTION  : This subroutine is called after the reading of the number of PFTS and the correspondence
!!                table defined by the user in the configuration files. \n
!!                With the correspondence table, the subroutine can search the default values for the parameter
!!                even if the PFTs are classified in a random order (except bare soil). \n
!!                With the correspondence table, the subroutine can initialize the pft parameters in function
!!                of the flags activated (ok_sechiba, ok_stomate, routing,...).\n
!!
!! RECENT CHANGE(S): Didier Solyga : Simplified PFT loops : use vector notation. 
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_init()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables

    !! 0.4 Local variables

    INTEGER(i_std)                :: jv            !! Index (unitless) 
    !_ ================================================================================================================================ 

    !
    ! 1. Correspondance between the PFTs values and thes MTCs values 
    !


    ! 1.1 For parameters used anytime

    PFT_name(:) = MTC_name(pft_to_mtc(:))
    !
    ! Vegetation structure 
    !
    veget_ori_fixed_test_1(:) = veget_ori_fixed_mtc(pft_to_mtc(:))
    llaimax(:) = llaimax_mtc(pft_to_mtc(:))
    llaimin(:) = llaimin_mtc(pft_to_mtc(:))
    height_presc(:) = height_presc_mtc(pft_to_mtc(:))
    z0_over_height(:) = z0_over_height_mtc(pft_to_mtc(:))
    ratio_z0m_z0h(:) = ratio_z0m_z0h_mtc(pft_to_mtc(:))
    type_of_lai(:) = type_of_lai_mtc(pft_to_mtc(:))
    natural(:) = natural_mtc(pft_to_mtc(:))
    !
    ! Water - sechiba
    !
    IF (zmaxh == 2.0) THEN
       IF (printlev>=2) WRITE(numout,*)'Initialize humcst using reference values for 2m soil depth'
       humcste(:) = humcste_ref2m(pft_to_mtc(:))  ! values for 2m soil depth
    ELSE IF (zmaxh == 4.0) THEN
       IF (printlev>=2) WRITE(numout,*)'Initialize humcst using reference values for 4m soil depth'
       humcste(:) = humcste_ref4m(pft_to_mtc(:))  ! values for 4m soil depth 
    ELSE
       IF (printlev>=2) WRITE(numout,*)'Note that humcste is initialized with values for 2m soil depth bur zmaxh=', zmaxh
       humcste(:) = humcste_ref2m(pft_to_mtc(:))  ! values for 2m soil depth
    END IF
    !
    ! Soil - vegetation
    !
    pref_soil_veg(:) = pref_soil_veg_mtc(pft_to_mtc(:))
    !
    ! Photosynthesis
    !
    is_c4(:) = is_c4_mtc(pft_to_mtc(:))
    vcmax_fix(:) = vcmax_fix_mtc(pft_to_mtc(:))
    downregulation_co2_coeff(:) = downregulation_co2_coeff_mtc(pft_to_mtc(:))
    downregulation_co2_coeff_new(:) = downregulation_co2_coeff_new_mtc(pft_to_mtc(:))
    E_KmC(:)      = E_KmC_mtc(pft_to_mtc(:))
    E_KmO(:)      = E_KmO_mtc(pft_to_mtc(:))
    E_Sco(:)      = E_Sco_mtc(pft_to_mtc(:))
    E_gamma_star(:) = E_gamma_star_mtc(pft_to_mtc(:))
    E_Vcmax(:)    = E_Vcmax_mtc(pft_to_mtc(:))
    E_Jmax(:)     = E_Jmax_mtc(pft_to_mtc(:))
    aSV(:)        = aSV_mtc(pft_to_mtc(:))
    bSV(:)        = bSV_mtc(pft_to_mtc(:))
    tphoto_min(:) = tphoto_min_mtc(pft_to_mtc(:))
    tphoto_max(:) = tphoto_max_mtc(pft_to_mtc(:))
    aSJ(:)        = aSJ_mtc(pft_to_mtc(:))
    bSJ(:)        = bSJ_mtc(pft_to_mtc(:))
    D_Vcmax(:)     = D_Vcmax_mtc(pft_to_mtc(:))
    D_Jmax(:)     = D_Jmax_mtc(pft_to_mtc(:))
    E_gm(:)       = E_gm_mtc(pft_to_mtc(:)) 
    S_gm(:)       = S_gm_mtc(pft_to_mtc(:)) 
    D_gm(:)       = D_gm_mtc(pft_to_mtc(:)) 
    E_Rd(:)       = E_Rd_mtc(pft_to_mtc(:))
    Vcmax25(:)    = Vcmax25_mtc(pft_to_mtc(:))
    arJV(:)       = arJV_mtc(pft_to_mtc(:))
    brJV(:)       = brJV_mtc(pft_to_mtc(:))
    KmC25(:)      = KmC25_mtc(pft_to_mtc(:))
    KmO25(:)      = KmO25_mtc(pft_to_mtc(:))
    Sco25(:)      = Sco25_mtc(pft_to_mtc(:))
    gm25(:)       = gm25_mtc(pft_to_mtc(:)) 
    gamma_star25(:)  = gamma_star25_mtc(pft_to_mtc(:))
    a1(:)         = a1_mtc(pft_to_mtc(:))
    b1(:)         = b1_mtc(pft_to_mtc(:))
    g0(:)         = g0_mtc(pft_to_mtc(:))
    h_protons(:)  = h_protons_mtc(pft_to_mtc(:))
    fpsir(:)      = fpsir_mtc(pft_to_mtc(:))
    fQ(:)         = fQ_mtc(pft_to_mtc(:))     
    fpseudo(:)    = fpseudo_mtc(pft_to_mtc(:))    
    kp(:)         = kp_mtc(pft_to_mtc(:))
    alpha(:)      = alpha_mtc(pft_to_mtc(:))
    gbs(:)        = gbs_mtc(pft_to_mtc(:))
    theta(:)      = theta_mtc(pft_to_mtc(:))        
    alpha_LL(:)   = alpha_LL_mtc(pft_to_mtc(:))
    stress_vcmax(:) = stress_vcmax_mtc(pft_to_mtc(:))
    stress_gs(:)    = stress_gs_mtc(pft_to_mtc(:))
    stress_gm(:)    = stress_gm_mtc(pft_to_mtc(:))
    ext_coeff(:) = ext_coeff_mtc(pft_to_mtc(:))
    ext_coeff_vegetfrac(:) = ext_coeff_vegetfrac_mtc(pft_to_mtc(:))
    !
    !! Define labels from physiologic characteristics 
    !
    leaf_tab(:) = leaf_tab_mtc(pft_to_mtc(:)) 
    pheno_model(:) = pheno_model_mtc(pft_to_mtc(:))   
    !
    is_tree(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) <= 2 ) is_tree(jv) = .TRUE.
    END DO
    !
    is_deciduous(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) /= "none") ) is_deciduous(jv) = .TRUE.
    END DO
    !
    is_evergreen(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) == "none") ) is_evergreen(jv) = .TRUE.
    END DO
    !
    is_needleleaf(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) == 2 ) is_needleleaf(jv) = .TRUE.
    END DO


    ! 1.2 For sechiba parameters

    IF (ok_sechiba) THEN
       !
       ! Vegetation structure - sechiba
       !
       rveg_pft(:) = rveg_mtc(pft_to_mtc(:))
       !
       ! Evapotranspiration -  sechiba
       !
       rstruct_const(:) = rstruct_const_mtc(pft_to_mtc(:))
       kzero(:) = kzero_mtc(pft_to_mtc(:))
       !
       ! Water - sechiba
       !
       wmax_veg(:) = wmax_veg_mtc(pft_to_mtc(:))
       IF ( OFF_LINE_MODE ) THEN
          throughfall_by_pft(:) = 0.
       ELSE
          throughfall_by_pft(:) = throughfall_by_mtc(pft_to_mtc(:))
       ENDIF
       !
       ! Albedo - sechiba
       !
       snowa_aged_vis(:) = snowa_aged_vis_mtc(pft_to_mtc(:))
       snowa_aged_nir(:) = snowa_aged_nir_mtc(pft_to_mtc(:))
       snowa_dec_vis(:) = snowa_dec_vis_mtc(pft_to_mtc(:)) 
       snowa_dec_nir(:) = snowa_dec_nir_mtc(pft_to_mtc(:)) 
       alb_leaf_vis(:) = alb_leaf_vis_mtc(pft_to_mtc(:))  
       alb_leaf_nir(:) = alb_leaf_nir_mtc(pft_to_mtc(:))
       !-
    ENDIF !(ok_sechiba)

    ! 1.3 For BVOC parameters

    IF (ok_bvoc) THEN
       !
       ! Biogenic Volatile Organic Compounds
       !
       em_factor_isoprene(:) = em_factor_isoprene_mtc(pft_to_mtc(:))
       em_factor_monoterpene(:) = em_factor_monoterpene_mtc(pft_to_mtc(:))
       LDF_mono = LDF_mono_mtc 
       LDF_sesq = LDF_sesq_mtc 
       LDF_meth = LDF_meth_mtc 
       LDF_acet = LDF_acet_mtc 

       em_factor_apinene(:) = em_factor_apinene_mtc(pft_to_mtc(:))
       em_factor_bpinene(:) = em_factor_bpinene_mtc(pft_to_mtc(:))
       em_factor_limonene(:) = em_factor_limonene_mtc(pft_to_mtc(:))
       em_factor_myrcene(:) = em_factor_myrcene_mtc(pft_to_mtc(:))
       em_factor_sabinene(:) = em_factor_sabinene_mtc(pft_to_mtc(:))
       em_factor_camphene(:) = em_factor_camphene_mtc(pft_to_mtc(:))
       em_factor_3carene(:) = em_factor_3carene_mtc(pft_to_mtc(:))
       em_factor_tbocimene(:) = em_factor_tbocimene_mtc(pft_to_mtc(:))
       em_factor_othermonot(:) = em_factor_othermonot_mtc(pft_to_mtc(:))
       em_factor_sesquiterp(:) = em_factor_sesquiterp_mtc(pft_to_mtc(:))

       beta_mono = beta_mono_mtc
       beta_sesq = beta_sesq_mtc
       beta_meth = beta_meth_mtc
       beta_acet = beta_acet_mtc
       beta_oxyVOC = beta_oxyVOC_mtc

       em_factor_ORVOC(:) = em_factor_ORVOC_mtc(pft_to_mtc(:)) 
       em_factor_OVOC(:) = em_factor_OVOC_mtc(pft_to_mtc(:))
       em_factor_MBO(:) = em_factor_MBO_mtc(pft_to_mtc(:))
       em_factor_methanol(:) = em_factor_methanol_mtc(pft_to_mtc(:))
       em_factor_acetone(:) = em_factor_acetone_mtc(pft_to_mtc(:)) 
       em_factor_acetal(:) = em_factor_acetal_mtc(pft_to_mtc(:))
       em_factor_formal(:) = em_factor_formal_mtc(pft_to_mtc(:))
       em_factor_acetic(:) = em_factor_acetic_mtc(pft_to_mtc(:))
       em_factor_formic(:) = em_factor_formic_mtc(pft_to_mtc(:))
       em_factor_no_wet(:) = em_factor_no_wet_mtc(pft_to_mtc(:))
       em_factor_no_dry(:) = em_factor_no_dry_mtc(pft_to_mtc(:))
       Larch(:) = Larch_mtc(pft_to_mtc(:)) 
       !-
    ENDIF !(ok_bvoc)

    ! 1.4 For stomate parameters

    IF (ok_stomate) THEN
       !
       ! Vegetation structure - stomate
       !
       sla(:) = sla_mtc(pft_to_mtc(:))
       availability_fact(:) = availability_fact_mtc(pft_to_mtc(:))
       !
       ! Allocation - stomate
       !
       R0(:) = R0_mtc(pft_to_mtc(:)) 
       S0(:) = S0_mtc(pft_to_mtc(:)) 
       !
       ! Respiration - stomate
       !
       frac_growthresp(:) = frac_growthresp_mtc(pft_to_mtc(:))  
       maint_resp_slope_c(:) = maint_resp_slope_c_mtc(pft_to_mtc(:))               
       maint_resp_slope_b(:) = maint_resp_slope_b_mtc(pft_to_mtc(:))
       maint_resp_slope_a(:) = maint_resp_slope_a_mtc(pft_to_mtc(:))
       cm_zero_leaf(:) = cm_zero_leaf_mtc(pft_to_mtc(:))
       cm_zero_sapabove(:) = cm_zero_sapabove_mtc(pft_to_mtc(:))
       cm_zero_sapbelow(:) = cm_zero_sapbelow_mtc(pft_to_mtc(:)) 
       cm_zero_heartabove(:) = cm_zero_heartabove_mtc(pft_to_mtc(:)) 
       cm_zero_heartbelow(:) = cm_zero_heartbelow_mtc(pft_to_mtc(:))
       cm_zero_root(:) = cm_zero_root_mtc(pft_to_mtc(:))
       cm_zero_fruit(:) = cm_zero_fruit_mtc(pft_to_mtc(:))
       cm_zero_carbres(:) = cm_zero_carbres_mtc(pft_to_mtc(:))
       !
       ! Fire - stomate
       !
       flam(:) = flam_mtc(pft_to_mtc(:))
       resist(:) = resist_mtc(pft_to_mtc(:))
       !
       ! Flux - LUC
       !
       coeff_lcchange_1(:) = coeff_lcchange_1_mtc(pft_to_mtc(:))
       coeff_lcchange_10(:) = coeff_lcchange_10_mtc(pft_to_mtc(:))
       coeff_lcchange_100(:) = coeff_lcchange_100_mtc(pft_to_mtc(:))
       !
       ! Phenology
       !
       !
       ! 1. Stomate
       !
       lai_max_to_happy(:) = lai_max_to_happy_mtc(pft_to_mtc(:))  
       lai_max(:) = lai_max_mtc(pft_to_mtc(:))
       pheno_type(:) = pheno_type_mtc(pft_to_mtc(:))
       !
       ! 2. Leaf Onset
       !
       pheno_gdd_crit_c(:) = pheno_gdd_crit_c_mtc(pft_to_mtc(:))
       pheno_gdd_crit_b(:) = pheno_gdd_crit_b_mtc(pft_to_mtc(:))         
       pheno_gdd_crit_a(:) = pheno_gdd_crit_a_mtc(pft_to_mtc(:))
       pheno_moigdd_t_crit(:) = pheno_moigdd_t_crit_mtc(pft_to_mtc(:))
       ngd_crit(:) =  ngd_crit_mtc(pft_to_mtc(:))
       ncdgdd_temp(:) = ncdgdd_temp_mtc(pft_to_mtc(:)) 
       hum_frac(:) = hum_frac_mtc(pft_to_mtc(:))
       hum_min_time(:) = hum_min_time_mtc(pft_to_mtc(:))
       tau_sap(:) = tau_sap_mtc(pft_to_mtc(:))
       tau_leafinit(:) = tau_leafinit_mtc(pft_to_mtc(:))  
       tau_fruit(:) = tau_fruit_mtc(pft_to_mtc(:))
       ecureuil(:) = ecureuil_mtc(pft_to_mtc(:))
       alloc_min(:) = alloc_min_mtc(pft_to_mtc(:))
       alloc_max(:) = alloc_max_mtc(pft_to_mtc(:))
       demi_alloc(:) = demi_alloc_mtc(pft_to_mtc(:))
       leaflife_tab(:) = leaflife_mtc(pft_to_mtc(:))
       !
       ! 3. Senescence
       !
       leaffall(:) = leaffall_mtc(pft_to_mtc(:))
       leafagecrit(:) = leafagecrit_mtc(pft_to_mtc(:))
       senescence_type(:) = senescence_type_mtc(pft_to_mtc(:)) 
       senescence_hum(:) = senescence_hum_mtc(pft_to_mtc(:)) 
       nosenescence_hum(:) = nosenescence_hum_mtc(pft_to_mtc(:)) 
       max_turnover_time(:) = max_turnover_time_mtc(pft_to_mtc(:))
       min_turnover_time(:) = min_turnover_time_mtc(pft_to_mtc(:))
       min_leaf_age_for_senescence(:) = min_leaf_age_for_senescence_mtc(pft_to_mtc(:))
       senescence_temp_c(:) = senescence_temp_c_mtc(pft_to_mtc(:))
       senescence_temp_b(:) = senescence_temp_b_mtc(pft_to_mtc(:))
       senescence_temp_a(:) = senescence_temp_a_mtc(pft_to_mtc(:))
       gdd_senescence(:) = gdd_senescence_mtc(pft_to_mtc(:))
       always_init(:) = always_init_mtc(pft_to_mtc(:))
       !
       ! DGVM
       !
       residence_time(:) = residence_time_mtc(pft_to_mtc(:))
       tmin_crit(:) = tmin_crit_mtc(pft_to_mtc(:))
       tcm_crit(:) = tcm_crit_mtc(pft_to_mtc(:))
       !-
    ENDIF !(ok_stomate)

  END SUBROUTINE pft_parameters_init


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_alloc
!!
!>\BRIEF         This subroutine allocates memory needed for the PFT parameters 
!! in function  of the flags activated.  
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_alloc()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables 

    !! 0.4 Local variables

    LOGICAL :: l_error                             !! Diagnostic boolean for error allocation (true/false) 
    INTEGER :: ier                                 !! Return value for memory allocation (0-N, unitless)

    !_ ================================================================================================================================


    !
    ! 1. Parameters used anytime
    !

    l_error = .FALSE.

    ALLOCATE(pft_to_mtc(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pft_to_mtc. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(PFT_name(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for PFT_name. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(height_presc(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for height_presc. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(z0_over_height(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for z0_over_height. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(ratio_z0m_z0h(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ratio_z0m_z0h. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_tree(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_tree. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(natural(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for natural. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_c4(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_c4. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(humcste(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for humcste. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(downregulation_co2_coeff(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for downregulation_co2_coeff. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(downregulation_co2_coeff_new(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for downregulation_co2_coeff_new. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_KmC(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_KmC. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_KmO(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_KmO. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_Sco(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Sco. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_gamma_star(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_gamma_star. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Vcmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_Jmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Jmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(aSV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for aSV. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(bSV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for bSV. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(tphoto_min(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for tphoto_min. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(tphoto_max(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for tphoto_max. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(aSJ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for aSJ. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(bSJ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for bSJ. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(D_Vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for D_Vcmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(D_Jmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for D_Jmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(E_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for E_gm. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF
    
    ALLOCATE(S_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for S_gm. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF
    
    ALLOCATE(D_gm(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for D_gm. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF
    
    ALLOCATE(E_Rd(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for E_Rd. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(Vcmax25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Vcmax25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(arJV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for arJV. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(brJV(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for brJV. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(KmC25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for KmC25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(KmO25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for KmO25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(Sco25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for Sco25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
    
    ALLOCATE(gm25(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN 
       WRITE(numout,*) ' Memory allocation error for gm25. We stop. We need nvm words = ',nvm 
       STOP 'pft_parameters_alloc' 
    END IF

    ALLOCATE(gamma_star25(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for gamma_star25. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(a1(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for a1. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(b1(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for b1. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(g0(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for g0. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(h_protons(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for h_protons. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(fpsir(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fpsir. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(fQ(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fQ. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(fpseudo(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for fpseudo. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(kp(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for kp. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(alpha(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for alpha. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(gbs(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for gbs. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(theta(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for theta. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(alpha_LL(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for alpha_LL. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(stress_vcmax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_vcmax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
    
    ALLOCATE(stress_gs(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_gs. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF
    
    ALLOCATE(stress_gm(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for stress_gm. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(ext_coeff(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ext_coeff. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(ext_coeff_vegetfrac(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for ext_coeff_vegetfrac. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(veget_ori_fixed_test_1(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for veget_ori_fixed_test_1. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(llaimax(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for llaimax. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(llaimin(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for llaimin. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(type_of_lai(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for type_of_lai. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(vcmax_fix(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for vcmax_fix. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(pref_soil_veg(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pref_soil_veg. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(leaf_tab(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for leaf_tab. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(pheno_model(nvm),stat=ier)
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for pheno_model. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_deciduous(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0) 
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_deciduous. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_evergreen(nvm),stat=ier) 
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_evergreen. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_needleleaf(nvm),stat=ier)  
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_needleleaf. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF

    ALLOCATE(is_tropical(nvm),stat=ier)   
    l_error = l_error .OR. (ier /= 0)
    IF (l_error) THEN
       WRITE(numout,*) ' Memory allocation error for is_tropical. We stop. We need nvm words = ',nvm
       STOP 'pft_parameters_alloc'
    END IF


    !
    ! 2. Parameters used if ok_sechiba only
    !
    IF ( ok_sechiba ) THEN

       l_error = .FALSE.

       ALLOCATE(rstruct_const(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rstruct_const. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(kzero(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for kzero. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(rveg_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for rveg_pft. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(wmax_veg(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for wmax_veg. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(throughfall_by_pft(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for throughfall_by_pft. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(snowa_aged_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_aged_vis. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(snowa_aged_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_aged_nir. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(snowa_dec_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_dec_vis. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(snowa_dec_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for snowa_dec_nir. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(alb_leaf_vis(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alb_leaf_vis. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(alb_leaf_nir(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alb_leaf_nir. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       IF( ok_bvoc ) THEN

          l_error = .FALSE.

          ALLOCATE(em_factor_isoprene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_isoprene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_monoterpene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_monoterpene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_apinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_apinene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_bpinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_bpinene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_limonene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_limonene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_myrcene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_myrcene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_sabinene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_sabinene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_camphene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_camphene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_3carene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_3carene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_tbocimene(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_tbocimene. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_othermonot(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_othermonot. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_sesquiterp(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_sesquiterp. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF


          ALLOCATE(em_factor_ORVOC(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_ORVOC. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_OVOC(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_OVOC. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_MBO(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_MBO. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_methanol(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_methanol. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_acetone(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetone. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_acetal(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetal. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_formal(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_formal. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_acetic(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_acetic. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_formic(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_formic. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_no_wet(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_no_wet. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(em_factor_no_dry(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0)       
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for em_factor_no_dry. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

          ALLOCATE(Larch(nvm),stat=ier)
          l_error = l_error .OR. (ier /= 0) 
          IF (l_error) THEN
             WRITE(numout,*) ' Memory allocation error for Larch. We stop. We need nvm words = ',nvm
             STOP 'pft_parameters_alloc'
          END IF

       ENDIF ! (ok_bvoc) 

    ENDIF !(ok_sechiba)

    !
    ! 3. Parameters used if ok_stomate only
    !
    IF ( ok_stomate ) THEN

       l_error = .FALSE.

       ALLOCATE(sla(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for sla. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(availability_fact(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for availability_fact. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(R0(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for R0. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(S0(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for S0. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(L0(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for L0. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_gdd_crit_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_c. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_gdd_crit_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_b. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_gdd_crit_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit_a. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_gdd_crit(nvm,3),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_gdd_crit. We stop. We need nvm words = ',nvm*3
          STOP 'pft_parameters_alloc'
       END IF
       pheno_gdd_crit(:,:) = zero

       ALLOCATE(pheno_moigdd_t_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_moigdd_t_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(ngd_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ngd_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(ncdgdd_temp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ncdgdd_temp. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(hum_frac(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for hum_frac. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(hum_min_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for hum_min_time. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tau_sap(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tau_sap. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tau_leafinit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tau_leafinit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tau_fruit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tau_fruit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(ecureuil(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for ecureuil. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(alloc_min(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alloc_min. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(alloc_max(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for alloc_max. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(demi_alloc(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(frac_growthresp(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for frac_growthresp. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(maint_resp_slope(nvm,3),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maint_resp_slope. We stop. We need nvm*3 words = ',nvm*3
          STOP 'pft_parameters_alloc'
       END IF
       maint_resp_slope(:,:) = zero

       ALLOCATE(maint_resp_slope_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maint_resp_slope_c. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(maint_resp_slope_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maint_resp_slope_b. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(maint_resp_slope_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maint_resp_slope_a. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_maint_zero(nvm,nparts),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_maint_zero. We stop. We need nvm*nparts words = ',nvm*nparts
          STOP 'pft_parameters_alloc'
       END IF
       coeff_maint_zero(:,:) = zero

       ALLOCATE(cm_zero_leaf(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_leaf. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_sapabove(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_sapabove. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_sapbelow(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_sapbelow. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_heartabove(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_heartabove. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_heartbelow(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_heartbelow. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_root(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_root. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_fruit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_fruit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cm_zero_carbres(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cm_zero_carbres. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(flam(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(resist(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for resist. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_lcchange_1(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_1. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_lcchange_10(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_10. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(coeff_lcchange_100(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for coeff_lcchange_100. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(lai_max_to_happy(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_max_to_happy. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(lai_max(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for lai_max. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(pheno_type(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for pheno_type. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(leaffall(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaffall. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(leafagecrit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leafagecrit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_type(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_hum(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_hum. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(nosenescence_hum(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for nosenescence_hum. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(max_turnover_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for max_turnover_time. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(min_turnover_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_turnover_time. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(min_leaf_age_for_senescence(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for min_leaf_age_for_senescence. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_temp_c(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_c. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_temp_b(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_b. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_temp_a(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp_a. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(senescence_temp(nvm,3),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for senescence_temp. We stop. We need nvm*3 words = ',nvm*3
          STOP 'pft_parameters_alloc'
       END IF
       senescence_temp(:,:) = zero

       ALLOCATE(gdd_senescence(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for gdd_senescence. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(always_init(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for always_init. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(residence_time(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for residence_time. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tmin_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tmin_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(tcm_crit(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for tcm_crit. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(lai_initmin(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for . We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(bm_sapl(nvm,nparts,nelements),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for bm_sapl. We stop. We need nvm*nparts*nelements words = ',& 
               &  nvm*nparts*nelements
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(migrate(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for migrate. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(maxdia(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for maxdia. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(cn_sapl(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for cn_sapl. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(leaf_timecst(nvm),stat=ier)
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaf_timecst. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

       ALLOCATE(leaflife_tab(nvm),stat=ier)   
       l_error = l_error .OR. (ier /= 0)
       IF (l_error) THEN
          WRITE(numout,*) ' Memory allocation error for leaflife_tab. We stop. We need nvm words = ',nvm
          STOP 'pft_parameters_alloc'
       END IF

    ENDIF ! (ok_stomate)

  END SUBROUTINE pft_parameters_alloc

!! ================================================================================================================================
!! SUBROUTINE   : config_pft_parameters 
!!
!>\BRIEF          This subroutine will read the imposed values for the global pft
!! parameters (sechiba + stomate). It is not called if IMPOSE_PARAM is set to NO.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_pft_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variable

    INTEGER(i_std) :: jv                   !! Index (untiless)

    !_ ================================================================================================================================ 


    !
    ! Vegetation structure
    !

    !Config Key   = LEAF_TAB
    !Config Desc  = leaf type : 1=broad leaved tree, 2=needle leaved tree, 3=grass 4=bare ground
    !Config if    = OK_STOMATE
    !Config Def   = 4, 1, 1, 2, 1, 1, 2, 1, 2, 3, 3, 3, 3, 4 
    !Config Help  = 
    !Config Units = [-] 
    CALL getin_p('LEAF_TAB',leaf_tab)

    !Config Key   = PHENO_MODEL
    !Config Desc  = which phenology model is used? (tabulated) 
    !Config if    = OK_STOMATE
    !Config Def   = none, none, moi, none, none, ncdgdd, none, ncdgdd, ngd, moigdd, moigdd, moigdd, moigdd, none
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('PHENO_MODEL',pheno_model)

    !! Redefine the values for is_tree, is_deciduous, is_needleleaf, is_evergreen if values have been modified
    !! in run.def

    is_tree(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) <= 2 ) is_tree(jv) = .TRUE.
    END DO
    !
    is_deciduous(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) /= "none") ) is_deciduous(jv) = .TRUE.
    END DO
    !
    is_evergreen(:) = .FALSE.
    DO jv = 1,nvm
       IF ( is_tree(jv) .AND. (pheno_model(jv) == "none") ) is_evergreen(jv) = .TRUE.
    END DO
    !
    is_needleleaf(:) = .FALSE.
    DO jv = 1,nvm
       IF ( leaf_tab(jv) == 2 ) is_needleleaf(jv) = .TRUE.
    END DO


    !Config Key   = SECHIBA_LAI
    !Config Desc  = laimax for maximum lai(see also type of lai interpolation)
    !Config if    = OK_SECHIBA or IMPOSE_VEG
    !Config Def   = 0., 8., 8., 4., 4.5, 4.5, 4., 4.5, 4., 2., 2., 2., 2., 0.
    !Config Help  = Maximum values of lai used for interpolation of the lai map
    !Config Units = [m^2/m^2]
    CALL getin_p('SECHIBA_LAI',llaimax)

    !Config Key   = LLAIMIN
    !Config Desc  = laimin for minimum lai(see also type of lai interpolation)
    !Config if    = OK_SECHIBA or IMPOSE_VEG
    !Config Def   = 0., 8., 0., 4., 4.5, 0., 4., 0., 0., 0., 0., 0., 0., 0.
    !Config Help  = Minimum values of lai used for interpolation of the lai map
    !Config Units = [m^2/m^2]
    CALL getin_p('LLAIMIN',llaimin)

    !Config Key   = SLOWPROC_HEIGHT
    !Config Desc  = prescribed height of vegetation 
    !Config if    = OK_SECHIBA
    !Config Def   = 0., 30., 30., 20., 20., 20., 15., 15., 15., .5, .6, 1., 1., 0.
    !Config Help  =
    !Config Units = [m] 
    CALL getin_p('SLOWPROC_HEIGHT',height_presc)

    !Config Key   = Z0_OVER_HEIGHT
    !Config Desc  = factor to calculate roughness height from height of canopy 
    !Config if    = OK_SECHIBA
    !Config Def   = 0., 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('Z0_OVER_HEIGHT',z0_over_height)

    !
    !Config Key   = RATIO_Z0M_Z0H
    !Config Desc  = Ratio between z0m and z0h
    !Config Def   = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 
    !Config if    = OK_SECHIBA
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('RATIO_Z0M_Z0H',ratio_z0m_z0h)


    !Config Key   = TYPE_OF_LAI
    !Config Desc  = Type of behaviour of the LAI evolution algorithm 
    !Config if    = OK_SECHIBA
    !Config Def   = inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter, inter
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TYPE_OF_LAI',type_of_lai)

    !Config Key   = NATURAL
    !Config Desc  = natural? 
    !Config if    = OK_SECHIBA, OK_STOMATE
    !Config Def   = y, y, y, y, y, y, y, y, y, y, y, n, n, y 
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('NATURAL',natural)


    !
    ! Photosynthesis
    !

    !Config Key   = IS_C4
    !Config Desc  = flag for C4 vegetation types
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = n, n, n, n, n, n, n, n, n, n, n, y, n, y, n
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('IS_C4',is_c4)

    !Config Key   = VCMAX_FIX
    !Config Desc  = values used for vcmax when STOMATE is not activated
    !Config if    = OK_SECHIBA and NOT(OK_STOMATE)
    !Config Def   = 0., 40., 50., 30., 35., 40.,30., 40., 35., 60., 60., 70., 70., 0.
    !Config Help  =
    !Config Units = [micromol/m^2/s] 
    CALL getin_p('VCMAX_FIX',vcmax_fix)

    !Config Key   = DOWNREG_CO2
    !Config Desc  = coefficient for CO2 downregulation (unitless)
    !Config if    = OK_CO2 and DOWNREGULATION_CO2
    !Config Def   = 0., 0.38, 0.38, 0.28, 0.28, 0.28, 0.22, 0.22, 0.22, 0.26, 0.26, 0.26, 0.26, 0.
    !Config Help  = 
    !Config Units = [-]
    CALL getin_p('DOWNREG_CO2',downregulation_co2_coeff)

    !Config Key   = DOWNREG_CO2_NEW
    !Config Desc  = coefficient for CO2 downregulation (unitless)
    !Config if    = OK_CO2 and DOWNREGULATION_CO2_NEW
    !Config Def   = 0., 0.35, 0.35, 0.26, 0.26, 0.26, 0.20, 0.20, 0.20, 0.24, 0.03, 0.24, 0.03, 0.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('DOWNREG_CO2_NEW',downregulation_co2_coeff_new)

    !Config Key   = E_KmC
    !Config Desc  = Energy of activation for KmC
    !Config if    = 
    !Config Def   = undef,  79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., 79430., undef
    !Config Help  = See Medlyn et al. (2002) 
    !Config Units = [J mol-1]
    CALL getin_p('E_KMC',E_KmC)

    !Config Key   = E_KmO
    !Config Desc  = Energy of activation for KmO
    !Config if    = 
    !Config Def   = undef, 36380.,  36380.,  36380.,  36380.,  36380., 36380., 36380., 36380., 36380., 36380., 36380., 36380., undef
    !Config Help  = See Medlyn et al. (2002) 
    !Config Units = [J mol-1]
    CALL getin_p('E_KMO',E_KmO)

    !Config Key   = E_Sco
    !Config Desc  = Energy of activation for Sco
    !Config if    = 
    !Config Def   = undef, -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., -24460., undef
    !Config Help  = See Table 2 of Yin et al. (2009) - Value for C4 plants is not mentioned - We use C3 for all plants
    !Config Units = [J mol-1]
    CALL getin_p('E_SCO',E_Sco)
    
    !Config Key   = E_gamma_star
    !Config Desc  = Energy of activation for gamma_star
    !Config if    = 
    !Config Def   = undef, 37830.,  37830.,  37830.,  37830.,  37830., 37830., 37830., 37830., 37830., 37830., 37830., 37830., undef
    !Config Help  = See Medlyn et al. (2002) from Bernacchi al. (2001) 
    !Config Units = [J mol-1]
    CALL getin_p('E_GAMMA_STAR',E_gamma_star)

    !Config Key   = E_Vcmax
    !Config Desc  = Energy of activation for Vcmax
    !Config if    = 
    !Config Def   = undef, 71513., 71513., 71513., 71513., 71513., 71513., 71513., 71513., 71513., 67300., 71513., 67300., undef
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Kattge & Knorr (2007) for C3 plants (table 3)
    !Config Units = [J mol-1]
    CALL getin_p('E_VCMAX',E_Vcmax)

    !Config Key   = E_Jmax
    !Config Desc  = Energy of activation for Jmax
    !Config if    = 
    !Config Def   = undef, 49884., 49884., 49884., 49884., 49884., 49884., 49884., 49884., 49884., 77900., 49884., 77900., undef 
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Kattge & Knorr (2007) for C3 plants (table 3)
    !Config Units = [J mol-1]
    CALL getin_p('E_JMAX',E_Jmax)

    !Config Key   = aSV
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax
    !Config if    = 
    !Config Def   = undef, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 668.39, 641.64, 668.39, 641.64, undef 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation and that at for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
    !Config Units = [J K-1 mol-1]
    CALL getin_p('ASV',aSV)

    !Config Key   = bSV
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Entropy term for Vcmax
    !Config if    = 
    !Config Def   = undef, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, -1.07, 0., -1.07, 0., undef 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation
    !Config Units = [J K-1 mol-1 °C-1]
    CALL getin_p('BSV',bSV)

    !Config Key   = TPHOTO_MIN
    !Config Desc  = minimum photosynthesis temperature (deg C)
    !Config if    = OK_STOMATE
    !Config Def   = undef,  -4., -4., -4., -4.,-4.,-4., -4., -4., -4., -4., -4., -4., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TPHOTO_MIN',tphoto_min)

    !Config Key   = TPHOTO_MAX
    !Config Desc  = maximum photosynthesis temperature (deg C)
    !Config if    = OK_STOMATE
    !Config Def   = undef, 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., 55., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('TPHOTO_MAX',tphoto_max)

    !Config Key   = aSJ
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Entropy term for Jmax
    !Config if    = 
    !Config Def   = undef, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 659.70, 630., 659.70, 630., undef 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - and Table 2 of Yin et al. (2009) for C4 plants
    !Config Units = [J K-1 mol-1]
    CALL getin_p('ASJ',aSJ)

    !Config Key   = bSJ
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Entropy term for Jmax
    !Config if    = 
    !Config Def   = undef, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, -0.75, 0., -0.75, 0., undef 
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation
    !Config Units = [J K-1 mol-1 °C-1]
    CALL getin_p('BSJ',bSJ)

    !Config Key   = D_Vcmax
    !Config Desc  = Energy of deactivation for Vcmax
    !Config if    = 
    !Config Def   = undef, 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 192000., 200000., 192000., undef
    !Config Help  = Medlyn et al. (2002) also uses 200000. for C3 plants (same value than D_Jmax). 'Consequently', we use the value of D_Jmax for C4 plants.
    !Config Units = [J mol-1]
    CALL getin_p('D_VCMAX',D_Vcmax)

    !Config Key   = D_Jmax
    !Config Desc  = Energy of deactivation for Jmax
    !Config if    = 
    !Config Def   = undef, 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 200000., 192000., 200000., 192000., undef
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [J mol-1]
    CALL getin_p('D_JMAX',D_Jmax)
    
    !Config Key   = E_gm 
    !Config Desc  = Energy of activation for gm 
    !Config if    =  
    !Config Def   = undef, 49600., 49600., 49600., 49600., 49600., 49600., 49600., 49600., 49600., undef, 49600., undef, undef
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J mol-1] 
    CALL getin_p('E_GM',E_gm) 
    
    !Config Key   = S_gm 
    !Config Desc  = Entropy term for gm 
    !Config if    =  
    !Config Def   = undef, 1400., 1400., 1400., 1400., 1400., 1400., 1400., 1400., 1400., undef, 1400., undef, undef 
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J K-1 mol-1] 
    CALL getin_p('S_GM',S_gm) 
    
    !Config Key   = D_gm 
    !Config Desc  = Energy of deactivation for gm 
    !Config if    =  
    !Config Def   = undef, 437400., 437400., 437400., 437400., 437400., 437400., 437400., 437400., 437400., undef, 437400., undef, undef 
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [J mol-1] 
    CALL getin_p('D_GM',D_gm) 
    
    !Config Key   = E_Rd
    !Config Desc  = Energy of activation for Rd
    !Config if    = 
    !Config Def   = undef, 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., 46390., undef
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [J mol-1]
    CALL getin_p('E_RD',E_Rd)

    !Config Key   = VCMAX25
    !Config Desc  = Maximum rate of Rubisco activity-limited carboxylation at 25°C
    !Config if    = OK_STOMATE
    !Config Def   = undef, 45.0, 45.0, 35.0, 40.0, 50.0, 45.0, 35.0, 35.0, 50.0, 50.0, 60.0, 60.0, undef
    !Config Help  =
    !Config Units = [micromol/m^2/s]
    CALL getin_p('VCMAX25',Vcmax25)

    !Config Key   = ARJV
    !Config Desc  = a coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio 
    !Config if    = OK_STOMATE
    !Config Def   = undef, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 2.59, 1.715, 2.59, 1.715, undef
    !Config Help  = See Table 3 of Kattge & Knorr (2007) - For C4 plants, we assume that there is no acclimation and that for a temperature of 25°C, aSV is the same for both C4 and C3 plants (no strong jusitification - need further parametrization)
    !Config Units = [mu mol e- (mu mol CO2)-1]
    CALL getin_p('ARJV',arJV)

    !Config Key   = BRJV
    !Config Desc  = b coefficient of the linear regression (a+bT) defining the Jmax25/Vcmax25 ratio 
    !Config if    = OK_STOMATE
    !Config Def   = undef, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, -0.035, 0., -0.035, 0., undef
    !Config Help  = See Table 3 of Kattge & Knorr (2007) -  We assume No acclimation term for C4 plants
    !Config Units = [(mu mol e- (mu mol CO2)-1) (°C)-1]
    CALL getin_p('BRJV',brJV)

    !Config Key   = KmC25
    !Config Desc  = Michaelis–Menten constant of Rubisco for CO2 at 25°C
    !Config if    = 
    !Config Def   = undef, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 404.9, 650., 404.9, 650., undef
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3 plants
    !Config Units = [ubar]
    CALL getin_p('KMC25',KmC25)

    !Config Key   = KmO25
    !Config Desc  = Michaelis–Menten constant of Rubisco for O2 at 25°C
    !Config if    = 
    !Config Def   = undef, 278400., 278400., 278400., 278400., 278400., 278400., 278400., 278400., 278400., 450000., 278400., 450000., undef
    !Config Help  = See Table 2 of Yin et al. (2009) for C4 plants and Medlyn et al. (2002) for C3 plants
    !Config Units = [ubar]
    CALL getin_p('KMO25',KmO25)

    !Config Key   = Sco25
    !Config Desc  = Relative CO2 /O2 specificity factor for Rubisco at 25Â°C
    !Config if    = 
    !Config Def   = undef, 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2800., 2590., 2800., 2590., undef
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [bar bar-1]
    CALL getin_p('SCO25',Sco25)
    
    !Config Key   = gm25 
    !Config Desc  = Mesophyll diffusion conductance at 25ÃÂ°C 
    !Config if    =  
    !Config Def   = undef, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, undef, 0.4, undef, undef 
    !Config Help  = See legend of Figure 6 of Yin et al. (2009) and review by Flexas et al. (2008) - gm is not used for C4 plants 
    !Config Units = [mol m-2 s-1 bar-1] 
    CALL getin_p('GM25',gm25) 
    
    !Config Key   = gamma_star25
    !Config Desc  = Ci-based CO2 compensation point in the absence of Rd at 25°C (ubar)
    !Config if    = 
    !Config Def   = undef, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, 42.75, undef
    !Config Help  = See Medlyn et al. (2002) for C3 plants - For C4 plants, we use the same value (probably uncorrect)
    !Config Units = [ubar]
    CALL getin_p('gamma_star25',gamma_star25)

    !Config Key   = a1
    !Config Desc  = Empirical factor involved in the calculation of fvpd
    !Config if    = 
    !Config Def   = undef, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.72, 0.85, 0.72, undef
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('A1',a1)

    !Config Key   = b1
    !Config Desc  = Empirical factor involved in the calculation of fvpd
    !Config if    = 
    !Config Def   = undef, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.20, 0.14, 0.20, undef
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('B1',b1)

    !Config Key   = g0
    !Config Desc  = Residual stomatal conductance when irradiance approaches zero 
    !Config if    = 
    !Config Def   = undef, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.00625, 0.01875, 0.00625, 0.01875, undef 
    !Config Help  = Value from ORCHIDEE - No other reference.
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('G0',g0)

    !Config Key   = h_protons
    !Config Desc  = Number of protons required to produce one ATP
    !Config if    = 
    !Config Def   = undef, 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., undef 
    !Config Help  = See Table 2 of Yin et al. (2009) - h parameter
    !Config Units = [mol mol-1]
    CALL getin_p('H_PROTONS',h_protons)

    !Config Key   = fpsir
    !Config Desc  = Fraction of PSII e− transport rate partitioned to the C4 cycle
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.4, undef, 0.4, undef 
    !Config Help  = See Table 2 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('FPSIR',fpsir)

    !Config Key   = fQ
    !Config Desc  = Fraction of electrons at reduced plastoquinone that follow the Q-cycle
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 1., undef, 1., undef
    !Config Help  = See Table 2 of Yin et al. (2009) - Values for C3 plants are not used
    !Config Units = [-]
    CALL getin_p('FQ',fQ)

    !Config Key   = fpseudo
    !Config Desc  = Fraction of electrons at PSI that follow pseudocyclic transport 
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.1, undef, 0.1, undef
    !Config Help  = See Table 2 of Yin et al. (2009) - Values for C3 plants are not used
    !Config Units = [-]
    CALL getin_p('FPSEUDO',fpseudo)

    !Config Key   = kp
    !Config Desc  = Initial carboxylation efficiency of the PEP carboxylase
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.7, undef, 0.7, undef
    !Config Help  = See Table 2 of Yin et al. (2009) 
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('KP',kp)

    !Config Key   = alpha
    !Config Desc  = Fraction of PSII activity in the bundle sheath
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.1, undef, 0.1, undef
    !Config Help  = See legend of Figure 6 of Yin et al. (2009)
    !Config Units = [-]
    CALL getin_p('ALPHA',alpha)

    !Config Key   = gbs
    !Config Desc  = Bundle-sheath conductance
    !Config if    = 
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.003, undef, 0.003, undef
    !Config Help  = See legend of Figure 6 of Yin et al. (2009)
    !Config Units = [mol m−2 s−1 bar−1]
    CALL getin_p('GBS',gbs)

    !Config Key   = theta
    !Config Desc  = Convexity factor for response of J to irradiance
    !Config if    = 
    !Config Def   = undef, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, undef
    !Config Help  = See Table 2 of Yin et al. (2009)   
    !Config Units = [−]
    CALL getin_p('THETA',theta)

    !Config Key   = alpha_LL
    !Config Desc  = Conversion efficiency of absorbed light into J at strictly limiting light
    !Config if    = 
    !Config Def   = undef, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372, undef
    !Config Help  = See comment from Yin et al. (2009) after eq. 4
    !Config Units = [mol e− (mol photon)−1]
    CALL getin_p('ALPHA_LL',alpha_LL)

    !Config Key   = STRESS_VCMAX
    !Config Desc  = Stress on vcmax
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_VCMAX', stress_vcmax)
    
    !Config Key   = STRESS_GS
    !Config Desc  = Stress on gs
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_GS', stress_gs)
    
    !Config Key   = STRESS_GM
    !Config Desc  = Stress on gm
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('STRESS_GM', stress_gm)

    !Config Key   = EXT_COEFF
    !Config Desc  = extinction coefficient of the Monsi&Seaki relationship (1953)
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('EXT_COEFF',ext_coeff)

    !Config Key   = EXT_COEFF_VEGETFRAC
    !Config Desc  = extinction coefficient used for the calculation of the bare soil fraction 
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('EXT_COEFF_VEGETFRAC',ext_coeff_vegetfrac)

    !
    ! Water-hydrology - sechiba
    !

    !Config Key   = HYDROL_HUMCSTE
    !Config Desc  = Root profile
    !Config Def   = humcste_ref2m or humcste_ref4m depending on zmaxh
    !Config if    = OK_SECHIBA
    !Config Help  = See module constantes_mtc for different default values
    !Config Units = [m]
    CALL getin_p('HYDROL_HUMCSTE',humcste)

    !
    ! Soil - vegetation
    !

    !Config Key   = PREF_SOIL_VEG
    !Config Desc  = The soil tile number for each vegetation
    !Config if    = OK_SECHIBA or OK_STOMATE
    !Config Def   = 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4
    !Config Help  = Gives the number of the soil tile on which we will
    !Config         put each vegetation. This allows to divide the hydrological column
    !Config Units = [-]        
    CALL getin_p('PREF_SOIL_VEG',pref_soil_veg)

  END SUBROUTINE config_pft_parameters


!! ================================================================================================================================
!! SUBROUTINE   : config_sechiba_pft_parameters
!!
!>\BRIEF        This subroutine will read the imposed values for the sechiba pft
!! parameters. It is not called if IMPOSE_PARAM is set to NO. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_sechiba_pft_parameters()

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.1 Input variables

    !! 0.4 Local variable

    !_ ================================================================================================================================ 

    !
    ! Evapotranspiration -  sechiba
    !

    !Config Key   = RSTRUCT_CONST
    !Config Desc  = Structural resistance 
    !Config if    = OK_SECHIBA
    !Config Def   = 0.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0,  2.5,  2.0,  2.0,  2.0, 0.0
    !Config Help  =
    !Config Units = [s/m]
    CALL getin_p('RSTRUCT_CONST',rstruct_const)

    !Config Key   = KZERO
    !Config Desc  = A vegetation dependent constant used in the calculation of the surface resistance.
    !Config if    = OK_SECHIBA
    !Config Def   = 0.0, 12.E-5, 12.E-5, 12.e-5, 12.e-5, 25.e-5, 12.e-5,25.e-5, 25.e-5, 30.e-5, 30.e-5, 30.e-5, 30.e-5, 0.0 
    !Config Help  =
    !Config Units = [kg/m^2/s]
    CALL getin_p('KZERO',kzero)

    !Config Key   = RVEG_PFT
    !Config Desc  = Artificial parameter to increase or decrease canopy resistance.
    !Config if    = OK_SECHIBA
    !Config Def   = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.
    !Config Help  = This parameter is set by PFT.
    !Config Units = [-]
    CALL getin_p('RVEG_PFT',rveg_pft)    

    !
    ! Water-hydrology - sechiba
    !

    !Config Key   = WMAX_VEG
    !Config Desc  = Maximum field capacity for each of the vegetations (Temporary): max quantity of water
    !Config if    = OK_SECHIBA
    !Config Def   = 150., 150., 150., 150., 150., 150., 150.,150., 150., 150., 150., 150., 150., 150.
    !Config Help  =
    !Config Units = [kg/m^3]
    CALL getin_p('WMAX_VEG',wmax_veg)

    !Config Key   = PERCENT_THROUGHFALL_PFT
    !Config Desc  = Percent by PFT of precip that is not intercepted by the canopy. Default value depend on run mode.
    !Config if    = OK_SECHIBA
    !Config Def   = Case offline [0. 0. 0....] else [30. 30. 30.....]
    !Config Help  = During one rainfall event, PERCENT_THROUGHFALL_PFT% of the incident rainfall
    !Config         will get directly to the ground without being intercepted, for each PFT.
    !Config Units = [%]
    CALL getin_p('PERCENT_THROUGHFALL_PFT',throughfall_by_pft)
    throughfall_by_pft(:) = throughfall_by_pft(:) / 100. 


    !
    ! Albedo - sechiba
    !

    !Config Key   = SNOWA_AGED_VIS
    !Config Desc  = Minimum snow albedo value for each vegetation type after aging (dirty old snow), visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.74, 0.0, 0.0, 0.08, 0.24, 0.07, 0.18, 0.18, 0.33, 0.57, 0.57, 0.57, 0.57, 0.74
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('SNOWA_AGED_VIS',snowa_aged_vis)

    !Config Key   = SNOWA_AGED_NIR
    !Config Desc  = Minimum snow albedo value for each vegetation type after aging (dirty old snow), near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.50, 0.0, 0.0, 0.10, 0.37, 0.08, 0.16, 0.17, 0.27, 0.44, 0.44, 0.44, 0.44, 0.50  
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('SNOWA_AGED_NIR',snowa_aged_nir)

    !Config Key   = SNOWA_DEC_VIS
    !Config Desc  = Decay rate of snow albedo value for each vegetation type as it will be used in condveg_snow, visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.21, 0.0, 0.0, 0.14, 0.08, 0.17, 0.05, 0.06, 0.09, 0.15, 0.15, 0.15, 0.15, 0.21 
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('SNOWA_DEC_VIS',snowa_dec_vis)

    !Config Key   = SNOWA_DEC_NIR
    !Config Desc  = Decay rate of snow albedo value for each vegetation type as it will be used in condveg_snow, near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.13, 0.0, 0.0, 0.10, 0.10, 0.16, 0.04, 0.07, 0.08, 0.12, 0.12, 0.12, 0.12, 0.13
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('SNOWA_DEC_NIR',snowa_dec_nir)

    !Config Key   = ALB_LEAF_VIS
    !Config Desc  = leaf albedo of vegetation type, visible albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.00, 0.04, 0.04, 0.04, 0.04, 0.03, 0.03, 0.03, 0.03, 0.06, 0.06, 0.06, 0.06, 0.00
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('ALB_LEAF_VIS',alb_leaf_vis)

    !Config Key   = ALB_LEAF_NIR
    !Config Desc  = leaf albedo of vegetation type, near infrared albedo
    !Config if    = OK_SECHIBA
    !Config Def   = 0.00, 0.23, 0.18, 0.18, 0.20, 0.24, 0.15, 0.26, 0.20, 0.24, 0.27, 0.28, 0.26, 0.00
    !Config Help  = Values optimized for ORCHIDEE2.0
    !Config Units = [-]
    CALL getin_p('ALB_LEAF_NIR',alb_leaf_nir)

    IF ( ok_bvoc ) THEN
       !
       ! BVOC
       !

       !Config Key   = ISO_ACTIVITY
       !Config Desc  = Biogenic activity for each age class : isoprene
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.5, 1.5, 1.5, 0.5
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('ISO_ACTIVITY',iso_activity)

       !Config Key   = METHANOL_ACTIVITY
       !Config Desc  = Isoprene emission factor for each age class : methanol
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 1., 1., 0.5, 0.5
       !Config Help  =
       !Config Units = [-]
       CALL getin_p('METHANOL_ACTIVITY',methanol_activity)

       !Config Key   = EM_FACTOR_ISOPRENE
       !Config Desc  = Isoprene emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 24., 24., 8., 16., 45., 8., 18., 0.5, 12., 18., 5., 5., 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_ISOPRENE',em_factor_isoprene)

       !Config Key   = EM_FACTOR_MONOTERPENE
       !Config Desc  = Monoterpene emission factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 2.0, 2.0, 1.8, 1.4, 1.6, 1.8, 1.4, 1.8, 0.8, 0.8,  0.22, 0.22, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_MONOTERPENE',em_factor_monoterpene)

       !Config Key   = C_LDF_MONO 
       !Config Desc  = Monoterpenes fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.6
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_MONO',LDF_mono)

       !Config Key   = C_LDF_SESQ 
       !Config Desc  = Sesquiterpenes fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.5
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_SESQ',LDF_sesq)

       !Config Key   = C_LDF_METH 
       !Config Desc  = Methanol fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.8
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_METH',LDF_meth)

       !Config Key   = C_LDF_ACET 
       !Config Desc  = Acetone fraction dependancy to light
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.2
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_LDF_ACET',LDF_acet)

       !Config Key   = EM_FACTOR_APINENE 
       !Config Desc  = Alfa pinene  emission factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 1.35, 1.35, 0.85, 0.95, 0.75, 0.85, 0.60, 1.98, 0.30, 0.30, 0.09, 0.09, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_APINENE',em_factor_apinene)

       !Config Key   = EM_FACTOR_BPINENE
       !Config Desc  = Beta pinene  emission factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.30, 0.30, 0.35, 0.25, 0.20, 0.35, 0.12, 0.45, 0.16, 0.12, 0.05, 0.05, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_BPINENE',em_factor_bpinene)

       !Config Key   = EM_FACTOR_LIMONENE
       !Config Desc  = Limonene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.25, 0.25, 0.20, 0.25, 0.14, 0.20, 0.135, 0.11, 0.19, 0.42, 0.03, 0.03, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_LIMONENE',em_factor_limonene)

       !Config Key   = EM_FACTOR_MYRCENE
       !Config Desc  = Myrcene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.20, 0.20, 0.12, 0.11, 0.065, 0.12, 0.036, 0.075, 0.08,  0.085, 0.015, 0.015, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_MYRCENE',em_factor_myrcene)

       !Config Key   = EM_FACTOR_SABINENE
       !Config Desc  = Sabinene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.20, 0.20, 0.12, 0.17, 0.70, 0.12, 0.50, 0.09, 0.085, 0.075, 0.02, 0.02, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_SABINENE',em_factor_sabinene)

       !Config Key   = EM_FACTOR_CAMPHENE 
       !Config Desc  = Camphene  emission factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.15, 0.15, 0.10, 0.10, 0.01, 0.10, 0.01, 0.07, 0.07, 0.08, 0.01, 0.01, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_CAMPHENE',em_factor_camphene)

       !Config Key   = EM_FACTOR_3CARENE 
       !Config Desc  = 3-Carene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.13, 0.13, 0.42, 0.02, 0.055, 0.42,0.025, 0.125, 0.085, 0.085, 0.065, 0.065, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_3CARENE',em_factor_3carene)

       !Config Key   = EM_FACTOR_TBOCIMENE
       !Config Desc  = T-beta-ocimene  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.25, 0.25, 0.13, 0.09, 0.26, 0.13, 0.20, 0.085, 0.18, 0.18, 0.01, 0.01, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_TBOCIMENE', em_factor_tbocimene)

       !Config Key   = EM_FACTOR_OTHERMONOT
       !Config Desc  = Other monoterpenes  emission factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.17, 0.17, 0.11, 0.11, 0.125, 0.11, 0.274, 0.01, 0.15, 0.155, 0.035, 0.035, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_OTHERMONOT',em_factor_othermonot)

       !Config Key   = EM_FACTOR_SESQUITERP 
       !Config Desc  = Sesquiterpenes  emission factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.45, 0.45, 0.13, 0.3, 0.36, 0.15, 0.3, 0.25, 0.6, 0.6, 0.08, 0.08, 0.
       !Config Help  =
       !Config Units = [ugC/g/h] 
       CALL getin_p('EM_FACTOR_SESQUITERP',em_factor_sesquiterp)



       !Config Key   = C_BETA_MONO 
       !Config Desc  = Monoterpenes temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.1
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_MONO',beta_mono)

       !Config Key   = C_BETA_SESQ 
       !Config Desc  = Sesquiterpenes temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.17
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_SESQ',beta_sesq)

       !Config Key   = C_BETA_METH 
       !Config Desc  = Methanol temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.08
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_METH',beta_meth)

       !Config Key   = C_BETA_ACET 
       !Config Desc  = Acetone temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.1
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_ACET',beta_acet)

       !Config Key   = C_BETA_OXYVOC 
       !Config Desc  = Other oxygenated BVOC temperature dependency coefficient
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0.13
       !Config Help  =
       !Config Units = []
       CALL getin_p('C_BETA_OXYVOC',beta_oxyVOC)

       !Config Key   = EM_FACTOR_ORVOC
       !Config Desc  = ORVOC emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0.
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ORVOC',em_factor_ORVOC)

       !Config Key   = EM_FACTOR_OVOC
       !Config Desc  = OVOC emissions factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0.
       !Config Help  =
       !Config Units = [ugC/g/h]        
       CALL getin_p('EM_FACTOR_OVOC',em_factor_OVOC)

       !Config Key   = EM_FACTOR_MBO
       !Config Desc  = MBO emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 2.e-5, 2.e-5, 1.4, 2.e-5, 2.e-5, 0.14, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 2.e-5, 0.
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_MBO',em_factor_MBO)

       !Config Key   = EM_FACTOR_METHANOL
       !Config Desc  = Methanol emissions factor 
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.8, 0.8, 1.8, 0.9, 1.9, 1.8, 1.8, 1.8, 0.7, 0.9, 2., 2., 0.
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_METHANOL',em_factor_methanol)

       !Config Key   = EM_FACTOR_ACETONE
       !Config Desc  = Acetone emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.25, 0.25, 0.3, 0.2, 0.33, 0.3, 0.25, 0.25, 0.2, 0.2, 0.08, 0.08, 0.
       !Config Help  =
       !Config Units = [ugC/g/h]     
       CALL getin_p('EM_FACTOR_ACETONE',em_factor_acetone)

       !Config Key   = EM_FACTOR_ACETAL
       !Config Desc  = Acetaldehyde emissions factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.2, 0.2, 0.2, 0.2, 0.25, 0.25, 0.16, 0.16, 0.12, 0.12, 0.035, 0.02, 0.
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ACETAL',em_factor_acetal)

       !Config Key   = EM_FACTOR_FORMAL
       !Config Desc  = Formaldehyde emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.04, 0.04, 0.08, 0.04, 0.04, 0.04, 0.04, 0.04, 0.025, 0.025, 0.013, 0.013, 0.
       !Config Help  = 
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_FORMAL',em_factor_formal)

       !Config Key   = EM_FACTOR_ACETIC
       !Config Desc  = Acetic Acid emissions factor
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.025, 0.025,0.025,0.022,0.08,0.025,0.022,0.013,0.012,0.012,0.008,0.008, 0.
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_ACETIC',em_factor_acetic)

       !Config Key   = EM_FACTOR_FORMIC
       !Config Desc  = Formic Acid emissions factor
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 0.015, 0.015, 0.02, 0.02, 0.025, 0.025, 0.015, 0.015,0.010,0.010,0.008,0.008, 0.
       !Config Help  =
       !Config Units = [ugC/g/h]  
       CALL getin_p('EM_FACTOR_FORMIC',em_factor_formic)

       !Config Key   = EM_FACTOR_NO_WET
       !Config Desc  = NOx emissions factor wet soil emissions and exponential dependancy factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 2.6, 0.06, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.36, 0.36, 0.36, 0.36, 0.
       !Config Help  =
       !Config Units = [ngN/m^2/s]
       CALL getin_p('EM_FACTOR_NO_WET',em_factor_no_wet)

       !Config Key   = EM_FACTOR_NO_DRY
       !Config Desc  = NOx emissions factor dry soil emissions and exponential dependancy factor 
       !Config if    = CHEMISTRY_BVOC
       !Config Def   = 0., 8.60, 0.40, 0.22, 0.22, 0.22, 0.22, 0.22, 0.22, 2.65, 2.65, 2.65, 2.65, 0.
       !Config Help  =
       !Config Units = [ngN/m^2/s] 
       CALL getin_p('EM_FACTOR_NO_DRY',em_factor_no_dry)

       !Config Key   = LARCH
       !Config Desc  = Larcher 1991 SAI/LAI ratio
       !Config if    = CHEMISTRY_BVOC 
       !Config Def   = 0., 0.015, 0.015, 0.003, 0.005, 0.005, 0.003, 0.005, 0.003, 0.005, 0.005, 0.008, 0.008, 0.
       !Config Help  =
       !Config Units = [-]  
       CALL getin_p('LARCH',Larch)

    ENDIF ! (ok_bvoc)

  END SUBROUTINE config_sechiba_pft_parameters


!! ================================================================================================================================
!! SUBROUTINE   : config_stomate_pft_parameters 
!!
!>\BRIEF         This subroutine will read the imposed values for the stomate pft
!! parameters. It is not called if IMPOSE_PARAM is set to NO.
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE config_stomate_pft_parameters

    IMPLICIT NONE

    !! 0. Variables and parameters declaration

    !! 0.4 Local variable

    !_ ================================================================================================================================

    !
    ! Vegetation structure
    !

    !Config Key   = SLA
    !Config Desc  = specif leaf area 
    !Config if    = OK_STOMATE
    !Config Def   = 1.5E-2, 1.53E-2, 2.6E-2, 9.26E-3, 2E-2, 2.6E-2, 9.26E-3, 2.6E-2, 1.9E-2, 2.6E-2, 2.6E-2, 2.6E-2, 2.6E-2, 1.5E-2
    !Config Help  =
    !Config Units = [m^2/gC]
    CALL getin_p('SLA',sla)


    !Config Key   = AVAILABILITY_FACT 
    !Config Desc  = Calculate dynamic mortality in lpj_gap, pft dependent parameter
    !Config If    = OK_STOMATE 
    !Config Def   = undef, 0.14, 0.14, 0.10, 0.10, 0.10, 0.05, 0.05, 0.05, undef, undef, undef, undef, undef 
    !Config Help  = 
    !Config Units = [-]   
    CALL getin_p('AVAILABILITY_FACT',availability_fact)

    !
    ! Allocation - stomate
    !
    !
    !Config Key   = R0 
    !Config Desc  = Standard root allocation 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, .30, .30, .30, .30, .30, .30, .30, .30, .30, .30, .30, .30, undef
    !Config Help  = 
    !Config Units = [-]    
    CALL getin_p('R0',R0)

    !Config Key   = S0 
    !Config Desc  = Standard sapwood allocation 
    !Config If    = OK_STOMATE 
    !Config Def   = undef, .25, .25, .30, .30, .30, .30, .30, .30, .30, .30, .30, .30, undef
    !Config Help  = 
    !Config Units = [-]    
    CALL getin_p('S0',S0)

    !
    ! Respiration - stomate
    !

    !Config Key   = FRAC_GROWTHRESP
    !Config Desc  = fraction of GPP which is lost as growth respiration
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.35, 0.35, 0.28, 0.28, 0.28, 0.35, 0.35, 0.35, 0.28, 0.28, 0.28, 0.28, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('FRAC_GROWTHRESP',frac_growthresp) 

    !Config Key   = MAINT_RESP_SLOPE_C
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant c of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.12, 0.12, 0.16, 0.16, 0.16, 0.25, 0.25, 0.25, 0.16, 0.12, 0.16, 0.12, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('MAINT_RESP_SLOPE_C',maint_resp_slope_c) 

    !Config Key   = MAINT_RESP_SLOPE_B
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant b of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, .0, .0, .0, .0, .0, .0, .0, .0, -.00133, .0, -.00133, .0, undef 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('MAINT_RESP_SLOPE_B',maint_resp_slope_b)

    !Config Key   = MAINT_RESP_SLOPE_A
    !Config Desc  = slope of maintenance respiration coefficient (1/K), constant a of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, .0, undef    
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('MAINT_RESP_SLOPE_A',maint_resp_slope_a)

    !Config Key   = CM_ZERO_LEAF
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for leaves, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 2.35E-3, 2.62E-3, 1.01E-3, 2.35E-3, 2.62E-3, 1.01E-3,2.62E-3, 2.05E-3, 2.62E-3, 2.62E-3, 2.62E-3, 2.62E-3, undef
    !Config Help  =
    !Config Units = [g/g/day]
    CALL getin_p('CM_ZERO_LEAF',cm_zero_leaf)

    !Config Key   = CM_ZERO_SAPABOVE
    !Config Desc  = maintenance respiration coefficient at 0 deg C,for sapwood above, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, undef
    !Config Help  =
    !Config Units = [g/g/day]
    CALL getin_p('CM_ZERO_SAPABOVE',cm_zero_sapabove)

    !Config Key   = CM_ZERO_SAPBELOW
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for sapwood below, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, undef 
    !Config Help  =
    !Config Units = [g/g/day]
    CALL getin_p('CM_ZERO_SAPBELOW',cm_zero_sapbelow)

    !Config Key   = CM_ZERO_HEARTABOVE
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for heartwood above, tabulated
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., undef
    !Config Help  =
    !Config Units = [g/g/day]
    CALL getin_p('CM_ZERO_HEARTABOVE',cm_zero_heartabove)

    !Config Key   = CM_ZERO_HEARTBELOW
    !Config Desc  = maintenance respiration coefficient at 0 deg C,for heartwood below, tabulated
    !Config if    = OK_STOMATE 
    !Config Def   = undef, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., undef
    !Config Help  =
    !Config Units = [g/g/day] 
    CALL getin_p('CM_ZERO_HEARTBELOW',cm_zero_heartbelow)

    !Config Key   = CM_ZERO_ROOT
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for roots, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef,1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3,1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, 1.67E-3, undef
    !Config Help  =
    !Config Units = [g/g/day] 
    CALL getin_p('CM_ZERO_ROOT',cm_zero_root)

    !Config Key   = CM_ZERO_FRUIT
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for fruits, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4,1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, undef    
    !Config Help  =
    !Config Units = [g/g/day] 
    CALL getin_p('CM_ZERO_FRUIT',cm_zero_fruit)

    !Config Key   = CM_ZERO_CARBRES
    !Config Desc  = maintenance respiration coefficient at 0 deg C, for carbohydrate reserve, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4,1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, 1.19E-4, undef
    !Config Help  =
    !Config Units = [g/g/day] 
    CALL getin_p('CM_ZERO_CARBRES',cm_zero_carbres)

    !
    ! Fire - stomate
    !

    !Config Key   = FLAM
    !Config Desc  = flamability: critical fraction of water holding capacity
    !Config if    = OK_STOMATE
    !Config Def   = undef, .15, .25, .25, .25, .25, .25, .25, .25, .25, .25, .35, .35, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('FLAM',flam)

    !Config Key   = RESIST
    !Config Desc  = fire resistance
    !Config if    = OK_STOMATE
    !Config Def   = undef, .95, .90, .12, .50, .12, .12, .12, .12, .0, .0, .0, .0, undef 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('RESIST',resist)

    !
    ! Flux - LUC
    !

    !Config Key   = COEFF_LCCHANGE_1
    !Config Desc  = Coeff of biomass export for the year
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.897, 0.897, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, 0.597, undef 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_1',coeff_lcchange_1)

    !Config Key   = COEFF_LCCHANGE_10
    !Config Desc  = Coeff of biomass export for the decade
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.103, 0.103, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.403, 0.299, 0.403, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_10',coeff_lcchange_10)

    !Config Key   = COEFF_LCCHANGE_100
    !Config Desc  = Coeff of biomass export for the century
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0., 0., 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0.104, 0., 0.104, 0., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('COEFF_LCCHANGE_100',coeff_lcchange_100)

    !
    ! Phenology
    !

    !Config Key   = LAI_MAX_TO_HAPPY
    !Config Desc  = threshold of LAI below which plant uses carbohydrate reserves
    !Config if    = OK_STOMATE
    !Config Def   = undef, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, .5, undef 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('LAI_MAX_TO_HAPPY',lai_max_to_happy) 

    !Config Key   = LAI_MAX
    !Config Desc  = maximum LAI, PFT-specific
    !Config if    = OK_STOMATE
    !Config Def   = undef, 7.0, 5.0, 5.0, 4.0, 5.0, 3.5, 4.0, 3.0, 2.5, 2.0, 5.0, 5.0, undef
    !Config Help  =
    !Config Units = [m^2/m^2]
    CALL getin_p('LAI_MAX',lai_max)

    !Config Key   = PHENO_TYPE
    !Config Desc  = type of phenology, 0=bare ground 1=evergreen,  2=summergreen,  3=raingreen,  4=perennial
    !Config if    = OK_STOMATE
    !Config Def   = 0, 1, 3, 1, 1, 2, 1, 2, 2, 4, 4, 2, 3, 0
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_TYPE',pheno_type)

    !
    ! Phenology : Leaf Onset
    !

    !Config Key   = PHENO_GDD_CRIT_C
    !Config Desc  = critical gdd, tabulated (C), constant c of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 270., 400., 125., 400., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_C',pheno_gdd_crit_c)

    !Config Key   = PHENO_GDD_CRIT_B
    !Config Desc  = critical gdd, tabulated (C), constant b of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef,undef, undef, 6.25, 0., 0., 0., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_B',pheno_gdd_crit_b)

    !Config Key   = PHENO_GDD_CRIT_A
    !Config Desc  = critical gdd, tabulated (C), constant a of aT^2+bT+c
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 0.03125,  0., 0., 0., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('PHENO_GDD_CRIT_A',pheno_gdd_crit_a)

    !Config Key   = PHENO_MOIGDD_T_CRIT
    !Config Desc  = Average temperature threashold for C4 grass used in pheno_moigdd
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 22.0, undef, undef, undef
    !Config Help  =
    !Config Units = [C]
    CALL getin_p('PHENO_MOIGDD_T_CRIT',pheno_moigdd_t_crit)

    !Config Key   = NGD_CRIT
    !Config Desc  = critical ngd, tabulated. Threshold -5 degrees
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, 0., undef, undef, undef, undef, undef, undef
    !Config Help  = NGD : Number of Growing Days.
    !Config Units = [days]
    CALL getin_p('NGD_CRIT',ngd_crit)

    !Config Key   = NCDGDD_TEMP
    !Config Desc  = critical temperature for the ncd vs. gdd function in phenology
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 5., undef, 0., undef, undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [C] 
    CALL getin_p('NCDGDD_TEMP',ncdgdd_temp)

    !Config Key   = HUM_FRAC
    !Config Desc  = critical humidity (relative to min/max) for phenology
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .5, undef, undef, undef, undef, undef,  undef, .5, .5, .5, .5, undef     
    !Config Help  =
    !Config Units = [%]
    CALL getin_p('HUM_FRAC',hum_frac)

    !Config Key   = HUM_MIN_TIME
    !Config Desc  = minimum time elapsed since moisture minimum
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 50., undef, undef, undef, undef, undef, undef, 35., 35., 75., 75., undef
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('HUM_MIN_TIME',hum_min_time)

    !Config Key   = TAU_SAP
    !Config Desc  = sapwood -> heartwood conversion time
    !Config if    = OK_STOMATE
    !Config Def   = undef, 730., 730., 730., 730., 730., 730., 730., 730., undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('TAU_SAP',tau_sap)

    !Config Key   = TAU_LEAFINIT
    !Config Desc  = time to attain the initial foliage using the carbohydrate reserve
    !Config if    = OK_STOMATE
    !Config Def   = undef, 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., undef
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('TAU_LEAFINIT',tau_leafinit) 

    !Config Key   = TAU_FRUIT
    !Config Desc  = fruit lifetime
    !Config if    = OK_STOMATE
    !Config Def   = undef, 90., 90., 90., 90., 90., 90., 90., 90., undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('TAU_FRUIT',tau_fruit)

    !Config Key   = ECUREUIL
    !Config Desc  = fraction of primary leaf and root allocation put into reserve
    !Config if    = OK_STOMATE
    !Config Def   = undef, .0, 1., .0, .0, 1., .0, 1., 1., 1., 1., 1., 1., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ECUREUIL',ecureuil)

    !Config Key   = ALLOC_MIN
    !Config Desc  = minimum allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, undef, undef, undef, undef, undef 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ALLOC_MIN',alloc_min)

    !Config Key   = ALLOC_MAX
    !Config Desc  = maximum allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('ALLOC_MAX',alloc_max)

    !Config Key   = DEMI_ALLOC 
    !Config Desc  = mean allocation above/below = f(age) - 30/01/04 NV/JO/PF
    !Config if    = OK_STOMATE
    !Config Def   = undef, 5., 5., 5., 5., 5., 5., 5., 5., undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('DEMI_ALLOC',demi_alloc)

    !Config Key   = LEAFLIFE_TAB
    !Config Desc  = leaf longevity
    !Config if    = OK_STOMATE
    !Config Def   = undef, .5, 2., .33, 1., 2., .33, 2., 2., 2., 2., 2., 2., undef 
    !Config Help  =
    !Config Units = [years]
    CALL getin_p('LEAFLIFE_TAB',leaflife_tab)

    !
    ! Phenology : Senescence
    !
    !
    !Config Key   = LEAFFALL
    !Config Desc  = length of death of leaves, tabulated 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 10., undef, undef, 10., undef, 10., 10., 10., 10., 10., 10., undef 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LEAFFALL',leaffall)

    !Config Key   = LEAFAGECRIT
    !Config Desc  = critical leaf age, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, 730., 180., 910., 730., 180., 910., 180., 180., 120., 120., 90., 90., undef  
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('LEAFAGECRIT',leafagecrit) 

    !Config Key   = SENESCENCE_TYPE
    !Config Desc  = type of senescence, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = none, none, dry, none, none, cold, none, cold, cold, mixed, mixed, mixed, mixed, none 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TYPE',senescence_type) 

    !Config Key   = SENESCENCE_HUM
    !Config Desc  = critical relative moisture availability for senescence
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .3, undef, undef, undef, undef, undef, undef, .2, .2, .3, .2, undef 
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('SENESCENCE_HUM',senescence_hum)

    !Config Key   = NOSENESCENCE_HUM
    !Config Desc  = relative moisture availability above which there is no humidity-related senescence
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, .8, undef, undef, undef, undef, undef, undef, .3, .3, .3, .3, undef 
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('NOSENESCENCE_HUM',nosenescence_hum) 

    !Config Key   = MAX_TURNOVER_TIME
    !Config Desc  = maximum turnover time for grasse
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef,  80.,  80., 80., 80., undef 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('MAX_TURNOVER_TIME',max_turnover_time)

    !Config Key   = MIN_TURNOVER_TIME
    !Config Desc  = minimum turnover time for grasse 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, 10., 10., 10., 10., undef 
    !Config Help  =
    !Config Units = [days]
    CALL getin_p('MIN_TURNOVER_TIME',min_turnover_time)

    !Config Key   = MIN_LEAF_AGE_FOR_SENESCENCE
    !Config Desc  = minimum leaf age to allow senescence g
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, 90., undef, undef, 90., undef, 60., 60., 30., 30., 30., 30., undef
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p('MIN_LEAF_AGE_FOR_SENESCENCE',min_leaf_age_for_senescence)

    !Config Key   = SENESCENCE_TEMP_C
    !Config Desc  = critical temperature for senescence (C), constant c of aT^2+bT+c, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 12., undef, 7., 2., -1.375, 5., 5., 10., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TEMP_C',senescence_temp_c)

    !Config Key   = SENESCENCE_TEMP_B
    !Config Desc  = critical temperature for senescence (C), constant b of aT^2+bT+c ,tabulated
    !Config if    = OK_STOMATE 
    !Config Def   = undef, undef, undef, undef, undef, 0., undef, 0., 0., .1, 0., 0., 0., undef
    !Config Help  =
    !Config Units = [-]
    CALL getin_p('SENESCENCE_TEMP_B',senescence_temp_b)

    !Config Key   = SENESCENCE_TEMP_A
    !Config Desc  = critical temperature for senescence (C), constant a of aT^2+bT+c , tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, 0., undef, 0., 0.,.00375, 0., 0., 0., undef 
    !Config Help  =
    !Config Units = [-] 
    CALL getin_p('SENESCENCE_TEMP_A',senescence_temp_a)

    !Config Key   = GDD_SENESCENCE
    !Config Desc  = minimum gdd to allow senescence of crops  
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, undef, 950., 4000., undef
    !Config Help  =
    !Config Units = [days] 
    CALL getin_p("GDD_SENESCENCE", gdd_senescence)

    !Config Key   = ALWAYS_INIT
    !Config Desc  = Take carbon from atmosphere if carbohydrate reserve too small
    !Config if    = OK_STOMATE
    !Config Def   = y, y, y, y, y, y, y, y, y, y, n, y, y, y
    !Config Help  =
    !Config Units = [BOOLEAN]
    CALL getin_p('ALWAYS_INIT',always_init)

    !
    ! DGVM
    !

    !Config Key   = RESIDENCE_TIME
    !Config Desc  = residence time of trees
    !Config if    = OK_DGVM and NOT(LPJ_GAP_CONST_MORT)
    !Config Def   = undef, 30.0, 30.0, 40.0, 40.0, 40.0, 80.0, 80.0, 80.0, 0.0, 0.0, 0.0, 0.0, undef 
    !Config Help  =
    !Config Units = [years]
    CALL getin_p('RESIDENCE_TIME',residence_time)

    !Config Key   = TMIN_CRIT
    !Config Desc  = critical tmin, tabulated
    !Config if    = OK_STOMATE
    !Config Def   = undef,  0.0, 0.0, -30.0, -14.0, -30.0, -45.0, -45.0, undef, undef, undef, undef, undef, undef
    !Config Help  = 
    !Config Units = [C]
    CALL getin_p('TMIN_CRIT',tmin_crit)

    !Config Key   = TCM_CRIT
    !Config Desc  = critical tcm, tabulated 
    !Config if    = OK_STOMATE
    !Config Def   = undef, undef, undef, 5.0, 15.5, 15.5, -8.0, -8.0, -8.0, undef, undef, undef, undef, undef
    !Config Help  =
    !Config Units = [C]
    CALL getin_p('TCM_CRIT',tcm_crit)

  END SUBROUTINE config_stomate_pft_parameters


!! ================================================================================================================================
!! SUBROUTINE   : pft_parameters_clear
!!
!>\BRIEF         This subroutine deallocates memory at the end of the simulation. 
!!
!! DESCRIPTION  : None
!!
!! RECENT CHANGE(S): None
!!
!! MAIN OUTPUT VARIABLE(S): None
!!
!! REFERENCE(S) : None
!!
!! FLOWCHART    : None
!! \n
!_ ================================================================================================================================

  SUBROUTINE pft_parameters_clear

    l_first_pft_parameters = .TRUE.

    IF (ALLOCATED(pft_to_mtc)) DEALLOCATE(pft_to_mtc)
    IF (ALLOCATED(PFT_name)) DEALLOCATE(PFT_name)
    IF (ALLOCATED(veget_ori_fixed_test_1)) DEALLOCATE(veget_ori_fixed_test_1)   
    IF (ALLOCATED(llaimax)) DEALLOCATE(llaimax)
    IF (ALLOCATED(llaimin)) DEALLOCATE(llaimin)
    IF (ALLOCATED(height_presc)) DEALLOCATE(height_presc)   
    IF (ALLOCATED(z0_over_height)) DEALLOCATE(z0_over_height)   
    IF (ALLOCATED(ratio_z0m_z0h)) DEALLOCATE(ratio_z0m_z0h)   
    IF (ALLOCATED(type_of_lai)) DEALLOCATE(type_of_lai)
    IF (ALLOCATED(is_tree)) DEALLOCATE(is_tree)
    IF (ALLOCATED(natural)) DEALLOCATE(natural)
    IF (ALLOCATED(is_deciduous)) DEALLOCATE(is_deciduous)
    IF (ALLOCATED(is_evergreen)) DEALLOCATE(is_evergreen)
    IF (ALLOCATED(is_needleleaf)) DEALLOCATE(is_needleleaf)
    IF (ALLOCATED(is_tropical)) DEALLOCATE(is_tropical)
    IF (ALLOCATED(humcste)) DEALLOCATE(humcste)
    IF (ALLOCATED(pref_soil_veg)) DEALLOCATE(pref_soil_veg)
    IF (ALLOCATED(is_c4)) DEALLOCATE(is_c4)  
    IF (ALLOCATED(vcmax_fix)) DEALLOCATE(vcmax_fix)
    IF (ALLOCATED(downregulation_co2_coeff)) DEALLOCATE(downregulation_co2_coeff) 
    IF (ALLOCATED(downregulation_co2_coeff_new)) DEALLOCATE(downregulation_co2_coeff_new) 
    IF (ALLOCATED(E_KmC)) DEALLOCATE(E_KmC)
    IF (ALLOCATED(E_KmO)) DEALLOCATE(E_KmO)
    IF (ALLOCATED(E_Sco)) DEALLOCATE(E_Sco)
    IF (ALLOCATED(E_gamma_star)) DEALLOCATE(E_gamma_star)
    IF (ALLOCATED(E_Vcmax)) DEALLOCATE(E_Vcmax)
    IF (ALLOCATED(E_Jmax)) DEALLOCATE(E_Jmax)
    IF (ALLOCATED(aSV)) DEALLOCATE(aSV)
    IF (ALLOCATED(bSV)) DEALLOCATE(bSV)
    IF (ALLOCATED(tphoto_min)) DEALLOCATE(tphoto_min)
    IF (ALLOCATED(tphoto_max)) DEALLOCATE(tphoto_max)
    IF (ALLOCATED(aSJ)) DEALLOCATE(aSJ)
    IF (ALLOCATED(bSJ)) DEALLOCATE(bSJ)
    IF (ALLOCATED(D_Vcmax)) DEALLOCATE(D_Vcmax)
    IF (ALLOCATED(D_Jmax)) DEALLOCATE(D_Jmax)
    IF (ALLOCATED(E_gm)) DEALLOCATE(E_gm) 
    IF (ALLOCATED(S_gm)) DEALLOCATE(S_gm) 
    IF (ALLOCATED(D_gm)) DEALLOCATE(D_gm) 
    IF (ALLOCATED(E_Rd)) DEALLOCATE(E_Rd)
    IF (ALLOCATED(Vcmax25)) DEALLOCATE(Vcmax25)
    IF (ALLOCATED(arJV)) DEALLOCATE(arJV)
    IF (ALLOCATED(brJV)) DEALLOCATE(brJV)
    IF (ALLOCATED(KmC25)) DEALLOCATE(KmC25)
    IF (ALLOCATED(KmO25)) DEALLOCATE(KmO25)
    IF (ALLOCATED(Sco25)) DEALLOCATE(Sco25)
    IF (ALLOCATED(gm25)) DEALLOCATE(gm25) 
    IF (ALLOCATED(gamma_star25)) DEALLOCATE(gamma_star25)
    IF (ALLOCATED(a1)) DEALLOCATE(a1)
    IF (ALLOCATED(b1)) DEALLOCATE(b1)
    IF (ALLOCATED(g0)) DEALLOCATE(g0)
    IF (ALLOCATED(h_protons)) DEALLOCATE(h_protons)
    IF (ALLOCATED(fpsir)) DEALLOCATE(fpsir)
    IF (ALLOCATED(fQ)) DEALLOCATE(fQ)
    IF (ALLOCATED(fpseudo)) DEALLOCATE(fpseudo)
    IF (ALLOCATED(kp)) DEALLOCATE(kp)
    IF (ALLOCATED(alpha)) DEALLOCATE(alpha)
    IF (ALLOCATED(gbs)) DEALLOCATE(gbs)
    IF (ALLOCATED(theta)) DEALLOCATE(theta)
    IF (ALLOCATED(alpha_LL)) DEALLOCATE(alpha_LL)
    IF (ALLOCATED(stress_vcmax)) DEALLOCATE(stress_vcmax)
    IF (ALLOCATED(stress_gs)) DEALLOCATE(stress_gs)
    IF (ALLOCATED(stress_gm)) DEALLOCATE(stress_gm)
    IF (ALLOCATED(ext_coeff)) DEALLOCATE(ext_coeff)
    IF (ALLOCATED(ext_coeff_vegetfrac)) DEALLOCATE(ext_coeff_vegetfrac)
    IF (ALLOCATED(rveg_pft)) DEALLOCATE(rveg_pft)
    IF (ALLOCATED(rstruct_const)) DEALLOCATE(rstruct_const)
    IF (ALLOCATED(kzero)) DEALLOCATE(kzero)
    IF (ALLOCATED(wmax_veg)) DEALLOCATE(wmax_veg)
    IF (ALLOCATED(throughfall_by_pft)) DEALLOCATE(throughfall_by_pft)
    IF (ALLOCATED(snowa_aged_vis)) DEALLOCATE(snowa_aged_vis)
    IF (ALLOCATED(snowa_aged_nir)) DEALLOCATE(snowa_aged_nir)
    IF (ALLOCATED(snowa_dec_vis)) DEALLOCATE(snowa_dec_vis)
    IF (ALLOCATED(snowa_dec_nir)) DEALLOCATE(snowa_dec_nir)
    IF (ALLOCATED(alb_leaf_vis)) DEALLOCATE(alb_leaf_vis)
    IF (ALLOCATED(alb_leaf_nir)) DEALLOCATE(alb_leaf_nir)   
    IF (ALLOCATED(em_factor_isoprene)) DEALLOCATE(em_factor_isoprene)
    IF (ALLOCATED(em_factor_monoterpene)) DEALLOCATE(em_factor_monoterpene)
    IF (ALLOCATED(em_factor_apinene)) DEALLOCATE(em_factor_apinene)
    IF (ALLOCATED(em_factor_bpinene)) DEALLOCATE(em_factor_bpinene)
    IF (ALLOCATED(em_factor_limonene)) DEALLOCATE(em_factor_limonene)
    IF (ALLOCATED(em_factor_myrcene)) DEALLOCATE(em_factor_myrcene)
    IF (ALLOCATED(em_factor_sabinene)) DEALLOCATE(em_factor_sabinene)
    IF (ALLOCATED(em_factor_camphene)) DEALLOCATE(em_factor_camphene)
    IF (ALLOCATED(em_factor_3carene)) DEALLOCATE(em_factor_3carene)
    IF (ALLOCATED(em_factor_tbocimene)) DEALLOCATE(em_factor_tbocimene)
    IF (ALLOCATED(em_factor_othermonot)) DEALLOCATE(em_factor_othermonot)
    IF (ALLOCATED(em_factor_sesquiterp)) DEALLOCATE(em_factor_sesquiterp)
    IF (ALLOCATED(em_factor_ORVOC)) DEALLOCATE(em_factor_ORVOC)
    IF (ALLOCATED(em_factor_OVOC)) DEALLOCATE(em_factor_OVOC)
    IF (ALLOCATED(em_factor_MBO)) DEALLOCATE(em_factor_MBO)
    IF (ALLOCATED(em_factor_methanol)) DEALLOCATE(em_factor_methanol)
    IF (ALLOCATED(em_factor_acetone)) DEALLOCATE(em_factor_acetone)
    IF (ALLOCATED(em_factor_acetal)) DEALLOCATE(em_factor_acetal)
    IF (ALLOCATED(em_factor_formal)) DEALLOCATE(em_factor_formal)
    IF (ALLOCATED(em_factor_acetic)) DEALLOCATE(em_factor_acetic)
    IF (ALLOCATED(em_factor_formic)) DEALLOCATE(em_factor_formic)
    IF (ALLOCATED(em_factor_no_wet)) DEALLOCATE(em_factor_no_wet)
    IF (ALLOCATED(em_factor_no_dry)) DEALLOCATE(em_factor_no_dry)
    IF (ALLOCATED(Larch)) DEALLOCATE(Larch)
    IF (ALLOCATED(leaf_tab)) DEALLOCATE(leaf_tab)
    IF (ALLOCATED(sla)) DEALLOCATE(sla)
    IF (ALLOCATED(availability_fact)) DEALLOCATE(availability_fact)
    IF (ALLOCATED(R0)) DEALLOCATE(R0)
    IF (ALLOCATED(S0)) DEALLOCATE(S0)
    IF (ALLOCATED(L0)) DEALLOCATE(L0)
    IF (ALLOCATED(frac_growthresp)) DEALLOCATE(frac_growthresp)
    IF (ALLOCATED(maint_resp_slope)) DEALLOCATE(maint_resp_slope)
    IF (ALLOCATED(maint_resp_slope_c)) DEALLOCATE(maint_resp_slope_c)
    IF (ALLOCATED(maint_resp_slope_b)) DEALLOCATE(maint_resp_slope_b)
    IF (ALLOCATED(maint_resp_slope_a)) DEALLOCATE(maint_resp_slope_a)
    IF (ALLOCATED(coeff_maint_zero)) DEALLOCATE(coeff_maint_zero)
    IF (ALLOCATED(cm_zero_leaf)) DEALLOCATE(cm_zero_leaf)
    IF (ALLOCATED(cm_zero_sapabove)) DEALLOCATE(cm_zero_sapabove)
    IF (ALLOCATED(cm_zero_sapbelow)) DEALLOCATE(cm_zero_sapbelow)
    IF (ALLOCATED(cm_zero_heartabove)) DEALLOCATE(cm_zero_heartabove)
    IF (ALLOCATED(cm_zero_heartbelow)) DEALLOCATE(cm_zero_heartbelow)
    IF (ALLOCATED(cm_zero_root)) DEALLOCATE(cm_zero_root)
    IF (ALLOCATED(cm_zero_fruit)) DEALLOCATE(cm_zero_fruit)
    IF (ALLOCATED(cm_zero_carbres)) DEALLOCATE(cm_zero_carbres)
    IF (ALLOCATED(flam)) DEALLOCATE(flam)
    IF (ALLOCATED(resist)) DEALLOCATE(resist)
    IF (ALLOCATED(coeff_lcchange_1)) DEALLOCATE(coeff_lcchange_1)
    IF (ALLOCATED(coeff_lcchange_10)) DEALLOCATE(coeff_lcchange_10)
    IF (ALLOCATED(coeff_lcchange_100)) DEALLOCATE(coeff_lcchange_100)
    IF (ALLOCATED(lai_max_to_happy)) DEALLOCATE(lai_max_to_happy)
    IF (ALLOCATED(lai_max)) DEALLOCATE(lai_max)
    IF (ALLOCATED(pheno_model)) DEALLOCATE(pheno_model)
    IF (ALLOCATED(pheno_type)) DEALLOCATE(pheno_type)
    IF (ALLOCATED(pheno_gdd_crit_c)) DEALLOCATE(pheno_gdd_crit_c)
    IF (ALLOCATED(pheno_gdd_crit_b)) DEALLOCATE(pheno_gdd_crit_b)
    IF (ALLOCATED(pheno_gdd_crit_a)) DEALLOCATE(pheno_gdd_crit_a)
    IF (ALLOCATED(pheno_gdd_crit)) DEALLOCATE(pheno_gdd_crit)
    IF (ALLOCATED(pheno_moigdd_t_crit)) DEALLOCATE(pheno_moigdd_t_crit)
    IF (ALLOCATED(ngd_crit)) DEALLOCATE(ngd_crit)
    IF (ALLOCATED(ncdgdd_temp)) DEALLOCATE(ncdgdd_temp)
    IF (ALLOCATED(hum_frac)) DEALLOCATE(hum_frac)
    IF (ALLOCATED(hum_min_time)) DEALLOCATE(hum_min_time)
    IF (ALLOCATED(tau_sap)) DEALLOCATE(tau_sap)
    IF (ALLOCATED(tau_leafinit)) DEALLOCATE(tau_leafinit)
    IF (ALLOCATED(tau_fruit)) DEALLOCATE(tau_fruit)
    IF (ALLOCATED(ecureuil)) DEALLOCATE(ecureuil)
    IF (ALLOCATED(alloc_min)) DEALLOCATE(alloc_min)
    IF (ALLOCATED(alloc_max)) DEALLOCATE(alloc_max)
    IF (ALLOCATED(demi_alloc)) DEALLOCATE(demi_alloc)
    IF (ALLOCATED(leaflife_tab)) DEALLOCATE(leaflife_tab)
    IF (ALLOCATED(leaffall)) DEALLOCATE(leaffall)
    IF (ALLOCATED(leafagecrit)) DEALLOCATE(leafagecrit)
    IF (ALLOCATED(senescence_type)) DEALLOCATE(senescence_type)
    IF (ALLOCATED(senescence_hum)) DEALLOCATE(senescence_hum)
    IF (ALLOCATED(nosenescence_hum)) DEALLOCATE(nosenescence_hum)
    IF (ALLOCATED(max_turnover_time)) DEALLOCATE(max_turnover_time)
    IF (ALLOCATED(min_turnover_time)) DEALLOCATE(min_turnover_time)
    IF (ALLOCATED(min_leaf_age_for_senescence)) DEALLOCATE(min_leaf_age_for_senescence)
    IF (ALLOCATED(senescence_temp_c)) DEALLOCATE(senescence_temp_c)
    IF (ALLOCATED(senescence_temp_b)) DEALLOCATE(senescence_temp_b)
    IF (ALLOCATED(senescence_temp_a)) DEALLOCATE(senescence_temp_a)
    IF (ALLOCATED(senescence_temp)) DEALLOCATE(senescence_temp)
    IF (ALLOCATED(gdd_senescence)) DEALLOCATE(gdd_senescence)
    IF (ALLOCATED(always_init)) DEALLOCATE(always_init)
    IF (ALLOCATED(residence_time)) DEALLOCATE(residence_time)
    IF (ALLOCATED(tmin_crit)) DEALLOCATE(tmin_crit)
    IF (ALLOCATED(tcm_crit)) DEALLOCATE(tcm_crit)
    IF (ALLOCATED(lai_initmin)) DEALLOCATE(lai_initmin)
    IF (ALLOCATED(bm_sapl)) DEALLOCATE(bm_sapl)
    IF (ALLOCATED(migrate)) DEALLOCATE(migrate)
    IF (ALLOCATED(maxdia)) DEALLOCATE(maxdia)
    IF (ALLOCATED(cn_sapl)) DEALLOCATE(cn_sapl)
    IF (ALLOCATED(leaf_timecst)) DEALLOCATE(leaf_timecst)

  END SUBROUTINE pft_parameters_clear

END MODULE pft_parameters
