! =================================================================================================================================
! PROGRAM       : forcesoil
!
! CONTACT	: orchidee-help _at_ listes.ipsl.fr
!
! LICENCE      	: IPSL (2006). This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF        This subroutine runs the soilcarbon submodel using specific initial conditions 
!! and driving variables in order to obtain soil carbon stocks closed to the steady-state values 
!! quicker than when using the ''full'' ORCHIDEE.  
!!	
!!\n DESCRIPTION: None
!! This subroutine computes the soil carbon stocks by calling the soilcarbon routine at each time step. \n
!! The aim is to obtain soil carbon stocks closed to the steady-state values and ultimately to create  
!! an updated stomate restart file for the stomate component. The state variables of the subsystem are the clay content 
!! (fixed value) and the soil carbon stocks. Initial conditions for the state variables are read in an  
!! input stomate restart file. Driving variables are Soil carbon input, Water and Temperature stresses on 
!! Organic Matter decomposition. Driving variables are read from a specific forcing file produced by a former run of ORCHIDEE
!! (SECHIBA+STOMATE). \n 
!! The FORCESOIL program first consists in reading a set of input files, allocating variables and 
!! preparing output stomate restart file. \n                                                             
!! Then, a loop over time is performed in which the soilcarbon routine is called at each time step. \n
!! Last, final values of the soil carbon stocks are written into the output stomate restart file. \n
!! No flag is associated with the use of the FORCESOIL program. \n
!!
!! RECENT CHANGE(S): None
!!
!! REFERENCE(S) : None    
!!
!! FLOWCHART    : None
!!
!! SVN		:
!! $HeadURL: $ 
!! $Date: $
!! $Revision: $
!! \n
!_ =================================================================================================================================

PROGRAM forcesoil
 
  USE netcdf
  !-
  USE defprec
  USE constantes
  USE constantes_mtc
  USE pft_parameters 
  USE stomate_data
  USE ioipsl_para
  USE mod_orchidee_para
  USE stomate_soilcarbon

  !-
  IMPLICIT NONE
  !-
  CHARACTER(LEN=80)                          :: sto_restname_in,sto_restname_out
  INTEGER(i_std)                             :: iim,jjm                !! Indices (unitless)

  INTEGER(i_std),PARAMETER                   :: llm = 1                !! Vertical Layers (requested by restini routine) (unitless)
  INTEGER(i_std)                             :: kjpindex               !! Domain size (unitless)

  INTEGER(i_std)                             :: itau_dep,itau_len      !! Time step read in the restart file (?) 
                                                                       !! and number of time steps of the simulation (unitless) 
  CHARACTER(LEN=30)                          :: time_str               !! Length of the simulation (year)
  REAL(r_std)                                :: dt_files               !! time step between two successive itaus (?) 
                                                                       !! (requested by restini routine) (seconds)
  REAL(r_std)                                :: date0                  !! Time at which itau = 0 (requested by restini routine) (?)
  INTEGER(i_std)                             :: rest_id_sto            !! ID of the input restart file (unitless)
  CHARACTER(LEN=20), SAVE                    :: thecalendar = 'noleap' !! Type of calendar defined in the input restart file 
                                                                       !! (unitless)
  !-
  CHARACTER(LEN=100)                         :: Cforcing_name          !! Name of the forcing file (unitless)
  INTEGER                                    :: Cforcing_id            !! ID of the forcing file (unitless)
  INTEGER                                    :: v_id                   !! ID of the variable 'Index' stored in the forcing file 
                                                                       !! (unitless)
  REAL(r_std)                                :: dt_forcesoil           !! Time step at which soilcarbon routine is called (days) 
  INTEGER                                    :: nparan                 !! Number of values stored per year in the forcing file 
                                                                       !! (unitless)
  INTEGER                                    :: nbyear
  INTEGER(i_std),DIMENSION(:),ALLOCATABLE    :: indices                !! Grid Point Index used per processor (unitless)
  INTEGER(i_std),DIMENSION(:),ALLOCATABLE    :: indices_g              !! Grid Point Index for all processor (unitless)
  REAL(r_std),DIMENSION(:),ALLOCATABLE       :: x_indices_g            !! Grid Point Index for all processor (unitless)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE     :: lon, lat               !! Longitude and Latitude of each grid point defined 
                                                                       !! in lat/lon (2D) (degrees)
  REAL(r_std),DIMENSION(llm)                 :: lev                    !! Number of level (requested by restini routine) (unitless)


  INTEGER                                    :: i,m,iatt,iv,iyear  !! counters (unitless)

  CHARACTER(LEN=80)                          :: var_name
  CHARACTER(LEN=800)                         :: taboo_vars         !! string used for storing the name of the variables 
                                                                   !! of the stomate restart file that are not automatically 
                                                                   !! duplicated from input to output restart file (unitless)
  REAL(r_std),DIMENSION(1)                   :: xtmp               !! scalar read/written in restget/restput routines (unitless)
  INTEGER(i_std),PARAMETER                   :: nbvarmax=300       !! maximum # of variables assumed in the stomate restart file 
                                                                   !! (unitless)
  INTEGER(i_std)                             :: nbvar              !! # of variables effectively present 
                                                                   !! in the stomate restart file (unitless)
  CHARACTER(LEN=50),DIMENSION(nbvarmax)      :: varnames           !! list of the names of the variables stored 
                                                                   !! in the stomate restart file (unitless)
  INTEGER(i_std)                             :: varnbdim           !! # of dimensions of a given variable 
                                                                   !! of the stomate restart file
  INTEGER(i_std),PARAMETER                   :: varnbdim_max=20    !! maximal # of dimensions assumed for any variable 
                                                                   !! of the stomate restart file 
  INTEGER,DIMENSION(varnbdim_max)            :: vardims            !! length of each dimension of a given variable 
                                                                   !! of the stomate restart file
  LOGICAL                                    :: l1d                !! boolean : TRUE if all dimensions of a given variable 
                                                                   !! of the stomate restart file are of length 1 (ie scalar) 
                                                                   !! (unitless)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE     :: var_3d             !! matrix read/written in restget/restput routines (unitless)
  REAL(r_std)                                :: x_tmp              !! temporary variable used to store return value 
                                                                   !! from nf90_get_att (unitless)
  CHARACTER(LEN=10)  :: part_str                                   !! string suffix indicating the index of a PFT 
  REAL(r_std),ALLOCATABLE,SAVE,DIMENSION(:)  :: clay_g             !! clay fraction (nbpglo) (unitless)
  REAL(r_std),DIMENSION(:,:,:,:),ALLOCATABLE :: soilcarbon_input_g !! soil carbon input (nbpglob,ncarb,nvm,time) 
                                                                   !! (\f$gC m^{-2} dt_forcesoil^{-1}\f$) 
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: control_temp_g     !! Temperature control (nbp_glo,above/below,time) on OM decomposition 
                                                                   !! (unitless)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: control_moist_g    !! Moisture control (nbp_glo,abo/below,time) on OM decomposition 
                                                                   !! ?? Should be defined per PFT as well (unitless)
  REAL(r_std),DIMENSION(:,:,:),ALLOCATABLE   :: carbon_g           !! Soil carbon stocks (nbp_glo,ncarb,nvm) (\f$gC m^{-2}\f$)

  REAL(r_std),ALLOCATABLE :: clay(:)                   !! clay fraction (nbp_loc) (unitless)
  REAL(r_std),ALLOCATABLE :: soilcarbon_input(:,:,:,:) !! soil carbon input (nbp_loc,ncarb,nvm,time) 
                                                       !! (\f$gC m^{-2} dt_forcesoil^{-1}\f$) 
  REAL(r_std),ALLOCATABLE :: control_temp(:,:,:)       !! Temperature control (nbp_loc,above/below,time) on OM decomposition 
                                                       !! (unitless)
  REAL(r_std),ALLOCATABLE :: control_moist(:,:,:)      !! Moisture control (nbp_loc,abo/below,time) on OM decomposition 
                                                       !! ?? Should be defined per PFT as well (unitless)
  REAL(r_std),ALLOCATABLE :: carbon(:,:,:)             !! Soil carbon stocks (nbp_loc,ncarb,nvm) (\f$gC m^{-2}\f$)
  REAL(r_std),DIMENSION(:,:),ALLOCATABLE     :: resp_hetero_soil  !! Heterotrophic respiration (\f$gC m^{-2} dt_forcesoil^{-1}\f$) 
                                                                  !! (requested by soilcarbon routine but not used here) 

  INTEGER(i_std)                             :: printlev_loc      !! Local write level
  INTEGER(i_std)                             :: ier,iret          !! Used for hangling errors 

  CHARACTER(LEN=30) :: temp_name 
  LOGICAL :: l_error                                              !! boolean for memory allocation
  REAL(r_std), ALLOCATABLE, DIMENSION(:,:,:,:) :: matrixA         !! Carbon fluxes matrix 

!_ =================================================================================================================================
 
  CALL Init_orchidee_para
  CALL init_timer

! Set specific write level to forcesoil using PRINTLEV_forcesoil=[0-4] in run.def. 
! The global printlev is used as default value. 
  printlev_loc=get_printlev('forcesoil')

!-
! Configure the number of PFTS 
!-
  
  ! 1. Read the number of PFTs
  !
  !Config Key   = NVM
  !Config Desc  = number of PFTs  
  !Config If    = OK_SECHIBA or OK_STOMATE
  !Config Def   = 14
  !Config Help  = The number of vegetation types define by the user
  !Config Units = [-]
  CALL getin_p('NVM',nvm)

  ! 2. Allocation
  l_error = .FALSE.
  ALLOCATE(pft_to_mtc(nvm),stat=ier)
  l_error = l_error .OR. (ier .NE. 0)
  IF (l_error) THEN
     STOP 'pft_to_mtc (forcesoil only) : error in memory allocation'
  ENDIF

  ! 3. Initialisation of the correspondance table
  pft_to_mtc(:) = undef_int
  
  ! 4.Reading of the conrrespondance table in the .def file
  !
  !Config Key   = PFT_TO_MTC
  !Config Desc  = correspondance array linking a PFT to MTC
  !Config if    = OK_SECHIBA or OK_STOMATE
  !Config Def   = 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
  !Config Help  =
  !Config Units = [-]
  CALL getin_p('PFT_TO_MTC',pft_to_mtc)

  ! 4.1 if nothing is found, we use the standard configuration
  IF(nvm == nvmc ) THEN
     IF(pft_to_mtc(1) == undef_int) THEN
        WRITE(numout,*) 'Note to the user : we will use ORCHIDEE to its standard configuration'
        pft_to_mtc(:) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 /)
     ENDIF
  ELSE   
     IF(pft_to_mtc(1) == undef_int) THEN
        WRITE(numout,*)' The array PFT_TO_MTC is empty : we stop'
     ENDIF
  ENDIF
  
  ! 4.2 What happened if pft_to_mtc(j) > nvmc (if the mtc doesn't exist)?
  DO i = 1, nvm
     IF(pft_to_mtc(i) > nvmc) THEN
        WRITE(numout,*) "the MTC you chose doesn't exist"
        STOP 'we stop reading pft_to_mtc'
     ENDIF
  ENDDO
  
  ! 4.3 Check if pft_to_mtc(1) = 1 
  IF(pft_to_mtc(1) /= 1) THEN
     WRITE(numout,*) 'the first pft has to be the bare soil'
     STOP 'we stop reading next values of pft_to_mtc'
  ELSE
     DO i = 2,nvm
        IF(pft_to_mtc(i) == 1) THEN
           WRITE(numout,*) 'only pft_to_mtc(1) has to be the bare soil'
           STOP 'we stop reading pft_to_mtc'
        ENDIF
     ENDDO
  ENDIF
  
  ! 5. Allocate and initialize natural ans is_c4
  
  ! 5.1 Memory allocation
  l_error = .FALSE.
  ALLOCATE(natural(nvm),stat=ier)
  l_error = l_error .OR. (ier .NE. 0)
  ALLOCATE(is_c4(nvm),stat=ier)

  IF (l_error) THEN
     STOP 'natural or is_c4 (forcesoil only) : error in memory allocation'
  ENDIF

  ! 5.2 Initialisation
  DO i = 1, nvm
     natural(i) = natural_mtc(pft_to_mtc(i))
     is_c4(i) = is_c4_mtc(pft_to_mtc(i))
  ENDDO

  !!- 
  !! 1. Initialisation stage
  !! Reading a set of input files, allocating variables and preparing output restart file.     
  !!-
  ! Define restart file name
  ! for reading initial conditions (sto_restname_in var) and for writting final conditions (sto_restname_out var). 
  ! User values are used if present in the .def file.
  ! If not present, default values (stomate_start.nc and stomate_rest_out.c) are used.
  !-
  IF (is_root_prc) THEN
     sto_restname_in = 'stomate_start.nc'
     CALL getin ('STOMATE_RESTART_FILEIN',sto_restname_in)
     WRITE(numout,*) 'STOMATE INPUT RESTART_FILE: ',TRIM(sto_restname_in)
     sto_restname_out = 'stomate_rest_out.nc'
     CALL getin ('STOMATE_RESTART_FILEOUT',sto_restname_out)
     WRITE(numout,*) 'STOMATE OUTPUT RESTART_FILE: ',TRIM(sto_restname_out)
     !-
     ! Open the input file and Get some Dimension and Attributes ID's 
     !-
     iret = NF90_OPEN (sto_restname_in, NF90_NOWRITE, rest_id_sto)
     iret = NF90_INQUIRE_DIMENSION (rest_id_sto,1,len=iim_g)
     iret = NF90_INQUIRE_DIMENSION (rest_id_sto,2,len=jjm_g)
     iret = NF90_INQ_VARID (rest_id_sto, "time", iv)
     iret = NF90_GET_ATT (rest_id_sto, iv, 'calendar',thecalendar)
     iret = NF90_CLOSE (rest_id_sto)
     i=INDEX(thecalendar,ACHAR(0))
     IF ( i > 0 ) THEN
        thecalendar(i:20)=' '
     ENDIF
     !-
     ! Allocate longitudes and latitudes
     !-
     ALLOCATE (lon(iim_g,jjm_g))
     ALLOCATE (lat(iim_g,jjm_g))
     lon(:,:) = zero
     lat(:,:) = zero
     lev(1)   = zero
     !-
     CALL restini &
          & (sto_restname_in, iim_g, jjm_g, lon, lat, llm, lev, &
          &  sto_restname_out, itau_dep, date0, dt_files, rest_id_sto)
  ENDIF

  CALL bcast(date0)
  CALL bcast(thecalendar)
  WRITE(numout,*) "calendar = ",thecalendar
  !-
  ! calendar
  !-
  CALL ioconf_calendar (thecalendar)
  CALL ioget_calendar  (one_year,one_day)
  CALL ioconf_startdate(date0)
  !
  !! For master process only
  !
  IF (is_root_prc) THEN
     !-
     ! define forcing file's name (Cforcing_name var)
     ! User value is used if present in the .def file
     ! If not, default (NONE) is used
     !-
     Cforcing_name = 'NONE'
     CALL getin ('STOMATE_CFORCING_NAME',Cforcing_name)
     !-
     ! Open FORCESOIL's forcing file to read some basic info (dimensions, variable ID's)
     ! and allocate variables.
     !-
     iret = NF90_OPEN (TRIM(Cforcing_name),NF90_NOWRITE,Cforcing_id)
     IF (iret /= NF90_NOERR) THEN
        CALL ipslerr (3,'forcesoil', &
             &        'Could not open file : ', &
             &          Cforcing_name,'(Do you have forget it ?)')
     ENDIF
     !-
     ! Total Domain size is stored in nbp_glo variable
     !-
     ier = NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'kjpindex',x_tmp)
     nbp_glo = NINT(x_tmp)
     !-
     ! Number of values stored per year in the forcing file is stored in nparan var.
     !-
     ier = NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'nparan',x_tmp)
     nparan = NINT(x_tmp)
     ier = NF90_GET_ATT (Cforcing_id,NF90_GLOBAL,'nbyear',x_tmp)
     nbyear = NINT(x_tmp)
     !-
     ALLOCATE (indices_g(nbp_glo))
     ALLOCATE (clay_g(nbp_glo))
     !-
     ALLOCATE (x_indices_g(nbp_glo),stat=ier)
     ier = NF90_INQ_VARID (Cforcing_id,'index',v_id)
     ier = NF90_GET_VAR   (Cforcing_id,v_id,x_indices_g)
     indices_g(:) = NINT(x_indices_g(:))
     WRITE(numout,*) mpi_rank,"indices globaux : ",indices_g
     DEALLOCATE (x_indices_g)
     !-
     ier = NF90_INQ_VARID (Cforcing_id,'clay',v_id)
     ier = NF90_GET_VAR   (Cforcing_id,v_id,clay_g)
     !-
     ! time step of forcesoil program (in days)
     !-
     dt_forcesoil = one_year / FLOAT(nparan)
     ! Initialize dt_sechiba in module constantes_var. This is needed in carbonsoil.
     dt_sechiba=dt_forcesoil*one_day
     WRITE(numout,*) 'time step (d): ',dt_forcesoil
     WRITE(numout,*) 'nparan: ',nparan
     WRITE(numout,*) 'nbyear: ',nbyear    
     !-
     ! read and write the variables in the output restart file we do not modify within the Forcesoil program
     ! ie all variables stored in the input restart file except those stored in taboo_vars
     !-
     taboo_vars ='$lon$ $lat$ $lev$ $nav_lon$ $nav_lat$ $nav_lev$ $time$ $time_steps$ '// &
          &             '$day_counter$ $dt_days$ $date$ '
     !-
     DO m = 1,nvm
        WRITE(part_str,'(I2)') m
        IF (m < 10) part_str(1:1) = '0'
        temp_name = '$carbon_'//part_str(1:LEN_TRIM(part_str))//'$'
        taboo_vars = TRIM(taboo_vars)//' '//TRIM(temp_name)
     ENDDO
     !-
     CALL ioget_vname(rest_id_sto, nbvar, varnames)
     !-
     ! read and write some special variables (1D or variables that we need)
     !-
     var_name = 'day_counter'
     CALL restget (rest_id_sto, var_name, 1, 1, 1, itau_dep, .TRUE., xtmp)
     CALL restput (rest_id_sto, var_name, 1, 1, 1, itau_dep, xtmp)
     !-
     var_name = 'dt_days'
     CALL restget (rest_id_sto, var_name, 1, 1, 1, itau_dep, .TRUE., xtmp)
     CALL restput (rest_id_sto, var_name, 1, 1, 1, itau_dep, xtmp)
     !-
     var_name = 'date'
     CALL restget (rest_id_sto, var_name, 1, 1, 1, itau_dep, .TRUE., xtmp)
     CALL restput (rest_id_sto, var_name, 1, 1, 1, itau_dep, xtmp)
     !-
     DO iv=1,nbvar
        !-- check if the variable is to be written here
        IF (INDEX(taboo_vars,'$'//TRIM(varnames(iv))//'$') == 0 ) THEN
           !---- get variable dimensions, especially 3rd dimension
           CALL ioget_vdim &
                &      (rest_id_sto, varnames(iv), varnbdim_max, varnbdim, vardims)
           l1d = ALL(vardims(1:varnbdim) == 1)
           !---- read it
           IF (l1d) THEN
              CALL restget &
                   &        (rest_id_sto, TRIM(varnames(iv)), 1, vardims(3), &
                   &         1, itau_dep, .TRUE., xtmp)
           ELSE
              ALLOCATE( var_3d(nbp_glo,vardims(3)), stat=ier)
              IF (ier /= 0) STOP 'ALLOCATION PROBLEM'
              !----
              CALL restget &
                   &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                   &         1, itau_dep, .TRUE., var_3d, "gather", nbp_glo, indices_g)
           ENDIF
           !---- write it
           IF (l1d) THEN
              CALL restput &
                   &        (rest_id_sto, TRIM(varnames(iv)), 1, vardims(3), &
                   &         1, itau_dep, xtmp)
           ELSE
              CALL restput &
                   &        (rest_id_sto, TRIM(varnames(iv)), nbp_glo, vardims(3), &
                   &         1, itau_dep, var_3d, 'scatter',  nbp_glo, indices_g)
              !----
              DEALLOCATE(var_3d)
           ENDIF
        ENDIF
     ENDDO
     !-
     ! read soil carbon stocks values stored in the input restart file
     !-
     ALLOCATE(carbon_g(nbp_glo,ncarb,nvm))
     carbon_g(:,:,:) = val_exp
     DO m = 1, nvm
        WRITE (part_str, '(I2)') m
        IF (m<10) part_str(1:1)='0'
        var_name = 'carbon_'//part_str(1:LEN_TRIM(part_str))
        CALL restget &
             &    (rest_id_sto, var_name, nbp_glo, ncarb , 1, itau_dep, &
             &     .TRUE., carbon_g(:,:,m), 'gather', nbp_glo, indices_g)
        IF (ALL(carbon_g(:,:,m) == val_exp)) carbon_g(:,:,m) = zero
        !-- do not write this variable: it will be modified.
     ENDDO
     WRITE(numout,*) "date0 : ",date0, itau_dep
     !-
     ! Analytical spinup is set to false
     !
     spinup_analytic = .FALSE.

     ! Length of the run (in Years)
     ! User value is used if present in the .def file
     ! If not, default value (10000 Years) is used
     !-
     WRITE(time_str,'(a)') '10000Y'
     CALL getin('TIME_LENGTH', time_str)
     write(numout,*) 'Number of years for carbon spinup : ',time_str
     ! transform into itau
     CALL tlen2itau(time_str, dt_forcesoil*one_day, date0, itau_len)
     write(numout,*) 'Number of time steps to do: ',itau_len
     !-
     ! read soil carbon inputs, water and temperature stresses on OM decomposition 
     ! into the forcing file - We read an average year.
     !-
     ALLOCATE(soilcarbon_input_g(nbp_glo,ncarb,nvm,nparan*nbyear))
     ALLOCATE(control_temp_g(nbp_glo,nlevs,nparan*nbyear))
     ALLOCATE(control_moist_g(nbp_glo,nlevs,nparan*nbyear))
     !-
     ier = NF90_INQ_VARID (Cforcing_id,'soilcarbon_input',v_id)
     ier = NF90_GET_VAR   (Cforcing_id,v_id,soilcarbon_input_g)
     ier = NF90_INQ_VARID (Cforcing_id,   'control_moist',v_id)
     ier = NF90_GET_VAR   (Cforcing_id,v_id,control_moist_g)
     ier = NF90_INQ_VARID (Cforcing_id,    'control_temp',v_id)
     ier = NF90_GET_VAR   (Cforcing_id,v_id,control_temp_g)
     !-
     ier = NF90_CLOSE (Cforcing_id)
     !-
  ENDIF
  CALL bcast(nparan)
  CALL bcast(nbyear)
  CALL bcast(dt_forcesoil)
  CALL bcast(iim_g)
  CALL bcast(jjm_g)
  CALL bcast(nbp_glo)
  CALL bcast(itau_dep)
  CALL bcast(itau_len)
  IF (.NOT. ALLOCATED(indices_g)) ALLOCATE (indices_g(nbp_glo))
  CALL bcast(indices_g)
 
  !
  ! We must initialize data_para :
  CALL init_orchidee_data_para_driver(nbp_glo,indices_g)

  kjpindex=nbp_loc
  jjm=jj_nb
  iim=iim_g
  IF (printlev_loc>=3) WRITE(numout,*) "Local grid : ",kjpindex,iim,jjm

  !---
  !--- Create the index table
  !---
  !--- This job returns a LOCAL kindex.
  !---
  ALLOCATE (indices(kjpindex),stat=ier)
  !
  !! scattering to all processes in parallel mode
  !
  CALL scatter(indices_g,indices)
  indices(1:kjpindex)=indices(1:kjpindex)-(jj_begin-1)*iim_g
  IF (printlev_loc>=3) WRITE(numout,*) mpi_rank,"indices locaux = ",indices(1:kjpindex)
  !-
  ! Allocation of the variables for a processor
  !-
  ALLOCATE(clay(kjpindex))
  ALLOCATE(soilcarbon_input(kjpindex,ncarb,nvm,nparan*nbyear))
  ALLOCATE(control_temp(kjpindex,nlevs,nparan*nbyear))
  ALLOCATE(control_moist(kjpindex,nlevs,nparan*nbyear))
  ALLOCATE(carbon(kjpindex,ncarb,nvm))
  ALLOCATE(resp_hetero_soil(kjpindex,nvm))
  ALLOCATE(matrixA(kjpindex,nvm,nbpools,nbpools))
  DO i = 1,nbpools
     matrixA(:,:,i,i) = un
  ENDDO
  iatt = 0

  !-
  ! Initialization of the variables for a processor
  !-
  CALL Scatter(clay_g,clay)
  CALL Scatter(soilcarbon_input_g,soilcarbon_input)
  CALL Scatter(control_temp_g,control_temp)
  CALL Scatter(control_moist_g,control_moist)
  CALL Scatter(carbon_g,carbon)  

!-
! Configuration of the parameters
!-

  !Config Key   = FRAC_CARB_AP
  !Config Desc  = frac carb coefficients from active pool: depends on clay content
  !Config if    = OK_STOMATE 
  !Config Def   = 0.004
  !Config Help  = fraction of the active pool going to the passive pool
  !Config Units = [-]
  CALL getin_p('FRAC_CARB_AP',frac_carb_ap)  
  !
  !Config Key   = FRAC_CARB_SA
  !Config Desc  = frac_carb_coefficients from slow pool
  !Config if    = OK_STOMATE 
  !Config Def   = 0.42
  !Config Help  = fraction of the slow pool going to the active pool
  !Config Units = [-] 
  CALL getin_p('FRAC_CARB_SA',frac_carb_sa)
  !
  !Config Key   = FRAC_CARB_SP
  !Config Desc  = frac_carb_coefficients from slow pool
  !Config if    = OK_STOMATE 
  !Config Def   = 0.03
  !Config Help  = fraction of the slow pool going to the passive pool
  !Config Units = [-] 
  CALL getin_p('FRAC_CARB_SP',frac_carb_sp)
  !
  !Config Key   = FRAC_CARB_PA
  !Config Desc  = frac_carb_coefficients from passive pool
  !Config if    = OK_STOMATE 
  !Config Def   = 0.45
  !Config Help  = fraction of the passive pool going to the passive pool
  !Config Units = [-] 
  CALL getin_p('FRAC_CARB_PA',frac_carb_pa)
  !
  !Config Key   = FRAC_CARB_PS
  !Config Desc  = frac_carb_coefficients from passive pool
  !Config if    = OK_STOMATE 
  !Config Def   = 0.0
  !Config Help  = fraction of the passive pool going to the passive pool
  !Config Units = [-]
  CALL getin_p('FRAC_CARB_PS',frac_carb_ps)
  !
  !Config Key   = ACTIVE_TO_PASS_CLAY_FRAC
  !Config Desc  = 
  !Config if    = OK_STOMATE 
  !Config Def   =  .68  
  !Config Help  =
  !Config Units = [-]
  CALL getin_p('ACTIVE_TO_PASS_CLAY_FRAC',active_to_pass_clay_frac)
  !
  !Config Key   = CARBON_TAU_IACTIVE
  !Config Desc  = residence times in carbon pools
  !Config if    = OK_STOMATE 
  !Config Def   = 0.149
  !Config Help  =
  !Config Units = [days] 
  CALL getin_p('CARBON_TAU_IACTIVE',carbon_tau_iactive)
  !
  !Config Key   = CARBON_TAU_ISLOW
  !Config Desc  = residence times in carbon pools
  !Config if    = OK_STOMATE 
  !Config Def   = 5.48
  !Config Help  =
  !Config Units = [days]
  CALL getin_p('CARBON_TAU_ISLOW',carbon_tau_islow)
  !
  !Config Key   = CARBON_TAU_IPASSIVE
  !Config Desc  = residence times in carbon pools
  !Config if    = OK_STOMATE 
  !Config Def   = 241.
  !Config Help  =
  !Config Units = [days]
  CALL getin_p('CARBON_TAU_IPASSIVE',carbon_tau_ipassive)
  !
  !Config Key   = FLUX_TOT_COEFF
  !Config Desc  =
  !Config if    = OK_STOMATE 
  !Config Def   = 1.2, 1.4,.75
  !Config Help  =
  !Config Units = [days]
  CALL getin_p('FLUX_TOT_COEFF',flux_tot_coeff)

  !!-
  !! 2. Computational step
  !! Loop over time - Call of soilcarbon routine at each time step 
  !! Updated soil carbon stocks are stored into carbon variable
  !! We only keep the last value of carbon variable (no time dimension).
  !!-
  iyear=1
  DO i=1,itau_len
     iatt = iatt+1
     IF (iatt > nparan*nbyear) THEN
        IF (printlev_loc>=3) WRITE(numout,*) iyear
        iatt = 1
        iyear=iyear+1
     ENDIF
     CALL soilcarbon &
          &    (kjpindex, clay, &
          &     soilcarbon_input(:,:,:,iatt), &
          &     control_temp(:,:,iatt), control_moist(:,:,iatt), &
          &     carbon, resp_hetero_soil, &
          &     matrixA)
  ENDDO
  WRITE(numout,*) "End of soilcarbon LOOP."

  !
  !! Gathering of variables towards main processor in parallel mode
  !
  CALL Gather(carbon,carbon_g)
  !!-
  !! 3. write new carbon stocks into the ouput restart file
  !!-
  IF (is_root_prc) THEN
     DO m=1,nvm
        WRITE (part_str, '(I2)') m
        IF (m<10) part_str(1:1)='0'
        var_name = 'carbon_'//part_str(1:LEN_TRIM(part_str))
        CALL restput &
             &    (rest_id_sto, var_name, nbp_glo, ncarb , 1, itau_dep, &
             &     carbon_g(:,:,m), 'scatter', nbp_glo, indices_g)
     ENDDO
     !-
     CALL getin_dump
     CALL restclo
  ENDIF
#ifdef CPP_PARA
  CALL MPI_FINALIZE(ier)
#endif
  WRITE(numout,*) "End of forcesoil."
  !--------------------
END PROGRAM forcesoil
