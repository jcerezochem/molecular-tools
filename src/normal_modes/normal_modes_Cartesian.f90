program normal_modes_cartesian


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS 
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to visualize vibrations obtained in internal coordinates.
    !
    !============================================================================    

    !*****************
    !   MODULE LOAD
    !*****************
    !============================================
    !   Generic
    !============================================
    use alerts
    use line_preprocess
    use constants 
    use verbosity
    use matrix
    use matrix_print
    use io !printing
    !============================================
    !   Structure types module
    !============================================
    use structure_types
    !============================================
    !   File readers
    !============================================
    use generic_io
    use generic_io_molec
    use xyz_manage_molec
    use gro_manage
    use gaussian_manage
    !============================================
    !  Structure-related modules
    !============================================
    use molecular_structure
    use atomic_geom
    use symmetry
    !============================================
    !  Vibrational 
    !============================================
    use internal_module
    use zmat_manage 
    use vibrational_analysis
    use vertical_model
    use thermochemistry
    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: use_symmetry=.false.,    &
               include_hbonds=.false.,  &
               vertical=.false.,        &
               analytic_Bder=.false.,   &
               check_symmetry=.true.,   &
               full_diagonalize=.false.,&
               animate=.true.,          &
               orthogonalize=.false.,   &
               Eckart_frame=.true.,     &
               modes_as_internals=.false., &
               original_internal=.false.,  &
               rm_gradcoord=.false.,       &
               apply_projection_matrix=.false., &
               print_modes=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule, molec_aux
    type(str_bonded) :: zmatgeom
    real(8),dimension(NDIM) :: X0, Y0, Z0
    integer,dimension(1:NDIM) :: isym
    integer,dimension(4,1:NDIM,1:NDIM) :: Osym
    integer :: nsym
    integer :: Nat, Nvib0, Nvib, Ns, Nrt, Nf
    character(len=5) :: PG
    real(8) :: Tthermo=0.d0
    !Job info
    character(len=20) :: calc, method, basis, method_
    character(len=150):: title
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=50) :: dummy_char
    real(8) :: dist
    !io flags
    integer :: error, info
    ! Auxiliar real arrays/scalars
    real(8),dimension(1:NDIM,1:NDIM) :: Aux, Aux2
    real(8),dimension(1:NDIM)        :: Vec, Vec1
    real(8),dimension(:),allocatable :: Vec_alloc
    real(8) :: Theta, Theta2
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii, jj, kk, iop
    !=============

    !====================== 
    ! PES topology and normal mode things
    real(8),dimension(1:NDIM,1:NDIM) :: LL, D, P, Lcartinv, Lmwc
    real(8),dimension(1:NDIM*NDIM)   :: Hlt
    real(8),dimension(:,:),allocatable :: Hess
    real(8),dimension(NDIM) :: Freq, Factor, Grad, Grad1
    !Moving normal modes
    character(len=50) :: selection="none"
    real(8) :: Amplitude = 2.d0, qcoord
    integer,dimension(1:NDIM) :: nm=0
    real(8) :: Qstep
    logical :: call_vmd = .false.
    character(len=10000) :: vmdcall
    integer :: Nsteps, Nsel=0, istep
    !MOVIE things
    logical :: movie_vmd = .false.
    integer :: movie_cycles=0,& !this means no movie
               movie_steps
    !====================== 

    !====================== 
    !INTERNAL CODE THINGS
    real(8),dimension(1:NDIM,1:NDIM) :: B, G, Asel, Aselinv
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    real(8),dimension(1:NDIM,1:NDIM) :: X,Xinv
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: S, Sref, S0
    integer,dimension(NDIM) :: S_sym
    ! Switches
    character(len=5) :: def_internal="ALL", def_internal_aux
    character(len=2) :: scan_type="NM"
    !Coordinate map
    integer,dimension(NDIM) :: Zmap
    ! Number of ic (Shortcuts)
    integer :: nbonds, nangles, ndihed, nimprop
    !====================== 

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_SYM=12,  &
               I_RMF=16,  &
               I_CNX=17,  &
               I_RMC=18,  &
               I_MAS=19,  &
               O_GRO=20,  &
               O_G09=21,  &
               O_G96=22,  &
               O_Q  =23,  &
               O_NUM=24,  &
               O_NUMD=27, &
               O_MOV=25,  &
               O_FCHK=26, &
               O_ARR=27,  &
               S_VMD=30

    !files
    character(len=10) :: ft ="guess",  ftg="guess",  fth="guess", ftn="guess"
    character(len=200):: inpfile  ="state1.fchk", &
                         gradfile ="same", &
                         hessfile ="same", &
                         nmfile   ="none", &
                         intfile  ="none", &
                         rmzfile  ="none", &
                         symm_file="none", &
                         cnx_file ="guess",&
                         mass_file="none", &
                         rm_custom_coord="none", &
                         rm_custom_mode ="none", &
                         outfchkfile='none'
    !Structure files to be created
    character(len=100) :: g09file,qfile, tmpfile, g96file, grofile,numfile,numfwfile,numbwfile,arrowfile
    !status
    integer :: IOstatus
    !===================

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
    integer :: lenght, N
    !====================== 

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

    call cpu_time(ti)

    !--------------------------
    ! Tune io
    !--------------------------
    ! Set unit for alert messages
    alert_unt=6
    ! Activate notes
    silent_notes = .false.
    !--------------------------

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(molecule)
    call allocate_atoms(molec_aux)
    !===========================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(&
                     ! input data
                     inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn,mass_file,&
                     ! Options (general)
                     Amplitude,call_vmd,include_hbonds,selection,vertical, &
                     outfchkfile,                                          &
                     ! Options (Cartesian)
                     full_diagonalize,                                     &
                     ! Remove coordinate along gradient
                     rm_gradcoord,rm_custom_coord,rm_custom_mode,          &
                     ! Animation and Movie
                     animate,movie_vmd, movie_cycles,print_modes,          &   
                     ! Options (internal)
                     use_symmetry,def_internal,intfile,rmzfile,            & !except scan_type
                     ! Additional vib options
                     Eckart_frame,orthogonalize,modes_as_internals,        &
                     original_internal, apply_projection_matrix,           &
                     ! connectivity file
                     cnx_file,                                             &
                     ! thermochemical analysis
                     Tthermo,                                              &
                     ! (hidden)
                     analytic_Bder)
    call set_word_upper_case(def_internal)

 
    ! ---------------------------------
    ! 1. READ DATA
    call heading(6,"INPUT DATA")
    ! ---------------------------------
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)
    if (fth == "guess") &
    call split_line_back(hessfile,".",null,fth)
    if (ftg == "guess") &
    call split_line_back(gradfile,".",null,ftg)
    if (ftn == "guess") &
    call split_line_back(nmfile,".",null,ftn)

    ! Manage special files (fcc) 
    if (adjustl(ft) == "fcc-state" .or. adjustl(ftn) == "fcc-state") then
        call alert_msg("note","fcc-state files needs fcc-input as -f and statefile as -fnm")
        ft ="fcc-state"
        ftn="fcc-state"
        ! inpfile has Nat, Nvib, and Masses          <= in inpfile
        ! statefile has coordinates and normal modes <= in hessfile
        ! Generic generic readers parse the state (not the inpfile)
        ! so we get the info here
        open(I_INP,file=inpfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
        read(I_INP,*) Nat 
        molecule%natoms = Nat
        read(I_INP,*) Nvib
        do i=1,Nat 
            read(I_INP,*) molecule%atom(i)%mass
            !Set atomnames from atommasses
            call atominfo_from_atmass(molecule%atom(i)%mass,  &
                                      molecule%atom(i)%AtNum, &
                                      molecule%atom(i)%name)
        enddo
        close(I_INP)
        ! Now put the statefile in the inpfile
        inpfile=nmfile
    elseif (adjustl(ftn) == "log") then
        ! Need to read the standard orientation, not from summary section
        ft = "log-stdori"
    endif
        
    ! STRUCTURE FILE
    call statement(6,"READING MOLECULE FILE (STRUCTURE)...")
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,molecule,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (State1)")
    ! Get job info if it is a Gaussian file
    if (ft == "log" .or. ft== "fchk") then
        rewind(I_INP) ! this should be a generic_job_rewind() call
        call read_gauss_job(I_INP,ft,calc,method,basis)
        ! Whichever, calc type was, se now need SP
        calc="SP"
    else
        calc="SP"
        method="B3LYP"
        basis="6-31G(d)"
    endif
    close(I_INP)
    ! Shortcuts
    Nat = molecule%natoms

    ! Read mass from file if given
    if (adjustl(mass_file) /= "none") then
        print'(/,X,A)', "Reading atomic masses from: "//trim(adjustl(mass_file))
        open(I_MAS,file=mass_file,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(mass_file)) )
        do i=1,Nat
            read(I_MAS,*,iostat=IOstatus) molecule%atom(i)%mass 
            if (IOstatus /= 0) call alert_msg( "fatal","While reading "//trim(adjustl(mass_file)) )
        enddo
        close(I_MAS)
    endif


    ! VIBRATIONAL ANALYSIS: either read from file or from diagonalization of Hessian
    if (adjustl(nmfile) /= "none") then
        call statement(6,"READING NORMAL MODES FROM FILE...")
        open(I_INP,file=nmfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(nmfile)) )
        call generic_nm_reader(I_INP,ftn,Nat,Nvib,Freq,LL)
        ! Show frequencies
        if (verbose>0) &
         call print_vector(6,Freq,Nvib,"Frequencies (cm-1)")
        ! The reader provide L in Normalized Cartesian. Need to Transform to Cartesian now
        call LcartNrm_to_Lmwc(Nat,Nvib,molecule%atom(:)%mass,LL,LL)
        call Lmwc_to_Lcart(Nat,Nvib,molecule%atom(:)%mass,LL,LL,error)
        close(I_INP)

    else
        ! ACTUALLY PERFORM THE ANALYSIS
        ! HESSIAN FILE
        call statement(6,"READING HESSIAN FILE...")
        open(I_INP,file=hessfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile)) )
        call generic_Hessian_reader(I_INP,fth,Nat,Hlt,error)
            i=3*Nat*(3*Nat+1)/2
        if (error /= 0) call alert_msg("fatal","Error reading Hessian (State1)")
        close(I_INP)
        
        ! GRADIENT FILE
        if (adjustl(gradfile) /= "none") then
            call statement(6,"READING GRADIENT FILE...")
            open(I_INP,file=gradfile,status='old',iostat=IOstatus)
            if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
            call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
            close(I_INP)
            if (error /= 0) then
                print*, "Error reading the Gradient. It will be set to zero"
                Grad(1:3*Nat) = 0.d0
            endif
        else
            Grad(1:3*Nat) = 0.d0
        endif


        call heading(6,"PREPARING VIBRATIONAL ANALYSIS")
        if (vertical) then
            call statement(6,"At non-stationary point.")
        else
            call statement(6,"At stationary point.")
        endif

        ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
        call subheading(6,"Preliminary vibrational analysis (Cartesian)")
        if (rm_gradcoord) then
            call statement(6,"Vibrations on the 3N-7 space")
            call vibrations_Cart(Nat,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,molecule%atom(:)%Mass,&
                             Hlt,Nvib,LL,Freq,error_flag=error,Dout=D,Grad=Grad)
            ! Store the number of vibrational degrees on freedom on Nvib0
            ! Nvib stores the reduced dimensionality
            Nvib0=Nvib+1
        elseif (adjustl(rm_custom_coord) /= "none") then
            call subheading(6,"Vibrations on the 3N-7 space",upper_case=.true.)
            call subheading(6,"Vibrational analysis removing one custom coordinate")
            ! Read the custom coordinate. Store in Grad1
            open(I_RMC,file=rm_custom_coord,status="old")
            do i=1,3*Nat
                ii = (i-1)/3+1
                read(I_RMC,*) Grad1(i)
                Grad1(i) = Grad1(i)*molecule%atom(ii)%mass
            enddo
            ! 
            call vibrations_Cart(Nat,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,molecule%atom(:)%Mass,&
                             Hlt,Nvib,LL,Freq,error_flag=error,Dout=D,Grad=Grad1)
            ! Store the number of vibrational degrees on freedom on Nvib0
            ! Nvib stores the reduced dimensionality
            Nvib0 = Nvib+1
            rm_gradcoord=.true.
        elseif (adjustl(rm_custom_mode) /= "none") then
            call subheading(6,"Vibrations on the 3N-7 space",upper_case=.true.)
            call subheading(6,"Vibrational analysis removing one custom mode")
            ! Read the custom coordinate. Store in Grad1
            open(I_RMC,file=rm_custom_mode,status="old")
            do i=1,3*Nat
                ii = (i-1)/3+1
                read(I_RMC,*) Grad1(i)
!                 Grad1(i) = Grad1(i)*molecule%atom(ii)%mass
            enddo
            ! 
            call vibrations_Cart(Nat,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,molecule%atom(:)%Mass,&
                             Hlt,Nvib,LL,Freq,error_flag=error,Dout=D,Grad=Grad1)
            ! Store the number of vibrational degrees on freedom on Nvib0
            ! Nvib stores the reduced dimensionality
            Nvib0 = Nvib+1
            rm_gradcoord=.true.
        else
            call statement(6,"Vibrations on the 3N-6 space")
            call vibrations_Cart(Nat,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,molecule%atom(:)%Mass,Hlt,&
                             Nvib,LL,Freq,error_flag=error,Dout=D)
            Nvib0=Nvib
        endif

        if (vertical.or.apply_projection_matrix) then
            call statement(6,"CORRECTIONS FOR NON-STATIONARY POINTS ACTIVATED",keep_case=.true.)
       
            ! MANAGE INTERNAL COORDS
            ! --------------------------------
            call subheading(6,"GETTING CONNECTIVITY")
            ! Get connectivity 
            if (cnx_file == "guess") then
                call statement(6,"Guess connectivity based on distance criteria")
                call guess_connect(molecule)
            else
                call statement(6,"Reading connectivity from file: "//trim(adjustl(cnx_file)))
                open(I_CNX,file=cnx_file,status='old')
                call read_connect(I_CNX,molecule)
                close(I_CNX)
            endif

            call subheading(6,"Managing internal coordinates")
            ! Manage symmetry
            if (.not.use_symmetry) then
                molecule%PG="C1"
            else if (trim(adjustl(symm_file)) /= "none") then
                call alert_msg("note","Using custom symmetry file: "//trim(adjustl(symm_file)) )
                open(I_SYM,file=symm_file)
                do i=1,molecule%natoms
                    read(I_SYM,*) j, isym(j)
                enddo
                close(I_SYM)
                !Set PG to CUStom
                molecule%PG="CUS"
            else
                molecule%PG="XX"
                call symm_atoms(molecule,isym)
            endif
            
            !Generate bonded info
            call gen_bonded(molecule)
            
            !From now on work in au
            call set_geom_units(molecule,"Bohr")

            ! Define internal set
            ! Get the set of internal coordinates
            call subheading(6,"Generating internal set for analysis")
            call define_internal_set(molecule,def_internal,intfile,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Aux2)
            !From now on, we'll use atomic units
!             call set_geom_units(molecule,"Bohr")
            if (Ns > Nvib0) then
                call internals_mapping(molecule%geom,zmatgeom,Zmap)
            elseif (def_internal=="ZMAT".and.rmzfile/="none") then
                ! We also get a Zmap
                call internals_mapping(molecule%geom,zmatgeom,Zmap)
                Nvib0=Ns
!             elseif (Ns < Nvib) then
                print*, "Ns", Ns
                print*, "Nvib0", Nvib0
                call alert_msg("warning","Reduced coordinates only produce animations with rmzfiles")
                ! Need to freeze unused coords to its input values
            endif

            ! If we use uncorrected modes as linear combinations for internals
            ! we now compute them. 
            ! This means using linearised internals (should be equivalent to Cart)
            if (modes_as_internals) then
                ! Get Hess matrix from H(lowertriangular)
                allocate (Hess(1:NDIM,1:NDIM))
                Hess(1:3*Nat,1:3*Nat) = Hlt_to_Hess(3*Nat,Hlt)
                call verbose_mute()
                call internal_Wilson(molecule,Ns,S,B,ModeDef)
                call internal_Gmetric(Nat,Ns,molecule%atom(:)%mass,B,G)
                call redundant2nonredundant(Ns,Nvib0,G,Asel)
                call verbose_continue()
                G(1:Nvib0,1:Nvib0) = matrix_basisrot(Nvib0,Ns,Asel(1:Ns,1:Nvib0),G,&
                                                   counter=.true.)
                B(1:Nvib0,1:3*Nat) = matrix_product(Nvib0,3*Nat,Ns,Aselinv,B)
                call HessianCart2int(Nat,Nvib0,Hess,molecule%atom(:)%mass,B,G)
                call gf_method(Nvib0,Nvib0,G,Hess,LL,Freq,X,Xinv)
                ! Save L(in LL) and Linv(in Aux)
                Aux(1:Nvib0,1:Nvib0)= inverse_realgen(Nvib0,LL)
                Aux(1:Nvib0,1:Ns)   = matrix_product(Nvib0,Ns,Nvib0,Aux,Asel,tB=.true.)
                LL(1:Ns,1:Nvib0)    = matrix_product(Ns,Nvib0,Nvib0,Asel,LL)
                deallocate(Hess)
            endif

            !------------------------------------------------
            ! Get B, G and Bders to perform the correction
            !------------------------------------------------
            print'(X,A,/)', "Prepare correction for non-stationary points"
            call internal_Wilson(molecule,Ns,S,B,ModeDef)
            call internal_Gmetric(Nat,Ns,molecule%atom(:)%mass,B,G)
            call calc_Bder(molecule,Ns,Bder,analytic_Bder)
            ! Select internal linear combinations
            call subheading(6,"Construct linear combination of original internal coordinates")
            if (modes_as_internals) then
                print*, "Using internal defined as uncorrected modes"
                Asel(1:Ns,1:Nvib0) = LL(1:Ns,1:Nvib0)
                Aselinv(1:Nvib0,1:Ns) = Aux(1:Nvib0,1:Ns)
            elseif (original_internal.and.Nvib0==Ns) then
                print*, "Using internal without linear combination"
                Asel(1:Ns,1:Nvib0) = 0.d0
                do i=1,Nvib0
                    Asel(i,i) = 1.d0 
                enddo
                Aselinv(1:Nvib0,1:Ns) = Asel(1:Ns,1:Nvib0)
            else
                call statement(6,"Getting internals from eigevector of G")
                call redundant2nonredundant(Ns,Nvib0,G,Asel)
                Aselinv(1:Nvib0,1:Ns) = transpose(Asel(1:Ns,1:Nvib0))
            endif
            ! Rotate Gmatrix
            G(1:Nvib0,1:Nvib0) = matrix_basisrot(Nvib0,Ns,Aselinv(1:Nvib0,1:Ns),G)
            ! Check if we want orthogonalization
            if (orthogonalize) then
                print*, "Orthogonalyzing internals..."
                Xinv(1:Nvib0,1:Nvib0) = 0.d0
                X(1:Nvib0,1:Nvib0)     = 0.d0
                do i=1,Nvib0
                    Xinv(i,i) = 1.d0/dsqrt(G(i,i))
                    X(i,i)     = dsqrt(G(i,i))
                enddo
                ! Rotate Gmatrix (again)
                G(1:Nvib0,1:Nvib0) = matrix_basisrot(Nvib0,Nvib0,Xinv(1:Nvib0,1:Nvib0),G)
                ! Update Asel(inv)
                Aselinv(1:Nvib0,1:Ns) = matrix_product(Nvib0,Ns,Nvib0,Xinv,Aselinv)
                Asel(1:Ns,1:Nvib0)    = matrix_product(Ns,Nvib0,Nvib0,Asel,X)
            endif
            ! Rotate Bmatrix
            B(1:Nvib0,1:3*Nat) = matrix_product(Nvib0,3*Nat,Ns,Aselinv,B)
            ! Rotate Bders
            do j=1,3*Nat
                Bder(1:Nvib0,j,1:3*Nat) =  matrix_product(Nvib0,3*Nat,Ns,Aselinv,Bder(1:Ns,j,1:3*Nat))
            enddo

            if (apply_projection_matrix) then
                ! Get projection matrix
                P(1:3*Nat,1:3*Nat) = projection_matrix(Nat,Nvib0,B,molecule%atom(1:Nat)%Mass)
                Aux(1:3*Nat,1:3*Nat) = matrix_product(3*Nat,3*Nat,3*Nat,P,P,tA=.true.)
                call MAT0(6,Aux,3*Nat,3*Nat,"Check prj (1.A)")
                Aux(1:3*Nat,1:3*Nat) = matrix_product(3*Nat,3*Nat,3*Nat,P,P,tB=.true.)
                call MAT0(6,Aux,3*Nat,3*Nat,"Check prj (1.B)")
                call MAT0(6,P,3*Nat,3*Nat,"P matrix (1)")

! print*, "New projection"
!                 P(1:3*Nat,1:3*Nat) = projection_matrix2(Nat,molecule%atom(1:Nat)%X, &
!                                                             molecule%atom(1:Nat)%Y, &
!                                                             molecule%atom(1:Nat)%Z, &
!                                                             molecule%atom(1:Nat)%Mass)
! 
!                 Aux(1:3*Nat,1:3*Nat) = matrix_product(3*Nat,3*Nat,3*Nat,P,P,tA=.true.)
!                 call MAT0(6,Aux,3*Nat,3*Nat,"Check prj (2.A)")
!                 Aux(1:3*Nat,1:3*Nat) = matrix_product(3*Nat,3*Nat,3*Nat,P,P,tB=.true.)
!                 call MAT0(6,Aux,3*Nat,3*Nat,"Check prj (2.B)")

!                 call MAT0(6,P,3*Nat,3*Nat,"P matrix (2)")
                ! And rotate gradient
                do i=1,3*Nat
                    Vec1(i) = 0.d0
                    do k=1,Nvib0
                        Vec1(i) = Vec1(i) + P(i,k)*Grad(k)
                    enddo
                enddo
                Grad(1:3*Nat) = Vec1(1:3*Nat)
            endif

            !Reset Angs
            call set_geom_units(molecule,"Angs")

        endif


        ! VIBRATIONAL ANALYSIS
        call heading(6,"FINAL VIBRATIONAL ANALYSIS",print_always=.true.)
        !-------------------------------------
        ! Vibrational analysis:
        !-------------------------------------
        ! With grad correction (vertical)
        !   Ht = D M^-1/2 B^t Hs B M^-1/2 D^t  (assuming G constant)
        ! Without grad correction: 
        !   Ht = D M^1/2 Hx M^1/2 D^t
        !-------------------------------------
        ! Get Hess matrix from H(lowertriangular)
        allocate (Hess(1:NDIM,1:NDIM))
        Hess(1:3*Nat,1:3*Nat) = Hlt_to_Hess(3*Nat,Hlt)
        if (vertical) then
            ! (Hess is already constructed)
            ! Hs (with the correction)
            ! First get: Hx' = Hx - gs^t\beta
            ! 1. Get gs from gx
            Vec(1:3*Nat) = Grad(1:3*Nat)
            call Gradcart2int(Nat,Nvib0,Vec,molecule%atom(:)%mass,B,G)
            ! 2. Multiply gs^t\beta and
            ! 3. Apply the correction
            ! Bder(i,j,K)^t * gq(K)
            do i=1,3*Nat
            do j=1,3*Nat
                Aux2(i,j) = 0.d0
                do k=1,Nvib0
                    Aux2(i,j) = Aux2(i,j) + Bder(k,i,j)*Vec(k)
                enddo
                Hess(i,j) = Hess(i,j) - Aux2(i,j)
            enddo
            enddo
            if (verbose>2) then
                print*, "Correction matrix to be applied on Hx:"
                call MAT0(6,Aux2,3*Nat,3*Nat,"gs*Bder matrix")
            endif
            
            if (check_symmetry) then
                print'(/,X,A)', "---------------------------------------"
                print'(X,A  )', " Check effect of symmetry operations"
                print'(X,A  )', " on the correction term gs^t\beta"
                print'(X,A  )', "---------------------------------------"
                molecule%PG="XX"
                call symm_atoms(molecule,isym,Osym,rotate=.false.,nsym_ops=nsym)
                ! Check the symmetry of the correction term
                ! Check all detected symmetry ops
                do iop=1,Nsym
                    Aux(1:3*Nat,1:3*Nat) = dfloat(Osym(iop,1:3*Nat,1:3*Nat))
                    Aux(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,Aux,Aux2,counter=.true.)
                    Theta=0.d0
                    do i=1,3*Nat 
                    do j=1,3*Nat 
                        if (Theta < abs(Aux(i,j)-Aux2(i,j))) then
                            Theta = abs(Aux(i,j)-Aux2(i,j))
                            Theta2=Aux2(i,j)
                        endif
                    enddo
                    enddo
                    print'(X,A,I0)', "Symmetry operation :   ", iop
                    print'(X,A,F10.6)',   " Max abs difference : ", Theta
                    print'(X,A,F10.6,/)', " Value before sym op: ", Theta2
                enddo
                print'(X,A,/)', "---------------------------------------"
            endif
            ! Get Hs
            call HessianCart2int(Nat,Nvib0,Hess,molecule%atom(:)%mass,B,G)
            ! B^t Hs B [~Hx]
            Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,Nvib0,B,Hess,counter=.true.)
            
            if (adjustl(outfchkfile) /= 'none') then
                
                !================
                !REWRITE FCHK
                !================
                
                print'(X,A)', "WRITTING FCHK..."
                print'(X,A)', "  Output file: "//trim(adjustl(outfchkfile))
                open(O_FCHK,file=outfchkfile,status="unknown",iostat=IOstatus)
!                 if (overwrite) then
!                     open(O_FCHK,file=outfchkfile,status="replace",iostat=IOstatus)
!                 else
!                     open(O_FCHK,file=outfchkfile,status="new",iostat=IOstatus)
!                     if (IOstatus /= 0) &
!                      call alert_msg("fatal","Cannot open output for writting. Use -ow to force overwrite")
!                 endif
                ! Title and job info
                write(O_FCHK,'(A)') "FCHK created with reform_fchk from "//trim(adjustl(inpfile))
!                 write(O_FCHK,'(A10,A60,A10)') adjustl(calc_type), adjustl(method), adjustl(basis)
                call write_fchk(O_FCHK,"Number of atoms","I",0,A,(/Nat/),error)
!                 call write_fchk(O_FCHK,"Charge","I",0,A,(/charge/),error)
!                 call write_fchk(O_FCHK,"Multiplicity","I",0,A,(/mult/),error)
                
                !Atomic Numbers and Nuclear charges
                N=molecule%natoms
                call write_fchk(O_FCHK,"Atomic numbers",'I',N,A,molecule%atom(1:N)%AtNum,error)
                N=molecule%natoms
                allocate(IA(1:1),A(1:N))
                A(1:N) = float(molecule%atom(1:N)%AtNum)
                call write_fchk(O_FCHK,"Nuclear charges",'R',N,A,IA,error)
                deallocate(A,IA)
                !Coordinates
                N=3*molecule%natoms
                allocate(IA(1:1),A(1:N))
                do i=1,N/3
                    j=3*i
                    A(j-2) = molecule%atom(i)%x/BOHRtoANGS
                    A(j-1) = molecule%atom(i)%y/BOHRtoANGS
                    A(j)   = molecule%atom(i)%z/BOHRtoANGS
                enddo
                call write_fchk(O_FCHK,"Current cartesian coordinates",'R',N,A,IA,error)
                deallocate(A,IA)
                !Atomic weights
                call write_fchk(O_FCHK,"Integer atomic weights",'I',3*Nat,A,int(molecule%atom(:)%mass),error)
                call write_fchk(O_FCHK,"Real atomic weights",'R',3*Nat,molecule%atom(:)%mass,IA,error)
!                 !Energy 
!                 if (have_SCF) &
!                     call write_fchk(O_FCHK,"SCF Energy",'R',0,(/E_scf/),IA,error)
!                 if (have_TD) &
!                     call write_fchk(O_FCHK,"CIS Energy",'R',0,(/E_td/),IA,error)
!                 if (have_TOT) &
!                     call write_fchk(O_FCHK,"Total Energy",'R',0,(/E_tot/),IA,error)
!                 !Gradient
!                 if (have_gradient) then
!                     call write_fchk(O_FCHK,"Cartesian Gradient",'R',3*Nat,Grad,IA,error)
!                 endif
                !Hessian
!                 if (have_hessian) then
                    ! Transform the Hess into Hlt 
                    N=3*Nat*(3*Nat+1)/2
                    Hlt(1:N) = Hess_to_Hlt(3*Nat,Hess)
                    call write_fchk(O_FCHK,"Cartesian Force Constants",'R',N,Hlt,IA,error)
!                 endif
                
                print'(X,A,/)', "Done"
                
                
                
                
                
            endif
        endif
!         if (apply_projection_matrix) then
!             ! Project out rotation and translation
!             Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,P,Hess)
!         endif

        ! M^-1/2 [Hx] M^-1/2
        do i=1,3*Nat
        do j=1,i
            ii = (i-1)/3+1
            jj = (j-1)/3+1
            Hess(i,j) = Hess(i,j) &
                        /dsqrt(molecule%atom(ii)%mass*AMUtoAU) &
                        /dsqrt(molecule%atom(jj)%mass*AMUtoAU)
            Hess(j,i) = Hess(i,j)
        enddo
        enddo

        if (apply_projection_matrix) then
            Nrt = 3*Nat-Nvib
            P(1:3*Nat,1:3*Nat) = matrix_product(3*Nat,3*Nat,Nrt,D(1:3*Nat,1:Nrt),D(1:3*Nat,1:Nrt),tB=.true.)
            Aux(1:3*Nat,1:3*Nat) = identity_matrix(3*Nat)
            P(1:3*Nat,1:3*Nat) = Aux(1:3*Nat,1:3*Nat) - P(1:3*Nat,1:3*Nat)
            ! Project out rotation and translation
            Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,P,Hess)
        endif

        ! Eckart traslation and rotation can be computed for testing
        if (.not.Eckart_frame.and..not.full_diagonalize) then
            call alert_msg("note","With -noEckart -fulldiag is always done")
            full_diagonalize=.true.
        endif

        ! Get number of Trans+Rot 
        Nrt = 3*Nat - Nvib

        ! Get normal modes: 
        !   1.- T, R and Vib (separated Eckart frames)
        !   2.- T+R+Vib (same frame)
        !   3.- Vib (separated Eckart frame)
        if (full_diagonalize.and.Eckart_frame) then
            LL(1:3*Nat,1:3*Nat) = 0.d0
            call statement(6,"Getting T, R and Vib in separated Eckart frames")
            ! 1a) Get T modes fromt he 3x3 block in the Eckart frame
            ! D [M^-1/2 [Hx] M^1/2] D^t
            Nrt = 3
            Aux(1:Nrt,1:Nrt) = matrix_basisrot(Nrt,3*Nat,D(1:3*Nat,1:Nrt),Hess,counter=.true.)
            ! Diagonalize and get data
            call diagonalize_full(Aux(1:Nrt,1:Nrt),Nrt,LL(1:Nrt,1:Nrt),Freq(1:Nrt),"lapack")
            ! 1b) Get T modes fromt he 3x3 block in the Eckart frame
            Nrt = 3*Nat - Nvib
            ! If there is an additional, do not sort it with rotations
            if (rm_custom_coord/="none".or.rm_gradcoord) &
              Nrt=Nrt-1
            ! D [M^-1/2 [Hx] M^1/2] D^t
            Aux(1:Nrt-3,1:Nrt-3) = matrix_basisrot(Nrt-3,3*Nat,D(1:3*Nat,4:Nrt),Hess,counter=.true.)
            ! Diagonalize and get data
            call diagonalize_full(Aux(1:Nrt-3,1:Nrt-3),Nrt-3,LL(4:Nrt,4:Nrt),Freq(4:Nrt),"lapack")
            ! If there is an additional removed coordinate, sort it out now
            if (rm_custom_coord/="none".or.rm_gradcoord) then
                Nrt=Nrt+1
                ! D [M^-1/2 [Hx] M^1/2] D^t
                Aux(1:1,1:1) = matrix_basisrot(1,3*Nat,D(1:3*Nat,Nrt:Nrt),Hess,counter=.true.)
                ! Diagonalize and get data
                call diagonalize_full(Aux(1:1,1:1),1,LL(Nrt:Nrt,Nrt:Nrt),Freq(Nrt:Nrt),"lapack")
            endif
            ! 2) Get Vib modes fromt he NvibxNvib Eckart frame
            ! D [M^-1/2 [Hx] M^1/2] D^t
            Hess(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,3*Nat,D(1:3*Nat,Nrt+1:3*Nat),Hess,counter=.true.)
            ! Diagonalize and get data
            call diagonalize_full(Hess(1:Nvib,1:Nvib),Nvib,LL(Nrt+1:3*Nat,Nrt+1:3*Nat),Freq(Nrt+1:3*Nat),"lapack")
            ! 3) Update "Nvib" to actual number of computed modes
            Nvib = 3*Nat
            Nrt  = 0
        else if (full_diagonalize) then
            call statement(6,"Getting T+R+Vib in the same frame")
            ! 1) Get T+R+Vib in the same frame
            ! D [M^-1/2 [Hx] M^1/2] D^t
            Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,D(1:3*Nat,1:3*Nat),Hess,counter=.true.)
            ! Check couplings between internal and external
            ! Translations
            Theta = 0.d0
            Theta2= 0.d0
            k=0
            do i=1,3
            do j=7,3*Nat 
                Theta =Theta+dabs(Hess(i,j))
                Theta2=max(Theta2,abs(Hess(i,j)))
                k=k+1
            enddo
            enddo
            print*, "Translation-Vibration couplings"
            print*, " Average: ", Theta/k 
            print*, " Max    : ", Theta2
            ! Rotations 
            Theta = 0.d0
            Theta2= 0.d0
            k=0
            do i=4,6
            do j=7,3*Nat 
                Theta =Theta+dabs(Hess(i,j))
                Theta2=max(Theta2,abs(Hess(i,j)))
                k=k+1
            enddo
            enddo
            print*, "Rotation-Vibration couplings"
            print*, " Average: ", Theta/k 
            print*, " Max    : ", Theta2 
            ! Check t7
            Theta = 0.d0
            Theta2= 0.d0
            k=0
            do i=7,7
            do j=8,3*Nat 
                Theta =Theta+dabs(Hess(i,j))
                Theta2=max(Theta2,abs(Hess(i,j)))
                k=k+1
            enddo
            enddo
            print*, "t7-Vibration-t7 couplings"
            print*, " Average: ", Theta/k 
            print*, " Max    : ", Theta2 
            Theta = 0.d0
            Theta2= 0.d0
            ! Diagonalize and get data
            call diagonalize_full(Hess(1:3*Nat,1:3*Nat),3*Nat,LL(1:3*Nat,1:3*Nat),Freq(1:3*Nat),"lapack")
            ! 2) Update "Nvib" to actual number of computed modes
            Nvib = 3*Nat
            Nrt  = 0
        else ! Only vibrations, in the Eckart frame
            call statement(6,"Getting only Vib in the Eckart frame")
            ! 1) Get Vib in the Eckart frame
            !! D [M^-1/2 [Hx] M^1/2] D^t
            Hess(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,3*Nat,D(1:3*Nat,Nrt+1:3*Nat),Hess,counter=.true.)
            ! Diagonalize and get data
            call diagonalize_full(Hess(1:Nvib,1:Nvib),Nvib,LL(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
        endif

        ! Rotate the modes to the Cartesian frame
        Nrt = 3*Nat - Nvib
        LL(1:3*Nat,1:Nvib) = matrix_product(3*Nat,Nvib,Nvib,D(1:3*Nat,Nrt+1:3*Nat),LL(1:Nvib,1:Nvib))

        !Check FC
        if (verbose>1) &
         call print_vector(6,Freq*1.d6,Nvib,"FORCE CONSTANTS x 10^6 (A.U.)")
        !Transform to FC to Freq
        do i=1,Nvib
              Freq(i) = sign(dsqrt(abs(Freq(i))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
                             Freq(i))
        enddo
        if (verbose>0) &
            call print_vector(6,Freq,Nvib,"Frequencies (cm-1)")

        ! If printing modes, take Lcar^-1
        if (print_modes) then
            Lmwc(1:3*Nat,1:Nvib) = LL(1:3*Nat,1:Nvib)
            do i=1,Nvib
            do j=1,3*Nat
                Lcartinv(j,i) = 0.d0
                jj = (j-1)/3+1
                Lcartinv(j,i) = Lcartinv(j,i) + LL(j,i) * dsqrt(molecule%atom(jj)%mass)
            enddo 
            enddo
        endif

        !Transform L to Cartesian (output in AU(mass) as with internal)
        call Lmwc_to_Lcart(Nat,Nvib,molecule%atom(1:Nat)%mass,LL,LL)

    endif ! Read normal modes OR do vibrationa analisis

    ! Get extra data (g_Q...)
    if (verbose>1) then
        ! Gradient
        ! gQ = L^t gx
        do i=1,Nvib
            Vec(i) = 0.d0
            do k=1,3*Nat
                Vec(i) = Vec(i) + LL(k,i) * Grad(k)
            enddo
        enddo
        call print_vector(6,Vec*1e5,Nvib,"Grad_Q (a.u.) x10^5")
    endif

    !Define the Factor to convert shift into addimensional displacements
    ! from the shift in SI units:
    Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib))*1.d2*clight*2.d0*PI/plankbar)
    ! but we need it from au not SI
    Factor(1:Nvib)=Factor(1:Nvib)*BOHRtoM*dsqrt(AUtoKG)

    ! Check if we should continue
    if (.not.animate) then
        call cpu_time(tf)
        write(6,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti
        stop
    endif

    ! Get the selection of normal modes to represent
    if (adjustl(selection) == "all") then
        Nsel=Nvib
        nm(1:Nvib) = (/(i, i=1,Nvib)/)
    else if (adjustl(selection) == "none") then
        Nsel=0
    else
        call selection2intlist(selection,nm,Nsel)
    endif

    if (Tthermo /= 0.d0) then
        ! Do thermochemical analysis
        call thermo(Nat,Nvib,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,molecule%atom(:)%Mass,Freq,Tthermo)
    endif


    !==========================================================0
    !  Normal mode displacements
    !==========================================================0
    call set_geom_units(molecule,"Angs")
    ! Initialization
    X0(1:Nat) = molecule%atom(1:Nat)%x
    Y0(1:Nat) = molecule%atom(1:Nat)%y
    Z0(1:Nat) = molecule%atom(1:Nat)%z
    Nsteps = 101
    if ( mod(Nsteps,2) == 0 ) Nsteps = Nsteps + 1 ! ensure odd number of steps (so we have same left and right)
    ! Qstep is dimless
    Qstep = Amplitude/float(Nsteps-1)*2.d0  ! Do the range (-A ... +A)
    molecule%atom(1:molecule%natoms)%resname = "RES" ! For printing
    ! Run over all selected modes/internals
    do jj=1,Nsel 
        k=0 ! equilibrium corresponds to k=0
        j = nm(jj)
        if (verbose>0.or.print_modes) &
         print'(X,A,I0,A)', "Generating Mode ", j, "..."
        if (print_modes) then
            print*, "   Lcart^-1              Lmwc^t"
            do i=1,3*Nat
            print'(2ES22.12)', Lcartinv(i,j), Lmwc(i,j)
            enddo 
        endif

        ! Set initial values for the scanned coordinate
        qcoord = 0.d0 

        ! Prepare and open files
        call prepare_files(j,grofile,g09file,g96file,numfile,numfwfile,numbwfile,qfile,arrowfile,title,full_diagonalize)
        open(O_GRO,file=grofile)
        open(O_G09,file=g09file)
        open(O_G96,file=g96file)
        open(O_Q  ,file=qfile)
        open(O_NUM,file=numfile)
        open(O_ARR,file=arrowfile)
        
        ! Print arrowfile
        do k=1,Nat
            write(O_ARR,*) k, LL(3*k-2:3*k,j)
        enddo
        close(O_ARR)

        !===========================
        !Start from equilibrium. 
        !===========================
        molecule%atom(1:Nat)%x = X0(1:Nat)
        molecule%atom(1:Nat)%y = Y0(1:Nat)
        molecule%atom(1:Nat)%z = Z0(1:Nat)
        if (verbose>1) &
         print'(/,A,I0)', "STEP:", k
        ! Update
        write(molecule%title,'(A,I0,A,2(X,F12.6))') &
         trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
        ! Print initial structure
        call write_gro(O_GRO,molecule)
        !===========================
        !Half Forward oscillation: from step "Eq + dQ" to "Eq + (N-1)/2 dQ"
        !===========================
        !Initialize distacen criterion for rmsd_fit_frame_brute SR
        dist=0.d0
        do istep = 1,(nsteps-1)/2
            k=k+1
            if (verbose>1) &
             print'(/,A,I0)', "STEP:", k
            ! Update values
            ! qcoord has AU
            qcoord = qcoord + Qstep/Factor(j)
            write(molecule%title,'(A,I0,A,2(X,F12.6))') &
             trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
            ! Displace in AA (displace always from reference to avoid error propagation)
            i=istep
            molecule%atom(1:Nat)%x = X0(1:Nat)
            molecule%atom(1:Nat)%y = Y0(1:Nat)
            molecule%atom(1:Nat)%z = Z0(1:Nat)
            call displace_Xcoord(LL(1:3*Nat,j),molecule%natoms,Qstep/Factor(j)*BOHRtoANGS*i,&
                                 molecule%atom(:)%x,molecule%atom(:)%y,molecule%atom(:)%z)
            ! PRINT
            ! Write GRO from the beginign and G96/G09 only when reach max amplitude
            call write_gro(O_GRO,molecule)
            if (k==(nsteps-1)/2) then
                call write_gcom(O_G09,molecule,g09file,calc,method,basis,molecule%title)
                call write_g96(O_G96,molecule)
                write(O_Q,*) qcoord, qcoord*Factor(j)
            endif
        enddo
        !=======================================
        ! Reached amplitude. Back oscillation: from step "MaxAmp + dQ" to "MaxAmp + (N-2) dQ" == -MaxAmp
        !=======================================
        ! This is the part reported in G09/G96 files
        do istep = 1,nsteps-1
            k=k+1
            if (verbose>1) &
             print'(/,A,I0)', "STEP:", k
            ! Update values
            qcoord = qcoord - Qstep/Factor(j)
            write(molecule%title,'(A,I0,A,2(X,F12.6))') &
             trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
            ! Displace in AA (displace always from reference to avoid error propagation)
            i=(nsteps-1)/2-istep
            molecule%atom(1:Nat)%x = X0(1:Nat)
            molecule%atom(1:Nat)%y = Y0(1:Nat)
            molecule%atom(1:Nat)%z = Z0(1:Nat)
            call displace_Xcoord(LL(1:3*Nat,j),molecule%natoms,Qstep/Factor(j)*BOHRtoANGS*i,&
                                 molecule%atom(:)%x,molecule%atom(:)%y,molecule%atom(:)%z)
            ! PRINT
            ! Write G96/GRO every step and G09 scan every 10 steps
            ! except the 5 poinst around minimum, which are all printed
            call write_gro(O_GRO,molecule)
            call write_g96(O_G96,molecule)
            if (mod(k,10) == 0) then
                call write_gcom(O_G09,molecule,g09file,calc,method,basis,molecule%title)
                write(O_Q,*) qcoord, qcoord*Factor(j)
            endif
            ! Write 5 poinst around minimum for numerical dierivatives
            if (k>=nsteps-3.and.k<=nsteps+1) then
                call write_gcom(O_NUM,molecule,numfile,calc,method,basis,molecule%title)
            endif
            if (k==nsteps-2) then
                open(O_NUMD,file=numfwfile)
                method_=adjustl(method)//" TD NoSym"
                call write_gcom(O_NUMD,molecule,numfwfile,calc,method_,basis,molecule%title)
                close(O_NUMD)
            endif
            if (k==nsteps) then
                open(O_NUMD,file=numbwfile)
                method_=adjustl(method)//" TD NoSym"
                call write_gcom(O_NUMD,molecule,numbwfile,calc,method_,basis,molecule%title)
                close(O_NUMD)
            endif
        enddo
        !=======================================
        ! Reached amplitude. Half Forward oscillation (till one step before equilibrium, so that we concatenate well)
        !=======================================
        do istep = 1,(nsteps-1)/2-1
            k=k+1
            if (verbose>1) &
             print'(/,A,I0)', "STEP:", k
            ! Update values
            qcoord = qcoord + Qstep/Factor(j)
            write(molecule%title,'(A,I0,A,2(X,F12.6))') &
             trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
            ! Displace in AA (displace always from reference to avoid error propagation)
            i=-(nsteps-1)/2+istep
            molecule%atom(1:Nat)%x = X0(1:Nat)
            molecule%atom(1:Nat)%y = Y0(1:Nat)
            molecule%atom(1:Nat)%z = Z0(1:Nat)
            call displace_Xcoord(LL(1:3*Nat,j),molecule%natoms,Qstep/Factor(j)*BOHRtoANGS*i,&
                                 molecule%atom(:)%x,molecule%atom(:)%y,molecule%atom(:)%z)
            ! PRINT
            ! Write only GRO 
            call write_gro(O_GRO,molecule)
        enddo
        open(O_GRO)
        open(O_G09)
        open(O_G96)
        open(O_Q  )
        open(O_NUM)
    enddo

    call summary_alerts

    call cpu_time(tf)
    write(6,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    ! CALL EXTERNAL PROGRAM TO RUN ANIMATIONS

    if (call_vmd) then
        open(S_VMD,file="vmd_conf.dat",status="replace")
        !Set general display settings (mimic gv)
        write(S_VMD,*) "color Display Background iceblue"
        write(S_VMD,*) "color Name {C} silver"
        write(S_VMD,*) "axes location off"
        !Set molecule representation
        do i=0,Nsel-1
            j = nm(i+1)
            write(S_VMD,*) "mol representation CPK"
!            write(S_VMD,*) "mol addrep 0"
            if (i==0) then
                write(S_VMD,*) "molinfo ", i, " set drawn 1"
            else
                write(S_VMD,*) "molinfo ", i, " set drawn 0"
            endif
            write(S_VMD,*) "mol addrep ", i
            write(dummy_char,'(A,I4,X,F8.2,A)') "{Mode",j, Freq(j),"cm-1}"
            dummy_char=trim(adjustl(dummy_char))
            write(S_VMD,*) "mol rename ", i, trim(dummy_char)
        enddo
        write(S_VMD,*) "display projection Orthographic"
        close(S_VMD)
        !Call vmd
        vmdcall = 'vmd -m '
        do i=1,Nsel
            j = nm(i)
            ! Get filenames (we want grofile name)
            call prepare_files(j,grofile,g09file,g96file,numfile,numfwfile,numbwfile,qfile,arrowfile,title,full_diagonalize)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile))
        enddo
        vmdcall = trim(adjustl(vmdcall))//" -e vmd_conf.dat"
        open(O_MOV,file="vmd_call.cmd")
        write(O_MOV,*) trim(adjustl(vmdcall))
        close(O_MOV)
        call system(vmdcall)
    endif

    if (movie_cycles > 0) then
        open(S_VMD,file="vmd_movie.dat",status="replace")
        !Set general display settings (mimic gv)
        write(S_VMD,*) "color Display Background white"
        write(S_VMD,*) "color Name {C} silver"
        write(S_VMD,*) "axes location off"
        write(S_VMD,*) "display projection Orthographic"
        !Set molecule representation
        do i=0,Nsel-1
            j = nm(i+1)
            ! Get filenames (we want grofile name)
            call prepare_files(j,grofile,g09file,g96file,numfile,numfwfile,numbwfile,qfile,arrowfile,title,full_diagonalize)
            write(S_VMD,*) "mol representation CPK"
            write(S_VMD,*) "molinfo ", i, " set drawn 0"
            write(S_VMD,*) "mol addrep ", i
            write(dummy_char,'(A,I4,X,F8.2,A)') "{Mode",j, Freq(j),"cm-1}"
            dummy_char=trim(adjustl(dummy_char))
            write(S_VMD,*) "mol rename ", i, trim(dummy_char)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile))
        enddo
        write(S_VMD,'(A)') "#====================="
        write(S_VMD,'(A)') "# Start movies"
        write(S_VMD,'(A)') "#====================="
        !Set length of the movie
        movie_steps = movie_cycles*20
        do i=0,Nsel-1
            j = nm(i+1)
            write(S_VMD,'(A,I4)') "# Mode", j
            write(tmpfile,*) j
            tmpfile="Mode"//trim(adjustl(tmpfile))
            write(S_VMD,*) "molinfo ", i, " set drawn 1"
            write(S_VMD,*) "set figfile "//trim(adjustl(tmpfile))
            write(S_VMD,'(A,I3,A)') "for {set xx 0} {$xx <=", movie_steps,&
                                    "} {incr xx} {"
            write(S_VMD,*) "set x [expr {($xx-($xx/20)*20)*10}]"
            write(S_VMD,*) 'echo "step $x"'
            write(S_VMD,*) "animate goto $x"
            write(S_VMD,*) "render Tachyon $figfile-$xx.dat"
            write(S_VMD,'(A)') '"/usr/local/lib/vmd/tachyon_LINUX" -aasamples 12 '//& 
                           '$figfile-$xx.dat -format TARGA -o $figfile-$xx.tga'
            write(S_VMD,'(A)') 'convert -font URW-Palladio-Roman -pointsize 30 -draw '//&
                           '"text 30,70 '//"'"//trim(adjustl(tmpfile))//"'"//&
                           '" $figfile-$xx.tga $figfile-$xx.jpg'
            write(S_VMD,'(A)') "}"
            !Updated ffmpeg call. The output is now loadable from ipynb
            write(S_VMD,'(A)') 'ffmpeg -i $figfile-%d.jpg -vcodec libx264 -s 640x360 $figfile.mp4'
            write(S_VMD,*) "molinfo ", i, " set drawn 0"
        enddo
        write(S_VMD,*) "exit"
        close(S_VMD)
        !Call vmd
        vmdcall = 'vmd -m '
        do i=1,Nsel
        vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile))
        enddo
        vmdcall = trim(adjustl(vmdcall))//" -e vmd_movie.dat -size 500 500"
        open(O_MOV,file="movie.cmd",status="replace")
        write(O_MOV,'(A)') trim(adjustl(vmdcall))
        write(O_MOV,'(A)') "rm Mode*jpg Mode*dat Mode*tga"
        close(O_MOV)
        print*, ""
        print*, "============================================================"
        print*, "TO GENERATE THE MOVIES (AVI) EXECUTE COMMANDS IN 'movie.cmd'"
        print*, "(you may want to edit 'vmd_movie.dat'  first)"
        print*, "============================================================"
        print*, ""
    endif

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(&
                           ! input data
                           inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn,mass_file,&
                           ! Options (general)
                           Amplitude,call_vmd,include_hbonds,selection,vertical, &
                           outfchkfile,                                          &
                           ! Options (Cartesian)
                           full_diagonalize,                                     &
                           ! Remove coordinate along gradient
                           rm_gradcoord,rm_custom_coord,rm_custom_mode,          &
                           ! Movie
                           animate,movie_vmd, movie_cycles,print_modes,          &
                           ! Options (internal)
                           use_symmetry,def_internal,intfile,rmzfile,            &
                           ! Additional vib options
                           Eckart_frame,orthogonalize,modes_as_internals,        &
                           original_internal, apply_projection_matrix,           &
                           ! connectivity file
                           cnx_file,                                             &
                           ! thermochemical analysis
                           Tthermo,                                              &
                           ! (hidden)
                           analytic_Bder)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn,selection, &
                                          !Internal
                                          def_internal,intfile,rmzfile,cnx_file, &
                                          rm_custom_coord,rm_custom_mode, mass_file, outfchkfile
        real(8),intent(inout)          :: Amplitude, Tthermo
        logical,intent(inout)          :: call_vmd, include_hbonds,vertical,movie_vmd,full_diagonalize,animate,&
                                          rm_gradcoord, &
                                          ! Internal
                                          use_symmetry,analytic_Bder, &
                                          ! Other
                                          Eckart_frame, orthogonalize, modes_as_internals, original_internal, &
                                          apply_projection_matrix, print_modes
        integer,intent(inout)          :: movie_cycles

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        character(len=500) :: input_command
        ! iargc type must be specified with implicit none (strict compilation)
        integer :: iargc
        
        !Initialize input_command
        call getarg(0,input_command)
        !Get input flags
        do i=1,iargc()
            call getarg(i,arg)
            input_command = trim(adjustl(input_command))//" "//trim(adjustl(arg))
        enddo

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-f") 
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.
                case ("-ft") 
                    call getarg(i+1, ft)
                    argument_retrieved=.true.

                case ("-fhess") 
                    call getarg(i+1, hessfile)
                    argument_retrieved=.true.
                case ("-fth") 
                    call getarg(i+1, fth)
                    argument_retrieved=.true.

                case ("-fgrad") 
                    call getarg(i+1, gradfile)
                    argument_retrieved=.true.
                case ("-ftg") 
                    call getarg(i+1, ftg)
                    argument_retrieved=.true.

                case ("-fnm") 
                    call getarg(i+1, nmfile)
                    argument_retrieved=.true.
                case ("-ftn") 
                    call getarg(i+1, ftn)
                    argument_retrieved=.true.

                case ("-fmass") 
                    call getarg(i+1, mass_file)
                    argument_retrieved=.true.

                case ("-rmgrad")
                    rm_gradcoord=.true.
                case ("-normgrad")
                    rm_gradcoord=.false.

                case ("-rmcoord") 
                    call getarg(i+1, rm_custom_coord)
                    argument_retrieved=.true.

                case ("-rmmode") 
                    call getarg(i+1, rm_custom_mode)
                    argument_retrieved=.true.

                case ("-Eckart")
                    Eckart_frame=.true.
                case ("-noEckart")
                    Eckart_frame=.false.

                case ("-origint")
                    original_internal=.true.
                case ("-noorigint")
                    original_internal=.false.

                case ("-fulldiag")
                    full_diagonalize=.true.
                case ("-nofulldiag")
                    full_diagonalize=.false.

                case ("-orth")
                    orthogonalize=.true.
                case ("-noorth")
                    orthogonalize=.false.

                case ("-modes2int")
                    modes_as_internals=.true.
                case ("-nomodes2int")
                    modes_as_internals=.false.

                case ("-nm") 
                    call getarg(i+1, selection)
                    argument_retrieved=.true.

                case ("-disp") 
                    call getarg(i+1, arg)
                    read(arg,*) Amplitude
                    argument_retrieved=.true.

                case ("-vmd")
                    call_vmd=.true.
                case ("-novmd")
                    call_vmd=.false.

                case ("-animate")
                    animate=.true.
                case ("-noanimate")
                    animate=.false.

                case ("-movie")
                    call getarg(i+1, arg)
                    read(arg,*) movie_cycles
                    movie_vmd=.true.
                    argument_retrieved=.true.

                case ("-include_hb")
                    include_hbonds=.true.

                ! Internal (for correction) options
                case ("-vert")
                    vertical=.true.
                case ("-novert")
                    vertical=.false.
                    
                case ("-writefchk") 
                    call getarg(i+1, outfchkfile)
                    argument_retrieved=.true.

                case ("-sym")
                    use_symmetry=.true.
                case ("-nosym")
                    use_symmetry=.false.

                case ("-cnx") 
                    call getarg(i+1, cnx_file)
                    argument_retrieved=.true.

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-rmzfile") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.

                case ("-intmode")
                    call getarg(i+1, def_internal)
                    argument_retrieved=.true.

                case ("-thermo")
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    read(arg,*) Tthermo

                case ("-prj-tr")
                    apply_projection_matrix=.true.
                case ("-noprj-tr")
                    apply_projection_matrix=.false.

                case ("-print")
                    print_modes=.true.
                case ("-noprint")
                    print_modes=.false.


                ! (HIDDEN FLAG)
                case ("-anaBder")
                    analytic_Bder=.true.
                case ("-noanaBder")
                    analytic_Bder=.false.

        
                case ("-h")
                    need_help=.true.

                ! Control verbosity
                case ("-quiet")
                    verbose=0
                    silent_notes = .true.
                case ("-concise")
                    verbose=1
                case ("-v")
                    verbose=2
                case ("-vv")
                    verbose=3
                    silent_notes=.false.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       ! Manage defaults
       ! If not declared, hessfile and gradfile are the same as inpfile
       ! unless we are using nm file
       if (adjustl(nmfile) == "none") then
           if (adjustl(hessfile) == "same") then
               hessfile=inpfile
               if (adjustl(fth) == "guess")  fth=ft
           endif
           if (adjustl(gradfile) == "same") then
               gradfile=inpfile
               if (adjustl(ftg) == "guess")  ftg=ft
           endif
           ftn="-"
       else
           if (adjustl(hessfile) /= "same") &
            call alert_msg("note","Using nm file, disabling Hessian file")
           hessfile="none"
           fth="-"
           if (adjustl(gradfile) /= "same") &
            call alert_msg("note","Using nm file, disabling gradient file")
           gradfile="none"
           ftg="-"
       endif


       !Print options (to stderr)
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '             N M    C A R T E S I A N '    
        write(6,'(/,A)') '      Perform vibrational analysis based on  '
        write(6,'(A)')   '             Cartesian coordinates '        
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag           Description                     Value'
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,*)       '-f             Input file (structure&default)  ', trim(adjustl(inpfile))
        write(6,*)       '-ft            \_ FileType                     ', trim(adjustl(ft))
        write(6,*)       '-fhess         Hessian file                    ', trim(adjustl(hessfile))
        write(6,*)       '-fth           \_ FileType                     ', trim(adjustl(fth))
        write(6,*)       '-fgrad         Hessian file                    ', trim(adjustl(gradfile))
        write(6,*)       '-ftg           \_ FileType                     ', trim(adjustl(ftg))
        write(6,*)       '-fnm           Gradient file                   ', trim(adjustl(nmfile))
        write(6,*)       '-ftn           \_ FileType                     ', trim(adjustl(ftn))
        write(6,*)       '-fmass         Mass file (optional)            ', trim(adjustl(mass_file))       
        write(6,*)       ''
        write(6,*)       ' ** Options for vibration analysis **'
        write(6,*)       '-[no]prj-tr    Project out tras+rot           ', apply_projection_matrix
        write(6,*)       '-[no]rmgrad    Remove coordinate along the    ', rm_gradcoord
        write(6,*)       '               grandient                      '
        write(6,*)       '-rmcoord       Remove custom coordinate       ', rm_custom_coord
        write(6,*)       '               (in this file)                      '
        write(6,*)       '-rmmode        Remove custom mode             ', rm_custom_mode
        write(6,*)       '               (in this file)                      '
        write(6,*)       '-[no]fulldiag  Diagonalize the 3Nx3N matrix   ',  full_diagonalize
        write(6,*)       '-[no]Eckart    Include translation & rotation ',  Eckart_frame
        write(6,*)       '               in the Eckart frame            '
        write(6,*)       '-[no]orth      Use orthogonalized internals   ',  orthogonalize
        write(6,*)       '-[no]modes2int Use uncorrected modes as the   ',  modes_as_internals
        write(6,*)       '               definiton of the internals     '
        write(6,*)       '-[no]origint   Use originally defined inter-  ',  original_internal
        write(6,*)       '               nal w\o linear combinations    '
        write(6,*)       '               (valid for non-redundant sets) '
        write(6,*)       ''
        write(6,*)       ' ** Options for themochemistry **'
        write(6,*)       '-thermo        Temp (K) for thermochemistry   ', Tthermo
        write(6,*)       '               (0.0 means no analysis)'
        write(6,*)       ''
        write(6,*)       ' ** Options for animation **'
        write(6,*)       '-[no]animate   Generate animation files       ',  animate
        write(6,*)       '-nm            Selection of normal modes to    ', trim(adjustl(selection))
        write(6,*)       '               generate animations            '
        write(6,'(X,A,F5.2)') &
                         '-disp          Mode displacements for animate ',  Amplitude
        write(6,*)       '               (dimensionless displacements)'
        write(6,*)       '-[no]vmd       Launch VMD after computing the ',  call_vmd
        write(6,*)       '               modes (needs VMD installed)'
        write(6,'(X,A,I0)') &
                         '-movie         Number of cycles to record on   ',  movie_cycles
        write(6,*)       '               a movie with the animation'
        write(6,*)       '-[no]print     Print modes in MWC              ',  print_modes
        write(6,*)       ''
        write(6,*)       ' ** Options for correction non-stationary points **'
        write(6,*)       '-[no]vert      Correct with B derivatives for ',  vertical
        write(6,*)       '               non-stationary points'
        write(6,*)       '-writefchk     Write "corrected" fchk file     ', trim(adjustl(outfchkfile))
        write(6,*)       '-cnx           Connectivity [filename|guess]   ', trim(adjustl(cnx_file))
        write(6,*)       '-intmode       Internal set:[zmat|sel|all]     ', trim(adjustl(def_internal))
        write(6,*)       '-intfile       File with ICs (for "sel")       ', trim(adjustl(intfile))
        write(6,*)       '-[no]sym       Use symmetry to form Zmat      ',  use_symmetry
!         write(6,*)       '-rmzfile       File deleting ICs from Zmat     ', trim(adjustl(rmzfile))
        write(6,*)       ''
        write(6,*)       '-h             Display this help              ',  need_help
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(A)') 'Input command:'
        write(6,'(A)') trim(adjustl(input_command))   
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                       'Verbose level:  ', verbose        
        write(6,'(A,/)') '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input

    subroutine prepare_files(icoord,grofile,g09file,g96file,numfile,numfwfile,numbwfile,qfile,arrowfile,title,fulldiag)

        integer,intent(in) :: icoord
        character(len=*),intent(out) :: grofile,g09file,g96file,numfile,numfwfile,numbwfile,qfile,arrowfile,title
        logical,intent(in) :: fulldiag

        !Local
        character(len=150) :: dummy_char
        character(len=5)   :: full_label

        if (full_diagonalize) then
            full_label="full-"
        else
            full_label=""
        endif

        write(dummy_char,"(I0,X,A)") icoord
        title   = "Animation of normal mode "//trim(adjustl(dummy_char))
        g09file = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart.com"
        g96file = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart.g96"
        qfile   = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart_steps.dat"
        grofile = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart.gro"
        numfile = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart_num.com"
        numfwfile= "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart_numder_fw.com"
        numbwfile= "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart_numder_bw.com"
        arrowfile= "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart_arrows.dat"

        return

    end subroutine prepare_files


    subroutine displace_Xcoord(Lc,Nat,Qstep,X,Y,Z)

        real(8),dimension(:),intent(in)   :: Lc
        real(8),intent(in)                :: Qstep 
        integer,intent(in)                :: Nat
        real(8),dimension(:),intent(inout):: X,Y,Z

        !Local
        integer :: k, kk
        real(8) :: magic_number

!         magic_number=1.8895d0

        do k=1,Nat
            kk = 3*(k-1)
            X(k) = X(k) + Lc(kk+1) * Qstep
            Y(k) = Y(k) + Lc(kk+2) * Qstep
            Z(k) = Z(k) + Lc(kk+3) * Qstep
        enddo

        return
       
    end subroutine displace_Xcoord

end program normal_modes_cartesian

