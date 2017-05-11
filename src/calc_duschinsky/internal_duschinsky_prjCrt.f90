program internal_duschinski


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Version: internal_duschinski_v13
    !
    ! Description:
    ! -----------
    ! Program to analyse vibrations in term of internal coordinates.
    !
    ! NOTES
    ! For vertical, the analysis of the nm displacement can be done:
    !  * in the S-space (internal coordinates). Activated with -vert
    !  * in the Q-space (final state nm). Activated with -vert2
    ! The same convention is used to compute Er
    ! For AH, we always work in the S-space
    !==============================================================


    !*****************
    !   MODULE LOAD
    !*****************
    !============================================
    !   Generic
    !============================================
    use io
    use alerts
    use line_preprocess
    use constants 
    use verbosity
    use matrix
    use matrix_print
    !============================================
    !   Structure types module
    !============================================
    use structure_types
    !============================================
    !   File readers
    !============================================
    use generic_io
    use generic_io_molec
    !============================================
    !  Structure-related modules
    !============================================
    use molecular_structure
    use atomic_geom
    use symmetry
    !============================================
    !  Internal thingies
    !============================================
    use internal_module
    use zmat_manage 
    use vibrational_analysis
    use vertical_model

    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: use_symmetry=.false.   ,&
               modred=.false.         ,&
               tswitch=.false.        ,&
               symaddapt=.false.      ,&
               vertical=.false.       ,&
               verticalQspace2=.false.,&
               verticalQspace1=.false.,&
               gradcorrectS1=.false.  ,&  
               gradcorrectS2=.false.  ,&
               same_red2nonred_rotation=.true., &
               analytic_Bder=.true., &
               check_symmetry=.true., &
               orthogonalize=.false., &
               original_internal=.false., &
               force_real=.false., &
               apply_projection_matrix=.false., &
               complementay_projection=.false., &
               move_to_min=.false.
    character(len=4) :: def_internal='all'
    character(len=4) :: def_internal0='defa' ! defa(ult) is "the same as working set"
    character(len=1) :: reference_frame='F'
    character(len=10):: model="adia"
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1, state2
    type(str_bonded) :: geomS, geom0
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns, NNvib, Nvib0, Ns0, NsS, Nz, Nf, NvibP, NvibP2, Nvib1, Nvib2
    character(len=5) :: PG
    !Bonded info
    integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !B and G matrices
    real(8),dimension(NDIM,NDIM) :: B1, B2, Bprj, Bprj2
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X, X1inv,X2inv, L1,L2,L1inv, &
                                        Asel1,Asel2,Asel1inv,Asel2inv, gBder,Asel, &
                                        G0, B0, P, Fltr
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    !Duschisky
    real(8),dimension(NDIM,NDIM) :: G1, G2, Jdus
    !T0 - switching effects
    real(8),dimension(3,3) :: T
    !AUXILIAR MATRICES
    real(8),dimension(NDIM,NDIM) :: Aux, Aux2
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: Freq1, Freq2, S1, S2, Vec1, Vec2, Q0, FC
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !Shifts
    real(8),dimension(NDIM) :: Delta
    real(8) :: Delta_p, Er
    ! Z-mat and redundant geoms
    integer,dimension(NDIM) :: Zmap
    type(str_bonded) :: zmatgeom, allgeom, inputgeom
    !====================== 

    !=========================
    ! Distance calculation stuff
    logical :: get_distances=.true.
    real(8) :: ff, f1, f0, time, dt, dist, area
    !====================== 

    !====================== 
    !
    real(8),dimension(:),allocatable :: A
    integer :: error
    character(len=200) :: msg
    character(len=5)   :: current_symm
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=16) :: dummy_char
    real(8) :: Theta, Theta2, Theta3
    ! RMZ things
    character :: rm_type
    integer :: rm_zline, Nrm, nbonds_rm, nangles_rm, ndiheds_rm
    integer,dimension(100) :: bond_rm, angle_rm, dihed_rm
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, k90,k95,k99, nn, imin, imax,&
               i1,i2,i3,i4, iop
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_ZMAT=11, &
               I_SYM=12,  &
               I_RED=13,  &
               I_ADD=14,  &
               I_AD2=15,  &
               I_RMF=16,  &
               I_CNX=17,  &
               I_MAS=18,  &
               O_DUS=20,  &
               O_DIS=21,  &
               O_DMAT=22, &
               O_DUS2=23, &
               O_DIS2=24, &
               O_STAT=25
    !files
    character(len=10) :: ft ="guess",  ftg="guess",  fth="guess", &
                         ft2 ="guess", ftg2="guess", fth2="guess"
    character(len=200):: inpfile  ="state1.fchk", &
                         gradfile ="same", &
                         hessfile ="same", &
                         inpfile2 ="none", &
                         gradfile2="same", &
                         hessfile2="same", &
                         intfile  ="none", &
                         intfile0 ="default", & ! default is "the same as working set"
                         rmzfile  ="none", &
                         symm_file="none", &
                         cnx_file="guess", &
                         mass_file="none"
    !status
    integer :: IOstatus
    !===================

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

! (End of variables declaration) 
!==================================================================================

    call cpu_time(ti)

    !--------------------------
    ! Tune io
    !--------------------------
    ! Set unit for alert messages
    alert_unt=6
    !--------------------------

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(state1)
    call allocate_atoms(state2)
    !===========================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,&
                     inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                     intfile,intfile0,rmzfile,def_internal,def_internal0,&
                     use_symmetry,cnx_file,mass_file,&
!                    tswitch,
                     symaddapt,same_red2nonred_rotation,analytic_Bder,&
                     model,vertical,verticalQspace2,verticalQspace1,&
                     gradcorrectS1,gradcorrectS2,&
                     orthogonalize,original_internal,force_real,reference_frame,&
                     apply_projection_matrix,complementay_projection,move_to_min)
    call set_word_upper_case(def_internal)
    call set_word_upper_case(reference_frame)
    call set_word_upper_case(model)

    ! 1. INTERNAL VIBRATIONAL ANALYSIS ON STATE1 AND STATE2

    !===========
    !State 1
    !===========
    if (verbose>0) then
        print*, ""
        print*, "=========="
        print*, " STATE 1"
        print*, "=========="
    endif
 
    ! READ DATA (each element from a different file is possible)
    ! ---------------------------------
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)
    if (fth == "guess") &
    call split_line_back(hessfile,".",null,fth)
    if (ftg == "guess") &
    call split_line_back(gradfile,".",null,ftg)

    ! STRUCTURE FILE
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,state1,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (State1)")
    close(I_INP)
    ! Shortcuts
    Nat = state1%natoms

    ! Read mass from file if given
    if (adjustl(mass_file) /= "none") then
        print'(/,X,A)', "Reading atomic masses from: "//trim(adjustl(mass_file))
        open(I_MAS,file=mass_file,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(mass_file)) )
        do i=1,Nat
            read(I_MAS,*,iostat=IOstatus) state1%atom(i)%mass 
            if (IOstatus /= 0) call alert_msg( "fatal","While reading "//trim(adjustl(mass_file)) )
        enddo
        close(I_MAS)
    endif

    ! HESSIAN FILE
    open(I_INP,file=hessfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile)) )
    allocate(A(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fth,Nat,A,error) 
    if (error /= 0) call alert_msg("fatal","Error reading Hessian (State1)")
    close(I_INP)
    ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
    call vibrations_Cart(Nat,state1%atom(:)%X,state1%atom(:)%Y,state1%atom(:)%Z,state1%atom(:)%Mass,A,&
                         Nvib,L1,Freq1,error_flag=error)
    ! And also initialize Nvib0
    Nvib0 = Nvib
    k=0
    do i=1,3*Nat
    do j=1,i
        k=k+1
        Hess(i,j) = A(k)
        Hess(j,i) = A(k)
    enddo 
    enddo
    deallocate(A)

    if (gradcorrectS1) then
        ! GRADIENT FILE
        open(I_INP,file=gradfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
        call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
        close(I_INP)
    endif

    !----------------------------------
    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! General Actions:
    !******************
    ! Get connectivity 
    if (cnx_file == "guess") then
        call guess_connect(state1)
    else
        print'(/,A,/)', "Reading connectivity from file: "//trim(adjustl(cnx_file))
        open(I_CNX,file=cnx_file,status='old')
        call read_connect(I_CNX,state1)
        close(I_CNX)
    endif

    ! Manage symmetry
    if (.not.use_symmetry) then
        state1%PG="C1"
    else if (trim(adjustl(symm_file)) /= "none") then
        call alert_msg("note","Using custom symmetry file: "//trim(adjustl(symm_file)) )
        open(I_SYM,file=symm_file)
        do i=1,state1%natoms
            read(I_SYM,*) j, isym(j)
        enddo
        close(I_SYM)
        !Set PG to CUStom
        state1%PG="CUS"
    else
        state1%PG="XX"
        call symm_atoms(state1,isym)
    endif

    call set_geom_units(state1,"BOHR")

    !--------------------------------
    ! Compute projection matrix
    !--------------------------------
    call heading(6,"COMPUTING PROJECTION MATRIX (CARTESIAN)")

    !---------------------------------------
    call subheading(6,"Internal Sets")

    call statement(6,"Set with ALL internal coordinates from connectivity",keep_case=.true.)
    call gen_bonded(state1)
    call define_internal_set(state1,"ALL",intfile,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Fltr)
    allgeom = state1%geom

    call statement(6,"Set indicated on INPUT: "//trim(def_internal),keep_case=.true.)
    call gen_bonded(state1)
    call define_internal_set(state1,def_internal,intfile,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Fltr)
    ! If not using combinations, the Filter need to be reset to the identity matrix
    if (Nf==0) then
        Nf=Ns
        Fltr(1:Nf,1:Ns) = identity_matrix(Nf)
    endif
    inputgeom = state1%geom
    !---------------------------------------


    ! Compute B matrix in the inputgeom for projection
    call subheading(6,"Compute B with INPUT set")

    ! Compute B and apply Filter
    call internal_Wilson(state1,Ns,S1,Bprj)
    Bprj(1:Nf,1:3*Nat) = matrix_product(Nf,3*Nat,Ns,Fltr,Bprj)
    Ns=Nf

    NvibP = Nvib
    !*********************************************
    ! Really need to do redundant2nonredundant?
    ! Compute G
    call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,Bprj,G1)
    call redundant2nonredundant(Ns,NvibP,G1,Asel1)
    Bprj(1:NvibP,1:3*Nat) = matrix_product(NvibP,3*Nat,Ns,Asel1,Bprj,tA=.true.)
    !*********************************************

    ! Compute projecton matrix. Method3 provide the matrix to be applied to
    ! Cartesian Hessian/Grad, but defined on MWC
    call statement(6,"Computing projection matrix")
    P(1:3*Nat,1:3*Nat) = projection_matrix3(Nat,NvibP,Bprj,state1%atom(:)%Mass)
    if (complementay_projection) then
        Aux(1:3*Nat,1:3*Nat) = identity_matrix(3*Nat)
        P(1:3*Nat,1:3*Nat) =  Aux(1:3*Nat,1:3*Nat)-P(1:3*Nat,1:3*Nat)
    endif


    ! Get back allgeom
    call heading(6,"Computing modes with ALL set")
    Ns = allgeom%nbonds+allgeom%nangles+allgeom%ndihed
    state1%geom = allgeom


    ! Rotate gradient and Hessian
    call subheading(6,"Modes projecting the Hess/Grad")
    call statement(6,"Rotating Gradient and Hessian")
!     ! Get full Cartesian Hessian from Hlt
!     Hess(1:3*Nat,1:3*Nat) = Hlt_to_Hess(3*Nat,Hlt)
!     Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,P,Hess,counter=.true.)
    ! Get Cartesian gradient
!     Grad(1:3*Nat) = grdx
!     Grad(1:3*Nat) = matrix_vector_product(3*Nat,3*Nat,P,Grad)

    ! Compute B, G and, if needed, Bder
    call statement(6,"Computing B, G (with ALL set)",keep_case=.true.)
    call internal_Wilson(state1,Ns,S1,B1,ModeDef)
    call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,B1,G1)
    if (gradcorrectS1) then
        call statement(6,"...and Bder for non-stationary point calculation")
        call calc_Bder(state1,Ns,Bder,analytic_Bder)
    endif

    ! Get non-redundant space
    call subsubheading(6,"Getting the actual vibrational space dimension")
    call redundant2nonredundant(Ns,Nvib,G1,Asel1)
    ! Get inverse rotation
    do i=1,Nvib
    do j=1,Ns
        Asel1inv(i,j) = Asel1(j,i)
    enddo
    enddo
    call statement(6,"Rotate B, G to non-redundant space", keep_case=.true.)
    ! Rotate Bmatrix
    B1(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel1,B1,tA=.true.)
    ! Rotate Gmatrix
    G1(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel1(1:Ns,1:Nvib),G1,counter=.true.)
    if (gradcorrectS1) then
        call statement(6,"Also rotate Bder to non-redundant space", keep_case=.true.)
        ! Rotate Bders
        do j=1,3*Nat
            Bder(1:Nvib,j,1:3*Nat) =  matrix_product(Nvib,3*Nat,Ns,Asel1,Bder(1:Ns,j,1:3*Nat),tA=.true.)
        enddo

        ! Get the Correction now
        call statement(6," Getting gs^t\beta term")
        ! The correction is applied with the Nvib0 SET
        ! Correct Hessian as
        ! Hx' = Hx - gs^t\beta
        ! 1. Get gs from gx
        call Gradcart2int(Nat,Nvib,Grad,state1%atom(:)%mass,B1,G1)
        ! 2. Multiply gs^t\beta and
        ! 3. Apply the correction
        ! Bder(i,j,K)^t * gq(K)
        do i=1,3*Nat
        do j=1,3*Nat
            gBder(i,j) = 0.d0
            do k=1,Nvib
                gBder(i,j) = gBder(i,j) + Bder(k,i,j)*Grad(k)
            enddo
        enddo
        enddo
        if (verbose>2) then
            print*, "Correction matrix to be applied on Hx:"
            call MAT0(6,gBder,3*Nat,3*Nat,"gs*Bder matrix")
        endif

        if (check_symmetry) then
            call check_symm_gsBder(state1,gBder)
        endif

        ! Apply term to Hess
        Hess(1:3*Nat,1:3*Nat) = Hess(1:3*Nat,1:3*Nat) - gBder(1:3*Nat,1:3*Nat)

    endif

    ! Get Hessian in internal coordinates and compute modes
    Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,P,Hess,counter=.true.)
    call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B1,G1)
    call gf_method(Nvib,Nvib1,G1,Hess,L1,Freq1,X,X1inv)

    ! If only one state is give, exit now
    if (adjustl(inpfile2) == "none") stop

    !===========
    !State 2
    !===========
    if (verbose>0) then
        print*, ""
        print*, "=========="
        print*, " STATE 2"
        print*, "=========="
    endif
 
    ! READ DATA (each element from a different file is possible)
    ! ---------------------------------
    !Guess filetypes
    if (ft2 == "guess") &
    call split_line_back(inpfile2,".",null,ft2)
    if (fth2 == "guess") &
    call split_line_back(hessfile2,".",null,fth2)
    if (ftg2 == "guess") &
    call split_line_back(gradfile2,".",null,ftg2)

    ! STRUCTURE FILE
    open(I_INP,file=inpfile2,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile2)) )
    call generic_strmol_reader(I_INP,ft2,state2,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (State2)")
    close(I_INP)
    ! Check that the structure is that of the State1 if vertical is used
    if (vertical) then
        call set_geom_units(state1,"Angs")
        Theta2=0.d0
        Theta3=0.d0
        do i=1,state1%natoms
            Theta  = calc_atm_dist(state1%atom(i),state2%atom(i))
            Theta2 = Theta2 + (Theta)**2
            Theta3 = max(Theta3,Theta)
        enddo
        Theta = sqrt(Theta2/state1%natoms)
        if (Theta3>1.d-4) then
            print'(X,A,X,F8.3)  ', "RMSD_struct (AA):", Theta
            print'(X,A,X,E10.3,/)', "Max Dist", Theta3
            call alert_msg("warning","vertical model is requested but State1 and State2 do not have the same structure")
        endif
        call set_geom_units(state1,"Bohr")
    endif
    ! Shortcuts
    if (Nat /= state2%natoms) call alert_msg("fatal","Initial and final states don't have the same number of atoms.")

    ! Read mass from file if given
    if (adjustl(mass_file) /= "none") then
        print'(/,X,A)', "Reading atomic masses from: "//trim(adjustl(mass_file))
        open(I_MAS,file=mass_file,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(mass_file)) )
        do i=1,Nat
            read(I_MAS,*,iostat=IOstatus) state2%atom(i)%mass 
            if (IOstatus /= 0) call alert_msg( "fatal","While reading "//trim(adjustl(mass_file)) )
        enddo
        close(I_MAS)
    endif

    ! HESSIAN FILE
    open(I_INP,file=hessfile2,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile2)) )
    allocate(A(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fth2,Nat,A,error) 
    if (error /= 0) call alert_msg("fatal","Error reading Hessian (State2)")
    close(I_INP)
    
    ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
    ! This is also used to check if the symmetry changed from the other state
    call vibrations_Cart(Nat,state2%atom(:)%X,state2%atom(:)%Y,state2%atom(:)%Z,state2%atom(:)%Mass,A,&
                         NNvib,L2,Freq2,error_flag=error)
    if (NNvib /= 3*Nat-6) call alert_msg("warning","Linear molecule (at state2). Things can go very bad.")
    if (NNvib > Nvib .and. .not.apply_projection_matrix) then
        print'(/,A,/)', "Using reduced space from State1"
        NNvib = Nvib
    endif
    k=0
    do i=1,3*Nat
    do j=1,i
        k=k+1
        Hess(i,j) = A(k)
        Hess(j,i) = A(k)
    enddo 
    enddo
    deallocate(A)

    if (vertical.or.gradcorrectS2) then
        ! GRADIENT FILE
        open(I_INP,file=gradfile2,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile2)) )
        call generic_gradient_reader(I_INP,ftg2,Nat,Grad,error)
        if (error /= 0) call alert_msg("fatal","Error reading gradient (State2)")
        close(I_INP)
    endif
    
    
    !----------------------------------
    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! General Actions:
    !******************
    ! Get connectivity 
    if (cnx_file == "guess") then
        call guess_connect(state2)
    else
        print'(/,A,/)', "Reading connectivity from file: "//trim(adjustl(cnx_file))
        open(I_CNX,file=cnx_file,status='old')
        call read_connect(I_CNX,state2)
        close(I_CNX)
    endif

    ! Manage symmetry
    if (.not.use_symmetry) then
        state2%PG="C1"
    else if (trim(adjustl(symm_file)) /= "none") then
        call alert_msg("note","Using custom symmetry file: "//trim(adjustl(symm_file)) )
        open(I_SYM,file=symm_file)
        do i=1,state2%natoms
            read(I_SYM,*) j, isym(j)
        enddo
        close(I_SYM)
        !Set PG to CUStom
        state2%PG="CUS"
    else
        state2%PG="XX"
        call symm_atoms(state2,isym)
    endif

    call set_geom_units(state2,"BOHR")

    !--------------------------------
    ! Compute projection matrix
    !--------------------------------
    call heading(6,"COMPUTING PROJECTION MATRIX (CARTESIAN)")

    !---------------------------------------
    call subheading(6,"Internal Sets")

    call statement(6,"Set with ALL internal coordinates from connectivity",keep_case=.true.)
    call gen_bonded(state2)
    call define_internal_set(state2,"ALL",intfile,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Fltr)
    allgeom = state2%geom

    call statement(6,"Set indicated on INPUT: "//trim(def_internal),keep_case=.true.)
    call gen_bonded(state2)
    call define_internal_set(state2,def_internal,intfile,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Fltr)
    ! If not using combinations, the Filter need to be reset to the identity matrix
    if (Nf==0) then
        Nf=Ns
        Fltr(1:Nf,1:Ns) = identity_matrix(Nf)
    endif
    inputgeom = state2%geom
    !---------------------------------------


    ! Compute B matrix in the inputgeom for projection
    call subheading(6,"Compute B with INPUT set")

    ! Compute B and apply Filter
    call internal_Wilson(state2,Ns,S2,Bprj)
    Bprj(1:Nf,1:3*Nat) = matrix_product(Nf,3*Nat,Ns,Fltr,Bprj)
    Ns=Nf

    NvibP = Nvib
    !*********************************************
    ! Really need to do redundant2nonredundant?
    ! Compute G
    call internal_Gmetric(Nat,Ns,state2%atom(:)%mass,Bprj,G2)
    call redundant2nonredundant(Ns,NvibP,G2,Asel2)
    Bprj(1:NvibP,1:3*Nat) = matrix_product(NvibP,3*Nat,Ns,Asel2,Bprj,tA=.true.)
    !*********************************************

    ! Compute projecton matrix. Method3 provide the matrix to be applied to
    ! Cartesian Hessian/Grad, but defined on MWC
    call statement(6,"Computing projection matrix")
    P(1:3*Nat,1:3*Nat) = projection_matrix3(Nat,NvibP,Bprj,state2%atom(:)%Mass)
    if (complementay_projection) then
        Aux(1:3*Nat,1:3*Nat) = identity_matrix(3*Nat)
        P(1:3*Nat,1:3*Nat) =  Aux(1:3*Nat,1:3*Nat)-P(1:3*Nat,1:3*Nat)
    endif


    ! Get back allgeom
    call heading(6,"Computing modes with ALL set")
    Ns = allgeom%nbonds+allgeom%nangles+allgeom%ndihed
    state2%geom = allgeom


    ! Rotate gradient and Hessian
    call subheading(6,"Modes projecting the Hess/Grad")
    call statement(6,"Rotating Gradient and Hessian")
!     ! Get full Cartesian Hessian from Hlt
!     Hess(1:3*Nat,1:3*Nat) = Hlt_to_Hess(3*Nat,Hlt)
!     Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,P,Hess,counter=.true.)
    ! Get Cartesian gradient
!     Grad(1:3*Nat) = grdx
!     Grad(1:3*Nat) = matrix_vector_product(3*Nat,3*Nat,P,Grad)

    ! Compute B, G and, if needed, Bder
    call statement(6,"Computing B, G (with ALL set)",keep_case=.true.)
    call internal_Wilson(state2,Ns,S2,B2,ModeDef)
    call internal_Gmetric(Nat,Ns,state2%atom(:)%mass,B2,G2)
    if (gradcorrectS2) then
        call statement(6,"...and Bder for non-stationary point calculation")
        call calc_Bder(state2,Ns,Bder,analytic_Bder)
    endif

    ! Get non-redundant space
    call subsubheading(6,"Getting the actual vibrational space dimension")
    call redundant2nonredundant(Ns,Nvib,G2,Asel2)
    call statement(6,"Rotate B, G to non-redundant space", keep_case=.true.)
    ! Rotate Bmatrix
    B2(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel2,B2,tA=.true.)
    ! Rotate Gmatrix
    G2(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel2(1:Ns,1:Nvib),G2,counter=.true.)
    if (gradcorrectS2) then
        call statement(6,"Also rotate Bder to non-redundant space", keep_case=.true.)
        ! Rotate Bders
        do j=1,3*Nat
            Bder(1:Nvib,j,1:3*Nat) =  matrix_product(Nvib,3*Nat,Ns,Asel2,Bder(1:Ns,j,1:3*Nat),tA=.true.)
        enddo

        ! Get the Correction now
        call statement(6," Getting gs^t\beta term")
        ! The correction is applied with the Nvib0 SET
        ! Correct Hessian as
        ! Hx' = Hx - gs^t\beta
        ! 1. Get gs from gx
        call Gradcart2int(Nat,Nvib,Grad,state2%atom(:)%mass,B2,G2)
        ! 2. Multiply gs^t\beta and
        ! 3. Apply the correction
        ! Bder(i,j,K)^t * gq(K)
        do i=1,3*Nat
        do j=1,3*Nat
            gBder(i,j) = 0.d0
            do k=1,Nvib
                gBder(i,j) = gBder(i,j) + Bder(k,i,j)*Grad(k)
            enddo
        enddo
        enddo
        if (verbose>2) then
            print*, "Correction matrix to be applied on Hx:"
            call MAT0(6,gBder,3*Nat,3*Nat,"gs*Bder matrix")
        endif

        if (check_symmetry) then
            call check_symm_gsBder(state2,gBder)
        endif

        ! Apply term to Hess
        Hess(1:3*Nat,1:3*Nat) = Hess(1:3*Nat,1:3*Nat) - gBder(1:3*Nat,1:3*Nat)

    endif

    ! Get Hessian in internal coordinates and compute modes
    Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,P,Hess,counter=.true.)
    call HessianCart2int(Nat,Nvib,Hess,state2%atom(:)%mass,B2,G2)
    call gf_method(Nvib,Nvib2,G2,Hess,L2,Freq2,X,X2inv)
    
    if (Nvib1 /= Nvib2) then
        call alert_msg("fatal","State1 and State2 have different number of normal modes")
    endif
    NvibP = Nvib1
!     NvibP = Nvib
    
    !==========================================
    ! CHECKS ON THE INTERNAL SETS
    !==========================================
    ! Evaluate orthogonality
    if (verbose>0) then
     print'(2/,X,A)', "============================================"
     print*,          " Internal Coordinates Orthogonality Checks  "
     print*,          "============================================"
     print*,          "Analysing: D = G1^-1/2 (A1^-1 A2) G2^1/2"
    endif
    
    same_red2nonred_rotation=.false.
    if (.not.same_red2nonred_rotation) then
        ! We need to include the fact that Asel1 /= Asel2, i.e.,
        ! rotate to the same redundant space L(red) = Asel * L(non-red)
        ! J = L1^-1 A1^-1 A2 L2
        ! so store in Aux the following part: [L1^-1 A1^-1 A2]
        Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Ns,Asel1inv,Asel2)
        Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,X1inv,Aux)
    else
        Aux(1:Nvib,1:Nvib) = X1inv(1:Nvib,1:Nvib)
    endif
    
    Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,Aux,X)
    open (O_DMAT,file="D_matrix_abs.dat",status="replace")
    do i=1,Nvib
        write(O_DMAT,'(600f8.2)') dabs(Aux(i,1:Nvib))
    enddo
    close(O_DMAT)
    open (O_DMAT,file="D_matrix.dat",status="replace")
    do i=1,Nvib
        write(O_DMAT,'(600f8.2)') Aux(i,1:Nvib)
    enddo
    close(O_DMAT)
    if (verbose>0) &
     print'(X,A,/)', "(D matrix written to files: D_matrix.dat and D_matrix_abs.dat)"
    
    if (verbose>0) then
        ! COMPUTE DETERMINANT AND TRACE OF D_matrix
        theta = 0.d0
        do i=1,Nvib
            theta = theta+Aux(i,i)
        enddo    
        print'(X,A,F8.2,A,I0,A)', "Trace", theta, "  (",Nvib,")"
        theta = -100.d0  !max
        theta2 = 100.d0  !min
        do i=1,Nvib
            if (Aux(i,i) > theta) then
                theta = Aux(i,i)
                imax = i
            endif
            if (Aux(i,i) < theta2) then
                theta2 = Aux(i,i)
                imin = i
            endif
        enddo 
        print*, "Min diagonal", theta2, trim(adjustl(ModeDef(imin)))
        print*, "Max diagonal", theta,  trim(adjustl(ModeDef(imax)))
        theta = -100.d0  !max
        theta2 = 100.d0  !min
        theta3 = 0.d0
        k=0
        do i=1,Nvib
        do j=1,Nvib
            if (i==j) cycle
            k=k+1
            theta3 = theta3 + dabs(Aux(i,j))
            theta  = max(theta, Aux(i,j))
            theta2 = min(theta2,Aux(i,j))
        enddo
        enddo 
        print*, "Min off-diagonal", theta2    
        print*, "Max off-diagonal", theta
        print*, "AbsDev off-diagonal sum", theta3
        print*, "AbsDev off-diagonal per element", theta3/dfloat(k)
        theta = determinant_realgen(Nvib,Aux)
        print*, "Determinant", theta
        print*, "-------------------------------------------------------"
        print*, ""
    endif
    
    
    ! ===========================================
    ! COMPUTE DUSCHINSKY MATRIX AND DISPLACEMENT VECTOR 
    ! ===========================================
    !--------------------------
    ! J-matrix (Duschinski)
    !--------------------------
    ! Orthogonal Duschinski (from orthogonalized ICs. Not used in this version)
    ! Get orthogonal modes:  L' = G^-1/2 L
    Aux(1:Nvib,1:NvibP)  = matrix_product(Nvib,NvibP,Nvib,X1inv,L1)
    Aux2(1:Nvib,1:NvibP) = matrix_product(Nvib,NvibP,Nvib,X2inv,L2)
    ! Duschinsky matrix (orth) stored in JdusO = L1'^t L2'
!     JdusO(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,Aux,Aux2,tA=.true.)
    !Store L1' in Aux2 to later be used to get the displacement
    Aux2(1:Nvib,1:NvibP)=Aux(1:Nvib,1:NvibP)
    
    ! Non-Orthogonal Duschinski (the one we use)
    if (verbose>0) &
     print*, "Calculating Duschisky..."
    !Inverse of L1 (and store in L1inv)
!     L1inv(1:Nvib,1:Nvib) = inverse_realgen(Nvib,L1(1:Nvib,1:Nvib))
    ! Generalized inverse: (L^t L)^-1 L^
    Aux(1:NvibP,1:NvibP) = matrix_product(NvibP,NvibP,Nvib,L1,L1,tA=.true.)
    Aux(1:NvibP,1:NvibP) = inverse_realgen(NvibP,Aux)
    L1inv(1:NvibP,1:Nvib)  = matrix_product(NvibP,Nvib,NvibP,Aux,L1,tB=.true.) 

    ! Account for different rotations to non-redundant set 
    ! but preserve the inverse L1 matrix in L1
! always do redundant2nonredundant
!     if (Nvib<Ns .and. .not.same_red2nonred_rotation) then
    if (.not.same_red2nonred_rotation) then
        print*, "  Using different A rotations"
        ! We need to include the fact that Asel1 /= Asel2, i.e.,
        ! rotate to the same redundant space L(red) = Asel * L(non-red)
        ! J = L1^-1 A1^-1 A2 L2
        ! so store in Aux the following part: [L1^-1 A1^-1 A2]
        Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Ns,Asel1inv,Asel2)
        Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,L1inv,Aux)
    else
        print*, "  Using the same A rotation for both states"
        Aux(1:NvibP,1:Nvib) = L1inv(1:NvibP,1:Nvib)
    endif
    !J = L1^-1 [A1^t A2] L2 (stored in J).
    Jdus(1:NvibP,1:NvibP) = matrix_product(NvibP,NvibP,Nvib,Aux,L2)


    !--------------------------
    ! ICs displacement
    !--------------------------
    if (verbose>0) then
        print*, ""
        print*, "=========================="
        print*, " SHIFTS (internal coord)"
        print*, "=========================="
    endif
    if (vertical) then
        print*, "Vertical model"
        ! GET MINIMUM IN INTERNAL COORDINATES
        ! At this point
        !  * Hess has the Hessian  of State2 in internal coords (output from HessianCart2int)
        !  * Grad has the gradient of State2 in internal coords (output from HessianCart2int)
        ! Inverse of the Hessian
        ! In order to support projected matrix where Nvib is reduced
!         Aux(1:Nvib0,1:Nvib0) = inverse_realgen(Nvib0,Hess)
        call generalized_inv(Nvib,Nvib1,Hess,Aux)
        if (Nvib1 /= NvibP) then
            call alert_msg("warning","Non-zero eigenvales from generalized inverse do not match Nvib")
        endif
        ! DeltaS0 = -Hs^1 * gs
        do i=1,Nvib
            Delta(i)=0.d0
            do k=1,Nvib0
                Delta(i) = Delta(i)-Aux(i,k) * Grad(k)
            enddo 
        enddo
! always do redundant2nonredundant
!         if (Nvib<Ns) then
            !Transform Delta' into Delta (for control purposes only)
            ! Delta = A Delta'
            do i=1,Ns
                Vec1(i) = 0.d0
                do k=1,Nvib
                    Vec1(i) = Vec1(i) + Asel1(i,k)*Delta(k)
                enddo
            enddo
!         else
!             Vec1(1:Nvib)=Delta(1:Nvib)
!         endif

        ! Get coordinates
        do i=1,Ns
            S2(i) = S1(i) + Vec1(i)
        enddo

    else ! Adiabatic case
        print*, "Adiabatic model"
        ! Bonds
        do i=1,state1%geom%nbonds
            Delta(i) = S2(i)-S1(i)
        enddo
        ! Angles
        do j=i,i+state1%geom%nangles-1
            Delta(j) = S2(j)-S1(j)
        enddo
        ! Dihedrals
        do k=j,j+state1%geom%ndihed-1
            Delta(k) = S2(k)-S1(k)
            Delta_p = S2(k)-S1(k)+2.d0*PI
            if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
            Delta_p = S2(k)-S1(k)-2.d0*PI
            if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
        enddo
        ! Impropers
        do i=k,k+state1%geom%nimprop-1
            Delta(i) = S2(i)-S1(i)
            Delta_p = S2(i)-S1(i)+2.d0*PI
            if (dabs(Delta_p) < dabs(Delta(i))) Delta(i)=Delta_p
            Delta_p = S2(i)-S1(i)-2.d0*PI
            if (dabs(Delta_p) < dabs(Delta(i))) Delta(i)=Delta_p
        enddo

        ! Store Deltas into the auxiliar vector Vec1
        Vec1(1:Ns) = Delta(1:Ns)
! always do redundant2nonredundant
!         if (Nvib<Ns) then
            ! Need to transform Deltas into the non-redundant coordinates set
            ! Delta' = A^-1 Delta
            do i=1,Nvib
                Vec2(i) = 0.d0
                do k=1,Ns
                    Vec2(i) = Vec2(i) + Asel1inv(i,k)*Delta(k)
                enddo
            enddo
            Delta(1:Nvib) = Vec2(1:Nvib)
!         endif

    endif

    if (verbose>0) then
        print*, "List of ICs and Deltas"
        if (state1%geom%nbonds /= 0) &
         print'(X,A,/,X,A)', "Bonds (Angs)",&
                             "  IC   Description      State2    State1     Delta"
        do i=1,state1%geom%nbonds
            print'(I5,X,A,X,3(F8.3,2X))', i,trim(ModeDef(i)),S2(i)*BOHRtoANGS,S1(i)*BOHRtoANGS,Vec1(i)*BOHRtoANGS
        enddo
        if (state1%geom%nangles /= 0) &
         print'(X,A,/,X,A)', "Angles (deg)",&
                             "  IC   Description                State2    State1     Delta"
        do j=i,i+state1%geom%nangles-1
            print'(I5,X,A,X,3(F8.2,2X))', j, trim(ModeDef(j)), S2(j)*180.d0/PI,S1(j)*180.d0/PI,Vec1(j)*180.d0/PI
        enddo
        if (state1%geom%ndihed /= 0) &
         print'(X,A,/,X,A)', "Dihedrals (deg)",&
                             "  IC   Description                          State2    State1     Delta"
        do k=j,j+state1%geom%ndihed-1
            print'(I5,X,A,X,3(F8.2,2X))', k, trim(ModeDef(k)), S2(k)*180.d0/PI,S1(k)*180.d0/PI,Vec1(k)*180.d0/PI
        enddo
        if (state1%geom%nimprop /= 0) &
         print'(X,A,/,X,A)', "Impropers (deg)",&
                             "  IC   Description                          State2    State1     Delta"
        do i=k,k+state1%geom%nimprop-1
            print'(I5,X,A,X,3(F8.2,2X))', i, trim(ModeDef(i)), S2(i)*180.d0/PI,S1(i)*180.d0/PI,Vec1(i)*180.d0/PI
        enddo
    endif


    !--------------------------
    ! K-vector (NM-shifts)
    !--------------------------
    if (verticalQspace2) then 
        ! Apply vertical model in the Qspace (get S2 equilibrium --Delta-- within the normal mode coordinate space)
        ! 
        ! K = -J * Lambda_f^-1 * L2^t * gs
        ! Convert Freq into FC. Store in FC for future use
        do i=1,NvibP
            FC(i) = sign((Freq2(i)*2.d0*pi*clight*1.d2)**2/HARTtoJ*BOHRtoM**2*AUtoKG,Freq2(i))
        enddo
        ! Lambda_f^-1 * L2^t
        do i=1,NvibP
            Aux(i,1:Nvib) = L2(1:Nvib,i) / FC(i)
        enddo
        ! -[Lambda_f^-1 * L2^t] * gs
        do i=1,NvibP
            Q0(i)=0.d0
            do k=1,Nvib
                Q0(i) = Q0(i) - Aux(i,k) * Grad(k)
            enddo
            ! Change imag to real if requested. This is done now, once the shift was computed
            ! So the displacement is anyway computed towards the stationary point of the quadratic PES
            ! It would be equivalent to also change the gradient if we did the change imag to real before 
            ! this point (hence the warning message)
            if (FC(i)<0) then
                print*, i, FC(i), Freq2(i)
                if (force_real) then 
                    FC(i)    = abs(FC(i))
                    Freq2(i) = abs(Freq2(i))
                    call alert_msg("warning","Negative FC turned real (Gradient also changed)")
                else
                    call alert_msg("warning","A negative FC found")
                endif
            endif
        enddo
        ! J * [-Lambda_f^-1 * L2^t * gs]
        do i=1,NvibP
            Vec1(i)=0.d0
            do k=1,NvibP
                Vec1(i) = Vec1(i) + Jdus(i,k) * Q0(k)
            enddo
        enddo

    elseif (verticalQspace1) then 
        ! Use the same strategy as in the Cartesian case: 1) normal modes at State1; 2) Rotate normal modes
        ! At this point
        !  * Hess has the Hessian  of State2 in internal coords (output from HessianCart2int)
        !  * Grad has the gradient of State2 in internal coords (output from HessianCart2int)
        !
        ! Rotate to Q1-space (IC)
        !
        ! HESIAN
        !  H_Q' = L^t Hs L (also H_Q' = L^-1 G Hs L)
        ! Note:
        !  * L1 contains the normal matrix (now the inverse is in L1inv)
        Hess(1:NvibP,1:NvibP) = matrix_basisrot(NvibP,Nvib,L1,Hess,counter=.true.)
        !
        ! GRADIENT
        !  g_Q' = L^t gs
        do i=1,NvibP
            Vec1(i) = 0.d0
            do k=1,Nvib
                Vec1(i) = Vec1(i) + Aux2(k,i) * Grad(k)
            enddo
        enddo
        Grad(1:NvibP) = Vec1

        ! Diagonalize Hessian in Q1-space to get State2 FC and Duschinski rotation
        call diagonalize_full(Hess(1:NvibP,1:NvibP),Nvib,Jdus(1:NvibP,1:NvibP),FC(1:NvibP),"lapack")
        Freq2(1:NvibP) = FC2Freq(NvibP,FC)
        call print_vector(6,Freq2,NvibP,"Frequencies (cm-1) -- from Q1-space")

        ! Get shift vector (also compute Qo'')
        ! First compute Qo''
        ! Q0 = - FC^-1 * J^t * gQ
        do i=1,NvibP
            Q0(i) = 0.d0
            do k=1,NvibP
                Q0(i) = Q0(i) - Jdus(k,i) * Grad(k)
            enddo
            Q0(i) = Q0(i) / FC(i)
            ! Change imag to real if requested. This is done now, once the shift was computed
            ! So the displacement is anyway computed towards the stationary point of the quadratic PES
            ! It would be equivalent to also change the gradient if we did the change imag to real before 
            ! this point (hence the warning message)
            if (FC(i)<0) then
                print*, i, FC(i), Freq2(i)
                if (force_real) then 
                    FC(i)    = abs(FC(i))
                    Freq2(i) = abs(Freq2(i))
                    call alert_msg("warning","Negative FC turned real (Gradient also changed)")
                else
                    call alert_msg("warning","A negative FC found")
                endif
             endif
        enddo
        if (verbose>2) then
            call print_vector(6,FC*1e5,Nvib,"FC - int")
            call print_vector(6,Grad*1e5,Nvib,"Grad - int")
            call print_vector(6,Q0,Nvib,"Q0 - int")
        endif
        ! K = J * Q0
        do i=1,NvibP
            Vec1(i) = 0.d0
            do k=1,NvibP
                Vec1(i) = Vec1(i) + Jdus(i,k) * Q0(k)
            enddo
        enddo

    else
        ! K = L1^-1 DeltaS (this is State 1 respect to state 2) . 
        ! Notes
        !   * L1inv stores the inverse
        !   * Delta in non-redundant IC set
        do i=1,NvibP
            Vec1(i) = 0.d0
            do k=1,Nvib
                Vec1(i) = Vec1(i) + L1inv(i,k)*Delta(k)
            enddo
        enddo
!         !Orthogonal: K=L1'^t DeltaS'
!         do i=1,Nvib
!             Vec2(i) = 0.d0
!             do k=1,Nvib
!                 Vec2(i) = Vec1(i) + Aux2(k,i)*Delta(k)
!             enddo
!         enddo
    endif

    if (verbose>2) then
        call MAT0(6,Jdus,NvibP,NvibP,"DUSCHINSKI MATRIX")
        call print_vector(6,Vec1,NvibP,"NORMAL MODE SHIFT")
    endif

    !Analyze Duschinsky matrix
    call analyze_duschinsky(6,NvibP,Jdus,Vec1,Freq1,Freq2)


    !===================================
    ! Reorganization energy
    !===================================
    print'(/,X,A)', "REORGANIZATION ENERGY"
    if (verticalQspace2) then
        print*, "Vertical model / Q2space"
        ! Normal-mode space
        ! Er = -L2^t gs * Q0 - 1/2 * Q0^t * Lambda_f * Q0
        ! At this point: 
        ! * Grad: in internal coords
        ! * Q0: DeltaQ in final state nm
        ! * FC: diagonal force constants for final state
        Er = 0.d0
        do i=1,NvibP
            ! Compute gQ(i) = L2^t * gs
            Theta = 0.d0
            do j=1,Nvib
                Theta =  Theta + L2(j,i)*Grad(j)
            enddo
            Er = Er - Theta * Q0(i) - 0.5d0 * FC(i) * Q0(i)**2
        enddo
    elseif (verticalQspace1) then
        print*, "Vertical model / Q1space"
        ! Normal-mode space
        ! Er = -L2^t gs * Q0 - 1/2 * Q0^t * Lambda_f * Q0
        ! At this point: 
        ! * Grad: in Q1-space
        ! * Q0: DeltaQ in final state nm
        ! * FC: diagonal force constants for final state
        Er = 0.d0
        do i=1,NvibP
            ! Get gradient in state2 Qspace
            ! gQ2 = J^t * gQ1
            Theta=0.d0
            do k=1,NvibP
                Theta = Theta + Jdus(k,i) * Grad(k)
            enddo
            Er = Er - Theta * Q0(i) - 0.5d0 * FC(i) * Q0(i)**2
        enddo
    elseif (vertical) then
        print*, "Vertical model / IC space"
        ! Internal-coordinates space
        ! Er = -gs * DeltaS - 1/2 DeltaS^t * Hs * DeltaS
        ! At this point (all data in the non-redundant IC space)
        ! * Grad: in internal coords
        ! * Delta: DeltaS 
        ! * Hess: Hessian in internal coords
        !
        ! Fisrt, compute DeltaS^t * Hs * DeltaS
        Theta=0.d0
        do j=1,Nvib
        do k=1,Nvib
            Theta = Theta + Delta(j)*Delta(k)*Hess(j,k)
        enddo
        enddo
        Er = -Theta*0.5d0
        do i=1,Nvib
            Er = Er - Grad(i)*Delta(i)
        enddo
    else ! Adiabatic
        print*, "Adiabatic model / IC space"
        ! Er = 1/2 DeltaS^t * Hs * DeltaS
        ! At this point (all data in the non-redundant IC space)
        ! * Grad: in internal coords
        ! * Delta: DeltaS 
        ! * Hess: Hessian in internal coords
        !
        ! Fisrt, compute DeltaS^t * Hs * DeltaS
        Theta=0.d0
        do j=1,Nvib
        do k=1,Nvib
            Theta = Theta + Delta(j)*Delta(k)*Hess(j,k)
        enddo
        enddo
        Er = Theta*0.5d0
    endif
    print'(X,A,F12.6)',   "Reorganization energy (AU) = ", Er
    print'(X,A,F12.6,/)', "Reorganization energy (eV) = ", Er*HtoeV


    ! ====================
    ! Print state files (better compute the dipole derivs directly in the Q-space) This also requires a change in FCclasses
    ! otherwise, the result will be approx, because FCclasses will first orthogonalize L1 and L2...
    ! ====================
    ! The L matrix obtained from Ls_to_Lcart is consistent with the  Cartesia geom of the state used. 
    ! In case the states are rotated, the above transformation should only be done on one of the states
    ! and the other must be obtained rotating the one transformed. Note that transformed one would be 
    ! consistent with the input geometries, and therefore, with the dipole derivatives
    if (reference_frame=="I".or.vertical) then
        print*, "Reference state to report statefiles: Initial"
        !HTi
        ! * L1: from Ls_to_Lcart
        call Ls_to_Lcart(Nat,Nvib,state1%atom(:)%mass,B1,G1,L1,L1,error)
        ! * L2: from rotation (Duschinski) of L1,   L2 = L1*J
        L2(1:3*Nat,1:Nvib) = matrix_product(3*Nat,Nvib,Nvib,L1,Jdus)
    else ! HTmode='F'
        print*, "Reference state to report statefiles: Final"
        ! HTf
        !  * L2: from Ls_to_Lcart
        call Ls_to_Lcart(Nat,Nvib,state2%atom(:)%mass,B2,G2,L2,L2,error)
        !  * L1: from rotation (Duschinski) of L2   L1 = L2*J^-1
        ! Need inverse of J
        Aux(1:Nvib,1:Nvib) = inverse_realgen(NvibP,Jdus)
        L1(1:3*Nat,1:Nvib) = matrix_product(3*Nat,Nvib,Nvib,L2,Aux)
    endif
    ! STATE1
    ! Compute new state_file
    ! T1(g09) = mu^1/2 m B^t G1^-1 L1
    !Print state
    open(O_STAT,file="state_file_1")
    call set_geom_units(state1,"Angs")
    do i=1,Nat
        write(O_STAT,'(E17.8)') state1%atom(i)%x
        write(O_STAT,'(E17.8)') state1%atom(i)%y
        write(O_STAT,'(E17.8)') state1%atom(i)%z
    enddo
    call Lcart_to_LcartNrm(Nat,Nvib,L1,Aux,error)
    do i=1,3*Nat
    do j=1,NvibP
        write(O_STAT,'(E17.8)') Aux(i,j)
    enddo
    enddo
    do j=1,NvibP
        if (force_real.and.Freq1(j)<0) then
            print*, Freq1(j)
            call alert_msg("warning","An imagainary frequency turned real (state1)")
            Freq1(j) = abs(Freq1(j))
        endif
        write(O_STAT,'(F10.4)') Freq1(j)
    enddo
    close(O_STAT)
    ! STATE2
    if (vertical) then
        ! Deactivate state2 coords to avoid confusion
        state2%atom(1:Nat)%x=0.d0
        state2%atom(1:Nat)%y=0.d0
        state2%atom(1:Nat)%z=0.d0
    endif
    !Print state
    open(O_STAT,file="state_file_2")
    call set_geom_units(state2,"Angs")
    ! Note that the geometry is the input one (not displaced for vertical)
    ! But it is ok for FCclasses (it is not using it AFIK) What about HT??
    do i=1,Nat
        write(O_STAT,'(E17.8)') state2%atom(i)%x
        write(O_STAT,'(E17.8)') state2%atom(i)%y
        write(O_STAT,'(E17.8)') state2%atom(i)%z
    enddo
    call Lcart_to_LcartNrm(Nat,Nvib,L2,Aux,error)
    do i=1,3*Nat
    do j=1,NvibP
        write(O_STAT,'(E17.8)') Aux(i,j)
    enddo
    enddo
    do j=1,NvibP
        if (force_real.and.Freq2(j)<0) then
            print*, Freq2(j)
            call alert_msg("warning","An imagainary frequency turned real (state2)")
            Freq2(j) = abs(Freq2(j))
        endif
        write(O_STAT,'(F10.4)') Freq2(j)
    enddo
    close(O_STAT)

    !============================================
    ! PRINT DUSCHINSKI AND DISPLACEMENT TO FILES
    !============================================
    print*, "Printing Duschinski matrix to 'duschinsky.dat'"
    open(O_DUS, file="duschinsky.dat")
!     open(O_DUS2,file="duschinsky_orth.dat")
    print'(X,A,/)', "Printing Shift vector to 'displacement.dat'"
    open(O_DIS, file="displacement.dat")
!     open(O_DIS2,file="displacement_orth.dat")
    do i=1,NvibP
    do j=1,NvibP
        write(O_DUS,*)  Jdus(i,j)
!         write(O_DUS2,*) JdusO(i,j)
    enddo 
        write(O_DIS,*)  Vec1(i)
!         write(O_DIS2,*) Vec2(i)
    enddo
    close(O_DUS)
!     close(O_DUS2)
    close(O_DIS)
!     close(O_DIS2)

    if (get_distances) then
        print*, "DISTANCES FROM STATE1 TO STATE2 IN Q1-SPACE"
        ! Euclidean distance
        dist=0.d0
        do i=1,NvibP
            dist = dist + Vec1(i)**2
        enddo
        dist=dsqrt(dist)
        print'(X,A,F10.4)', "Euclidean distance in nm space", dist/dsqrt(AMUtoAU)
        
        ! RC path distance
        FC(1:NvibP) = Freq2FC(NvibP,Freq1)
        dist=0.d0
        area=1.d0
        dt=5.d2
        time = 0.d0
        Nvib0=NvibP
        k=0
        do while (dabs(area) > 1d-10 .and. Nvib0>0)
            f0=0.d0
            f1=0.d0
            Nvib0=NvibP
            j = 0
            do i=1,NvibP
                ff = FC(i)**2*Vec1(i)**2*dexp(-2.d0*dabs(FC(i))*time)
                f0 = f0 + ff
                f1 = f1 + FC(i)**2*Vec1(i)**2*dexp(-2.d0*dabs(FC(i))*(time+dt))
                ! Discard modes that reached the baseline
                if (ff < 5e-24) then
                    Nvib0=Nvib0-1
                else
                    j = j+1
                    Vec2(j) = Vec1(i)
                endif
            enddo
            NvibP = Nvib0
            Vec1(1:NvibP) = Vec2(1:NvibP)
            f0 = dsqrt(f0)
            f1 = dsqrt(f1)
            area = 0.5d0*(f0+f1)*dt
            dist = dist + area
            time=time+dt
            k=k+1
            if (k==100000000) exit
        enddo
        if (k<100000000) then
            print'(X,A,F10.4,/)', "Contour distance in IRC space ", dist/dsqrt(AMUtoAU)
        else
            print'(X,A,F10.4,/)', "Contour distance in IRC space (integration failed)"
        endif
    endif

    call summary_alerts

    call cpu_time(tf)
    write(6,'(/,A,X,F12.3,/)') "CPU (s) for internal vib analysis: ", tf-ti

    stop

    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,&
                           inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                           intfile,intfile0,rmzfile,def_internal,def_internal0,&
                           use_symmetry,cnx_file,mass_file,&
!                          tswitch,
                           symaddapt,same_red2nonred_rotation,analytic_Bder,&
                           model,vertical,verticalQspace2,verticalQspace1,&
                           gradcorrectS1,gradcorrectS2,&
                           orthogonalize,original_internal,force_real,reference_frame,&
                           apply_projection_matrix,complementay_projection,move_to_min)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,&
                                          inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                                          intfile,intfile0,rmzfile,def_internal,def_internal0,&
                                          cnx_file,reference_frame,model,mass_file !, symfile
        logical,intent(inout)          :: use_symmetry, vertical, verticalQspace2, &
                                          verticalQspace1, &
                                          gradcorrectS1, gradcorrectS2, symaddapt, &
                                          same_red2nonred_rotation,analytic_Bder, &
                                          orthogonalize,original_internal,force_real, &
                                          apply_projection_matrix,move_to_min,complementay_projection
!         logical,intent(inout) :: tswitch

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        character(len=500) :: input_command
        character(len=10)  :: MODEL_UPPER
        ! iargc type must be specified with implicit none (strict compilation)
        integer :: iargc

        ! Tune defaults
        logical :: gradcorrectS1_default=.true., &
                   gradcorrectS2_default=.true.

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

                case ("-f2") 
                    call getarg(i+1, inpfile2)
                    argument_retrieved=.true.
                case ("-ft2") 
                    call getarg(i+1, ft2)
                    argument_retrieved=.true.

                case ("-fhess2") 
                    call getarg(i+1, hessfile2)
                    argument_retrieved=.true.
                case ("-fth2") 
                    call getarg(i+1, fth2)
                    argument_retrieved=.true.

                case ("-fgrad2") 
                    call getarg(i+1, gradfile2)
                    argument_retrieved=.true.
                case ("-ftg2") 
                    call getarg(i+1, ftg2)
                    argument_retrieved=.true.

                case ("-fmass") 
                    call getarg(i+1, mass_file)
                    argument_retrieved=.true.

                case ("-cnx") 
                    call getarg(i+1, cnx_file)
                    argument_retrieved=.true.

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-intfile2") 
                    call getarg(i+1, intfile0)
                    argument_retrieved=.true.

                case ("-rmzfile") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.
                ! Kept for backward compatibility (but replaced by -rmzfile)
                case ("-rmz") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.

                case ("-intmode0")
                    call getarg(i+1, def_internal0)
                    argument_retrieved=.true.

                case ("-intmode")
                    call getarg(i+1, def_internal)
                    argument_retrieved=.true.
                ! Kept for backward compatibility (but replaced by -intmode)
                case ("-intset")
                    call getarg(i+1, def_internal)
                    argument_retrieved=.true.

                case ("-sym")
                    use_symmetry=.true.
                case ("-nosym")
                    use_symmetry=.false.

                case ("-symaddapt")
                    symaddapt=.true.

                case("-samerot")
                    same_red2nonred_rotation=.true.
                case("-nosamerot")
                    same_red2nonred_rotation=.false.

                ! Options to tune the model
                !================================================================
                ! This is the new (and now standard way to get the model)
                case ("-model")
                    call getarg(i+1, model)
                    argument_retrieved=.true.
                !The others are kept for backward compatibility
!                 case ("-vertQ1")
!                     vertical=.true.
!                     verticalQspace2=.false.
!                     verticalQspace1=.true.
!                 case ("-vertQ2")
!                     vertical=.true.
!                     verticalQspace2=.true.
!                     verticalQspace1=.false.
!                 case ("-vert")
!                     vertical=.true.
!                     verticalQspace2=.false.
!                 case ("-novert")
!                     vertical=.false.
!                     verticalQspace2=.false.
                !================================================================

                case ("-ref") 
                    call getarg(i+1, reference_frame)
                    argument_retrieved=.true.

                case ("-force-real")
                    force_real=.true.
                case ("-noforce-real")
                    force_real=.false.

                case ("-prjS")
                    apply_projection_matrix=.true.
                    complementay_projection=.false.
                case ("-noprjS")
                    apply_projection_matrix=.false.
                case ("-prjS-c")
                    complementay_projection=.true.
                    apply_projection_matrix=.true.
                case ("-noprjS-c")
                    complementay_projection=.false.

                case ("-orth")
                    orthogonalize=.true.
                case ("-noorth")
                    orthogonalize=.false.
                    
                case ("-origint")
                    original_internal=.true.
                case ("-noorigint")
                    original_internal=.false.
        
                case ("-h")
                    need_help=.true.

                case ("-corrS2")
                    gradcorrectS2=.true.
                    gradcorrectS2_default=.false.
                case ("-nocorrS2")
                    gradcorrectS2=.false.
                    gradcorrectS2_default=.false.
                case ("-corrS1")
                    gradcorrectS1=.true.
                    gradcorrectS1_default=.false.
                case ("-nocorrS1")
                    gradcorrectS1=.false.
                    gradcorrectS1_default=.false.

                !HIDDEN

                case ("-anaBder")
                    analytic_Bder=.true.
                case ("-noanaBder")
                    analytic_Bder=.false.

                ! Control verbosity
                case ("-quiet")
                    verbose=0
                case ("-concise")
                    verbose=1
                case ("-v")
                    verbose=2
                case ("-vv")
                    verbose=3

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       ! Manage defaults
       ! If not declared, hessfile and gradfile are the same as inpfile
       if (adjustl(hessfile) == "same") then
           hessfile=inpfile
           if (adjustl(fth) == "guess")  fth=ft
       endif
       if (adjustl(gradfile) == "same") then
           gradfile=inpfile
           if (adjustl(ftg) == "guess")  ftg=ft
       endif
       if (adjustl(hessfile2) == "same") then
           hessfile2=inpfile2
           if (adjustl(fth2) == "guess")  fth2=ft2
       endif
       if (adjustl(gradfile2) == "same") then
           gradfile2=inpfile2
           if (adjustl(ftg2) == "guess")  ftg2=ft2
       endif

       ! Set old options for the model with the new input key
       MODEL_UPPER = adjustl(model)
       call set_word_upper_case(MODEL_UPPER)
       select case (adjustl(MODEL_UPPER))
           case ("ADIA") 
               vertical=.false.
               verticalQspace1=.false.
               verticalQspace2=.false.
           case ("AS") 
               vertical=.false.
               verticalQspace1=.false.
               verticalQspace2=.false.
           case ("VERT") 
               vertical=.true.
               verticalQspace1=.false.
               verticalQspace2=.false.
           case ("VERT-A") 
               vertical=.true.
               verticalQspace1=.false.
               verticalQspace2=.false.
               move_to_min=.true.
           case ("VERTQ1") 
               vertical=.true.
               verticalQspace1=.true.
               verticalQspace2=.false.
           case ("VERTQ2") 
               vertical=.true.
               verticalQspace1=.false.
               verticalQspace2=.true.
           case default
               call alert_msg("warning","Unkown model to describe PESs: "//adjustl(model))
               need_help=.true.
       end select

       if (gradcorrectS2_default) then
           ! All vertical models include the correction
           if (vertical) then
               gradcorrectS2=.true.
           else
               gradcorrectS2=.false.
           endif
       endif
       if (gradcorrectS1_default) then
           gradcorrectS1=.false.
       endif

       ! Take defaults for the internal set for correction only
       if (def_internal0 == "defa") def_internal0=def_internal
       if (adjustl(intfile0)=="default") intfile0=intfile

       !Print options (to stdout)
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '        I N T E R N A L   D U S C H I N S K Y '    
        write(6,'(/,A)') '         Duschinski analysis for Adiabatic and    '
        write(6,'(A,/)') '        Vertical model in internal coordinates          '   
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag         Description                   Value'
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,*) '-f           Input file (State1)           ', trim(adjustl(inpfile))
        write(6,*) '-ft          \_ FileType                   ', trim(adjustl(ft))
        write(6,*) '-fhess       Hessian(S1) file              ', trim(adjustl(hessfile))
        write(6,*) '-fth         \_ FileType                   ', trim(adjustl(fth))
        write(6,*) '-fgrad       Gradient(S1) file             ', trim(adjustl(gradfile))
        write(6,*) '-ftg         \_ FileType                   ', trim(adjustl(ftg))
        write(6,*) '-f2          Input file (State2)           ', trim(adjustl(inpfile2))
        write(6,*) '-ft2         \_ FileType                   ', trim(adjustl(ft2))
        write(6,*) '-fhess2      Hessian(S2) file              ', trim(adjustl(hessfile2))
        write(6,*) '-fth2        \_ FileType                   ', trim(adjustl(fth2))
        write(6,*) '-fgrad2      Gradient(S2) file             ', trim(adjustl(gradfile2))
        write(6,*) '-ftg2        \_ FileType                   ', trim(adjustl(ftg2))
        write(6,*) '-fmass       Mass file                     ', trim(adjustl(mass_file))
        write(6,*) ''
        write(6,*) ' ** Options for state_files ** '
        write(6,*) '-ref         Reference state to output the ', reference_frame
        write(6,*) '             L matrices in its Cartesian '
        write(6,*) '             frame [I|F]'
        write(6,*) '-[no]force-real Turn imaginary frequences ', force_real
        write(6,*) '              to real (also affects Er)'
        write(6,*) '               '                       
        write(6,*) ' ** Options Internal Coordinates **           '
        write(6,*) '-[no]prjS    Apply projection matrix to    ', apply_projection_matrix
        write(6,*) '             rotate Grad and Hess'
        write(6,*) '-[no]prjS-c  Apply complementay projection ', complementay_projection
        write(6,*) '-cnx         Connectivity [filename|guess] ', trim(adjustl(cnx_file))
        write(6,*) '-intmode0    Internal set:[zmat|sel|all]   ', trim(adjustl(def_internal0))
        write(6,*) '-intfile2    File with ICs (for "sel")     ', trim(adjustl(intfile0))
        write(6,*) '-intmode     Internal set:[zmat|sel|all]   ', trim(adjustl(def_internal))
        write(6,*) '-intfile     File with ICs (for "sel")     ', trim(adjustl(intfile))
        write(6,*) '-rmzfile     File deleting ICs from Zmat   ', trim(adjustl(rmzfile))
        write(6,*) '-[no]sym     Use symmetry to form Zmat    ',  use_symmetry
        write(6,*) '-[no]samerot Use S1 red->non-red rotation ',  same_red2nonred_rotation
        write(6,*) '             for S2'
        write(6,*) '-[no]orth    Use orthogonalized internals ',  orthogonalize
        write(6,*) '-[no]origint Use originally defined inter-',  original_internal
        write(6,*) '             nal w\o linear combinations  '
        write(6,*) '             (valid for non-redundant sets)'
        write(6,*) '               '
        write(6,*) ' ** Options Vertical Model **'
        write(6,*) '-model       Model for harmonic PESs       ', trim(adjustl(model))
        write(6,*) '             [vert|vertQ1|vertQ2|vert-a|adia]     '
        write(6,*) '-[no]corrS1  Correct S1 for non-stationary ', gradcorrectS1
        write(6,*) '-[no]corrS2  Correct S2 for non-stationary ', gradcorrectS2
        write(6,*) '               '
        write(6,*) '-h           Display this help            ',  need_help
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(A)') 'Input command:'
        write(6,'(A)') trim(adjustl(input_command))   
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                       'Verbose level:  ', verbose        
        write(6,'(A)') '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


end program internal_duschinski

