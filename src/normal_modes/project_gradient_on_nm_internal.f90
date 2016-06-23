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
    use ff_build
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
               analytic_Bder=.false., &
               check_symmetry=.true., &
               orthogonalize=.false., &
               original_internal=.false., &
               force_real=.false., &
               apply_projection_matrix=.false.
    character(len=4) :: def_internal='all'
    character(len=4) :: def_internal0='all'
    character(len=1) :: reference_frame='F'
    character(len=2) :: grad_state="S2"
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1, state2
    type(str_bonded) :: geomS, geom0
    integer,dimension(NDIM) :: isym
    integer :: Nat, Nvib, Ns, NNvib, Nvib0, Ns0, NsS, Nf
    character(len=5) :: PG
    !Bonded info
    integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !B and G matrices
    real(8),dimension(NDIM,NDIM) :: B1, B2
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad, Grad1
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X, X1inv,X2inv, L1,L2,L1inv, &
                                        Asel1,Asel2,Asel1inv,Asel2inv, gBder, &
                                        G0, B0, P
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
    real(8) :: Theta, Theta2, Theta3, Theta4
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
                         intfile0 ="none", &
                         rmzfile  ="none", &
                         symm_file="none", &
                         cnx_file="guess"
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

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,&
                     inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                     intfile,intfile0,rmzfile,def_internal,def_internal0,&
                     use_symmetry,cnx_file, &
!                    tswitch,
                     symaddapt,same_red2nonred_rotation,analytic_Bder,&
                     vertical,verticalQspace2,verticalQspace1,&
                     gradcorrectS1,gradcorrectS2,&
                     orthogonalize,original_internal,force_real,reference_frame,&
                     apply_projection_matrix,grad_state)
    call set_word_upper_case(def_internal)
    call set_word_upper_case(reference_frame)
    call set_word_upper_case(grad_state)

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

    ! GRADIENT FILE
    if (grad_state /= "NO") then
        open(I_INP,file=gradfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
        call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
        close(I_INP)
        ! Use Grad1 if corrections are needed
        Grad1(1:Nvib) = Grad(1:Nvib)
    endif

    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
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
        msg = "Using custom symmetry file: "//trim(adjustl(symm_file)) 
        call alert_msg("note",msg)
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

    !Specific actions (for Nvib and Nvib0 sets)
    !*****************
    if (gradcorrectS1) then
        ! Nvib0 SET
        !-----------
        print'(/,X,A)', "COMPUTING CORRECTION FOR NON-STATIONARY   "
        print*, "------------------------------------------"
        ! The internal set for the correction does not need to be the same 
        ! as the one to represent the normal modes
        !Generate bonded info
        call gen_bonded(state1)

        !---------------------------------------
        ! NOW, GET THE ACTUAL WORKING INTERNAL SET
        call define_internal_set(state1,def_internal0,intfile0,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Aux2)
        !---------------------------------------
        ! Save the geom for the state2
        geom0=state1%geom
        Ns0=Ns

        ! Get G, B, and Bder 
        call internal_Wilson(state1,Ns,S1,B0,ModeDef)
        call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,B0,G0)
        call calc_Bder(state1,Ns,Bder,analytic_Bder)

        ! The diagonalization of the G matrix can be donne with all sets
        ! (either redundant or non-redundant), and it is the best way to 
        ! set the number of vibrations. The following call also set Nvib0
        call redundant2nonredundant(Ns,Nvib0,G0,Asel1)
        ! Rotate Bmatrix
        B0(1:Nvib0,1:3*Nat) = matrix_product(Nvib0,3*Nat,Ns,Asel1,B0,tA=.true.)
        ! Rotate Gmatrix
        G0(1:Nvib0,1:Nvib0) = matrix_basisrot(Nvib0,Ns,Asel1(1:Ns,1:Nvib0),G0,counter=.true.)
        ! Rotate Bders
        do j=1,3*Nat
            Bder(1:Nvib0,j,1:3*Nat) =  matrix_product(Nvib0,3*Nat,Ns,Asel1,Bder(1:Ns,j,1:3*Nat),tA=.true.)
        enddo

        if (apply_projection_matrix) then
            ! Get projection matrix
            P(1:3*Nat,1:3*Nat) = projection_matrix(Nat,Nvib0,B0)
            ! And rotate gradient
            do i=1,3*Nat
                Vec1(i) = 0.d0
                do k=1,Nvib0
                    Vec1(i) = Vec1(i) + P(i,k)*Grad1(k)
                enddo
            enddo
            Grad1(1:3*Nat) = Vec1(1:3*Nat)
        endif

        ! Get the Correction now
        print*, " Getting the correction term: gs^t\beta"
        ! The correction is applied with the Nvib0 SET
        ! Correct Hessian as
        ! Hx' = Hx - gs^t\beta
        ! 1. Get gs from gx
        Vec1(1:3*Nat) = Grad1(1:3*Nat)
        call Gradcart2int(Nat,Nvib0,Vec1,state1%atom(:)%mass,B0,G0)
        ! 2. Multiply gs^t\beta and
        ! 3. Apply the correction
        ! Bder(i,j,K)^t * gq(K)
        do i=1,3*Nat
        do j=1,3*Nat
            gBder(i,j) = 0.d0
            do k=1,Nvib0
                gBder(i,j) = gBder(i,j) + Bder(k,i,j)*Vec1(k)
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
    endif


    !Nvib SET
    !---------
    print'(/,X,A)', "GETTING INTERNAL SET TO DESCRIBE MODES   "
    print*, "------------------------------------------"
    ! Refress connectivity
    call gen_bonded(state1)

    ! Define internal set
    call define_internal_set(state1,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns,Nf,Aux2)
    ! Save the geom for the state2
    geomS=state1%geom
    NsS=Ns

    !From now on, we'll use atomic units
    call set_geom_units(state1,"Bohr")


    ! INTERNAL COORDINATES

    !SOLVE GF METHOD TO GET NM AND FREQ
    call internal_Wilson(state1,Ns,S1,B1,ModeDef)
    call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,B1,G1)

    ! SET REDUNDANT/SYMETRIZED/CUSTOM INTERNAL SETS
!     if (symaddapt) then (implement in an analogous way as compared with the transformation from red to non-red
    if (original_internal.and.Nvib==Ns) then
        print*, "Using internal without linear combination"
        Asel1(1:Ns,1:Nvib) = 0.d0
        do i=1,Nvib
            Asel1(i,i) = 1.d0 
        enddo
        Asel1inv(1:Nvib,1:Ns) = Asel1(1:Ns,1:Nvib)
    else
        print*, "Using internal from eigevectors of G"
        call redundant2nonredundant(Ns,Nvib,G1,Asel1)
        ! Rotate Gmatrix
        G1(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel1(1:Ns,1:Nvib),G1,&
                                            counter=.true.)
        ! Check if we want orthogonalization
        if (orthogonalize) then
            print*, "Orthogonalyzing state1 internals..."
            X1inv(1:Nvib,1:Nvib) = 0.d0
            X(1:Nvib,1:Nvib)     = 0.d0
            do i=1,Nvib
                X1inv(i,i) = 1.d0/dsqrt(G1(i,i))
                X(i,i)     = dsqrt(G1(i,i))
            enddo
            ! Rotate Gmatrix (again)
            G1(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Nvib,X1inv(1:Nvib,1:Nvib),G1)
            ! Update Asel(inv)
            Asel1inv(1:Nvib,1:Ns) = matrix_product(Nvib,Ns,Nvib,X1inv,Asel1,tB=.true.)
            Asel1(1:Ns,1:Nvib)    = matrix_product(Ns,Nvib,Nvib,Asel1,X)
        else
            Asel1inv(1:Nvib,1:Ns) = transpose(Asel1(1:Ns,1:Nvib))
        endif
        ! Rotate Bmatrix
        B1(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel1inv,B1)
    endif

    if (apply_projection_matrix) then
        ! Get projection matrix
        P(1:3*Nat,1:3*Nat) = projection_matrix(Nat,Nvib,B1)
        ! And rotate gradient
        do i=1,3*Nat
            Vec1(i) = 0.d0
            do k=1,Nvib
                Vec1(i) = Vec1(i) + P(i,k)*Grad(k)
            enddo
        enddo
        Grad(1:3*Nat) = Vec1(1:3*Nat)
    endif

    if (gradcorrectS1) then
        do i=1,3*Nat
        do j=1,3*Nat
            ! Apply correction to the Hessian term
            Hess(i,j) = Hess(i,j) - gBder(i,j)
        enddo
        enddo
    endif

    if (apply_projection_matrix) then
        ! Project out rotation and translation
        Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,P,Hess)
    endif

    call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B1,G1)
    call gf_method(Nvib,Nvib,G1,Hess,L1,Freq1,X,X1inv)
    if (verbose>0) then
        ! Analyze normal modes in terms of the redundant set
        Aux(1:Ns,1:Nvib) = matrix_product(Ns,Nvib,Nvib,Asel1,L1)
        if (use_symmetry) then
            call analyze_internal(Nvib,Ns,Aux,Freq1,ModeDef,S_sym)
        else
            call analyze_internal(Nvib,Ns,Aux,Freq1,ModeDef)
        endif
    endif

    ! Convert also the gradient to internal (for future use)
    call Gradcart2int(Nat,Nvib,Grad,state1%atom(:)%mass,B1,G1)


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

    ! HESSIAN FILE: Not read

    ! GRADIENT FILE
    if (grad_state=="S2") then
        print*, "Reading gradient from State2"
        open(I_INP,file=gradfile2,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile2)) )
        call generic_gradient_reader(I_INP,ftg2,Nat,Grad,error)
        if (error /= 0) call alert_msg("fatal","Error reading gradient (State2)")
        close(I_INP)
    endif


    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
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
        msg = "Using custom symmetry file: "//trim(adjustl(symm_file)) 
        call alert_msg("note",msg)
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
    if (state1%PG /= state2%PG) then
        print*, "PG(State1): ", state1%PG
        print*, "PG(State2): ", state2%PG
        if (.not.use_symmetry) then
            call alert_msg("note","Initial and final state have different symmetry")
        else
            call alert_msg("warning","Initial and final state have different symmetry")
        endif
    endif

    !Nvib SET
    !---------
    print'(/,X,A)', "GETTING INTERNAL SET TO DESCRIBE MODES   "
    print*, "------------------------------------------"
    ! Refress connectivity
    call gen_bonded(state2)

    ! Define internal set => taken from state1
    state2%geom = geomS
    Ns = NsS

    !From now on, we'll use atomic units
    call set_geom_units(state2,"Bohr")


    ! INTERNAL COORDINATES

    !SOLVE GF METHOD TO GET NM AND FREQ
    call internal_Wilson(state2,Ns,S2,B2,ModeDef)
    call internal_Gmetric(Nat,Ns,state2%atom(:)%mass,B2,G2)
    ! Handle redundant/symtrized sets
!     if (symaddapt) then (implement in an analogous way as compared with the transformation from red to non-red
    if (original_internal.and.Nvib==Ns) then
        print*, "Using internal without linear combination"
        Asel2(1:Ns,1:Nvib) = 0.d0
        do i=1,Nvib
            Asel2(i,i) = 1.d0 
        enddo
        Asel2inv(1:Nvib,1:Ns) = Asel2(1:Ns,1:Nvib)
    else
        if (same_red2nonred_rotation) then
            print*, "Using internal from eigevectors of G (state1)"
            ! Using Asel1 (from state1)
            Asel2(1:Ns,1:Nvib) = Asel1(1:Ns,1:Nvib)
            Asel2inv(1:Nvib,1:Ns) = Asel1inv(1:Nvib,1:Ns)
            ! Rotate Gmatrix
            G2(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel2inv(1:Nvib,1:Ns),G2)
        else
            print*, "Using internal from eigevectors of G (state2)"
            call redundant2nonredundant(Ns,Nvib,G2,Asel2)
            ! Rotate Gmatrix
            G2(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel2(1:Ns,1:Nvib),G2,&
                                                counter=.true.)
            if (orthogonalize) then
                print*, "Orthogonalyzing state2 internals..."
                X2inv(1:Nvib,1:Nvib) = 0.d0
                X(1:Nvib,1:Nvib)     = 0.d0
                do i=1,Nvib
                    X(i,i)     = dsqrt(G2(i,i))
                    X2inv(i,i) = 1.d0/X(i,i)
                enddo
                ! Rotate Gmatrix (again)
                G2(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Nvib,X2inv(1:Nvib,1:Nvib),G2)
                ! Update Asel(inv)
                Asel2inv(1:Nvib,1:Ns) = matrix_product(Nvib,Ns,Nvib,X2inv,Asel2,tB=.true.)
                Asel2(1:Ns,1:Nvib)    = matrix_product(Ns,Nvib,Nvib,Asel2,X)
            else
                Asel2inv(1:Nvib,1:Ns) = transpose(Asel2(1:Ns,1:Nvib))
            endif
        endif
        ! Rotate Bmatrix
        B2(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel2inv,B2)
    endif

    if (grad_state=="S2") then
        if (apply_projection_matrix) then
            ! Get projection matrix
            P(1:3*Nat,1:3*Nat) = projection_matrix(Nat,Nvib,B2)
            ! And rotate gradient
            do i=1,3*Nat
                Vec1(i) = 0.d0
                do k=1,Nvib
                    Vec1(i) = Vec1(i) + P(i,k)*Grad(k)
                enddo
            enddo
            Grad(1:3*Nat) = Vec1(1:3*Nat)
        endif
        ! Convert also the gradient to internal (for future use)
        call Gradcart2int(Nat,Nvib,Grad,state2%atom(:)%mass,B2,G2)
    endif

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
    ! ICs displacement
    !--------------------------
    if (verbose>0) then
        print*, ""
        print*, "=========================="
        print*, " SHIFTS (internal coord)"
        print*, "=========================="
    endif
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
    ! Need to transform Deltas into the non-redundant coordinates set
    ! Delta' = A^-1 Delta
    do i=1,Nvib
        Vec2(i) = 0.d0
        do k=1,Ns
            Vec2(i) = Vec2(i) + Asel1inv(i,k)*Delta(k)
        enddo
    enddo
    Delta(1:Nvib) = Vec2(1:Nvib)


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


    !Inverse of L1 (and store in L1inv)
    L1inv(1:Nvib,1:Nvib) = inverse_realgen(Nvib,L1(1:Nvib,1:Nvib))
    !--------------------------
    ! K-vector (NM-shifts)
    !--------------------------
    ! K = L1^-1 DeltaS (this is State 1 respect to state 2) . 
    ! Notes
    !   * L1inv stores the inverse
    !   * Delta in non-redundant IC set
    do i=1,Nvib
        Vec1(i) = 0.d0
        do k=1,Nvib
            Vec1(i) = Vec1(i) + L1inv(i,k)*Delta(k)
        enddo
    enddo


    !==================================
    ! At this point
    !---------------
    ! Vec1: K vector
    ! Hess: Hessian  (S2) in internal coordinates
    ! Grad: Gradient (S2) in internal coordinates

    ! Get coordinate in the direction of the gradient in the nm space
    !  g_Q = L1^t gs
    do i=1,Nvib
        Vec2(i) = 0.d0
        do k=1,Nvib
            Vec2(i) = Vec2(i) + L1(k,i) * Grad(k)
        enddo
    enddo
    call print_vector(6,Vec1,Nvib,"Normal mode shift")

    call print_vector(6,Grad,Nvib,"Gradient in NM")

    ! Compute harmonic gradient
    FC(1:Nvib) = Freq2FC(Nvib,Freq1)
    do i=1,Nvib
        Grad(i) = Vec1(i)*FC(i)
    enddo

    call print_vector(6,Grad,Nvib,"Harmonic gradient in NM")

    print*, ""

    Theta2 = vector_dot_product(Nvib,Vec1,Vec1)
    Vec1(1:Nvib) = Vec1(1:Nvib)/dsqrt(Theta2)
    Theta = vector_dot_product(Nvib,Vec2,Vec2)
    Vec2(1:Nvib) = Vec2(1:Nvib)/dsqrt(Theta)
    Theta4 = vector_dot_product(Nvib,Grad,Grad)
    Grad(1:Nvib) = Grad(1:Nvib)/dsqrt(Theta4)

    call print_vector(6,Vec1,Nvib,"Normalized straight direction")
    call print_vector(6,Vec2,Nvib,"Normalized gradient direction")

    Theta3 = vector_dot_product(Nvib,Vec1,Vec2)
!     print*, Theta, dacos(Theta)*180.d0/pi
    print'(X,A18,X,F12.4)', "Distance =", dsqrt(Theta2/AMUtoAU)
    print'(X,A18,X,ES12.4)', "|Grad| =", dsqrt(Theta)
    print'(X,A18,X,ES12.4)', "|Grad(H)| =", dsqrt(Theta4)
    print'(X,A18,X,F12.4)', "Angle w/ grad =", dacos(Theta3)*180.d0/pi
    Theta3 = vector_dot_product(Nvib,Grad,Vec1)
    print'(X,A18,X,F12.4)', "Angle w/ grad(H) =", dacos(Theta3)*180.d0/pi
    Theta3 = vector_dot_product(Nvib,Grad,Vec2)
    print'(X,A18,X,F12.4)', "Angle btw grads =", dacos(Theta3)*180.d0/pi
    print*, ""

    call summary_alerts

    call cpu_time(tf)
    write(6,'(A,F12.3)') "CPU (s) for internal vib analysis: ", tf-ti

    stop

    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,&
                           inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                           intfile,intfile0,rmzfile,def_internal,def_internal0,&
                           use_symmetry,cnx_file,&
!                          tswitch,
                           symaddapt,same_red2nonred_rotation,analytic_Bder,&
                           vertical,verticalQspace2,verticalQspace1,&
                           gradcorrectS1,gradcorrectS2,&
                           orthogonalize,original_internal,force_real,reference_frame,&
                           apply_projection_matrix,grad_state)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,&
                                          inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                                          intfile,intfile0,rmzfile,def_internal,def_internal0,&
                                          cnx_file,reference_frame, & !, symfile
                                          grad_state
        logical,intent(inout)          :: use_symmetry, vertical, verticalQspace2, &
                                          verticalQspace1, &
                                          gradcorrectS1, gradcorrectS2, symaddapt, &
                                          same_red2nonred_rotation,analytic_Bder, &
                                          orthogonalize,original_internal,force_real, &
                                          apply_projection_matrix
!         logical,intent(inout) :: tswitch

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        character(len=500) :: input_command
        character(len=10)  :: model="adia", MODEL_UPPER
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

                case ("-cnx") 
                    call getarg(i+1, cnx_file)
                    argument_retrieved=.true.

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-intfile0") 
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
                case ("-vertQ1")
                    vertical=.true.
                    verticalQspace2=.false.
                    verticalQspace1=.true.
                case ("-vertQ2")
                    vertical=.true.
                    verticalQspace2=.true.
                    verticalQspace1=.false.
                case ("-vert")
                    vertical=.true.
                    verticalQspace2=.false.
                case ("-novert")
                    vertical=.false.
                    verticalQspace2=.false.
                !================================================================

                case ("-ref") 
                    call getarg(i+1, reference_frame)
                    argument_retrieved=.true.

                case ("-force-real")
                    force_real=.true.
                case ("-noforce-real")
                    force_real=.false.

                case ("-prj-tr")
                    force_real=.true.
                case ("-noprj-tr")
                    force_real=.false.

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

                case ("-gradS")
                    call getarg(i+1, grad_state)
                    argument_retrieved=.true.
              

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
           case ("VERT") 
               vertical=.true.
               verticalQspace1=.false.
               verticalQspace2=.false.
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
        write(6,*) ''
        write(6,*) '-gradS       State from which Grad is      ', grad_state
        write(6,*) '             taken [S1|S2]'
!         write(6,*) ' ** Options for state_files ** '
!         write(6,*) '-ref         Reference state to output the ', reference_frame
!         write(6,*) '             L matrices in its Cartesian '
!         write(6,*) '             frame [I|F]'
!         write(6,*) '-[no]force-real Turn imaginary frequences ', force_real
!         write(6,*) '              to real (also affects Er)'
        write(6,*) '               '                       
        write(6,*) ' ** Options Internal Coordinates **           '
        write(6,*) '-[no]prj-tr  Apply projection matrix to   ', apply_projection_matrix
        write(6,*) '             rotate Grad and Hess'
        write(6,*) '-cnx         Connectivity [filename|guess] ', trim(adjustl(cnx_file))
        write(6,*) '-intmode0    Internal set:[zmat|sel|all]   ', trim(adjustl(def_internal0))
        write(6,*) '-intfile0    File with ICs (for "sel")     ', trim(adjustl(intfile0))
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
!         write(6,*) '-model       Model for harmonic PESs       ', trim(adjustl(model))
!         write(6,*) '             [vert|vertQ1|vertQ2|adia]     '
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

