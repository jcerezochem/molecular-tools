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

    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: use_symmetry=.false. ,&
               modred=.false.       ,&
               tswitch=.false.      ,&
               symaddapt=.false.    ,&
               vertical=.false.     ,&
               vertical_method2=.false. ,&
               only_state2=.true.   ,&  
               gradcorrect=.true.
    character(len=4) :: def_internal='zmat'
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1, state2
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns
    character(len=5) :: PG
    !Bonded info
    integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !B and G matrices
    real(8),dimension(NDIM,NDIM) :: B
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X, X1inv,X2inv, L1,L2, Asel1,Asel2
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    !Duschisky
    real(8),dimension(NDIM,NDIM) :: G1, G2
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
               i1,i2,i3,i4
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
                         intfile  ="none",       &
                         rmzfile  ="none",       &
                         symm_file="none"
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

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,&
                     inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                     intfile,rmzfile,def_internal,use_symmetry, &
!                    tswitch,symaddapt, &
                     vertical,vertical_method2,only_state2,gradcorrect)
    call set_word_upper_case(def_internal)

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
    open(I_INP,file=gradfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
    call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
    close(I_INP)

    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! Get connectivity 
    call guess_connect(state1)

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

    !Generate bonded info
    call gen_bonded(state1)

    ! Define internal set
    call define_internal_set(state1,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns)

    !From now on, we'll use atomic units
    call set_geom_units(state1,"Bohr")


    ! INTERNAL COORDINATES

    !SOLVE GF METHOD TO GET NM AND FREQ
    call internal_Wilson(state1,Nvib,S1,B,ModeDef)
    call internal_Gmetric(Nat,Nvib,state1%atom(:)%mass,B,G1)
    if (vertical.and..not.only_state2.and.gradcorrect) then
        call NumBder(state1,Nvib,Bder)
        call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B,G1,Grad=Grad,Bder=Bder)
    else
        call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B,G1)
    endif
    call gf_method(Nvib,G1,Hess,L1,Freq1,X,X1inv)
    if (verbose>0) then
        ! Analyze normal modes
        if (use_symmetry) then
            call analyze_internal(Nvib,L1,Freq1,ModeDef,S_sym)
        else
            call analyze_internal(Nvib,L1,Freq1,ModeDef)
        endif
    endif

    ! Compute new state_file
    ! T1(g09) = mu^1/2 m B^t G1^-1 L1
    call Ls_to_Lcart(Nat,Nvib,state1%atom(:)%mass,B,G1,L1,Aux,error)
    call Lcart_to_LcartNrm(Nat,Nvib,Aux,Aux2,error)
    !Checking normalization
    if (verbose>0) then
        print'(/,X,A)', "Checking normalization of Tcart (G09)"
        do j=1,Nvib
            theta=0.d0
            do i=1,3*Nat
                theta=theta+Aux2(i,j)**2
            enddo
            print'(F15.8)', theta
        enddo
        print*, ""
    endif
    !Print state
    open(O_STAT,file="state_file_1")
    call set_geom_units(state1,"Angs")
    do i=1,Nat
        write(O_STAT,*) state1%atom(i)%x
        write(O_STAT,*) state1%atom(i)%y
        write(O_STAT,*) state1%atom(i)%z
    enddo
    do i=1,3*Nat
    do j=1,Nvib
        write(O_STAT,*) Aux2(i,j)
    enddo
    enddo
    do j=1,Nvib
        write(O_STAT,'(F12.5)') Freq1(j)
    enddo
    close(O_STAT)

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
    ! Shortcuts
    if (Nat /= state2%natoms) call alert_msg("fatal","Initial and final states don't have the same number of atoms.")

    ! HESSIAN FILE
    open(I_INP,file=hessfile2,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile2)) )
    allocate(A(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fth2,Nat,A,error) 
    if (error /= 0) call alert_msg("fatal","Error reading Hessian (State2)")
    close(I_INP)
    ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
    call vibrations_Cart(Nat,state2%atom(:)%X,state2%atom(:)%Y,state2%atom(:)%Z,state2%atom(:)%Mass,A,&
                         Nvib,L2,Freq2,error_flag=error)
    if (Nvib /= 3*Nat-6) call alert_msg("warning","Linear molecule (at state2). Things can go very bad.")
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
    open(I_INP,file=gradfile2,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile2)) )
    call generic_gradient_reader(I_INP,ftg2,Nat,Grad,error)
    if (error /= 0) call alert_msg("fatal","Error reading gradient (State2)")
    close(I_INP)

    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! Get connectivity 
    call guess_connect(state2)

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
        if (.not.use_symmetry) then
            call alert_msg("note","Initial and final state have different symmetry")
        else
            call alert_msg("warning","Initial and final state have different symmetry")
        endif
    endif

    !Generate bonded info
    call gen_bonded(state2)

    ! Define internal set => taken from state1
    state2%geom = state1%geom

    !From now on, we'll use atomic units
    call set_geom_units(state2,"Bohr")


    ! INTERNAL COORDINATES

    !SOLVE GF METHOD TO GET NM AND FREQ
    call internal_Wilson(state2,Nvib,S2,B,ModeDef)
    call internal_Gmetric(Nat,Nvib,state2%atom(:)%mass,B,G2)
    if (vertical.and.gradcorrect) then
        call NumBder(state2,Nvib,Bder)
        call HessianCart2int(Nat,Nvib,Hess,state2%atom(:)%mass,B,G2,Grad=Grad,Bder=Bder)
    else
        call HessianCart2int(Nat,Nvib,Hess,state2%atom(:)%mass,B,G2)
    endif
    call gf_method(Nvib,G2,Hess,L2,Freq2,X,X2inv)
    if (verbose>0) then
        ! Analyze normal modes
        if (use_symmetry) then
            call analyze_internal(Nvib,L2,Freq2,ModeDef,S_sym)
        else
            call analyze_internal(Nvib,L2,Freq2,ModeDef)
        endif
    endif

    ! Compute new state_file
    call Ls_to_Lcart(Nat,Nvib,state2%atom(:)%mass,B,G2,L2,Aux,error)
    call Lcart_to_LcartNrm(Nat,Nvib,Aux,Aux2,error)
    !Checking normalization
    if (verbose>0) then
        print'(/,X,A)', "Checking normalization of Tcart (G09)"
        do j=1,Nvib
            theta=0.d0
            do i=1,3*Nat
                theta=theta+Aux2(i,j)**2
            enddo
            print'(F15.8)', theta
        enddo
        print*, ""
    endif
    !Print state
    open(O_STAT,file="state_file_2")
    call set_geom_units(state2,"Angs")
    do i=1,Nat
        write(O_STAT,*) state2%atom(i)%x
        write(O_STAT,*) state2%atom(i)%y
        write(O_STAT,*) state2%atom(i)%z
    enddo
    do i=1,3*Nat
    do j=1,Nvib
        write(O_STAT,*) Aux2(i,j)
    enddo
    enddo
    do j=1,Nvib
        write(O_STAT,'(F12.5)') Freq2(j)
    enddo
    close(O_STAT)


    !==========================================
    ! CHECKS ON THE INTERNAL SETS
    !==========================================
    ! Evaluate orthogonality
    if (verbose>0) then
     print'(/,X,A)', "Checking simultaneous orthogonality: D = G1^1/2 G2^1/2"
     print*,         "-------------------------------------------------------"
     print*,         "Analysing: D = G1^1/2 G2^1/2"
    endif
    Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,X1inv,X)
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
     print'(X,A,/)', "(D matrix has been written to files: D_matrix.dat and D_matrix_abs.dat)"
    
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
    ! Orthogonal Duschinski
    !--------------------------
    ! Get orthogonal modes:  L' = G^-1/2 L
    Aux(1:Nvib,1:Nvib)  = matrix_product(Nvib,Nvib,Nvib,X1inv,L1)
    Aux2(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,X2inv,L2)
    ! Duschinsky matrix (orth) stored in G2 = L1'^t L2'
    G2(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,Aux,Aux2,tA=.true.)
    !Store L1' in Aux2 to later be used to get the displacement
    Aux2(1:Nvib,1:Nvib)=Aux(1:Nvib,1:Nvib)


   !-- IF REDUNDANT
    if (Ns > Nvib) then
        if (verbose>0) &
         print'(/,X,A,/)', "Working with redundant coordianates. Computing A matrices"
        !----------------------------------------------------
        ! Stop here for now
        print*, " *** UNDER DEVELOPMENT *** "
        call alert_msg("fatal","Internal analysis with redundant coordianates is not yet ready")
        !----------------------------------------------------
        !Prepare A matrix = A2^T * A1
        Aux2(1:Ns,1:Ns) = matrix_product(Ns,Ns,Ns,Asel2,Asel1,tA=.true.)
        !Compute inverse as A1^T * A2, eventually store in Asel1
        Asel1(1:Ns,1:Ns) = matrix_product(Ns,Ns,Ns,Asel1,Asel2,tA=.true.)
        if (verbose>0) then
            Aux(1:Ns,1:Ns) = matrix_product(Ns,Ns,Ns,Asel1,Aux2)
            call MAT0(6,Aux,Ns,Ns,"Asel test")
        endif
    endif
   !-- ENDOF IF REDUNDANT

    !--------------------------
    ! Non-Orthogonal Duschinski
    !--------------------------
    if (verbose>0) &
     print*, "Calculating Duschisky..."
    !Inverse of L1 (and store in L1)
    L1(1:Nvib,1:Nvib) = inverse_realgen(Nvib,L1(1:Nvib,1:Nvib))
    !J = L1^-1 L2 (stored in G1).
    G1(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,L1,L2)

    if (vertical.and..not.vertical_method2) then
        ! GET MINIMUM IN INTERNAL COORDINATES
        ! At this point
        !  * Hess has the Hessian  of State2 in internal coords (output from HessianCart2int)
        !  * Grad has the gradient of State2 in internal coords (output from HessianCart2int)
        ! Inverse of the Hessian
        Aux(1:Nvib,1:Nvib) = inverse_realgen(Nvib,Hess)
        ! DeltaS0 = -Hs^1 * gs
        do i=1,Nvib
            Delta(i)=0.d0
            do k=1,Nvib
                Delta(i) = Delta(i)-Aux(i,k) * Grad(k)
            enddo 
        enddo
        if (verbose>2) then
            k=0
            print*, "DELTA BONDS"
            do i=1,state1%geom%nbonds
                k = k+1
                print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Delta(k), Delta(k)*BOHRtoANGS
            enddo
            print*, "DELTA ANGLES"
            do i=1,state1%geom%nangles
                k = k+1
                print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Delta(k), Delta(k)*180.d0/PI
            enddo
            print*, "DELTA DIHEDRALS"
            do i=1,state1%geom%ndihed
                k = k+1
                print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Delta(k), Delta(k)*180.d0/PI
            enddo
        endif
        ! Get coordinates
        do i=1,Nvib
            S2(i) = S1(i) + Delta(i)
        enddo
    endif


    if (.not.vertical_method2) then
        if (verbose>0) then
            print*, ""
            print*, "=========================="
            print*, " SHIFTS (internal coord)"
            print*, "=========================="
        endif
        if (verbose>0) &
         print*, "Bonds"
        do i=1,state1%geom%nbonds
            Delta(i) = S2(i)-S1(i)
            if (verbose>0) &
             print'(I5,3(F8.2,2X))', i, S2(i),S1(i), Delta(i)
        enddo
        if (verbose>0) &
         print*, "Angles"
        do j=i,i+state1%geom%nangles-1
            Delta(j) = S2(j)-S1(j)
            if (verbose>0) &
             print'(I5,3(F8.2,2X))', j, S2(j)*180.d0/PI,S1(j)*180.d0/PI,Delta(j)*180.d0/PI
        enddo
        if (verbose>0) &
         print*, "Dihedrals"
        do k=j,j+state1%geom%ndihed-1
            Delta(k) = S2(k)-S1(k)
            Delta_p = S2(k)-S1(k)+2.d0*PI
            if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
            Delta_p = S2(k)-S1(k)-2.d0*PI
            if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
            if (verbose>0) &
             print'(I5,3(F8.2,2X))', k, S2(k)*180.d0/PI,S1(k)*180.d0/PI,Delta(k)*180.d0/PI
        enddo

    else ! vertical with method2 
        ! K = -J * Lambda_f^-1 * L2^t * gs
        ! Convert Freq into FC. Store in FC for future use
        do i=1,Nvib
            FC(i) = sign((Freq2(i)*2.d0*pi*clight*1.d2)**2/HARTtoJ*BOHRtoM**2*AUtoKG,Freq2(i))
            if (FC(i)<0) then
                print*, i, FC(i)
!                 FC(i) = -FC(i)
                call alert_msg("warning","A negative FC found")
            endif
        enddo
        ! Lambda_f^-1 * L2^t
        do i=1,Nvib
            Aux(i,1:Nvib) = L2(1:Nvib,i) / FC(i)
        enddo
        ! -[Lambda_f^-1 * L2^t] * gs
        do i=1,Nvib
            Q0(i)=0.d0
            do k=1,Nvib
                Q0(i) = Q0(i) - Aux(i,k) * Grad(k)
            enddo
        enddo
        ! J * [-Lambda_f^-1 * L2^t * gs]
        do i=1,Nvib
            Vec1(i)=0.d0
            do k=1,Nvib
                Vec1(i) = Vec1(i) + G1(i,k) * Q0(k)
            enddo
        enddo
    endif

!     print*, ""
!     print*, "=========================="
!     print*, " SHIFTS (orthogonal internal coord)"
!     print*, "=========================="
!     print*, "Setting orthogonal internal coordinates..."
!     do i=1,Nvib
!         Vec1(i) = 0.d0
!         do j=1,Nvib
!             Vec1(i) = Vec1(i) + X1inv(
!         enddo
!     enddo

    if (symaddapt) then
        if (vertical.or.vertical_method2) call alert_msg("fatal","Symmetry conflicts with vertical")
        print*, ""
        print*, "Using symmetry addapted coordinates"
        print*, "                     Coord                          Displacement"
        print*, "  ---------------------------------------------------------------"
        do i=1,Nvib
            if (S_sym(i) <= i) cycle
            j=S_sym(i)
            Vec1(1) = Delta(i)+Delta(j)
            Vec1(2) = Delta(i)-Delta(j)
            Delta(i) = Vec1(1)
            Delta(j) = Vec1(2)
        enddo
        do i=1,Nvib
            print'(X,A45,3X,F12.5)', trim(adjustl(ModeDef(i))), Delta(i)
        enddo
        print*, "  ---------------------------------------------------------------"

    elseif (Ns /= Nvib) then
        if (vertical.or.vertical_method2) call alert_msg("fatal","Reduced/Augmented spaces conflict with vertical")
        do i=1,Ns
            Vec1(i) = 0.d0
            do k=1,Ns
                Vec1(i) = Vec1(i) + Asel1(i,k)*Delta(k)
            enddo
        enddo
        Delta(1:Nvib) = Vec1(1:Nvib)
    endif


    if (.not.vertical_method2) then
        ! K = L1^-1 DeltaS (this is State 1 respect to state 2) . L1 already stores the inverse!
        do i=1,Nvib
            Vec1(i) = 0.d0
            do k=1,Nvib
                Vec1(i) = Vec1(i) + L1(i,k)*Delta(k)
            enddo
        enddo
        !Orthogonal: K=L1'^t DeltaS'
        do i=1,Nvib
            Vec2(i) = 0.d0
            do k=1,Nvib
                Vec2(i) = Vec1(i) + Aux2(k,i)*Delta(k)
            enddo
        enddo
    endif
    if (verbose>1) &
     call MAT0(6,G1,Nvib,Nvib,"DUSCHINSKI MATRIX")

    !Analyze Duschinsky matrix
    call analyze_duschinsky(6,Nvib,G1,Vec1,Freq1,Freq2)


    !===================================
    ! Reorganization energy
    !===================================
    if (vertical_method2) then
        ! Normal-mode space
        ! Er = -L2^t gs * Q0 - 1/2 * Q0^t * Lambda_f * Q0
        ! At this point: 
        ! * Grad: in internal coords
        ! * Q0: DeltaQ in final state nm
        ! * FC: diagonal force constants for final state
        Er = 0.d0
        do i=1,Nvib
            ! Compute gQ(i) = L2^t * gs
            Theta = 0.d0
            do j=1,Nvib
                Theta =  Theta + L2(j,i)*Grad(j)
            enddo
            Er = Er - Theta * Q0(i) - 0.5d0 * FC(i) * Q0(i)**2
        enddo
    else
        ! Internal-coordinates space
        ! Er = -gs * DeltaS - 1/2 DeltaS^t * Hs * DeltaS
        ! At this point: 
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
    endif
    print'(X,A,F12.6)',   "Reorganization energy (AU) = ", Er
    print'(X,A,F12.6,/)', "Reorganization energy (eV) = ", Er*HtoeV


    !============================================
    ! PRINT DUSCHINSKI AND DISPLACEMENT TO FILES
    !============================================
    !Prints Duschisky matrix and displacements to files
    open(O_DUS, file="duschinsky.dat")
    open(O_DUS2,file="duschinsky_orth.dat")
    open(O_DIS, file="displacement.dat")
    open(O_DIS2,file="displacement_orth.dat")
    do i=1,Nvib
    do j=1,Nvib
        write(O_DUS,*)  G1(i,j)
        write(O_DUS2,*) G2(i,j)
    enddo 
        write(O_DIS,*)  Vec1(i)
        write(O_DIS2,*) Vec2(i)
    enddo
    close(O_DUS)
    close(O_DUS2)
    close(O_DIS)
    close(O_DIS2)

    call summary_alerts

    call cpu_time(tf)
    write(6,'(A,F12.3)') "CPU (s) for internal vib analysis: ", tf-ti

    stop

    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,&
                           inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                           intfile,rmzfile,def_internal,use_symmetry, &
!                          tswitch,symaddapt, & symfile
                           vertical,vertical_method2,only_state2,gradcorrect)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,&
                                          inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                                          intfile,rmzfile,def_internal !, symfile
        logical,intent(inout)          :: use_symmetry, vertical, vertical_method2, &
                                          only_state2, gradcorrect
!         logical,intent(inout) :: tswitch, symaddapt

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg

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

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-rmzfile") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.
                ! Kept for backward compatibility (but replaced by -rmzfile)
                case ("-rmz") 
                    call getarg(i+1, rmzfile)
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

                case ("-vert2")
                    vertical=.true.
                    vertical_method2=.true.
                case ("-vert")
                    vertical=.true.
                    vertical_method2=.false.
                case ("-novert")
                    vertical=.false.
                    vertical_method2=.false.

                case ("-onlys2")
                    only_state2=.true.
                case ("-noonlys2")
                    only_state2=.false.
        
                case ("-h")
                    need_help=.true.

                !HIDDEN

                case ("-correct")
                    gradcorrect=.true.
                case ("-nocorrect")
                    gradcorrect=.false.

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

       !Print options (to stderr)
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,'(/,A)') '        I N T E R N A L   A N A L Y S I S '    
        write(6,'(/,A)') '      Perform vibrational analysis based on  '
        write(6,'(/,A)') '            internal coordinates (D-V9.1)'        
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,*) '-f              ', trim(adjustl(inpfile))
        write(6,*) '-ft             ', trim(adjustl(ft))
        write(6,*) '-fhess          ', trim(adjustl(hessfile))
        write(6,*) '-fth            ', trim(adjustl(fth))
        write(6,*) '-fgrad          ', trim(adjustl(gradfile))
        write(6,*) '-ftg            ', trim(adjustl(ftg))
        write(6,*) '-f2             ', trim(adjustl(inpfile2))
        write(6,*) '-ft2            ', trim(adjustl(ft2))
        write(6,*) '-fhess2         ', trim(adjustl(hessfile2))
        write(6,*) '-fth2           ', trim(adjustl(fth2))
        write(6,*) '-fgrad2         ', trim(adjustl(gradfile2))
        write(6,*) '-ftg2           ', trim(adjustl(ftg2))
        write(6,*) '-intmode        ', trim(adjustl(def_internal))
        write(6,*) '-intfile        ', trim(adjustl(intfile))
        write(6,*) '-rmzfile        ', trim(adjustl(rmzfile))
        write(6,*) '-[no]sym       ',  use_symmetry
        write(6,*) '-[no]vert      ',  vertical
        write(6,*) '-[no]correct   ',  gradcorrect
        write(6,*) '-vert2         ',  vertical_method2
        write(6,*) '-[no]onlys2    ',  only_state2
        write(6,*) '-h             ',  need_help
        write(6,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


end program internal_duschinski

