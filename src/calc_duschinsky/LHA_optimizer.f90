program optimizer


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
    use xyz_manage_molec
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
    logical :: use_symmetry=.false. ,&
               modred=.false.       ,&
               tswitch=.false.      ,&
               symaddapt=.false.    ,&
               vertical=.true.      ,&
               apply_projection_matrix=.false.
    character(len=4) :: def_internal='zmat', def_internal_aux
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1,state2
    type(str_bonded) :: zmatgeom, allgeom
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !B and G matrices
    real(8),dimension(NDIM,NDIM) :: B1,B2, B, G1,G2,Ginv
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad, FC, Q0
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv,X2,X2inv, L1,L2, Asel1, Asel2, Asel, P
    real(8),dimension(3,3) :: IM, Xrot1, Xrot2
    real(8),dimension(3)   :: Rtras
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    !Duschisky
    real(8),dimension(NDIM,NDIM) :: G
    !T0 - switching effects
    real(8),dimension(3,3) :: T
    !AUXILIAR MATRICES
    real(8),dimension(NDIM,NDIM) :: Aux, Aux2
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: Freq, S1,S1old, S2, Vec, Vec2, mu, Factor, Vec1
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !Shifts
    real(8),dimension(NDIM) :: Delta
    real(8) :: Delta_p, Er_int, Er_crt, Er_qcrt, Er_qint, rmsd, dx,dy,dz
    !Coordinate map
    integer,dimension(NDIM) :: Zmap
    !====================== 

    !====================== 
    !Hessian
    real(8),dimension(:),allocatable :: Hlt

    integer :: error, N
    !====================== 

    !====================== 
    !Auxiliars for LAPACK matrix nversion
    integer :: info
    integer,dimension(NDIM) :: ipiv
    real(8),dimension(NDIM,NDIM) :: work
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=16) :: dummy_char
    real(8) :: Theta, Theta2, Theta3
    ! Messages
    character(len=200) :: msg
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
               I_CNX=17,  &
               O_DUS=20,  &
               O_DIS=21,  &
               O_DMAT=22, &
               O_DUS2=23, &
               O_DIS2=24, &
               O_STAT=25, &
               O_STR =26
    !files
    character(len=10) :: ft ="guess", ftg="guess", fth="guess"
    character(len=200):: inpfile  ="input.fchk", &
                         gradfile ="same", &
                         hessfile ="same", &
                         intfile  ="none", &
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
!     call generic_input_parser(inpfile, "-f" ,"c",&
!                               filetype,"-ft","c",&
!                               )
    call parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,cnx_file,&
                     intfile,rmzfile,def_internal,use_symmetry,vertical,apply_projection_matrix)
    call set_word_upper_case(def_internal)

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
    call generic_strmol_reader(I_INP,ft,state1)
    close(I_INP)
    ! Shortcuts
    Nat = state1%natoms

    ! HESSIAN FILE
    open(I_INP,file=hessfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile)) )
    allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fth,Nat,Hlt,error) 
    close(I_INP)
    ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
    call vibrations_Cart(Nat,state1%atom(:)%X,state1%atom(:)%Y,state1%atom(:)%Z,state1%atom(:)%Mass,Hlt,&
                         Nvib,L1,Freq,error)
    k=0
    do i=1,3*Nat
    do j=1,i
        k=k+1
        Hess(i,j) = Hlt(k)
        Hess(j,i) = Hlt(k)
    enddo 
    enddo
!     deallocate(Hlt)

    ! GRADIENT FILE
    open(I_INP,file=gradfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
    call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
    close(I_INP)

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
    else if (trim(adjustl(symm_file)) /= "NONE") then
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
    if (def_internal/="ZMAT") then 
        print*, "Preliminary Zmat analysis"
        ! Get Zmat first
        def_internal_aux="ZMAT"
        call define_internal_set(state1,def_internal_aux,"none","none",use_symmetry,isym,S_sym,Ns)
        ! Get only the geom, and reuse molecule
        zmatgeom=state1%geom
        ! And reset bonded parameters
        call gen_bonded(state1)
    endif
    call define_internal_set(state1,def_internal,intfile,rmzfile,use_symmetry,isym,S_sym,Ns)
    if (Ns > Nvib) then
        call internals_mapping(state1%geom,zmatgeom,Zmap)
    elseif (Ns < Nvib) then
        print*, "Ns", Ns
        print*, "Nvib", Nvib
        call alert_msg("fatal","Reduced coordinates cases still not implemented")
        ! Need to freeze unused coords to its input values
    endif

    !From now on, we'll use atomic units
    call set_geom_units(state1,"bohr")


    !=================================
    ! INTERNAL COORDINATES
    !=================================
    print'(/,A)', "=============================="
    print'(X,A)', " INTERNAL COORDINATES"
    print'(A)',   "=============================="
    call set_geom_units(state1,"Bohr")
    !SOLVE GF METHOD TO GET NM AND FREQ
    call internal_Wilson(state1,Ns,S1,B1,ModeDef)
    call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,B1,G1)
    if (vertical) then
        call calc_Bder(state1,Ns,Bder,.true.)
    endif
    state2=state1

    ! Get generalized inverse of G
    print*, "Compute generalized inverse of G"
    call generalized_inv(Ns,Nvib,G1,Ginv)
    print*, " Ns  ", Ns
    print*, " Nvib", Nvib
    print*, "Done"
    ! Projection matrix
    P(1:Ns,1:Ns) = matrix_product(Ns,Ns,Ns,G1,Ginv)
    if (verbose>1) call MAT0(6,P,Ns,Ns,"P'")
    print*, "Diag(P')"
    do i=1,Ns
        print'(I4,F8.2)', i, P(i,i)
    enddo
    print*, ""

    if (vertical) then
        ! (Hess is already constructed)
        ! Hs (with the correction)
        ! First get: Hx' = Hx - gs^t\beta
        ! 1. Get gs from gx
        print*, "Get gradient in internal coords"
        call Gradcart2intRed(Nat,Ns,Grad,state1%atom(:)%mass,B1,Ginv)
        print*, "Done"
        ! Project 
        Grad(1:Ns) = matrix_vector_product(Ns,Ns,P,Grad)
        ! 2. Multiply gs^t\beta and
        ! 3. Apply the correction
        ! Bder(i,j,K)^t * gq(K)
        do i=1,3*Nat
        do j=1,3*Nat
            Aux2(i,j) = 0.d0
            do k=1,Ns
                Aux2(i,j) = Aux2(i,j) + Bder(k,i,j)*Grad(k)
            enddo
            Hess(i,j) = Hess(i,j) - Aux2(i,j)
        enddo
        enddo
    endif
    print*, "Get Hessian in internal coords"
    call HessianCart2intRed(Nat,Ns,Hess,state1%atom(:)%mass,B1,Ginv)
    ! Project Hessian
    Hess(1:Ns,1:Ns) = matrix_basisrot(Ns,Ns,P,Hess)
    print*, "Done"

    print*, "Get minimum in internal coords"
    ! Get minimum in internal coordinates
    ! At this point
    !  * Hess has the Hessian  of State2 in internal coords (output from HessianCart2int)
    !  * Grad has the gradient of State2 in internal coords (output from HessianCart2int)
    ! Inverse of the Hessian
    print*, "Compute (generalized) inverse of H"
!     Aux(1:Ns,1:Ns) = inverse_realgen(Ns,Hess)
    i=Nvib
    call generalized_inv(Ns,i,Hess,Aux)
    print*, " Ns  ", Ns
    print*, " Nvib", i
    print*, "Done"
    ! Project inverse Hessian
    Aux(1:Ns,1:Ns) = matrix_basisrot(Ns,Ns,P,Aux)
    ! DeltaS0 = -Hs^1 * gs
    do i=1,Ns
        Delta(i)=0.d0
        do k=1,Ns
            Delta(i) = Delta(i)-Aux(i,k) * Grad(k)
        enddo 
    enddo
    Vec(1:Ns) = Delta(1:Ns)

    ! Print
    k=0
    print*, "DELTA BONDS"
    do i=1,state1%geom%nbonds
        k = k+1
        print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Vec(k), Vec(k)*BOHRtoANGS
    enddo
    print*, "DELTA ANGLES"
    do i=1,state1%geom%nangles
        k = k+1
        print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Vec(k), Vec(k)*180.d0/PI
    enddo
    print*, "DELTA DIHEDRALeS"
    do i=1,state1%geom%ndihed
        k = k+1
        print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Vec(k), Vec(k)*180.d0/PI
    enddo
    do i=1,Ns
        S1(i) = S1(i) + Vec(i)
    enddo
    S2(1:Ns) = S1(1:Ns)

    ! Map to Zmat if needed
    if (adjustl(def_internal) /= "ZMAT") then
        ! From now on, we use the zmatgeom 
        allgeom = state1%geom
        state1%geom = zmatgeom
        S1(1:Nvib) = map_Zmatrix(Nvib,S1,Zmap)
    endif

    ! Get Coordinate by using a Z-matrix selection:
    call zmat2cart(state1,S1)

    !Transform to AA and export coords and put back into BOHR
    call set_geom_units(state1,"Angs")
    open(70,file="minim_harmonic_Inter_LHA.xyz")
    call write_xyz(70,state1)
    close(70)
    print*, ""
    if (adjustl(def_internal) /= "ZMAT") then
        state1%geom = allgeom
        call internal_Wilson(state1,Ns,S1,B1,ModeDef)
        print'(/,X,A)', "Check redundant to non-redundant transformation"
        k=0
        print*, "BONDS", allgeom%nbonds
        do i=1,allgeom%nbonds
            k=k+1
            print'(A,X,3F10.4)', trim(adjustl(ModeDef(k))), &
                                 S2(k), S1(k), S2(k)-S1(k)
        enddo
        print*, "ANGLES", allgeom%nangles
        do i=1,allgeom%nangles
            k=k+1
            print'(A,X,3F10.3)', trim(adjustl(ModeDef(k))), & 
                                 S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
        enddo
        print*, "DIHEDRALS", allgeom%ndihed
        do i=1,allgeom%ndihed
            k=k+1
            print'(A,X,3F10.3)', trim(adjustl(ModeDef(k))), &
                                 S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
        enddo
    endif

! stop
!     if (Ns==Nvib) stop

    ! From internal to Cartesian with the iterative method
    print'(/,X,A)', "Getting Cartesian with iterative method"
    call set_geom_units(state2,"Bohr")
    ! A=M^-1 B^t G^-
    do i=1,Ns
    do j=1,3*Nat
        jj=(j-1)/3+1
        Aux(i,j) = B1(i,j) / state2%atom(jj)%mass / AMUtoAU
    enddo
    enddo
    Aux2(1:3*Nat,1:Ns) = matrix_product(3*Nat,Ns,Ns,Aux,Ginv,tA=.true.)
    rmsd=1.d0
    k=0
    do while (rmsd > 1.d-6)
        k=k+1
        ii=0
        rmsd=0.d0
        do i=1,Nat
            ii=3*(i-1) + 1
            dx = 0.d0
            dy = 0.d0
            dz = 0.d0
            do j=1,Ns
                dx = dx + Aux2(ii+0,j) * Vec(j)
                dy = dy + Aux2(ii+1,j) * Vec(j)
                dz = dz + Aux2(ii+2,j) * Vec(j)
            enddo
            state2%atom(i)%x = state2%atom(i)%x + dx
            state2%atom(i)%y = state2%atom(i)%y + dy
            state2%atom(i)%z = state2%atom(i)%z + dz
            rmsd=rmsd+dx**2+dy**2+dz**2
        enddo
        rmsd=dsqrt(rmsd/Nat)
        print*, k, "rmsd (a.u.) = ", rmsd
        S1old(1:Ns) = S1(1:Ns)
        call verbose_mute()
        call compute_internal(state2,Ns,S1)
        call verbose_continue()
        jj=0
        do j=1,state2%geom%nbonds
            jj=jj+1
            Vec1(jj) = S1(jj) - S1old(jj)
        enddo
        do j=1,state2%geom%nangles
            jj=jj+1
            Vec1(jj) = S1(jj) - S1old(jj)
            if (Vec1(jj) >= 2*pi) then
                Vec1(jj)=Vec1(jj)-2*pi
            elseif (Vec1(jj) <= -2*pi) then
                Vec1(jj)=Vec1(jj)+2*pi
            endif
        enddo
        do j=1,state2%geom%ndihed
            jj=jj+1
            Vec1(jj) = S1(jj) - S1old(jj)
            if (Vec1(jj) >= 2*pi) then
                Vec1(jj)=Vec1(jj)-2*pi
            elseif (Vec1(jj) <= -2*pi) then
                Vec1(jj)=Vec1(jj)+2*pi
            endif
        enddo


        jj=0
        do j=1,state2%geom%nbonds
            jj=jj+1
            print*, Vec(jj), Vec1(jj)
        enddo
print*, ""
        do j=1,state2%geom%nangles
            jj=jj+1
            print*, Vec(jj)*180.d0/PI, Vec1(jj)*180.d0/PI
        enddo
print*, ""
        do j=1,state2%geom%ndihed
            jj=jj+1
            print*, Vec(jj)*180.d0/PI, Vec1(jj)*180.d0/PI
        enddo

        Vec(1:Ns) = Vec(1:Ns) - Vec1(1:Ns)
rmsd=0.d0


    enddo
    print'(X,A,I0,A)', "Iterations converged in ", k, " cyles"

    call set_geom_units(state2,"Angs")
    open(70,file="minim_harmonic_Inter_LHAit.xyz")
    call write_xyz(70,state2)
    close(70)
stop

    print*, ""
    if (Ns /= Nvib) then
        state2%geom = allgeom
        call internal_Wilson(state2,Ns,S1,B1,ModeDef)
        print'(/,X,A)', "Check redundant to non-redundant transformation"
        k=0
        print*, "BONDS", allgeom%nbonds
        do i=1,allgeom%nbonds
            k=k+1
            print'(A,X,3F10.4)', trim(adjustl(ModeDef(k))), &
                                 S2(k), S1(k), S2(k)-S1(k)
        enddo
        print*, "ANGLES", allgeom%nangles
        do i=1,allgeom%nangles
            k=k+1
            print'(A,X,3F10.3)', trim(adjustl(ModeDef(k))), & 
                                 S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
        enddo
        print*, "DIHEDRALS", allgeom%ndihed
        do i=1,allgeom%ndihed
            k=k+1
            print'(A,X,3F10.3)', trim(adjustl(ModeDef(k))), &
                                 S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
        enddo
    endif

    stop

    ! Compute G matrix at the minimum
    !=================================
    ! INTERNAL COORDINATES
    !=================================
    print'(/,A)', "=============================="
    print'(X,A)', " INTERNAL COORDINATES"
    print'(A)',   "=============================="
    call set_geom_units(state1,"Bohr")
    !SOLVE GF METHOD TO GET NM AND FREQ
    call internal_Wilson(state1,Ns,S1,B1,ModeDef)
    call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,B1,G2)

    ! SET REDUNDANT/SYMETRIZED/CUSTOM INTERNAL SETS
!     if (symaddapt) then (implement in an analogous way as compared with the transformation from red to non-red
!     if (Ns > Nvib) then ! Redundant
!         call redundant2nonredundant(Ns,Nvib,G2,Asel1)
        ! Rotate Bmatrix
!         B1(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel1,B1,tA=.true.)
        ! Rotate Gmatrix
        G2(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel1(1:Ns,1:Nvib),G2,counter=.true.)
!     endif

    ! G1 vs G2
    Aux(1:Nvib,1:Nvib) = inverse_realgen(Nvib,G2)
    Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,G1,Aux)
    call MAT0(6,Aux,Nvib,Nvib,"G1 G2^-1")

    print*, "Redoing the Freq analysis"
    call gf_method(Nvib,G2,Hess,L1,Freq,X1,X1inv)


    !===================================
    ! Reorganization energy
    !===================================
    ! Internal-coordinates space
    ! Er = -gs * DeltaS - 1/2 DeltaS^t * Hs * DeltaS
    ! At this point: 
    ! * Grad: in internal coords
    ! * Delta: DeltaS (non-redundant)
    ! * Hess: Hessian in internal coords
    !
    ! Fisrt, compute DeltaS^t * Hs * DeltaS
    Theta=0.d0
    do j=1,Nvib
    do k=1,Nvib
        Theta = Theta + Delta(j)*Delta(k)*Hess(j,k)
    enddo
    enddo
    Er_int = -Theta*0.5d0
    do i=1,Nvib
        Er_int = Er_int - Grad(i)*Delta(i)
    enddo

    ! In Qint-space
    ! Minimization
    ! Q0 = -Lambda^-1 * L^t gs
    ! 
    ! Convert Freq into FC. Store in FC for future use
    do i=1,Nvib
        FC(i) = sign((Freq(i)*2.d0*pi*clight*1.d2)**2/HARTtoJ*BOHRtoM**2*AUtoKG,Freq(i))
        if (FC(i)<0) then
            print*, i, FC(i)
!             FC(i) = -FC(i)
            call alert_msg("warning","A negative FC found")
        endif
    enddo
    ! Lambda^-1 * L1^t
    do i=1,Nvib
        Aux(i,1:Nvib) = L1(1:Nvib,i) / FC(i)
    enddo
    ! -[Lambda^-1 * L1^t] * gs
    do i=1,Nvib
        Q0(i)=0.d0
        do k=1,Nvib
            Q0(i) = Q0(i) - Aux(i,k) * Grad(k)
        enddo
    enddo
    !===================================
    ! Reorganization energy
    !===================================
    ! Normal-mode space
    ! Er = -L1^t gs * Q0 - 1/2 * Q0^t * Lambda * Q0
    ! At this point: 
    ! * Grad: in internal coords
    ! * Q0: DeltaQint 
    ! * FC: diagonal force constants
    Er_qint = 0.d0
    do i=1,Nvib
        ! Compute gQ(i) = L1^t * gs
        Theta = 0.d0
        do j=1,Nvib
            Theta =  Theta + L1(j,i)*Grad(j)
        enddo
        Er_qint = Er_qint - Theta * Q0(i) - 0.5d0 * FC(i) * Q0(i)**2
    enddo


    ! PRINT
    print*, "CARTESIAN COORDINATES"
    print'(X,A,F12.6)',   "Reorganization energy (AU) = ", Er_crt
    print'(X,A,F12.6,/)', "Reorganization energy (eV) = ", Er_crt*HtoeV
    print*, "NORMAL-MODE COORDINATES (derived in Cartesian)"
    print'(X,A,F12.6)',   "Reorganization energy (AU) = ", Er_qcrt
    print'(X,A,F12.6,/)', "Reorganization energy (eV) = ", Er_qcrt*HtoeV
    print*, "INTERNAL COORDINATES"
    print'(X,A,F12.6)',   "Reorganization energy (AU) = ", Er_int
    print'(X,A,F12.6,/)', "Reorganization energy (eV) = ", Er_int*HtoeV
    print*, "NORMAL-MODE COORDINATES (derived in internal)"
    print'(X,A,F12.6)',   "Reorganization energy (AU) = ", Er_qint
    print'(X,A,F12.6,/)', "Reorganization energy (eV) = ", Er_qint*HtoeV


    call summary_alerts

    call cpu_time(tf)
    write(6,'(A,F12.3)') "CPU (s) for internal vib analysis: ", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,cnx_file,& 
                           intfile,rmzfile,def_internal,use_symmetry,vertical,apply_projection_matrix)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,&
                                          intfile,rmzfile,def_internal,cnx_file
        logical,intent(inout)          :: use_symmetry, vertical,apply_projection_matrix
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

                case ("-vert")
                    vertical=.true.
                case ("-novert")
                    vertical=.false.

                case ("-cnx") 
                    call getarg(i+1, cnx_file)
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

                case ("-prj-tr")
                    apply_projection_matrix=.true.
                case ("-noprj-tr")
                    apply_projection_matrix=.false.
        
                case ("-h")
                    need_help=.true.

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


       !Print options (to stderr)
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '           L H A   o p t i mi z e r '    
        write(6,'(/,A)') '  Displace structure from vertical to adiabatic geoms  '   
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag           Description                   Value'
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,*)       '-f              Input (vertical) file        ', trim(adjustl(inpfile))
        write(6,*)       '-ft             \_ FileType                  ', trim(adjustl(ft))
        write(6,*)       '-fhess          Hessian file                 ', trim(adjustl(hessfile))
        write(6,*)       '-fth            \_ FileType                  ', trim(adjustl(fth))
        write(6,*)       '-fgrad          Gradient File                ', trim(adjustl(gradfile))
        write(6,*)       '-ftg            \_ FileType                  ', trim(adjustl(ftg))
        write(6,*)       '-intmode        Internal set:[zmat|sel|all]  ', trim(adjustl(def_internal))
        write(6,*)       '-cnx           Connectivity [filename|guess] ', trim(adjustl(cnx_file))
        write(6,*)       '-intfile        File with ICs (for "sel")    ', trim(adjustl(intfile))
        write(6,*)       '-rmzfile        File deleting ICs from Zmat  ', trim(adjustl(rmzfile))
        write(6,*)       '-[no]vert       Vertical model              ',  vertical
        write(6,*)       '-[no]prj-tr     Apply projection matrix     ',  apply_projection_matrix
        write(6,*)       '-h              Display this help           ',  need_help
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                   'Verbose level:  ', verbose        
        write(6,'(A)')   '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program optimizer

