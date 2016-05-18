program NewtonRaphson_step


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
               analytic_Bder=.true. ,&
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
    real(8),dimension(NDIM,NDIM) :: B1,B2, B, G1,G2, Ginv
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad, FC, Q0, Grad0
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv,X2,X2inv, L1,L2, Asel1, Asel2, Asel, P, Hess0, Pc
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
    real(8),dimension(NDIM) :: Freq, S1, S2, Vec, Vec2, mu, Factor, Vec1, Vec3
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !Shifts
    real(8),dimension(NDIM) :: Delta
    real(8) :: Delta_p, Er_int, Er_crt, Er_qcrt, Er_qint
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
    real(8) :: Theta, Theta2, Theta3, Theta4, thr
    ! Messages
    character(len=200) :: msg
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, k90,k95,k99, nn, imin, imax,&
               i1,i2,i3,i4, iter
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
                     intfile,rmzfile,def_internal,use_symmetry,vertical,apply_projection_matrix,analytic_Bder)
    call set_word_upper_case(def_internal)
    call set_word_upper_case(symm_file)

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
    Hess0(1:3*Nat,1:3*Nat)=Hess(1:3*Nat,1:3*Nat)
!     deallocate(Hlt)

    ! GRADIENT FILE
    open(I_INP,file=gradfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
    call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
    Grad0(1:3*Nat)=Grad(1:3*Nat)
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


    !==============================
    ! CARTESIAN COORDINATES
    !==============================
    print'(/,A)', "=============================="
    print'(X,A)', " CARTESIAN COORDINATES"
    print'(A)',   "=============================="
    ! Get minimum in Cartesian coordinates: x0 = - F^-1 grad
    Aux(1:3*Nat,1:3*Nat) = inverse_realsym(3*Nat,Hess)
    ! Matrix x vector 
    do i=1, 3*Nat
        Vec(i)=0.d0
        do k=1,3*Nat
            Vec(i) = Vec(i) - Aux(i,k) * Grad(k)
        enddo 
    enddo
    print*, "DELTA R (x,y,z), Angstrong"
    do i=1,Nat 
        j=3*i-2
        print'(I3,3X, 3F10.3)', i, Vec(j)*BOHRtoANGS, Vec(j+1)*BOHRtoANGS, Vec(j+2)*BOHRtoANGS
    enddo
    print*, ""

    state2=state1
    do i=1,Nat 
        j=3*i-2
        state2%atom(i)%x = (state1%atom(i)%x + Vec(j+0))*BOHRtoANGS
        state2%atom(i)%y = (state1%atom(i)%y + Vec(j+1))*BOHRtoANGS
        state2%atom(i)%z = (state1%atom(i)%z + Vec(j+2))*BOHRtoANGS
    enddo
    open(70,file="minim_harmonic_Cart.xyz")
    call write_xyz(70,state2)
    close(70)

    ! Check the rotation of the Ekart frame
    print'(/,A)', "------------------------------------------------------------"
    print'(X,A)', "ESTIMATION OF THE MOLECULAR TRASLATION/ROTATION (CARTESIAN) "
    print'(A)',   "------------------------------------------------------------"

    call set_geom_units(state1,"Angs")
    call set_geom_units(state2,"Angs")

    !Traslation:
    call get_com(state1)
    call get_com(state2)
    Rtras(1) = state1%comX - state2%comX
    Rtras(2) = state1%comY - state2%comY
    Rtras(3) = state1%comZ - state2%comZ
    call print_vector(6,Rtras,3,"Traslation between Vertical and Adiabatic")

    ! The rotation can be computed from the diagonalization of the matrix
    ! of moment of inertia for each geometry
    call inertia(state1,IM)
    call diagonalize_full(IM(1:3,1:3),3,Xrot1(1:3,1:3),Vec2(1:3),"lapack")
    if (verbose>1) &
     call MAT1(6,Xrot1,Vec2,3,3,"Xrot (state1)")
    call inertia(state2,IM)
    call diagonalize_full(IM(1:3,1:3),3,Xrot2(1:3,1:3),Vec2(1:3),"lapack")
    if (verbose>1) &
     call MAT1(6,Xrot2,Vec2,3,3,"Xrot (state2)")
    !
    ! The rotation from one geometry to the other is then:
    ! Rot = Xrot1^t  Xrot2
    Xrot1(1:3,1:3) = matrix_product(3,3,3,Xrot1,Xrot2,tA=.true.)

    call MAT0(6,Xrot1,3,3,"Rotation between Vertical and Adiabatic")

    ! Check the vibrational analysis at the state2 estimated geom
    print'(/,A)', "-------------------------------------------------"
    print'(X,A)', "VIBRATIONAL ANALYSIS WITH STATE2 GEOM (ESTIMATED)"
    print'(A)',   "-------------------------------------------------"
    call vibrations_Cart(Nat,state2%atom(:)%X,state2%atom(:)%Y,state2%atom(:)%Z,state2%atom(:)%Mass,Hlt,&
                         Nvib,L1,Vec2,error)
    deallocate(Hlt)

    !-------------------------------
    ! Reorganization energy
    !-------------------------------
    ! Cartesian-coordinates space
    ! Er = -gx * DeltaX - 1/2 DeltaX^t * Hx * DeltaX
    ! At this point: 
    ! * Grad: in Cartesian coords
    ! * Vec: DeltaX 
    ! * Hess: Hessian in Cartesian coords
    !
    ! Fisrt, compute DeltaX^t * Hs * DeltaX
    Theta=0.d0
    do j=1,3*Nat
    do k=1,3*Nat
        Theta = Theta + Vec(j)*Vec(k)*Hess(j,k)
    enddo
    enddo
    Er_crt = -Theta*0.5d0
    do i=1,3*Nat
        Er_crt = Er_crt - Grad(i)*Vec(i)
    enddo

    ! In Qcart-space
    ! First convert L to Lcart
    call Lmwc_to_Lcart(Nat,Nvib,state1%atom(:)%mass,L1,L1,error)
    ! Minimization
    ! Q0 = -Lambda^-1 * L^t gx
    ! 
    ! Convert Freq into FC. Store in FC for future use
    do i=1,Nvib
        FC(i) = sign((Freq(i)*2.d0*pi*clight*1.d2)**2/HARTtoJ*BOHRtoM**2*AUtoKG,Freq(i))
        if (FC(i)<0) then
            print*, i, FC(i)
            call alert_msg("warning","A negative FC found")
        endif
    enddo
    ! Lambda^-1 * L1^t
    do i=1,Nvib
        Aux(i,1:3*Nat) = L1(1:3*Nat,i) / FC(i)
    enddo
    ! -[Lambda^-1 * L1^t] * gx
    do i=1,Nvib
        Q0(i)=0.d0
        do k=1,3*Nat
            Q0(i) = Q0(i) - Aux(i,k) * Grad(k)
        enddo
    enddo
    !-------------------------
    ! Reorganization energy
    !-------------------------
    ! Normal-mode space
    ! Er = -L1^t gx * Q0 - 1/2 * Q0^t * Lambda * Q0
    ! At this point: 
    ! * Grad: in Cartesian coords
    ! * Q0: DeltaQ 
    ! * FC: diagonal force constants
    Er_qcrt = 0.d0
    do i=1,Nvib
        ! Compute gQ(i) = L1^t * gx
        Theta = 0.d0
        do j=1,3*Nat
            Theta =  Theta + L1(j,i)*Grad(j)
        enddo
        Er_qcrt = Er_qcrt - Theta * Q0(i) - 0.5d0 * FC(i) * Q0(i)**2
    enddo



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
        call calc_Bder(state1,Ns,Bder,analytic_Bder)
    endif

    ! SET REDUNDANT/SYMETRIZED/CUSTOM INTERNAL SETS
!     if (symaddapt) then (implement in an analogous way as compared with the transformation from red to non-red
!     if (Ns > Nvib) then ! Redundant
        call redundant2nonredundant(Ns,Nvib,G1,Asel1)
        ! Rotate Bmatrix
        B1(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel1,B1,tA=.true.)
        ! Rotate Gmatrix
        G1(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel1(1:Ns,1:Nvib),G1,counter=.true.)
        ! Rotate Bders
        if (vertical) then
            do j=1,3*Nat
                Bder(1:Nvib,j,1:3*Nat) =  matrix_product(Nvib,3*Nat,Ns,Asel1,Bder(1:Ns,j,1:3*Nat),tA=.true.)
            enddo
        endif

    P(1:Ns,1:Ns) = matrix_product(Ns,Ns,Nvib,Asel1,Asel1,tB=.true.)
    call MAT0(6,P,Ns,Ns,"Prj(1)")
!     endif

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

    if (vertical) then
        ! (Hess is already constructed)
        ! Hs (with the correction)
        ! First get: Hx' = Hx - gs^t\beta
        ! 1. Get gs from gx
        call Gradcart2int(Nat,Nvib,Grad,state1%atom(:)%mass,B1,G1)
        ! 2. Multiply gs^t\beta and
        ! 3. Apply the correction
        ! Bder(i,j,K)^t * gq(K)
        do i=1,3*Nat
        do j=1,3*Nat
            Aux2(i,j) = 0.d0
            do k=1,Nvib
                Aux2(i,j) = Aux2(i,j) + Bder(k,i,j)*Grad(k)
            enddo
            Hess(i,j) = Hess(i,j) - Aux2(i,j)
        enddo
        enddo
    endif
    if (apply_projection_matrix) then
        ! Project out rotation and translation
        Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,P,Hess)
    endif
    call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B1,G1)
    call gf_method(Nvib,G1,Hess,L1,Freq,X1,X1inv)

    ! Get minimum in internal coordinates
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
!     if (Nvib<Ns) then
        !Transform Delta' into Delta 
        ! Delta = A Delta'
        do i=1,Ns
            Vec(i) = 0.d0
            do k=1,Nvib
                Vec(i) = Vec(i) + Asel1(i,k)*Delta(k)
            enddo
        enddo

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
    print*, "DELTA DIHEDRALS"
    do i=1,state1%geom%ndihed
        k = k+1
        print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Vec(k), Vec(k)*180.d0/PI
    enddo
    do i=1,Ns
        S1(i) = S1(i) + Vec(i)
    enddo
    S2(1:Ns) = S1(1:Ns)

    ! Map to Zmat if needed
    state2=state1
    ! From now on, we use the zmatgeom 
    allgeom = state1%geom
    if (Ns /= Nvib) then
        state1%geom = zmatgeom
        S1(1:Nvib) = map_Zmatrix(Nvib,S1,Zmap)
    endif
    call zmat2cart(state1,S1)

    !Transform to AA and export coords and put back into BOHR
    call set_geom_units(state1,"Angs")
    open(70,file="minim_harmonic_Int_Zmat.xyz")
    call write_xyz(70,state1)
    close(70)
    call set_geom_units(state1,"BOHR")
    print*, ""
!     if (Ns /= Nvib) then
        state1%geom = allgeom
        call verbose_mute()
        call compute_internal(state1,Ns,S1,ModeDef)
        call verbose_continue()
        print'(/,X,A)', "Check redundant to non-redundant transformation"
        k=0
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "BONDS", allgeom%nbonds
        do i=1,allgeom%nbonds
            k=k+1
            print'(3F10.4)', S2(k), S1(k), S2(k)-S1(k)
            Theta = (S2(k)-S1(k)) * BOHRtoANGS
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
        enddo
        Theta3 = dsqrt(Theta3/allgeom%nbonds)
        print'(X,5(A,F12.6))',  "(1)MaxDev-bond(AA) = ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "ANGLES", allgeom%nangles
        do i=1,allgeom%nangles
            k=k+1
            print'(3F10.3)', S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
            Theta = (S2(k)-S1(k)) * 180.d0/PI
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
        enddo
        Theta3 = dsqrt(Theta3/allgeom%nangles)
        print'(X,5(A,F12.6))',  "(1)MaxDev-angl(deg)= ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "DIHEDRALS", allgeom%ndihed
        do i=1,allgeom%ndihed
            k=k+1
            print'(3F10.3)', S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
            Theta = (S2(k)-S1(k)) * 180.d0/PI
            if (abs(Theta)>abs(Theta+360.d0)) then
                Theta = Theta+360.d0
            else if (abs(Theta)>abs(Theta-360.d0)) then
                Theta = Theta-360.d0
            endif
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
        enddo
        Theta3 = dsqrt(Theta3/allgeom%ndihed)
        print'(X,5(A,F12.6))',  "(1)MaxDev-dihe(deg)= ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3
!     endif



    ! Perform the step with the iterative procedure
    call set_geom_units(state2,"BOHR")
    ! A = M^-1 B^t G^-
    ! Compute Ginv(red) from non-red one
    Ginv(1:Nvib,1:Nvib) = 0.d0
    do i=1,Nvib
        Ginv(i,i) = 1.d0/G1(i,i)
    enddo
    Ginv(1:Ns,1:Ns) = matrix_basisrot(Ns,Nvib,Asel1(1:Ns,1:Nvib),Ginv,counter=.false.)
    ! Rotate back Bmatrix
    B1(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Nvib,Asel1,B1)

    ! From scratch
    state2%geom = allgeom
    state1=state2

    call intshif2cart(state2,Vec)



    call set_geom_units(state2,"Angs")
    open(70,file="minim_harmonic_Int_it.xyz")
    call write_xyz(70,state2)
    close(70)

    ! Export fch with the Hessian transfomed back at the minimum
    open(70,file="minim_harmonic_Int_it.fchk")
    call generic_structure_writer(70,'fchk',Nat,state2%atom(1:Nat)%X,&
                                                state2%atom(1:Nat)%Y,&
                                                state2%atom(1:Nat)%Z,&
                                                state2%atom(1:Nat)%Mass,&
                                                state2%atom(1:Nat)%Name,&
                                                error,"NewtonRaphson step")
    
    call internal_Wilson(state2,Ns,S2,B2,ModeDef)
    Hess(1:Ns,1:Ns) = matrix_basisrot(Ns,Nvib,Asel1(1:Ns,1:Nvib),Hess,counter=.false.)
    Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,Ns,B2,Hess,counter=.true.)
    N = 3*Nat*(3*Nat+1)/2
    allocate(Hlt(1:N))
    Hlt(1:N) = Hess_to_Hlt(3*Nat,Hess)
    call write_fchk(70,"Cartesian Force Constants","R",N,Hlt,(/0/),error)
    deallocate(Hlt)
    close(70)
    call set_geom_units(state2,"BOHR")

!     if (Ns /= Nvib) then
        state2%geom = allgeom
        call verbose_mute()
        call compute_internal(state2,Ns,S1,ModeDef)
        call verbose_continue()
        print'(/,X,A)', "Check redundant to non-redundant transformation"
        k=0
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "BONDS", allgeom%nbonds
        do i=1,allgeom%nbonds
            k=k+1
            print'(3F10.4)', S2(k), S1(k), S2(k)-S1(k)
            Theta = (S2(k)-S1(k)) * BOHRtoANGS
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
        enddo
        Theta3 = dsqrt(Theta3/allgeom%nbonds)
        print'(X,5(A,F12.6))',  "(2)MaxDev-bond(AA) = ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "ANGLES", allgeom%nangles
        do i=1,allgeom%nangles
            k=k+1
            print'(3F10.3)', S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
            Theta = (S2(k)-S1(k)) * 180.d0/PI
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
        enddo
        Theta3 = dsqrt(Theta3/allgeom%nangles)
        print'(X,5(A,F12.6))',  "(2)MaxDev-angl(deg)= ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "DIHEDRALS", allgeom%ndihed
        do i=1,allgeom%ndihed
            k=k+1
            print'(3F10.3)', S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
            Theta = (S2(k)-S1(k)) * 180.d0/PI
            if (abs(Theta)>abs(Theta+360.d0)) then
                Theta = Theta+360.d0
            else if (abs(Theta)>abs(Theta-360.d0)) then
                Theta = Theta-360.d0
            endif
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
        enddo
        Theta3 = dsqrt(Theta3/allgeom%ndihed)
        print'(X,5(A,F12.6))',  "(2)MaxDev-dihe(deg)= ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3
!     endif








    ! Compute G matrix at the minimum
    !=================================
    ! INTERNAL COORDINATES
    !=================================
!     print'(/,A)', "=============================="
!     print'(X,A)', " INTERNAL COORDINATES"
!     print'(A)',   "=============================="
!     call set_geom_units(state1,"Bohr")
!     !SOLVE GF METHOD TO GET NM AND FREQ
!     call verbose_mute()
!     call internal_Wilson(state1,Ns,S1,B1,ModeDef)
!     call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,B1,G2)

    ! SET REDUNDANT/SYMETRIZED/CUSTOM INTERNAL SETS
!     if (symaddapt) then (implement in an analogous way as compared with the transformation from red to non-red
!     if (Ns > Nvib) then ! Redundant
!         call redundant2nonredundant(Ns,Nvib,G2,Asel1)
        ! Rotate Bmatrix
!         B1(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel1,B1,tA=.true.)
        ! Rotate Gmatrix
!         G2(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel1(1:Ns,1:Nvib),G2,counter=.true.)
!     endif

!     ! G1 vs G2
!     Aux(1:Nvib,1:Nvib) = inverse_realgen(Nvib,G2)
!     Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,G1,Aux)
!     call MAT0(6,Aux,Nvib,Nvib,"G1 G2^-1")
! 
!     print*, "Redoing the Freq analysis"
!     call gf_method(Nvib,G2,Hess,L1,Freq,X1,X1inv)


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
    print*, ""
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


!********************************************
        !Transform back to Delta'
        do i=1,Ns
            Vec1(i) = 0.d0
            do k=1,Ns
                Vec1(i) = Vec1(i) + Asel1(k,i)*Vec(k)
            enddo
        enddo

        print*, "Non-zero"
        do i=1,Nvib
            print'(3(ES14.6,2X))', Delta(i), Vec1(i), Delta(i)-Vec1(i)
        enddo
        print*, "Zero"
        do i=Nvib+1,Ns
            print'(3(ES14.6,2X))', Vec1(i)
        enddo

!----------- From scratch using Ginv as Pulay and Peng -------------
    Hess(1:3*Nat,1:3*Nat) = Hess0(1:3*Nat,1:3*Nat)
    Grad(1:3*Nat) = Grad0(1:3*Nat)
    call internal_Wilson(state1,Ns,S1,B1,ModeDef)
    call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,B1,G1)
    if (vertical) then
        call calc_Bder(state1,Ns,Bder,analytic_Bder)
    endif
    print*, "Get Ginv and P and alpha(1-P)"
    call generalized_inv(Ns,Nvib,G1,Ginv)
    P(1:Ns,1:Ns) = matrix_product(Ns,Ns,Ns,G1,Ginv)
!     Aux(1:Ns,1:Ns) = matrix_product(Ns,Ns,Nvib,Asel1,Asel1,tB=.true.)
!     call MAT0(6,P,Ns,Ns,"P(1)")
!     call MAT0(6,Aux,Ns,Ns,"P(2)")
!     call MAT0(6,(Aux-P)*1e5,Ns,Ns,"[P(1)-P(2)] x1e5")
!     Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Ns,Asel1,Asel1,tA=.true.)
!     call MAT0(6,Aux,Nvib,Nvib,"A^t A")
! stop

    Pc(1:Ns,1:Ns) = 0.d0
    do i=1,Ns
        Pc(i,i) = 1d0
    enddo
    Pc(1:Ns,1:Ns) = Pc(1:Ns,1:Ns)-P(1:Ns,1:Ns)
    if (vertical) then
        ! (Hess is already constructed)
        ! Hs (with the correction)
        ! First get: Hx' = Hx - gs^t\beta
        ! 1. Get gs from gx
        call Gradcart2int_red(Nat,Ns,Grad,state1%atom(:)%mass,B1,G1)
        do i=1,Ns
            Vec2(i) = 0.d0
            do k=1,Ns
                Vec2(i) = Vec2(i) + P(i,k)*Grad(k)
            enddo
        enddo
!         Vec3(1:Ns) = Grad(1:Ns) 

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
    print'(/,X,A)', "Convert Hessian from Cart to Int(red)"
    call HessianCart2int_red(Nat,Ns,Hess,state1%atom(:)%mass,B1,G1)
!     Aux(1:Ns,1:Ns) = matrix_product(Ns,Ns,Ns,P,P)
!     call MAT0(6,P,Ns,Ns,"P")
!     call MAT0(6,Aux,Ns,Ns,"P^2")
!     print'(/,X,A)', "Project Hessian"
!     Aux(1:Ns,1:Ns) = matrix_basisrot(Ns,Ns,P,Hess)
!     call MAT0(6,Hess*1e5,Ns,Ns,"Hess x1e5")
!     call MAT0(6,Aux*1e5,Ns,Ns,"PHessP x1e5")
!     call MAT0(6,(Aux-Hess)*1e5,Ns,Ns,"(PHessP-Hess) x1e5")


!     Hess(1:Ns,1:Ns) = Hess(1:Ns,1:Ns) + 1.d3*Pc(1:Ns,1:Ns)
    print'(/,X,A)', "Inverse projected Hessian"
    call generalized_inv(Ns,Nvib,Hess,Aux)
!     Aux(1:Ns,1:Ns) = inverse_realgen(Ns,Hess)
!     print'(/,X,A)', "Project back inverse Hessian"
!     Aux(1:Ns,1:Ns) = matrix_basisrot(Ns,Ns,P,Aux)
! print*, "Nvib", Nvib



    ! DeltaS0 = -Hs^1 * gs
    do i=1,Ns
        Vec2(i)=0.d0
        do k=1,Ns
            Vec2(i) = Vec2(i)-Aux(i,k) * Grad(k)
        enddo 
    enddo

        do i=1,Ns
            Vec3(i) = 0.d0
            do k=1,Ns
                Vec3(i) = Vec3(i) + Asel1(k,i)*Vec2(k)
            enddo
        enddo

        print*, "Non-zero"
        do i=1,Nvib
            print'(3(ES14.6,2X))', Delta(i), Vec1(i), Vec3(i)
        enddo
        print*, "Zero"
        do i=Nvib+1,Ns
            print'(3(ES14.6,2X))', Vec1(i), Vec3(i)
        enddo

        print*, "Redundant"
        do i=1,Ns
            print'(3(ES14.6,2X))', Vec(i), Vec2(i)
        enddo
         

!         stop

!     else
!         Vec(1:Nvib)=Delta(1:Nvib)
!     endif
    Vec(1:Ns) = Vec2(1:Ns)

    ! From scratch
    state1%geom = allgeom

    call intshif2cart(state1,Vec,maxiter_set=100)



    call set_geom_units(state1,"Angs")
    open(70,file="minim_harmonic_Int_it2.xyz")
    call write_xyz(70,state1)
    close(70)
    call set_geom_units(state1,"BOHR")
!     if (Ns /= Nvib) then
        state1%geom = allgeom
        call verbose_mute()
        call compute_internal(state1,Ns,S1,ModeDef)
        call verbose_continue()
        print'(/,X,A)', "Check redundant to non-redundant transformation"
        k=0
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "BONDS", allgeom%nbonds
        do i=1,allgeom%nbonds
            k=k+1
            print'(3F10.4)', S2(k), S1(k), S2(k)-S1(k)
            Theta = (S2(k)-S1(k)) * BOHRtoANGS
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
        enddo
        Theta3 = dsqrt(Theta3/allgeom%nbonds)
        print'(X,5(A,F12.6))',  "(3)MaxDev-bond(AA) = ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "ANGLES", allgeom%nangles
        do i=1,allgeom%nangles
            k=k+1
            print'(3F10.3)', S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
            Theta = (S2(k)-S1(k)) * 180.d0/PI
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
        enddo
        Theta3 = dsqrt(Theta3/allgeom%nangles)
        print'(X,5(A,F12.6))',  "(3)MaxDev-angl(deg)= ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3
        Theta=0.d0
        Theta2 = 0.d0
        Theta3 = 0.d0
        Theta4 = 1.d10
        print*, "DIHEDRALS", allgeom%ndihed
        do i=1,allgeom%ndihed
            k=k+1
!             print'(3F10.3)', S2(k)*180.d0/PI, S1(k)*180.d0/PI, (S2(k)-S1(k))*180.d0/PI
            Theta = (S2(k)-S1(k)) * 180.d0/PI
            if (abs(Theta)>abs(Theta+360.d0)) then
                Theta = Theta+360.d0
            else if (abs(Theta)>abs(Theta-360.d0)) then
                Theta = Theta-360.d0
            endif
            Theta2 = max(Theta2,abs(Theta))
            Theta4 = min(Theta4,abs(Theta))
            Theta3 = Theta3 + Theta**2
            print'(A,3F10.3)', ModeDef(k), S2(k)*180.d0/PI, S1(k)*180.d0/PI, Theta
        enddo
        Theta3 = dsqrt(Theta3/allgeom%ndihed)
        print'(X,5(A,F12.6))',  "(3)MaxDev-dihe(deg)= ", Theta2, "   MinDev", Theta4, "   RMSD", Theta3

!********************************************

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,cnx_file,& 
                           intfile,rmzfile,def_internal,use_symmetry,vertical,apply_projection_matrix,analytic_Bder)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,&
                                          intfile,rmzfile,def_internal,cnx_file
        logical,intent(inout)          :: use_symmetry, vertical,apply_projection_matrix
        logical,intent(inout)          :: analytic_Bder
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        ! iargc type must be specified with implicit none (strict compilation)
        integer :: iargc

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


       !Print options (to stderr)
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '        V E R T I C A L 2 A D I A B A T I C '    
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
        write(6,*)       '-sym            Use symmetry                ',  use_symmetry
        write(6,*)       '-h              Display this help           ',  need_help
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                   'Verbose level:  ', verbose        
        write(6,'(A)')   '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program NewtonRaphson_step

