program vertical2adiabatic


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
    use xyz_manage
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
               vertical=.true.
    character(len=4) :: def_internal='zmat'
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1,state2
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !B and G matrices
    real(8),dimension(NDIM,NDIM) :: B1,B2, B, G1,G2
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad, FC, Q0
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv,X2,X2inv, L1,L2, Asel1, Asel2, Asel
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
    real(8),dimension(NDIM) :: Freq, S1, S2, Vec, Vec2, mu, Factor
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !Shifts
    real(8),dimension(NDIM) :: Delta
    real(8) :: Delta_p, Er_int, Er_crt, Er_qcrt, Er_qint
    !====================== 

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
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
!     call generic_input_parser(inpfile, "-f" ,"c",&
!                               filetype,"-ft","c",&
!                               )
    call parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,intfile,rmzfile,def_internal,use_symmetry,vertical)
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
    allocate(A(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fth,Nat,A,error) 
    close(I_INP)
    ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
    call vibrations_Cart(Nat,state1%atom(:)%X,state1%atom(:)%Y,state1%atom(:)%Z,state1%atom(:)%Mass,A,&
                         Nvib,L1,Freq,error)
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
    call define_internal_set(state1,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns)

    !From now on, we'll use atomic units
    call set_geom_units(state1,"bohr")


    ! GET MINIMUM IN CARTESIAN COORDINATES: x0 = F^-1 grad
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
        state2%atom(i)%x = (state1%atom(i)%x - Vec(j+0))*BOHRtoANGS
        state2%atom(i)%y = (state1%atom(i)%y - Vec(j+1))*BOHRtoANGS
        state2%atom(i)%z = (state1%atom(i)%z - Vec(j+2))*BOHRtoANGS
    enddo
    open(70,file="minim_harmonic_Cart.xyz")
    call write_xyz(70,state2)
    close(70)

    !===================================
    ! Reorganization energy
    !===================================
    ! Cartesian-coordinates space
    ! Er = -gx * DeltaX - 1/2 DeltaX^t * Hx * DeltaX
    ! At this point: 
    ! * Grad: in Cartesian coords
    ! * Vec: DeltaX 
    ! * Hess: Hessian in Cartesian coords
    !
    ! Fisrt, compute DeltaS^t * Hs * DeltaS
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
!             FC(i) = -FC(i)
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
    !===================================
    ! Reorganization energy
    !===================================
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


    ! INTERNAL COORDINATES

    !SOLVE GF METHOD TO GET NM AND FREQ
    call internal_Wilson(state1,Nvib,S1,B1,ModeDef)
    call internal_Gmetric(Nat,Nvib,state1%atom(:)%mass,B1,G1)
    if (vertical) then
        call calc_Bder(state1,Nvib,Bder)
        call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B1,G1,Grad=Grad,Bder=Bder)
    else
        call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B1,G1)
    endif
    call gf_method(Nvib,G1,Hess,L1,Freq,X1,X1inv)

    !Compute new state_file for 2
    ! T2(g09) = mu^1/2 m B^t G2^-1 L2
    ! Compute G1^-1 (it is X1inv * X1inv
    Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,X1inv,X1inv)
    ! Compute B1^t G1^-1
    do i=1,3*Nat
    do j=1,Nvib
        Aux2(i,j) = 0.d0
        do k=1,Nvib
         Aux2(i,j)=Aux2(i,j)+B1(k,i)*Aux(k,j)
        enddo
    enddo
    enddo
    ! Compute [B1^t G2^-1] L2
    Aux(1:3*Nat,1:Nvib) = matmul(Aux2(1:3*Nat,1:Nvib),L1(1:Nvib,1:Nvib))
    ! Compute mu^1/2 m [B^t G2^-1 L2] (masses are in UMA in the fchk)
!     print*, state2%atom(1)%name, state2%atom(1)%mass
!     print*, mu(1), mu(Nvib)
    i=0
    do k=1,Nat
    do kk=1,3
    i=i+1
    do j=1,Nvib
        Aux2(i,j) = 1.d0/&!dsqrt(mu(j)*UMAtoAU)        / &
                    state1%atom(k)%mass/UMAtoAU * &
                    Aux(i,j)
    enddo
    enddo
    enddo
    !Compute reduced masses
    do j=1,Nvib
    mu(j)=0.d0
    do i=1,3*Nat
        mu(j)=mu(j)+Aux2(i,j)**2
    enddo
    mu(j) = 1.d0/mu(j)
    enddo
    !Normalize with mu
    i=0
    do k=1,Nat
    do kk=1,3
    i=i+1
    do j=1,Nvib
        Aux2(i,j) = Aux2(i,j)*dsqrt(mu(j))
    enddo
    enddo
    enddo
    !Print state
    open(O_STAT,file="state_file_1")
    do i=1,Nat
        write(O_STAT,*) state1%atom(i)%x*BOHRtoANGS
        write(O_STAT,*) state1%atom(i)%y*BOHRtoANGS
        write(O_STAT,*) state1%atom(i)%z*BOHRtoANGS
    enddo
    do i=1,3*Nat
    do j=1,Nvib
        write(O_STAT,*) Aux2(i,j)
    enddo
    enddo
    do j=1,Nvib
        write(O_STAT,'(F12.5)') Freq(j)
    enddo
    close(O_STAT)

    if (verbose>1) then
    print*, "B1=", B1(1,1)
    do i=1,Nvib
        print'(100(F8.3,2X))', B1(i,1:Nvib)
    enddo
    endif


    ! GET MINIMUM IN INTERNAL COORDINATES
    if (verbose>1) then
        Vec(1:Nvib) = (/(Hess(i,i), i=1,Nvib)/)
        call print_vector(6,Vec,Nvib,"Diagonal FC (internal)")
        call print_vector(6,Grad,Nvib,"Gradient (internal)")
    endif
    Aux(1:Nvib,1:Nvib) = inverse_realgen(Nvib,Hess)
    ! Matrix x vector 
    do i=1,Nvib
        Vec(i)=0.d0
        do k=1,Nvib
            Vec(i) = Vec(i)-Aux(i,k) * Grad(k)
        enddo 
    enddo
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
    do i=1,Nvib
        S1(i) = S1(i) + Vec(i)
    enddo
    call zmat2cart(state1,S1)
    !Transform to AA and export coords and put back into BOHR
    state1%atom(1:Nat)%x = state1%atom(1:Nat)%x*BOHRtoANGS
    state1%atom(1:Nat)%y = state1%atom(1:Nat)%y*BOHRtoANGS
    state1%atom(1:Nat)%z = state1%atom(1:Nat)%z*BOHRtoANGS
    open(70,file="minim_harmonic_Inter.xyz")
    call write_xyz(70,state1)
    close(70)
    print*, ""

    !===================================
    ! Reorganization energy
    !===================================
    ! Internal-coordinates space
    ! Er = -gs * DeltaS - 1/2 DeltaS^t * Hs * DeltaS
    ! At this point: 
    ! * Grad: in internal coords
    ! * Vec: DeltaS 
    ! * Hess: Hessian in internal coords
    !
    ! Fisrt, compute DeltaS^t * Hs * DeltaS
    Theta=0.d0
    do j=1,Nvib
    do k=1,Nvib
        Theta = Theta + Vec(j)*Vec(k)*Hess(j,k)
    enddo
    enddo
    Er_int = -Theta*0.5d0
    do i=1,Nvib
        Er_int = Er_int - Grad(i)*Vec(i)
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

    subroutine parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,intfile,rmzfile,def_internal,use_symmetry,vertical)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,&
                                          intfile,rmzfile,def_internal
        logical,intent(inout)          :: use_symmetry, vertical
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
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,'(/,A)') '        V E R T I C A L 2 A D I A B A T I C '    
        write(6,'(/,A)') '  Displace structure from vertical to adiabatic geoms  '
        write(6,'(/,A)') '           '        
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,*) '-f              ', trim(adjustl(inpfile))
        write(6,*) '-ft             ', trim(adjustl(ft))
        write(6,*) '-fhess          ', trim(adjustl(hessfile))
        write(6,*) '-fth            ', trim(adjustl(fth))
        write(6,*) '-fgrad          ', trim(adjustl(gradfile))
        write(6,*) '-ftg            ', trim(adjustl(ftg))
        write(6,*) '-intmode        ', trim(adjustl(def_internal))
        write(6,*) '-intfile        ', trim(adjustl(intfile))
        write(6,*) '-rmzfile        ', trim(adjustl(rmzfile))
        write(6,*) '-[no]vert      ',  vertical
        write(6,*) '-h             ',  need_help
        write(6,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program vertical2adiabatic

