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
               vertical=.true.      ,&
               do_correct_num=.false., &
               do_correct_int=.false., &
               gradcorrectS1=.false.
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
    real(8),dimension(NDIM,NDIM) :: B, G1
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv, L1,L2, Asel1, Asel2, L1int
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
    real(8),dimension(NDIM) :: Freq1, Freq2, S1, S2, Vec, Vec1, mu, Q0, FC
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !Shifts
    real(8),dimension(NDIM) :: Delta
    real(8) :: Delta_p, Er
    !====================== 

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: Hlt
    integer :: error
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
               I_DER=17,  &
               O_DUS=20,  &
               O_DIS=21,  &
               O_DMAT=22, &
               O_DUS2=23, &
               O_DIS2=24, &
               O_STAT=25, &
               O_STR =26
    !files
    character(len=10) :: ft ="guess", fth="guess", ftg="guess", ftgv="guess", fthv="guess" 
    character(len=200):: inpfile  ="input.fchk", &
                         hessfile ="same", &
                         gradfile="same",&
                         hessfile_v ="vertical.fchk", &
                         gradfile_v ="same", &
                         intfile  ="none",       &
                         rmzfile  ="none",       &
                         symm_file="none",     & 
                         derfile="base", derfile_base, &
                         tmpfile
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
    call parse_input(inpfile,ft,gradfile,ftg,hessfile,fth,gradfile_v,ftgv,hessfile_v,fthv,&
                     intfile,rmzfile,def_internal,use_symmetry,derfile,do_correct_num,do_correct_int,&
                     gradcorrectS1,vertical)
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
    if (ftgv == "guess") &
    call split_line_back(gradfile_v,".",null,ftgv)
    if (fthv == "guess") &
    call split_line_back(hessfile_v,".",null,fthv)

    ! STRUCTURE FILE
    print'(X,A)', "READING STATE1 FILE (STRUCTURE)..."
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,state1)
    close(I_INP)
    ! Shortcuts
    Nat = state1%natoms
    print'(X,A,/)', "Sucess!"

    ! HESSIAN FILE (State1)
    print'(X,A)', "READING STATE1 FILE (HESSIAN)..."
    open(I_INP,file=hessfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile)) )
    allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fth,Nat,Hlt,error) 
    close(I_INP)
    print'(X,A,/)', "Sucess!"
    ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
    print'(X,A)', "Preliminar vibrational analysis (Cartesian coordinates)..."
    call vibrations_Cart(Nat,state1%atom(:)%X,state1%atom(:)%Y,state1%atom(:)%Z,state1%atom(:)%Mass,Hlt,&
                         Nvib,L1,Freq1,error_flag=error)
    ! Get Lcart = M^1/2 Lmwc 
    call Lmwc_to_Lcart(Nat,Nvib,state1%atom(:)%Mass,L1,L1,error)
    if (verbose>2) &
        call MAT0(6,L1,3*Nat,Nvib,"Lcart matrix")

    ! GRADIENT FILE (State1)
    if (gradcorrectS1) then
        print'(X,A)', "READING STATE1 FILE (GRADIENT)..."
        open(I_INP,file=gradfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
        call generic_gradient_reader(I_INP,ftg,Nat,Grad,error) 
        close(I_INP)
        print'(X,A,/)', "Sucess!"
    else
        print'(X,A,/)', "Assuming gradient for State1 equal to zero"
        Grad(1:3*Nat) = 0.d0
    endif

    if (do_correct_int) then
        !***************************************************************
        ! Compute normal modes in internal coordinates (only needed in do_correct_int=.true.)
        print'(/,X,A)', "VIBRATIONAL ANALYSIS IN INTERNAL COORDINATES..."
    
        ! Reconstruct full hessian
        k=0
        do i=1,3*Nat
        do j=1,i
            k=k+1
            Hess(i,j) = Hlt(k)
            Hess(j,i) = Hlt(k)
        enddo 
        enddo
    
        !Generate bonded info
        !From now on, we'll use atomic units
        call guess_connect(state1)
        call gen_bonded(state1)
    
        ! Define internal set
        call define_internal_set(state1,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns)
    
        !From now on, we'll use atomic units
        call set_geom_units(state1,"Bohr")
    
        ! INTERNAL COORDINATES
    
        !SOLVE GF METHOD TO GET NM AND FREQ
        call internal_Wilson_new(state1,Ns,S1,B,ModeDef)
        call internal_Gmetric(Nat,Ns,state1%atom(:)%mass,B,G1)
        call calc_BDer(state1,Ns,Bder)
    
        ! SET REDUNDANT/SYMETRIZED/CUSTOM INTERNAL SETS
    !     if (symaddapt) then (implement in an analogous way as compared with the transformation from red to non-red
        if (Ns > Nvib) then ! Redundant
            call redundant2nonredundant(Ns,Nvib,G1,Asel1)
            ! Rotate Bmatrix
            B(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel1,B,tA=.true.)
            ! Rotate Gmatrix
            G1(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel1(1:Ns,1:Nvib),G1,counter=.true.)
            ! Rotate Bders
            if (vertical) then
                do j=1,3*Nat
                    Bder(1:Nvib,j,1:3*Nat) =  matrix_product(Nvib,3*Nat,Ns,Asel1,Bder(1:Ns,j,1:3*Nat),tA=.true.)
                enddo
            endif
        endif
    
        if (vertical) then
            call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B,G1,Grad=Grad,Bder=Bder)
        else
            call HessianCart2int(Nat,Nvib,Hess,state1%atom(:)%mass,B,G1)
        endif
        ! Do not overwrite Freq1 (although should be identical)
        call gf_method(Nvib,G1,Hess,L1int,Vec,X1,X1inv)
        if (verbose>1) then
            ! Analyze normal modes
            if (use_symmetry) then
                call analyze_internal(Nvib,Ns,L1int,Vec,ModeDef,S_sym)
            else
                call analyze_internal(Nvib,Ns,L1int,Vec,ModeDef)
            endif
        endif
        !***************************************************************
    endif
    deallocate(Hlt)

    ! HESSIAN FILE (State2)
    print'(/,X,A)', "READING STATE2 FILE (HESSIAN)..."
    open(I_INP,file=hessfile_v,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile_v)) )
    allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fthv,Nat,Hlt,error) 
    close(I_INP)
    k=0
    do i=1,3*Nat
    do j=1,i
        k=k+1
        Hess(i,j) = Hlt(k)
        Hess(j,i) = Hlt(k)
    enddo 
    enddo
    deallocate(Hlt)

    ! GRADIENT FILE
    print'(/,X,A)', "READING STATE2 FILE (GRADIENT)..."
    open(I_INP,file=gradfile_v,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile_v)) )
    call generic_gradient_reader(I_INP,ftgv,Nat,Grad,error)
    close(I_INP)

    !*****************************************************    
    ! Apply matrix derivative if the option is enabled
    !*****************************************************
    if (do_correct_int) then
        print'(X,A,/)', "Apply correction for vertical case based on internal vibrational analysis..."
        Aux(1:3*Nat,1:3*Nat) = Hess(1:3*Nat,1:3*Nat)
        ! Compute gQ
        ! Convert Gradient to normal mode coordinates in state1 Qspace.
        ! We use the internal normal modes and not the Cartesian
        ! ones to get the sign consistent with the internal mode
        ! definition. Note that, apart from the sign, both should be
        ! equivalent in state1 Qspace
        call HessianCart2int(Nat,Nvib,Aux,state1%atom(:)%mass,B,G1,Grad=Grad,Bder=Bder)
        do i=1,Nvib
            Vec(i) = 0.d0
            do k=1,Nvib !3*Nat
                Vec(i) = Vec(i) + L1int(k,i) * Grad(k)
            enddo
        enddo
        Grad(1:Nvib) = Vec(1:Nvib)

        ! Compute LLL^Q = Lint^-1 Lder
        ! And directly: gQ * LLL^Q
        L1int(1:Nvib,1:Nvib) = inverse_realgen(Nvib,L1int)
        do j=1,3*Nat 
        do k=1,3*Nat
            Aux(j,k) = 0.d0
            do i=1,Nvib
                Theta = 0.d0
                do l=1,Nvib
                    Theta = Theta + L1int(i,l) * Bder(l,j,k)
                enddo
                Aux(j,k) = Aux(j,k) + Grad(i) * Theta
            enddo        
        enddo
        enddo

        ! Compute (Hess - gQ LLL^Q)
        Hess(1:3*Nat,1:3*Nat) = Hess(1:3*Nat,1:3*Nat) - Aux(1:3*Nat,1:3*Nat)

        !Compute H_Q = L1^t (Hess - gQ LLL^Q) L1
        Hess(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,3*Nat,L1,Hess,counter=.true.)

    elseif (do_correct_num) then
        print'(X,A,/)', "Apply correction for vertical case based on numerical derivates (not-tested)..."
        ! Correct with numerical derivatives of Lcart matrix
        ! The derivatives are computed externally and fed through
        ! files. This is not working for the moment.
        !
        !Compute H_Q = L1^t Hess L1  +  gx LLL^x
        Hess(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,3*Nat,L1,Hess,counter=.true.)
        ! Fill Lder tensor
        derfile_base=derfile
        do j=1,Nvib
            write(derfile,'(A,I0,A)') trim(adjustl(derfile_base)), j, ".dat"
            open(I_DER,file=derfile,status='old',iostat=IOstatus)
            if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(derfile)) )
            do i=1,3*Nat
                ! Use the symbol Bder, but it is Lder!
                read(I_DER,*) Bder(i,j,1:Nvib)
            enddo
            close(I_DER)
        enddo

        if (verbose>2) then
            do i=1,3*Nat
                write(tmpfile,'(A,I0,A)') "Lder *10^6, Cart=",i
                call MAT0(6,Bder(i,:,:)*1.e6,Nvib,Nvib,trim(tmpfile))
            enddo
        endif

        do i=1,Nvib
        do j=1,Nvib
            Aux(i,j) = 0.d0
            do l=1,3*Nat
                Aux(i,j) = Aux(i,j) + Grad(l) * Bder(l,j,i) 
            enddo
            Hess(i,j) = Hess(i,j) + Aux(i,j)
        enddo
        enddo

        ! We need gQ...

    else
        print'(X,A,/)', "Do not apply any correction for vertical"
        ! Do not apply any correction (original PCCP2011 implementation)
        !Compute H_Q = L1^t Hess L1 
        Hess(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,3*Nat,L1,Hess,counter=.true.)

        ! We need gQ, in state1 normal modes
        ! We could use either internal or Cartesian coordinates to get it.
        ! Lets try Cartesia modes (we're using Lcart, not Lmwc)
        ! gQ = L^t * gx
        do i=1,Nvib 
            Vec(i) = 0.d0
            do k=1,3*Nat 
                Vec(i) = Vec(i) + L1(k,i) * Grad(k)
            enddo
        enddo
        Grad(1:Nvib) = Vec(1:Nvib)
    endif


    !-------------------
    ! Duschisky matrix
    !-------------------
    print'(/,X,A,/)', "DIAGONALIZE HESSIAN IN Q1-SPACE..."
    ! The matrix that diagonalizes the Hessian in Q1 modes is the Duschisky matrix
    call diagonalize_full(Hess(1:Nvib,1:Nvib),Nvib,G1(1:Nvib,1:Nvib),FC(1:Nvib),"lapack")

    !---------
    !Check FC
    !---------
    if (verbose>1) &
        call print_vector(6,FC*1.d6,Nvib,"FORCE CONSTANTS x 10^6 (A.U.)")
    !Transform FC to Freq
    do i=1,Nvib
        Freq2(i) = sign(dsqrt(abs(FC(i))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
                         FC(i))
        if (FC(i)<0) then
            print*, i, FC(i)
            call alert_msg("warning","A negative FC found")
        endif
    enddo
    if (verbose>0) &
        call print_vector(6,Freq2,Nvib,"Frequencies (cm-1)")

    !--------------------------------------------------------
    ! Shift vector
    !--------------------------------------------------------
    ! K = J^t * Q0 = J^t * [-FC^-1 * J * gQ]
    ! Where
    ! * J is the Duschisky matrix for the VH model
    !   STORED in G1
    ! * Q0   is the equilibrium geometry in terms of state2 Qspace
    !        note that the displacement is actually -Qo, but in state1 Qspace
    !   TO BE STORED in Q0 
    ! * FC are the diagonal force constant matrix for state2 (diagonal in the state2 Qspace)
    !   STORED in FC
    ! * gQ is the gradient of state2 in state1 Qspace 
    !   STORED in Grad
    !--------------------------------------------------------
    ! Q0 = - FC^-1 * J^t * gQ
    print'(X,A,/)', "COMPUTE SHIFT VECTOR..."
    do i=1,Nvib
        Q0(i) = 0.d0
        do k=1,Nvib
            Q0(i) = Q0(i) - G1(k,i) * Grad(k) / FC(i)
        enddo
    enddo
   ! K = J^t * Q0
    do i=1,Nvib
        Vec1(i) = 0.d0
        do k=1,Nvib
            Vec1(i) = Vec1(i) - G1(i,k) * Q0(k)
        enddo
    enddo
! This is simpler, but do not store Q0...
!     ! J^t * FC^-1 * J
!     Aux(1:Nvib,1:Nvib) = diag_basisrot(Nvib,Nvib,G1,1.d0/Freq2(1:Nvib),counter=.true.)
!     ! -[J^t * FC^-1 * J] * gQ
!     do i=1,Nvib
!         Vec1(i) = 0.d0
!         do k=1,Nvib
!             Vec1(i) = Vec1(i) - Aux(i,k) * Grad(k)
!         enddo
!     enddo

    !Analyze Duschinsky matrix
    call analyze_duschinsky(6,Nvib,G1,Vec1,Freq1,Freq2)


    !=======================
    ! REORGANIZATION ENERGY
    !=======================
    print*, "REORGANIZATION ENERGY"
    ! Normal-mode space
    ! Er = -gQ * Q0 - 1/2 * Q0^t * Lambda_f * Q0
    ! At this point: 
    ! * Grad: gradient in state1 normal modes (Qspace)
    ! * Q0  : DeltaQ in state2 Qspace
    ! * FC  : diagonal force constants for final state
    Er = 0.d0
    do i=1,Nvib
        ! Get gradient in state2 Qspace
        ! gQ2 = J^t * gQ1
        Theta=0.d0
        do k=1,Nvib
            Theta = Theta + G1(k,i) * Grad(k)
        enddo
        Er = Er - Theta * Q0(i) - 0.5d0 * FC(i) * Q0(i)**2
    enddo
    print'(X,A,F12.6)',   "Reorganization energy (AU) = ", Er
    print'(X,A,F12.6,/)', "Reorganization energy (eV) = ", Er*HtoeV

    !============================================
    ! PRINT DUSCHINSKI AND DISPLACEMENT TO FILES
    !============================================
    print*, "Printing Duschinski matrix to 'duschinsky.dat'"
    open(O_DUS, file="duschinsky.dat")
    print'(X,A,/)', "Printing Shift vector to 'displacement.dat'"
    open(O_DIS, file="displacement.dat")
    do i=1,Nvib
    do j=1,Nvib
        write(O_DUS,*)  G1(i,j)
    enddo 
        write(O_DIS,*)  Vec1(i)
    enddo
    close(O_DUS)
    close(O_DIS)

    !====================
    ! Print state files
    !====================
    ! State1
    call Lcart_to_LcartNrm(Nat,Nvib,L1,Aux,error)
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
        write(O_STAT,*) Aux(i,j)
    enddo
    enddo
    do j=1,Nvib
        write(O_STAT,'(F12.5)') Freq1(j)
    enddo
    close(O_STAT)
    ! State2
    ! L2 = L1 * J
    L2(1:3*Nat,1:Nvib) = matrix_product(3*Nat,Nvib,Nvib,L1,G1)
    call Lcart_to_LcartNrm(Nat,Nvib,L2,Aux,error)
    !Print state
    ! Note that the geometry is that of state1 (not displaced for vertical)
    ! But it is ok for FCclasses (it is not using it AFIK) What about HT??
    open(O_STAT,file="state_file_2")
    do i=1,Nat
        write(O_STAT,*) state1%atom(i)%x
        write(O_STAT,*) state1%atom(i)%y
        write(O_STAT,*) state1%atom(i)%z
    enddo
    do i=1,3*Nat
    do j=1,Nvib
        write(O_STAT,*) Aux(i,j)
    enddo
    enddo
    do j=1,Nvib
        write(O_STAT,'(F12.5)') Freq2(j)
    enddo
    close(O_STAT)



    call summary_alerts

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,gradfile,ftg,hessfile,fth,gradfile_v,ftgv,hessfile_v,fthv,intfile,&
                           rmzfile,def_internal,use_symmetry,derfile,do_correct_num,do_correct_int,&
                           gradcorrectS1,vertical)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,gradfile,ftg,hessfile,fth,gradfile_v,ftgv,hessfile_v,fthv,&
                                          intfile,rmzfile,def_internal,derfile
        logical,intent(inout)          :: use_symmetry,do_correct_num,do_correct_int,gradcorrectS1,vertical
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

                case ("-fder") 
                    call getarg(i+1, derfile)
                    argument_retrieved=.true.

                case ("-fhessv") 
                    call getarg(i+1, hessfile_v)
                    argument_retrieved=.true.
                case ("-fthv") 
                    call getarg(i+1, fthv)
                    argument_retrieved=.true.
                case ("-f2") 
                    call getarg(i+1, hessfile_v)
                    argument_retrieved=.true.
                case ("-ft2") 
                    call getarg(i+1, fthv)
                    argument_retrieved=.true.

                case ("-fgradv") 
                    call getarg(i+1, gradfile_v)
                    argument_retrieved=.true.
                case ("-ftgv") 
                    call getarg(i+1, ftgv)
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

                case("-vert")
                    vertical=.true.
                case("-novert")
                    vertical=.false.

                case ("-sym")
                    use_symmetry=.true.
                case ("-nosym")
                    use_symmetry=.false.

                case ("-correct-num")
                    do_correct_num=.true.
                case ("-nocorrect-num")
                    do_correct_num=.false.
        
                case ("-correct-int")
                    do_correct_int=.true.
                case ("-nocorrect-int")
                    do_correct_int=.false.
                    gradcorrectS1=.false.

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
           gradfile=hessfile
           if (adjustl(ftg) == "guess")  ftg=fth
       endif
       if (adjustl(gradfile_v) == "same") then
           gradfile_v=hessfile_v
           if (adjustl(ftgv) == "guess")  ftgv=fthv
       endif


       !Print options (to stderr)
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,'(/,A)') '      C A R T E S I A N  D U S C H I N S K Y        '
        write(6,'(/,A)') '         Duschinski analysis for Vertical         '
        write(6,'(A,/)') '          model in Cartesian coordinates          '        
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,*) '-f                Input file State1(S1)      ', trim(adjustl(inpfile))
        write(6,*) '-ft               \_FileType                 ', trim(adjustl(ft))
        write(6,*) '-fhess            S1 Hessian file name       ', trim(adjustl(hessfile))
        write(6,*) '-fth              \_FileType                 ', trim(adjustl(fth))
        write(6,*) '-fgrad            S1 gradient file name      ', trim(adjustl(gradfile))
        write(6,*) '-ftg              \_FileType                 ', trim(adjustl(ftg))
        write(6,*) '-fhessv           S2(vertical) Hessian file  ', trim(adjustl(hessfile_v))
        write(6,*) '                  (-f2 is a synonim)         '
        write(6,*) '-fthv             \_FileType                 ', trim(adjustl(fthv))
        write(6,*) '                  (-ft2 is a synonim)        '
        write(6,*) '-fgradv           S2(vertical) gradient file ', trim(adjustl(gradfile_v))
        write(6,*) '-ftgv             \_FileType                 ', trim(adjustl(ftgv))
        write(6,*) '-vert             Vertical model              ', vertical
        write(6,*) ''
        write(6,*) '** Options correction method (vertical) **'
        write(6,*) '-[no]correct-num  Correction with numerical   ', do_correct_num
        write(6,*) '                  derivatives of L1 (Cart)    '
        write(6,*) '-[no]correct-int  Correction with analytical  ', do_correct_int
        write(6,*) '                  L1 ders based on internal   '
        write(6,*) '                  analysis                    '
        write(6,*) '-fder             Numerical derivative file  ', trim(adjustl(derfile))
        write(6,*) '                  basename (-correct-num)    '
        write(6,*) '-intmode          Internal set [zmat|sel|all]', trim(adjustl(def_internal))
        write(6,*) '                  (-correct-int)'
        write(6,*) '-intfile          File with internal set def.', trim(adjustl(intfile))
        write(6,*) '                  (-correct-int -intmode sel)'
        write(6,*) '-[no]corrS1       Correct S1 at vib-in       ', gradcorrectS1
!         write(6,*) '-rmzfile        ', trim(adjustl(rmzfile))
        write(6,*) ''
        write(6,*) '-h               ',  need_help
        write(6,*) '--------------------------------------------------'
        if (.not.vertical) call alert_msg("fatal", 'Adiabatic model not yet implemented' )
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program vertical2adiabatic

