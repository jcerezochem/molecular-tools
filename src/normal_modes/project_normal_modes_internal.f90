program projection_normal_modes_int


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
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
    logical :: use_symmetry=.false.,   &
               include_hbonds=.false., &
               vertical=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule, molec_aux
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns, Nf
    character(len=5) :: PG
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=50) :: dummy_char
    real(8) :: dist
    !io flags
    integer :: error, info
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, jj
    !=============

    !====================== 
    ! PES topology and normal mode things
    real(8),dimension(:),allocatable :: A
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, LL
    real(8),dimension(NDIM) :: Freq, Factor, Grad
    !Coord
    real(8) :: qcoord
    !Delta
    real(8),dimension(1:NDIM) :: Delta
    real(8)                   :: Delta_p
    ! Integral
    real(8),dimension(1:NDIM) :: FC, Q0, Q0b
    real(8) :: f0, f1, area, t, dt, ff
    integer :: Nvib0
    !====================== 

    !======================
    ! Auxiliar
    real(8),dimension(NDIM,NDIM) :: Aux2
    !======================

    !====================== 
    !INTERNAL CODE THINGS
    real(8),dimension(1:NDIM,1:NDIM) :: B,G
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    real(8),dimension(1:NDIM,1:NDIM) :: X,Xinv
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: S1,S2
    integer,dimension(NDIM) :: S_sym
    ! Switches
    character(len=5) :: def_internal="ZMAT"
    !Coordinate map
    integer,dimension(NDIM) :: Zmap
    !====================== 

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_STR=11,  &
               I_SYM=12,  &
               I_RMF=16,  &
               O_DIS=20

    !files
    character(len=10) :: ft ="guess",  ftg="guess",  fth="guess", fts="guess"
    character(len=200):: inpfile  ="equil.fchk", &
                         gradfile ="same", &
                         hessfile ="same", &
                         intfile  ="none", &
                         rmzfile  ="none", &
                         symm_file="none", &
                         strfile  ="structure.xyz", &
                         dispfile ="disp.dat"
    !status
    integer :: IOstatus
    !===================

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

    call cpu_time(ti)

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(molecule)
    call allocate_atoms(molec_aux)
    !===========================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(&
                     ! input data
                     inpfile,ft,hessfile,fth,gradfile,ftg,strfile,fts,  &
                     ! Options (general)
                     include_hbonds,vertical,                           &
                     ! Options (internal)
                     use_symmetry,def_internal,intfile,rmzfile)


    ! INTERNAL VIBRATIONAL ANALYSIS
 
    ! 1. READ DATA
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
    call generic_strmol_reader(I_INP,ft,molecule,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (State1)")
    close(I_INP)
    ! Shortcuts
    Nat = molecule%natoms

    ! HESSIAN FILE
    open(I_INP,file=hessfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile)) )
    allocate(A(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fth,Nat,A,error) 
    if (error /= 0) call alert_msg("fatal","Error reading Hessian (State1)")
    close(I_INP)
    ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
    call vibrations_Cart(Nat,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,molecule%atom(:)%Mass,A,&
                         Nvib,LL,Freq,error_flag=error)
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
    if (vertical) then
        open(I_INP,file=gradfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
        call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
        close(I_INP)
    endif


    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! Get connectivity 
    call guess_connect(molecule)

    ! Manage symmetry
    if (.not.use_symmetry) then
        molecule%PG="C1"
    else if (trim(adjustl(symm_file)) /= "NONE") then
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

    ! Define internal set
    call define_internal_set(molecule,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns,Nf,Aux2)
    if (Ns > Nvib) then
        call alert_msg("fatal","Non-redundan coordinate set needs mapping (still on dev)")
        ! Need mapping from whole set to Zmat
    elseif (Ns > Nvib) then
        call alert_msg("fatal","Reduced coordinates cases still not implemented")
        ! Need to freeze unused coords to its input values
    endif

    !From now on, we'll use atomic units
    call set_geom_units(molecule,"Bohr")

    !SOLVE GF METHOD TO GET NM AND FREQ
    call internal_Wilson(molecule,Nvib,S1,B,ModeDef)
    call internal_Gmetric(Nat,Nvib,molecule%atom(:)%mass,B,G)
    if (vertical) then
        call calc_Bder(molecule,Nvib,Bder)
        ! (Hess is already constructed)
        ! Hs (with the correction)
        ! First get: Hx' = Hx - gs^t\beta
        ! 1. Get gs from gx
        call Gradcart2int(Nat,Nvib,Grad,molecule%atom(:)%mass,B,G)
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
    call HessianCart2int(Nat,Nvib,Hess,molecule%atom(:)%mass,B,G)
    call gf_method(Nvib,Nvib,G,Hess,LL,Freq,X,Xinv)
    if (verbose>0) then
        ! Analyze normal modes
        if (use_symmetry) then
            call analyze_internal(Nvib,Ns,LL,Freq,ModeDef,S_sym)
        else
            call analyze_internal(Nvib,Ns,LL,Freq,ModeDef)
        endif
    endif
    !Define the Factor to convert shift into addimensional displacements
    ! from the shift in SI units:
    Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib))*1.d2*clight*2.d0*PI/plankbar)
    ! but we need it from au not SI
    Factor(1:Nvib)=Factor(1:Nvib)*BOHRtoM*dsqrt(AUtoKG)


    !==========================================================0
    !  Projection on Normal modes
    !==========================================================0
    ! 1. Structure to be analyzed
    if (fts == "guess") &
    call split_line_back(strfile,".",null,fts)
    open(I_STR,file=strfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(strfile)) )
    call generic_strmol_reader(I_STR,fts,molec_aux,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (Strfile)")
    close(I_STR)

    ! Set to A.U.
    call set_geom_units(molec_aux,"Bohr")

    ! 2. Get internal coordinates
    molec_aux%geom = molecule%geom
    call internal_Wilson(molec_aux,Nvib,S2,B,ModeDef)

    ! 3. Compute DeltaS
    if (verbose>0) then
        print*, ""
        print*, "=========================="
        print*, " SHIFTS (internal coord)"
        print*, "=========================="
    endif
    if (verbose>0) &
     print*, "Bonds"
    do i=1,molecule%geom%nbonds
        Delta(i) = S2(i)-S1(i)
        if (verbose>0) &
         print'(I5,3(F8.2,2X))', i, S2(i),S1(i), Delta(i)
    enddo
    if (verbose>0) &
     print*, "Angles"
    do j=i,i+molecule%geom%nangles-1
        Delta(j) = S2(j)-S1(j)
        if (verbose>0) &
         print'(I5,3(F8.2,2X))', j, S2(j)*180.d0/PI,S1(j)*180.d0/PI,Delta(j)*180.d0/PI
    enddo
    if (verbose>0) &
     print*, "Dihedrals"
    do k=j,j+molecule%geom%ndihed-1
        Delta(k) = S2(k)-S1(k)
        Delta_p = S2(k)-S1(k)+2.d0*PI
        if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
        Delta_p = S2(k)-S1(k)-2.d0*PI
        if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
        if (verbose>0) &
         print'(I5,3(F8.2,2X))', k, S2(k)*180.d0/PI,S1(k)*180.d0/PI,Delta(k)*180.d0/PI
    enddo

    ! 4. Compute DeltaQ
    LL(1:Nvib,1:Nvib) = inverse_realgen(Nvib,LL)
    open(O_DIS,file=dispfile)
    ! dispfile: Column 1: AU units
    !           Column 2: dimesionless
    do j=1,Nvib
        qcoord = 0.d0
        do i=1,Nvib
            qcoord = qcoord + LL(j,i) * Delta(i)
        enddo

        write(O_DIS,'(X,I6,X,2(F12.6,2X))') j,      &
                                            qcoord, &        ! AU
                                            qcoord*Factor(j) ! Dimless
        Q0(j) = qcoord
    enddo
    close(O_DIS)

    ! Euclidean distance
    print*, ""
    dist=0.d0
    do i=1,Nvib
        dist = dist + Q0(i)**2
    enddo
    dist=dsqrt(dist)
    print'(X,A,F10.4)', "Euclidean distance in nm space", dist/dsqrt(AMUtoAU)

    ! RC path distance
    FC(1:Nvib) = Freq2FC(Nvib,Freq)
    dist=0.d0
    area=1.d0
    dt=5.d1
    t = 0.d0
    Nvib0=Nvib
    do while (dabs(area) > 1d-10 .and. Nvib0>0)
        f0=0.d0
        f1=0.d0
        Nvib0=Nvib
        j = 0
        do i=1,Nvib
            ff = FC(i)**2*Q0(i)**2*dexp(-2.d0*FC(i)*t)
            f0 = f0 + ff
            f1 = f1 + FC(i)**2*Q0(i)**2*dexp(-2.d0*FC(i)*(t+dt))
            ! Discard modes that reached the baseline
            if (ff < 5e-24) then
                Nvib0=Nvib0-1
            else
                j = j+1
                Q0b(j) = Q0(i)
            endif
        enddo
        Nvib = Nvib0
        Q0(1:Nvib) = Q0b(1:Nvib)
        f0 = dsqrt(f0)
        f1 = dsqrt(f1)
        area = 0.5d0*(f0+f1)*dt
        dist = dist + area
        t=t+dt
    enddo
    print'(X,A,F10.4)', "Contour distance in IRC space ", dist/dsqrt(AMUtoAU)

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(&
                           ! input data
                           inpfile,ft,hessfile,fth,gradfile,ftg,strfile,fts,  &
                           ! Options (general)
                           include_hbonds,vertical,                           &
                           ! Options (internal)
                           use_symmetry,def_internal,intfile,rmzfile)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,strfile,fts, &
                                          intfile,rmzfile,def_internal
        logical,intent(inout)          :: include_hbonds,vertical, use_symmetry

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

                case ("-fs") 
                    call getarg(i+1, strfile)
                    argument_retrieved=.true.
                case ("-fts") 
                    call getarg(i+1, fts)
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

                case ("-sym")
                    use_symmetry=.true.
                case ("-nosym")
                    use_symmetry=.false.

                case ("-vert")
                    vertical=.true.
                case ("-novert")
                    vertical=.false.

                case ("-include_hb")
                    include_hbonds=.true.
        
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
           fth=ft
       endif
       if (adjustl(gradfile) == "same") then
           gradfile=inpfile
           ftg=ft
       endif


       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '          PROJECT NORMAL MODES INTERNAL '           
        write(0,'(/,A)') '                                        '
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(ft))
        write(0,*) '-fhess          ', trim(adjustl(hessfile))
        write(0,*) '-fth            ', trim(adjustl(fth))
        write(0,*) '-fgrad          ', trim(adjustl(gradfile))
        write(0,*) '-ftg            ', trim(adjustl(ftg))
        write(0,*) '-fs             ', trim(adjustl(strfile))
        write(0,*) '-fts            ', trim(adjustl(fts))
        write(0,*) '-intmode        ', trim(adjustl(def_internal))
        write(0,*) '-intfile        ', trim(adjustl(intfile))
        write(0,*) '-rmzfile        ', trim(adjustl(rmzfile))
        write(0,*) '-[no]sym       ',  use_symmetry
        write(0,*) '-[no]vert      ',  vertical
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input




end program projection_normal_modes_int

