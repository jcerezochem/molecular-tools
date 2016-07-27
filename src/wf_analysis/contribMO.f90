program contribMO

    !==============================================================
    ! This code uses MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    ! Description
    !  Utility to compute the contribution of atom (label=1) in the MO
    ! Compilation intructions:
    !  gfortran ../modules/constants_mod.f90 ../modules/alerts.f90 ../modules/structure_types_v2_ALLOC.f90 ../modules/allocation_mod.f90 ../modules/line_preprocess.f90 ../modules/gaussian_manage_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/atomic_geom_v2.f90 ../modules/MatrixMod.f90 contribMO_v4.f90 -o contribMO_v4.exe -cpp -llapack
    !
    ! History
    ! v2 (20/07/2012):
    !       -at changed to -frg: the file may contain now several fragments
    !       NTO specific feature: contributions for all NTO and average.
    ! v3 (28/06/2013)
    !       Enlarge memory
    ! v4 (6/11/2013)
    !       Disable NTO averaging. This is not conceptually clear. If there is not one predominant NTO,
    !        NTO analysis is just meaningless (or is it?)
    !        Now it gives a warning if NTO coef is < 0.9
    !***************************
    ! v4.4 (Adapt to distributed modules)
    !       Change real(4) to real(8)
    !
    ! The calculation is done following the proceduce in GaussSum
    ! File: gausssum/popanalysis.py 
    ! Line 495
    ! Code (python):
    !     contrib = [x * numpy.dot(x,overlap) for x in MOCoeff]
    ! where contrib is the contains the MO**2 "corrected"
    !==============================================================

    use matrix
    use matrix_print
    use gaussian_manage

    implicit none

    integer,parameter :: BASIS_SIZE = 2000

    real(8),dimension(:,:),allocatable :: MO, S, MOc
    integer,dimension(:),allocatable :: AOmap
    
    integer :: IFCHK =11, &
               IMO   =12, &
               IAT   =13
    character(len=50) :: fchkfile,      &
                         mo_file="all", &
                         at_file="ndx"

    integer :: N,M, i, j, k, ii, frg, Nb, No

    real(8) :: coeff_total, coeff_metal

    integer :: Metal_Index = 1
    character(len=5) :: mo_type
    character(len=1) :: dtype
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA

    !OM selection stuff
    integer :: Nselect, Nselect_alpha, Nselect_beta, Nfrg
    integer,dimension(1:10) :: Nat
    integer,dimension(1:10,1:500) :: atom_index
    integer,dimension(:),allocatable :: MO_select_alpha, MO_select_beta

    !NTO averaging
    integer :: Nel
    real(8) :: Norm
    real(8),dimension(1:10) :: hole, particle
    real(8),dimension(1:BASIS_SIZE) :: W, W2
    logical :: nto_av = .false., do_nto=.false.

    call parse_input(fchkfile,mo_file,at_file,Metal_Index,do_nto)

    open(IFCHK,file=fchkfile,status="old")

!     call read_MO(IFCHK,MO,N,"alpha",-1)
    call read_fchk(IFCHK,"Number of basis functions",dtype,N,A,IA,j)
    rewind(IFCHK)
    Nb = IA(1)
    deallocate(IA)
    call read_fchk(IFCHK,"Number of independent functions",dtype,N,A,IA,j)
    rewind(IFCHK)
    No = IA(1)
    deallocate(IA)
    call read_fchk(IFCHK,"Alpha MO coefficients",dtype,N,A,IA,j)
    rewind(IFCHK)
    allocate(MO(1:Nb,1:No),S(1:Nb,1:Nb),MOc(1:Nb,1:No))
    k=0
    do i=1,No
    do j=1,Nb
        k=k+1
        MO(j,i) = A(k)
    enddo
    enddo
    deallocate(A)
    call overlap_matrix(Nb,No,MO,S)
    call MAT0(99,S,Nb,Nb,"S")
    call MAT0(98,MO,Nb,No,"C")
    !Compute overlap corrected MOs 
    MOc(1:Nb,1:No) = matrix_product(Nb,No,Nb,S,MO)
    deallocate(S)
    allocate(AOmap(1:Nb))
    call AO_atom_map(IFCHK,AOmap(1:Nb))
    
    !Atom selection
    if (adjustl(at_file) == "ndx") then
        Nfrg   = 1
        Nat(1) = 1
        atom_index(1,1) = Metal_Index
    else
        open(IAT,file=at_file,status="old")
        read(IAT,*) Nfrg
        do i=1,Nfrg
            read(IAT,*) Nat(i)
            read(IAT,*) atom_index(i,1:Nat(i))
        enddo
    endif
    
    allocate(MO_select_alpha(1:No))
    !Enable OM selection (unless do_nto is enabled -- in this case only the first one is analyzed)
    if (adjustl(mo_file) == "all") then ! .or. nto_av) then
        Nselect_alpha = No
        Nselect_beta  = No
        MO_select_alpha(1:No) = (/(i,i=1,No)/)
!         MO_select_beta(1:No)  = (/(i,i=1,No)/)
    else if (adjustl(mo_file) == "alpha") then
        Nselect_alpha = No
        Nselect_beta  = 0
        MO_select_alpha(1:No) = (/(i,i=1,No)/)
    else
        open(IMO,file=mo_file,status="old")
        read(IMO,*) Nselect
        Nselect_alpha = 0
        Nselect_beta  = 0
        do i=1,Nselect
            read(IMO,*) mo_type, ii
            if (adjustl(mo_type) == "alpha") then 
                Nselect_alpha = Nselect_alpha + 1
                MO_select_alpha(Nselect_alpha) = ii
            else if (adjustl(mo_type) == "beta") then
                Nselect_beta  = Nselect_beta + 1
!                 MO_select_beta(Nselect_beta)   = ii
            endif
        enddo
    endif

    !Alpha
    !=====
    !NTO averaging: get NTO coefficients (W) (stored in MO energies section)
    ! note that after this call M=N --- disabled
    if (do_nto) then
        print*, "NTO is disable for the moment"
        stop
        write(0,*) "=============="
        write(0,*) " ANALISIS NTO"
        write(0,*) "=============="
        write(0,*) "NOTE: NTO averaging has been disabled, the analysis will be done only"
        write(0,*) "      on the first NTO orbital. This is conceptually not clear. Think" 
        write(0,*) "      about it and report back, thanks."

        call read_NrEl(IFCHK,Nel,"alpha")
        call read_MOenergy(IFCHK,W,No,"alpha",-1)
!         !Normalize
        Norm = 0.0
        do i=1,Nel
            Norm = Norm + W(i)
        enddo
        do i=1,Nel
            W2(i) = W(i)/Norm
            !The particle is equal
            j = 2*Nel + 1 - i
            W2(j) = W2(j)/Norm
        enddo
! 
!         !Set the first NTO orbital
!         Nselect_alpha = 2
!         MO_select_alpha(1) = Nel
!         MO_select_alpha(2) = Nel+1
    endif


    do frg = 1, Nfrg
        hole(frg) = 0.0
        particle(frg) = 0.0

        do ii=1,Nselect_alpha
            i = MO_select_alpha(ii)
            ! Compute contributions
            coeff_total = 0.0
            coeff_metal = 0.0
            do j=1,Nb
                coeff_total = coeff_total + MO(j,i)*MOc(j,i)
                !Contribution of fragments
                do k=1,Nat(frg)
                    if ( AOmap(j) == atom_index(frg,k) ) &
                     coeff_metal = coeff_metal + MO(j,i)*MOc(j,i)
                enddo
            enddo
            if (.not.do_nto) &
            print'(X,A,I3,A,I3,A,I3,A,I4,A, F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to MO alpha", i,":",coeff_metal/coeff_total * 100
            !NTO av: store in array
            if (do_nto) then
                if (i == Nel) then
                    hole(frg)     = hole(frg) + coeff_metal/coeff_total * 100.0 !* W(i) -- averaging disabled
                else if (i == Nel+1) then
                    particle(frg) = particle(frg) + coeff_metal/coeff_total * 100.0 !* W(i) -- averaging disabled
                endif
            endif

        enddo

    enddo

    !NTO:
    if (do_nto) then
        print*, "Sum of Alpha Coefficients: ", Norm
        print*, "First NTO coef           : ", W(Nel)!, W(Nel+1) 
        do frg = 1, Nfrg
            print'(X,A,I3,A,I3,A,I3,A,F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to alpha HOLE    :",hole(frg)
            print'(X,A,I3,A,I3,A,I3,A,F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to alpha PARTICLE:",particle(frg)
        enddo
    endif

    
    !If no beta requested, stop here!
    deallocate(MO,MOc,MO_select_alpha)
    if ( Nselect_beta == 0 ) stop
    print*, "Unrestricted calculation not yet supported"
    stop
!    print*, ""
    
    !Beta
    !=====
    call read_MO(IFCHK,MO,N,"beta",-1)
    !NTO averaging: get NTO coefficients (W) (stored in MO energies section)
    ! note that afeter this call M=N
    if (do_nto) then
        call read_NrEl(IFCHK,Nel,"beta")
        call read_MOenergy(IFCHK,W,N,"beta",-1)
        !Normalize
        Norm = 0.0
        do i=1,Nel
            Norm = Norm + W(i)
        enddo
        do i=1,Nel
            W2(i) = W2(i)/Norm
            !The particle is equal
            j = 2*Nel + 1 - i
            W2(j) = W2(j)/Norm
        enddo
    endif


    do frg = 1, Nfrg
        hole(frg) = 0.0
        particle(frg) = 0.0

        do ii=1,Nselect_beta
            i = MO_select_beta(ii)
            ! Compute contributions
            coeff_total = 0.0
            coeff_metal = 0.0
            do j=1,M
                coeff_total = coeff_total + abs(MO(j,i))
                !Contribution of fragments
                do k=1,Nat(frg)
                    if ( AOmap(j) == atom_index(frg,k) ) &
                     coeff_metal = coeff_metal + abs(MO(j,i))
                enddo
            enddo
            if (.not.do_nto) &
            print'(X,A,I3,A,I3,A,I3,A,I4,A, F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to MO beta ", i,":",coeff_metal/coeff_total * 100
            !NTO av: store in array
            if (do_nto) then
                if (i == Nel) then
                    hole(frg)     = hole(frg) + coeff_metal/coeff_total * 100.0 !* W(i)
                else if (i == Nel+1) then
                    particle(frg) = particle(frg) + coeff_metal/coeff_total * 100.0 !* W(i)
                endif
            endif

        enddo

    enddo

    !NTO av:
    if (do_nto) then
        print*, "Sum of Beta Coefficients: ", Norm
        print*, "First NTO coef           : ", W(Nel)!, W(Nel+1) 
        do frg = 1, Nfrg
            print'(X,A,I3,A,I3,A,I3,A,F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to beta  HOLE    :",hole(frg)
            print'(X,A,I3,A,I3,A,I3,A,F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to beta  PARTICLE:",particle(frg)
        enddo
    endif
     
    stop

    contains

    subroutine parse_input(inpfile,mo_file,at_file,Metal_Index,do_nto)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile, mo_file, at_file
        integer,intent(inout) :: Metal_Index
        logical,intent(inout) :: do_nto
        
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg, dummy_char

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

                case ("-mo")
                    call getarg(i+1, mo_file)
                    argument_retrieved=.true.

                case ("-frg")
                    call getarg(i+1, at_file)
                    argument_retrieved=.true.

                case ("-ndx")
                    call getarg(i+1, dummy_char)
                    read(dummy_char,*) Metal_Index
                    argument_retrieved=.true.

                case ("-do_nto")
!                     nto_av = .true.
                    do_nto = .true.

                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo


       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '              c o n t r i b  M O '
        write(0,'(/,A)') '     Compute the contribution of an atom '
        write(0,'(A)')   '            to a given set of MOs'
        write(0,'(/,A)') '                    V4'
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-mo             ', trim(adjustl(mo_file))
        write(0,*) '-frg            ', trim(adjustl(at_file))
        write(dummy_char,*) Metal_Index
        write(0,*) '-ndx            ', trim(adjustl(dummy_char))
        write(0,*) '-do_nto        ',  do_nto
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


    subroutine AO_atom_map(unt,AOmap)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the mapping between atomic orbitals and atoms from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! AOmap (int,dim(:);out): mapping array
        !==============================================================

        integer,dimension(:),intent(out) :: AOmap

        integer,dimension(5000) :: ShellType,AO_shells_map

        !Reading stuff
        character(len=50) :: header_of_section
        character(len=240) :: line=""

        !Auxiliar variables and Dummies
!         character(len=30) :: dummy_char
!         integer :: dummy_int

        !=============
        !Counters
        integer :: i,k
        integer :: Nshells, Shell_size
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================



        !Read shell types
        header_of_section="Shell types"
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching shell types")
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') Nshells
                    exit
                endif
        enddo

        read(unt,*) ShellType(1:Nshells)

        !Read mapping
        header_of_section="Shell to atom map"
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching AO mappings")
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') Nshells
                    exit
                endif
        enddo

        read(unt,*) AO_shells_map(1:Nshells)

        k=0
        do i=1,Nshells
            !Translate shell type to shell size
            select case (ShellType(i))
             case (0)
                 ! S
                 Shell_size = 1
             case (1)
                 ! P
                 Shell_size = 3
             case (-1)
                 ! SP
                 Shell_size = 4
             case (2)
                 ! D (cartesian)
                 Shell_size = 6
             case (-2)
                 ! D (pure)
                 Shell_size = 5
             case (3)
                 ! F (cartesian)
                 Shell_size = 10
             case (-3)
                 ! F (pure)
                 Shell_size = 7
             case default
                 Shell_size = 2*iabs(ShellType(i)) + 1
            end select

            AOmap(k+1:k+Shell_size) = AO_shells_map(i)
            k = k + Shell_size

        enddo

        rewind(unt)

        return   

    end subroutine AO_atom_map


    subroutine read_MO(unt,MO,N,dtype,isel)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read density matrix (in a.o. basis) from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! MO(real,dim(:,:);out): array which columns correspond to MO coefficients
        !                        in rows: MO(1:M,i) ar MO_i coefficients
        ! N(int;out): number of elements read
        ! dtype(char;out) : type of orbitals (alpha/beta)
        ! isel(char;in): target a given MO. If -1, retrieve all
        !==============================================================


        !Arguments
#ifdef DOUBLE
        double precision, dimension(:,:), intent(out) :: MO
#else
        real, dimension(:,:), intent(out) :: MO
#endif
        integer, intent(in)  :: isel
        integer, intent(out) :: N
        character(len=*), intent(in) :: dtype

        !Local auxiliar array (will store ALL coefficients in the SR)
#ifdef DOUBLE
        double precision, dimension(:,:),allocatable :: MO_L
#else
        real, dimension(:,:),allocatable :: MO_L
#endif
        integer :: M

        !Reading stuff
        character(len=50)  :: header_of_section
        character(len=240) :: line=""
        character(len=1)   :: data_type

        !Auxiliar variables and Dummies
!         character(len=30) :: dummy_char
!         integer :: dummy_int

        !=============
        !Counters
        integer :: i,j
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================


        !Set variables
        select case ( adjustl(dtype) )
            case ( "alpha" )
             header_of_section="Alpha MO coefficients"
            case ( "beta" )
             header_of_section="Beta MO coefficients"
            case default
             call alert_msg( "fatal","Incorrect usage of read_MO subroutine. dtype cannot be: "//trim(adjustl(dtype)) )
        end select

        ! Search MO type
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching density: "//trim(adjustl(dtype)))
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') N
                    exit
                endif
        enddo

        !Calculate the number of states as the sqrt of N
        M = int ( sqrt(float(N)) )

        !Allocate MO_L
        allocate(MO_L(1:M,1:M))

        ! Read density
         read(unt,*) MO_L(1:M,1:M)

        if (isel == -1) then
            MO(1:M,1:M) = MO_L(1:M,1:M)
        else
            MO(1,1:M)   = MO_L(isel,1:M)
        endif
 

        rewind(unt)
        return

    end subroutine read_MO

    subroutine read_MOenergy(unt,E_MO,N,dtype,isel)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read density matrix (in a.o. basis) from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! E_MO(real,dim(:);out): vector with MO energies
        ! N(int;out): number of elements read
        ! dtype(char;out) : type of orbitals (alpha/beta)
        ! isel(char;in): target a given MO. If -1, retrieve all
        !==============================================================


        !Arguments
#ifdef DOUBLE
        double precision, dimension(:), intent(out) :: E_MO
#else
        real, dimension(:), intent(out) :: E_MO
#endif
        integer, intent(in)  :: isel
        integer, intent(out) :: N
        character(len=*), intent(in) :: dtype

        !Local auxiliar array (will store ALL coefficients in the SR)
#ifdef DOUBLE
        double precision, dimension(:),allocatable :: E_MO_L
#else
        real, dimension(:),allocatable :: E_MO_L
#endif
        integer :: M

        !Reading stuff
        character(len=50) :: header_of_section
        character(len=240) :: line=""

        !Auxiliar variables and Dummies
!         character(len=30) :: dummy_char
!         integer :: dummy_int

        !=============
        !Counters
        integer :: i,j
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================


        !Set variables
        select case ( adjustl(dtype) )
            case ( "alpha" )
             header_of_section="Alpha Orbital Energies"
            case ( "beta" )
             header_of_section="Beta Orbital Energies"
            case default
             call alert_msg( "fatal","Incorrect usage of read_MOenergy subroutine. dtype cannot be: "//trim(adjustl(dtype)) )
        end select

        ! Search MO type
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching density: "//trim(adjustl(dtype)))
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') N
                    exit
                endif
        enddo

        !Allocate E_MO_L
        allocate(E_MO_L(1:N))

        ! Read density
         read(unt,*) E_MO_L(1:N)

        if (isel == -1) then
            E_MO(1:N) = E_MO_L(1:N)
        else
            E_MO(1)   = E_MO_L(isel)
        endif
 

        rewind(unt)
        return

    end subroutine read_MOenergy


    subroutine read_NrEl(unt,N,dtype)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Get nr of electrons from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! Nel(int;out): number of electrons
        ! dtype(char;out) : type of orbitals (alpha/beta)
        !==============================================================


        !Arguments
        integer, intent(out) :: N
        character(len=*), intent(in) :: dtype

        !Reading stuff
        character(len=50) :: header_of_section
        character(len=240) :: line=""

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================


        !Set variables
        select case ( adjustl(dtype) )
            case ( "total" )
             header_of_section="Number of electrons"
            case ( "alpha" )
             header_of_section="Number of alpha electrons"
            case ( "beta" )
             header_of_section="Number of beta electrons"
            case default
             call alert_msg( "fatal","Incorrect usage of read_NrEl subroutine. dtype cannot be: "//trim(adjustl(dtype)) )
        end select

        ! Search MO type
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while getting NrEl: "//trim(adjustl(dtype)))
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') N
                    exit
                endif
        enddo

        rewind(unt)
        return

    end subroutine read_NrEl

    subroutine overlap_matrix(Nb,No,MO,S)

        integer,intent(in)                 :: Nb, No
        real(8),dimension(:,:),intent(in)  :: MO
        real(8),dimension(:,:),intent(out) :: S

        ! S = (C^t)^-1 C^-1 = (C^-1)^t C^-1
        ! In the case of non-independent basis, Norb<Nbasis
        ! we use the generalized inverse
        ! S = C (C^tC)^-1 (C^tC)^-1 C^t
        S(1:No,1:No) = matrix_product(No,No,Nb,MO,MO,tA=.true.)
        S(1:No,1:No) = inverse_realgen(No,S)
        S(1:No,1:No) = matrix_product(No,No,No,S,S)
        S(1:Nb,1:No) = matrix_product(Nb,No,No,MO,S)
        S(1:Nb,1:Nb) = matrix_product(Nb,Nb,No,S,MO,tB=.true.)
        
        return

    end subroutine overlap_matrix


end program contribMO

    

