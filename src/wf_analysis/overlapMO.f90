program printMO

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
    use verbosity

    implicit none

    integer,parameter :: BASIS_SIZE = 2000

    real(8),dimension(:,:),allocatable :: MO1,MO2, S1,S2, MOc
    integer,dimension(:),allocatable :: AOmap
    
    integer :: IFCHK =11, &
               IMO   =12, &
               IAT   =13
    character(len=50) :: fchkfile1="file1.fchk",&
                         fchkfile2="none",      &
                         mo_sel="all"

    integer :: N,M, i, j, k, ii, frg, Nb, No

    real(8) :: coeff_total, coeff_metal

    integer :: Metal_Index = 1
    character(len=5) :: mo_type
    character(len=1) :: dtype
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    
    !Auxiliar
    real(8),dimension(:,:),allocatable :: Aux
    real(8),dimension(:),allocatable   :: Vec

    !OM selection stuff
    integer :: Nselect
    integer,dimension(:),allocatable :: MO_select
    integer,dimension(1:10) :: Nat
    real(8) :: OV, OVmax
    integer :: kmax
    logical :: do_fast=.false.

    !NTO averaging
    integer :: Nel
    real(8) :: Norm
    real(8),dimension(1:10) :: hole, particle
    real(8),dimension(1:BASIS_SIZE) :: W, W2
    
    ! Set verbosity
    verbose = 1

    call parse_input(fchkfile1,fchkfile2,mo_sel,do_fast)
    
    ! Open first file 
    open(IFCHK,file=fchkfile1,status="old")
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
    allocate(MO1(1:Nb,1:No),S1(1:Nb,1:Nb),MOc(1:Nb,1:No))
    k=0
    do i=1,No
    do j=1,Nb
        k=k+1
        MO1(j,i) = A(k)
    enddo
    enddo
    deallocate(A)
    call overlap_matrix(Nb,No,MO1,S1)
!     call MAT0(99,S1,Nb,Nb,"S1")
!     call MAT0(98,MO1,Nb,No,"C")
    deallocate(MOc)
    close(IFCHK)
    
    if (fchkfile2/="none") then
        ! Open second file 
        open(IFCHK,file=fchkfile2,status="old")
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
        allocate(MO2(1:Nb,1:No),S2(1:Nb,1:Nb),MOc(1:Nb,1:No))
        k=0
        do i=1,No
        do j=1,Nb
            k=k+1
            MO2(j,i) = A(k)
        enddo
        enddo
        deallocate(A)
        call overlap_matrix(Nb,No,MO2,S2)
    !     call MAT0(99,S,Nb,Nb,"S")
    !     call MAT0(98,MO2,Nb,No,"C")
        deallocate(MOc)
        close(IFCHK)
    endif
    
    
    ! Get selection
    allocate(MO_select(1:No))
    if (adjustl(mo_sel) == "all") then
        Nselect = No
        MO_select  = (/(i,i=1,No)/)
    else
        call selection2intlist(mo_sel,MO_select,Nselect)
    endif

    ! Allocate auxiliars for overlap 
    allocate(Aux(nb,nb),Vec(nb))
    if (fchkfile2/="none" .and. .not.do_fast) then
        ! Compute S1^1/2 and S1^1/2
        ! S^1/2 = U (Sd)^1/2 U^t
        !
        ! S1^1/2 (store result in S1)
        call diagonalize_full(S1(1:nb,1:nb),nb,Aux(1:nb,1:nb),Vec(1:nb),"lapack")
        do i=1,nb
            do j=1,nb
                S1(i,j) = 0.d0
                do k=1,nb
!                     if (dabs(Vec(k))<1d-10) cycle
                    S1(i,j) = S1(i,j) + Aux(i,k)*dsqrt(Vec(k))*Aux(j,k)
                enddo
            enddo
        enddo
        ! S2^1/2 (store result in S2)
        call diagonalize_full(S2(1:nb,1:nb),nb,Aux(1:nb,1:nb),Vec(1:nb),"lapack")
        do i=1,nb
            do j=1,nb
                S2(i,j) = 0.d0
                do k=1,nb
!                     if (dabs(Vec(k))<1d-10) cycle
                    S2(i,j) = S2(i,j) + Aux(i,k)*dsqrt(Vec(k))*Aux(j,k)
                enddo
            enddo
        enddo
        !
        ! Get S1^1/2 S2^1/2
        Aux=matrix_product(nb,nb,nb,S1,S2)
    else
        Aux=S1
    endif

    do i=1,Nselect
        ii=MO_select(i)
        print*, "MO:", ii
        if (verbose>1) then
        print*, verbose
            print*, "  Bf#    MO1    MO2"
            print*, "----------------------------------"
            OV = 0.d0
            do j=1,Nb
                if (fchkfile2/="none") then
                    OV = OV + MO1(j,ii)*MO2(j,ii)
                else
                    do k=1,Nb
                    OV = OV + MO1(j,ii)*MO1(j,ii)
                    enddo
                endif
                if (dabs(MO1(j,ii)) > 0.1d0) then
                    if (fchkfile2=="none") then
                        print'(I5,2X,F10.4)', j, MO1(j,ii)
                    else
                        print'(I5,2X,2F10.4)', j, MO1(j,ii), MO2(j,ii)
                    endif
                endif
            enddo
            print*, "----------------------------------"
!             print'(X,A,F8.3)', "Overlap(Uncorrected)=", OV
        endif
        ! Compute overlap as 
        !  OV = MO1^t S1^1/2 S2^1/2 MO2
        Vec = matrix_vector_product(nb,nb,Aux,MO1(:,ii),tA=.true.)
        OVmax=0.d0
        do k=1,No
            if (fchkfile2/="none") then
                OV = vector_dot_product(nb,Vec,MO2(:,k))
            else
                OV = vector_dot_product(nb,Vec,MO1(:,k))
            endif
            if (dabs(OV)>dabs(OVmax)) then
                OVmax=OV
                kmax=k
            endif
            if (dabs(OV)>0.8d0) exit
        enddo
        print'(X,A,F8.3,A,I0,A)', "Overlap=", OVmax, "(with OM2:", kmax, ")"
        
        if (fchkfile2/="none") then
            ! Inverse coeff if negative overlap
            if (OV<0) MO2(:,ii)=-MO2(:,ii)
        endif
    enddo

    if (fchkfile2/="none") then
        ! Rewrite MO2 section
        allocate(A(nb*no))
        k=0
        do i=1,No
        do j=1,Nb
            k=k+1
            A(k) = MO2(j,i)
        enddo
        enddo
        call write_fchk(99,"Alpha MO coefficients",'R',nb*no,A,(/0/),j)
    endif

     
    stop

    contains

    subroutine parse_input(inpfile,inpfile2,mo_sel,do_fast)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,inpfile2, mo_sel
        logical,intent(inout)          :: do_fast
        
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
                case ("-f2")
                    call getarg(i+1, inpfile2)
                    argument_retrieved=.true.

                case ("-mo")
                    call getarg(i+1, mo_sel)
                    argument_retrieved=.true.
                    
                case ("-fast")
                    do_fast=.true.
                case ("-nofast")
                    do_fast=.false.

                case ("-h")
                    need_help=.true.
                    
                case ("-v")
                    verbose=2

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo


       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '              p r i n t  M O '
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-f2             ', trim(adjustl(inpfile2))
        write(0,*) '-mo             ', trim(adjustl(mo_sel))
        write(0,*) '-[no]fast      ',  do_fast
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


end program printMO

    

