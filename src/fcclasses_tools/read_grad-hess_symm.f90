program read_grad_hess_symm

    !Compilation
    ! $FC ../modules/alerts.f90 ../modules/line_preprocess.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/gaussian_manage_vNoTypes.f90 read_state_fchk.f90 -o read_state_fchk.exe -cpp -DDOUBLE 
    ! THIS IS A SELF CONTAINED VERSION (no needs molecular tools)

    !NOTES
    ! Normal modes matrix indeces: T(1:Nvib,1:3Nat)
    ! This might be the contrary as the usual convention
    ! (anywaym who cares, since we use a Tvector)

!     use gaussian_fchk_manage
    use alerts
    use line_preprocess
    use symmetry_mod_notypes

    implicit none

    real(8),parameter :: BOHRtoAMS = 5.2917720859D-1
    !CONVERSIONS
    double precision, parameter :: uma_to_Kg=1.66053873d-27,  &
                                   H_to_J=4.3597482d-18,      &
                                   bohr_to_m=5.291772083d-11, &
                                   A_to_m=1d-10

    !CONSTANTS
    double precision, parameter :: pi=4.d0*datan(1.d0),       &
                                   c_vac=2.99792458d8,        &
                                   plank=6.62606957d-34   

    !Interesting info..
    integer :: Nat
    real(8),dimension(1:1000) :: GEOM, GRAD
    real(8),dimension(1:1000,1:1000) :: HESS
    integer,dimension(1:1000) :: isym
    

    !Auxiliars
    character(len=100) :: section
    integer :: N, N_T
    character :: dtype, cnull
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: I
    integer :: error
    real(8) :: Norm
    !Counters
    integer :: ii,j,k,l

    !I/O
    character(len=100) :: fchkfile, outfile
    integer :: ios
    integer :: I_FCHK=10,&
               O_STA =20

    ! Get filenames from commandline or read standard input/output
    if ( iargc() == 0 ) then
        I_FCHK = 5
        O_STA  = 6
    else
        call getarg(1, fchkfile)
        open(I_FCHK,file=fchkfile,status="old",iostat=ios)
        if ( ios /= 0 ) call alert_msg("fatal","could not open the file "//trim(adjustl(fchkfile)))
    endif

    !INFORMATION IN THE FCHK FILE
    ! Natoms
    call read_fchk(I_FCHK,"Number of atoms",dtype,N,A,I,error)
    if (error == 0) then
        Nat = I(1)
        deallocate(I)
    endif

    ! Geom
    call read_fchk(I_FCHK,"Current cartesian coordinates",dtype,N,A,I,error)
    if (error == 0) then
        GEOM(1:N) = A*BOHRtoAMS
        deallocate(A)
    endif


    ! Grad
    call read_fchk(I_FCHK,"Cartesian Gradient",dtype,N,A,I,error)
    if (error == 0) then
        GRAD(1:N) = A
        deallocate(A)
    endif


    ! Hessian
    call read_fchk(I_FCHK,"Cartesian Force Constants",dtype,N,A,I,error)
    if (error == 0) then
        k=0
        do ii=1,3*Nat
            do j=1,ii
                k=k+1
                Hess(ii,j) = A(k) 
                Hess(j,ii) = Hess(ii,j)
            enddo
        enddo
        deallocate(A)
    endif

    call symm_atoms_nt(geom,Nat,isym)
    do ii=1,Nat
        k=ii*3-2
        print*, "Atom", ii, "Sym", isym(ii)
        print'(A,3(F8.3))', "  ", GRAD(k:k+2) 
        k=isym(ii)*3-2
        print'(A,3(F8.3))', "  ", GRAD(k:k+2)
        print*, "==========" 
    enddo

    do ii=1,Nat
    do j=1,ii
        k=ii*3-2
        l=j*3-2
        print*, "Atoms", ii, j, "Syms", isym(ii), isym(j)
        print'(A,3(G15.3))', "  ", Hess(k,l), Hess(k+1,l+1), Hess(k+2,l+2)
        k=isym(ii)*3-2
        l=isym(j) *3-2
        print'(A,3(G15.3))', "  ", Hess(k,l), Hess(k+1,l+1), Hess(k+2,l+2)
        print*, "==========" 
    enddo
    enddo


    close(O_STA)


    stop

    contains

    subroutine read_fchk(unt,section,data_type,N,A,I,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 1.0/March 2012)
        !==============================================================
        !Description
        ! Generic SR to read any section of the checkpoint
        ! Enter deallocated arrays
        !Arguments
        ! unt (int;in): unit number of the log file
        ! section(char,in): name of the section to be read
        ! data_type(char,out): Integer (I) or Real (R) data read
        ! A(real,dimension(:)): Real array to store real data
        ! I(integer,dimension(:)): Int array to store int data
        ! error_flag(integer,out): 0: success
        !                          1: section not found
        !==============================================================

        integer,intent(in) :: unt
        character(len=*),intent(in) :: section
        character(len=1),intent(out) :: data_type
        integer,intent(out) :: N
! #ifdef DOUBLE
        double precision, dimension(:), allocatable, intent(out) :: A
! #else
!         real, dimension(:), allocatable, intent(out) :: A
! #endif
        integer,dimension(:), allocatable, intent(out) :: I
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=42) :: section_full
        character(len=1) :: is_array
        character(len=40) :: cdata
        !I/O
        integer :: IOstatus
        
        
        ! Search section
        error_flag = 0
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    call alert_msg("note","Section '"//trim(adjustl(section))//"' not present in the FCHK file.")
                    error_flag=1
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                    read(line,'(A42)') section_full
                    if (adjustl(section_full) == adjustl(section)) exit
                endif
        enddo

        !Get info about section from line just read
        read(line,'(A42,X,A1,3X,A1,X,A)') section_full, data_type, is_array, cdata
        if (is_array /= "N") then
            !Is not an array
            N=1
            if ( data_type == "R" ) then 
                allocate( A(1:1) )
                read(cdata,*) A(1)
            elseif ( data_type == "I" ) then
                allocate( I(1:1) )
                read(cdata,*) I(1) 
            endif
        else
            read(cdata,*) N
            if ( data_type == "R" ) then
                allocate( A(1:N) )
                read(unt,*) A(1:N)
            elseif ( data_type == "I" ) then
                allocate( I(1:N) )
                read(unt,*) I(1:N)
            endif
        endif 

        rewind(unt)
        return

    end subroutine read_fchk

end program read_grad_hess_symm

