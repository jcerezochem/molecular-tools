program get_adiabaticE

    !Compilation
    ! $FC ../modules/alerts.f90 ../modules/gaussian_fchk_manage_v2.f90 get_adiabaticE.f90 -o get_adiabaticE.exe -cpp -DDOUBLE 

    !NOTES
    ! Normal modes matrix indeces: T(1:Nvib,1:3Nat)
    ! This might be the contrary as the usual convention
    ! (anywaym who cares, since we use a Tvector)

    use gaussian_fchk_manage
    use alerts
    use line_preprocess

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
    integer :: Nat, Nvib
    real(8) :: E0, E1, DeltaE
    

    !Auxiliars
    character(len=100) :: section
    integer :: N
    character :: dtype, cnull
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    integer :: error
    !Counters
    integer :: j,k

    !I/O
    character(len=100) :: fchkfile1, fchkfile2
    integer :: I_FCHK=10

    ! Get filenames from commandline or read standard input
    if ( iargc() /= 2 ) then
        call alert_msg("fatal","FCHK files fo optimized GS and Exc required (in that order)")
    else
        call getarg(1, fchkfile1)
        call getarg(2, fchkfile2)
    endif

    open(I_FCHK,file=fchkfile1,status="old",iostat=ios)
    if ( ios /= 0 ) call alert_msg("fatal","could not open the file "//trim(adjustl(fchkfile1)))

    !INFORMATION IN THE FCHK FILE
    ! GS energy
!     call read_fchk(I_FCHK,"Total Energy",dtype,N,A,IA,error)
    call read_fchk(I_FCHK,"SCF Energy",dtype,N,A,IA,error)
    if (error == 0) then
        E0 = A(1)
        deallocate(A)
    endif

    close(I_FCHK)

    open(I_FCHK,file=fchkfile2,status="old",iostat=ios)
    if ( ios /= 0 ) call alert_msg("fatal","could not open the file "//trim(adjustl(fchkfile2)))

    !INFORMATION IN THE FCHK FILE
    ! TD energy
!     call read_fchk(I_FCHK,"ETran state values",dtype,N,A,IA,error)
    call read_fchk(I_FCHK,"CIS Energy",dtype,N,A,IA,error)
    if (error == 0) then
        E1 = A(1)
        deallocate(A)
    endif

    DeltaE = (E1 - E0)*27.2116d0

    print*, DeltaE


    stop

end program get_adiabaticE

