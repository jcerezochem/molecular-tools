program test_read_log

    !Compilation
    ! $FC ../modules/alerts.f90 ../modules/line_preprocess.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/gaussian_manage_vNoTypes.f90 read_state_log.f90 -o read_state_log.exe -cpp -DDOUBLE

    !NOTES
    ! Normal modes matrix indeces: T(1:Nvib,1:3Nat)
    ! This might be the contrary as the usual convention
    ! (anywaym who cares, since we use a Tvector)

    use gaussian_manage_notypes
    use alerts
    use line_preprocess
    use constants

    !Interesting info..
    integer :: Nat, Nvib
    real(8),dimension(1:1000) :: GEOM, Gvector, mu
    real(8),dimension(:),allocatable :: Tvector
    real(8),dimension(1:1000,1:1000) :: T
    character(len=5) :: PG

    !Auxiliars
    character(len=100) :: section
    integer :: N, N_T
    character :: dtype, cnull
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: I
    integer :: error
    real(8) :: Norm
    !Counters
    integer :: j,k

    !I/O
    character(len=100) :: logfile, outfile
    integer :: ios
    integer :: I_LOG =11,&
               O_STA =20

    ! Get filenames from commandline or read standard input/output
    if ( iargc() == 0 ) then
        I_LOG = 5
        O_STA  = 6
    else
        call getarg(1, logfile)
        open(I_LOG,file=logfile,status="old",iostat=ios)
        if ( ios /= 0 ) call alert_msg("fatal","could not open the file "//trim(adjustl(logfile)))
    endif

    !INFORMATION IN THE LOG FILE
    ! Natoms
    call get_Natoms(I_LOG,Nat)

    ! Nvib
    call get_PG(I_LOG,PG)
    !Determine linearity according to point group (* means \infty)
    if ( INDEX(PG,"*") /= 0 ) then
        !Linear molecule 
        Nvib = 3*Nat - 5
    else
        !Non Linear molecule
        Nvib = 3*Nat - 6
    endif

    ! Geom
    call get_ori_geom(I_LOG,GEOM(1:3*Nat),"Standard orientation",error)
    if ( error /= 0 ) call get_ori_geom(I_LOG,GEOM(1:3*Nat),"Input orientation",error)

    ! Freq and Normal modes
    N_T = 3*Nat * Nvib
    allocate( Tvector(1:N_T) )
    call read_freq_NT(I_LOG,Nvib,Nat,Gvector(1:Nvib),mu(1:Nvib),Tvector(1:N_T),error)
    if (error == 1) then
        call alert_msg("warning", "Normal modes (T) matrix in low precision. This will lead to poor results")
    elseif (error == -1) then
        call alert_msg("warning","No frequency information in the log file: only geom will be written")
    endif
    close(I_LOG)

    ! Reconstruct Lcart (non-symmetric)
    l=0
    do j=1,Nvib
        do k=1,3*Nat
            l=l+1
            T(j,k) = Tvector(l)
         enddo
    enddo
    ! Rebuild Tvector in FCclasses compliant way
    l=0
    do j=1,3*Nat
        do k=1,Nvib
            l=l+1
            Tvector(l) = T(k,j)
         enddo
    enddo


    !WRITE STATE FILE
    if (O_STA /= 6) then
        call split_line(logfile,".",outfile,cnull)
        outfile = "state_"//trim(adjustl(outfile))//"_log"
        open(O_STA,file=outfile,status="replace")
    endif
    do j=1,3*Nat
        write(O_STA,'(E17.8)') GEOM(j)
    enddo
    if (error == -1) stop
    do j=1,N_T
        write(O_STA,'(E17.8)') Tvector(j)
    enddo
    do j=1,Nvib
        write(O_STA,'(F10.4)') Gvector(j)
    enddo
    close(O_STA)


    stop

end program test_read_log

