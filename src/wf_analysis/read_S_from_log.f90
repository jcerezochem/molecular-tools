program read_S

    use line_preprocess

    implicit none

    character(len=200) :: line
    real(8),dimension(:),allocatable :: vec
    real(8),dimension(:,:),allocatable :: S
    !IO
    integer :: I_LOG=5
    integer :: IOstatus
    ! Counters
    integer :: i, j, Nbasis, N, Nc


    do 
        read(I_LOG,'(A)') line
        if (INDEX(line,'basis functions,') /= 0) exit
    enddo
    read(line,*) Nbasis
    allocate( vec(1:Nbasis+1), S(1:Nbasis,1:Nbasis) )

    do 
        read(I_LOG,'(A)') line
        if (INDEX(line,'SSO for IR=          1') /= 0) exit
    enddo

    i=0
    read(I_LOG,'(A)') line
    call string2vector(line,vec,Nc,' ')
    do 
        do j=1,Nbasis
            read(I_LOG,'(A)',iostat=IOstatus) line
            call string2vector(line,vec,N,' ')
            S(j,i+1:i+Nc) = vec(2:N)
        enddo
        i = i + Nc
        if (i == Nbasis) exit
        read(I_LOG,'(A)') line
        call string2vector(line,vec,Nc,' ')
    enddo
            
    do i=1,Nbasis
        print'(10000G15.6)', S(i,1:Nbasis)
    enddo

    stop


end program read_S



