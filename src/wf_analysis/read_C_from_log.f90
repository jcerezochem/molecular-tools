program read_C

    use line_preprocess

    implicit none

    character(len=200) :: line
    real(8),dimension(:),allocatable :: vec
    real(8),dimension(:,:),allocatable :: C
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
    allocate( vec(1:Nbasis+1), C(1:Nbasis,1:Nbasis) )

    do 
        read(I_LOG,'(A)') line
        if (INDEX(line,'Molecular Orbital Coefficients:') /= 0) exit
    enddo

    i=0
    read(I_LOG,'(A)') line
    call string2vector(line,vec,Nc,' ')
    read(I_LOG,'(A)') line ! O/V (type)
    read(I_LOG,'(A)') line ! Eigenvalues
    do 
        do j=1,Nbasis
            read(I_LOG,'(A)',iostat=IOstatus) line
            read(line,'(21X,1000F10.5)') C(j,i+1:i+Nc)
        enddo
        i = i + Nc
        if (i == Nbasis-1) exit
        read(I_LOG,'(A)') line
        call string2vector(line,vec,Nc,' ')
        read(I_LOG,'(A)') line ! O/V (type)
        read(I_LOG,'(A)') line ! Eigenvalues
    enddo
            
    do i=1,Nbasis
        print'(10000G15.6)', C(i,1:Nbasis)
    enddo

    stop


end program read_C



