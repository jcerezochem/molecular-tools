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
    integer :: i, j, Nbasis, N


    ! Get Nbasis first
    do 
        read(I_LOG,'(A)') line
        if (INDEX(line,'[MO]') /= 0) exit
    enddo
    read(I_LOG,*) line ! Ene=
    read(I_LOG,*) line ! Spin=
    read(I_LOG,*) line ! Occup=
    do 
        read(I_LOG,'(A)') line
        if (INDEX(line,'Ene=') /= 0) exit
        read(line,*) Nbasis
    enddo
    rewind(I_LOG)
    allocate( C(1:Nbasis,1:Nbasis) )
    write(0,*) "Basis", Nbasis

    do 
        read(I_LOG,'(A)') line
        if (INDEX(line,'[MO]') /= 0) exit
    enddo
    do i=1,Nbasis !(=Nmo)
        read(I_LOG,*) line ! Ene=
        read(I_LOG,*) line ! Spin=
        read(I_LOG,*) line ! Occup=
        do j=1,Nbasis
            read(I_LOG,*) N, C(j,i)
        enddo
    enddo

    do i=1,Nbasis
        print'(10000G15.6)', C(i,1:Nbasis)
    enddo

    stop


end program read_C



