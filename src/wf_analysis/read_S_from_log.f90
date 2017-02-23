program read_S

    use line_preprocess

    implicit none

    character(len=200) :: line
    real(8),dimension(:),allocatable :: vec
    real(8),dimension(:,:),allocatable :: S
    real(8),dimension(:), allocatable :: Slt
    !IO
    integer :: I_LOG=5
    integer :: IOstatus
    ! Counters
    integer :: i, j, Nbasis, N, Nc, k
    integer :: icols, nblocks, istart, isum, imax, imin, ii, ib


    do 
        read(I_LOG,'(A)') line
        if (INDEX(line,'basis functions,') /= 0) exit
    enddo
    read(line,*) Nbasis
    N = Nbasis*(Nbasis+1)/2
    allocate( vec(1:Nbasis+1), S(1:Nbasis,1:Nbasis), Slt(1:N) )

    do 
        read(I_LOG,'(A)') line
        if (INDEX(line,'*** Overlap ***') /= 0) exit
    enddo

    !Organized in blocks of 6 cols
    icols = 5
    nblocks = Nbasis/icols
    if (mod(Nbasis,icols) /= 0) nblocks = nblocks+1
    do ib = 1,nblocks
        ! Place the tip on the correct line
        read(I_LOG,'(A)') line ! header
        ! Initialize auxiliar counters
        istart = (ib-1)*icols
        isum = 0
        do i=1,Nbasis
            ! Accumunalte terms (isum) but cycle if the block is not read
            isum = isum + i
            if (i<=istart) cycle
            ! Determine the parts we are goint to read
            ! The j-th row is from 1:sum(1···j)
            ! But at the current round we get for the j-th col 
            !   istart:min(istart+i,istart+6)
            imin = isum-i + istart+1
            ii = min(i-istart,icols)
            imax = isum-i + istart+ii
            read(I_LOG,'(7X,1000(X,D13.6))') Slt(imin:imax)
        enddo
    enddo

    k=0
    do i=1,Nbasis
    do j=1,i
        k=k+1
        S(i,j) = Slt(k)
        S(j,i) = S(i,j)
    enddo
    enddo
    deallocate(Slt)
            
    do i=1,Nbasis
        print'(10000G15.6)', S(i,1:Nbasis)
    enddo

    stop


end program read_S



