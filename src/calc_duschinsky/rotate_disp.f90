program rotate_disp

!     use matrix !-- TO BE DONE
    implicit none

    integer,parameter :: NDIM=600

    integer :: Nvib=0
    real(8),dimension(:),allocatable   :: GK, GK2
    real(8),dimension(:,:),allocatable :: GM
    
    ! Files IO
    character(len=200) :: dispfile="displacement.dat",&
                          duschfile="duschinsky.dat",&
                          outdisp="displacement_rot.dat"
    integer :: I_DIS=10, &
               I_DUS=11, &
               O_DIS=20
    integer :: IOstatus

    !Counters
    integer :: i,j


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(Nvib,dispfile,duschfile,outdisp)


    ! Allocate arrays
    if (Nvib/=0) then
        allocate(GK(1:Nvib),GK2(1:Nvib),GM(1:Nvib,1:Nvib))
    else
        allocate(GK(1:NDIM),GK2(1:NDIM),GM(1:NDIM,1:NDIM))
    endif

    ! Read Duschinsky and displacement vector
    open(I_DUS,file=duschfile)
    open(I_DIS,file=dispfile)
    if (Nvib/=0) then
        do i=1,Nvib
            do j=1,Nvib
            read(I_DUS,*) GM(i,j)
            enddo
            read(I_DIS,*) GK(i)
        enddo
    else ! guess Nvib from dispfile
        i=0
        do
            i=i+1
            read(I_DIS,*,iostat=IOstatus) GK(i)
            if (IOstatus /= 0) exit
        enddo
        Nvib=i-1
        do i=1,Nvib
            do j=1,Nvib
            read(I_DUS,*) GM(i,j)
            enddo
        enddo
    endif
    close(I_DIS)
    close(I_DUS)

    ! The inverse rotation requires the inverse Duschinsky matrix
    ! (do the inverse,as with internal it is not exactly orthogonal)
!     GM(1:Nvib,1:Nvib) = inverse_realgen(Nvib,GM) !--TO BE DONE

    ! Rotate displacement to the final state with Duschinsky matrix
    ! K2 = J^-1 * K1 (but keep the sign)
    ! For now, we use J^t
    do i=1,Nvib
        GK2(i) = 0.d0
        do j=1,Nvib
            GK2(i) = GK2(i) + GM(j,i)*GK(j) 
        enddo
    enddo

    ! Print rotated displacement
    open(O_DIS,file=outdisp)
    do i=1,Nvib
        write(O_DIS,*) GK2(i)
    enddo
    close(O_DIS)

    stop

    contains

    subroutine parse_input(Nvib,dispfile,duschfile,outdisp)

        ! Arguments
        integer,intent(inout) :: Nvib
        character(len=*),intent(inout) :: dispfile,duschfile,outdisp

        !Local
       ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg


        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-disp") 
                    call getarg(i+1, dispfile)
                    argument_retrieved=.true.
                case ("-dus") 
                    call getarg(i+1, duschfile)
                    argument_retrieved=.true.
                case ("-o") 
                    call getarg(i+1, outdisp)
                    argument_retrieved=.true.
                case ("-nv") 
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    read(arg,*) Nvib
                case ("-h")
                    need_help=.true.

                case default
                    print*, "ERROR: Unkown command line argument: "//trim(adjustl(arg))
                    stop
            end select
        enddo 

       !Print options (to stdx)
        write(0,'(/,A)')    '========================================================'
        write(0,'(/,A)')    '             R O T A T E   D I S P '    
        write(0,'(/,A)')    '  Rotate displacement vector to the final state'        
!         call print_version()
        write(0,'(/,A)')    '========================================================'
        write(0,'(/,A)')    '-------------------------------------------------------------------'
        write(0,'(A)')      ' Flag         Description                      Value'
        write(0,'(A)')      '-------------------------------------------------------------------'
        write(0,'(X,A,I0)') '-nv          Number of vibrational coods      ', Nvib
        write(0,*)          '-disp        File with displacement vector    ', trim(adjustl(dispfile))
        write(0,*)          '-dus         File with Duschinsky matrix      ', trim(adjustl(duschfile))
        write(0,*)          '-o           Output file with displacement    ', trim(adjustl(outdisp))
        write(0,*)          '-h           This help                       ',  need_help
        write(0,*)          '-------------------------------------------------------------------'
        if (need_help) stop

        return


    end subroutine parse_input

end program rotate_disp