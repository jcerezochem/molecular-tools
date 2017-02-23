program substract_disps

!     use matrix !-- TO BE DONE
    implicit none

    integer,parameter :: NDIM=600

    integer :: Nvib=0
    real(8),dimension(:),allocatable   :: GK, GK2
    
    ! Files IO
    character(len=200) :: dispfile1="displacement1.dat",&
                          dispfile2="displacement2.dat",&
                          outdisp="delta_disp.dat"
    integer :: I_DIS=10, &
               O_DIS=20
    integer :: IOstatus

    !Counters
    integer :: i,j


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(Nvib,dispfile1,dispfile2,outdisp)


    ! Allocate arrays
    if (Nvib/=0) then
        allocate(GK(1:Nvib),GK2(1:Nvib))
    else
        allocate(GK(1:NDIM),GK2(1:NDIM))
    endif

    ! Read displacement vectors
    open(I_DIS,file=dispfile1,status='old')
    if (Nvib/=0) then
        do i=1,Nvib
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
    endif
    close(I_DIS)
    open(I_DIS,file=dispfile2,status='old')
    do i=1,Nvib
        read(I_DIS,*) GK2(i)
    enddo
    close(I_DIS)


    !Substract 1-2
    do i=1,Nvib
        GK2(i) = GK(i) - GK2(i)
    enddo

    ! Print rotated displacement
    open(O_DIS,file=outdisp)
    do i=1,Nvib
        write(O_DIS,*) GK2(i)
    enddo
    close(O_DIS)

    stop

    contains

    subroutine parse_input(Nvib,dispfile1,dispfile2,outdisp)

        ! Arguments
        integer,intent(inout) :: Nvib
        character(len=*),intent(inout) :: dispfile1,dispfile2,outdisp

        !Local
       ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        ! iargc type must be specified with implicit none (strict compilation)
        integer :: iargc


        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-disp1") 
                    call getarg(i+1, dispfile1)
                    argument_retrieved=.true.
                case ("-disp2") 
                    call getarg(i+1, dispfile2)
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
        write(0,'(/,A)')    '          S U B S T R A C T   D I S P S'    
        write(0,'(/,A)')    '  Rotate displacement vector to the final state'
        write(0,'(/,A)')    '               (disp1 - disp2)'        
!         call print_version()
        write(0,'(/,A)')    '========================================================'
        write(0,'(/,A)')    '-------------------------------------------------------------------'
        write(0,'(A)')      ' Flag         Description                      Value'
        write(0,'(A)')      '-------------------------------------------------------------------'
        write(0,'(X,A,I0)') '-nv          Number of vibrational coods      ', Nvib
        write(0,*)          '-disp1       File with displacement vector 1  ', trim(adjustl(dispfile1))
        write(0,*)          '-disp2       File with displacement vector 2  ', trim(adjustl(dispfile2))
        write(0,*)          '-o           Output file with displacement    ', trim(adjustl(outdisp))
        write(0,*)          '-h           This help                       ',  need_help
        write(0,*)          '-------------------------------------------------------------------'
        if (need_help) stop

        return


    end subroutine parse_input

end program substract_disps
