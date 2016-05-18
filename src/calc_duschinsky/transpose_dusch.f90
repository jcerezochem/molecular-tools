program transpose_dusch

!     use matrix !-- TO BE DONE
    implicit none

    integer,parameter :: NDIM=600

    integer :: Nvib=0
    real(8),dimension(:),allocatable   :: GK, GK2
    real(8),dimension(:,:),allocatable :: GM
    
    ! Files IO
    character(len=200) :: duschfile="duschinsky.dat",&
                          outdusch="duschinsky_inv.dat"
    integer :: I_DUS=11, &
               O_DUS=20
    integer :: IOstatus
    real(8) :: r

    !Counters
    integer :: i,j


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(Nvib,duschfile,outdusch)


    ! Allocate arrays
    if (Nvib/=0) then
        allocate(GM(1:Nvib,1:Nvib))
    else
        allocate(GM(1:NDIM,1:NDIM))
    endif

    ! Read Duschinsky and displacement vector
    open(I_DUS,file=duschfile)
    if (Nvib/=0) then
        do i=1,Nvib
            do j=1,Nvib
            read(I_DUS,*) GM(i,j)
            enddo
        enddo
    else ! guess Nvib from dispfile
        i=0
        do
            i=i+1
            read(I_DUS,*,iostat=IOstatus) r
            if (IOstatus /= 0) exit
        enddo
        Nvib=int(sqrt(dfloat(i-1)))
        rewind(I_DUS)
        do i=1,Nvib
            do j=1,Nvib
            read(I_DUS,*) GM(i,j)
            enddo
        enddo
    endif
    close(I_DUS)

    ! The inverse rotation requires the inverse Duschinsky matrix
    ! (do the inverse,as with internal it is not exactly orthogonal)
!     GM(1:Nvib,1:Nvib) = inverse_realgen(Nvib,GM) !--TO BE DONE

    ! Print traspose displacement
    open(O_DUS,file=outdusch)
    do i=1,Nvib
        do j=1,Nvib
        write(O_DUS,*) GM(j,i)
        enddo
    enddo
    close(O_DUS)

    stop

    contains

    subroutine parse_input(Nvib,duschfile,outdusch)

        ! Arguments
        integer,intent(inout) :: Nvib
        character(len=*),intent(inout) :: duschfile,outdusch

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
                case ("-dus") 
                    call getarg(i+1, duschfile)
                    argument_retrieved=.true.
                case ("-o") 
                    call getarg(i+1, outdusch)
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
        write(0,'(/,A)')    '          T R A N S P O S E    D U S C H '    
        write(0,'(/,A)')    '  Transpose Duschinsky (to inverse the transformation)'        
!         call print_version()
        write(0,'(/,A)')    '========================================================'
        write(0,'(/,A)')    '-------------------------------------------------------------------'
        write(0,'(A)')      ' Flag         Description                      Value'
        write(0,'(A)')      '-------------------------------------------------------------------'
        write(0,'(X,A,I0)') '-nv          Number of vibrational coods      ', Nvib
        write(0,*)          '-dus         File with Duschinsky matrix      ', trim(adjustl(duschfile))
        write(0,*)          '-o           Output file with Duschinsky      ', trim(adjustl(outdusch))
        write(0,*)          '-h           This help                       ',  need_help
        write(0,*)          '-------------------------------------------------------------------'
        if (need_help) stop

        return


    end subroutine parse_input

end program transpose_dusch
