program get_MO_energy

    !==============================================================
    ! This code uses MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    ! Description
    !  Utility to retrieve NTO coeficients
    ! Compilation intructions:
    !  gfortran ../modules/alerts.f90 ../modules/line_preprocess.f90 ../modules/gaussian_fchk_manage_v2.f90 get_MO_energy.f90 -o get_MO_energy.exe -cpp 
    !
    !==============================================================

    use gaussian_manage
    use line_preprocess


    implicit none

    

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
    integer :: error, N
    !====================== 

    character(len=150) :: fileinp 
    integer :: I_FCHK=10

    real(8),dimension(1000) :: MO_E
    integer,dimension(1000) :: alpha, beta
    integer :: nalpha,nbeta
    integer :: i


    call parse_input(fileinp,alpha,nalpha,beta,nbeta)

    open(I_FCHK,file=fileinp,status="old")

    !Alpha
    N=0
    call read_fchk(I_FCHK,"Alpha Orbital Energies",dtype,N,A,IA,error)
    if (error == 0) then
        MO_E(1:N) = A(1:N)
        deallocate(A)
    endif

    if (N/=0 .and. nalpha/=0) then
        print*, "Alpha orbital energies (a.u)"
        do i=1,nalpha
            print*,  alpha(i), MO_E(alpha(i))
        enddo
    endif
    print*, ""

    !Beta
    N=0
    call read_fchk(I_FCHK,"Beta Orbital Energies",dtype,N,A,IA,error)
    if (error == 0) then
        MO_E(1:N) = A(1:N)
        deallocate(A)
    endif

    if (N/=0 .and. nbeta/=0) then
        print*, "Beta orbital energies (a.u)"
        do i=1,nbeta
            print*,  beta(i), MO_E(beta(i))
        enddo
    endif
    print*, ""

    stop

    contains

    subroutine sort_vec_max(V,N)

        implicit none

#ifdef DOUBLE
        double precision,dimension(:),intent(inout) :: V
        double precision :: aux
#else
        real,dimension(:),intent(inout) :: V
        real :: aux
#endif  
!         integer,dimension(:),intent(inout) :: IORD
        integer,intent(in) :: N
        integer :: i,j, iaux

!         !Intialize IORD
!         do i=1,N
!             IORD(i) = i
!         enddo

        do i=1,N-1
            do j=i+1,N
                if (V(j)>V(i)) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                    !Track the index permutations in IORD
!                     iaux = IORD(i)
!                     IORD(i) = IORD(j)
!                     IORD(j) = iaux
                endif
            enddo
        enddo

        return

    end subroutine sort_vec_max    

    subroutine parse_input(inpfile,alpha,nalpha,beta,nbeta)

    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile
        integer,dimension(:),intent(inout) :: alpha, beta
        integer,intent(inout) :: nalpha, nbeta
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg, alpha_list="", beta_list=""

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-f") 
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.
!                 case ("-fti") 
!                     call getarg(i+1, filetype_inp)
!                     argument_retrieved=.true.

                case ("-a")
                    call getarg(i+1, alpha_list)
                    argument_retrieved=.true.

                case ("-b")
                    call getarg(i+1, beta_list)
                    argument_retrieved=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 


        ! Define alpha and beta orbital lists
        !--------------------------------------
        call string2vectorINT(alpha_list,alpha,nalpha)
        call string2vectorINT(beta_list,beta,nbeta)

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '              M O   E N E R G Y '    
        write(0,'(/,A)') '        Get MO energies from fchk files '       
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-a              ', trim(adjustl(alpha_list))
        write(0,*) '-b              ', trim(adjustl(beta_list))
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return

    end subroutine parse_input  

    subroutine string2vectorINT(raw_vector,array_vector,n_elem)

        !Description
        ! Tranforms a string of comma sepparated values into an
        ! array of such real vaues 

        character(len=*),intent(inout) :: raw_vector
        integer,dimension(:),intent(out) :: array_vector
        integer,intent(out) :: n_elem

        !Local
        character(len=len_trim(raw_vector)) :: copy_vector
        character(len=240) :: auxchar
        integer :: i,i_max,i_min, ierr
    
        n_elem = 0
        
        !Read unknown length vector two ways
        
        copy_vector=adjustl(raw_vector)

        !1. comma sepparated
        if ( INDEX(copy_vector,',') /= 0 ) then
            i=0
            do 
                i=i+1
                if ( INDEX(copy_vector,',') /= 0 ) then
                    call split_line(copy_vector,',',auxchar,copy_vector)
                    read(auxchar,*) array_vector(i)
                else 
                    read(copy_vector,*) array_vector(i)
                    exit
                endif
            enddo  
            n_elem=i
        !2. i-j range
        else if ( INDEX(copy_vector,'-') /= 0 ) then
            call split_line(copy_vector,'-',auxchar,copy_vector)
            read(auxchar,*) i_max 
            read(copy_vector,*) i
            i_min = min(i_max,i)
            i_max = max(i_max,i)
            n_elem = i_max-i_min+1
            array_vector(1:n_elem) = (/ (i,i=i_min,i_max) /)
        !Single orbital 
        else if (len_trim(copy_vector) /= 0) then
            n_elem = 1
            read(copy_vector,*,iostat=ierr) array_vector(1)
        endif

        return

    end subroutine string2vectorINT

end program get_MO_energy

    

