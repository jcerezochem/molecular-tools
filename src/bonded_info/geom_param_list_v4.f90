program geom_param_list

    !==============================================================
    ! This code uses MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    ! Description
    !  Input a list of geometric parameters as comma separated atom
    !  numbers. The parameter to be calculated is decided from the
    !  the number of itmes: 2=bond, 3=angle, 4=dihedral angle
    !==============================================================

    !Compilation instructions: now using automake (v4)

    !History
    !v4: adapted to v4 modules (disabled allocation)

    use alerts
    use structure_types
    use line_preprocess
    use gaussian_manage
    use pdb_manage
    use gro_manage
    use atomic_geom

    implicit none

    type(str_resmol) :: molecule
    integer :: i, nlines, nitems
    character(len=100) :: inpfile="test.log",listfile="list.dat"
    character(len=30) :: list_item
    integer,dimension(1:4) :: atlabels
    character(len=10) :: filetype="guess"
    character(len=1) :: null
    logical :: labels = .true.
    !I/O
    integer :: IGeom=10 , &
               IList=11
    integer :: IOstatus
#ifdef DOUBLE
    double precision :: param
#else
    real :: param
#endif


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,listfile,filetype,labels)

    ! 1. READ DATA
    ! ---------------------------------
    !Structure
    open(IGeom,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line(filetype,".",null,filetype)
    endif

    call generic_strfile_read(IGeom,filetype,molecule)
         
    close(IGeom)

    !Read input list
    open(IList,file=listfile,status="old")
    read(IList,*) nlines
    do i=1,nlines
        read(IList,'(A)') list_item
        call string2vector_int(list_item,atlabels,nitems)
        if (nitems == 2) then
            param=calc_dist(molecule%atom(atlabels(1)),molecule%atom(atlabels(2)))
            if (.not.labels) print'(F9.4,2X)', param
            if (labels)      print'(2I5,11X,F9.4)', atlabels(1), atlabels(2), param
        elseif (nitems == 3) then
            param=calc_angle(molecule%atom(atlabels(1)),molecule%atom(atlabels(2)),molecule%atom(atlabels(3)))
            if (.not.labels) print'(F9.2,2X)', param
            if (labels)      print'(3I5,6X,F9.2)', atlabels(1), atlabels(2), atlabels(3), param
        elseif (nitems == 4) then
            param=calc_dihed(molecule%atom(atlabels(1)),molecule%atom(atlabels(2)),molecule%atom(atlabels(3)),&
                             molecule%atom(atlabels(4)))
            if (.not.labels) print'(F9.2,2X)', param
            if (labels)      print'(4I5,X,F9.2)', atlabels(1), atlabels(2), atlabels(3), atlabels(4), param
        endif
        
    enddo

    stop 

    contains

    subroutine parse_input(inpfile,listfile,filetype,labels)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,listfile,filetype
        logical,intent(inout)          :: labels
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false., &
                   do_refine_charges
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
                case ("-f") 
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.

                case ("-ft") 
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.

                case ("-l") 
                    call getarg(i+1, listfile)
                    argument_retrieved=.true.

                case ("-labels")
                    labels=.true.

                case ("-nolabels")
                    labels=.false.

                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 


       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '                GEOM_PARAM '    
        write(0,'(/,A)') '     Automatic geometric parameter generator '        
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(filetype))
        write(0,*) '-l              ', trim(adjustl(listfile))
        write(0,*) '-[no]labels    ',  labels
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(inout) :: filetype
        type(str_resmol),intent(inout) :: molec

        !local axu
        !Read gaussian log auxiliars
        type(str_molprops),allocatable :: props
        character :: null

        select case (adjustl(filetype))
            case("gro")
             call read_gro(unt,molec)
            case("pdb")
             call read_pdb_new(unt,molec)
            case("log")
             allocate(props)
             call parse_summary(unt,molec,props,"struct_only")
             deallocate(props)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

!         call atname2element(molec)

        return

    end subroutine generic_strfile_read


end program geom_param_list

    

