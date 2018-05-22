program geom_param_list

    !==============================================================
    ! This code uses MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    ! Description
    !  Input a list of geometric parameters as comma separated atom
    !  numbers. The parameter to be calculated is decided from the
    !  the number of itmes: 2=bond, 3=angle, 4=dihedral angle
    !==============================================================

    !*****************
    !   MODULE LOAD
    !*****************
    !============================================
    !   Generic
    !============================================
    use alerts
    use line_preprocess
    use constants 
    use verbosity
    use matrix
    use matrix_print
    use io, only: uout
    !============================================
    !   Structure types module
    !============================================
    use structure_types
    !============================================
    !   File readers
    !============================================
    use generic_io
    use generic_io_molec
    use xyz_manage
    use gaussian_manage
    !============================================
    !  Structure-related modules
    !============================================
    use molecular_structure
    use metrics
    use atomic_geom
    use symmetry

    implicit none

    type(str_resmol) :: molecule
    integer :: i, nlines, nitems
    character(len=100) :: inpfile="test.log",listfile="list.dat"
    character(len=30) :: list_item
    integer,dimension(1:4) :: atlabels
    character(len=10) :: filetype="guess"
    character(len=1) :: dihed_mod
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

    !===========================
    ! Set output unit to stderr
    uout = 0
    !===========================

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(molecule)
    !===========================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,listfile,filetype,labels)

    ! 1. READ DATA
    ! ---------------------------------
    !Structure
    open(IGeom,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line_back(inpfile,".",null,filetype)
    endif

    call generic_strmol_reader(IGeom,filetype,molecule)
    close(IGeom)

    !Read input list
    open(IList,file=listfile,status="old")
    read(IList,*) nlines
    do i=1,nlines
        read(IList,'(A)') list_item
        ! Remove comments
        call split_line(list_item,';',list_item,null)
        ! Get modifies (for improper dihedrals)
        call split_line(list_item,'(',list_item,dihed_mod)
        call selection2intlist(list_item,atlabels,nitems)
        if (nitems == 2) then
            param=calc_atm_dist(molecule%atom(atlabels(1)),molecule%atom(atlabels(2)))
            if (.not.labels) print'(F11.6,2X)', param
            if (labels)      print'(2I5,11X,F12.6)', atlabels(1), atlabels(2), param
        elseif (nitems == 3) then
            param=calc_atm_angle(molecule%atom(atlabels(1)),molecule%atom(atlabels(2)),molecule%atom(atlabels(3)))
            if (.not.labels) print'(F9.2,2X)', param*180.d0/pi
            if (labels)      print'(3I5,6X,F9.2)', atlabels(1), atlabels(2), atlabels(3), param*180.d0/pi
        elseif (nitems == 4 .and. dihed_mod /= 'I') then
            param=calc_atm_dihed_new(molecule%atom(atlabels(1)),molecule%atom(atlabels(2)),molecule%atom(atlabels(3)),&
                             molecule%atom(atlabels(4)))
            if (.not.labels) print'(F9.2,2X)', param*180.d0/pi
            if (labels)      print'(4I5,X,F9.2)', atlabels(1), atlabels(2), atlabels(3), atlabels(4), param*180.d0/pi
        elseif (nitems == 4  .and. dihed_mod == 'I') then
            param=calc_atm_improper(molecule%atom(atlabels(1)),molecule%atom(atlabels(2)),molecule%atom(atlabels(3)),&
                             molecule%atom(atlabels(4)))
            if (.not.labels) print'(F9.2,2X)', param*180.d0/pi
            if (labels)      print'(4I5,A3,X,F9.2)', atlabels(1), atlabels(2), atlabels(3), atlabels(4),'(I)', param*180.d0/pi
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
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '                G E O M   P A R A M '    
        write(0,'(/,A)') '     Automatic geometric parameter generator '     
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input file                       ', trim(adjustl(inpfile))
        write(0,*)       '-ft          \_ FileTyep                      ', trim(adjustl(filetype))
        write(0,*)       '-l           File with coordinates(IC) list   ', trim(adjustl(listfile))
        write(0,*)       '-[no]labels  Show IC labels                  ',  labels
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


end program geom_param_list

    

