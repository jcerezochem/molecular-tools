program get_isotopologues


    !*****************
    !   MODULE LOAD
    !*****************
    !============================================
    !   Generic
    !============================================
    use io
    use alerts
    use line_preprocess
    use constants 
    use verbosity
    use matrix
    use matrix_print
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
    !============================================
    !  Structure-related modules
    !============================================
    use molecular_structure
    use ff_build
    use atomic_geom
    use symmetry

    implicit none

    integer,parameter :: NDIM = 600

    ! System info
    type(str_resmol) :: molecule
    integer :: nsym
    integer,dimension(NDIM,10) :: isym
    integer,dimension(NDIM) :: nsym_uniq

    ! IO
    character(len=10)  :: ft="guess"
    character(len=200) :: inpfile="input.xyz"
    integer :: I_INP=10
    integer :: IOstatus

    ! Counters
    integer :: i,j,k

    ! Auxiliar
    character :: null
    
    

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,ft)

 
    ! 1. READ INPUT
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)

    ! STRUCTURE FILE
    print'(X,A)', "READING STATE1 FILE (STRUCTURE)..."
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,molecule)
    close(I_INP)

    ! 2. GUESS SYMMETRY AND WRITE
    ! -------------------------------
    call symm_atoms_v2(molecule,isym)
    do i=1,molec%natoms
        print*, i, isym(i,1:nsym)
    enddo

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft
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
                case ("-f") 
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.
                case ("-ft") 
                    call getarg(i+1, ft)
                    argument_retrieved=.true.        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        ! Some checks on the input
        !----------------------------

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '        G E T    I S O T O P O L O G U E S '    
        write(0,'(/,A)') '   Get the number of isotopologues of a molecule'         
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(ft))
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input

end program get_isotopologues