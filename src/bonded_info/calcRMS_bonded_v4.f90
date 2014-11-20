program calcRMS_bonded


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to compute the RMSD for all bonded parameters in a 
    ! strucure file, compared with a reference. Conectivity is guessed
    ! by distances.
    !
    ! Compilation instructions (for mymake script): now using automake (v4)
    !
    ! Change log:
    ! v2: use slightly modified modules. calc_dihed and calc_angle now returns angle in rad, so changes are made accordingly
    ! v4: using v4 modules (disable allocation)
    ! Added RMSD of the whole structure (atomic positions)
    !============================================================================    

    use structure_types
    use line_preprocess
    use ff_build
    use gro_manage
    use pdb_manage
    use alerts
    use constants
    use atomic_geom
    use gaussian_manage
    use gaussian_fchk_manage
    use molecular_structure
    use xyz_manage

    implicit none

    common /ALERT_COUNT/ n_notes, n_errors
    common /GLOBAL_OPTS/ do_refine_charges, I_DB, rename_atoms

    !====================== 
    !Number of error/notes
    integer :: n_notes=0, n_errors=0
    !====================== 

    !====================== 
    !Options 
    logical :: debug=.false., nonH=.false.
    ! Those related to ff_type
    logical :: do_refine_charges=.false.,rename_atoms=.false.
    integer :: I_DB = -1 !This is to use hybridization database
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molec, ref_molec
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: ierr
    character(1) :: null
    character(len=16) :: dummy_char
    !====================== 

    !=============
    !Counters
    integer :: i,j,k
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_REF=11,  &
               O_STR=20  
    !files
    character(len=10) :: ft="guess", ft_ref="guess"
    character(len=200):: inpfile="input.pdb",          &
                         reffile="ref.pdb"
    !status
    integer :: IOstatus
    !===================

    !New things for bonds
    integer,dimension(500,2) :: bond
    integer :: nbonds
    real :: calc, ref, dev, rmsd, dif

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,reffile,ft,ft_ref,debug,nonH)

    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(ft) == "guess") then
        ! Guess file type
        call split_line(inpfile,".",null,ft)
    endif

    call generic_strfile_read(I_INP,ft,molec)
         
    close(I_INP)

    !Read reference structure
    open(I_INP,file=reffile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(reffile)) )

    if (adjustl(ft_ref) == "guess") then
        ! Guess file type
        call split_line(reffile,".",null,ft_ref)
    endif

    call generic_strfile_read(I_INP,ft_ref,ref_molec)

    close(I_INP)


    ! Get connectivity from the residue

    call guess_connect(ref_molec)
    call gen_bonded(ref_molec)

    dev = 0.0
    k = 0
    if (debug) print*, "LIST OF BONDS"
    do i=1,ref_molec%geom%nbonds
        if (nonH) then
            if (adjustl(ref_molec%atom(ref_molec%geom%bond(i,1))%name) == "H" .or. &
                adjustl(ref_molec%atom(ref_molec%geom%bond(i,2))%name) == "H") cycle
        endif
        !Using an external counter in case nonH is used
        k=k+1
        ref  = calc_dist(ref_molec%atom(ref_molec%geom%bond(i,1)),ref_molec%atom(ref_molec%geom%bond(i,2))) 
        calc = calc_dist(molec%atom(ref_molec%geom%bond(i,1)),molec%atom(ref_molec%geom%bond(i,2)))
        dif = abs(calc - ref)
        if (debug) &
        print'(A2,A1,I2,A5,A2,A1,I2,A1,X,2(F7.2,X),F11.6)', &
              ref_molec%atom(ref_molec%geom%bond(i,1))%name, "(", ref_molec%geom%bond(i,1), ") -- ",&
              ref_molec%atom(ref_molec%geom%bond(i,2))%name, "(", ref_molec%geom%bond(i,2), ")", calc, ref, dif
        dev = dev + (calc - ref)**2
    enddo
    rmsd = sqrt(dev/k)
    print'(A,/)', '---------------------' 

    print'(X,A,X,F8.3,/)', "RMSD-bonds (AA):", rmsd

    dev = 0.0
    k = 0
    if (debug) print*, "LIST OF ANGLES"
    do i=1,ref_molec%geom%nangles
        if (nonH) then
            if (adjustl(ref_molec%atom(ref_molec%geom%angle(i,1))%name) == "H" .or. &
                adjustl(ref_molec%atom(ref_molec%geom%angle(i,2))%name) == "H" .or. &
                adjustl(ref_molec%atom(ref_molec%geom%angle(i,3))%name) == "H") cycle
        endif
        !Using an external counter in case nonH is used
        k=k+1
        ref  = calc_angle(ref_molec%atom(ref_molec%geom%angle(i,1)),&
                          ref_molec%atom(ref_molec%geom%angle(i,2)),&
                          ref_molec%atom(ref_molec%geom%angle(i,3)))
        calc = calc_angle(molec%atom(ref_molec%geom%angle(i,1)),&
                          molec%atom(ref_molec%geom%angle(i,2)),&
                          molec%atom(ref_molec%geom%angle(i,3)))
        calc = calc*180.d0/PI
        ref  = ref*180.d0/PI
        dif = abs(calc - ref)
        if (debug) &
        print'(2(A2,A1,I2,A5),A2,A1,I2,A1,X,2(F7.2,X),F11.6)', &
              ref_molec%atom(ref_molec%geom%angle(i,1))%name, "(", ref_molec%geom%angle(i,1), ") -- ",&
              ref_molec%atom(ref_molec%geom%angle(i,2))%name, "(", ref_molec%geom%angle(i,2), ") -- ",&
              ref_molec%atom(ref_molec%geom%angle(i,3))%name, "(", ref_molec%geom%angle(i,3), ")"    ,&
              calc, ref, dif
        dev = dev + (calc - ref)**2
    enddo
    rmsd = sqrt(dev/k)

    print'(X,A,X,F8.3,/)', "RMSD-angles (deg):", rmsd


    dev = 0.0
    k = 0
    if (debug) print*, "LIST OF DIHEDRALS"
    do i=1,ref_molec%geom%ndihed
        if (nonH) then
            if (adjustl(ref_molec%atom(ref_molec%geom%dihed(i,1))%name) == "H" .or. &
                adjustl(ref_molec%atom(ref_molec%geom%dihed(i,2))%name) == "H" .or. &
                adjustl(ref_molec%atom(ref_molec%geom%dihed(i,3))%name) == "H" .or. &
                adjustl(ref_molec%atom(ref_molec%geom%dihed(i,4))%name) == "H") cycle
        endif
        !Using an external counter in case nonH is used
        k=k+1
        ref  = calc_dihed(ref_molec%atom(ref_molec%geom%dihed(i,1)),&
                          ref_molec%atom(ref_molec%geom%dihed(i,2)),&
                          ref_molec%atom(ref_molec%geom%dihed(i,3)),&
                          ref_molec%atom(ref_molec%geom%dihed(i,4)))
        calc = calc_dihed(molec%atom(ref_molec%geom%dihed(i,1)),&
                          molec%atom(ref_molec%geom%dihed(i,2)),&
                          molec%atom(ref_molec%geom%dihed(i,3)),&
                          molec%atom(ref_molec%geom%dihed(i,4)))
        calc = calc*180.d0/PI
        ref  = ref*180.d0/PI
        dif = abs(calc - ref)
        dif = min(dif,abs(dif-360.))
        if (debug) &
        print'(3(A2,A1,I2,A5),A2,A1,I2,A1,X,2(F7.2,X),F11.6)', &
              ref_molec%atom(ref_molec%geom%dihed(i,1))%name, "(", ref_molec%geom%dihed(i,1), ") -- ",&
              ref_molec%atom(ref_molec%geom%dihed(i,2))%name, "(", ref_molec%geom%dihed(i,2), ") -- ",&
              ref_molec%atom(ref_molec%geom%dihed(i,3))%name, "(", ref_molec%geom%dihed(i,3), ") -- ",&
              ref_molec%atom(ref_molec%geom%dihed(i,4))%name, "(", ref_molec%geom%dihed(i,4), ")"    ,&
              calc, ref, dif
        dev = dev + (dif)**2
    enddo
    rmsd = sqrt(dev/k)

    print'(X,A,X,F8.3,/)', "RMSD-dihedrals (deg):", rmsd


    dev = 0.0
    k = 0
    if (debug) print*, "LIST OF ATOMS"
    do i=1,ref_molec%natoms
        if (nonH) then
            if (adjustl(ref_molec%atom(i)%name) == "H") cycle
        endif
        !Using an external counter in case nonH is used
        k=k+1
        dif = calc_dist(molec%atom(i),ref_molec%atom(i))
        if (debug) &
        print'(A2,A1,I2,A1,X,F11.6)', &
              ref_molec%atom(i)%name, "(", i, ")", dif
        dev = dev + (dif)**2
    enddo
    rmsd = sqrt(dev/k)
    print'(A,/)', '---------------------' 

    print'(X,A,X,F8.3,/)', "RMSD-struct (AA):", rmsd
   

    ! 9999. CHECK ERROR/NOTES
    ! -------------------------------------------------
    if (n_notes > 0) then 
        write(dummy_char,*) n_notes
        write(6,'(/,A,/)') "There were "//trim(adjustl(dummy_char))//" note(s) in this run"
    endif
    if (n_errors > 0) then 
        write(dummy_char,*) n_errors
        write(6,'(/,A,/)') "There were "//trim(adjustl(dummy_char))//" warning(s) in this run"
        call alert_msg("warning", "Files generated, but with too many warnings")
    endif


    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,reffile,ft,ft_ref,debug,nonH)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,ft_ref,reffile
        logical,intent(inout) :: debug, nonH
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

                case ("-r") 
                    call getarg(i+1, reffile)
                    argument_retrieved=.true.
                case ("-ftr") 
                    call getarg(i+1, ft_ref)
                    argument_retrieved=.true.

                case ("-dbg")
                    debug=.true.

                case ("-nonH")
                    nonH=.true.
        
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
        write(0,'(/,A)') '          R M S D - C A L C U L A T O R '    
        write(0,'(/,A)') '      Compare two structures using rmsd values '        
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(ft))
        write(0,*) '-r              ', trim(adjustl(reffile))
        write(0,*) '-ftr            ', trim(adjustl(ft_ref))
        write(0,*) '-dbg           ',  debug
        write(0,*) '-nonH          ',  nonH
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
            case("g96")
             call read_g96(unt,molec)
            case("xyz")
             call read_xyz(unt,molec)
            case("gro")
             call read_gro(unt,molec)
            case("pdb")
             call read_pdb_new(unt,molec)
            case("log")
             allocate(props)
             call parse_summary(unt,molec,props,"struct_only")
             deallocate(props)
            case("fchk")
             call read_fchk_geom(unt,molec)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

        call atname2element(molec)


        return


    end subroutine generic_strfile_read

end program calcRMS_bonded

