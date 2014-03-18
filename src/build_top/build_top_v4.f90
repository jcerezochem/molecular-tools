program build_top

    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to generate a topology file from gaussian log files.
    ! Both coordinates and atomic charges are read. Atom types are
    ! guessed from the (also guessed) connectivity [but it works...]
    !
    ! Compilation instructions (for mymake script):
    !make$ gfortran ../modules/alerts.f90 ../modules/structure_types_v2_ALLOC.f90 ../modules/allocation_mod.f90 ../modules/line_preprocess.f90 ../modules/ff_build_module_v2.f90 ../modules/gro_manage_v2.f90 ../modules/gaussian_manage_v2.f90 ../modules/pdb_manage_v2.f90 ../modules/molpro_manage_v2.f90 BUILD_TOP_v2.f90 -o BUILD_TOP_v2.exe -cpp
    !
    ! Change log:
    ! Feb 2012: ff_build and related subroutine now work on residue not on system (beta)
    ! May 2012: Uses version 2 of molecular tools (v2)
    !           Implement allocation_mod
    ! Jun 2012 (v2) Some features are broken (un-reverted tests?) -- SOLVED
    ! Version 4: Uses v4 modules (allocation disabled)
    !============================================================================    

    use structure_types
    use line_preprocess
    use ff_build
    use molecular_structure
    use gro_manage
    use gaussian_manage
    use pdb_manage
    use molpro_manage
    use alerts

    implicit none

    common /ALERT_COUNT/ n_notes, n_errors
    common /GLOBAL_OPTS/ do_refine_charges, I_DB, rename_atoms

    !====================== 
    !Number of error/notes
    integer :: n_notes=0, n_errors=0
    !====================== 

    !====================== 
    !Options 
    logical :: call_vmd=.false.,          & 
               do_refine_charges=.false., & !.true.
               united_atom=.false.,       &
               rename_atoms=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molec !, molecUA
    type(str_resmol),dimension(10) :: residue, residueUA
    character(len=6) :: resname="UNK"
    character(len=8) :: q_type="mulliken"
    real :: charge
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: id, ih
    character(len=1) :: null
    character(len=4) :: atname
    character(len=16) :: dummy_char
    integer :: dummy_int
    !====================== 

    !=============
    !Counters
    integer :: i,j
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_MUL=11,  & 
               I_DB=0,    & !zero means no external DB
               S_SPDB=30, &
               O_GRO=20,  &
               O_TOP=21
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile="input.log",          &
                         outgro="out.gro",             &
                         outtop="out.top",             &
                         database="default"               
    !status
    integer :: IOstatus
    !===================


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype,outgro,outtop,q_type,resname,call_vmd,database,united_atom)

    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Load atomic coordinates on "molecule"
    if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line(inpfile,".",null,filetype)
    endif

    call generic_strfile_read(I_INP,filetype,molec)
         
    close(I_INP)


    ! 2. ASSOCIATE ATOM TYPES AND BONDED TERMS (if not disabled)
    ! -----------------------------------------
    if ( adjustl(outtop) /= "none" )  then
        !Attypes DataBase
        if ( adjustl(database) /= "default" ) then
            I_DB=12
            open(I_DB, file=database, status="old", IOSTAT=IOstatus)
            if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(database)) )
        endif

        call sist2res(molec,residue)
!         call sist_nres(molec)

        !Only one residue supported (for the momenfor the moment)
        if (molec%nres /= 1) call alert_msg("fatal","Only structures with one residue are currently supported")
        call guess_connect(residue(1))
        call build_ff(residue(1))
        call gen_bonded(residue(1))

        if (I_DB /= 0) close(I_DB)
    endif

    !United atoms model:
    if (united_atom) then
        call aa2ua(residue(1),residueUA(1))
        residue(1) = residueUA(1)
        !Recalculate connections 
        call guess_connect(residue(1))
        call gen_bonded(residue(1))
        call res2sist(residue,molec%nres,molec%BoxX,molec%BoxY,molec%BoxZ,molec)
    endif


    ! 3. WRITE COORDINATES (.gro) AND TOPOLOGY (.top)
    ! -------------------------------------------------
    !GROFILE:
    open(O_GRO,file=outgro,status='replace')
    call write_gro(O_GRO,molec)
    close(O_GRO)

    !TOPFILE (if not disabled)
    if ( trim(adjustl(outtop)) /= "none" ) then
        !The topology builder uses residue as input:
!         call sist2res(molec,residue)
        open(O_TOP,file=outtop,status='replace')
        call split_line(outtop,".",null,filetype)
        select case (adjustl(filetype))
            case("itp")
             call write_top_new(O_TOP,residue(1))
            case("top")
             call write_top_new(O_TOP,residue(1))
            case("rtp")
             call write_rtp(O_TOP,residue(1))
            case default
             !By default a top format is use
             call write_top_new(O_TOP,residue(1))
             call alert_msg("note","Topology file written in *.top format")
        end select
        close(O_TOP)
    endif


    ! 4. REPRESENT THE MOLECULE (using VMD)
    !----------------------------------------------
    if (call_vmd) &
        call system( "vmd "//outgro )


    ! 9999. CHECK ERROR/NOTES
    ! -------------------------------------------------
    if (n_notes > 0) then 
        write(dummy_char,*) n_notes
        write(6,'(/,A,/)') "There were "//trim(adjustl(dummy_char))//" note(s) in this run"
    endif
    if (n_errors > 0) then 
        write(dummy_char,*) n_errors
        write(6,'(/,A,/)') "There were "//trim(adjustl(dummy_char))//" warning(s) in this run"
        call alert_msg("fatal", "Files generated, but with too many warnings")
    endif

    ! Summary output files
    if ( trim(adjustl(outtop)) /= "none" ) then
        write(6,'(/,A,/)') "2 output files generated: "//trim(adjustl(outgro))//",  "//trim(adjustl(outtop))
    else
        write(6,'(/,A,/)') "1 output file generated: "//trim(adjustl(outgro))
    endif


    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,outfile,outfile2,q_type,resname,call_vmd,database,united_atom)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        common /GLOBAL_OPTS/ do_refine_charges

        character(len=*),intent(inout) :: inpfile,filetype,outfile,outfile2,resname,q_type,database
        logical,intent(inout) :: call_vmd, united_atom
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

                case ("-o") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.

                case ("-p") 
                    call getarg(i+1, outfile2)
                    argument_retrieved=.true.

                case ("-chrg") 
                    call getarg(i+1, q_type)
                    argument_retrieved=.true.

                case ("-res") 
                    call getarg(i+1, resname)
                    argument_retrieved=.true.

                case ("-vmd")
                    call_vmd=.true.
        
                case ("-db")
                    call getarg(i+1, database)
                    argument_retrieved=.true.

                case ("-ua")
                     united_atom = .true.

                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        ! Some checks on the input
        !----------------------------

        !-p disabling accepts different inputs
        if ( trim(adjustl(outtop)) == "none" .or.   &
             trim(adjustl(outtop)) == "NONE" .or.   &
             trim(adjustl(outtop)) == "no"   .or.   &
             trim(adjustl(outtop)) == "NO"        ) then 
            outtop = "none"
        endif

        !-chrg options
        if ( q_type /= "mulliken" .and.  &
             q_type /= "elec_pot" .and.   &
             q_type /= "all_zero" ) then
            call alert_msg("fatal", 'Options for charges (-chrg) are "mulliken", "elec_pot" and "all_zero", select one of these')
        endif

        if ( q_type == "all_zero" ) do_refine_charges=.false.
          

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '              B U I L D _ T O P '    
        write(0,'(/,A)') '     A connectivity based topology builder '        
        write(0,'(A)')   '                for GROMACS'    
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(filetype))
        write(0,*) '-o              ', trim(adjustl(outgro))
        write(0,*) '-p              ', trim(adjustl(outtop))
        write(0,*) '-chrg           ', trim(adjustl(q_type))
        write(0,*) '-res            ', trim(adjustl(resname))
        write(0,*) '-db             ', trim(adjustl(database))
        write(0,*) '-vmd           ',  call_vmd
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
            case("molpro")
             call get_log_xyz(I_INP,molec)
             call read_charge(I_INP,molec,"all_zero")
             call alert_msg("note","No charges available in a MOLRPO log file. Setting all of them to zero")
             rename_atoms = .true.
            case("log")
             allocate(props)
             call parse_summary(I_INP,molec,props,"struct_only")
             deallocate(props)
             call read_charge(I_INP,molec,q_type)
             rename_atoms = .true.
            case("gro")
             call read_gro(I_INP,molec)
             call read_charge(I_INP,molec,"all_zero")
             call alert_msg("note","No charges available in a GRO file. Setting all of them to zero")
            case("pdb")
             call read_pdb_new(I_INP,molec)
             call read_charge(I_INP,molec,"all_zero")
             call alert_msg("note","No charges available in a PDB file. Setting all of them to zero")
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

        call atname2element(molec)

        return

    end subroutine generic_strfile_read

end program build_top

