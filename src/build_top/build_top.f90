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
    !============================================
    !  Internal thingies
    !============================================
    use zmat_manage 

    implicit none

!     common /GLOBAL_OPTS/ do_refine_charges, I_DB, rename_atoms

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
    integer :: nat_alloc=-1
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
    call parse_input(inpfile,filetype,outgro,outtop,q_type,resname,database,united_atom,nat_alloc)

    ! Allocate atoms
    if (nat_alloc<0) then
        call allocate_atoms(molec)
        do i=1,10
            call allocate_atoms(residue(i))
            call allocate_atoms(residueUA(i))
        enddo
    else 
        call allocate_atoms(molec,nat_alloc)
        do i=1,10
            call allocate_atoms(residue(i),nat_alloc)
            call allocate_atoms(residueUA(i),nat_alloc)
        enddo
    endif
    
    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Load atomic coordinates on "molecule"
    if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line(inpfile,".",null,filetype)
    endif
    
    call generic_strmol_reader(I_INP,filetype,molec)

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
             call write_top(O_TOP,residue(1))
            case("top")
             call write_top(O_TOP,residue(1))
            case("rtp")
             call write_rtp(O_TOP,residue(1))
            case default
             !By default a top format is use
             call write_top(O_TOP,residue(1))
             call alert_msg("note","Topology file written in *.top format")
        end select
        close(O_TOP)
    endif


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

    subroutine parse_input(inpfile,filetype,outfile,outfile2,q_type,resname,database,united_atom,nat_alloc)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype,outfile,outfile2,resname,q_type,database
        logical,intent(inout) :: united_atom
        integer :: nat_alloc
        ! Local
        character(len=500) :: input_command
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        
        call getarg(0,input_command)
        !Get input flags
        do i=1,iargc()
            call getarg(i,arg)
            input_command = trim(adjustl(input_command))//" "//trim(adjustl(arg))
        enddo

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
        
                case ("-db")
                    call getarg(i+1, database)
                    argument_retrieved=.true.

                case ("-ua")
                     united_atom = .true.
                     
                case ("-alloc-atm")
                    call getarg(i+1,arg)
                    read(arg,*) nat_alloc
                    argument_retrieved=.true.
                     

                case ("-h")
                    need_help=.true.
                    
                    
                ! Control verbosity
                case ("-quiet")
                    verbose=0
                case ("-concise")
                    verbose=1
                case ("-v")
                    verbose=2
                case ("-vv")
                    verbose=3
                    

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
! 
!         if ( q_type == "all_zero" ) do_refine_charges=.false.
          

       !Print options (to stdout)    
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '              B U I L D _ T O P '    
        write(6,'(/,A)') '     A connectivity based topology builder '        
        write(6,'(A)')   '                for GROMACS'     
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag         Description                   Value'
        write(6,'(A)')   '-------------------------------------------------------------------'


        write(0,*) '-h             ',  need_help
        
        write(6,*) '-f           Input file (structure)        ', trim(adjustl(inpfile))
        write(6,*) '-ft          \_ FileType                   ', trim(adjustl(filetype))
        write(6,*) '-o           Output (com) file             ', trim(adjustl(outgro))
        write(6,*) '-p           Topology file                 ', trim(adjustl(outtop))
        write(6,*) '-chrg        Charge type (mulliken...)     ', trim(adjustl(q_type))
        
        write(6,*) '-res         Reside name                   ', trim(adjustl(resname))
        write(6,*) '-db          Database to set atom types    ', trim(adjustl(database))
        write(6,'(X,A,I0)') &
                   '-alloc-atm   Atoms to allocate(-1:default) ', nat_alloc
        write(6,*) '-h           Show this help and quit       ',  need_help
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(A)') 'Input command:'
        write(6,'(A)') trim(adjustl(input_command))   
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                       'Verbose level:  ', verbose        
        write(6,'(A)') '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input

end program build_top

