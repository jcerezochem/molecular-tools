program get_reordering


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to analyse vibrations in term of internal coordinates.
    !
    ! Compilation instructions (for mymake script):
    !make$ echo "COMPILER: $FC"; sleep 1; $FC ../modules/alerts.f90 ../modules/structure_types_v3.f90 ../modules/line_preprocess.f90 ../modules/ff_build_module_v3.f90 ../modules/gro_manage_v2.f90 ../modules/pdb_manage_v2.f90 ../modules/constants_mod.f90 ../modules/atomic_geom_v2.f90 ../modules/gaussian_manage_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/symmetry_mod.f90 ../modules/MatrixMod.f90 internal_SR_v6.f90 internal_duschinski_v5.f90 -llapack -o internal_duschinski_v5.exe -cpp -DDOUBLE
    !
    ! Change log:
    !
    ! TODO:
    ! ------
    !
    ! History
    ! V3: Internal analysis is based on internal_duschinski_v5 (never finished...)
    ! V4: Internal analysis is based on internal_duschinski_v7
    !  V4b (not in the main streamline!): includes the generation of scans calc. for Gaussian.
    !  V4c: the same as 4b. Bug fixes on guess_connect (ff_build module_v3 was buggy)
    !
    !Addapted to v4 release (distribution upgrade). Feb '14
    !v4 releases:
    !v4.0.1:
    ! -use redundant coordinates
    !v4.0.1.1:
    ! - include UnSym (MOLCAS) file as input (freqs and nm)
    !v4.0.1.2:
    ! - if the calculation is only a internal scan, do not need the hessian (so do not try to read it)
    !
    !============================================================================    


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
    use atomic_geom
    use symmetry


    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !System variables
    type(str_resmol)  :: molecule, molecule2
    integer,dimension(1:NDIM) :: isym
    integer :: Nat
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: error
    character(1) :: null
    real(8) :: dist_thr=0.1d0, dist
    !====================== 

    !=============
    !Counters
    integer :: i,j,k
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10, &
               I_REF=11
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile ="input.fchk"
    character(len=10) :: reffiletype="guess"
    character(len=200):: reffile ="ref.fchk"
    !status
    integer :: IOstatus
    !===================

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

! (End of variables declaration) 
!==================================================================================
    call cpu_time(ti)

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(molecule)
    call allocate_atoms(molecule2)
    !===========================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype,reffile,reffiletype)
 
    !================
    ! READ DATA
    !================
    if (adjustl(filetype) == "guess") &
    call split_line_back(inpfile,".",null,filetype)
    if (adjustl(reffiletype) == "guess") &
    call split_line_back(reffile,".",null,reffiletype)
    
    ! STRUCTURE FILE
    print'(X,A)', "READING STRUCTURE..."
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,filetype,molecule)
    close(I_INP)
    
    print'(X,A)', "READING REFERENCE STRUCTURE..."
    open(I_REF,file=reffile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_REF,reffiletype,molecule2)
    close(I_REF)
    
    !Shortcuts
    Nat = molecule%natoms
    print'(X,A,/)', "Done"

    print*, "Structure  --> Reference"
    do i=1,Nat 
        do j=1,Nat 
            dist = calc_atm_dist(molecule%atom(i), molecule2%atom(j))
            if (dist < dist_thr) then
                print*, i,j,'  --  ', molecule%atom(i)%name, molecule2%atom(j)%name
            endif
        enddo
    enddo



    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,reffile,reffiletype)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype, reffile, reffiletype
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
                case ("-f") 
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.
                case ("-ft") 
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.  
                case ("-r") 
                    call getarg(i+1, reffile)
                    argument_retrieved=.true.
                case ("-frt") 
                    call getarg(i+1, reffiletype)
                    argument_retrieved=.true. 
                case ("-thr") 
                    call getarg(i+1, sym_thr)
                    argument_retrieved=.true.
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       !Print options (to stdx)
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '             G E T    R E O R D E R I N G '    
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input file                       ', trim(adjustl(inpfile))
        write(0,*)       '-ft          \_ FileTyep                      ', trim(adjustl(filetype))
        write(0,*)       '-r           Input file                       ', trim(adjustl(reffile))
        write(0,*)       '-frt         \_ FileTyep                      ', trim(adjustl(reffiletype))
        write(0,*)       '-thr         Symm. thr [tight|normal|loose]   ', trim(adjustl(sym_thr))
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program get_reordering

