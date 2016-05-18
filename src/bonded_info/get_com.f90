program get_centers_mg


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Forms a gro structure with the one present in the fchk
    !
    ! Compilation instructions (for mymake script):
    !make$ echo "COMPILER: $FC"; sleep 1; $FC ../modules/alerts.f90 ../modules/structure_types_v2.f90 ../modules/line_preprocess.f90  ../modules/gro_manage_v2.f90 ../modules/pdb_manage_v2.f90 ../modules/constants_mod.f90 ../modules/gaussian_manage_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/xyz_manage.f90 ../modules/ff_build_module_v3.f90 fchk2gro_v2.1.f90 -o fchk2gro_v2.1.exe -cpp 
    !
    ! Change log:
    !
    ! TODO:
    ! ------
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
    use ff_build
    use metrics
    use atomic_geom
    use symmetry


    implicit none

    !====================== 
    !System variables
    type(str_resmol) :: molec, &
                        molec_aux
    !====================== 

    !=============
    !Counters and dummies
    integer :: i,j,k,l, jj,kk, iat
    real(8),dimension(3,3) :: T, X
    real(8),dimension(3)   :: Ip
    character(len=1) :: null
    logical :: overwrite=.false.,&
               make_connect=.false.
    !Swap related counters
    integer :: iat_new, iat_orig, nswap
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10, &
               I_SWP=11, &
               O_OUT=20  
    !files
    character(len=5) :: resname="read"
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile="input.fchk"
    logical :: do_inertia=.false.
    !status
    integer :: IOstatus
    character(len=7) :: stat="new" !do not overwrite when writting
    !===================


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype,do_inertia)

 
    ! 1. READ INPUT
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(filetype) == "guess") call split_line_back(inpfile,".",null,filetype)
    call generic_strmol_reader(I_INP,filetype,molec)
    close(I_INP)

    ! 2. MAKE THE THING
    ! ------------------------------
    call get_cog(molec)
    call get_com(molec)

    print*, "CENTER OF GEOMETRY"
    print*, molec%cogX, molec%cogY, molec%cogZ
    print*, "-"
    print*, "CENTER OF MASS"
    print*, molec%comX, molec%comY, molec%comZ

    if (do_inertia) then
        call translate_molec(molec,(/-molec%comX,-molec%comY,-molec%comZ/))
        ! Compute inertia axes
        call inertia(molec,T)
        call diagonalize_full(T,3,X,Ip,"lapack")
        
        ! Principal moment of inertia
        Ip(1) = Ip(1) * ANGStoM**2 * AMUtoKG
        Ip(2) = Ip(2) * ANGStoM**2 * AMUtoKG
        Ip(3) = Ip(3) * ANGStoM**2 * AMUtoKG
        ! Rotational constants
        Ip(1) = plank/8.d0/pi**2/Ip(1)
        Ip(2) = plank/8.d0/pi**2/Ip(2)
        Ip(3) = plank/8.d0/pi**2/Ip(3)
        
        print'(2/,X,A)', "ROTATIONAN CONSNTANTS (GHz)"
        print'(3F12.5)', Ip(1)*1.d-9, Ip(2)*1.d-9, Ip(3)*1.d-9

        call MAT0(6,X,3,3,"AXIS OF INERTIA")

    endif

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,do_inertia)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype
        logical,intent(inout)          :: do_inertia
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

                case ("-iner")
                    do_inertia=.true.
                case ("-noiner")
                    do_inertia=.false.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        ! Some checks on the input
        !----------------------------

       !Print options (to stderr)
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '            G  E T    C E N T E R S '      
        write(0,'(/,A)') '     Computes COM, COG and, optionally, Inertia'     
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input file                       ', trim(adjustl(inpfile))
        write(0,*)       '-ft          \_ FileTyep                      ', trim(adjustl(filetype))
        write(0,*)       '-[no]iner    Compute inertia axis            ',  do_inertia
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )


        return
    end subroutine parse_input


end program get_centers_mg

