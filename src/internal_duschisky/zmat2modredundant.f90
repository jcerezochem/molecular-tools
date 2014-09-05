program internal_duschinski


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Version: internal_duschinski_v13
    !
    ! Description:
    ! -----------
    ! Program to analyse vibrations in term of internal coordinates.
    !
    ! Compilation instructions (for mymake script):
    !make$ echo "COMPILER: $FC"; sleep 1; $FC modules/alerts.f90 modules/structure_types_v3.2.f90 modules/line_preprocess.f90 modules/ff_build_module_v3.f90 modules/gro_manage_v2.f90 modules/pdb_manage_v2.f90 modules/constants_mod.f90 modules/atomic_geom_v2.f90 modules/gaussian_manage_v2.f90 modules/gaussian_fchk_manage_v2.f90 modules/symmetry_mod.f90 modules/MatrixMod.f90 modules/internal_SR_v8.4_corrected.f90 internal_duschinski_v13.f90 -llapack -o internal_duschinski_v13.exe -cpp -DDOUBLE
    !
    ! Change log:
    !
    ! TODO:
    ! ------
    !
    ! History
    ! V4: distributable version
    !     - only system types (molec) are used (needs V3 of structure_types and ff_build)
    !     - Allocation disabled
    !     - Clean up
    ! V5: detect symmetry of vibrations (only for Ci symmetry). Needs internal_SR_v6.f90
    ! V6: NDIM increased to 600 (recommended internal_SR_v7.f90, but also works with v6)
    !     Added quality checks
    ! V7: Proper treatment of redundants. Involves changes in SR (need version SR_v8)
    !     Use of symmetry addapted internal coordinates
    !     Added readZ optionç
    ! V8: Review and refine v7. Delete/clarify incongruent parts. Results not changed from v7 
    !     Add approxmiate orthogonal J=L1'^T L2'^T option and K=L1^T DeltaS
    ! V9: Add the possibility to extract the Hessian from gaussian log files
    ! V9.1: Add symm_file to specify a custom atom symmetry
    !       Added "split_line_back" to allow different relative PATHS to files
    ! V9.2: Add generation of state files (also for log files: no need to read mu)
    !       This version REQUIRES internal_SR_v8_corrected to solve a serious bug
    !...................................................................................................
    ! OUT OF TRACK:
    ! v10: Fixing Displacement in orthogonal internal coords (buggy version, only works partially)
    !  >> Ongoing works are left on an alternative track! (to be merged!)
    ! v11: Include derivatives for B to allow VH model in internal coordianates
    !      REQUIRES: internal_SR_v9
    !...................................................................................................
    ! v13: merge v9.2 with advances in v11 (B derivatives)
    !*********************************
    !
    ! v13_v4: addapt to molecular tools distribution
    !
    !============================================================================    

!*****************
!   MODULE LOAD
!*****************
!============================================
!   Generic (structure_types independent)
!============================================
    use alerts
    use line_preprocess
    use constants
!   Matrix manipulation (i.e. rotation matrices)
    use MatrixMod
!============================================
!   Structure types module
!============================================
    use structure_types
!============================================
!   Structure dependent modules
!============================================
    use gro_manage
    use pdb_manage
    use gaussian_manage
    use gaussian_fchk_manage
    use xyz_manage
!   Structural parameters
    use molecular_structure
    use ff_build
!   Bond/angle/dihed meassurement
    use atomic_geom
!   Symmetry support
    use symmetry_mod
!   For internal thingies
    use internal_module
    use zmat_manage

    implicit none

    integer,parameter :: NDIM = 400

    !====================== 
    !Options 
    logical :: nosym=.true.   ,&
               zmat=.true.    ,&
               tswitch=.false.,&
               symaddapt=.false.,&
               vertical=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Nred
    character(len=5) :: PG
    !Bonded info
    integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !AUXILIAR MATRICES
    real(8),dimension(NDIM,NDIM) :: Aux, Aux2
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !====================== 

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
    integer :: error, N
    !====================== 

    !====================== 
    !Auxiliars for LAPACK matrix nversion
    integer :: info
    integer,dimension(NDIM) :: ipiv
    real(8),dimension(NDIM,NDIM) :: work
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=16) :: dummy_char
    real(8) :: Theta, Theta2, Theta3
    !Read gaussian log auxiliars
    type(str_molprops),allocatable :: props
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, k90,k95,k99, nn, imin, imax,&
               i1,i2,i3,i4
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_ZMAT=11, &
               I_SYM=12,  &
               I_RED=13,  &
               O_DUS=20,  &
               O_DIS=21,  &
               O_DMAT=22, &
               O_DUS2=23, &
               O_DIS2=24, &
               O_STAT=25
    !files
    character(len=10) :: filetype="guess", ft
    character(len=200):: inpfile ="input.fchk",  &
                         inpfile2="input2.fchk", &
                         zmatfile="NO", &
                         symm_file="NO"
    !Control of stdout
    logical :: verbose=.false.
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

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype,nosym,zmat,verbose,zmatfile,symm_file)
 
    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Read structure
    ft=filetype
    call generic_strfile_read(I_INP,ft,state1)
    !Shortcuts
    Nat = state1%natoms
    Nvib = 3*Nat-6

    ! Get connectivity from the residue (needs to be in ANGS)
    call guess_connect(state1)
    if (nosym) then
        PG="C1"
    else if (trim(adjustl(symm_file)) /= "NO") then
        write(*,*) ""
        write(*,*) "Using custom symmetry file: "//trim(adjustl(symm_file)) 
        write(*,*) ""
        open(I_SYM,file=symm_file)
        do i=1,state1%natoms
            read(I_SYM,*) j, isym(j)
        enddo
        close(I_SYM)
        !Set PG to CUStom
        PG="CUS"
    else
        PG="XX"
        call symm_atoms(state1,isym)
        PG=state1%PG
    endif
    !From now on, we'll use atomic units
    state1%atom(:)%x = state1%atom(:)%x/BOHRtoANGS
    state1%atom(:)%y = state1%atom(:)%y/BOHRtoANGS
    state1%atom(:)%z = state1%atom(:)%z/BOHRtoANGS

    !Generate bonded info
    call gen_bonded(state1)
    !...Also get bond array (to be done in gen_bonded, would be cleaner)
    k = 0
    do i=1,state1%natoms-1
        do j=1,state1%atom(i)%nbonds
            if ( state1%atom(i)%connect(j) > i )  then
                k = k + 1
                state1%geom%bond(k,1) = i
                state1%geom%bond(k,2) = state1%atom(i)%connect(j)
            endif 
        enddo
    enddo
    state1%geom%nbonds = k
   

    !GEN BONDED SET FOR INTERNAL COORD
!     if (.not.zmat) then
!         print*, "Custom internal coordianates"
!         open(I_RED,file="modred.dat") 
!         call modredundant(I_RED,state1)
!         close(I_RED)
!         zmat=.false.
    if (zmat) then
        if (adjustl(zmatfile) == "NO") then
            call build_Z(state1,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
        else
            open(I_ZMAT,file=zmatfile,status="old")
            print*, "Z-matrix read from "//trim(adjustl(zmatfile))
            call read_Z(state1,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym,I_ZMAT)
            close(I_ZMAT)
            !Deactivate symaddapt (for the moment)
            PG = "C1"
        endif
        !Z-mat
        state1%geom%bond(1:Nat-1,1:2) = bond_s(2:Nat,1:2)
        state1%geom%angle(1:Nat-2,1:3) = angle_s(3:Nat,1:3)
        state1%geom%dihed(1:Nat-3,1:4) = dihed_s(4:Nat,1:4)
        state1%geom%nbonds  = Nat-1
        state1%geom%nangles = Nat-2
        state1%geom%ndihed  = Nat-3
    endif !otherwise all parameters are used

    !WRITE MODEREDUNDANT FILE
    open(I_RED,file="modred.dat") 
    do i=1,state1%geom%nbonds
        write(I_RED,'(A,100I10)') "B", state1%geom%bond(i,1:2)
    enddo
    do i=1,state1%geom%nangles
        write(I_RED,'(A,100I10)') "A", state1%geom%angle(i,1:3)
    enddo
    do i=1,state1%geom%ndihed
        write(I_RED,'(A,100I10)') "D", state1%geom%dihed(i,1:4)
    enddo
    close(I_RED)

    call cpu_time(tf)
    write(6,'(A,F12.3)') "CPU (s) for internal vib analysis: ", tf-ti

    stop


    !==============================================
    contains
    !=============================================
!         call parse_input(inpfile,filetype,nosym,zmat,verbose,zmatfile,symm_file)
    subroutine parse_input(inpfile,filetype,nosym,zmat,verbose,zmatfile,symm_file)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype,zmatfile,symm_file
        logical,intent(inout) :: nosym, verbose, zmat
        ! Localconsole in kate
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
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.

                case ("-nosym")
                    nosym=.true.
                case ("-sym")
                    nosym=.false.

                case ("-symfile")
                    nosym=.false.
                    call getarg(i+1, symm_file)
                    argument_retrieved=.true.

                case ("-readz") 
                    call getarg(i+1, zmatfile)
                    argument_retrieved=.true.

                case ("-zmat")
                    zmat=.true.
                case ("-nozmat")
                    zmat=.false.

                case ("-v")
                    verbose=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        ! Some checks on the input
        !----------------------------
        if (symaddapt.and.nosym) then
            print*, ""
            print*, "Symmetry addapted internal coordintes implies -sym. Turning on..."
            print*, ""
            nosym=.false.
        endif

       !Print options (to stderr)
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,'(/,A)') '        I N T E R N A L   A N A L Y S I S '    
        write(6,'(/,A)') '      Perform vibrational analysis based on  '
        write(6,'(/,A)') '            internal coordinates (D-V9.1)'        
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,*) '-f              ', trim(adjustl(inpfile))
        if (nosym) dummy_char="NO "
        if (.not.nosym) dummy_char="YES"
        write(6,*) '-[no]sym        ', dummy_char
        if (zmat) dummy_char="YES"
        if (.not.zmat) dummy_char="NO "
        write(6,*) '-symfile        ', trim(adjustl(symm_file))
        write(6,*) '-[no]zmat       ', dummy_char
        write(6,*) '-readz          ', trim(adjustl(zmatfile))
        write(6,*) '-v             ', verbose
        write(6,*) '-h             ',  need_help
        write(6,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(inout) :: filetype
        type(str_resmol),intent(inout) :: molec

        !local
        type(str_molprops) :: props

        if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line(inpfile,".",null,filetype)
        select case (adjustl(filetype))
            case("gro")
             call read_gro(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("pdb")
             call read_pdb_new(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("log")
             call parse_summary(I_INP,molec,props,"struct_only")
             call atname2element(molec)
             call assign_masses(molec)
            case("fchk")
             call read_fchk_geom(I_INP,molec)
             call atname2element(molec)
!              call assign_masses(molec) !read_fchk_geom includes the fchk masses
            case default
             call alert_msg("fatal","Trying to guess, but file type but not known: "//adjustl(trim(filetype))&
                        //". Try forcing the filetype with -ft flag (available: log, fchk)")
        end select

        else
        ! Predefined filetypes
        select case (adjustl(filetype))
            case("gro")
             call read_gro(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("pdb")
             call read_pdb_new(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("log")
             call parse_summary(I_INP,molec,props,"struct_only")
             call atname2element(molec)
             call assign_masses(molec)
            case("fchk")
             call read_fchk_geom(I_INP,molec)
             call atname2element(molec)
!              call assign_masses(molec) !read_fchk_geom includes the fchk masses
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select
        endif


        return


    end subroutine generic_strfile_read
       

end program internal_duschinski

