program read_opt


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
    use ff_build
    use atomic_geom
    use symmetry
    !============================================
    !  Vibrational info
    !============================================
    use vibrational_analysis, only: Hlt_to_Hess,Hess_to_Hlt


    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: overwrite=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol)  :: molecule, molec_aux
    character(len=50) :: calc_type,basis
    character(len=60) :: method
    integer           :: charge, mult
    integer,dimension(1:NDIM) :: iord
    integer :: Nat, Nvib
    character(len=5) :: PG
    ! Energies
    real(8) :: E_scf, E_td, E_tot
    !======================

    !=====================
    ! Available data swithes
    logical :: have_gradient, have_hessian,&
               ! Energies
               have_SCF,have_TD,have_TOT
    !===================== 

    !====================== 
    !MATRICES
    !Other matrices
    real(8),dimension(:),allocatable :: R, E
    !====================== 

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
    integer :: error, N, lenght
    !====================== 

    ! IRC vars
    integer :: ninfo, nsteps

    !====================== 
    !Auxiliars for LAPACK matrix nversion
    integer :: info
    integer,dimension(NDIM) :: ipiv
    real(8),dimension(NDIM,NDIM) :: work
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=50) :: dummy_char
    character(len=240):: line
    integer :: iat_new, iat_orig, j_orig, j_new, nswap
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, nn, imin, imax, iii, i_scan
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               O_STR=20,  &
               O_GAU=21
    !files
    character(len=10) :: filetype="guess", filetype_out="guess"
    character(len=200):: inpfile ="input.fchk", &
                         outfile="default"
    character(len=200) :: basefile, title
    character(len=10)  :: extension, label
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
    call parse_input(inpfile,filetype,outfile,filetype_out,overwrite)

    ! Output names are labeled with steps
    call split_line_back(outfile,".",basefile,extension)
 
    !================
    ! READ DATA
    !================
    if (adjustl(filetype) == "guess") &
    call split_line_back(inpfile,".",null,filetype)
    if (filetype /= "fchk") call alert_msg("fatal","Only FCHK are currently supported")
    if (adjustl(filetype_out) == "guess") &
    call split_line_back(outfile,".",null,filetype_out)

    ! STRUCTURE FILE
    print'(X,A)', "READING STRUCTURE..."
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,filetype,molecule)
    !Shortcuts
    Nat = molecule%natoms
    Nvib = 3*Nat-6
    print'(X,A,/)', "Done"

    ! JOB INFO
    print'(X,A)', "READING JOB INFO..."
    rewind(I_INP)
    call read_gauss_job(I_INP,filetype,calc_type,method,basis)
    print'(X,A)', " Job type: "//trim(adjustl(calc_type))
    print'(X,A)', " Method  : "//trim(adjustl(method))
    print'(X,A)', " Basis   : "//trim(adjustl(basis))
    rewind(I_INP)
    call read_gauss_chargemult(I_INP,filetype,charge,mult)
    print'(X,A,I0)', " Charge  : ", charge
    print'(X,A,I0)', " Mult.   : ", mult
    print'(X,A,/)', "Done"

    ! Opt iNFO
    !-----------
    ! Data info
    call read_fchk(I_INP,"Optimization Num results per geometry",dtype,N,A,IA,error)
    if (error /= 0) call alert_msg("fatal","Looking for Scan info")
    ninfo=IA(1)
    deallocate(IA)
    ! 3Nat
    call read_fchk(I_INP,"Optimization Num geometry variables",dtype,N,A,IA,error)
    if (error /= 0) call alert_msg("fatal","Looking for Scan info")
    if (IA(1) /= 3*Nat) call alert_msg("fatal","Atoms in Scan are not consisntent")
    deallocate(IA)
    ! Run over the only Opt step (Kown number of steps)
    i_scan = 0
    do ii=1,1
        write(line,'(A9,I8,X,A22)') "Opt point",i_scan+1,"Results for each geome"
        ! Energies and Distances
        call read_fchk(I_INP,adjustl(line),dtype,N,A,IA,error)
        if (error /= 0) exit
        i_scan = i_scan + 1
        nsteps = N/ninfo
        ! Read all points
        allocate(E(1:nsteps),R(1:nsteps))
        j=0
        do i=1,nsteps
            j=j+1
            E(i) = A(j)
            j=j+1
            R(i) = A(j)
        enddo
        deallocate(A)
        
        ! Read all geoms and write to files
        write(line,'(A9,I8,X,A10)') "Opt point",i_scan,"Geometries"
        call read_fchk(I_INP,adjustl(line),dtype,N,A,IA,error)
        if (error /= 0) call alert_msg("fatal","Looking for Scan info (Geometries)")
        if (nsteps /= N/(3*Nat)) call alert_msg("fatal","Inconsistency in the number of steps")

        do i=1,nsteps
            k=(i-1)*Nat*3
            do j=1,Nat 
                k=k+1
                molecule%atom(j)%x = A(k)*BOHRtoANGS
                k=k+1
                molecule%atom(j)%y = A(k)*BOHRtoANGS
                k=k+1
                molecule%atom(j)%z = A(k)*BOHRtoANGS
            enddo
            label = int20char(i,2)
            write(title,'(A,I0,5X,A,F15.6,X,A,F10.4)') "Opt step ", i, "E=", E(i)
            outfile=trim(adjustl(basefile))//trim(adjustl(label))//"."//trim(adjustl(extension))
            open(O_STR,file=outfile)
            call generic_strmol_writer(O_STR,filetype_out,molecule,title=title)
            close(O_STR)
            ! Write g09 input
            outfile=trim(adjustl(basefile))//trim(adjustl(label))//".com"
            open(O_GAU,file=outfile)
            call write_gcom(O_GAU,molecule,&
                                  !Optional args
                                  chkname=outfile,&! 
                                  calc=calc_type,  &! e.g. SP, Freq, Opt...
                                  method=method,   &!
                                  basis=basis,     &!
                                  title=title      )
            close(O_GAU)
        enddo
        
        deallocate(A)
    enddo

    write(0,*) 'Number of opt points', i_scan

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,outfile,filetype_out,overwrite)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype,outfile,filetype_out
        logical,intent(inout)          :: overwrite
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
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.  
                case("-o")
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                case ("-fto") 
                    call getarg(i+1, filetype_out)
                    argument_retrieved=.true.  
                case("-ow")
                    overwrite=.true.
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        !The default output is now xyz
        if (adjustl(outfile) == "default") then
            call split_line_back(inpfile,".",outfile,null)
            outfile=trim(adjustl(outfile))//".xyz"
        endif

       !Print options (to stdx)
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '             R E A D    O P T '    
        write(0,'(/,A)') '       Read Optimization steps from FCHK' 
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input file                       ', trim(adjustl(inpfile))
        write(0,*)       '-ft          \_ FileTyep                      ', trim(adjustl(filetype))
        write(0,*)       '             (the molecule is rotated FIRST)'
        write(0,*)       '-o           Output file                      ', trim(adjustl(outfile))
        write(0,*)       '-fto         \_ FileTyep                      ', trim(adjustl(filetype_out))
        write(0,*)       '-ow          Force overwrite output          ',  overwrite
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program read_opt

