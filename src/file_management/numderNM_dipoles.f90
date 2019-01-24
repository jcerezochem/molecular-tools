program numderNM_dipoles


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
    character(len=50) :: calc_type,basis
    character(len=60) :: method
    integer           :: charge, mult
    integer,dimension(1:NDIM) :: iord
    integer :: Nat, Nvib
    character(len=5) :: PG
    ! Energies
    real(8) :: E_scf, E_td, E_tot
    ! FCHK data
    real(8),dimension(1:NDIM) :: DipFW, DipBW, Dip0
    real(8),dimension(1)      :: DD ! fake
    !
    integer :: iflag1, iflag2
    logical :: is_ddip
    ! DELTA
    real(8) :: delta, delta2
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
    real(8) :: R, E
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
    integer :: I_INP =10,  &
               O_STR =20,  &
               O_GAU =21,  &
               O_FCHK=22,  &
               O_OUT =23
    !files
    character(len=10)  :: filetype="guess", filetype_out="guess"
    character(len=200) :: inpfile ="input.fchk", &
                          outfile="default"
    character(len=200) :: basefile="basename", title, path
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
    call parse_input(basefile,filetype,outfile,filetype_out,overwrite)
 
    !================
    ! READ DATA
    !================
    ! We assume that a base file called basename.fchk with TD SP (at minimum) exists
    ! And take trdips from this fileinp
    call split_line_back(basefile,'.fchk',line,null)
    basefile=adjustl(line)
    inpfile=trim(adjustl(basefile))//".fchk"
    if (adjustl(filetype) == "guess") &
    call split_line_back(inpfile,".",null,filetype)
    if (filetype /= "fchk") call alert_msg("fatal","Only FCHK are currently supported")
    ! STRUCTURE FILE
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    ! Only Nat is needed
    call generic_natoms_reader(I_INP,filetype,Nat)
    Nvib=3*Nat-6
    ! JOB INFO
    print'(X,A)', "READING JOB INFO..."
    call read_gauss_job(I_INP,filetype,calc_type,method,basis)
    if (adjustl(calc_type) /= "SP") then
        print'(X,A)', " Job type: "//trim(adjustl(calc_type))
        call alert_msg('warning','The calculation type is not Single Point')
    endif
    ! TRDIPS (at ref geom)
    iflag1=-1
    iflag2=-1
    is_ddip=.false.
    call read_gaussfchk_dip(I_INP,iflag1,iflag2,is_ddip,'eldip',Dip0,DD,error)
    close(I_INP)
    
    ! Compute derivatives
    allocate(A(1:9*Nat))
    is_ddip=.false.
    k=0
    do i=1,Nvib
        ! ** Q **
        ! Delta changes for each mode
        delta=0.d0
        ! Bw
        iflag1=-1
        iflag2=-1
        write(inpfile,'(A,I0,A)') "Mode", i, "_Cart_numder_bw.fchk"
        open(I_INP,file=inpfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
        call read_gaussfchk_dip(I_INP,iflag1,iflag2,is_ddip,'eldip',DipBW,DD,error)
        call read_disp(I_INP,delta)
        delta2 = -delta
        close(I_INP)
        ! Fw
        iflag1=-1
        iflag2=-1
        write(inpfile,'(A,I0,A)') "Mode", i, "_Cart_numder_fw.fchk"
        open(I_INP,file=inpfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
        call read_gaussfchk_dip(I_INP,iflag1,iflag2,is_ddip,'eldip',DipFW,DD,error)
        rewind(I_INP)
        call read_disp(I_INP,delta)
        delta2 = delta2+delta
        close(I_INP)
        ! Compute derivative (in AU)
        k=k+1 !mu_x
        A(k) = (dipFW(1)-dipBW(1))/delta2
        k=k+1 !mu_y
        A(k) = (dipFW(2)-dipBW(2))/delta2
        k=k+1 !mu_z
        A(k) = (dipFW(3)-dipBW(3))/delta2
    
    enddo    
    
    ! Open out-eldip file
    outfile='eldip_nm_'//trim(adjustl(basefile))
    open(O_OUT,file=outfile,status='unknown',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(outfile)) )
    
    write(O_OUT,'(3(X,E18.9))') Dip0(1:3)
    write(O_OUT,'(3(X,E18.9))') Dip0(1:3)
    do i=1,Nvib
        ii=3*(i-1)
        write(O_OUT,'(3(X,E18.9))') A(ii+1:ii+3)
    enddo
    close(O_OUT)
    print'(/,X,A)', "eldip written to file: "//trim(adjustl(outfile))
    
    

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
        write(0,'(/,A)') '             P R E P A R E   N U M D E R  C O M S '    
        write(0,'(/,A)') '       Prepare input files (.com) for numerical ' 
        write(0,'(A)')   '       derivatives displacing the initila structure'
        write(0,'(A)')   '               +/- 0.001 Angs for each atom'
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input file                       ', trim(adjustl(inpfile))
        write(0,*)       '-ft          \_ FileTyep                      ', trim(adjustl(filetype))
        write(0,*)       '-o           Output file                      ', trim(adjustl(outfile))
        write(0,*)       '-fto         \_ FileTyep                      ', trim(adjustl(filetype_out))
        write(0,*)       '-ow          Force overwrite output          ',  overwrite
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       
    subroutine read_disp(unt,disp)
    
        integer,intent(in) :: unt
        real(8),intent(out):: disp
        
        !Local
        character(len=200) :: line, values
        character          :: null
        integer :: ioflag
        
        rewind(unt)
        do
            read(unt,'(A)',iostat=ioflag) line
            if (ioflag/=0) call alert_msg('fatal','Searching disp')
            if (index(line,"Animation of normal mode") /= 0) then
                call split_line(line,"=",null,values)
                read(values,*) disp
                exit
            endif
            
        enddo
        
        return
        
    end subroutine read_disp

end program numderNM_dipoles

