program read_scan


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to extract individuated information from the FCHK of
    ! a Scan job. For each scan point, it generates:
    !  * xyz file with the last geom in the scan point (contraint opt)
    !  * com file with an input for further jobs on each point
    !  * fchk file with geom, energy and gradient (at each contraint opt point)
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
    type(str_resmol)  :: molecule, molec_aux
    character(len=50) :: calc_type,basis
    character(len=60) :: method
    integer           :: charge, mult
    integer,dimension(1:NDIM) :: iord
    integer :: Nat, Nvib
    character(len=5) :: PG
    ! Energies
    real(8) :: E_scf, E_td, E_tot
    ! FCHK data
    real(8),dimension(1:NDIM) :: Grad
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
               O_GAU=21,  &
               O_FCHK=22
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

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(molecule)
    call allocate_atoms(molec_aux)
    !===========================

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

    ! Scan iNFO
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
    ! Run over all scan geoms (Unkown number of steps)
    ! Allocate for at most 100 points
    allocate(E(1:100),R(1:100))
    i_scan = 0
    do 
        write(line,'(A9,I8,X,A22)') "Opt point",i_scan+1,"Results for each geome"
        ! Energies and Distances
        call read_fchk(I_INP,adjustl(line),dtype,N,A,IA,error)
        if (error /= 0) exit
        i_scan = i_scan + 1
        nsteps = N/ninfo
        ! Read last point
        E(i_scan) = A(N-1)
        R(i_scan) = A(N)
        deallocate(A)
        
        ! Read all geoms at the scan point, and get the last point (constraint minimum)
        write(line,'(A9,I8,X,A10)') "Opt point",i_scan,"Geometries"
        call read_fchk(I_INP,adjustl(line),dtype,N,A,IA,error)
        if (error /= 0) call alert_msg("fatal","Looking for Scan info (Geometries)")
        if (nsteps /= N/(3*Nat)) call alert_msg("fatal","Inconsistency in the number of steps")

        k=(nsteps-1)*Nat*3
        do j=1,Nat 
            k=k+1
            molecule%atom(j)%x = A(k)*BOHRtoANGS
            k=k+1
            molecule%atom(j)%y = A(k)*BOHRtoANGS
            k=k+1
            molecule%atom(j)%z = A(k)*BOHRtoANGS
        enddo
        
        ! Read all Grads at the scan point, and get the last point (constraint minimum)
        write(line,'(A9,I8,X,A22)') "Opt point",i_scan,"Gradient at each geome"
        call read_fchk(I_INP,adjustl(line),dtype,N,A,IA,error)
        if (error /= 0) call alert_msg("fatal","Looking for Scan info (Gradient)")
        if (nsteps /= N/(3*Nat)) call alert_msg("fatal","Inconsistency in the number of steps")

        k=(nsteps-1)*Nat*3
        do j=1,3*Nat
            Grad(j) = A(j)
        enddo
        deallocate(A)
        
        ! Print to files
        label = int20char(i_scan,2)
        write(title,'(A,I0,5X,A,F15.6,X,A,F10.4)') "Scan step ", i_scan, "E=", E(i_scan)
        outfile=trim(adjustl(basefile))//"_"//trim(adjustl(label))//"."//trim(adjustl(extension))
        open(O_STR,file=outfile)
        call generic_strmol_writer(O_STR,filetype_out,molecule,title=title)
        close(O_STR)
        ! Write g09 input
        outfile=trim(adjustl(basefile))//"_"//trim(adjustl(label))//".com"
        open(O_GAU,file=outfile)
        call write_gcom(O_GAU,molecule,&
                              !Optional args
                              chkname=outfile,&! 
                              calc=calc_type,  &! e.g. SP, Freq, Opt...
                              method=method,   &!
                              basis=basis,     &!
                              title=title      )
        close(O_GAU)
        
        ! Write FCHK
        outfile=trim(adjustl(basefile))//"_"//trim(adjustl(label))//".fchk"
        open(O_FCHK,file=outfile)
        ! Title and job info
        write(O_FCHK,'(A)') "FCHK created with read_scan from "//trim(adjustl(inpfile))
        write(O_FCHK,'(A10,A60,A10)') adjustl(calc_type), adjustl(method), adjustl(basis)
        call write_fchk(O_FCHK,"Number of atoms","I",0,A,(/Nat/),error)
        call write_fchk(O_FCHK,"Charge","I",0,A,(/charge/),error)
        call write_fchk(O_FCHK,"Multiplicity","I",0,A,(/mult/),error)
        N=molecule%natoms
        call write_fchk(O_FCHK,"Atomic numbers",'I',N,A,molecule%atom(1:N)%AtNum,error)
        N=molecule%natoms
        allocate(IA(1:1),A(1:N))
        A(1:N) = float(molecule%atom(1:N)%AtNum)
        call write_fchk(O_FCHK,"Nuclear charges",'R',N,A,IA,error)
        deallocate(A,IA)
        !Coordinates
        N=3*molecule%natoms
        allocate(IA(1:1),A(1:N))
        do i=1,N/3
            j=3*i
            A(j-2) = molecule%atom(i)%x/BOHRtoANGS
            A(j-1) = molecule%atom(i)%y/BOHRtoANGS
            A(j)   = molecule%atom(i)%z/BOHRtoANGS
        enddo
        call write_fchk(O_FCHK,"Current cartesian coordinates",'R',N,A,IA,error)
        deallocate(A,IA)
        !Atomic weights
        call write_fchk(O_FCHK,"Integer atomic weights",'I',3*Nat,A,int(molecule%atom(:)%mass),error)
        call write_fchk(O_FCHK,"Real atomic weights",'R',3*Nat,molecule%atom(:)%mass,IA,error)
        !Energy 
        call write_fchk(O_FCHK,"Total Energy",'R',0,(/E(i_scan)/),IA,error)
        !Gradient
        call write_fchk(O_FCHK,"Cartesian Gradient",'R',3*Nat,Grad,IA,error)
        close(O_FCHK)
    enddo

    write(0,*) 'Number of scan points', i_scan

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
        write(0,'(/,A)') '             R E A D    S C A N '    
        write(0,'(/,A)') '         Read RlxScan steps from FCHK' 
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
       

end program read_scan

