program modify_fchk


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
    use molcas_unsym_manage
!   Structural parameters
    use molecular_structure

    implicit none

    integer,parameter :: NDIM = 800

    !====================== 
    !Options 
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule, molec_aux
    type(str_job)    :: job
    integer,dimension(1:NDIM) :: iord
    integer :: Nat, Nvib
    character(len=5) :: PG
    logical :: hess_present=.true.
    !====================== 

    !====================== 
    !MATRICES
    !Other matrices
    real(8),dimension(1:NDIM,1:NDIM) :: Hess1, Hess2
    !VECTORS
!     real(8),dimension(NDIM) :: Grad1, Grad2, Grad_aux
    !ENERGY
    real(8) :: E_mod
    character(len=10) :: energy_units
    !====================== 

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
    integer :: error, N
    !Read gaussian log auxiliars
    type(str_molprops),allocatable :: props
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
    character(len=50) :: dummy_char
    integer :: iat_new, iat_orig, j_orig, j_new
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, nn, imin, imax, iii
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_MOD=12,  &
               I_ENE=13,  &
               O_FCHK=20
    !files
    character(len=10) :: filetype= "guess", &
                         filetype_inp
    character(len=200):: inpfile = "input.fchk",  &
                         modhess = "hess.gmx",  &
                         enefile = "energy.dat", &
                         outfile = "output.fchk"
    !Control of stdout
    logical :: verbose=.false.
    !status
    integer :: IOstatus
    !===================

    !===========================
    ! Advanced file reading
    character(len=240) :: line
    character(len=50),dimension(5)  :: args
    integer :: nargs
    !============================

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

! (End of variables declaration) 
!==================================================================================
    call cpu_time(ti)

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,modhess,filetype,enefile,outfile)
 
    ! 1. READ DATA
    ! ---------------------------------
    ! 1a. REFERENCE HESSIAN AND STRUCTURE (fchk only)
    !Read structure
    call split_line_back(inpfile,".",null,filetype_inp)
    if (adjustl(filetype_inp) /= "fchk") &
        call alert_msg( "fatal","Expected input format is fchk but got "//trim(adjustl(filetype_inp)) )
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    call generic_strfile_read(I_INP,filetype_inp,molecule)
    !Shortcuts
    Nat = molecule%natoms
    Nvib = 3*Nat-6

    !FCHK file    
    call read_fchk(I_INP,"Cartesian Force Constants",dtype,N,A,IA,error)
    if (error == 0) then
        ! RECONSTRUCT THE FULL HESSIAN
        k=0
        do i=1,3*Nat
            do j=1,i
                k=k+1
                Hess1(i,j) = A(k) 
                Hess1(j,i) = Hess1(i,j)
            enddo
        enddo
        deallocate(A)
    else
        call alert_msg( "note","Hessian not found in the fchk")
        hess_present=.false.
    endif


    if (hess_present) then
        ! 1b. ADDITIONAL HESSIAN
        open(I_MOD,file=modhess,status='old',iostat=IOstatus)
        if (adjustl(filetype) == "guess") call split_line_back(modhess,".",null,filetype)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(modhess)) )
    
        !Read the and Gradient Hessian: only on supported files
        if (adjustl(filetype) == "fchk") then
            !FCHK file    
            call read_fchk(I_MOD,"Cartesian Force Constants",dtype,N,A,IA,error)
            if (error == 0) then
                ! RECONSTRUCT THE FULL HESSIAN
                k=0
                do i=1,3*Nat
                    do j=1,i
                        k=k+1
                        Hess2(i,j) = A(k) 
                        Hess2(j,i) = Hess2(i,j)
                    enddo
                enddo
                deallocate(A)
            endif
        else if (adjustl(filetype) == "gmx") then
            open(I_MOD,file=modhess,status="old")
            call read_gro_hess(I_MOD,N,Hess2,error)
            close(I_MOD)
        else
            call alert_msg("fatal","Filet type "//filetype//" does not support Gradient/Hessian r/w")
        endif
    
    
        ! MODIFY HESSIAN
print*, Hess1(1,1), Hess2(1,1)
print*, Hess1(3*Nat,3*Nat), Hess2(3*Nat,3*Nat)
        Hess1(1:3*Nat,1:3*Nat) = Hess1(1:3*Nat,1:3*Nat) - Hess2(1:3*Nat,1:3*Nat)
print*, Hess1(1,1)
print*, Hess1(3*Nat,3*Nat)
    endif

    !REWRITE FCHK
    open(O_FCHK,file=outfile)
    !Copy lines till "Atomic Numbers"
    read (I_INP,'(A)') line
    line="Combined: "//trim(adjustl(inpfile))//" + "//trim(adjustl(modhess))
    do while (index(line,"Atomic numbers")==0)
        write(O_FCHK,'(A)') trim(adjustl(line))
        read(I_INP,'(A)',iostat=IOstatus) line
    enddo
    rewind(I_INP)
    !Atomic Numbers and Nuclear charges
    dtype="I"
    N=molecule%natoms
    allocate(IA(1:N),A(1:1))
    IA(1:N) = molecule%atom(1:N)%AtNum
    call write_fchk(O_FCHK,"Atomic numbers",dtype,N,A,IA,error)
    deallocate(A,IA)
    dtype="R"
    N=molecule%natoms
    allocate(IA(1:1),A(1:N))
    A(1:N) = float(molecule%atom(1:N)%AtNum)
    call write_fchk(O_FCHK,"Nuclear charges",dtype,N,A,IA,error)
    deallocate(A,IA)
    !Coordinates
    dtype="R"
    N=3*molecule%natoms
    allocate(IA(1:1),A(1:N))
    do i=1,N/3
        j=3*i
        A(j-2) = molecule%atom(i)%x/BOHRtoANGS
        A(j-1) = molecule%atom(i)%y/BOHRtoANGS
        A(j)   = molecule%atom(i)%z/BOHRtoANGS
    enddo
    call write_fchk(O_FCHK,"Current cartesian coordinates",dtype,N,A,IA,error)
    deallocate(A,IA)
    !Atomic weights
    call read_fchk(I_INP,"Integer atomic weights",dtype,N,A,IA,error)
    if (error == 0) then
        deallocate(IA)
        dtype="I"
        N=molecule%natoms
        allocate(IA(1:N),A(1:1))
        IA(1:N) = int(molecule%atom(1:N)%mass)
        call write_fchk(O_FCHK,"Integer atomic weights",dtype,N,A,IA,error)
        deallocate(A,IA)
        dtype="R"
        N=molecule%natoms
        allocate(IA(1:1),A(1:N))
        A(1:N) = molecule%atom(1:N)%mass
        call write_fchk(O_FCHK,"Real atomic weights",dtype,N,A,IA,error)
        deallocate(A,IA)
    endif
    !Energy 
    call read_fchk(I_INP,"SCF Energy",dtype,N,A,IA,error)
    if (error == 0) then
        N=0
        call write_fchk(O_FCHK,"SCF Energy",dtype,N,A,IA,error)
        deallocate(A)
    endif
    ! Read Energy from energy file
    open(I_ENE,file=enefile,status='old',iostat=IOstatus)
    read(I_ENE,'(A)') line
    close(I_ENE)

    call string2vector_char_new(line,args,nargs," ")
    read(args(1),*) E_mod

    ! Convert energy into AU
    if (nargs > 1 ) then
        read(args(1),*) energy_units
        call set_word_upper_case(energy_units)
        select case (adjustl(energy_units))
            case ("KJ/MOL") 
                E_mod = E_mod / CALtoJ / HtoKCALM
            case default
                call alert_msg("fatal","Unkown units of energy: "//&
                                trim(adjustl(energy_units)))
        end select
    else
        !By default, we expect kJ/mol   
        E_mod = E_mod / CALtoJ / HtoKCALM
    endif   


    allocate(A(1))
    A(1) = E_mod
    call write_fchk(O_FCHK,"Additional Energy",'R',0,A,IA,error)
    deallocate(A)
    call read_fchk(I_INP,"Total Energy",dtype,N,A,IA,error)
    if (error == 0) then
        A(1) = A(1) - E_mod
        N=0
        call write_fchk(O_FCHK,"Total Energy",dtype,N,A,IA,error)
        deallocate(A)
    endif
!     !Gradient
!     call read_fchk(I_INP,"Cartesian Gradient",dtype,N,A,IA,error)
!     if (error == 0) then
!         deallocate(A)
!         dtype="R"
!         N=3*molecule%natoms
!         allocate(IA(1:1),A(1:N))
!         A(1:N) = Grad(1:N)
!         call write_fchk(O_FCHK,"Cartesian Gradient",dtype,N,A,IA,error)
!         deallocate(A,IA)
!     endif
    if (hess_present) then
        !Hessian
        call read_fchk(I_INP,"Cartesian Force Constants",dtype,N,A,IA,error)
        if (error == 0) then
            deallocate(A)
            dtype="R"
            N=3*molecule%natoms*(3*molecule%natoms+1)/2
            allocate(IA(1:1),A(1:N))
            k=0
            do i=1,3*Nat
                do j=1,i
                    k=k+1
                    A(k) = Hess1(i,j)
                enddo
            enddo
            call write_fchk(O_FCHK,"Cartesian Force Constants",dtype,N,A,IA,error)
            deallocate(A,IA)
        endif
    endif

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    close(I_INP)
    close(O_FCHK)

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,modhess,filetype,enefile,outfile)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype,modhess,enefile,outfile
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
                case ("-f1") 
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.
                case ("-ft2") 
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.  
                case ("-f2") 
                    call getarg(i+1, modhess)
                    argument_retrieved=.true.      
                case ("-e") 
                    call getarg(i+1, enefile)
                    argument_retrieved=.true.      
                case("-o")
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '                MODIFY FCHK                       '    
        write(0,'(/,A)') '    Modify energy and Hessian in a FCHK file      '        
        write(0,'(/,A)') '        Revision: modify_fchk-v140103-1           '
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f1             ', trim(adjustl(inpfile))
        write(0,*) '-f2             ', trim(adjustl(modhess))
        write(0,*) '-ft2            ', trim(adjustl(filetype))
        write(0,*) '-e              ', trim(adjustl(enefile))
        write(0,*) '-o              ', trim(adjustl(outfile))
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
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
        call split_line_back(inpfile,".",null,filetype)
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
            case("UnSym")
             call read_molcas_geom(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
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
            case("UnSym")
             call read_molcas_geom(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select
        endif


        return


    end subroutine generic_strfile_read
       

end program modify_fchk

