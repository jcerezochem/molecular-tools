program reorder_fchk


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
    real(8),dimension(:),allocatable :: Hlt
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, Hess_aux
    !VECTORS
    real(8),dimension(NDIM) :: Grad, Grad_aux
    ! Rotation matrix 
    real(8),dimension(3,3) :: R
    !====================== 

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
    integer :: error, N, lenght
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
    character(len=240):: line
    integer :: iat_new, iat_orig, j_orig, j_new, nswap
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, nn, imin, imax, iii
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_ORD=11,  &
               I_ROT=12,  &
               O_FCHK=20
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile ="input.fchk",  &
                         orderfile = "none", &
                         outfile="output.fchk", &
                         rotfile="none"
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
    call parse_input(inpfile,filetype,orderfile,rotfile,outfile)
 
    !================
    ! READ DATA
    !================
    if (adjustl(filetype) == "guess") &
    call split_line_back(inpfile,".",null,filetype)

    if (adjustl(filetype) /= "log" .and. adjustl(filetype) /= "fchk") &
        call alert_msg("fatal","Only Gaussian log and fchk files supported")

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

    ! HESSIAN FILE
    have_hessian=.false.
    if (adjustl(calc_type) == "Freq") then
        rewind(I_INP)
        print'(X,A)', "READING HESSIAN..."
        allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
        call generic_Hessian_reader(I_INP,filetype,Nat,Hlt,error) 
        Hess(1:3*Nat,1:3*Nat) = Hlt_to_Hess(3*Nat,Hlt)
        print'(X,A,/)', "Done"
        have_hessian=.true.
    endif

    ! GRADIENT FILE
    have_gradient=.false.
    if (adjustl(calc_type) == "Freq") then
        ! Freq calcs have gradient in log and fchk
        have_gradient=.true.
    elseif (adjustl(calc_type) == "Opt"   .or. &
            adjustl(calc_type) == "FOpt"  .or. &
            adjustl(calc_type) == "Force") then
        ! Other calcs only have gradient on fchk
        if (adjustl(filetype) == "fchk") then
            have_gradient=.true.
        else
            have_gradient=.false.
        endif
    endif
    if (have_gradient) then
        rewind(I_INP)
        print'(X,A)', "READING GRADIENT..."
        call generic_gradient_reader(I_INP,filetype,Nat,Grad,error) 
        print'(X,A,/)', "Done"
    endif

    ! ENERGIES
    print'(X,A)', "READING ENERGIES..."
    have_SCF=.false.
    have_TD =.false.
    have_TOT=.false.
    if (adjustl(filetype) == "fchk") then
        call read_fchk(I_INP,"SCF Energy",dtype,N,A,IA,error)
        if (error == 0) then
            E_scf = A(1)
            have_SCF = .true.
            deallocate(A)
        endif
        call read_fchk(I_INP,"CIS Energy",dtype,N,A,IA,error)
        if (error == 0) then
            E_td = A(1)
            have_TD = .true.
            deallocate(A)
        endif
        call read_fchk(I_INP,"Total Energy",dtype,N,A,IA,error)
        if (error == 0) then
            E_tot = A(1)
            have_TOT = .true.
            deallocate(A)
        endif
    elseif (adjustl(filetype) == "log") then
        rewind(I_INP)
        ! we need to first compute the section lenght where properties are 
        call estimate_section_length(I_INP,5,lenght)
        rewind(I_INP)
        error=lenght ! we give the section lenght in the input error_flag
        call read_gausslog_property(I_INP,"HF",line,error)
        if (error == 0) then
            read(line,*) E_scf
            have_SCF = .true.
            !Update Total Energy with this value
            E_tot = E_scf
            have_TOT = .true.
        endif
        rewind(I_INP)
        call read_gauslog_tdenergy(I_INP,E_td,error)
        if (error == 0) then
            have_TD  = .true.
            !Update Total Energy with this value
            E_tot = E_td
            have_TOT = .true.
        endif
    endif
    if (have_SCF) print*, "  SCF Energy"
    if (have_TD ) print*, "  CIS Energy"
    if (have_TOT) print*, "  Total Energy"
    print'(X,A,/)', "Done"
    close(I_INP)


    !================
    !CHANGE ORDER (if needed)
    !================
    if (adjustl(orderfile)/="none") then
        print'(X,A)', "REORDERING:"
        print'(X,A)', "  STRUCTURE..."
        open(I_ORD,file=orderfile,status="old")
        ! Only change the selected ones
        read(I_ORD,*) nswap
        ! initialize iord
        iord(1:3*molecule%natoms) = [ (i, i=1,3*molecule%natoms) ]
        molec_aux = molecule
        do i=1,nswap !molecule%natoms
            read(I_ORD,*) iat_orig, iat_new
            molec_aux%atom(iat_new)=molecule%atom(iat_orig)
            !Build an ordering array with the 3N elements: (x1,y1,z1,x2,y2...) -> (x1',y1',z1',x2',y2'...)
            j_orig = 3*iat_orig
            j_new  = 3*iat_new
            iord(j_new-2) = j_orig-2
            iord(j_new-1) = j_orig-1
            iord(j_new)   = j_orig
        enddo
        close(I_ORD)
        molecule = molec_aux

        ! Reorder Gradient and Hessian if present
        if (have_gradient) then
            print'(X,A)', "  GRADIENT..."
            Grad_aux(1:3*Nat) = Grad(1:3*Nat)
            do i=1,3*molecule%natoms
                Grad(i) = Grad_aux(iord(i))
            enddo
        endif
        if (have_hessian) then
            print'(X,A)', "  HESSIAN..."
            Hess_aux(1:3*Nat,1:3*Nat) = Hess(1:3*Nat,1:3*Nat)
            do i=1,3*molecule%natoms
                do j=1,i
                    Hess(i,j) = Hess_aux(iord(i),iord(j))
                    Hess(j,i) = Hess(i,j)
                enddo
            enddo
        endif
        print'(X,A,/)', "Done"
    endif

    !================
    !ROTATE (if needed)
    !================
    if (adjustl(rotfile)/="none") then
        print'(X,A)', "ROTATING:"
        open(I_ROT,file=rotfile,status="old")
        ! Read rotation matrix.
        do i=1,3
            read(I_ROT,*) R(i,1:3)
        enddo
        close(I_ROT)
        print'(X,A)', "  STRUCTURE..."
        call rotate_molec(molecule,R)

        ! Reorder Gradient and Hessian if present
        if (have_gradient) then
            print'(X,A)', "  GRADIENT..."
            ! Grad' = R Grad
            Grad(1:3*Nat) = rotate3D_vector(3*Nat,Grad,R)
        endif
        if (have_hessian) then
            print'(X,A)', "  HESSIAN..."
            ! Rotate useing rotate3D_matrix
            ! Hess' = R Hess R^t
            !  [R Hess]
            Hess(1:3*Nat,1:3*Nat) = rotate3D_matrix(3*Nat,3*Nat,Hess,R)
            ! and in the second step
            ! (Hess')^t = R [R Hess]^t = Hess' (symmetric Hessian)
            Hess(1:3*Nat,1:3*Nat) = rotate3D_matrix(3*Nat,3*Nat,Hess,R,tA=.true.)
        endif
        print'(X,A,/)', "Done"
    endif


    !================
    !REWRITE FCHK
    !================
    print'(X,A)', "WRITTING FCHK..."
    print'(X,A)', "  Output file: "//trim(adjustl(outfile))
    open(O_FCHK,file=outfile)
    ! Title and job info
    write(O_FCHK,'(A)') "FCHK created with reform_fchk from "//trim(adjustl(inpfile))
    write(O_FCHK,'(A10,A60,A10)') adjustl(calc_type), adjustl(method), adjustl(basis)
    call write_fchk(O_FCHK,"Number of atoms","I",0,A,(/Nat/),error)
    call write_fchk(O_FCHK,"Charge","I",0,A,(/charge/),error)
    call write_fchk(O_FCHK,"Multiplicity","I",0,A,(/mult/),error)

!   THIS IS NOT SAFE. IF THE STARTING FCHK IS NOT STANDARD, IT MAY NEVER REACH "Atomic numbers" SECTION
!     !If fchk, copy first lines
!     if (adjustl(filetype) == "fchk") then
!         open(I_INP,file=inpfile,status='old',iostat=IOstatus)
!         !Copy lines till "Atomic Numbers" if this is a fchk
!         ! skipping the following lines:
!         read (I_INP,'(A)') line ! Title is changed
!         read (I_INP,'(A)') line ! Job info
!         read (I_INP,'(A)') line ! Number of atoms
!         !--
!         ! Info1-1 will be in a different place wrt original
!         read (I_INP,'(A)') line ! Info1-9 (1/2)
!         write(O_FCHK,'(A)') trim(adjustl(line))
!         read (I_INP,'(A)') line ! Info1-9 (2/2)
!         write(O_FCHK,'(A)') trim(adjustl(line))
!         ! skipping the following lines:
!         read (I_INP,'(A)') line ! Charge
!         read (I_INP,'(A)') line ! Multiplicity
!         !--
!         read(I_INP,'(A)',iostat=IOstatus) line
!         do while (index(line,"Atomic numbers")==0)
!             write(O_FCHK,'(A)') trim(adjustl(line))
!             read(I_INP,'(A)',iostat=IOstatus) line
!         enddo
!         close(I_INP)
!     endif

    !Atomic Numbers and Nuclear charges
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
    if (have_SCF) &
        call write_fchk(O_FCHK,"SCF Energy",'R',0,(/E_scf/),IA,error)
    if (have_TD) &
        call write_fchk(O_FCHK,"CIS Energy",'R',0,(/E_td/),IA,error)
    if (have_TOT) &
        call write_fchk(O_FCHK,"Total Energy",'R',0,(/E_tot/),IA,error)
    !Gradient
    if (have_gradient) then
        call write_fchk(O_FCHK,"Cartesian Gradient",'R',3*Nat,Grad,IA,error)
    endif
    !Hessian
    if (have_hessian) then
        ! Transform the Hess into Hlt 
        N=3*Nat*(3*Nat+1)/2
        Hlt(1:N) = Hess_to_Hlt(3*Nat,Hess)
        call write_fchk(O_FCHK,"Cartesian Force Constants",'R',N,Hlt,IA,error)
    endif

    print'(X,A,/)', "Done"

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    close(O_FCHK)

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,orderfile,rotfile,outfile)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype,orderfile,outfile,rotfile
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
                case ("-reor") 
                    call getarg(i+1, orderfile)
                    argument_retrieved=.true.      
                case ("-rot") 
                    call getarg(i+1, rotfile)
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

       !Print options (to stdx)
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '             R E F O R M   F C H K '    
        write(0,'(/,A)') '      Generate fchk file from log or fchk'        
        write(0,'(A)')   '        reordering or rotating the data' 
        write(0,'(A)')   '        (geometry, gradient and Hessian)' 
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input file                       ', trim(adjustl(inpfile))
        write(0,*)       '-ft          \_ FileTyep                      ', trim(adjustl(filetype))
        write(0,*)       '-reor        File with reordering instruction ', trim(adjustl(orderfile))
        write(0,*)       '-rot         File with 3x3 rotation matrix    ', trim(adjustl(rotfile))
        write(0,*)       '-f           Output file (fchk)               ', trim(adjustl(outfile))
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program reorder_fchk

