program normal_modes_animation


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
    use ff_build
!   Bond/angle/dihed meassurement
    use atomic_geom
!   Symmetry support
    use symmetry_mod
!   For internal thingies
    use internal_module
    use zmat_manage

    implicit none

    integer,parameter :: NDIM = 800

    !====================== 
    !System variables
    type(str_resmol) :: molecule
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Nred
    character(len=5) :: PG
    !Bonded info
    integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !B and G matrices
    real(8),dimension(1:NDIM,1:NDIM) :: B1, G1
    !Other matrices
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv, L1, Asel1, Asel
    !Auxiliar matrices
    real(8),dimension(1:NDIM,1:NDIM) :: Aux, Aux2
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: Freq, S1, Vec, Smap, mu
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !Shifts
    real(8),dimension(NDIM) :: Delta
    real(8) :: Delta_p
    !Coordinate map
    integer,dimension(NDIM) :: Zmap
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
    real(8) :: Theta, Theta2, dist
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, nn, imin, imax, iii
    integer :: nbonds, nangles, ndihed
    !=============

     !orientation things
    real(8),dimension(3,3) :: ori

    !=================================
    !NM stuff
    !=================================
    real(8) :: Amplitude = 1.d0, qcoord
    !Moving normal modes
    integer,dimension(1:1000) :: nm=0
    real(8) :: Qstep, d, rmsd1, rmsd2
    logical :: call_vmd = .false.
    character(len=10000) :: vmdcall
    integer :: Nsteps, Nsel, istep

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_ZMAT=11, &
               O_STAT=22, &
               S_VMD=30
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile ="input.fchk",  &
                         outfile   
    character(len=1)  :: cnull       
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

    ! 0. GET COMMAND LINE ARGUMENT: input file name
    call getarg(1, inpfile) 

    ! 1. INTERNAL VIBRATIONAL ANALYSIS 
 
    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Read structure
    call generic_strfile_read(I_INP,filetype,molecule)
    !Shortcuts
    Nat = molecule%natoms
    Nvib = 3*Nat-6
    !Read the Hessian: only two possibilities supported
    if (adjustl(filetype) == "log") then
        !Gaussian logfile
        allocate(props)
        call parse_summary(I_INP,molecule,props,"read_hess")
        !Caution: we NEED to read the Freq summary section
        if (adjustl(molecule%job%type) /= "Freq") &
          call alert_msg( "fatal","Section from the logfile is not a Freq calculation")
        ! RECONSTRUCT THE FULL HESSIAN
        k=0
        do i=1,3*Nat
            do j=1,i
                k=k+1
                Hess(i,j) = props%H(k) 
                Hess(j,i) = Hess(i,j)
            enddo
        enddo
        deallocate(props)
        !No job info read for the moment. Use "sensible" defaults
        molecule%job%title = ""
        molecule%job%type= "SP"
    else if (adjustl(filetype) == "fchk") then
        !FCHK file    
        call read_fchk(I_INP,"Cartesian Force Constants",dtype,N,A,IA,error)
        ! RECONSTRUCT THE FULL HESSIAN
        k=0
        do i=1,3*Nat
            do j=1,i
                k=k+1
                Hess(i,j) = A(k) 
                Hess(j,i) = Hess(i,j)
            enddo
        enddo
        deallocate(A)
        !Read job info
        call get_jobtype_fchk(I_INP,molecule%job,error)
        molecule%job%type= "SP"

    else if (adjustl(filetype) == "UnSym") then
        call read_molcas_hess(I_INP,N,Hess,error)
    endif
    close(I_INP)


    !====================================
    !INTERNAL COORDINATES MANAGEMENT
    !====================================
    ! Get connectivity from the residue (needs to be in ANGS, as it is -- default coord. output)
    ! Setting element from atom names is mandatory to use guess_connect
    call guess_connect(molecule)
    !From now on, we'll use atomic units
    molecule%atom(1:Nat)%x = molecule%atom(1:Nat)%x/BOHRtoANGS
    molecule%atom(1:Nat)%y = molecule%atom(1:Nat)%y/BOHRtoANGS
    molecule%atom(1:Nat)%z = molecule%atom(1:Nat)%z/BOHRtoANGS

    !Generate bonded info
    call gen_bonded(molecule)
! print*, "Nbonds", molecule%geom%nbonds
! print*, "Nangls", molecule%geom%nangles
! print*, "Ndihed", molecule%geom%ndihed

    !GENERATE SET FOR INTERNAL COORDINATES FROM Z-MATRIX
    call build_Z(molecule,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
    !Z-mat (always Z-mat is used)
    molecule%geom%bond(1:Nat-1,1:2) = bond_s(2:Nat,1:2)
    molecule%geom%angle(1:Nat-2,1:3) = angle_s(3:Nat,1:3)
    molecule%geom%dihed(1:Nat-3,1:4) = dihed_s(4:Nat,1:4)
    molecule%geom%nbonds  = Nat-1
    molecule%geom%nangles = Nat-2
    molecule%geom%ndihed  = Nat-3


    !INTERNAL COORDINATES ANALYSIS
    Asel1(1,1) = 99.d0 !out-of-range, as Asel is normalized -- this option is not tested
    S_sym(3*Nat) = 0
    call internal_Wilson(molecule,S1,S_sym,ModeDef,B1,G1,Asel1,verbose)


    !SOLVE GF METHOD TO GET NM AND FREQ
    !For redundant coordinates a non-redundant set is formed as a combination of
    !the redundant ones. The coefficients for the combination are stored in Asel
    !as they must be used for state 2 (not rederived!).
    call gf_method(Hess,molecule,S_sym,ModeDef,L1,B1,G1,Freq,Asel1,X1,X1inv,verbose) 


    !Compute new state_file for 2
    ! T2(g09) = mu^1/2 m B^t G2^-1 L2
    ! Compute G1^-1 (it is X1inv * X1inv
    Aux(1:Nvib,1:Nvib) = matmul(X1inv(1:Nvib,1:Nvib),X1inv(1:Nvib,1:Nvib))
    ! Compute B1^t G1^-1
    do i=1,3*Nat
    do j=1,Nvib
        Aux2(i,j) = 0.d0
        do k=1,Nvib
         Aux2(i,j)=Aux2(i,j)+B1(k,i)*Aux(k,j)
        enddo
    enddo
    enddo
    ! Compute [B1^t G2^-1] L2
    Aux(1:3*Nat,1:Nvib) = matmul(Aux2(1:3*Nat,1:Nvib),L1(1:Nvib,1:Nvib))
    ! Compute mu^1/2 m [B^t G2^-1 L2] (masses are in UMA in the fchk)
!     print*, state2%atom(1)%name, state2%atom(1)%mass
!     print*, mu(1), mu(Nvib)
    i=0
    do k=1,Nat
    do kk=1,3
    i=i+1
    do j=1,Nvib
        Aux2(i,j) = 1.d0/&!dsqrt(mu(j)*UMAtoAU)        / &
                    molecule%atom(k)%mass/UMAtoAU * &
                    Aux(i,j)
    enddo
    enddo
    enddo
    !Compute reduced masses
    do j=1,Nvib
    mu(j)=0.d0
    do i=1,3*Nat
        mu(j)=mu(j)+Aux2(i,j)**2
    enddo
    mu(j) = 1.d0/mu(j)
    enddo
    !Normalize with mu
    i=0
    do k=1,Nat
    do kk=1,3
    i=i+1
    do j=1,Nvib
        Aux2(i,j) = Aux2(i,j)*dsqrt(mu(j))
    enddo
    enddo
    enddo
    !Checking normalization
    print*, ""
    print*, "Checking normalization of Tcart (G09)"
    do j=1,Nvib
    theta=0.d0
    do i=1,3*Nat
        theta=theta+Aux2(i,j)**2
    enddo
    print*, theta
    enddo
    print*, ""
    !Print state
    call split_line(inpfile,".",outfile,cnull)
    outfile = "state_"//trim(adjustl(outfile))//"_hess"
    open(O_STAT,file=outfile)
    do i=1,Nat
        write(O_STAT,*) molecule%atom(i)%x*BOHRtoANGS
        write(O_STAT,*) molecule%atom(i)%y*BOHRtoANGS
        write(O_STAT,*) molecule%atom(i)%z*BOHRtoANGS
    enddo
    do i=1,3*Nat
    do j=1,Nvib
        write(O_STAT,*) Aux2(i,j)
    enddo
    enddo
    do j=1,Nvib
        write(O_STAT,'(F12.5)') Freq(j)
    enddo
    close(O_STAT)


    
    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================


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
       

end program normal_modes_animation

