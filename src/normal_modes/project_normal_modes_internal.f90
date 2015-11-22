program projection_normal_modes_int


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to analyse vibrations in term of internal coordinates.
    !
    ! Compilation instructions (for mymake script):
    !
    ! Change log:
    !
    ! TODO:
    ! ------
    !
    ! History
    ! v0: Adapted from nm_internal (v4.0.1.1)
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
    use psi4_manage
!   Structural parameters
    use molecular_structure
    use ff_build
!   Bond/angle/dihed meassurement
    use atomic_geom
!   Symmetry support
    use symmetry_mod
!   For internal thingies
    use internal_module,only:internal_Wilson,gf_method
    use zmat_manage

    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: nosym=.true.   ,&
               zmat=.true.    ,&
               tswitch=.false.,&
               symaddapt=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule, molec_aux
    type(str_job)    :: job
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
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv, L1, Asel
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(1:NDIM) :: Freq, S1, S2, Vec, Smap, Factor
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
    integer :: i,j,k,l, ii,jj,kk, iat, nn, imin, imax, iii, &
               i1, i2, i3, i4
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

    !=================================
    ! Internal scan stuff
    !=================================
    integer :: icoord
    logical :: scan_internal=.false.,&
               showZ=.false.

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_ZMAT=11, & 
               I_STR=12,  &
               I_RED=13,  &
               I_ADD=14,  &
               O_GRO=20,  &
               O_G09=21,  &
               O_Q  =22,  &
               O_NUM=23,  &
               S_VMD=30
    !files
    character(len=10) :: filetype="guess", &
                         filetype_str="guess"
    character(len=200):: inpfile ="input.fchk",  &
                         addfile ="additional.input", &
                         strfile ="struct.pdb",  &
                         zmatfile="NO",          &
                         numfile       
    character(len=100),dimension(1:1000) :: grofile
    character(len=100) :: g09file,qfile
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
    nm(1) = 0
    icoord=-1
    call parse_input(inpfile,strfile,nm,Nsel,Amplitude,filetype,filetype_str,nosym,zmat,verbose,tswitch,symaddapt,zmatfile,&
                     addfile,icoord,showZ,call_vmd)
    if (icoord /= -1) scan_internal=.true.
    if (showZ) scan_internal=.false.

    ! 1. INTERNAL VIBRATIONAL ANALYSIS 
 
    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Read structure
    if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line_back(inpfile,".",null,filetype)
    endif
    call generic_strfile_read(I_INP,filetype,molecule)
    !Shortcuts
    Nat = molecule%natoms
    Nvib = 3*Nat-6
    !Read the Hessian
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
    else if (adjustl(filetype) == "psi4") then
        N=molecule%natoms*3
        call read_psi_hess(I_INP,N,Hess,error)
    else if (adjustl(filetype) == "g96") then
        !The hessian should be given as additional input
        if (adjustl(addfile) == "additional.input") &
         call alert_msg("fatal","With a g96, and additional file should be provided with the Hessian")
        open(I_ADD,file=addfile,status="old")
        call read_gro_hess(I_ADD,N,Hess,error)
        close(I_ADD)
    endif
    close(I_INP)

    !NORMAL MODES SELECTION SWITCH
    if (Nsel == 0) then
        !The select them all
        Nsel = Nvib
        do i=1,Nsel
            nm(i) = i
        enddo
    endif
    if (Nsel > 1000) call alert_msg("fatal", "Too many normal modes. Dying")


    !====================================
    !INTERNAL COORDINATES MANAGEMENT
    !====================================
    ! Get connectivity from the residue (needs to be in ANGS, as it is -- default coord. output)
    ! Setting element from atom names is mandatory to use guess_connect
    call guess_connect(molecule)
    if (nosym) then
        PG="C1"
    else
        call symm_atoms(molecule,isym)
        PG=molecule%PG
    endif
    !From now on, we'll use atomic units
    molecule%atom(1:Nat)%x = molecule%atom(1:Nat)%x/BOHRtoANGS
    molecule%atom(1:Nat)%y = molecule%atom(1:Nat)%y/BOHRtoANGS
    molecule%atom(1:Nat)%z = molecule%atom(1:Nat)%z/BOHRtoANGS

    !Generate bonded info
    call gen_bonded(molecule)

    !GEN BONDED SET FOR INTERNAL COORD
    if (.not.zmat) then
!         print*, "Custom internal coordianates"
!         open(I_RED,file="modred.dat") 
!         call modredundant(I_RED,molecule)
!         close(I_RED)
!         zmat=.false.
        print*, "Using all internals"
! print*, "Review"
!         do i=1,molecule%geom%ndihed
!             i1=molecule%geom%dihed(i,1)
!             i2=molecule%geom%dihed(i,2)
!             i3=molecule%geom%dihed(i,3)
!             i4=molecule%geom%dihed(i,4)
!             print*, "Dihedral", i1, i2, i3, i4
!             print*, calc_improper(molecule%atom(i1),molecule%atom(i2),molecule%atom(i3),molecule%atom(i4))*180.d0/PI
!         enddo
    elseif (zmat) then
        if (adjustl(zmatfile) == "NO") then
            call build_Z(molecule,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
        else
            open(I_ZMAT,file=zmatfile,status="old")
            print*, "Z-matrix read from "//trim(adjustl(zmatfile))
            call read_Z(molecule,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym,I_ZMAT)
            close(I_ZMAT)
            !Deactivate symaddapt (for the moment)
            PG = "C1"
        endif
        !Z-mat
        molecule%geom%bond(1:Nat-1,1:2) = bond_s(2:Nat,1:2)
        molecule%geom%angle(1:Nat-2,1:3) = angle_s(3:Nat,1:3)
        molecule%geom%dihed(1:Nat-3,1:4) = dihed_s(4:Nat,1:4)
        molecule%geom%nbonds  = Nat-1
        molecule%geom%nangles = Nat-2
        molecule%geom%ndihed  = Nat-3
    endif !otherwise all parameters are used
    !Set number of redundant
    Nred = molecule%geom%nbonds  + &
           molecule%geom%nangles + &
           molecule%geom%ndihed

    !Set symmetry of internal (only if symmetry is detected)
!     if (adjustl(PG) == "C1") then
!         S_sym(3*Nat) = 1
!     else
!         do i=1,Nat-1
!             S_sym(i) = bond_sym(i+1)-1
!         enddo
!         do i=1,Nat-2
!             S_sym(i+Nat-1) = angle_sym(i+2)+Nat-3
!         enddo
!         do i=1,Nat-3
!             S_sym(i+2*Nat-3) = dihed_sym(i+3)+2*Nat-6
!         enddo
!     endif

!     !We send the option -sa within S_sym (confflict with redundant coord!!)
!     if (symaddapt) then
!         S_sym(3*Nat) = 1
!     else
!         S_sym(3*Nat) = 0
!     endif

!     print*, "Internal symm"
!     do i=1,Nvib
!         print*, i,S_sym(i)
!     enddo


    !INTERNAL COORDINATES ANALYSIS
    Asel(1,1) = 99.d0 !out-of-range, as Asel is normalized -- this option is not tested
    call internal_Wilson(molecule,S1,S_sym,ModeDef,B1,G1,Asel,verbose)
    if (showZ) stop

    !SOLVE GF METHOD TO GET NM AND FREQ
        !For redundant coordinates a non-redundant set is formed as a combination of
        !the redundant ones. The coefficients for the combination are stored in Asel
        !as they must be used for state 2 (not rederived!).
        call gf_method(Hess,molecule,S_sym,ModeDef,L1,B1,G1,Freq,Asel,X1,X1inv,verbose) 

        !Compute the inverse (Linv), stored in X1inv
        X1inv(1:Nvib,1:Nvib)=L1(1:Nvib,1:Nvib)
        call dgetrf(Nvib,Nvib, X1inv, NDIM, ipiv, info)
        call dgetri(Nvib, X1inv, NDIM, ipiv, work, NDIM, info)

        X1(1:Nvib,1:Nvib)=matmul(X1inv(1:Nvib,1:Nvib),L1(1:Nvib,1:Nvib))
!         print*, "L1 L1^-1="
!         do i=1,Nvib
!             print'(100(F8.3,2X))', X1(i,1:Nvib)
!         enddo
!         print*, ""

!         ii=0
!         do i=1,Nred
!             if (Freq(i) > 1.d-2) then
!                 ii=ii+1
!                 Freq(ii) = Freq(i)
!                 L1(1:Nred,ii) = L1(1:Nred,i)
!              endif
!         enddo

       !Use freqs. to make displacements equivalent in dimensionless units
!         Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib)))/5.d3
        !Factor to convert to adim units (from SI)
        Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib))*1.d2*clight*2.d0*PI/plankbar)

!     endif

    if (verbose) then
        print*, "L1="
        do i=1,Nred
            print'(100(F8.3,2X))', L1(i,1:Nred)
        enddo
        print*, "Selected Freq="
        do i=1,Nvib
            print'(100(F8.3,2X))', Freq(i)
        enddo
    endif


    !==========================================================0
    !  Projection on Normal modes
    !==========================================================0
    ! 1. Read structure
    open(I_INP,file=strfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(strfile)) )
    if (adjustl(filetype_str) == "guess") then
        ! Guess file type
        call split_line_back(strfile,".",null,filetype_str)
    endif
    call generic_strfile_read(I_INP,filetype_str,molec_aux)
    !From now on, we'll use atomic units
    molec_aux%atom(1:Nat)%x = molec_aux%atom(1:Nat)%x/BOHRtoANGS
    molec_aux%atom(1:Nat)%y = molec_aux%atom(1:Nat)%y/BOHRtoANGS
    molec_aux%atom(1:Nat)%z = molec_aux%atom(1:Nat)%z/BOHRtoANGS
    ! 2. Get internal coordinates
    molec_aux%geom = molecule%geom
    call internal_Wilson(molec_aux,S2,S_sym,ModeDef,B1,G1,Asel,verbose)
    ! 3. Compute DeltaS
    print*, ""
    print*, "=========================="
    print*, " SHIFTS (internal coord)"
    print*, "=========================="
    print*, "Bonds"
    do i=1,molecule%geom%nbonds
        Delta(i) = S2(i)-S1(i)
        print'(I5,3(F8.2,2X))', i, S2(i),S1(i), Delta(i)
    enddo
    print*, "Angles"
    do j=i,i+molecule%geom%nangles-1
        Delta(j) = S2(j)-S1(j)
        print'(I5,3(F8.2,2X))', j, S2(j)*180.d0/PI,S1(j)*180.d0/PI,Delta(j)*180.d0/PI
    enddo
    print*, "Dihedrals"
    do k=j,j+molecule%geom%ndihed-1
          if (S1(k) > S2(k)) then
              Delta_p = S2(k)+2.d0*PI - S1(k) 
          else
              Delta_p = S2(k) - S1(k)
          endif

          if (Delta_p < -PI) then
              Delta(k) = Delta_p + 2*PI
          elseif (Delta_p > PI) then
              Delta(k) = Delta_p - 2*PI
          else
              Delta(k) = Delta_p
          endif
!         Delta(k) = S2(k)-S1(k)
!         Delta_p = S2(k)-S1(k)+2.d0*PI
!         if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
!         Delta_p = S2(k)-S1(k)-2.d0*PI
!         if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
        print'(I5,3(F8.2,2X))', k, S2(k)*180.d0/PI,S1(k)*180.d0/PI,Delta(k)*180.d0/PI
    enddo
    ! 4. Compute DeltaQ
    do j=1,Nvib
        qcoord = 0.d0
        do i=1,Nvib
            qcoord = qcoord + X1inv(j,i) * Delta(i)
        enddo
        !fort.98 in atomic units
        write(98,*) j, qcoord
        qcoord = qcoord*BOHRtoM*dsqrt(AUtoKG) ! SI units
        !Adim
        qcoord = qcoord*Factor(j)
        !fort.99 in dimesionless disp
        write(99,*) j, qcoord
    enddo

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,strfile,nm,Nsel,Amplitude,filetype,filetype_str,nosym,zmat,verbose,tswitch,symaddapt,zmatfile,&
                           addfile,icoord,showZ,call_vmd)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,strfile,filetype,filetype_str,zmatfile,addfile
        logical,intent(inout) :: nosym, verbose, zmat, tswitch, symaddapt,showZ,call_vmd
        integer,dimension(:),intent(inout) :: nm
        integer,intent(inout) :: icoord
        integer,intent(out) :: Nsel
        real(8),intent(out) :: Amplitude
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg

        !Prelimirary defaults
        Nsel = 0

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

                case ("-add") 
                    call getarg(i+1, addfile)
                    argument_retrieved=.true.

                case ("-nm") 
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    call string2vector_int(arg,nm,Nsel)

                case ("-fs") 
                    call getarg(i+1, strfile)
                    argument_retrieved=.true.
                case ("-fts") 
                    call getarg(i+1, filetype_str)
                    argument_retrieved=.true.

                case ("-maxd") 
                    call getarg(i+1, arg)
                    read(arg,*) Amplitude
                    argument_retrieved=.true.

                case ("-vmd")
                    call_vmd=.true.

                case ("-nosym")
                    nosym=.true.
                case ("-sym")
                    nosym=.false.

                case ("-sa")
                    symaddapt=.true.
                case ("-nosa")
                    symaddapt=.false.

                case ("-readz") 
                    call getarg(i+1, zmatfile)
                    argument_retrieved=.true.

                case("-showz")
                    showZ=.true.

                case ("-zmat")
                    zmat=.true.
                case ("-nozmat")
                    zmat=.false.

                case ("-tswitch")
                    tswitch=.true.

                case ("-int")
                    call getarg(i+1, arg)
                    read(arg,*) icoord
                    argument_retrieved=.true.

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
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '          INTERNAL MODES ANIMATION '    
        write(0,'(/,A)') '      Perform vibrational analysis based on  '
        write(0,'(/,A)') '            internal coordinates (D-V7)'        
        write(0,'(/,A)') '         Revision: nm_internal-140320-1'
       write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(filetype))
        write(0,*) '-add            ', trim(adjustl(addfile))
!         write(0,*) '-nm            ', nm(1),"-",nm(Nsel)
        write(0,*) '-fs             ', trim(adjustl(strfile))
        write(0,*) '-fts            ', trim(adjustl(filetype_str))
!         write(0,*) '-vmd           ',  call_vmd
!         write(0,*) '-maxd          ',  Amplitude
        if (nosym) dummy_char="NO "
        if (.not.nosym) dummy_char="YES"
        write(0,*) '-[no]sym        ', dummy_char
        if (zmat) dummy_char="YES"
        if (.not.zmat) dummy_char="NO "
        write(0,*) '-[no]zmat       ', dummy_char
        write(0,*) '-readz          ', trim(adjustl(zmatfile))
        write(0,*) '-showz         ', showZ
!         if (tswitch) dummy_char="YES"
!         if (.not.tswitch) dummy_char="NO "
!         write(0,*) '-tswitch        ', dummy_char
!         if (symaddapt) dummy_char="YES"
!         if (.not.symaddapt) dummy_char="NO "
!         write(0,*) '-sa             ', dummy_char
!         write(0,*) '-int           ', icoord
        write(0,*) '-v             ', verbose
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

        select case (adjustl(filetype))
            case("xyz")
             call read_xyz(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("gro")
             call read_gro(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("g96")
             call read_g96(I_INP,molec)
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
                        //". Try forcing the filetype with -ft flag")
        end select



        return


    end subroutine generic_strfile_read
       

end program projection_normal_modes_int

