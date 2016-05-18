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
    use psi4_manage
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

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: nosym=.true.   ,&
               zmat=.true.    ,&
               tswitch=.false.,&
               symaddapt=.false., &
               include_hbonds=.false.
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
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv, L1, Asel1, Asel
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: Freq, S1, Vec, Smap, Factor
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
    logical :: call_vmd = .false., &
               movie_vmd = .false.
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
               I_ADD=12,  &
               O_GRO=20,  &
               O_G09=21,  &
               O_Q  =22,  &
               O_NUM=23,  &
               O_LIS=24,  &
               S_VMD=30
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile ="input.fchk",  &
                         addfile = "additional.input", &
                         zmatfile="NO",          &
                         numfile       
    character(len=100),dimension(1:1000) :: grofile
    character(len=100) :: nmfile='none', g09file,qfile, tmpfile
    !Control of stdout
    logical :: verbose=.false.
    !status
    integer :: IOstatus
    !===================

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

    !===================
    !MOVIE things
    integer :: movie_cycles=0,& !this means no movie
               movie_steps
    !===================

! (End of variables declaration) 
!==================================================================================
    call cpu_time(ti)

    ! 0. GET COMMAND LINE ARGUMENTS
    nm(1) = 0
    icoord=-1
    call parse_input(inpfile,addfile,nmfile,nm,Nsel,Amplitude,filetype,nosym,zmat,verbose,tswitch,symaddapt,&
                     zmatfile,icoord,showZ,call_vmd,movie_cycles,include_hbonds)
    if (icoord /= -1) scan_internal=.true.
! Por qué estaba este switch????
!    if (showZ) scan_internal=.false.
! Más bien debería ser al contrario, para que le deje hacer el internal analisys aunque no haya Hessiana
   if (showZ) scan_internal=.true.

    ! 1. INTERNAL VIBRATIONAL ANALYSIS 
 
    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Read structure
    if (adjustl(filetype) == "guess") call split_line_back(inpfile,".",null,filetype)
    call generic_strfile_read(I_INP,filetype,molecule)
    !Shortcuts
    Nat = molecule%natoms
    Nvib = 3*Nat-6

    !ONLY READ HESSIAN IF NM ANALYSIS IS REQUIRED
    if (.not.scan_internal) then
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
    call guess_connect(molecule,include_hbonds)
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

    !GENERATE SET FOR INTERNAL COORDINATES FROM Z-MATRIX
    ! Store this info in molec%geom module
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
    if (zmat) then
        !Z-mat
        molecule%geom%bond(1:Nat-1,1:2) = bond_s(2:Nat,1:2)
        molecule%geom%angle(1:Nat-2,1:3) = angle_s(3:Nat,1:3)
        molecule%geom%dihed(1:Nat-3,1:4) = dihed_s(4:Nat,1:4)
        molecule%geom%nbonds  = Nat-1
        molecule%geom%nangles = Nat-2
        molecule%geom%ndihed  = Nat-3
    else
        !otherwise all parameters from molec%geom are used
        nbonds = molecule%geom%nbonds
        nangles= molecule%geom%nangles
        ndihed = molecule%geom%ndihed
        print*, "Bonds Map"
        do j=2,Nat
        do i=1,nbonds
            if (bond_s(j,1)==molecule%geom%bond(i,1).and.&
                bond_s(j,2)==molecule%geom%bond(i,2)) then
                print*, "Zmat<->Map Redundant:", j-1,i
                Zmap(j-1) = i
            endif
            if (bond_s(j,2)==molecule%geom%bond(i,1).and.&
                bond_s(j,1)==molecule%geom%bond(i,2)) then
                print*, "Zmat<->Map Redundant:", j-1,i
                Zmap(j-1) = i
            endif
        enddo
        enddo
        print*, "Angles Map"
        do j=3,Nat
        do i=1,nangles
            if (angle_s(j,1)==molecule%geom%angle(i,1).and.&
                angle_s(j,2)==molecule%geom%angle(i,2).and.&
                angle_s(j,3)==molecule%geom%angle(i,3)) then
                print*, "Zmat<->Map Redundant:", j-3+Nat,i+nbonds
                Zmap(j-3+Nat) = i+nbonds
            endif
            if (angle_s(j,3)==molecule%geom%angle(i,1).and.&
                angle_s(j,2)==molecule%geom%angle(i,2).and.&
                angle_s(j,1)==molecule%geom%angle(i,3)) then
                print*, "Zmat<->Map Redundant:", j-3+Nat,i+nbonds
                Zmap(j-3+Nat) = i+nbonds
            endif
        enddo
        enddo
        print*, "Dihedral Map"
        do j=4,Nat
        do i=1,ndihed
            if (dihed_s(j,1)==molecule%geom%dihed(i,1).and.&
                dihed_s(j,2)==molecule%geom%dihed(i,2).and.&
                dihed_s(j,3)==molecule%geom%dihed(i,3).and.&
                dihed_s(j,4)==molecule%geom%dihed(i,4)) then
                print*, "Zmat<->Map Redundant:", j-6+2*Nat,i+nbonds+nangles
                Zmap(j-6+2*Nat) = i+nbonds+nangles
            endif
            if (dihed_s(j,4)==molecule%geom%dihed(i,1).and.&
                dihed_s(j,3)==molecule%geom%dihed(i,2).and.&
                dihed_s(j,2)==molecule%geom%dihed(i,3).and.&
                dihed_s(j,1)==molecule%geom%dihed(i,4)) then
                print*, "Zmat<->Map Redundant:", j-6+2*Nat,i+nbonds+nangles
                Zmap(j-6+2*Nat) = i+nbonds+nangles
            endif
        enddo
        enddo
    endif 
    !Set number of redundant
    Nred = molecule%geom%nbonds  + &
           molecule%geom%nangles + &
           molecule%geom%ndihed

    ! SYMMETRIC 
    !Set symmetry of internal (only if symmetry is detected)
    if (adjustl(PG) == "C1") then
        S_sym(3*Nat) = 1
    else
        do i=1,Nat-1
            S_sym(i) = bond_sym(i+1)-1
        enddo
        do i=1,Nat-2
            S_sym(i+Nat-1) = angle_sym(i+2)+Nat-3
        enddo
        do i=1,Nat-3
            S_sym(i+2*Nat-3) = dihed_sym(i+3)+2*Nat-6
        enddo
    endif

    !We send the option -sa within S_sym (confflict with redundant coord!!)
    ! S_sym(3*Nat) =  1 (sa=true) / 0 (sa=false)
    if (symaddapt) then
        S_sym(3*Nat) = 1
    else
        S_sym(3*Nat) = 0
    endif

    !INTERNAL COORDINATES ANALYSIS
    Asel1(1,1) = 99.d0 !out-of-range, as Asel is normalized -- this option is not tested
    call internal_Wilson(molecule,S1,S_sym,ModeDef,B1,G1,Asel1,verbose)
    if (showZ) stop

    !Two possible jobs: scan normal modes or scan a given internal coordinates
    if (scan_internal) then
        print'(A,X,I3,X,A)', "Scanning coordinate", icoord, ModeDef(icoord)
        Nsel=1
        nm(1) = icoord
        L1=0.d0
        L1(icoord,icoord) = 1.d0
        Freq = 1.d0
        Factor(1:Nvib) = 1.d0

     else

        !SOLVE GF METHOD TO GET NM AND FREQ
        !For redundant coordinates a non-redundant set is formed as a combination of
        !the redundant ones. The coefficients for the combination are stored in Asel
        !as they must be used for state 2 (not rederived!).
        call gf_method(Hess,molecule,S_sym,ModeDef,L1,B1,G1,Freq,Asel1,X1,X1inv,verbose) 

        ii=0
        do i=1,Nred
            if (Freq(i) > 1.d-2) then
                ii=ii+1
                Freq(ii) = Freq(i)
                L1(1:Nred,ii) = L1(1:Nred,i)
             endif
        enddo

       !Use freqs. to make displacements equivalent in dimensionless units
        Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib)))/5.d3

    endif

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
    !  Normal mode displacements
    !==========================================================0
    !Initialize auxiliar states
    Nsteps = 100
    Qstep = Amplitude/float(Nsteps)
    if ( mod(Nsteps,2) /= 0 ) Nsteps = Nsteps + 1
    do jj=1,Nsel 
        k=0
        molecule%atom(1:molecule%natoms)%resname = "RES"
        j = nm(jj)
        write(dummy_char,*) j
        if (scan_internal) then
            qcoord=S1(icoord)
!             qcoord=0.d0
            molecule%title = "Animation of internal coordinate "//trim(adjustl(ModeDef(j)))
            g09file="Coord"//trim(adjustl(dummy_char))//"_int.com"
            qfile="Coord"//trim(adjustl(dummy_char))//"_int_steps.dat"
            grofile(jj) = "Coord"//trim(adjustl(dummy_char))//"_int.gro" 
            numfile="Coord"//trim(adjustl(dummy_char))//"_int_num.com"
        else
            qcoord=0.d0
            molecule%title = "Animation of normal mode "//trim(adjustl(grofile(jj)))
            g09file="Mode"//trim(adjustl(dummy_char))//"_int.com"
            qfile="Mode"//trim(adjustl(dummy_char))//"_int_steps.dat"
            grofile(jj) = "Mode"//trim(adjustl(dummy_char))//"_int.gro"
            numfile="Mode"//trim(adjustl(dummy_char))//"_int_num.com"
        endif
        print*, "Writting results to: "//trim(adjustl(grofile(jj)))
        open(O_GRO,file=grofile(jj))
        open(O_G09,file=g09file)
        open(O_Q  ,file=qfile)

        !===========================
        !Start from equilibrium. 
        !===========================
        print*, "STEP:", k
        if (zmat) then
            Smap=S1
        else
            do i=1,Nvib
                ii=Zmap(i)
                Smap(i) = S1(ii)
            enddo
        endif
        call zmat2cart(molecule,bond_s,angle_s,dihed_s,Smap,verbose)
        !Transform to AA and export coords and put back into BOHR
        molecule%atom(1:Nat)%x = molecule%atom(1:Nat)%x*BOHRtoANGS
        molecule%atom(1:Nat)%y = molecule%atom(1:Nat)%y*BOHRtoANGS
        molecule%atom(1:Nat)%z = molecule%atom(1:Nat)%z*BOHRtoANGS
        call write_gro(O_GRO,molecule)
        molec_aux=molecule
        !===========================
        !Half Forward oscillation
        !===========================
        !Initialize distacen criterion for new check_ori2b SR
        dist=0.d0
        do istep = 1,nsteps/2
            k=k+1
            qcoord = qcoord + Qstep/Factor(j)
            print*, "STEP:", k
            write(dummy_char,*) k
            molecule%title = "File generated for "//trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
!             print*, "Bonds"
            do i=1,molecule%geom%nbonds
                S1(i) = S1(i) + L1(i,j) * Qstep/Factor(j)
            enddo
!             print*, "Angles"
            do ii=i,i+molecule%geom%nangles-1
                S1(ii) = S1(ii) + L1(ii,j) * Qstep/Factor(j)
            enddo
!             print*, "Dihedrals"
            do iii=ii,ii+molecule%geom%ndihed-1
               S1(iii) = S1(iii) + L1(iii,j) * Qstep/Factor(j)
               if (S1(iii) > PI) S1(iii)=S1(iii)-2.d0*PI
               if (S1(iii) < -PI) S1(iii)=S1(iii)+2.d0*PI
            enddo
            if (zmat) then
                Smap=S1
            else
                do i=1,Nvib
                    ii=Zmap(i)
                    Smap(i) = S1(ii)
                enddo
            endif
            call zmat2cart(molecule,bond_s,angle_s,dihed_s,Smap,verbose)
            !Transform to AA and comparae with last step (stored in state) -- comparison in AA
            molecule%atom(1:Nat)%x = molecule%atom(1:Nat)%x*BOHRtoANGS
            molecule%atom(1:Nat)%y = molecule%atom(1:Nat)%y*BOHRtoANGS
            molecule%atom(1:Nat)%z = molecule%atom(1:Nat)%z*BOHRtoANGS
            !call check_ori3(state,ref): efficient but not always works. If so, it uses check_ori2(state,ref)
            call check_ori3(molecule,molec_aux,info)
            if (info /= 0 .or. istep==1) then
                call check_ori2b(molecule,molec_aux,dist)
                !The threshold in 5% avobe the last meassured distance
                dist = dist*1.05
            endif
            call write_gro(O_GRO,molecule)
            !Write the max amplitude step to G09 scan
            if (k==nsteps/2) then
                molecule%job%title=trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
                molecule%title=trim(adjustl(g09file))
                call write_gcom(O_G09,molecule)
                write(O_Q,*) qcoord
            endif
            !Save last step in molec_aux (in AA)
            molec_aux=molecule
        enddo
        !=======================================
        ! Reached amplitude. Back oscillation
        !=======================================
        do istep = 1,nsteps/2-3
            k=k+1
            qcoord = qcoord - Qstep/Factor(j)
            print*, "STEP:", k
            write(dummy_char,*) k
            molecule%title = trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
!             print*, "Bonds"
            do i=1,molecule%geom%nbonds
                S1(i) = S1(i) - L1(i,j) * Qstep/Factor(j)
            enddo
!             print*, "Angles"
            do ii=i,i+molecule%geom%nangles-1
                S1(ii) = S1(ii) - L1(ii,j) * Qstep/Factor(j)
            enddo
!             print*, "Dihedrals"
            do iii=ii,ii+molecule%geom%ndihed-1
                S1(iii) = S1(iii) - L1(iii,j) * Qstep/Factor(j)
                if (S1(iii) > PI) S1(iii)=S1(iii)-2.d0*PI
                if (S1(iii) < -PI) S1(iii)=S1(iii)+2.d0*PI
            enddo
            if (zmat) then
                Smap=S1
            else
                do i=1,Nvib
                    ii=Zmap(i)
                    Smap(i) = S1(ii)
                enddo
            endif
            call zmat2cart(molecule,bond_s,angle_s,dihed_s,Smap,verbose)
            molecule%atom(1:Nat)%x = molecule%atom(1:Nat)%x*BOHRtoANGS
            molecule%atom(1:Nat)%y = molecule%atom(1:Nat)%y*BOHRtoANGS
            molecule%atom(1:Nat)%z = molecule%atom(1:Nat)%z*BOHRtoANGS

            call check_ori3(molecule,molec_aux,info)
            if (info /= 0) then
                call check_ori2b(molecule,molec_aux,dist)
                !The threshold in 5% avobe the last meassured distance
                dist = dist*1.05
            endif
            call write_gro(O_GRO,molecule)
            if (mod(k,10) == 0) then
                molecule%job%title=trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
                molecule%title=trim(adjustl(g09file))
                call write_gcom(O_G09,molecule)
                write(O_Q,*) qcoord
            endif
            !Save last step in molec_aux (in AA)
            molec_aux=molecule
        enddo
        !=======================================
        ! Reached equilibrium again (5 points for numerical second derivatives)
        !=======================================
        open(O_NUM,file=numfile,status="replace")
        do istep = 1,5
            k=k+1
            qcoord = qcoord - Qstep/Factor(j)
            print*, "STEP:", k
            write(dummy_char,*) k
            molecule%title=trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
!             print*, "Bonds"
            do i=1,molecule%geom%nbonds
                S1(i) = S1(i) - L1(i,j) * Qstep/Factor(j)
            enddo
!             print*, "Angles"
            do ii=i,i+molecule%geom%nangles-1
                S1(ii) = S1(ii) - L1(ii,j) * Qstep/Factor(j)
            enddo
!             print*, "Dihedrals"
            do iii=ii,ii+molecule%geom%ndihed-1
                S1(iii) = S1(iii) - L1(iii,j) * Qstep/Factor(j)
               if (S1(iii) > PI) S1(iii)=S1(iii)-2.d0*PI
                if (S1(iii) < -PI) S1(iii)=S1(iii)+2.d0*PI
            enddo
            if (zmat) then
                Smap=S1
            else
                do i=1,Nvib
                    ii=Zmap(i)
                    Smap(i) = S1(ii)
                enddo
            endif
            call zmat2cart(molecule,bond_s,angle_s,dihed_s,Smap,verbose)
            molecule%atom(1:Nat)%x = molecule%atom(1:Nat)%x*BOHRtoANGS
            molecule%atom(1:Nat)%y = molecule%atom(1:Nat)%y*BOHRtoANGS
            molecule%atom(1:Nat)%z = molecule%atom(1:Nat)%z*BOHRtoANGS

            call check_ori3(molecule,molec_aux,info)
            if (info /= 0) then
                call check_ori2b(molecule,molec_aux,dist)
                !The threshold in 5% avobe the last meassured distance
                dist = dist*1.05
            endif
            call write_gro(O_GRO,molecule)
            !This time write all five numbers
            molecule%job%title=trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
            molecule%title=trim(adjustl(g09file))
            call write_gcom(O_G09,molecule)
            write(dummy_char,*) qcoord
            molecule%job%title = "Displacement = "//trim(adjustl(dummy_char))
            molecule%title=trim(adjustl(numfile))
            call write_gcom(O_NUM,molecule)
            write(O_Q,*) qcoord
            molec_aux=molecule
        enddo
!         ! Numerical stimation of the second derivative
!         aux_der(2) = (aux_der(3) - aux_der(1))/2.d0/Qstep/Factor(j)
!         aux_der(4) = (aux_der(5) - aux_der(3))/2.d0/Qstep/Factor(j)
!         Freq_num = (aux_der(4) - aux_der(2))/2.d0/Qstep/Factor(j)
!         print*, Freq(j), &
!                 sign(dsqrt(abs(Freq_num)*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
!                      Freq_num)
        !=======================================
        ! Continue Back oscillation
        !=======================================
        do istep = 1,nsteps/2-2
            k=k+1
            qcoord = qcoord - Qstep/Factor(j)
            print*, "STEP:", k
            write(dummy_char,*) k
            molecule%title = trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
!             print*, "Bonds"
            do i=1,molecule%geom%nbonds
                S1(i) = S1(i) - L1(i,j) * Qstep/Factor(j)
            enddo
!             print*, "Angles"
            do ii=i,i+molecule%geom%nangles-1
                S1(ii) = S1(ii) - L1(ii,j) * Qstep/Factor(j)
            enddo
!             print*, "Dihedrals"
            do iii=ii,ii+molecule%geom%ndihed-1
                S1(iii) = S1(iii) - L1(iii,j) * Qstep/Factor(j)
               if (S1(iii) > PI) S1(iii)=S1(iii)-2.d0*PI
                if (S1(iii) < -PI) S1(iii)=S1(iii)+2.d0*PI
            enddo
            if (zmat) then
                Smap=S1
            else
                do i=1,Nvib
                    ii=Zmap(i)
                    Smap(i) = S1(ii)
                enddo
            endif
            call zmat2cart(molecule,bond_s,angle_s,dihed_s,Smap,verbose)
            molecule%atom(1:Nat)%x = molecule%atom(1:Nat)%x*BOHRtoANGS
            molecule%atom(1:Nat)%y = molecule%atom(1:Nat)%y*BOHRtoANGS
            molecule%atom(1:Nat)%z = molecule%atom(1:Nat)%z*BOHRtoANGS

            call check_ori3(molecule,molec_aux,info)
            if (info /= 0) then
                call check_ori2b(molecule,molec_aux,dist)
                !The threshold in 5% avobe the last meassured distance
                dist = dist*1.05
            endif
            call write_gro(O_GRO,molecule)
            if (mod(k,10) == 0) then
                molecule%job%title=trim(adjustl(grofile(jj)))//".step "//trim(adjustl(dummy_char))
                molecule%title=trim(adjustl(g09file))
                call write_gcom(O_G09,molecule)
                write(O_Q,*) qcoord
            endif
            molec_aux=molecule
        enddo
        !=======================================
        ! Reached amplitude. Half Forward oscillation (till equilibrium)
        !=======================================
        do istep = 1,nsteps/2-1
            k=k+1
            print*, "STEP:", k
!             print*, "Bonds"
            do i=1,molecule%geom%nbonds
                S1(i) = S1(i) + L1(i,j) * Qstep/Factor(j)
            enddo
!             print*, "Angles"
            do ii=i,i+molecule%geom%nangles-1
                S1(ii) = S1(ii) + L1(ii,j) * Qstep/Factor(j)
            enddo
!             print*, "Dihedrals"
            do iii=ii,ii+molecule%geom%ndihed-1
                S1(iii) = S1(iii) + L1(iii,j) * Qstep/Factor(j)
               if (S1(iii) > PI) S1(iii)=S1(iii)-2.d0*PI
                if (S1(iii) < -PI) S1(iii)=S1(iii)+2.d0*PI
            enddo
            if (zmat) then
                Smap=S1
            else
                do i=1,Nvib
                    ii=Zmap(i)
                    Smap(i) = S1(ii)
                enddo
            endif
            call zmat2cart(molecule,bond_s,angle_s,dihed_s,Smap,verbose)
            molecule%atom(1:Nat)%x = molecule%atom(1:Nat)%x*BOHRtoANGS
            molecule%atom(1:Nat)%y = molecule%atom(1:Nat)%y*BOHRtoANGS
            molecule%atom(1:Nat)%z = molecule%atom(1:Nat)%z*BOHRtoANGS

            call check_ori3(molecule,molec_aux,info)
            if (info /= 0) then
                call check_ori2b(molecule,molec_aux,dist)
                !The threshold in 5% avobe the last meassured distance
                dist = dist*1.05
            endif
            call write_gro(O_GRO,molecule)
            molec_aux=molecule
        enddo
        close(O_GRO)
        close(O_G09)
        close(O_Q)
    enddo

    if (call_vmd) then
        open(S_VMD,file="vmd_conf.dat",status="replace")
        !Set general display settings (mimic gv)
        write(S_VMD,*) "color Display Background iceblue"
        write(S_VMD,*) "color Name {C} silver"
        write(S_VMD,*) "axes location off"
        !Set molecule representation
        do i=0,Nsel-1
            j = nm(i+1)
            write(S_VMD,*) "mol representation CPK"
!            write(S_VMD,*) "mol addrep 0"
            if (i==0) then
                write(S_VMD,*) "molinfo ", i, " set drawn 1"
            else
                write(S_VMD,*) "molinfo ", i, " set drawn 0"
            endif
            write(S_VMD,*) "mol addrep ", i
            write(dummy_char,'(A,I4,X,F8.2,A)') "{Mode",j, Freq(j),"cm-1}"
            dummy_char=trim(adjustl(dummy_char))
            write(S_VMD,*) "mol rename ", i, trim(dummy_char)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile(i+1)))
        enddo
        write(S_VMD,*) "display projection Orthographic"
        close(S_VMD)
        !Call vmd
        vmdcall = 'vmd -m '
        do i=1,Nsel
        vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile(i)))
        enddo
        vmdcall = trim(adjustl(vmdcall))//" -e vmd_conf.dat"
        call system(vmdcall)
    endif

    if (movie_cycles > 0) then
        open(S_VMD,file="vmd_movie.dat",status="replace")
        !Set general display settings (mimic gv)
        write(S_VMD,*) "color Display Background white"
        write(S_VMD,*) "color Name {C} silver"
        write(S_VMD,*) "axes location off"
        write(S_VMD,*) "display projection Orthographic"
        !Set molecule representation
        do i=0,Nsel-1
            j = nm(i+1)
            write(S_VMD,*) "mol representation CPK"
            write(S_VMD,*) "molinfo ", i, " set drawn 0"
            write(S_VMD,*) "mol addrep ", i
            write(dummy_char,'(A,I4,X,F8.2,A)') "{Mode",j, Freq(j),"cm-1}"
            dummy_char=trim(adjustl(dummy_char))
            write(S_VMD,*) "mol rename ", i, trim(dummy_char)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile(i+1)))
        enddo
        write(S_VMD,'(A)') "#====================="
        write(S_VMD,'(A)') "# Start movies"
        write(S_VMD,'(A)') "#====================="
        !Set length of the movie
        movie_steps = movie_cycles*20
        do i=0,Nsel-1
            j = nm(i+1)
            write(S_VMD,'(A,I4)') "# Mode", j
            write(tmpfile,*) j
            tmpfile="Mode"//trim(adjustl(tmpfile))
            write(S_VMD,*) "molinfo ", i, " set drawn 1"
            write(S_VMD,*) "set figfile "//trim(adjustl(tmpfile))
            write(S_VMD,'(A,I3,A)') "for {set xx 0} {$xx <=", movie_steps,&
                                    "} {incr xx} {"
            write(S_VMD,*) "set x [expr {($xx-($xx/20)*20)*10}]"
            write(S_VMD,*) 'echo "step $x"'
            write(S_VMD,*) "animate goto $x"
            write(S_VMD,*) "render Tachyon $figfile-$xx.dat"
            write(S_VMD,'(A)') '"/usr/local/lib/vmd/tachyon_LINUX" -aasamples 12 '//& 
                           '$figfile-$xx.dat -format TARGA -o $figfile-$xx.tga'
            write(S_VMD,'(A)') 'convert -font URW-Palladio-Roman -pointsize 30 -draw '//&
                           '"text 30,70 '//"'"//trim(adjustl(tmpfile))//"'"//&
                           '" $figfile-$xx.tga $figfile-$xx.jpg'
            write(S_VMD,'(A)') "}"
            !Updated ffmpeg call. The output is now loadable from ipynb
            write(S_VMD,'(A)') 'ffmpeg -i $figfile-%d.jpg -vcodec libx264 -s 640x360 $figfile.mp4'
            write(S_VMD,*) "molinfo ", i, " set drawn 0"
        enddo
        write(S_VMD,*) "exit"
        close(S_VMD)
        !Call vmd
        vmdcall = 'vmd -m '
        do i=1,Nsel
        vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile(i)))
        enddo
        vmdcall = trim(adjustl(vmdcall))//" -e vmd_movie.dat -size 500 500"
        open(O_LIS,file="movie.cmd",status="replace")
        write(O_LIS,'(A)') trim(adjustl(vmdcall))
        write(O_LIS,'(A)') "rm Mode*jpg Mode*dat Mode*tga"
        close(O_LIS)
        print*, ""
        print*, "============================================================"
        print*, "TO GENERATE THE MOVIES (AVI) EXECUTE COMMANDS IN 'movie.cmd'"
        print*, "(you may want to edit 'vmd_movie.dat'  first)"
        print*, "============================================================"
        print*, ""
    endif
    
    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,addfile,nmfile,nm,Nsel,Amplitude,filetype,nosym,zmat,verbose,tswitch,symaddapt,&
                           zmatfile,icoord,showZ,call_vmd,movie_cycles,include_hbonds)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,addfile,nmfile,filetype,zmatfile
        logical,intent(inout) :: nosym, verbose, zmat, tswitch, symaddapt,showZ,call_vmd,include_hbonds
        integer,dimension(:),intent(inout) :: nm
        integer,intent(inout) :: icoord, movie_cycles
        integer,intent(out) :: Nsel
        real(8),intent(out) :: Amplitude
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        ! iargc type must be specified with implicit none (strict compilation)
        integer :: iargc

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

                case ("-nmf") 
                    call getarg(i+1, nmfile)
                    argument_retrieved=.true.

                case ("-maxd") 
                    call getarg(i+1, arg)
                    read(arg,*) Amplitude
                    argument_retrieved=.true.

                case ("-vmd")
                    call_vmd=.true.

                case ("-movie")
                    call getarg(i+1, arg)
                    read(arg,*) movie_cycles
                    argument_retrieved=.true.

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

                case ("-include_hb")
                    include_hbonds=.true.

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
        write(0,*) '-nm            ', nm(1),"-",nm(Nsel)
!         write(0,*) '-nmf           ', nm(1:Nsel)
        write(0,*) '-vmd           ',  call_vmd
        write(0,*) '-movie (cycles)',  movie_cycles
        write(0,*) '-maxd          ',  Amplitude
        if (nosym) dummy_char="NO "
        if (.not.nosym) dummy_char="YES"
        write(0,*) '-[no]sym        ', dummy_char
        if (zmat) dummy_char="YES"
        if (.not.zmat) dummy_char="NO "
        write(0,*) '-[no]zmat       ', dummy_char
        write(0,*) '-readz          ', trim(adjustl(zmatfile))
        write(0,*) '-showz         ', showZ
        if (tswitch) dummy_char="YES"
        if (.not.tswitch) dummy_char="NO "
        write(0,*) '-tswitch        ', dummy_char
        write(0,*) '-include_hb    ',  include_hbonds
        if (symaddapt) dummy_char="YES"
        if (.not.symaddapt) dummy_char="NO "
        write(0,*) '-sa             ', dummy_char
        write(0,*) '-int           ', icoord
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
            case("psi4")
             call read_psi_geom(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select


        return


    end subroutine generic_strfile_read
       

end program normal_modes_animation

