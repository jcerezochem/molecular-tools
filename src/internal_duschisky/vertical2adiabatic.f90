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
    use matrix
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

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: nosym=.true.      ,&
               modred=.false.    ,&
               tswitch=.false.   ,&
               symaddapt=.false. ,&
               vertical=.true.
    character(len=4) :: internal_set='zmat'
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1, state2
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
    real(8),dimension(NDIM,NDIM) :: B1,B2, B, G1,G2
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv,X2,X2inv, L1,L2, Asel1, Asel2, Asel
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    !Duschisky
    real(8),dimension(NDIM,NDIM) :: G
    !T0 - switching effects
    real(8),dimension(3,3) :: T
    !AUXILIAR MATRICES
    real(8),dimension(NDIM,NDIM) :: Aux, Aux2
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: Freq, S1, S2, Vec, Vec2, mu, Factor
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !Shifts
    real(8),dimension(NDIM) :: Delta
    real(8) :: Delta_p
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
    ! RMZ things
    character :: rm_type
    integer :: rm_zline, Nrm, nbonds_rm, nangles_rm, ndiheds_rm
    integer,dimension(100) :: bond_rm, angle_rm, dihed_rm
    ! Yes or no
    logical :: True=.true., False=.false.
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
               I_ADD=14,  &
               I_AD2=15,  &
               I_RMF=16,  &
               O_DUS=20,  &
               O_DIS=21,  &
               O_DMAT=22, &
               O_DUS2=23, &
               O_DIS2=24, &
               O_STAT=25
    !files
    character(len=10) :: filetype="guess", ft, &
                         filetype2="guess"
    character(len=200):: inpfile ="input.fchk",  &
                         inpfile2="none", &
                         addfile="no", &
                         addfile2="no", &
                         zmatfile="NO", &
                         rmzfile="NO", &
                         symm_file="NO", &
                         selfile="modred.dat"
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
    call parse_input(inpfile,inpfile2,filetype,filetype2,nosym,verbose,tswitch,symaddapt,vertical,symm_file,&
                     addfile,addfile2,zmatfile,rmzfile,internal_set,selfile)
    call set_word_upper_case(internal_set)

    ! 1. INTERNAL VIBRATIONAL ANALYSIS ON STATE1 AND STATE2
    

    !===========
    !State 1
    !===========
    print*, ""
    print*, "=========="
    print*, " STATE 1"
    print*, "=========="
 
    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Read structure
    if (adjustl(filetype) == "guess") then
        call split_line_back(inpfile,".",null,filetype)
    endif
    ft=filetype
    call generic_strfile_read(I_INP,ft,state1)
    !Shortcuts
    Nat = state1%natoms
    Nvib = 3*Nat-6

    ! READ HESSIAN

    if (adjustl(ft) == "log") then
        !Gaussian logfile
        allocate(props)
        call parse_summary(I_INP,state1,props,"read_hess")
        !Caution: we NEED to read the Freq summary section
        if (adjustl(state1%job%type) /= "Freq") &
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

    else if (adjustl(ft) == "fchk") then
        !Read Hessian from fchk
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
!         !Read reduced masses
!         call read_fchk(I_INP,"Vib-E2",dtype,N,A,IA,error)
!         mu(1:Nvib) = A(Nvib+1:Nvib+Nvib)
!         deallocate(A)
    else if (adjustl(ft) == "g96") then
        !The hessian should be given as additional input
        if (adjustl(addfile) == "no") &
         call alert_msg("fatal","With a g96, and additional file should be provided with the Hessian")
        open(I_ADD,file=addfile,status="old")
        call read_gro_hess(I_ADD,N,Hess,error)
        close(I_ADD)
    else
        call alert_msg("fatal","Unkown file to get the Hessian: "//ft)
    endif


!     ! READ GRADIENT
    if (trim(adjustl(ft)) == "log") then

          call alert_msg("fatal","GRADIENT from log files still under constructions (perdone las molestias)")

    else if (adjustl(ft) == "fchk") then
        !Read gradient from fchk
        call read_fchk(I_INP,"Cartesian Gradient",dtype,N,A,IA,error)
        Grad(1:N) = A(1:N)
        deallocate(A)
    else
        call alert_msg("note","Gradient could not be read")
    endif
    close(I_INP)

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
    if (adjustl(state1%PG) /= "XX") then
        print'(/,X,A,/)', "Symmetry "//state1%PG 
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
    if (adjustl(internal_set) == "SEL") then
        print*, "Custom internal coordianates"
        open(I_RED,file=selfile,iostat=IOstatus) 
        if (IOstatus /= 0) call alert_msg("fatal","Cannot open file: "//trim(adjustl(selfile)))
        call modredundant(I_RED,state1)
        close(I_RED)
! print*, "Review"
!         do i=1,state1%geom%nimprop
!             i1=state1%geom%improp(i,1)
!             i2=state1%geom%improp(i,2)
!             i3=state1%geom%improp(i,3)
!             i4=state1%geom%improp(i,4)
!             print*, "Improper", i1, i2, i3, i4
!             print*, calc_improper(state1%atom(i1),state1%atom(i2),state1%atom(i3),state1%atom(i4))*180.d0/PI
!         enddo

        ! If the number bonds+angles+dihed+.. is less than Nvib, we are in a reduced space
        ! if so, set Nvib to the new dimension
        Nred = state1%geom%nbonds  + &
               state1%geom%nangles + &
               state1%geom%ndihed
        if (Nvib > Nred) then
            Nvib = Nred
            print*, "Working with a reduced space"
        endif
    elseif (adjustl(internal_set) == "ZMAT") then
        if (adjustl(zmatfile) == "NO") then
            call build_Z(state1,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
        else
            open(I_ZMAT,file=zmatfile,status="old")
            print*, "Z-matrix read from "//trim(adjustl(zmatfile))
            call read_Z(state1,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym,I_ZMAT)
            close(I_ZMAT)
            !Deactivate symaddapt (for the moment)
!             PG = "C1"
        endif
        !Z-mat
        state1%geom%bond(1:Nat-1,1:2) = bond_s(2:Nat,1:2)
        state1%geom%angle(1:Nat-2,1:3) = angle_s(3:Nat,1:3)
        state1%geom%dihed(1:Nat-3,1:4) = dihed_s(4:Nat,1:4)
        state1%geom%nbonds  = Nat-1
        state1%geom%nangles = Nat-2
        state1%geom%ndihed  = Nat-3
        if (rmzfile /= "NO") then
            state1 = state1
            open(I_RMF,file=rmzfile,status='old')
            read(I_RMF,*) Nrm
            ! First, read lines to remove
            nbonds_rm  = 0
            nangles_rm = 0
            ndiheds_rm = 0
            do i=1,Nrm
                read(I_RMF,*) rm_zline, rm_type
                if (rm_type == "B") then
                    nbonds_rm = nbonds_rm + 1
                    bond_rm(nbonds_rm) = rm_zline
                elseif (rm_type == "A") then
                    nangles_rm = nangles_rm + 1
                    angle_rm(nangles_rm) = rm_zline
                elseif (rm_type == "D") then
                    ndiheds_rm = ndiheds_rm + 1
                    dihed_rm(ndiheds_rm) = rm_zline
                endif
            enddo
            ! Short remove lists
            call sort_vec_int(bond_rm(1:nbonds_rm),nbonds_rm,.true.)
            call sort_vec_int(angle_rm(1:nangles_rm),nangles_rm,True)
            call sort_vec_int(dihed_rm(1:ndiheds_rm),ndiheds_rm,True)
            ! The remove
            print*, "Removing bonds from lines:"
            do i=1,nbonds_rm
                print*, bond_rm(i)
                rm_zline = bond_rm(i) - 1
                do ii = rm_zline, state1%geom%nbonds-1
                    state1%geom%bond(ii,1:2) = state1%geom%bond(ii+1,1:2)
                enddo
                state1%geom%nbonds  = state1%geom%nbonds  - 1
                Nvib = Nvib - 1
            enddo
            print*, "Removing angles from lines:"
            do i=1,nangles_rm
                print*, angle_rm(i)
                rm_zline = angle_rm(i) - 2
                do ii = rm_zline, state1%geom%nangles-1
                    state1%geom%angle(ii,1:3) = state1%geom%angle(ii+1,1:3)
                enddo
                state1%geom%nangles  = state1%geom%nangles  - 1
                Nvib = Nvib - 1
            enddo
            print*, "Removing dihedrals from lines:"
            do i=1,ndiheds_rm
                print*, dihed_rm(i)
                rm_zline = dihed_rm(i) - 3
                do ii = rm_zline, state1%geom%ndihed-1
                    state1%geom%dihed(ii,1:4) = state1%geom%dihed(ii+1,1:4)
                enddo
                state1%geom%ndihed  = state1%geom%ndihed  - 1
                Nvib = Nvib - 1
            enddo
            print*, " "
        endif
    else if (adjustl(internal_set) == "ALL") then!otherwise all parameters are used
        print*, "Using all internal coordinates", state1%geom%nangles
        Nred = state1%geom%nbonds  + &
               state1%geom%nangles + &
               state1%geom%ndihed
    else
        call alert_msg("fatal","Unkownn option for -intset. Valid options are 'zmat', 'sel', 'all'")
    endif 

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
    if (symaddapt) then
        S_sym(3*Nat) = 1
    else
        S_sym(3*Nat) = 0
    endif


    ! GET MINIMUM IN CARTESIAN COORDINATES
print*, ""
print*, "GRAD (Cart)"
print*, ""
    print'(5ES16.8)', Grad(1:3*Nat)
print*, ""
    Aux(1:3*Nat,1:3*Nat) = inverse_realsym(3*Nat,Hess)
    Aux2(1:3*Nat,1:3*Nat) = matmul(Hess(1:3*Nat,1:3*Nat),Aux(1:3*Nat,1:3*Nat))
    do i=1,3*Nat
        print'(100F4.1)', Aux2(i,1:3*Nat) 
    enddo
    print*, ""
    ! Matrix x vector 
    do i=1, 3*Nat
        Vec(i)=0.d0
        do k=1,3*Nat
            Vec(i) = Vec(i) - Aux(i,k) * Grad(k)
        enddo 
    enddo
print*, "DELTA R"
do i=1,Nat 
    j=3*i-2
    print'(I3,3X, 3G10.3)', i, Vec(j)*BOHRtoANGS, Vec(j+1)*BOHRtoANGS, Vec(j+2)*BOHRtoANGS
enddo

    state2=state1
    do i=1,Nat 
        j=3*i-2
        state2%atom(i)%x = (state1%atom(i)%x - Vec(j+0))*BOHRtoANGS
        state2%atom(i)%y = (state1%atom(i)%y - Vec(j+1))*BOHRtoANGS
        state2%atom(i)%z = (state1%atom(i)%z - Vec(j+2))*BOHRtoANGS
    enddo
    open(70,file="minim_harmonic_Cart.xyz")
    call write_xyz(70,state2)
    close(70)

    ! INTERNAL COORDINATES


    !SOLVE GF METHOD TO GET NM AND FREQ
    !For redundant coordinates a non-redundant set is formed as a combination of
    !the redundant ones. The coefficients for the combination are stored in Asel
    !as they must be used for state 2 (not rederived!).
    Asel1(1,1) = 99.d0 !out-of-range, as Asel is normalized -- this option is not tested
    call internal_Wilson(state1,Nvib,S1,S_sym,ModeDef,B1,G1,Asel1,verbose)
!     if (vertical) then
        call NumBder(state1,Nvib,S_sym,Bder)
        call gf_method_V(Hess,state1,Nvib,S_sym,ModeDef,L1,B1,G1,Freq,Asel1,X1,X1inv,verbose,Grad=Grad,Bder=Bder) 
!     else
!         call gf_method_V(Hess,state1,Nvib,S_sym,ModeDef,L1,B1,G1,Freq,Asel1,X1,X1inv,verbose)
        !Define the Factor to convert into shift into addimensional displacements
        ! for the shift in SI units:
        Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib))*1.d2*clight*2.d0*PI/plankbar)
        ! but we have it in au
        Factor(1:Nvib)=Factor(1:Nvib)*BOHRtoM*dsqrt(AUtoKG)
!     endif 

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
                    state1%atom(k)%mass/UMAtoAU * &
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
    print'(F15.8)', theta
    enddo
    print*, ""
    !Print state
    open(O_STAT,file="state_file_1")
    do i=1,Nat
        write(O_STAT,*) state1%atom(i)%x*BOHRtoANGS
        write(O_STAT,*) state1%atom(i)%y*BOHRtoANGS
        write(O_STAT,*) state1%atom(i)%z*BOHRtoANGS
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

    if (verbose) then
    print*, "B1=", B1(1,1)
    do i=1,Nvib
        print'(100(F8.3,2X))', B1(i,1:Nvib)
    enddo
    endif


    ! GET MINIMUM IN INTERNAL COORDINATES
    print*, ""
    print*, "Diagonal FC (internal)"
    print*, ""
    do i=1,Nvib
        print'(1ES16.8)', Hess(i,i)
    enddo
    Aux(1:Nvib,1:Nvib) = inverse_realgen(Nvib,Hess)
!     Aux2(1:Nvib,1:Nvib) = matmul(Hess(1:Nvib,1:Nvib),Aux(1:Nvib,1:Nvib))
!     do i=1,Nvib
!         print'(100F5.2)', Aux2(i,1:Nvib) 
!     enddo
!     print*, ""
!     do i=1,Nvib
!         print'(100F5.2)', Hess(i,1:Nvib) 
!     enddo
!     print*, ""
    ! Matrix x vector 
    print*, "" 
    print*, "GRAD (internal)"
    print*, ""
    print'(5ES16.8)', Grad(1:Nvib)
    print*, ""
    print*, ""
    print*, "Diagonal FC (internal)"
    print*, ""
    do i=1,Nvib
        print'(1ES16.8)', Hess(i,i)
    enddo
    print*, ""
    do i=1,Nvib
        Vec(i)=0.d0
        do k=1,Nvib
            Vec(i) = Vec(i)-Aux(i,k) * Grad(k)
        enddo 
    enddo
    k=0
    print*, "DELTA BONDS"
    do i=1,state1%geom%nbonds
        k = k+1
        print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Vec(k), Vec(k)*BOHRtoANGS, Grad(k)
    enddo
    print*, "DELTA ANGLES"
    do i=1,state1%geom%nangles
        k = k+1
        print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Vec(k), Vec(k)*180.d0/PI, Grad(k)
    enddo
    print*, "DELTA DIHEDRALS"
    do i=1,state1%geom%ndihed
        k = k+1
        print'(A,3X,2(F8.3,3X),G10.3)', trim(adjustl(ModeDef(k))), Vec(k), Vec(k)*180.d0/PI, Grad(k)
    enddo
    do i=1,Nvib
        S1(i) = S1(i) + Vec(i)
    enddo
    call zmat2cart(state1,bond_s,angle_s,dihed_s,S1,verbose)
    !Transform to AA and export coords and put back into BOHR
    state1%atom(1:Nat)%x = state1%atom(1:Nat)%x*BOHRtoANGS
    state1%atom(1:Nat)%y = state1%atom(1:Nat)%y*BOHRtoANGS
    state1%atom(1:Nat)%z = state1%atom(1:Nat)%z*BOHRtoANGS
    open(70,file="minim_harmonic_Inter.xyz")
    call write_xyz(70,state1)
    close(70)

    print*, ""


    call cpu_time(tf)
    write(6,'(A,F12.3)') "CPU (s) for internal vib analysis: ", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,inpfile2,filetype,filetype2,nosym,verbose,tswitch,symaddapt,vertical,symm_file,&
                           addfile,addfile2,zmatfile,rmzfile,internal_set,selfile)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,inpfile2,filetype,filetype2,zmatfile,symm_file,&
                                          addfile, addfile2,rmzfile, internal_set, selfile
        logical,intent(inout) :: nosym, verbose, tswitch, symaddapt, vertical
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

                case ("-f2") 
                    call getarg(i+1, inpfile2)
                    argument_retrieved=.true.
                case ("-ft2") 
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.

                case ("-add") 
                    call getarg(i+1, addfile)
                    argument_retrieved=.true.

                case ("-add2") 
                    call getarg(i+1, addfile2)
                    argument_retrieved=.true.

                case ("-nosym")
                    nosym=.true.
                case ("-sym")
                    nosym=.false.

                case ("-symfile")
                    nosym=.false.
                    call getarg(i+1, symm_file)
                    argument_retrieved=.true.

                case ("-sa")
                    symaddapt=.true.
                case ("-nosa")
                    symaddapt=.false.

                case ("-zfile") 
                    call getarg(i+1, zmatfile)
                    argument_retrieved=.true.

                case ("-sfile") 
                    call getarg(i+1, selfile)
                    argument_retrieved=.true.

                case ("-rmz") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.

                case ("-intset")
                    call getarg(i+1, internal_set)
                    argument_retrieved=.true.

                case ("-tswitch")
                    tswitch=.true.

                case ("-vert")
                     vertical=.true.

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
        write(6,*) '-add            ', trim(adjustl(addfile))
        write(6,*) '-f2             ', trim(adjustl(inpfile2))
        write(6,*) '-add2           ', trim(adjustl(addfile2))
        if (nosym) dummy_char="NO "
        if (.not.nosym) dummy_char="YES"
        write(6,*) '-[no]sym        ', dummy_char
        write(6,*) '-symfile        ', trim(adjustl(symm_file))
        write(6,*) '-intset         ', trim(adjustl(internal_set))
        write(6,*) '-zfile          ', trim(adjustl(zmatfile))
        write(6,*) '-sfile          ', trim(adjustl(selfile))
        write(6,*) '-rmz            ', trim(adjustl(rmzfile))
        if (tswitch) dummy_char="YES"
        if (.not.tswitch) dummy_char="NO "
        write(6,*) '-tswitch        ', dummy_char
        if (symaddapt) dummy_char="YES"
        if (.not.symaddapt) dummy_char="NO "
        write(6,*) '-sa             ', dummy_char
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

        ! Predefined filetypes
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
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

        return


    end subroutine generic_strfile_read
       

end program internal_duschinski

