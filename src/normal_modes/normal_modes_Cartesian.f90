program normal_modes_cartesian


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS 
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to visualize vibrations obtained in internal coordinates.
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
    use gro_manage
    use gaussian_manage
    !============================================
    !  Structure-related modules
    !============================================
    use molecular_structure
    use ff_build
    use atomic_geom
    use symmetry
    !============================================
    !  Vibrational 
    !============================================
    use internal_module
    use zmat_manage 
    use vibrational_analysis
    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: use_symmetry=.false.,    &
               include_hbonds=.false.,  &
               vertical=.false.,        &
               analytic_Bder=.false.,   &
               check_symmetry=.true.,   &
               full_diagonalize=.false.,&
               animate=.true.
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule, molec_aux
    type(str_bonded) :: zmatgeom
    real(8),dimension(NDIM) :: X0, Y0, Z0
    integer,dimension(1:NDIM) :: isym
    integer,dimension(4,1:NDIM,1:NDIM) :: Osym
    integer :: nsym
    integer :: Nat, Nvib, Ns
    character(len=5) :: PG
    !Job info
    character(len=20) :: calc, method, basis
    character(len=150):: title
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=50) :: dummy_char
    real(8) :: dist
    !io flags
    integer :: error, info
    ! other flags
    integer :: verbose_current
    ! Auxiliar real arrays/scalars
    real(8),dimension(1:NDIM,1:NDIM) :: Aux, Aux2
    real(8),dimension(:),allocatable :: Vec_alloc
    real(8) :: Theta, Theta2
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii, jj, kk, iop
    !=============

    !====================== 
    ! PES topology and normal mode things
    real(8),dimension(1:NDIM,1:NDIM) :: LL
    real(8),dimension(1:NDIM*NDIM)   :: Hlt
    real(8),dimension(:,:),allocatable :: Hess
    real(8),dimension(NDIM) :: Freq, Factor, Grad
    !Moving normal modes
    character(len=50) :: selection="none"
    real(8) :: Amplitude = 2.d0, qcoord
    integer,dimension(1:NDIM) :: nm=0
    real(8) :: Qstep
    logical :: call_vmd = .false.
    character(len=10000) :: vmdcall
    integer :: Nsteps, Nsel=0, istep
    !MOVIE things
    logical :: movie_vmd = .false.
    integer :: movie_cycles=0,& !this means no movie
               movie_steps
    !====================== 

    !====================== 
    !INTERNAL CODE THINGS
    real(8),dimension(1:NDIM,1:NDIM) :: B, G, Asel
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    real(8),dimension(1:NDIM,1:NDIM) :: X,Xinv
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: S, Sref, S0
    integer,dimension(NDIM) :: S_sym
    ! Switches
    character(len=5) :: def_internal="ZMAT", def_internal_aux
    character(len=2) :: scan_type="NM"
    !Coordinate map
    integer,dimension(NDIM) :: Zmap
    ! Number of ic (Shortcuts)
    integer :: nbonds, nangles, ndihed, nimprop
    !====================== 

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_SYM=12,  &
               I_RMF=16,  &
               I_CNX=17,  &
               O_GRO=20,  &
               O_G09=21,  &
               O_G96=22,  &
               O_Q  =23,  &
               O_NUM=24,  &
               O_MOV=25,  &
               S_VMD=30

    !files
    character(len=10) :: ft ="guess",  ftg="guess",  fth="guess", ftn="guess"
    character(len=200):: inpfile  ="state1.fchk", &
                         gradfile ="same", &
                         hessfile ="same", &
                         nmfile   ="none", &
                         intfile  ="none", &
                         rmzfile  ="none", &
                         symm_file="none", &
                         cnx_file="guess"
    !Structure files to be created
    character(len=100) :: g09file,qfile, tmpfile, g96file, grofile,numfile
    !status
    integer :: IOstatus
    !===================

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

    call cpu_time(ti)

    !--------------------------
    ! Tune io
    !--------------------------
    ! Set unit for alert messages
    alert_unt=6
    ! Activate notes
    silent_notes = .false.
    !--------------------------

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(&
                     ! input data
                     inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn,      &
                     ! Options (general)
                     Amplitude,call_vmd,include_hbonds,selection,vertical, &
                     ! Options (Cartesian)
                     full_diagonalize,                                     &
                     ! Animation and Movie
                     animate,movie_vmd, movie_cycles ,                     &   
                     ! Options (internal)
                     use_symmetry,def_internal,intfile,rmzfile,            & !except scan_type
                     ! connectivity file
                     cnx_file,                                             &
                     ! (hidden)
                     analytic_Bder)

 
    ! 1. READ DATA
    ! ---------------------------------
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)
    if (fth == "guess") &
    call split_line_back(hessfile,".",null,fth)
    if (ftg == "guess") &
    call split_line_back(gradfile,".",null,ftg)
    if (ftn == "guess") &
    call split_line_back(nmfile,".",null,ftn)

    ! Manage special files (fcc) 
    if (adjustl(ft) == "fcc" .or. adjustl(ftn) == "fcc") then
        call alert_msg("note","fcc files needs fcc-input as -f and statefile as -ftn")
        ft ="fcc"
        ftn="fcc"
        ! inpfile has Nat, Nvib, and Masses          <= in inpfile
        ! statefile has coordinates and normal modes <= in hessfile
        ! Generic generic readers parse the state (not the inpfile)
        ! so we get the info here
        open(I_INP,file=inpfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
        read(I_INP,*) Nat 
        molecule%natoms = Nat
        read(I_INP,*) Nvib
        do i=1,Nat 
            read(I_INP,*) molecule%atom(i)%mass
            !Set atomnames from atommasses
            call atominfo_from_atmass(molecule%atom(i)%mass,  &
                                      molecule%atom(i)%AtNum, &
                                      molecule%atom(i)%name)
        enddo
        close(I_INP)
        ! Now put the statefile in the inpfile
        inpfile=nmfile
    elseif (adjustl(ftn) == "log") then
        ! Need to read the standard orientation, not from summary section
        ft = "log-stdori"
    endif
        
    ! STRUCTURE FILE
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,molecule,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (State1)")
    ! Get job info if it is a Gaussian file
    if (ft == "log" .or. ft== "fchk") then
        rewind(I_INP) ! this should be a generic_job_rewind() call
        call read_gauss_job(I_INP,ft,calc,method,basis)
        ! Whichever, calc type was, se now need SP
        calc="SP"
    else
        calc="SP"
        method="B3LYP"
        basis="6-31G(d)"
    endif
    close(I_INP)
    ! Shortcuts
    Nat = molecule%natoms


    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! Get connectivity 
    if (cnx_file == "guess") then
        call guess_connect(molecule)
    else
        print'(/,A,/)', "Reading connectivity from file: "//trim(adjustl(cnx_file))
        open(I_CNX,file=cnx_file,status='old')
        call read_connect(I_CNX,molecule)
        close(I_CNX)
    endif


    ! VIBRATIONAL ANALYSIS: either read from file or from diagonalization of Hessian
    if (adjustl(nmfile) /= "none") then
        open(I_INP,file=nmfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(nmfile)) )
        call generic_nm_reader(I_INP,ftn,Nat,Nvib,Freq,LL)
        ! Show frequencies
        if (verbose>0) &
         call print_vector(6,Freq,Nvib,"Frequencies (cm-1)")
        ! The reader provide L in Normalized Cartesian. Need to Transform to Cartesian now
        call LcartNrm_to_Lmwc(Nat,Nvib,molecule%atom(:)%mass,LL,LL)
        call Lmwc_to_Lcart(Nat,Nvib,molecule%atom(:)%mass,LL,LL,error)
        close(I_INP)

    else
        ! ACTUALLY PERFORM THE ANALYSIS
        ! HESSIAN FILE
        open(I_INP,file=hessfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile)) )
        call generic_Hessian_reader(I_INP,fth,Nat,Hlt,error)
        if (error /= 0) call alert_msg("fatal","Error reading Hessian (State1)")
        close(I_INP)
        
        ! GRADIENT FILE
        if (vertical) then
            open(I_INP,file=gradfile,status='old',iostat=IOstatus)
            if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
            call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
            close(I_INP)
        endif


        if (vertical) then
            ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
            print*, "Preliminary vibrational analysis"
            call vibrations_Cart(Nat,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,molecule%atom(:)%Mass,Hlt,&
                                 Nvib,LL,Freq,error_flag=error)
            ! Manage symmetry
            if (.not.use_symmetry) then
                molecule%PG="C1"
            else if (trim(adjustl(symm_file)) /= "none") then
                call alert_msg("note","Using custom symmetry file: "//trim(adjustl(symm_file)) )
                open(I_SYM,file=symm_file)
                do i=1,molecule%natoms
                    read(I_SYM,*) j, isym(j)
                enddo
                close(I_SYM)
                !Set PG to CUStom
                molecule%PG="CUS"
            else
                molecule%PG="XX"
                call symm_atoms(molecule,isym)
            endif
            
            !Generate bonded info
            call gen_bonded(molecule)
            
            !From now on work in au
            call set_geom_units(molecule,"Bohr")

            ! Define internal set
            if (def_internal=="SEL".or.def_internal=="ALL".or.(def_internal=="ZMAT".and.rmzfile/="none")) then 
                !Only if def_internal="all", we can print the animation mapping the Zmat
                !but NOT for "sel"
                print*, "Preliminary Zmat analysis"
                ! Get Zmat first
                def_internal_aux="ZMAT"
                call define_internal_set(molecule,def_internal_aux,"none","none",use_symmetry,isym,S_sym,Ns)
                ! Get only the geom, and reuse molecule
                zmatgeom=molecule%geom
                ! Compute S values (to get initial values for frozen (rmzfile) ICs)
                verbose_current=verbose
!                 verbose=0
                call set_geom_units(molecule,"Bohr")
                call internal_Wilson_new(molecule,Ns,S0,B)
                call set_geom_units(molecule,"Angs")
                verbose=verbose_current
                ! And reset bonded parameters
                call gen_bonded(molecule)
            endif
            call define_internal_set(molecule,def_internal,intfile,rmzfile,use_symmetry,isym,S_sym,Ns)
            if (Ns > Nvib) then
                call red2zmat_mapping(molecule,zmatgeom,Zmap)
            elseif (def_internal=="ZMAT".and.rmzfile/="none") then
                ! We also get a Zmap
                call red2zmat_mapping(molecule,zmatgeom,Zmap)
                Nvib=Ns
            elseif (Ns < Nvib) then
                print*, "Ns", Ns
                print*, "Nvib", Nvib
                call alert_msg("warning","Reduced coordinates only produce animations with rmzfiles")
                ! Need to freeze unused coords to its input values
            endif

            print'(X,A,/)', "Apply correction for non-stationary points"
            ! Get metric matrix and Bders
            call internal_Wilson_new(molecule,Ns,S,B,ModeDef)
            call internal_Gmetric(Nat,Ns,molecule%atom(:)%mass,B,G)
            call calc_Bder(molecule,Ns,Bder,analytic_Bder)
            if (Ns > Nvib) then
                call redundant2nonredundant(Ns,Nvib,G,Asel)
                ! Rotate Bmatrix
                B(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel,B,tA=.true.)
                ! Rotate Gmatrix
                G(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel(1:Ns,1:Nvib),G,counter=.true.)
                ! Rotate Bders
                do j=1,3*Nat
                    Bder(1:Nvib,j,1:3*Nat) =  matrix_product(Nvib,3*Nat,Ns,Asel,Bder(1:Ns,j,1:3*Nat),tA=.true.)
                enddo
            endif
            ! Get gradient in internal coordinates
            call Gradcart2int(Nat,Nvib,Grad,molecule%atom(:)%mass,B,G)

            ! Compute gs^t * Lder
            do i=1,3*Nat
            do j=1,3*Nat
                Aux(i,j) = 0.d0
                do k=1,Nvib
                    Aux(i,j) = Aux(i,j) + Bder(k,i,j)*Grad(k)
                enddo
            enddo
            enddo

            ! Compute modified Hessian (Hess - gs^t beta) directly on Hlt
            allocate(Vec_alloc(1:(3*Nat*(3*Nat+1))/2))
            Vec_alloc(1:(3*Nat*(3*Nat+1))/2) = Hess_to_Hlt(3*Nat,Aux)
            Hlt(1:(3*Nat*(3*Nat+1))/2) = Hlt(1:(3*Nat*(3*Nat+1))/2) - Vec_alloc(1:(3*Nat*(3*Nat+1))/2)
            deallocate(Vec_alloc)

            ! Also check the symmetry of the correction term
            if (check_symmetry) then
                print'(/,X,A)', "---------------------------------------"
                print'(X,A  )', " Check effect of symmetry operations"
                print'(X,A  )', " on the correction term gs^t\beta"
                print'(X,A  )', "---------------------------------------"
                molecule%PG="XX"
                call symm_atoms(molecule,isym,Osym,rotate=.false.,nsym_ops=nsym)
                ! Check the symmetry of the correction term
                ! First compute the correction term (was already computed)
                Aux2(1:3*Nat,1:3*Nat) = Aux(1:3*Nat,1:3*Nat)
                ! Print if verbose level is high
                if (verbose>2) &
                    call MAT0(6,Aux2,3*Nat,3*Nat,"gs*Bder matrix")
                ! Check all detected symmetry ops
                do iop=1,Nsym
                    Aux(1:3*Nat,1:3*Nat) = dfloat(Osym(iop,1:3*Nat,1:3*Nat))
                    Aux(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,Aux,Aux2,counter=.true.)
                    Theta=0.d0
                    do i=1,3*Nat 
                    do j=1,3*Nat
                        if (Theta < abs(Aux(i,j)-Aux2(i,j))) then
                            Theta = abs(Aux(i,j)-Aux2(i,j))
                            Theta2=Aux2(i,j)
                        endif
                    enddo
                    enddo
                    print'(X,A,I0)', "Symmetry operation :   ", iop
                    print'(X,A,F10.6)',   " Max abs difference : ", Theta
                    print'(X,A,F10.6,/)', " Value before sym op: ", Theta2
                enddo
                print'(X,A,/)', "---------------------------------------"
            endif

            !Reset Angs
            call set_geom_units(molecule,"Angs")

        endif

        ! VIBRATIONAL ANALYSIS
        ! Using either the full 3Natx3Nat space of internal frame (default)
        if (full_diagonalize) then
            ! Diag in 3Natx3Nat
            allocate(Hess(1:3*Nat,1:3*Nat))
            !Massweight the Hessian
            k=0
            do i=1,3*Nat
            do j=1,i
                k=k+1
                ii = (i-1)/3+1
                jj = (j-1)/3+1
                Hess(i,j) = Hlt(k)/sqrt(molecule%atom(ii)%mass * &
                                        molecule%atom(jj)%mass) / AMUtoAU
                Hess(j,i) = Hess(i,j)
            enddo 
            enddo
            ! Make Nvib equal to full space dimension
            Nvib = 3*Nat
            call diagonalize_full(Hess(1:Nvib,1:Nvib),Nvib,LL(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
            deallocate(Hess)
            !Check FC
            if (verbose>1) &
             call print_vector(6,Freq*1.d6,Nvib,"FORCE CONSTANTS x 10^6 (A.U.)")
            !Transform to FC to Freq
            do i=1,Nvib
                  Freq(i) = sign(dsqrt(abs(Freq(i))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
                                 Freq(i))
            enddo
            if (verbose>0) &
                call print_vector(6,Freq,Nvib,"Frequencies (cm-1)")
        else
            ! Run final vibrations_Cart 
            call vibrations_Cart(Nat,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,molecule%atom(:)%Mass,Hlt,&
                                 Nvib,LL,Freq,error_flag=error)
            if (error /= 0) then
                call alert_msg("fatal","Error in the diagonalization")
                stop
            endif
        endif
        !Transform L to Cartesian (output in AU(mass) as with internal)
        call Lmwc_to_Lcart(Nat,Nvib,molecule%atom(1:Nat)%mass,LL,LL)
    endif

    !Define the Factor to convert shift into addimensional displacements
    ! from the shift in SI units:
    Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib))*1.d2*clight*2.d0*PI/plankbar)
    ! but we need it from au not SI
    Factor(1:Nvib)=Factor(1:Nvib)*BOHRtoM*dsqrt(AUtoKG)

    ! Check if we should continue
    if (.not.animate) then
        call cpu_time(tf)
        write(6,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti
        stop
    endif

    ! Get the selection of normal modes to represent
    if (adjustl(selection) == "all") then
        Nsel=Nvib
        nm(1:Nvib) = (/(i, i=1,Nvib)/)
    else if (adjustl(selection) == "none") then
        Nsel=0
    else
        call selection2intlist(selection,nm,Nsel)
    endif

    !==========================================================0
    !  Normal mode displacements
    !==========================================================0
    ! Initialization
    X0(1:Nat) = molecule%atom(1:Nat)%x
    Y0(1:Nat) = molecule%atom(1:Nat)%y
    Z0(1:Nat) = molecule%atom(1:Nat)%z
    Nsteps = 101
    if ( mod(Nsteps,2) == 0 ) Nsteps = Nsteps + 1 ! ensure odd number of steps (so we have same left and right)
    ! Qstep is dimless
    Qstep = Amplitude/float(Nsteps-1)*2.d0  ! Do the range (-A ... +A)
    molecule%atom(1:molecule%natoms)%resname = "RES" ! For printing
    ! Run over all selected modes/internals
    do jj=1,Nsel 
        k=0 ! equilibrium corresponds to k=0
        j = nm(jj)
        if (verbose>0) &
         print'(X,A,I0,A)', "Generating Mode ", j, "..."

        ! Set initial values for the scanned coordinate
        qcoord = 0.d0 

        ! Prepare and open files
        call prepare_files(j,grofile,g09file,g96file,numfile,qfile,title,full_diagonalize)
        open(O_GRO,file=grofile)
        open(O_G09,file=g09file)
        open(O_G96,file=g96file)
        open(O_Q  ,file=qfile)
        open(O_NUM,file=numfile)

        !===========================
        !Start from equilibrium. 
        !===========================
        molecule%atom(1:Nat)%x = X0(1:Nat)
        molecule%atom(1:Nat)%y = Y0(1:Nat)
        molecule%atom(1:Nat)%z = Z0(1:Nat)
        if (verbose>1) &
         print'(/,A,I0)', "STEP:", k
        ! Update
        write(molecule%title,'(A,I0,A,2(X,F12.6))') &
         trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
        ! Print initial structure
        call write_gro(O_GRO,molecule)
        !===========================
        !Half Forward oscillation: from step "Eq + dQ" to "Eq + (N-1)/2 dQ"
        !===========================
        !Initialize distacen criterion for rmsd_fit_frame_brute SR
        dist=0.d0
        do istep = 1,(nsteps-1)/2
            k=k+1
            if (verbose>1) &
             print'(/,A,I0)', "STEP:", k
            ! Update values
            ! qcoord has AU
            qcoord = qcoord + Qstep/Factor(j)
            write(molecule%title,'(A,I0,A,2(X,F12.6))') &
             trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
            ! Displace in AA (displace always from reference to avoid error propagation)
            i=istep
            molecule%atom(1:Nat)%x = X0(1:Nat)
            molecule%atom(1:Nat)%y = Y0(1:Nat)
            molecule%atom(1:Nat)%z = Z0(1:Nat)
            call displace_Xcoord(LL(1:3*Nat,j),molecule%natoms,Qstep/Factor(j)*BOHRtoANGS*i,&
                                 molecule%atom(:)%x,molecule%atom(:)%y,molecule%atom(:)%z)
            ! PRINT
            ! Write GRO from the beginign and G96/G09 only when reach max amplitude
            call write_gro(O_GRO,molecule)
            if (k==(nsteps-1)/2) then
                call write_gcom(O_G09,molecule,g09file,calc,method,basis,molecule%title)
                call write_g96(O_G96,molecule)
                write(O_Q,*) qcoord, qcoord*Factor(j)
            endif
        enddo
        !=======================================
        ! Reached amplitude. Back oscillation: from step "MaxAmp + dQ" to "MaxAmp + (N-2) dQ" == -MaxAmp
        !=======================================
        ! This is the part reported in G09/G96 files
        do istep = 1,nsteps-1
            k=k+1
            if (verbose>1) &
             print'(/,A,I0)', "STEP:", k
            ! Update values
            qcoord = qcoord - Qstep/Factor(j)
            write(molecule%title,'(A,I0,A,2(X,F12.6))') &
             trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
            ! Displace in AA (displace always from reference to avoid error propagation)
            i=(nsteps-1)/2-istep
            molecule%atom(1:Nat)%x = X0(1:Nat)
            molecule%atom(1:Nat)%y = Y0(1:Nat)
            molecule%atom(1:Nat)%z = Z0(1:Nat)
            call displace_Xcoord(LL(1:3*Nat,j),molecule%natoms,Qstep/Factor(j)*BOHRtoANGS*i,&
                                 molecule%atom(:)%x,molecule%atom(:)%y,molecule%atom(:)%z)
            ! PRINT
            ! Write G96/GRO every step and G09 scan every 10 steps
            ! except the 5 poinst around minimum, which are all printed
            call write_gro(O_GRO,molecule)
            call write_g96(O_G96,molecule)
            if (mod(k,10) == 0) then
                call write_gcom(O_G09,molecule,g09file,calc,method,basis,molecule%title)
                write(O_Q,*) qcoord, qcoord*Factor(j)
            endif
            ! Write 5 poinst around minimum for numerical dierivatives
            if (k>=nsteps-3.and.k<=nsteps+1) then
                call write_gcom(O_NUM,molecule,numfile,calc,method,basis,molecule%title)
            endif
        enddo
        !=======================================
        ! Reached amplitude. Half Forward oscillation (till one step before equilibrium, so that we concatenate well)
        !=======================================
        do istep = 1,(nsteps-1)/2-1
            k=k+1
            if (verbose>1) &
             print'(/,A,I0)', "STEP:", k
            ! Update values
            qcoord = qcoord + Qstep/Factor(j)
            write(molecule%title,'(A,I0,A,2(X,F12.6))') &
             trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
            ! Displace in AA (displace always from reference to avoid error propagation)
            i=-(nsteps-1)/2+istep
            molecule%atom(1:Nat)%x = X0(1:Nat)
            molecule%atom(1:Nat)%y = Y0(1:Nat)
            molecule%atom(1:Nat)%z = Z0(1:Nat)
            call displace_Xcoord(LL(1:3*Nat,j),molecule%natoms,Qstep/Factor(j)*BOHRtoANGS*i,&
                                 molecule%atom(:)%x,molecule%atom(:)%y,molecule%atom(:)%z)
            ! PRINT
            ! Write only GRO 
            call write_gro(O_GRO,molecule)
        enddo
        open(O_GRO)
        open(O_G09)
        open(O_G96)
        open(O_Q  )
        open(O_NUM)
    enddo

    call summary_alerts

    call cpu_time(tf)
    write(6,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    ! CALL EXTERNAL PROGRAM TO RUN ANIMATIONS

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
        enddo
        write(S_VMD,*) "display projection Orthographic"
        close(S_VMD)
        !Call vmd
        vmdcall = 'vmd -m '
        do i=1,Nsel
            j = nm(i)
            ! Get filenames (we want grofile name)
            call prepare_files(j,grofile,g09file,g96file,numfile,qfile,title,full_diagonalize)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile))
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
            ! Get filenames (we want grofile name)
            call prepare_files(j,grofile,g09file,g96file,numfile,qfile,title,full_diagonalize)
            write(S_VMD,*) "mol representation CPK"
            write(S_VMD,*) "molinfo ", i, " set drawn 0"
            write(S_VMD,*) "mol addrep ", i
            write(dummy_char,'(A,I4,X,F8.2,A)') "{Mode",j, Freq(j),"cm-1}"
            dummy_char=trim(adjustl(dummy_char))
            write(S_VMD,*) "mol rename ", i, trim(dummy_char)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile))
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
        vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile))
        enddo
        vmdcall = trim(adjustl(vmdcall))//" -e vmd_movie.dat -size 500 500"
        open(O_MOV,file="movie.cmd",status="replace")
        write(O_MOV,'(A)') trim(adjustl(vmdcall))
        write(O_MOV,'(A)') "rm Mode*jpg Mode*dat Mode*tga"
        close(O_MOV)
        print*, ""
        print*, "============================================================"
        print*, "TO GENERATE THE MOVIES (AVI) EXECUTE COMMANDS IN 'movie.cmd'"
        print*, "(you may want to edit 'vmd_movie.dat'  first)"
        print*, "============================================================"
        print*, ""
    endif

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(&
                           ! input data
                           inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn,      &
                           ! Options (general)
                           Amplitude,call_vmd,include_hbonds,selection,vertical, &
                           ! Options (Cartesian)
                           full_diagonalize,                                     &
                           ! Movie
                           animate,movie_vmd, movie_cycles,                      &
                           ! Options (internal)
                           use_symmetry,def_internal,intfile,rmzfile,            &
                           ! connectivity file
                           cnx_file,                                             &
                           ! (hidden)
                           analytic_Bder)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn,selection, &
                                          !Internal
                                          def_internal,intfile,rmzfile,cnx_file
        real(8),intent(inout)          :: Amplitude
        logical,intent(inout)          :: call_vmd, include_hbonds,vertical,movie_vmd,full_diagonalize,animate,&
                                          ! Internal
                                          use_symmetry,analytic_Bder
        integer,intent(inout)          :: movie_cycles

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
                    call getarg(i+1, ft)
                    argument_retrieved=.true.

                case ("-fhess") 
                    call getarg(i+1, hessfile)
                    argument_retrieved=.true.
                case ("-fth") 
                    call getarg(i+1, fth)
                    argument_retrieved=.true.

                case ("-fgrad") 
                    call getarg(i+1, gradfile)
                    argument_retrieved=.true.
                case ("-ftg") 
                    call getarg(i+1, ftg)
                    argument_retrieved=.true.

                case ("-fnm") 
                    call getarg(i+1, nmfile)
                    argument_retrieved=.true.
                case ("-ftn") 
                    call getarg(i+1, ftn)
                    argument_retrieved=.true.

                case ("-fulldiag")
                    full_diagonalize=.true.
                case ("-nofulldiag")
                    full_diagonalize=.false.

                case ("-nm") 
                    call getarg(i+1, selection)
                    argument_retrieved=.true.

                case ("-disp") 
                    call getarg(i+1, arg)
                    read(arg,*) Amplitude
                    argument_retrieved=.true.

                case ("-vmd")
                    call_vmd=.true.
                case ("-novmd")
                    call_vmd=.false.

                case ("-animate")
                    animate=.true.
                case ("-noanimate")
                    animate=.false.

                case ("-movie")
                    call getarg(i+1, arg)
                    read(arg,*) movie_cycles
                    movie_vmd=.true.
                    argument_retrieved=.true.

                case ("-include_hb")
                    include_hbonds=.true.

                ! Internal (for correction) options
                case ("-vert")
                    vertical=.true.
                case ("-novert")
                    vertical=.false.

                case ("-sym")
                    use_symmetry=.true.
                case ("-nosym")
                    use_symmetry=.false.

                case ("-cnx") 
                    call getarg(i+1, cnx_file)
                    argument_retrieved=.true.

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-rmzfile") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.

                case ("-intmode")
                    call getarg(i+1, def_internal)
                    argument_retrieved=.true.

                ! (HIDDEN FLAG)
                case ("-anaBder")
                    analytic_Bder=.true.
                case ("-noanaBder")
                    analytic_Bder=.false.

        
                case ("-h")
                    need_help=.true.

                ! Control verbosity
                case ("-quiet")
                    verbose=0
                    silent_notes = .true.
                case ("-concise")
                    verbose=1
                case ("-v")
                    verbose=2
                case ("-vv")
                    verbose=3
                    silent_notes=.false.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       ! Manage defaults
       ! If not declared, hessfile and gradfile are the same as inpfile
       ! unless we are using nm file
       if (adjustl(nmfile) == "none") then
           if (adjustl(hessfile) == "same") then
               hessfile=inpfile
               if (adjustl(fth) == "guess")  fth=ft
           endif
           if (adjustl(gradfile) == "same") then
               gradfile=inpfile
               if (adjustl(ftg) == "guess")  ftg=ft
           endif
           ftn="-"
       else
           if (adjustl(hessfile) /= "same") &
            call alert_msg("note","Using nm file, disabling Hessian file")
           hessfile="none"
           fth="-"
           if (adjustl(gradfile) /= "same") &
            call alert_msg("note","Using nm file, disabling gradient file")
           gradfile="none"
           ftg="-"
       endif


       !Print options (to stderr)
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '             N M    C A R T E S I A N '    
        write(6,'(/,A)') '      Perform vibrational analysis based on  '
        write(6,'(A)')   '             Cartesian coordinates '        
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag           Description                     Value'
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,*)       '-f             Input file (structure&default)  ', trim(adjustl(inpfile))
        write(6,*)       '-ft            \_ FileType                     ', trim(adjustl(ft))
        write(6,*)       '-fhess         Hessian file                    ', trim(adjustl(hessfile))
        write(6,*)       '-fth           \_ FileType                     ', trim(adjustl(fth))
        write(6,*)       '-fgrad         Hessian file                    ', trim(adjustl(gradfile))
        write(6,*)       '-ftg           \_ FileType                     ', trim(adjustl(ftg))
        write(6,*)       '-fnm           Gradient file                   ', trim(adjustl(nmfile))
        write(6,*)       '-ftn           \_ FileType                     ', trim(adjustl(ftn))
        write(6,*)       '-[no]fulldiag  Diagonalize the 3Nx3N matrix   ',  full_diagonalize
        write(6,*)       ''
        write(6,*)       ' ** Options for animation **'
        write(6,*)       '-[no]animate   Generate animation files       ',  animate
        write(6,*)       '-nm            Selection of normal modes to    ', trim(adjustl(selection))
        write(6,*)       '               generate animations            '
        write(6,'(X,A,F5.2)') &
                         '-disp          Mode displacements for animate ',  Amplitude
        write(6,*)       '               (dimensionless displacements)'
        write(6,*)       '-[no]vmd       Launch VMD after computing the ',  call_vmd
        write(6,*)       '               modes (needs VMD installed)'
        write(6,'(X,A,I0)') &
                         '-movie         Number of cycles to record on   ',  movie_cycles
        write(6,*)       '               a movie with the animation'
        write(6,*)       ''
        write(6,*)       ' ** Options for correction non-stationary points **'
        write(6,*)       '-[no]vert      Correct with B derivatives for ',  vertical
        write(6,*)       '               non-stationary points'
        write(6,*)       '-cnx           Connectivity [filename|guess]   ', trim(adjustl(cnx_file))
        write(6,*)       '-intmode       Internal set:[zmat|sel|all]     ', trim(adjustl(def_internal))
        write(6,*)       '-intfile       File with ICs (for "sel")       ', trim(adjustl(intfile))
        write(6,*)       '-[no]sym       Use symmetry to form Zmat      ',  use_symmetry
!         write(6,*)       '-rmzfile       File deleting ICs from Zmat     ', trim(adjustl(rmzfile))
        write(6,*)       ''
        write(6,*)       '-h             Display this help              ',  need_help
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                       'Verbose level:  ', verbose        
        write(6,'(A)') '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input

    subroutine prepare_files(icoord,grofile,g09file,g96file,numfile,qfile,title,fulldiag)

        integer,intent(in) :: icoord
        character(len=*),intent(out) :: grofile,g09file,g96file,numfile,qfile,title
        logical,intent(in) :: fulldiag

        !Local
        character(len=150) :: dummy_char
        character(len=5)   :: full_label

        if (full_diagonalize) then
            full_label="full-"
        else
            full_label=""
        endif

        write(dummy_char,"(I0,X,A)") icoord
        title   = "Animation of normal mode "//trim(adjustl(dummy_char))
        g09file = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart.com"
        g96file = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart.g96"
        qfile   = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart_steps.dat"
        grofile = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart.gro"
        numfile = "Mode"//trim(adjustl(dummy_char))//"_"//trim(full_label)//"Cart_num.com"

        return

    end subroutine prepare_files


    subroutine displace_Xcoord(Lc,Nat,Qstep,X,Y,Z)

        real(8),dimension(:),intent(in)   :: Lc
        real(8),intent(in)                :: Qstep 
        integer,intent(in)                :: Nat
        real(8),dimension(:),intent(inout):: X,Y,Z

        !Local
        integer :: k, kk
        real(8) :: magic_number

!         magic_number=1.8895d0

        do k=1,Nat
            kk = 3*(k-1)
            X(k) = X(k) + Lc(kk+1) * Qstep
            Y(k) = Y(k) + Lc(kk+2) * Qstep
            Z(k) = Z(k) + Lc(kk+3) * Qstep
        enddo

        return
       
    end subroutine displace_Xcoord

end program normal_modes_cartesian

