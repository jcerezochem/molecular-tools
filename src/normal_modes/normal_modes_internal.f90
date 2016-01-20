program normal_modes_animation


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
    use xyz_manage_molec
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
    !  Internal thingies
    !============================================
    use internal_module
    use zmat_manage 
    use vibrational_analysis

    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: use_symmetry=.false.,   &
               include_hbonds=.false., &
               vertical=.false.,       &
               analytic_Bder=.false.,  &
               check_symmetry=.true.,  &
               animate=.true.,         &
               project_on_all=.true.,  &
               do_zmap
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule, molec_aux
    type(str_bonded) :: zmatgeom, allgeom, currentgeom
    integer,dimension(1:NDIM) :: isym
    integer,dimension(4,1:NDIM,1:NDIM) :: Osym
    integer :: Nsym
    integer :: Nat, Nvib, Ns
    integer :: Ns_zmat, Ns_all, Nvib_all
    character(len=5) :: PG
    !Job info
    character(len=20) :: calc, method, basis
    character(len=150):: title

    real(8) :: val1,val2,val3,val4,val5,val6
    integer :: imax,jmax,kmax, kk, jj1,jj2,jj3, kk1,kk2,kk3, kkmax, jjmax
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=50) :: dummy_char
    real(8) :: dist
    !io flags
    integer :: error, info
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, jj, iop
    !=============

    !====================== 
    ! PES topology and normal mode things
    real(8),dimension(:),allocatable :: Hlt
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, Hess_all, LL, LL_all
    real(8),dimension(NDIM) :: Freq, Freq_all, Factor, Grad, Grad_all
    !Moving normal modes
    character(len=50) :: selection="none"
    real(8) :: Amplitude = 2.d0, qcoord
    integer,dimension(1:NDIM) :: nm=0
    real(8) :: Qstep
    logical :: call_vmd = .false.
    character(len=10000) :: vmdcall
    integer :: Nsteps, Nsel, istep
    !MOVIE things
    logical :: movie_vmd = .false.
    integer :: movie_cycles=0,& !this means no movie
               movie_steps
    !====================== 

    !====================== 
    !INTERNAL CODE THINGS
    real(8),dimension(1:NDIM,1:NDIM) :: B, G, Asel, Asel_all
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    real(8),dimension(1:24,1:30,1:30) :: Bder1, Bder2, Bder3
    real(8),dimension(1:NDIM,1:NDIM) :: X,Xinv
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: S, Sref, Szmat, Sall
    integer,dimension(NDIM) :: S_sym
    ! Switches
    character(len=5) :: def_internal="ZMAT", def_internal_aux
    character(len=2) :: scan_type="NM"
    !Coordinate map
    integer,dimension(NDIM) :: Zmap, IntMap
    ! Number of ic (Shortcuts)
    integer :: nbonds, nangles, ndihed, nimprop
    !====================== 

    !====================== 
    !Auxiliar
    real(8),dimension(1:NDIM,1:NDIM) :: Aux, Aux2, Aux3
    real(8) :: Theta, Theta2
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
               O_PRJ=26,  &
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
    integer :: IOstatus, verbose_current
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
                     ! Movie
                     animate,movie_vmd, movie_cycles,                      &
                     ! Options (internal)
                     use_symmetry,def_internal,intfile,rmzfile,scan_type,  &
                     ! connectivity file
                     cnx_file,                                             &
                     ! (hidden)
                     analytic_Bder)
    call set_word_upper_case(def_internal)


    ! INTERNAL VIBRATIONAL ANALYSIS
 
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

    !Only read Grad/Hess or nm if we want to scan norma modes
    if (scan_type == "NM") then
        ! Vibrational analysis: either read from file (fcc) or from diagonalization of Hessian
        if (adjustl(nmfile) /= "none") then
            open(I_INP,file=nmfile,status='old',iostat=IOstatus)
            if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(nmfile)))
            call generic_nm_reader(I_INP,ftn,Nat,Nvib,Freq,LL)
            ! Show frequencies
            if (verbose>0) &
             call print_vector(6,Freq,Nvib,"Frequencies (cm-1)")
            ! The reader provide L in Normalized Cartesian. Need to Transform to Cartesian now
            call LcartNrm_to_Lmwc(Nat,Nvib,molecule%atom(:)%mass,LL,LL)
            call Lmwc_to_Lcart(Nat,Nvib,molecule%atom(:)%mass,LL,LL,error)
            close(I_INP)
        else
            ! HESSIAN FILE
            open(I_INP,file=hessfile,status='old',iostat=IOstatus)
            if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile)))
            allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
            call generic_Hessian_reader(I_INP,fth,Nat,Hlt,error)
            if (error /= 0) call alert_msg("fatal","Error reading Hessian (State1)")
            close(I_INP)
            ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
            call vibrations_Cart(Nat,molecule%atom(:)%X,molecule%atom(:)%Y,molecule%atom(:)%Z,&
                                 molecule%atom(:)%Mass,Hlt,Nvib,LL,Freq,error_flag=error)
            Nvib_all=Nvib
            k=0
            do i=1,3*Nat
            do j=1,i
                k=k+1
                Hess(i,j) = Hlt(k)
                Hess(j,i) = Hlt(k)
            enddo 
            enddo
            deallocate(Hlt)
           
            ! GRADIENT FILE
            if (vertical) then
                open(I_INP,file=gradfile,status='old',iostat=IOstatus)
                if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)))
                call generic_gradient_reader(I_INP,ftg,Nat,Grad,error)
                close(I_INP)
            endif
        endif
    else
        !We need to provide a value for Nvib. Lets assume non-liear molecules
        Nvib = 3*Nat-6
        Nvib_all=Nvib
    endif

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

    ! Define internal set
    Ns_zmat=0
    ! Shortcut logical: whether or not do zmap 
    !  1/ If -intmode zmat but we are removing coordinates
    do_zmap = def_internal=="ZMAT".and.rmzfile/="none"
    !  2/ If -intmode not equal zmat
    do_zmap = do_zmap.or.def_internal/="ZMAT"
    ! 3/ But only do it if we want the animation
    do_zmap = do_zmap.and.animate
    if (do_zmap) then 
        print*, "Preliminary Zmat analysis"
        ! Get Zmat first
        def_internal_aux="ZMAT"
        call define_internal_set(molecule,def_internal_aux,"none","none",use_symmetry,isym,S_sym,Ns)
        ! Get Ns at this point to latter form the proper Zmat
        Ns_zmat = Ns
        ! Get only the geom, and reuse molecule
        zmatgeom=molecule%geom
        ! Compute S values (to get initial values for frozen (rmzfile) ICs)
        verbose_current=verbose
!         verbose=0
        call set_geom_units(molecule,"Bohr")
        call internal_Wilson(molecule,Ns,Szmat,B)
        call set_geom_units(molecule,"Angs")
        verbose=verbose_current
        ! And reset bonded parameters
        call gen_bonded(molecule)
    endif
    if (project_on_all.and.def_internal/="ALL") then
        print'(/X,A,/)', "Preparing to make the projection on 'all' set"
        ! Prepare to compute selected frame on all coords
        ! Get all coordinates
        def_internal_aux="ALL"
        call define_internal_set(molecule,def_internal_aux,"none","none",use_symmetry,isym,S_sym,Ns)
        ! Get Ns at this point to latter form the proper Zmat
        Ns_all = Ns
        ! Get only the geom, and reuse molecule
        allgeom=molecule%geom
        ! And reset bonded parameters
        call gen_bonded(molecule)
        print'(X,A,/)', "Projection on 'all' prepared"
    elseif (def_internal=="ALL") then
        call alert_msg("note","The projection on the 'all' set is not done with -intmode all")
        project_on_all=.false.
    endif
    ! NOW, GET THE ACTUAL WORKING INTERNAL SET
    call define_internal_set(molecule,def_internal,intfile,rmzfile,use_symmetry,isym,S_sym,Ns)

    ! Set variables based on the working internal set 
    if (molecule%geom%nimprop/=0) then
        if (project_on_all) then
            call alert_msg("note","Cannot perform the projection with impropers")
            project_on_all=.false.
        endif
        if (do_zmap) then
            call alert_msg("note","Cannot produce animations with impropers")
            do_zmap=.false.
        endif
    endif
    if (do_zmap) then
        call internals_mapping(molecule%geom,zmatgeom,Zmap)
    endif

    !From now on, we'll use atomic units
    call set_geom_units(molecule,"Bohr")

    ! SCAN JOBS
    ! Get the selection of normal modes to represent
    if (adjustl(selection) == "all") then
        Nsel=Nvib
        nm(1:Nvib) = (/(i, i=1,Nvib)/)
    else if (adjustl(selection) == "none") then
        Nsel=0
    else
        call selection2intlist(selection,nm,Nsel)
    endif
    !Two possible jobs:
    if (scan_type=="IN") then
        !--------------------------------------
        ! 1. Internal Coordinates SCAN
        !--------------------------------------
        ! We need to call Wilson to get ModeDef
        call internal_Wilson(molecule,Nvib,S,B,ModeDef)
        LL(1:Ns,1:Ns)=0.d0
        do i=1,Nsel
            LL(nm(i),nm(i)) = 1.d0
            Freq(nm(i))     = 1.d0
            Factor(nm(i))   = 1.d0
        enddo
     else
        !--------------------------------------
        ! 2. Normal Mode SCAN
        !--------------------------------------
        ! If requested, get the vibrations with the 'all' set
        if (project_on_all) then
            print'(/,X,A)', "Vibrational anlysis with 'all' set"
            ! First get the geom for current state
            currentgeom=molecule%geom
            
            ! Then make the whole analysis for 'All' set
            molecule%geom=allgeom
            Hess_all(1:3*Nat,1:3*Nat) = Hess(1:3*Nat,1:3*Nat)
            call internal_Wilson(molecule,Ns_all,Sall,B,ModeDef)
            !SOLVE GF METHOD TO GET NM AND FREQ
            call internal_Gmetric(Nat,Ns_all,molecule%atom(:)%mass,B,G)

            if (vertical) then
                call calc_Bder(molecule,Ns_all,Bder,analytic_Bder)
            endif

            call redundant2nonredundant(Ns_all,Nvib_all,G,Asel_all)
            ! Rotate Bmatrix
            B(1:Nvib_all,1:3*Nat) = matrix_product(Nvib_all,3*Nat,Ns_all,Asel_all,B,tA=.true.)
            ! Rotate Gmatrix
            G(1:Nvib_all,1:Nvib_all) = matrix_basisrot(Nvib_all,Ns_all,Asel_all(1:Ns_all,1:Nvib_all),G,counter=.true.)
            ! Rotate Bders
            if (vertical) then
                do j=1,3*Nat
                    Bder(1:Nvib_all,j,1:3*Nat) =  &
                     matrix_product(Nvib_all,3*Nat,Ns_all,Asel_all,Bder(1:Ns_all,j,1:3*Nat),tA=.true.)
                enddo
            endif

            if (vertical) then
                ! We save the original gradient
                Grad_all(1:3*Nat) = Grad(1:3*Nat)
                call HessianCart2int(Nat,Nvib_all,Hess_all,molecule%atom(:)%mass,B,G,Grad_all,Bder)
                if (verbose>2) then
                    do i=1,Nvib_all
                        write(tmpfile,'(A,I0,A)') "Bder, ic=",i
                        call MAT0(6,Bder(i,:,:),3*Nat,3*Nat,trim(tmpfile))
                    enddo
                endif

                if (check_symmetry) then
                    print'(/,X,A)', "---------------------------------------"
                    print'(X,A  )', " Check effect of symmetry operations"
                    print'(X,A  )', " on the correction term gs^t\beta"
                    print'(X,A  )', "---------------------------------------"
                    molecule%PG="XX"
                    call symm_atoms(molecule,isym,Osym,rotate=.false.,nsym_ops=nsym)
                    ! Check the symmetry of the correction term
                    ! First compute the correction term
                    do i=1,3*Nat
                    do j=1,3*Nat
                        Aux2(i,j) = 0.d0
                        do k=1,Nvib_all
                            Aux2(i,j) = Aux2(i,j) + Bder(k,i,j)*Grad_all(k)
                        enddo
                    enddo
                    enddo
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

            else
                call HessianCart2int(Nat,Nvib_all,Hess_all,molecule%atom(:)%mass,B,G)
            endif
            call gf_method(Nvib_all,G,Hess_all,LL_all,Freq_all,X,Xinv)

           ! Restore original geom
           molecule%geom=currentgeom
            
            print'(X,A,/)', "Vibrational anlysis done. Ready for projection"

        endif

        call internal_Wilson(molecule,Ns,S,B,ModeDef)
        if (nmfile == "none") then
            !SOLVE GF METHOD TO GET NM AND FREQ
            call internal_Gmetric(Nat,Ns,molecule%atom(:)%mass,B,G)

            if (vertical) then
                call calc_Bder(molecule,Ns,Bder,analytic_Bder)
            endif

            ! The diagonalization of the G matrix can be donne with all sets
            ! (either redundant or non-redundant), and it is the best way to 
            ! set the number of vibrations. The following call also set Nvib
            call redundant2nonredundant(Ns,Nvib,G,Asel)
            ! Rotate Bmatrix
            B(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel,B,tA=.true.)
            ! Rotate Gmatrix
            G(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel(1:Ns,1:Nvib),G,counter=.true.)
            ! Rotate Bders
            if (vertical) then
                do j=1,3*Nat
                    Bder(1:Nvib,j,1:3*Nat) =  matrix_product(Nvib,3*Nat,Ns,Asel,Bder(1:Ns,j,1:3*Nat),tA=.true.)
                enddo
            endif

            if (vertical) then
                call HessianCart2int(Nat,Nvib,Hess,molecule%atom(:)%mass,B,G,Grad,Bder)
                if (verbose>2) then
                    do i=1,Nvib
                        write(tmpfile,'(A,I0,A)') "Bder, ic=",i
                        call MAT0(6,Bder(i,:,:),3*Nat,3*Nat,trim(tmpfile))
                    enddo
                endif

                if (check_symmetry) then
                    print'(/,X,A)', "---------------------------------------"
                    print'(X,A  )', " Check effect of symmetry operations"
                    print'(X,A  )', " on the correction term gs^t\beta"
                    print'(X,A  )', "---------------------------------------"
                    molecule%PG="XX"
                    call symm_atoms(molecule,isym,Osym,rotate=.false.,nsym_ops=nsym)
                    ! Check the symmetry of the correction term
                    ! First compute the correction term
                    do i=1,3*Nat
                    do j=1,3*Nat
                        Aux2(i,j) = 0.d0
                        do k=1,Nvib
                            Aux2(i,j) = Aux2(i,j) + Bder(k,i,j)*Grad(k)
                        enddo
                    enddo
                    enddo
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

            else
                call HessianCart2int(Nat,Nvib,Hess,molecule%atom(:)%mass,B,G)
            endif
            call gf_method(Nvib,G,Hess,LL,Freq,X,Xinv)
            if (verbose>1) then
                ! Analyze normal modes in terms of the redundant set
                if (Nvib<Ns) then
                    Aux(1:Ns,1:Nvib) = matrix_product(Ns,Nvib,Nvib,Asel,LL)
                else
                    Aux(1:Ns,1:Nvib) = LL(1:Ns,1:Nvib)
                endif
                if (use_symmetry) then
                    call analyze_internal(Nvib,Ns,Aux,Freq,ModeDef,S_sym)
                else
                    call analyze_internal(Nvib,Ns,Aux,Freq,ModeDef)
                endif
            endif
        else
            ! Transform LcartNrm read from file to Ls (this is approximate if -vert!)
            call Lcart_to_Ls(Nat,Nvib,B,LL,LL,error)
        endif
        !Define the Factor to convert shift into addimensional displacements
        ! from the shift in SI units:
        Factor(1:Nvib) = dsqrt(dabs(Freq(1:Nvib))*1.d2*clight*2.d0*PI/plankbar)
        ! but we need it from au not SI
        Factor(1:Nvib)=Factor(1:Nvib)*BOHRtoM*dsqrt(AUtoKG)
    endif

    if (project_on_all) then
        print*, "Projecting reduced normal mode set on whole set"
        ! Map 'all' modes with reduced space
        call internals_mapping(currentgeom,allgeom,IntMap)
   
        do i=1,Ns_all
            j=IntMap(i) 
            if (j/=0) then
                Aux2(i,1:Nvib) = Asel(j,1:Nvib)
            else
                Aux2(i,1:Nvib) = 0.d0
            endif
        enddo
        !
        ! Prj1 = (LL)^-1 * (Asel)^t * Asel_all * LL_all
        !
        ! (Asel)^t * Asel_all
        Aux(1:Nvib,1:Nvib_all) = matrix_product(Nvib,Nvib_all,Ns_all,Aux2,Asel_all,tA=.true.)
        ! (LL)^-1 * [(Asel)^t * Asel_all]
        Aux3(1:Nvib,1:Nvib) = inverse_realgen(Nvib,LL)
        Aux(1:Nvib,1:Nvib_all) = matrix_product(Nvib,Nvib_all,Nvib,Aux3,Aux)
        ! [(LL)^-1 * (Asel)^t * Asel_all] * LL_all
        Aux(1:Nvib,1:Nvib_all) = matrix_product(Nvib,Nvib_all,Nvib_all,Aux,LL_all)

        print*, "Projection an modes with 'all' internal set: "
        print'(X,A,I0,X,I0,/)', " Prj = (LL_current)^-1 * (LL_all). Size: ", Nvib, Nvib_all
        open(O_PRJ,file="Prj1_modes.dat")
        do i=1,Nvib 
        do j=1,Nvib_all
            write(O_PRJ,*) Aux(i,j)
        enddo
        enddo
        close(O_PRJ)

        !Make the assignememts
        print'(/,X,A)', "--------------------------------------------------------"
        print'(X,A)',   " Assignments of reduced set based on Prj"
        print'(X,A)',   "--------------------------------------------------------"
        print'(X,A)',   " Reduced/Freq     Complete/Freq      Coef^2     Norm"
        do i=1,Nvib
            Theta=0.d0
            Theta2=0.d0
            do j=1,Nvib_all
                if (Theta<=abs(Aux(i,j))) then
                    Theta=abs(Aux(i,j))
                    k=j
                endif 
                Theta2=Theta2+Aux(i,j)**2
            enddo
            print'(I4,2X,F8.2,4X,I4,2X,F8.2,3X,F8.3,3X,F8.3)', &
                          i, Freq(i), k, Freq_all(k), Theta**2, Theta2
        enddo
        print'(X,A,/)',   "--------------------------------------------------------"
        
    endif

    ! Check if we should continue
    if (molecule%geom%nimprop/=0) then
        if (Nsel/=0) call alert_msg("warning","Animations not yet possible "//&
                                              "when impropers are selected")
        Nsel=0
    endif
    if (def_internal=="SEL") then
        if (Nsel/=0) call alert_msg("warning","Animations may missbehave "//&
                                              "with -intmode sel")
    endif
    if (.not.animate) then
        Nsel=0
    endif

    ! If redundant set, transform from orthogonal non-redundant
!     if (Nvib < Ns) then
        print'(/,X,A,/)', "Transform from non-redundant orthogonal to original redundant set"
        LL(1:Ns,1:Nvib) = matrix_product(Ns,Nvib,Nvib,Asel,LL)
!     endif


    !==========================================================0
    !  Normal mode displacements
    !==========================================================0
    if (Nsel > 0) then
        ! Tune Ns to activate switches
        if (rmzfile/="none") Ns=0
        ! Take number of ICs as Shortcuts
        nbonds = molecule%geom%nbonds
        nangles= molecule%geom%nangles
        ndihed = molecule%geom%ndihed
        ! Initialization
        Sref = S
        ! To ensure that we always have the same orientation, we stablish the reference here
        ! this can be used to use the input structure as reference (this might need also a 
        ! traslation if not at COM -> not necesary, the L matrices are not dependent on the 
        ! COM position, only on the orientation)
        if (Ns_zmat /= 0 .and. scan_type == "NM") then
            ! From now on, we use the zmatgeom 
            molecule%geom = zmatgeom
            S(1:Ns_zmat) = map_Zmatrix(Ns_zmat,S,Zmap,Szmat)
        endif
        call zmat2cart(molecule,S)
        ! Save state as reference frame for RMSD fit (in AA)
        molec_aux=molecule
        ! Default steps (to be set by the user..)
        Nsteps = 101
        if ( mod(Nsteps,2) == 0 ) Nsteps = Nsteps + 1 ! ensure odd number of steps (so we have same left and right)
        ! Qstep is dimless
        Qstep = Amplitude/float(Nsteps-1)*2.d0  ! Do the range (-A ... +A)
        molecule%atom(1:molecule%natoms)%resname = "RES" ! For printing
    endif
    
    !---------------------------------------------
    ! Run over all selected modes/internals
    !---------------------------------------------
    do jj=1,Nsel 
        k=0 ! equilibrium corresponds to k=0
        j = nm(jj)
        if (scan_type == "NM") then
            if (verbose>0) &
             print'(X,A,I0,A)', "Generating Mode ", j, "..."
        else 
            if (verbose>0) &
             print'(X,A,I0,A)', "Generating Scan for IC ", j, ": "//trim(adjustl(ModeDef(j)))//"..."
        endif

        ! Set initial values for the scanned coordinate
        if (scan_type =="IN") then
            qcoord = Sref(j)
        else
            qcoord = 0.d0
        endif 

        ! Prepare and open files
        call prepare_files(j,ModeDef(j),scan_type,&
                           grofile,g09file,g96file,numfile,qfile,title)
        open(O_GRO,file=grofile)
        open(O_G09,file=g09file)
        open(O_G96,file=g96file)
        open(O_Q  ,file=qfile)
        open(O_NUM,file=numfile)

        !===========================
        !Start from equilibrium. 
        !===========================
        if (verbose>1) &
         print'(/,A,I0)', "STEP:", k
        ! Update
        write(molecule%title,'(A,I0,A,2(X,F12.6))') &
         trim(adjustl((title)))//" Step ",k," Disp = ", qcoord, qcoord*Factor(j)
!         call zmat2cart(molecule,S)
        !call rmsd_fit_frame(state,ref): efficient but not always works. If so, it uses rmsd_fit_frame_brute(state,ref)
        call rmsd_fit_frame(molecule,molec_aux,info)
        if (info /= 0) then
            print'(X,A,I0)', "RMSD fit failed at Step: ", k
            call rmsd_fit_frame_brute(molecule,molec_aux,dist)
        endif
        !Transform to AA and export coords and put back into BOHR
        call set_geom_units(molecule,"Angs")
        call write_gro(O_GRO,molecule)
        ! Save state as reference frame for RMSD fit (in AA)
        molec_aux=molecule
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
            ! Displace (displace always from reference to avoid error propagation)
            i=istep
            S=Sref
            call displace_Scoord(LL(:,j),nbonds,nangles,ndihed,Qstep/Factor(j)*i,S)
            ! Get Cart coordinates
            if (Ns /= Ns_zmat .and. scan_type == "NM") then
                S(1:Ns_zmat) = map_Zmatrix(Ns_zmat,S,Zmap,Szmat)
            endif
            call zmat2cart(molecule,S)
            !call rmsd_fit_frame(state,ref): efficient but not always works. If so, it uses rmsd_fit_frame_brute(state,ref)
            call rmsd_fit_frame(molecule,molec_aux,info)
            if (info /= 0) then
                print'(X,A,I0)', "RMSD fit failed at Step: ", k
                call rmsd_fit_frame_brute(molecule,molec_aux,dist)
            endif
            ! PRINT
            !Transform to AA and comparae with last step (stored in state)  -- this should be detected and fix by the subroutines
            call set_geom_units(molecule,"Angs")
            ! Write GRO from the beginign and G96/G09 only when reach max amplitude
            call write_gro(O_GRO,molecule)
            if (k==(nsteps-1)/2) then
                call write_gcom(O_G09,molecule,g09file,calc,method,basis,molecule%title)
                call write_g96(O_G96,molecule)
                write(O_Q,*) qcoord, qcoord*Factor(j)
            endif
            ! Save state as reference frame for RMSD fit (in AA)
            molec_aux=molecule
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
            ! Displace (displace always from reference to avoid error propagation)
            i=(nsteps-1)/2-istep
            S=Sref
            call displace_Scoord(LL(:,j),nbonds,nangles,ndihed,Qstep/Factor(j)*i,S)
            ! Get Cart coordinates
            if (Ns /= Ns_zmat .and. scan_type == "NM") then
                S(1:Ns_zmat) = map_Zmatrix(Ns_zmat,S,Zmap,Szmat)
            endif
            call zmat2cart(molecule,S)
            !Transform to AA and comparae with last step (stored in state) -- comparison in AA
            call set_geom_units(molecule,"Angs")
            !call rmsd_fit_frame(state,ref): efficient but not always works. If so, it uses rmsd_fit_frame_brute(state,ref)
            call rmsd_fit_frame(molecule,molec_aux,info)
            if (info /= 0) then
                print'(X,A,I0)', "RMSD fit failed at Step: ", k
                call rmsd_fit_frame_brute(molecule,molec_aux,dist)
            endif
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
            ! Save state as reference frame for RMSD fit (in AA)
            molec_aux=molecule
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
            ! Displace (displace always from reference to avoid error propagation)
            i=-(nsteps-1)/2+istep
            S=Sref
            call displace_Scoord(LL(:,j),nbonds,nangles,ndihed,Qstep/Factor(j)*i,S)
            ! Get Cart coordinates
            if (Ns /= Ns_zmat .and. scan_type == "NM") then
                S(1:Ns_zmat) = map_Zmatrix(Ns_zmat,S,Zmap,Szmat)
            endif
            call zmat2cart(molecule,S)
            !call rmsd_fit_frame(state,ref): efficient but not always works. If so, it uses rmsd_fit_frame_brute(state,ref)
            call rmsd_fit_frame(molecule,molec_aux,info)
            if (info /= 0) then
                print'(X,A,I0)', "RMSD fit failed at Step: ", k
                call rmsd_fit_frame_brute(molecule,molec_aux,dist)
            endif
            ! PRINT
            !Transform to AA and comparae with last step (stored in state)  -- this should be detected and fix by the subroutines
            call set_geom_units(molecule,"Angs")
            ! Write only GRO 
            call write_gro(O_GRO,molecule)
            ! Save state as reference frame for RMSD fit (in AA)
            molec_aux=molecule
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
            call prepare_files(j,ModeDef(j),scan_type,&
                              grofile,g09file,g96file,numfile,qfile,title)
            vmdcall = trim(adjustl(vmdcall))//" "//trim(adjustl(grofile))
        enddo
        vmdcall = trim(adjustl(vmdcall))//" -e vmd_conf.dat"
        open(O_MOV,file="vmd_call.cmd")
        write(O_MOV,*) trim(adjustl(vmdcall))
        close(O_MOV)
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
            call prepare_files(j,ModeDef(j),scan_type,&
                              grofile,g09file,g96file,numfile,qfile,title)
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
                           ! Movie
                           animate,movie_vmd, movie_cycles,                      &
                           ! Options (internal)
                           use_symmetry,def_internal,intfile,rmzfile,scan_type,  &
                           ! connectivity file
                           cnx_file,                                             &
                           ! (hidden)
                           analytic_Bder)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn, &
                                          intfile,rmzfile,scan_type,def_internal,selection, &
                                          cnx_file
        real(8),intent(inout)          :: Amplitude
        logical,intent(inout)          :: call_vmd, include_hbonds,vertical, use_symmetry,movie_vmd,animate,&
                                          analytic_Bder
        integer,intent(inout)          :: movie_cycles

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        character(len=200) :: int_selection, nm_selection

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

                case ("-sym")
                    use_symmetry=.true.
                case ("-nosym")
                    use_symmetry=.false.

                case ("-vert")
                    vertical=.true.
                case ("-novert")
                    vertical=.false.

                case ("-nm") 
                    scan_type="NM"
                    call getarg(i+1, selection)
                    argument_retrieved=.true.

                case ("-int")
                    scan_type="IN"
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
        
                case ("-h")
                    need_help=.true.

                ! HIDDEN FLAGS

                case ("-anaBder")
                    analytic_Bder=.true.
                case ("-noanaBder")
                    analytic_Bder=.false.

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

       ! Select internal or normal modes
       if (scan_type == "NM") then
           int_selection="-"
           nm_selection =selection
       elseif (scan_type == "IN") then
           nm_selection ="-"
           int_selection=selection
       endif


       !Print options (to stderr)
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '             N M    I N T E R N A L '    
        write(6,'(/,A)') '      Perform vibrational analysis based on  '
        write(6,'(A)')   '             internal coordinates '        
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
        write(6,*)       '-cnx           Connectivity [filename|guess]   ', trim(adjustl(cnx_file))
!         write(6,*)       '-fnm           Gradient file                   ', trim(adjustl(nmfile))
!         write(6,*)       '-ftn           \_ FileType                     ', trim(adjustl(ftn))
        write(6,*)       '-intmode       Internal set:[zmat|sel|all]     ', trim(adjustl(def_internal))
        write(6,*)       '-intfile       File with ICs (for "sel")       ', trim(adjustl(intfile))
        write(6,*)       '-rmzfile       File deleting ICs from Zmat     ', trim(adjustl(rmzfile))
        write(6,*)       '-[no]sym       Use symmetry to form Zmat      ',  use_symmetry
        write(6,*)       '-[no]vert      Correct with B derivatives for ',  vertical
        write(6,*)       '               non-stationary points'
        write(6,*)       ''
        write(6,*)       ' ** Options for animation **'
        write(6,*)       '-[no]animate   Generate animation files       ',  animate
        write(6,*)       '-nm            Selection of normal modes to    ', trim(adjustl(nm_selection))
        write(6,*)       '               generate animations             '
        write(6,*)       '-int           Selection of internal coords    ', trim(adjustl(int_selection))
        write(6,*)       '               to generate animations          '
        write(6,'(X,A,F5.2)') &
                         '-disp          Mode displacements for animate ',  Amplitude
        write(6,*)       '               (dimensionless displacements)'
        write(6,*)       '-[no]vmd       Launch VMD after computing the ',  call_vmd
        write(6,*)       '               modes (needs VMD installed)'
        write(6,'(X,A,I0)') &
                         '-movie         Number of cycles to record on   ',  movie_cycles
        write(6,*)       '               a movie with the animation'
        write(6,*)       ''
        write(6,*)       '-h             Display this help              ',  need_help
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                       'Verbose level:  ', verbose        
        write(6,'(A)') '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input

    subroutine prepare_files(icoord,label,scan_type,grofile,g09file,g96file,numfile,qfile,title)

        integer,intent(in) :: icoord
        character(len=*),intent(out) :: label, scan_type
        character(len=*),intent(out) :: grofile,g09file,g96file,numfile,qfile,title

        !Local
        character(len=150) :: dummy_char

        if (scan_type=="IN") then
            write(dummy_char,"(I0,X,A)") icoord
            title   = "Animation of IC "//trim(adjustl(dummy_char))//"("//trim(adjustl(label))//")"
            g09file = "Coord"//trim(adjustl(dummy_char))//"_int.com"
            g96file = "Coord"//trim(adjustl(dummy_char))//"_int.g96"
            qfile   = "Coord"//trim(adjustl(dummy_char))//"_int_steps.dat"
            grofile = "Coord"//trim(adjustl(dummy_char))//"_int.gro" 
            numfile = "Coord"//trim(adjustl(dummy_char))//"_int_num.com"
        else
            write(dummy_char,"(I0,X,A)") icoord
            title   = "Animation of normal mode "//trim(adjustl(dummy_char))
            g09file = "Mode"//trim(adjustl(dummy_char))//"_int.com"
            g96file = "Mode"//trim(adjustl(dummy_char))//"_int.g96"
            qfile   = "Mode"//trim(adjustl(dummy_char))//"_int_steps.dat"
            grofile = "Mode"//trim(adjustl(dummy_char))//"_int.gro"
            numfile = "Mode"//trim(adjustl(dummy_char))//"_int_num.com"
        endif

        return

    end subroutine prepare_files


    subroutine displace_Scoord(Lc,nbonds,nangles,ndihed,Qstep,S)

        real(8),dimension(:),intent(in)   :: Lc
        real(8),intent(in)                :: Qstep 
        integer,intent(in)                :: nbonds,nangles,ndihed
        real(8),dimension(:),intent(inout):: S

        !Local
        integer :: i, k

        k=0
        if (verbose>1) &
         print*, "Bonds"
        do i=1,nbonds
            k=k+1
            S(k) = S(k) + Lc(k) * Qstep
        enddo
        if (verbose>1) &
         print*, "Angles"
        do i=1,nangles
            k=k+1
            S(k) = S(k) + Lc(k) * Qstep
        enddo
        if (verbose>1) &
         print*, "Dihedrals"
        do i=1,ndihed
            k=k+1
            S(k) = S(k) + Lc(k) * Qstep
            if (S(k) >  PI) S(k)=S(k)-2.d0*PI
            if (S(k) < -PI) S(k)=S(k)+2.d0*PI
        enddo

        return
       
    end subroutine displace_Scoord

end program normal_modes_animation

