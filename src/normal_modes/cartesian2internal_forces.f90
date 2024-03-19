program normal_modes_internal


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
    use io
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
    use atomic_geom
    use symmetry
    !============================================
    !  Internal thingies
    !============================================
    use internal_module
    use zmat_manage 
    use vibrational_analysis
    use thermochemistry
    use vertical_model

    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: use_symmetry=.false.,   &
               include_hbonds=.false., &
               vertical=.false.,       &
               analytic_Bder=.true.,  &
               check_symmetry=.true.,  &
               animate=.true.,         &
               project_on_all=.false.,  &
               apply_projection_matrix=.false., &
               complementay_projection=.false., &
               rm_gradcoord=.false.,            &
               complementay_gradient = .false., &
               do_zmap
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule, molec_aux
    type(str_bonded) :: zmatgeom, allgeom, currentgeom, iputgeom
    integer,dimension(1:NDIM) :: isym
    integer,dimension(4,1:NDIM,1:NDIM) :: Osym
    integer :: Nsym
    integer :: Nat, Nvib, Ns, Nvib0, NvibP, NvibP2, Nf
    integer :: Ns_zmat, Ns_all, Nvib_all
    character(len=5) :: PG
    real(8) :: Tthermo=0.d0
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
    integer :: i,j,k,l, ii,jj, iop
    !=============

    !====================== 
    ! PES topology and normal mode things
    real(8),dimension(:),allocatable :: Hlt, grdx
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, Hess_all, LL, LL_all, gBder, P, Fltr
    real(8),dimension(NDIM) :: Freq, Freq_all, Factor, Grad, Vec1, Vec
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
    real(8),dimension(1:NDIM,1:NDIM) :: B, G, Asel, Asel_all, B0, G0, Asel0, Bprj, Bprj2
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    real(8),dimension(1:NDIM,1:NDIM) :: X,Xinv
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    character(len=400)                 :: CombDef
    !VECTORS
    real(8),dimension(NDIM) :: S, Sref, Szmat, Sall, DeltaS
    integer,dimension(NDIM) :: S_sym
    ! Switches
    character(len=4) :: def_internal="ALL",  & ! To do the vibrational analysis
                        def_internal0='defa',& ! defa(ult) is "the same as working set"
                        conversion_i2c="ZMAT"
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
    character(len=5) :: current_symm
    !====================== 

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_SYM=12,  &
               I_RMF=16,  &
               I_CNX=17,  &
               I_MAS=19,  &
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
                         intfile0 ="default", & ! default is "the same as working set"
                         rmzfile  ="none", &
                         symm_file="none", &
                         mass_file="none", &
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

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(molecule)
    call allocate_atoms(molec_aux)
    !===========================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(&
                     ! input data
                     inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn,mass_file,&
                     ! Options (general)
                     Amplitude,call_vmd,include_hbonds,selection,vertical, &
                     ! Movie
                     animate,movie_vmd, movie_cycles,conversion_i2c,       &
                     ! Options (internal)
                     use_symmetry,def_internal,def_internal0,intfile,intfile0,&
                     apply_projection_matrix,complementay_projection,      &
                     rmzfile,scan_type,  &
                     project_on_all,rm_gradcoord,complementay_gradient,    &
                     ! connectivity file
                     cnx_file,                                             &
                     ! thermochemical analysis
                     Tthermo,                                              &
                     ! (hidden)
                     analytic_Bder)
    call set_word_upper_case(def_internal)
    call set_word_upper_case(conversion_i2c)


    ! INTERNAL VIBRATIONAL ANALYSIS
 
    ! 1. READ DATA
    ! ---------------------------------
    call heading(6,"INPUT DATA")
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
    call statement(6,"READING STRUCTURE...")
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,molecule,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry")
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

    ! Read mass from file if given
    if (adjustl(mass_file) /= "none") then
        print'(/,X,A)', "Reading atomic masses from: "//trim(adjustl(mass_file))
        open(I_MAS,file=mass_file,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(mass_file)) )
        do i=1,Nat
            read(I_MAS,*,iostat=IOstatus) molecule%atom(i)%mass 
            if (IOstatus /= 0) call alert_msg( "fatal","While reading "//trim(adjustl(mass_file)) )
        enddo
        close(I_MAS)
    endif


           
    ! GRADIENT FILE
    if (adjustl(gradfile) /= "none") then
        call statement(6,"READING GRADIENT FILE...")
        open(I_INP,file=gradfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile)) )
        allocate(grdx(1:3*Nat))
        call generic_gradient_reader(I_INP,ftg,Nat,grdx,error)
        close(I_INP)
        if (error /= 0) then
            print*, "Error reading the Gradient. It will be set to zero"
            grdx(1:3*Nat) = 0.d0
        endif
    else
        grdx(1:3*Nat) = 0.d0
    endif


    !We need to provide a value for Nvib. Lets assume non-liear molecules
    Nvib = 3*Nat-6
    Nvib_all=Nvib

    
    !----------------------------------
    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! General Actions:
    !******************
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

    call set_geom_units(molecule,"BOHR")


    Grad(1:3*Nat) = grdx
    
    
    !Generate bonded info
    call gen_bonded(molecule)

    !---------------------------------------
    call subheading(6,"Getting the internal set")
    ! NOW, GET THE ACTUAL WORKING INTERNAL SET
    call define_internal_set(molecule,"ALL",intfile,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Fltr)
    allgeom = molecule%geom
    call statement(6,"And now get the actual working set")
    call define_internal_set(molecule,def_internal0,intfile0,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Fltr)
    !---------------------------------------    


    ! Compute B, G and, if needed, Bder
    call statement(6,"Computing B, G (with ALL set)",keep_case=.true.)
    call internal_Wilson(molecule,Ns,S,B,ModeDef)
    call internal_Gmetric(Nat,Ns,molecule%atom(:)%mass,B,G)

    ! Get non-redundant space
    call subsubheading(6,"Getting the actual vibrational space dimension")
    call redundant2nonredundant(Ns,Nvib,G,Asel)
    call statement(6,"Rotate B, G to non-redundant space", keep_case=.true.)
    ! Rotate Bmatrix
    B(1:Nvib,1:3*Nat) = matrix_product(Nvib,3*Nat,Ns,Asel,B,tA=.true.)
    ! Rotate Gmatrix
    G(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Ns,Asel(1:Ns,1:Nvib),G,counter=.true.)

    call Gradcart2int(Nat,Nvib,Grad,molecule%atom(:)%mass,B,G)
    
    ! Transform to redundant set 
    do i=1,Ns
        Vec1(i) = 0.0
        do k=1,Nvib
            Vec1(i) = Vec1(i) + Asel(i,k)*Grad(k)
        enddo
    enddo
    
    print*, "Internal Gradient (kJ/mol/[XX]), where [XX]=nm(bonds) or rad(ange/dihed)"
    k=0
    do i=1,molecule%geom%nbonds
        k=k+1
        print*, k, ModeDef(k), Vec1(k)*HtoKCALM*CALtoJ/BOHRtoNM
    enddo
    do i=1,molecule%geom%nangles
        k=k+1
        print*, k, ModeDef(k), Vec1(k)*HtoKCALM*CALtoJ
    enddo
    do i=1,molecule%geom%ndihed
        k=k+1
        print*, k, ModeDef(k), Vec1(k)*HtoKCALM*CALtoJ
    enddo
    
    print*, ""

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(&
                           ! input data
                           inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn,mass_file,&
                           ! Options (general)
                           Amplitude,call_vmd,include_hbonds,selection,vertical, &
                           ! Movie
                           animate,movie_vmd, movie_cycles,conversion_i2c,        &
                           ! Options (internal)
                           use_symmetry,def_internal,def_internal0,intfile,intfile0,&
                           apply_projection_matrix,complementay_projection,  &
                           rmzfile,scan_type,  &
                           project_on_all, rm_gradcoord,complementay_gradient,   &
                           ! connectivity file
                           cnx_file,                                             &
                           ! thermochemical analysis
                           Tthermo,                                              &
                           ! (hidden)
                           analytic_Bder)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,nmfile,ftn, &
                                          intfile,intfile0,rmzfile,scan_type,def_internal, &
                                          selection,cnx_file,def_internal0,conversion_i2c,mass_file
        real(8),intent(inout)          :: Amplitude,Tthermo
        logical,intent(inout)          :: call_vmd, include_hbonds,vertical, use_symmetry,movie_vmd,animate,&
                                          analytic_Bder,project_on_all,apply_projection_matrix,complementay_projection,&
                                          rm_gradcoord,complementay_gradient
        integer,intent(inout)          :: movie_cycles

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        character(len=200) :: int_selection, nm_selection
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

                case ("-fmass") 
                    call getarg(i+1, mass_file)
                    argument_retrieved=.true.

                case ("-cnx") 
                    call getarg(i+1, cnx_file)
                    argument_retrieved=.true.

                case ("-rmgrad")
                    rm_gradcoord=.true.
                case ("-normgrad")
                    rm_gradcoord=.false.
                case ("-rmgrad-c")
                    complementay_gradient=.true.
                    rm_gradcoord=.true.
                case ("-normgrad-c")
                    complementay_gradient=.false.

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-intfile2") 
                    call getarg(i+1, intfile0)
                    argument_retrieved=.true.

                case ("-rmzfile") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.

                case ("-intmode")
                    call getarg(i+1, def_internal)
                    argument_retrieved=.true.

                case ("-alg")
                    call getarg(i+1, conversion_i2c)
                    argument_retrieved=.true.

                case ("-intmode0")
                    call getarg(i+1, def_internal0)
                    argument_retrieved=.true.

                case ("-prjS")
                    apply_projection_matrix=.true.
                    complementay_projection=.false.
                case ("-noprjS")
                    apply_projection_matrix=.false.

                case ("-prjS-c")
                    apply_projection_matrix=.true.
                    complementay_projection=.true.
                case ("-noprjS-c")
                    apply_projection_matrix=.false.

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

                case ("-prjall")
                    project_on_all=.true.
                case ("-noprjall")
                    project_on_all=.false.


                case ("-movie")
                    call getarg(i+1, arg)
                    read(arg,*) movie_cycles
                    movie_vmd=.true.
                    argument_retrieved=.true.

                case ("-include_hb")
                    include_hbonds=.true.

                case ("-thermo")
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    read(arg,*) Tthermo
        
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

       ! Take defaults for the internal set for correction only
       if (def_internal0 == "defa") def_internal0=def_internal
       if (adjustl(intfile0)=="default") intfile0=intfile


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
        write(6,*)       '-fmass         Mass file (optional)            ', trim(adjustl(mass_file)) 
!         write(6,*)       '-fnm           Gradient file                   ', trim(adjustl(nmfile))
!         write(6,*)       '-ftn           \_ FileType                     ', trim(adjustl(ftn))
        write(6,*)       '-[no]prjS      Apply projection matrix to     ', apply_projection_matrix
        write(6,*)       '               rotate Grad and Hess.'
        write(6,*)       '               Projection P=B^+B, where the'
        write(6,*)       '-[no]prjS-c    Use the complentary projection ', complementay_projection
        write(6,*)       '               P=I - B^+B'
        write(6,*)       '-[no]rmgrad    Remove coordinate along the    ', rm_gradcoord
        write(6,*)       '               grandient                      '
        write(6,*)       '-[no]rmgrad-c  Use complementary projection   ', rm_gradcoord
        write(6,*)       '               (i.e., use project on grad)    '
        write(6,*)       '-intmode0      Internal set:[zmat|sel|all]    ', trim(adjustl(def_internal0))
        write(6,*)       '               to compute additional terms    '
        write(6,*)       '               (vertical method)              '
        write(6,*)       '-intmode       Internal set:[zmat|sel|all]     ', trim(adjustl(def_internal))
        write(6,*)       '-intfile       File with ICs (for "sel")       ', trim(adjustl(intfile))
        write(6,*)       '-intfile2      File with ICs (for "sel")       ', trim(adjustl(intfile0))
        write(6,*)       '               second file for double projeciton'
        write(6,*)       '               (only meaningful with -prjS[-c])'
        write(6,*)       '-rmzfile       File deleting ICs from Zmat     ', trim(adjustl(rmzfile))
        write(6,*)       '-[no]sym       Use symmetry to form Zmat      ',  use_symmetry
        write(6,*)       '-[no]vert      Correct with B derivatives for ',  vertical
        write(6,*)       '               non-stationary points'
        write(6,*)       '-[no]prjall    Project modes with current     ', project_on_all
        write(6,*)       '               internal set on those computed '
        write(6,*)       '               with the "-intmode all" set    '
        write(6,*)       ''
        write(6,*)       ' ** Options for themochemistry **'
        write(6,*)       '-thermo        Temp (K) for thermochemistry   ', Tthermo 
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
        write(6,*)       '-alg           Algorith to convert from Cart. ',  conversion_i2c
        write(6,*)       '               to internal [zmat|iter]'
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
        character(len=*),intent(in) :: label, scan_type
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
        ! "Bonds"
        do i=1,nbonds
            k=k+1
            S(k) = S(k) + Lc(k) * Qstep
        enddo
        ! "Angles"
        do i=1,nangles
            k=k+1
            S(k) = S(k) + Lc(k) * Qstep
        enddo
        ! "Dihedrals"
        do i=1,ndihed
            k=k+1
            S(k) = S(k) + Lc(k) * Qstep
            if (S(k) >  PI) S(k)=S(k)-2.d0*PI
            if (S(k) < -PI) S(k)=S(k)+2.d0*PI
        enddo

        return
       
    end subroutine displace_Scoord


end program normal_modes_internal

