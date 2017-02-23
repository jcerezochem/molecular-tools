program get_deltaX


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
    ! NOTES
    ! For vertical, the analysis of the nm displacement can be done:
    !  * in the S-space (internal coordinates). Activated with -vert
    !  * in the Q-space (final state nm). Activated with -vert2
    ! The same convention is used to compute Er
    ! For AH, we always work in the S-space
    !==============================================================


    !*****************
    !   MODULE LOAD
    !*****************
    !============================================
    !   Generic
    !============================================
    use io
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
    use vertical_model

    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: use_symmetry=.false.   ,&
               modred=.false.         ,&
               tswitch=.false.        ,&
               symaddapt=.false.      ,&
               vertical=.false.       ,&
               verticalQspace2=.false.,&
               verticalQspace1=.false.,&
               gradcorrectS1=.false.  ,&  
               gradcorrectS2=.false.  ,&
               same_red2nonred_rotation=.true., &
               analytic_Bder=.true., &
               check_symmetry=.true., &
               orthogonalize=.false., &
               original_internal=.false., &
               force_real=.false., &
               apply_projection_matrix=.false., &
               move_to_min=.false.
    character(len=4) :: def_internal='all'
    character(len=4) :: def_internal0='defa' ! defa(ult) is "the same as working set"
    character(len=1) :: reference_frame='F'
    character(len=10):: model="adia"
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1, state2
    type(str_bonded) :: geomS, geom0
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns, NNvib, Nvib0, Ns0, NsS, Nz, Nf
    character(len=5) :: PG
    !Bonded info
    integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !B and G matrices
    real(8),dimension(NDIM,NDIM) :: B1, B2
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X, X1inv,X2inv, L1,L2,L1inv, &
                                        Asel1,Asel2,Asel1inv,Asel2inv, gBder, &
                                        G0, B0, P
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    !Duschisky
    real(8),dimension(NDIM,NDIM) :: G1, G2, Jdus
    !T0 - switching effects
    real(8),dimension(3,3) :: T
    !AUXILIAR MATRICES
    real(8),dimension(NDIM,NDIM) :: Aux, Aux2
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: Freq1, Freq2, S1, S2, Vec1, Vec2, Q0, FC
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym
    !Shifts
    real(8),dimension(NDIM) :: Delta
    real(8) :: Delta_p, Er
    ! Z-mat and redundant geoms
    integer,dimension(NDIM) :: Zmap
    type(str_bonded) :: zmatgeom, allgeom
    !====================== 

    !=========================
    ! Distance calculation stuff
    logical :: get_distances=.true.
    real(8) :: ff, f1, f0, time, dt, dist, area
    !====================== 

    !====================== 
    !
    real(8),dimension(:),allocatable :: A
    integer :: error
    character(len=200) :: msg
    character(len=5)   :: current_symm
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=16) :: dummy_char
    real(8) :: Theta, Theta2, Theta3
    ! RMZ things
    character :: rm_type
    integer :: rm_zline, Nrm, nbonds_rm, nangles_rm, ndiheds_rm
    integer,dimension(100) :: bond_rm, angle_rm, dihed_rm
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, k90,k95,k99, nn, imin, imax,&
               i1,i2,i3,i4, iop
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
               I_CNX=17,  &
               I_MAS=18,  &
               O_DUS=20,  &
               O_DIS=21,  &
               O_DMAT=22, &
               O_DUS2=23, &
               O_DIS2=24, &
               O_STAT=25
    !files
    character(len=10) :: ft ="guess",  ftg="guess",  fth="guess", &
                         ft2 ="guess", ftg2="guess", fth2="guess"
    character(len=200):: inpfile  ="state1.fchk", &
                         gradfile ="same", &
                         hessfile ="same", &
                         inpfile2 ="none", &
                         gradfile2="same", &
                         hessfile2="same", &
                         intfile  ="none", &
                         intfile0 ="default", & ! default is "the same as working set"
                         rmzfile  ="none", &
                         symm_file="none", &
                         cnx_file="guess", &
                         mass_file="none"
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

    !--------------------------
    ! Tune io
    !--------------------------
    ! Set unit for alert messages
    alert_unt=6
    !--------------------------

    !===========================
    ! Allocate atoms (default)
    call allocate_atoms(state1)
    call allocate_atoms(state2)
    !===========================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,&
                     inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                     intfile,intfile0,rmzfile,def_internal,def_internal0,&
                     use_symmetry,cnx_file,mass_file,&
!                    tswitch,
                     symaddapt,same_red2nonred_rotation,analytic_Bder,&
                     model,vertical,verticalQspace2,verticalQspace1,&
                     gradcorrectS1,gradcorrectS2,&
                     orthogonalize,original_internal,force_real,reference_frame,&
                     apply_projection_matrix,move_to_min)
    call set_word_upper_case(def_internal)
    call set_word_upper_case(reference_frame)
    call set_word_upper_case(model)

    ! 1. INTERNAL VIBRATIONAL ANALYSIS ON STATE1 AND STATE2

    !===========
    !State 1
    !===========
    if (verbose>0) then
        print*, ""
        print*, "=========="
        print*, " STATE 1"
        print*, "=========="
    endif
 
    ! READ DATA (each element from a different file is possible)
    ! ---------------------------------
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)


    ! STRUCTURE FILE
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,state1,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (State1)")
    close(I_INP)
    ! Shortcuts
    Nat = state1%natoms
    Nvib = 3*Nat-6

    ! Read mass from file if given
    if (adjustl(mass_file) /= "none") then
        print'(/,X,A)', "Reading atomic masses from: "//trim(adjustl(mass_file))
        open(I_MAS,file=mass_file,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(mass_file)) )
        do i=1,Nat
            read(I_MAS,*,iostat=IOstatus) state1%atom(i)%mass 
            if (IOstatus /= 0) call alert_msg( "fatal","While reading "//trim(adjustl(mass_file)) )
        enddo
        close(I_MAS)
    endif


    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! Get connectivity 
    if (cnx_file == "guess") then
        call guess_connect(state1)
    else
        print'(/,A,/)', "Reading connectivity from file: "//trim(adjustl(cnx_file))
        open(I_CNX,file=cnx_file,status='old')
        call read_connect(I_CNX,state1)
        close(I_CNX)
    endif

    ! Manage symmetry
    if (.not.use_symmetry) then
        state1%PG="C1"
    else if (trim(adjustl(symm_file)) /= "none") then
        msg = "Using custom symmetry file: "//trim(adjustl(symm_file)) 
        call alert_msg("note",msg)
        open(I_SYM,file=symm_file)
        do i=1,state1%natoms
            read(I_SYM,*) j, isym(j)
        enddo
        close(I_SYM)
        !Set PG to CUStom
        state1%PG="CUS"
    else
        state1%PG="XX"
        call symm_atoms(state1,isym)
    endif

    !---------
    print'(/,X,A)', "GETTING INTERNAL SET "
    print*, "------------------------------------------"
    ! Refress connectivity
    call gen_bonded(state1)
    ! Define internal set
    call define_internal_set(state1,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns,Nf,Aux2)
    !From now on, we'll use atomic units
    call set_geom_units(state1,"Bohr")
    call compute_internal(state1,Ns,S1,ModeDef)

    !===========
    !State 2
    !===========
    if (verbose>0) then
        print*, ""
        print*, "=========="
        print*, " STATE 2"
        print*, "=========="
    endif

 
    ! READ DATA (each element from a different file is possible)
    ! ---------------------------------
    !Guess filetypes
    if (ft2 == "guess") &
    call split_line_back(inpfile2,".",null,ft2)


    ! STRUCTURE FILE
    open(I_INP,file=inpfile2,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile2)) )
    call generic_strmol_reader(I_INP,ft2,state2,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (State2)")
    close(I_INP)
    if (Nat /= state2%natoms) call alert_msg("fatal","Initial and final states don't have the same number of atoms.")

    ! Read mass from file if given
    if (adjustl(mass_file) /= "none") then
        print'(/,X,A)', "Reading atomic masses from: "//trim(adjustl(mass_file))
        open(I_MAS,file=mass_file,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(mass_file)) )
        do i=1,Nat
            read(I_MAS,*,iostat=IOstatus) state2%atom(i)%mass 
            if (IOstatus /= 0) call alert_msg( "fatal","While reading "//trim(adjustl(mass_file)) )
        enddo
        close(I_MAS)
    endif
    
    
    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! Get connectivity 
    if (cnx_file == "guess") then
        call guess_connect(state2)
    else
        print'(/,A,/)', "Reading connectivity from file: "//trim(adjustl(cnx_file))
        open(I_CNX,file=cnx_file,status='old')
        call read_connect(I_CNX,state2)
        close(I_CNX)
    endif
    
    ! Manage symmetry
    if (.not.use_symmetry) then
        state2%PG="C1"
    else if (trim(adjustl(symm_file)) /= "none") then
        msg = "Using custom symmetry file: "//trim(adjustl(symm_file)) 
        call alert_msg("note",msg)
        open(I_SYM,file=symm_file)
        do i=1,state2%natoms
            read(I_SYM,*) j, isym(j)
        enddo
        close(I_SYM)
        !Set PG to CUStom
        state2%PG="CUS"
    else
        state2%PG="XX"
        call symm_atoms(state2,isym)
    endif
    if (state1%PG /= state2%PG) then
        print*, "PG(State1): ", state1%PG
        print*, "PG(State2): ", state2%PG
        if (.not.use_symmetry) then
            call alert_msg("note","Initial and final state have different symmetry")
        else
            call alert_msg("warning","Initial and final state have different symmetry")
        endif
    endif
    
    
    !---------
    print'(/,X,A)', "GETTING INTERNAL SET "
    print*, "------------------------------------------"
    ! Refress connectivity
    call gen_bonded(state2)
    ! Define internal set
    call define_internal_set(state2,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns,Nf,Aux2)
    !From now on, we'll use atomic units
    call set_geom_units(state2,"Bohr")
    call compute_internal(state2,Ns,S2,ModeDef)
    
        
    !--------------------------
    ! ICs displacement
    !--------------------------
    if (verbose>0) then
        print*, ""
        print*, "=========================="
        print*, " SHIFTS (internal coord)"
        print*, "=========================="
    endif

    ! Bonds
    do i=1,state1%geom%nbonds
        Delta(i) = S2(i)-S1(i)
    enddo
    ! Angles
    do j=i,i+state1%geom%nangles-1
        Delta(j) = S2(j)-S1(j)
    enddo
    ! Dihedrals
    do k=j,j+state1%geom%ndihed-1
        Delta(k) = S2(k)-S1(k)
        Delta_p = S2(k)-S1(k)+2.d0*PI
        if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
        Delta_p = S2(k)-S1(k)-2.d0*PI
        if (dabs(Delta_p) < dabs(Delta(k))) Delta(k)=Delta_p
    enddo
    ! Impropers
    do i=k,k+state1%geom%nimprop-1
        Delta(i) = S2(i)-S1(i)
        Delta_p = S2(i)-S1(i)+2.d0*PI
        if (dabs(Delta_p) < dabs(Delta(i))) Delta(i)=Delta_p
        Delta_p = S2(i)-S1(i)-2.d0*PI
        if (dabs(Delta_p) < dabs(Delta(i))) Delta(i)=Delta_p
    enddo


    print*, "List of ICs and Deltas"
    if (state1%geom%nbonds /= 0) &
     print'(X,A,/,X,A)', "Bonds (Angs)",&
                         "  IC   Description      State2    State1     Delta"
    do i=1,state1%geom%nbonds
        print'(I5,X,A,X,3(F8.3,2X))', i,trim(ModeDef(i)),S2(i)*BOHRtoANGS,S1(i)*BOHRtoANGS,Delta(i)*BOHRtoANGS
    enddo
    if (state1%geom%nangles /= 0) &
     print'(X,A,/,X,A)', "Angles (deg)",&
                         "  IC   Description                State2    State1     Delta"
    do j=i,i+state1%geom%nangles-1
        print'(I5,X,A,X,3(F8.2,2X))', j, trim(ModeDef(j)), S2(j)*180.d0/PI,S1(j)*180.d0/PI,Delta(j)*180.d0/PI
    enddo
    if (state1%geom%ndihed /= 0) &
     print'(X,A,/,X,A)', "Dihedrals (deg)",&
                         "  IC   Description                          State2    State1     Delta"
    do k=j,j+state1%geom%ndihed-1
        print'(I5,X,A,X,3(F8.2,2X))', k, trim(ModeDef(k)), S2(k)*180.d0/PI,S1(k)*180.d0/PI,Delta(k)*180.d0/PI
    enddo
    if (state1%geom%nimprop /= 0) &
     print'(X,A,/,X,A)', "Impropers (deg)",&
                         "  IC   Description                          State2    State1     Delta"
    do i=k,k+state1%geom%nimprop-1
        print'(I5,X,A,X,3(F8.2,2X))', i, trim(ModeDef(i)), S2(i)*180.d0/PI,S1(i)*180.d0/PI,Delta(i)*180.d0/PI
    enddo

    ! Original DeltaX
    call set_geom_units(state1,"Bohr")
    call set_geom_units(state2,"Bohr")
    print*, ""
    print*, "ORIGINAL DELTA X (Bohr)"
    do i=1,Nat 
        print*, state2%atom(i)%x - state1%atom(i)%x
        print*, state2%atom(i)%y - state1%atom(i)%y
        print*, state2%atom(i)%z - state1%atom(i)%z
    enddo
    print*, ""

    ! DeltaX from DeltaS
    state2=state1
    call intshif2cart(state2,Delta,thr_set=1.d-12)
    print*, ""
    print*, "DELTA X FROM DELTA S (Bohr)"
    do i=1,Nat 
        print*, state2%atom(i)%x - state1%atom(i)%x
        print*, state2%atom(i)%y - state1%atom(i)%y
        print*, state2%atom(i)%z - state1%atom(i)%z
    enddo

    call summary_alerts

    call cpu_time(tf)
    write(6,'(/,A,X,F12.3,/)') "CPU (s) for internal vib analysis: ", tf-ti

    stop

    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,hessfile,fth,gradfile,ftg,&
                           inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                           intfile,intfile0,rmzfile,def_internal,def_internal0,&
                           use_symmetry,cnx_file,mass_file,&
!                          tswitch,
                           symaddapt,same_red2nonred_rotation,analytic_Bder,&
                           model,vertical,verticalQspace2,verticalQspace1,&
                           gradcorrectS1,gradcorrectS2,&
                           orthogonalize,original_internal,force_real,reference_frame,&
                           apply_projection_matrix,move_to_min)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile,ftg,&
                                          inpfile2,ft2,hessfile2,fth2,gradfile2,ftg2,&
                                          intfile,intfile0,rmzfile,def_internal,def_internal0,&
                                          cnx_file,reference_frame,model,mass_file !, symfile
        logical,intent(inout)          :: use_symmetry, vertical, verticalQspace2, &
                                          verticalQspace1, &
                                          gradcorrectS1, gradcorrectS2, symaddapt, &
                                          same_red2nonred_rotation,analytic_Bder, &
                                          orthogonalize,original_internal,force_real, &
                                          apply_projection_matrix,move_to_min
!         logical,intent(inout) :: tswitch

        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        character(len=500) :: input_command
        character(len=10)  :: MODEL_UPPER
        ! iargc type must be specified with implicit none (strict compilation)
        integer :: iargc

        ! Tune defaults
        logical :: gradcorrectS1_default=.true., &
                   gradcorrectS2_default=.true.

        !Initialize input_command
        call getarg(0,input_command)
        !Get input flags
        do i=1,iargc()
            call getarg(i,arg)
            input_command = trim(adjustl(input_command))//" "//trim(adjustl(arg))
        enddo

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

                case ("-f2") 
                    call getarg(i+1, inpfile2)
                    argument_retrieved=.true.
                case ("-ft2") 
                    call getarg(i+1, ft2)
                    argument_retrieved=.true.

                case ("-fhess2") 
                    call getarg(i+1, hessfile2)
                    argument_retrieved=.true.
                case ("-fth2") 
                    call getarg(i+1, fth2)
                    argument_retrieved=.true.

                case ("-fgrad2") 
                    call getarg(i+1, gradfile2)
                    argument_retrieved=.true.
                case ("-ftg2") 
                    call getarg(i+1, ftg2)
                    argument_retrieved=.true.

                case ("-fmass") 
                    call getarg(i+1, mass_file)
                    argument_retrieved=.true.

                case ("-cnx") 
                    call getarg(i+1, cnx_file)
                    argument_retrieved=.true.

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-intfile0") 
                    call getarg(i+1, intfile0)
                    argument_retrieved=.true.

                case ("-rmzfile") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.
                ! Kept for backward compatibility (but replaced by -rmzfile)
                case ("-rmz") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.

                case ("-intmode0")
                    call getarg(i+1, def_internal0)
                    argument_retrieved=.true.

                case ("-intmode")
                    call getarg(i+1, def_internal)
                    argument_retrieved=.true.
                ! Kept for backward compatibility (but replaced by -intmode)
                case ("-intset")
                    call getarg(i+1, def_internal)
                    argument_retrieved=.true.

                case ("-sym")
                    use_symmetry=.true.
                case ("-nosym")
                    use_symmetry=.false.

                case ("-symaddapt")
                    symaddapt=.true.

                case("-samerot")
                    same_red2nonred_rotation=.true.
                case("-nosamerot")
                    same_red2nonred_rotation=.false.

                ! Options to tune the model
                !================================================================
                ! This is the new (and now standard way to get the model)
                case ("-model")
                    call getarg(i+1, model)
                    argument_retrieved=.true.
                !The others are kept for backward compatibility
!                 case ("-vertQ1")
!                     vertical=.true.
!                     verticalQspace2=.false.
!                     verticalQspace1=.true.
!                 case ("-vertQ2")
!                     vertical=.true.
!                     verticalQspace2=.true.
!                     verticalQspace1=.false.
!                 case ("-vert")
!                     vertical=.true.
!                     verticalQspace2=.false.
!                 case ("-novert")
!                     vertical=.false.
!                     verticalQspace2=.false.
                !================================================================

                case ("-ref") 
                    call getarg(i+1, reference_frame)
                    argument_retrieved=.true.

                case ("-force-real")
                    force_real=.true.
                case ("-noforce-real")
                    force_real=.false.

                case ("-prj-tr")
                    apply_projection_matrix=.true.
                case ("-noprj-tr")
                    apply_projection_matrix=.false.

                case ("-orth")
                    orthogonalize=.true.
                case ("-noorth")
                    orthogonalize=.false.
                    
                case ("-origint")
                    original_internal=.true.
                case ("-noorigint")
                    original_internal=.false.
        
                case ("-h")
                    need_help=.true.

                case ("-corrS2")
                    gradcorrectS2=.true.
                    gradcorrectS2_default=.false.
                case ("-nocorrS2")
                    gradcorrectS2=.false.
                    gradcorrectS2_default=.false.
                case ("-corrS1")
                    gradcorrectS1=.true.
                    gradcorrectS1_default=.false.
                case ("-nocorrS1")
                    gradcorrectS1=.false.
                    gradcorrectS1_default=.false.

                !HIDDEN

                case ("-anaBder")
                    analytic_Bder=.true.
                case ("-noanaBder")
                    analytic_Bder=.false.

                ! Control verbosity
                case ("-quiet")
                    verbose=0
                case ("-concise")
                    verbose=1
                case ("-v")
                    verbose=2
                case ("-vv")
                    verbose=3

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       ! Manage defaults
       ! If not declared, hessfile and gradfile are the same as inpfile
       if (adjustl(hessfile) == "same") then
           hessfile=inpfile
           if (adjustl(fth) == "guess")  fth=ft
       endif
       if (adjustl(gradfile) == "same") then
           gradfile=inpfile
           if (adjustl(ftg) == "guess")  ftg=ft
       endif
       if (adjustl(hessfile2) == "same") then
           hessfile2=inpfile2
           if (adjustl(fth2) == "guess")  fth2=ft2
       endif
       if (adjustl(gradfile2) == "same") then
           gradfile2=inpfile2
           if (adjustl(ftg2) == "guess")  ftg2=ft2
       endif

       ! Set old options for the model with the new input key
       MODEL_UPPER = adjustl(model)
       call set_word_upper_case(MODEL_UPPER)
       select case (adjustl(MODEL_UPPER))
           case ("ADIA") 
               vertical=.false.
               verticalQspace1=.false.
               verticalQspace2=.false.
           case ("AS") 
               vertical=.false.
               verticalQspace1=.false.
               verticalQspace2=.false.
           case ("VERT") 
               vertical=.true.
               verticalQspace1=.false.
               verticalQspace2=.false.
           case ("VERT-A") 
               vertical=.true.
               verticalQspace1=.false.
               verticalQspace2=.false.
               move_to_min=.true.
           case ("VERTQ1") 
               vertical=.true.
               verticalQspace1=.true.
               verticalQspace2=.false.
           case ("VERTQ2") 
               vertical=.true.
               verticalQspace1=.false.
               verticalQspace2=.true.
           case default
               call alert_msg("warning","Unkown model to describe PESs: "//adjustl(model))
               need_help=.true.
       end select

       if (gradcorrectS2_default) then
           ! All vertical models include the correction
           if (vertical) then
               gradcorrectS2=.true.
           else
               gradcorrectS2=.false.
           endif
       endif
       if (gradcorrectS1_default) then
           gradcorrectS1=.false.
       endif

       ! Take defaults for the internal set for correction only
       if (def_internal0 == "defa") def_internal0=def_internal
       if (adjustl(intfile0)=="default") intfile0=intfile

       !Print options (to stdout)
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '        I N T E R N A L   D U S C H I N S K Y '    
        write(6,'(/,A)') '         Duschinski analysis for Adiabatic and    '
        write(6,'(A,/)') '        Vertical model in internal coordinates          '   
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag         Description                   Value'
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,*) '-f           Input file (State1)           ', trim(adjustl(inpfile))
        write(6,*) '-ft          \_ FileType                   ', trim(adjustl(ft))
        write(6,*) '-fhess       Hessian(S1) file              ', trim(adjustl(hessfile))
        write(6,*) '-fth         \_ FileType                   ', trim(adjustl(fth))
        write(6,*) '-fgrad       Gradient(S1) file             ', trim(adjustl(gradfile))
        write(6,*) '-ftg         \_ FileType                   ', trim(adjustl(ftg))
        write(6,*) '-f2          Input file (State2)           ', trim(adjustl(inpfile2))
        write(6,*) '-ft2         \_ FileType                   ', trim(adjustl(ft2))
        write(6,*) '-fhess2      Hessian(S2) file              ', trim(adjustl(hessfile2))
        write(6,*) '-fth2        \_ FileType                   ', trim(adjustl(fth2))
        write(6,*) '-fgrad2      Gradient(S2) file             ', trim(adjustl(gradfile2))
        write(6,*) '-ftg2        \_ FileType                   ', trim(adjustl(ftg2))
        write(6,*) '-fmass       Mass file                     ', trim(adjustl(mass_file))
        write(6,*) ''
        write(6,*) ' ** Options for state_files ** '
        write(6,*) '-ref         Reference state to output the ', reference_frame
        write(6,*) '             L matrices in its Cartesian '
        write(6,*) '             frame [I|F]'
        write(6,*) '-[no]force-real Turn imaginary frequences ', force_real
        write(6,*) '              to real (also affects Er)'
        write(6,*) '               '                       
        write(6,*) ' ** Options Internal Coordinates **           '
        write(6,*) '-[no]prj-tr  Apply projection matrix to   ', apply_projection_matrix
        write(6,*) '             rotate Grad and Hess'
        write(6,*) '-cnx         Connectivity [filename|guess] ', trim(adjustl(cnx_file))
        write(6,*) '-intmode0    Internal set:[zmat|sel|all]   ', trim(adjustl(def_internal0))
        write(6,*) '-intfile0    File with ICs (for "sel")     ', trim(adjustl(intfile0))
        write(6,*) '-intmode     Internal set:[zmat|sel|all]   ', trim(adjustl(def_internal))
        write(6,*) '-intfile     File with ICs (for "sel")     ', trim(adjustl(intfile))
        write(6,*) '-rmzfile     File deleting ICs from Zmat   ', trim(adjustl(rmzfile))
        write(6,*) '-[no]sym     Use symmetry to form Zmat    ',  use_symmetry
        write(6,*) '-[no]samerot Use S1 red->non-red rotation ',  same_red2nonred_rotation
        write(6,*) '             for S2'
        write(6,*) '-[no]orth    Use orthogonalized internals ',  orthogonalize
        write(6,*) '-[no]origint Use originally defined inter-',  original_internal
        write(6,*) '             nal w\o linear combinations  '
        write(6,*) '             (valid for non-redundant sets)'
        write(6,*) '               '
        write(6,*) ' ** Options Vertical Model **'
        write(6,*) '-model       Model for harmonic PESs       ', trim(adjustl(model))
        write(6,*) '             [vert|vertQ1|vertQ2|vert-a|adia]     '
        write(6,*) '-[no]corrS1  Correct S1 for non-stationary ', gradcorrectS1
        write(6,*) '-[no]corrS2  Correct S2 for non-stationary ', gradcorrectS2
        write(6,*) '               '
        write(6,*) '-h           Display this help            ',  need_help
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(A)') 'Input command:'
        write(6,'(A)') trim(adjustl(input_command))   
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                       'Verbose level:  ', verbose        
        write(6,'(A)') '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


end program get_deltaX

