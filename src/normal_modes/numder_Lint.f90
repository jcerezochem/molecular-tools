program numder_fc


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to compute the numerical freqs from G09 scan
    !
    ! Compilation instructions (for mymake script):
    !
    ! Change log:
    !
    ! TODO:
    ! ------
    !
    ! History

    !
    !============================================================================    

!*****************
!   MODULE LOAD
!*****************
!============================================
!   Generic (structure_types independent)
!============================================
    use alerts
    use structure_types
    use line_preprocess
    use constants
    use generic_io
    use generic_io_molec
    use vibrational_analysis
    use verbosity
    use molecular_structure
    use ff_build
    use symmetry
    use internal_module

    implicit none

    integer,parameter :: NDIM=600

    !====================== 
    !Options 
    logical :: use_symmetry=.false.,   &
               include_hbonds=.false., &
               vertical=.false.,       &
               analytic_Bder=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule, molec_aux
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns
    character(len=5) :: PG
    !Job info
    character(len=20) :: calc, method, basis
    character(len=150):: title

    !====================== 
    ! PES topology and normal mode things
    real(8),dimension(1:5,1:NDIM,1:NDIM) :: LL
    real(8),dimension(1:NDIM,1:NDIM) :: Lref
    real(8),dimension(1:NDIM,1:NDIM) :: LderQj, Drot
    real(8),dimension(1:NDIM*NDIM)   :: Hlt
    real(8),dimension(1:NDIM,1:NDIM) :: Hess
    real(8),dimension(NDIM) :: Freq, Grad
    character(len=2),dimension(NDIM) :: AtName
    integer :: error
    logical :: do_Lt=.false.
    !====================== 

    !====================== 
    !INTERNAL CODE THINGS
    real(8),dimension(1:NDIM,1:NDIM) :: B, G
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Bder
    real(8),dimension(1:NDIM,1:NDIM) :: X,Xinv
    !Save definitio of the modes in character
    character(len=100),dimension(NDIM) :: ModeDef
    !VECTORS
    real(8),dimension(NDIM) :: S
    integer,dimension(NDIM) :: S_sym
    ! Switches
    character(len=5) :: def_internal="ZMAT"
    !====================== 

    real(8),dimension(5) :: R
    real(8) :: factor
    logical :: corrected
    character(len=200) :: line, dummy_char
    character(len=13)  :: marker
    character(len=3)   :: density="SCF"

    integer :: Na, Nb
    integer :: i,j,k, istep, iread, IOstatus

    !================
    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_SYM=12,  &
               I_RMF=16,  & 
               O_DER=20

    !files
    character(len=10) :: ft ="guess",  ftg="guess",  fth="guess", ftn="guess"
    character(len=200):: inpfile  ="state1.fchk", &
                         gradfile ="same", &
                         hessfile ="same", &
                         nmfile   ="none", &
                         intfile  ="none", &
                         rmzfile  ="none", &
                         symm_file="none", &
                         reffile ="reference.log"




! (End of variables declaration) 
!==================================================================================

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,reffile,density,do_Lt)
    call set_word_upper_case(density)

    if (adjustl(density) == "SCF") then
        marker="SCF Done:  E("
    elseif (adjustl(density) == "TD") then
        marker="E(TD-HF/TD-KS"
    else
        print*, "Unkown density: "//density
        stop
    endif

    ! MANAGE INTERNAL COORDS
    ! ---------------------------------
    ! Initially read structure to get the internal set
    open(I_INP,file=inpfile,status="old") 
    call generic_strmol_reader(I_INP,'log',molecule)
    Nat =molecule%natoms
    Nvib=3*Nat-6
    close(I_INP)
    ! Get connectivity 
    call guess_connect(molecule)

    ! Manage symmetry
    if (.not.use_symmetry) then
        molecule%PG="C1"
    else if (trim(adjustl(symm_file)) /= "NONE") then
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
    call define_internal_set(molecule,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns)
    if (Ns > Nvib) then
        print*, "Ns > Nvib", Ns, Nvib
        call alert_msg("fatal","Non-redundan coordinate set needs mapping (still on dev)")
        ! Need mapping from whole set to Zmat
    elseif (Ns > Nvib) then
        call alert_msg("fatal","Reduced coordinates cases still not implemented")
        ! Need to freeze unused coords to its input values
    endif

    ! Analyze numders
    open(I_INP,file=inpfile,status="old")

    istep=0
    verbose=0
    do while ( IOstatus == 0 )
        istep=istep+1
        call summary_parser(I_INP,3,line,IOstatus)
        if (IOstatus /= 0) exit
        call split_line(line,'=',dummy_char,line)
        read(line,*) R(istep)
        call rewind_summary(I_INP) 
        call generic_strmol_reader(I_INP,'log',molecule)
        call gen_bonded(molecule)
        call define_internal_set(molecule,def_internal,intfile,rmzfile,use_symmetry,isym, S_sym,Ns)
        call rewind_summary(I_INP)
        call generic_Hessian_reader(I_INP,'log',Nat,Hlt,error)
        call rewind_summary(I_INP)
        k=0
        do i=1,3*Nat
        do j=1,i
            k=k+1
            Hess(i,j) = Hlt(k)
            Hess(j,i) = Hlt(k)
        enddo 
        enddo
        call generic_gradient_reader(I_INP,'log',Nat,Grad,error)
        print*, "Reading step...", istep, R(istep)

!         call vibrations_Cart(Nat,X,Y,Z,Mass,Hlt,Nvib,LL(istep,:,:),Freq,error_flag=error)
        call internal_Wilson_new(molecule,Nvib,S,B,ModeDef)
        !SOLVE GF METHOD TO GET NM AND FREQ
        call internal_Gmetric(Nat,Nvib,molecule%atom(:)%mass,B,G)
        call calc_Bder(molecule,Nvib,Bder,.true.)
        call HessianCart2int(Nat,Nvib,Hess,molecule%atom(:)%mass,B,G,Grad=Grad,Bder=Bder)
        call gf_method(Nvib,G,Hess,LL(istep,:,:),Freq,X,Xinv)
!         call analyze_internal(Nvib,LL(istep,:,:),Freq,ModeDef)
        LL(istep,1:Nvib,1:Nvib) = inverse_realgen(Nvib,LL(istep,:,:))
        call analyze_internal(Nvib,LL(istep,:,:),Freq,ModeDef)
    enddo
    close(I_INP)
    if (istep-1/=5) call alert_msg("fatal","logfile has not the the required number of steps. "//&
                                           "Was it generated by nm_internal or nm_cartesian?") 
    print*, ""

!     call analyze_internal(Nvib,LL(1,:,:),Freq,ModeDef)

    !Assume always the sign at equil. Correct using an element for the whole row
    do istep=1,5
        do k=1,Nvib
            do i=1,Nvib
                !use an element significantly large
                corrected=.false.
                if (abs(LL(istep,i,k)) > 0.001d0) then
                    if ( LL(istep,i,k) == sign(LL(istep,i,k),LL(3,i,k)) ) then
                        factor=1.d0
                    else
                        factor=-1.d0
                    endif
                    corrected=.true.
                    exit
                endif
            enddo
            if (.not.corrected) call alert_msg("fatal","Could not correct the signs")
            LL(istep,1:Nvib,k) = LL(istep,1:Nvib,k) * factor
        enddo
    enddo


    !    L(Na,   Nb)
    ! Lder(Na,Nb,Nb)
    if (do_Lt) then
        Na = Nvib
        Nb = Nvib
        do istep=1,5
            LL(istep,:,:) = transpose(LL(istep,:,:))
        enddo
    else 
        Na = Nvib
        Nb = Nvib
    endif


    ! Small disp R(4) - R(2)
    do i=1,Na
    do k=1,Nb
        LderQj(i,k) = ( LL(4,i,k) - LL(2,i,k) ) / ( R(4) - R(2) )
    enddo
    enddo
    !print derivates
    open(O_DER,file='Ders_s_disp.dat') 
    do i=1,Na
        write(O_DER,'(600(E15.8,X))') LderQj(i,1:Nb)
    enddo
    close(O_DER)

    ! Large disp R(5) - R(1)
    do i=1,Na
    do k=1,Nb
        LderQj(i,k) = ( LL(5,i,k) - LL(1,i,k) ) / ( R(5) - R(1) )
    enddo
    enddo
    !print derivates
    open(O_DER,file='Ders_l_disp.dat') 
    do i=1,Na
        write(O_DER,'(600(E15.8,X))') LderQj(i,1:Nb)
    enddo
    close(O_DER)

    ! And write the L at equilibrium for further checking
    open(O_DER,file='Lmatrix_ref.dat') 
    do i=1,Na
        write(O_DER,'(600(E15.8,X))') LL(3,i,1:Nb)
    enddo
    close(O_DER)

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,reffile,density,do_Lt)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile, density, reffile
        logical                        :: do_Lt
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

                case ("-ref") 
                    call getarg(i+1, reffile)
                    argument_retrieved=.true.

                case ("-Lt") 
                    do_Lt=.true.
                case ("-noLt") 
                    do_Lt=.false.

                case ("-dens") 
                    call getarg(i+1, density)
                    argument_retrieved=.true.
        
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
                    print*, "Unkown command line argument: "//adjustl(arg)
                    stop
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '          NUM DERIVATIVES FROM G09log FILE '          
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-dens           ', trim(adjustl(density))
        write(0,*) '-[no]Lt        ', do_Lt
        write(0,*) '-ref           ', trim(adjustl(reffile))
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program numder_fc

