program normal_modes_Cartesian

    !Compilation
    ! $FC ../modules/constants_mod.f90 ../modules/alerts.f90 ../modules/line_preprocess.f90 ../modules/structure_types_v2.f90 ../modules/gro_manage_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/geom_meter_v2.f90 ../modules/xyz_manage.f90 normal_modes_Cartesian_v1b.f90 -o normal_modes_Cartesian_v1b.exe -cpp -DDOUBLE 

    !NOTES
    ! Normal modes matrix indeces: T(1:Nvib,1:3Nat)
    ! This might be the contrary as the usual convention
    ! (anywaym who cares, since we use a Tvector)

    !V1b: from normal_mode_animation_v1b. Naming systematic in compliance with internal version
    !
    !Addapted to v4 release (distribution upgrade). 25/02/14
    !Additional Changes
    ! -Added filetype support with generic_strfile_read 


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
    use gaussian_manage_notypes
    use gaussian_fchk_manage
    use gaussian_manage_lowlevel
    use xyz_manage
    use fcc_manage
!   Structural parameters
    use molecular_structure
    use ff_build
!   Bond/angle/dihed meassurement
    use atomic_geom
!   vibrational analysis
    use vibrational_analysis

    type(str_resmol) :: molec

    !Interesting info..
    integer :: Nat, Nvib
    real(8),dimension(1:600) :: GEOM, RedMass, Freq, Grad
    real(8),dimension(1:600,1:600) :: L,Lt, H, Aux, Lplus, Lminus
    real(8),dimension(1:600,1:600,1:600) :: LL
    real(8),dimension(1:180300) :: Hlt
    character(len=2) :: atname 
    integer          :: Zat

    real(8) :: Amplitude, qcoord

    logical :: is_g09_modetype
    

    !Auxiliars
    character(1) :: null
    character(len=50) :: dummy_char
    real(8) :: Tmwc_ij
    ! Reading FCHK part
    character(len=10000) :: section
    integer :: N, N_T
    character :: dtype, cnull
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    integer :: error

    !Counters
    integer :: i,j,k, jdh,idh,istep, iat, kk

    !I/O
    integer :: I_INP=10,  &
               O_GRO=20,  &
               O_G09=21,  & 
               O_Q  =22,  &
               O_NUM=23,  &
               O_LIS=24,  &
               O_G96=25,  &
               S_VMD=30
    !files
    character(len=20) :: filetype="guess", ft
    character(len=200):: numderfile="none"
    character(len=200):: state1 ="state_initial.fchk", state2="state_vert.fchk"
    !Control of stdout
    logical :: verbose=.false.

    !MOVIE things
    integer :: movie_cycles, movie_steps

    !===============================================

    !Get options from command line
    call parse_input(state1,state2,numderfile,filetype,verbose)


    ! 1. CARTESIAN VIBRATIONAL ANALYSIS -- done by G09 (just read it) 
 
    ! 1. READ DATA
    ! ---------------------------------
    ! A. Initial state vibrational analysis
    open(I_INP,file=state1,status='old',iostat=IOstatus)
    ! Structure
    if (adjustl(filetype) == "guess") &
    call split_line_back(state1,".",null,filetype)
    call generic_strfile_read(I_INP,filetype,molec)

    ! Hessian (FCHK file)    
    call read_fchk(I_INP,"Cartesian Force Constants",dtype,N,A,IA,error)
    if (error == 0) then
        Hlt(1:N) = A(1:N)
        deallocate(A)
    else
        print*, "ERROR: 'Cartesian Force Constants' not found in fchk"
        stop     
    endif     

    Nat=molec%natoms
    call diag_int(Nat,                &
                  molec%atom(:)%x,    &
                  molec%atom(:)%y,    &
                  molec%atom(:)%z,    &
                  molec%atom(:)%mass, &
                  Hlt,                &
                  Nvib,               &
                  L,                  &
                  Freq,               &
                  error)
    close(I_INP)

    print*, ""
    print*, "FREQUENCIES (cm-1) - initial state"
    do i=1,Nvib
        print'(F12.4)', dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
    enddo
    print*, ""


    ! B. Final state vertical Hessian and gradient
    open(I_INP,file=state2,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(state2)) )

    ! Hessian (FCHK file)    
    call read_fchk(I_INP,"Cartesian Force Constants",dtype,N,A,IA,error)
    if (error == 0) then
        Hlt(1:N) = A(1:N)
        deallocate(A)
    else
        print*, "ERROR: 'Cartesian Force Constants' not found in fchk"
        stop     
    endif
    ! Massweight the Hessian (we keep on using UMA not UA, consistent as far as we don't mixup things)
    ! and store in a matrix array
    k=0
    do i=1,3*Nat
    do j=1,i
        k=k+1
        ii = (i-1)/3+1
        jj = (j-1)/3+1
        H(i,j) = Hlt(k)/sqrt(molec%atom(ii)%mass*molec%atom(jj)%mass) 
        H(j,i) = H(i,j)
    enddo 
    enddo


    ! Gradient (FCHK file)    
    call read_fchk(I_INP,"Cartesian Gradient",dtype,N,A,IA,error)
    if (error == 0) then
        Grad(1:N) = A(1:N)
        deallocate(A)
    else
        print*, "ERROR: 'Cartesian Gradient' not found in fchk"
        stop     
    endif

    close(I_INP)


    ! C. Numerical derivates of L
    open(I_INP,file=numderfile,status='old',iostat=IOstatus)
    do i=1,3*Nat
        ii = (i-1)/3+1
        ! Read Hessian from log (x+ / y+ / z+)
        call summary_parser(I_INP,6,section,error)
        read(section,*) Hlt(1:3*Nat*(3*Nat+1)/2) 
        call diag_int(Nat,                &
                      molec%atom(:)%x,    &
                      molec%atom(:)%y,    &
                      molec%atom(:)%z,    &
                      molec%atom(:)%mass, &
                      Hlt,Nvib,Lplus,Freq,error)
        print*, ""
        print*, "Freq. +IC", i
        do k=1,1
            print'(F12.4)', dsign(dsqrt(dabs(Freq(k))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(k))
        enddo
        print*, ""
        ! Read Hessian from log (x- / y- / z-)
        call summary_parser(I_INP,6,section,error)
        read(section,*) Hlt(1:3*Nat*(3*Nat+1)/2) 
        call diag_int(Nat,                &
                      molec%atom(:)%x,    &
                      molec%atom(:)%y,    &
                      molec%atom(:)%z,    &
                      molec%atom(:)%mass, &
                      Hlt,Nvib,Lminus,Freq,error)
        print*, ""
        print*, "Freq. -IC", i
        do k=1,1
            print'(F12.4)', dsign(dsqrt(dabs(Freq(k))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(k))
        enddo
        print*, ""
        ! LL = LLcart m^-1/2
        do j=1,3*Nat
        do k=1,Nvib
            LL(i,j,k) = (Lplus(j,k) - Lminus(j,k))/2.d-3*BOHRtoANGS/sqrt(molec%atom(ii)%mass)  ! dist in mwc
        enddo
        enddo
!         print*, ""
!         print'(X,A,3(I0,X),A,G12.4)', "LL(", i,i,1, ")=", LL(i,j,k)
!         print*, "---------------"
!         print*, ""
         
    enddo


    !---------------------------------------------
    ! Hessian in initial-state nm coordinates
    !---------------------------------------------
    if (vertical) then
        ! 1. Gradient:   g^Q = L^t m^-1/2 g^x
        do i=1,Nvib
            Freq(i) = 0.0d0
            do k=1,3*Nat
                kk = (k-1)/3+1
                Freq(i) = Freq(i) + L(k,i) / dsqrt(molec%atom(kk)%mass) * Grad(k)
            enddo
        enddo
        Grad(1:Nvib) = Freq(1:Nvib)
        print*, ""
        print*, "GRAD"
        print'(F12.4)', Grad(1:Nvib)
        print*, ""
        
        ! 2. Hessian
        ! * Correct for gradient: H^q - LL g^Q
        do i=1,3*Nat
        do j=1,3*Nat
            Aux(i,j) = 0.d0
            do k=1,Nvib
                Aux(i,j) = Aux(i,j) + Grad(k) * LL(j,i,k)
            enddo
!             Aux(i,j) = 0.d0
            H(i,j) = H(i,j) - Aux(i,j)
        enddo
        enddo
    endif
    ! * Hessian in initial-state nm:  H^Q = Lt (H^q - LL g^Q) L
    Aux(1:3*Nat,1:Nvib) = matmul(H(1:3*Nat,1:3*Nat),L(1:3*Nat,1:Nvib))
    Lt = transpose(L)
    H(1:Nvib,1:Nvib)    = matmul(Lt(1:Nvib,1:3*Nat),Aux(1:3*Nat,1:Nvib))
    ! * Diagonalize Hessian: autovectors form Duschinsky matrix, autovalues are the freqs
    call diagonalize_full(H(1:Nvib,1:Nvib),Nvib,L(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
    
    print*, ""
    print*, "FREQUENCIES (cm-1) - vertical state"
    do i=1,Nvib
        print'(F12.4)', dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
    enddo
    print*, ""


    stop

    !==============================================
    contains
    !=============================================

    subroutine parse_input(state1,state2,numderfile,filetype,verbose)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: state1,state2,numderfile,filetype
        logical,intent(inout) ::  verbose
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        real(8) :: maxd

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-f") 
                    call getarg(i+1, state1)
                    argument_retrieved=.true.
                case ("-ft") 
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.

                case ("-vf") 
                    call getarg(i+1, state2)
                    argument_retrieved=.true.

                case ("-numder") 
                    call getarg(i+1, numderfile)
                    argument_retrieved=.true.

                case ("-v")
                    verbose=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '               VIBRATIONS '    
        write(0,'(/,A)') '       ' 
        write(0,'(/,A)') '         Revision: **'
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(state1))
        write(0,*) '-ft             ', trim(adjustl(filetype))
        write(0,*) '-vf             ', trim(adjustl(state2))
        write(0,*) '-numder         ', trim(adjustl(numderfile))
        write(0,*) '-v             ', verbose
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(in) :: filetype
        type(str_resmol),intent(inout) :: molec

        select case (adjustl(filetype))
            case("fchk")
             call read_fchk_geom(I_INP,molec)
             call atname2element(molec)
!              call assign_masses(molec) !read_fchk_geom includes the fchk masses
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select


        return


    end subroutine generic_strfile_read


end program normal_modes_Cartesian

