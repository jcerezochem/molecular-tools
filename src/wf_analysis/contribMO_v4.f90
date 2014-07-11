program contribMO

    !==============================================================
    ! This code uses MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    ! Description
    !  Utility to compute the contribution of atom (label=1) in the MO
    ! Compilation intructions:
    !  gfortran ../modules/constants_mod.f90 ../modules/alerts.f90 ../modules/structure_types_v2_ALLOC.f90 ../modules/allocation_mod.f90 ../modules/line_preprocess.f90 ../modules/gaussian_manage_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/atomic_geom_v2.f90 ../modules/MatrixMod.f90 contribMO_v4.f90 -o contribMO_v4.exe -cpp -llapack
    !
    ! History
    ! v2 (20/07/2012):
    !       -at changed to -frg: the file may contain now several fragments
    !       NTO specific feature: contributions for all NTO and average.
    ! v3 (28/06/2013)
    !       Enlarge memory
    ! v4 (6/11/2013)
    !       Disable NTO averaging. This is not conceptually clear. If there is not one predominant NTO,
    !        NTO analysis is just meaningless (or is it?)
    !        Now it gives a warning if NTO coef is < 0.9
    !***************************
    ! v4.4 (Adapt to distributed modules)
    !       Change real(4) to real(8)
    !==============================================================

    use structure_types
    use gaussian_manage
    use gaussian_fchk_manage !, only : read_MO
    use atomic_geom

    implicit none

    integer,parameter :: BASIS_SIZE = 2000

    real(8),dimension(1:BASIS_SIZE,1:BASIS_SIZE) :: MO
    integer,dimension(1:BASIS_SIZE) :: AOmap
    
    integer :: IFCHK =11, &
               IMO   =12, &
               IAT   =13
    character(len=50) :: fchkfile,      &
                         mo_file="all", &
                         at_file="ndx"

    integer :: N,M, i, j, k, ii, frg

    real(8) :: coeff_total, coeff_metal

    integer :: Metal_Index = 1
    character(len=5) :: mo_type

    !OM selection stuff
    integer :: Nselect, Nselect_alpha, Nselect_beta, Nfrg
    integer,dimension(1:10) :: Nat
    integer,dimension(1:10,1:500) :: atom_index
    integer,dimension(1:BASIS_SIZE) :: MO_select_alpha, MO_select_beta

    !NTO averaging
    integer :: Nel
    real(8) :: Norm
    real(8),dimension(1:10) :: hole, particle
    real(8),dimension(1:BASIS_SIZE) :: W, W2
    logical :: nto_av = .false., do_nto=.false.

    call parse_input(fchkfile,mo_file,at_file,Metal_Index,do_nto)

    open(IFCHK,file=fchkfile,status="old")

    call read_MO(IFCHK,MO,N,"alpha",-1)
    M = int ( sqrt(float(N)) )
    call AO_atom_map(IFCHK,AOmap(1:M))
    
    !Atom selection
    if (adjustl(at_file) == "ndx") then
        Nfrg   = 1
        Nat(1) = 1
        atom_index(1,1) = Metal_Index
    else
        open(IAT,file=at_file,status="old")
        read(IAT,*) Nfrg
        do i=1,Nfrg
            read(IAT,*) Nat(i)
            read(IAT,*) atom_index(i,1:Nat(i))
        enddo
    endif
    
    
    !Enable OM selection (unless do_nto is enabled -- in this case only the first one is analyzed)
    if (adjustl(mo_file) == "all") then ! .or. nto_av) then
        Nselect_alpha = M
        Nselect_beta  = M
        MO_select_alpha(1:M) = (/(i,i=1,M)/)
        MO_select_beta(1:M)  = (/(i,i=1,M)/)
    else if (adjustl(mo_file) == "alpha") then
        Nselect_alpha = M
        Nselect_beta  = 0
        MO_select_alpha(1:M) = (/(i,i=1,M)/)
    else
        open(IMO,file=mo_file,status="old")
        read(IMO,*) Nselect
        Nselect_alpha = 0
        Nselect_beta  = 0
        do i=1,Nselect
            read(IMO,*) mo_type, ii
            if (adjustl(mo_type) == "alpha") then 
                Nselect_alpha = Nselect_alpha + 1
                MO_select_alpha(Nselect_alpha) = ii
            else if (adjustl(mo_type) == "beta") then
                Nselect_beta  = Nselect_beta + 1
                MO_select_beta(Nselect_beta)   = ii
            endif
        enddo
    endif

    !Alpha
    !=====
    !NTO averaging: get NTO coefficients (W) (stored in MO energies section)
    ! note that after this call M=N --- disabled
    if (do_nto) then
        write(0,*) "=============="
        write(0,*) " ANALISIS NTO"
        write(0,*) "=============="
        write(0,*) "NOTE: NTO averaging has been disabled, the analysis will be done only"
        write(0,*) "      on the first NTO orbital. This is conceptually not clear. Think" 
        write(0,*) "      about it and report back, thanks."

        call read_NrEl(IFCHK,Nel,"alpha")
        call read_MOenergy(IFCHK,W,N,"alpha",-1)
!         !Normalize
        Norm = 0.0
        do i=1,Nel
            Norm = Norm + W(i)
        enddo
        do i=1,Nel
            W2(i) = W(i)/Norm
            !The particle is equal
            j = 2*Nel + 1 - i
            W2(j) = W2(j)/Norm
        enddo
! 
!         !Set the first NTO orbital
!         Nselect_alpha = 2
!         MO_select_alpha(1) = Nel
!         MO_select_alpha(2) = Nel+1
    endif


    do frg = 1, Nfrg
        hole(frg) = 0.0
        particle(frg) = 0.0

        do ii=1,Nselect_alpha
            i = MO_select_alpha(ii)
            ! Compute contributions
            coeff_total = 0.0
            coeff_metal = 0.0
            do j=1,M
                coeff_total = coeff_total + abs(MO(j,i))
                !Contribution of fragments
                do k=1,Nat(frg)
                    if ( AOmap(j) == atom_index(frg,k) ) &
                     coeff_metal = coeff_metal + abs(MO(j,i))
                enddo
            enddo
            if (.not.do_nto) &
            print'(X,A,I3,A,I3,A,I3,A,I4,A, F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to MO alpha", i,":",coeff_metal/coeff_total * 100
            !NTO av: store in array
            if (do_nto) then
                if (i == Nel) then
                    hole(frg)     = hole(frg) + coeff_metal/coeff_total * 100.0 !* W(i) -- averaging disabled
                else if (i == Nel+1) then
                    particle(frg) = particle(frg) + coeff_metal/coeff_total * 100.0 !* W(i) -- averaging disabled
                endif
            endif

        enddo

    enddo

    !NTO:
    if (do_nto) then
        print*, "Sum of Alpha Coefficients: ", Norm
        print*, "First NTO coef           : ", W(Nel)!, W(Nel+1) 
        do frg = 1, Nfrg
            print'(X,A,I3,A,I3,A,I3,A,F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to alpha HOLE    :",hole(frg)
            print'(X,A,I3,A,I3,A,I3,A,F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to alpha PARTICLE:",particle(frg)
        enddo
    endif

    
    !If no beta requested, stop here!
    if ( Nselect_beta == 0 ) stop
!    print*, ""
    
    !Beta
    !=====
    call read_MO(IFCHK,MO,N,"beta",-1)
    !NTO averaging: get NTO coefficients (W) (stored in MO energies section)
    ! note that afeter this call M=N
    if (do_nto) then
        call read_NrEl(IFCHK,Nel,"beta")
        call read_MOenergy(IFCHK,W,N,"beta",-1)
        !Normalize
        Norm = 0.0
        do i=1,Nel
            Norm = Norm + W(i)
        enddo
        do i=1,Nel
            W2(i) = W2(i)/Norm
            !The particle is equal
            j = 2*Nel + 1 - i
            W2(j) = W2(j)/Norm
        enddo
    endif


    do frg = 1, Nfrg
        hole(frg) = 0.0
        particle(frg) = 0.0

        do ii=1,Nselect_beta
            i = MO_select_beta(ii)
            ! Compute contributions
            coeff_total = 0.0
            coeff_metal = 0.0
            do j=1,M
                coeff_total = coeff_total + abs(MO(j,i))
                !Contribution of fragments
                do k=1,Nat(frg)
                    if ( AOmap(j) == atom_index(frg,k) ) &
                     coeff_metal = coeff_metal + abs(MO(j,i))
                enddo
            enddo
            if (.not.do_nto) &
            print'(X,A,I3,A,I3,A,I3,A,I4,A, F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to MO beta ", i,":",coeff_metal/coeff_total * 100
            !NTO av: store in array
            if (do_nto) then
                if (i == Nel) then
                    hole(frg)     = hole(frg) + coeff_metal/coeff_total * 100.0 !* W(i)
                else if (i == Nel+1) then
                    particle(frg) = particle(frg) + coeff_metal/coeff_total * 100.0 !* W(i)
                endif
            endif

        enddo

    enddo

    !NTO av:
    if (do_nto) then
        print*, "Sum of Beta Coefficients: ", Norm
        print*, "First NTO coef           : ", W(Nel)!, W(Nel+1) 
        do frg = 1, Nfrg
            print'(X,A,I3,A,I3,A,I3,A,F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to beta  HOLE    :",hole(frg)
            print'(X,A,I3,A,I3,A,I3,A,F6.1)', "Fragment ",frg," with atoms ",atom_index(frg,1),'-',atom_index(frg,Nat(frg)),&
                                                    " contribution (%) to beta  PARTICLE:",particle(frg)
        enddo
    endif
     
    stop

    contains

    subroutine parse_input(inpfile,mo_file,at_file,Metal_Index,do_nto)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile, mo_file, at_file
        integer,intent(inout) :: Metal_Index
        logical,intent(inout) :: do_nto
        
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg, dummy_char

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

                case ("-mo")
                    call getarg(i+1, mo_file)
                    argument_retrieved=.true.

                case ("-frg")
                    call getarg(i+1, at_file)
                    argument_retrieved=.true.

                case ("-ndx")
                    call getarg(i+1, dummy_char)
                    read(dummy_char,*) Metal_Index
                    argument_retrieved=.true.

                case ("-do_nto")
!                     nto_av = .true.
                    do_nto = .true.

                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo


       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '              c o n t r i b  M O '
        write(0,'(/,A)') '     Compute the contribution of an atom '
        write(0,'(A)')   '            to a given set of MOs'
        write(0,'(/,A)') '                    V4'
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-mo             ', trim(adjustl(mo_file))
        write(0,*) '-frg            ', trim(adjustl(at_file))
        write(dummy_char,*) Metal_Index
        write(0,*) '-ndx            ', trim(adjustl(dummy_char))
        write(0,*) '-do_nto        ',  do_nto
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


end program contribMO

    

