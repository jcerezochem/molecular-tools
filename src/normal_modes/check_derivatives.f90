program check_derivatives


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
    logical :: use_symmetry=.false. ,&
               modred=.false.       ,&
               tswitch=.false.      ,&
               symaddapt=.false.    ,&
               vertical=.true.      ,&
               do_correct_vert=.false.
    character(len=4) :: def_internal='zmat'
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molecule
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns, Ni, Nj, Nk
    !====================== 

    !====================== 
    !INTERNAL VIBRATIONAL ANALYSIS
    !MATRICES
    !B and G matrices
    real(8),dimension(NDIM,NDIM) :: B1,B2, B, G1,G2
    !Other arrays
    real(8),dimension(1:NDIM) :: Grad
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, X1,X1inv,X2,X2inv, L1,L2, Asel1, Asel2, Asel
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM) :: Lder
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
    real(8),dimension(:),allocatable :: Hlt
    integer :: error
    real(8) :: scl
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=16) :: dummy_char
    real(8) :: Theta, Theta2, Theta3
    ! Messages
    character(len=200) :: msg
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
               I_DER=17,  &
               O_DUS=20,  &
               O_DIS=21,  &
               O_DMAT=22, &
               O_DUS2=23, &
               O_DIS2=24, &
               O_STAT=25, &
               O_STR =26
    !files
    character(len=10) :: ft ="guess", fth="guess", ftgv="guess", fthv="guess" 
    character(len=200):: inpfile  ="input.fchk", &
                         hessfile ="same", &
                         hessfile_v ="vertical.fchk", &
                         gradfile_v ="same", &
                         intfile  ="none",       &
                         rmzfile  ="none",       &
                         symm_file="none",     & 
                         derfile="base", derfile_base, &
                         tmpfile
    character(len=10) :: dertype="Q"
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
!     call generic_input_parser(inpfile, "-f" ,"c",&
!                               filetype,"-ft","c",&
!                               )
    call parse_input(inpfile,ft,derfile,dertype)
    call set_word_upper_case(dertype)


    ! 1. READ DATA
    ! ---------------------------------
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)
        
    ! STRUCTURE FILE
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,molecule,error)
    if (error /= 0) call alert_msg("fatal","Error reading geometry (State1)")
    Nat = molecule%natoms
    Nvib = 3*Nat-6


    if (dertype == "Q") then
        Ni = 3*Nat
        Nj = Nvib
        Nk = Nvib
        scl = 1.d2
    elseif (dertype == "X") then
        Ni = Nvib
        Nj = 3*Nat
        Nk = 3*Nat
        scl = 1.d0
    elseif (dertype == "QINT") then
        Ni = Nvib
        Nj = Nvib
        Nk = Nvib
        scl = 1.d0
    elseif (dertype == "S") then
        Ni = Nvib
        Nj = Nvib
        Nk = Nvib
        scl = 1.d0
    else

        call alert_msg("fatal","Unkown derivative type")
    endif

    ! Fill Lder tensor
    derfile_base=derfile
    do j=1,Nj
        write(derfile,'(A,I0,A)') trim(adjustl(derfile_base)), j, ".dat"
        print*, "Reading ", trim(adjustl(derfile))
        open(I_DER,file=derfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(derfile)) )
        do i=1,Ni
            read(I_DER,*) Lder(i,j,1:Nk)
        enddo
        close(I_DER)
    enddo

    if (verbose>0) then
        do i=1,Ni
            write(tmpfile,'(A,I0,A)') "Lder *scl, Cart=",i
            call MAT0(6,Lder(i,:,:)*scl,Nj,Nk,trim(tmpfile))
        enddo
    endif

!         do j=1,3*Nat
!             write(derfile,'(A,I0,A)') trim(adjustl(derfile_base)), j, ".dat"
!             open(I_DER,file=derfile,status='old',iostat=IOstatus)
!             if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(derfile)) )
!             do i=1,Nvib
!                 read(I_DER,*) Lder(i,j,1:Nvib)
!             enddo
!             close(I_DER)
!         enddo
! 
!         if (verbose>2) then
!             do i=1,Nvib
!                 write(tmpfile,'(A,I0,A)') "Lder, Q=",i
!                 call MAT0(6,Lder(i,:,:),3*Nat,3*Nat,trim(tmpfile))
!             enddo
!         endif
    
    call summary_alerts

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,derfile,dertype)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile, ft, derfile, dertype
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
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

                case ("-fder") 
                    call getarg(i+1, derfile)
                    argument_retrieved=.true.

                case ("-type") 
                    call getarg(i+1, dertype)
                    argument_retrieved=.true.
        
                case ("-h")
                    need_help=.true.

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



       !Print options (to stderr)
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,'(/,A)') '        Check Numerical Drivatives'    
        write(6,'(/,A)') '           '
        write(6,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(ft))
        write(6,*) '-fder           ', trim(adjustl(derfile))
        write(6,*) '-type           ', trim(adjustl(dertype))
        write(6,*) '-h             ',  need_help
        write(6,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program check_derivatives

