program state2hessian


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
    use xyz_manage
    !============================================
    !  Vibrationsl
    !============================================
    use vibrational_analysis

    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !System variables
    integer :: Nat, Nvib, N, AtNum
    character(len=100) :: title
    !MATRICES
    real(8),dimension(1:NDIM,1:NDIM) :: Hess, L, D
    !VECTORS
    real(8),dimension(NDIM) :: Freq, FC, X, Y, Z, Mass
    character(len=2),dimension(NDIM) :: AtName
    integer,dimension(NDIM) :: IA
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: Hlt
    integer :: error
    !====================== 

    !=============
    !Counters
    integer :: i,j,k, ii
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               O_FCHK=20
    !files 
    character(len=200):: inpfile  ="fcc.inp", &
                         fccfile  ="state_file", &
                         outfile  ="output.fchk"
    !status
    integer :: IOstatus
    integer :: current_verbose
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
!     alert_unt=6
    !--------------------------

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,fccfile,outfile)


    ! READ INPUT DATA 
    ! Input file
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    read(I_INP,*) Nat
    read(I_INP,*) Nvib
    do i=1,Nat 
        read(I_INP,*) Mass(i)
    enddo    
    close(I_INP)

    ! State file
    open(I_INP,file=fccfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(fccfile)) )

    ! Read geom
    do i=1,Nat 
       read(I_INP,*) X(i)
       read(I_INP,*) Y(i)
       read(I_INP,*) Z(i)
    enddo
    ! Read Lcart
    do j=1,3*Nat
    do i=1,Nvib
        read(I_INP,*) L(j,i)
    enddo
    enddo
    ! Read Freq
    do i=1,Nvib
        read(I_INP,*) Freq(i)
    enddo
    close(I_INP)
    if (verbose>1) &
     call print_vector(6,Freq,Nvib,"Frequencies -- Read")

    ! ROTATE DIAG FC MATRIX TO CARTESIAN HESSIAN
    ! Transform Lcart(normalized) to Lmwc 
    call LcartNrm_to_Lmwc(Nat,Nvib,Mass,L,L)
    ! Transform Freq to FC (diagonal matrix)
    FC(1:Nvib) = Freq2FC(Nvib,Freq)
    Hess(1:Nvib,1:Nvib) = 0.d0
    do i=1,Nvib
        Hess(i,i) = FC(i)
    enddo

    ! Hx = M^1/2 Lmwc FC Lmwc^t M^1/2
    ! First M^1/2 Lmwc:
    do j=1,Nvib 
    do i=1,3*Nat
        ii = (i-1)/3+1
        L(i,j) = L(i,j)*dsqrt(Mass(ii)*AMUtoAU)
    enddo
    enddo
    ! Then rotate FC to Hx
    Hess(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,Nvib,L,Hess,counter=.false.)

    ! Check the Hessian
    N=3*Nat*(3*Nat+1)/2
    allocate(Hlt(1:N))
    Hlt(1:N) = Hess_to_Hlt(3*Nat,Hess)
    if (verbose>1) then 
        print'(2/,X,A)', "Recalculating frequencies from new Hessian"
        print*,         "============================================"
        current_verbose=verbose
        verbose=1
        call vibrations_Cart(Nat,X,Y,Z,Mass,Hlt,Nvib,L,Freq)
        verbose=current_verbose
    endif


    ! Prepare to write the outfile 
    do i=1,Nat
        call atominfo_from_atmass(Mass(i),AtNum,AtName(i))
    enddo
    title="FCHK from FCclasses state file: "//trim(adjustl(fccfile)) 

    ! Write Geom Hessian in FCHK format
    open(O_FCHK,file=outfile)
    ! Geom
    call generic_structure_writer(O_FCHK,'fchk',Nat,X,Y,Z,Mass,AtName,title=trim(adjustl(title)))
    ! Hessian
    call write_fchk(O_FCHK,"Cartesian Force Constants",'R',N,Hlt,IA,error)
    close(O_FCHK)
    deallocate(Hlt)

    print'(/,X,A)', "File generated: "//trim(adjustl(outfile))    

    call summary_alerts

    call cpu_time(tf)
    write(6,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,fccfile,outfile)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,fccfile,outfile
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        character(len=500) :: input_command
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
                case ("-input") 
                    call getarg(i+1, inpfile)
                    argument_retrieved=.true.

                case ("-state") 
                    call getarg(i+1, fccfile)
                    argument_retrieved=.true.

                case ("-o") 
                    call getarg(i+1, outfile)
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

       !Print options (to stdout)    
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '           S T A T E  2  H E S S I A N            '
        write(6,'(/,A)') '   Transform FCclasses state file to Hessian (FCHK)    '      
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag         Description                   Value'
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,*)       '-input       FCclasses input file          ', trim(adjustl(inpfile))
        write(6,*)       '-state       FClasses state file           ', trim(adjustl(fccfile))
        write(6,*)       '-o           Output file (fchk)            ', trim(adjustl(outfile))
        write(6,*)       '-h               ',  need_help
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                       'Verbose level:  ', verbose        
        write(6,'(A)') '-------------------------------------------------------------------'     
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program state2hessian

