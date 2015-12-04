program vertical2adiabatic


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
    type(str_resmol) :: state1,state2
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns
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
    call parse_input(inpfile,ft,hessfile,fth,gradfile_v,ftgv,hessfile_v,fthv,&
                     intfile,rmzfile,def_internal,use_symmetry,derfile,do_correct_vert)
    call set_word_upper_case(def_internal)

    ! READ DATA (each element from a different file is possible)
    ! ---------------------------------
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)
    if (fth == "guess") &
    call split_line_back(hessfile,".",null,fth)
    if (ftgv == "guess") &
    call split_line_back(gradfile_v,".",null,ftgv)
    if (fthv == "guess") &
    call split_line_back(hessfile_v,".",null,fthv)

    ! STRUCTURE FILE
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,state1)
    close(I_INP)
    ! Shortcuts
    Nat = state1%natoms

    ! HESSIAN FILE (State1)
    open(I_INP,file=hessfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile)) )
    allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fth,Nat,Hlt,error) 
    close(I_INP)
    ! Run vibrations_Cart to get the number of Nvib (to detect linear molecules)
    call vibrations_Cart(Nat,state1%atom(:)%X,state1%atom(:)%Y,state1%atom(:)%Z,state1%atom(:)%Mass,Hlt,&
                         Nvib,L1,Freq,error)
    deallocate(Hlt)
    ! Get Lcart = M^1/2 Lmwc 
    call Lmwc_to_Lcart(Nat,Nvib,state1%atom(:)%Mass,L1,L1,error)
    if (verbose>2) &
        call MAT0(6,L1,3*Nat,Nvib,"Lcart matrix")

    ! HESSIAN FILE (State2)
    open(I_INP,file=hessfile_v,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(hessfile_v)) )
    allocate(Hlt(1:3*Nat*(3*Nat+1)/2))
    call generic_Hessian_reader(I_INP,fthv,Nat,Hlt,error) 
    close(I_INP)
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
    open(I_INP,file=gradfile_v,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(gradfile_v)) )
    call generic_gradient_reader(I_INP,ftgv,Nat,Grad,error)
    close(I_INP)

    !Compute H_Q = L1^t Hess L1  +  gx LLL
    Hess(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,3*Nat,L1,Hess,counter=.true.)

    ! Apply matrix derivative if the option is enabled
    if (do_correct_vert) then
        ! Fill Lder tensor
        derfile_base=derfile
        do j=1,Nvib
            write(derfile,'(A,I0,A)') trim(adjustl(derfile_base)), j, ".dat"
            open(I_DER,file=derfile,status='old',iostat=IOstatus)
            if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(derfile)) )
            do i=1,3*Nat
                read(I_DER,*) Lder(i,j,1:Nvib)
            enddo
            close(I_DER)
        enddo

        if (verbose>2) then
            do i=1,3*Nat
                write(tmpfile,'(A,I0,A)') "Lder *10^6, Cart=",i
                call MAT0(6,Lder(i,:,:)*1.e6,Nvib,Nvib,trim(tmpfile))
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
    
        do i=1,Nvib
        do j=1,Nvib
            Aux(i,j) = 0.d0
            do l=1,3*Nat
                Aux(i,j) = Aux(i,j) + Grad(l) * Lder(l,j,i) 
            enddo
            Hess(i,j) = Hess(i,j) + Aux(i,j)
        enddo
        enddo
    endif


    call diagonalize_full(Hess(1:Nvib,1:Nvib),Nvib,L2(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
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

    call summary_alerts

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,hessfile,fth,gradfile_v,ftgv,hessfile_v,fthv,intfile,&
                           rmzfile,def_internal,use_symmetry,derfile,do_correct_vert)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,hessfile,fth,gradfile_v,ftgv,hessfile_v,fthv,&
                                          intfile,rmzfile,def_internal,derfile
        logical,intent(inout)          :: use_symmetry,do_correct_vert
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

                case ("-fder") 
                    call getarg(i+1, derfile)
                    argument_retrieved=.true.

                case ("-fhessv") 
                    call getarg(i+1, hessfile_v)
                    argument_retrieved=.true.
                case ("-fthv") 
                    call getarg(i+1, fthv)
                    argument_retrieved=.true.

                case ("-fgradv") 
                    call getarg(i+1, gradfile_v)
                    argument_retrieved=.true.
                case ("-ftgv") 
                    call getarg(i+1, ftgv)
                    argument_retrieved=.true.

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-rmzfile") 
                    call getarg(i+1, rmzfile)
                    argument_retrieved=.true.
                ! Kept for backward compatibility (but replaced by -rmzfile)
                case ("-rmz") 
                    call getarg(i+1, rmzfile)
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

                case ("-correct")
                    do_correct_vert=.true.
                case ("-nocorrect")
                    do_correct_vert=.false.
        
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

       ! Manage defaults
       ! If not declared, hessfile and gradfile are the same as inpfile
       if (adjustl(hessfile) == "same") then
           hessfile=inpfile
           if (adjustl(fth) == "guess")  fth=ft
       endif
       if (adjustl(gradfile_v) == "same") then
           gradfile_v=hessfile_v
           if (adjustl(ftgv) == "guess")  ftgv=fthv
       endif


       !Print options (to stderr)
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,'(/,A)') '        V E R T I C A L  --PCCP-- '    
        write(6,'(/,A)') '  Displace structure from vertical to adiabatic geoms  '
        write(6,'(/,A)') '           '        
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,*) '-f              ', trim(adjustl(inpfile))
        write(6,*) '-ft             ', trim(adjustl(ft))
        write(6,*) '-fhess          ', trim(adjustl(hessfile))
        write(6,*) '-fth            ', trim(adjustl(fth))
        write(6,*) '-fhessv         ', trim(adjustl(hessfile_v))
        write(6,*) '-fthv           ', trim(adjustl(fthv))
        write(6,*) '-fgradv         ', trim(adjustl(gradfile_v))
        write(6,*) '-ftgv           ', trim(adjustl(ftgv))
        write(6,*) '-fder           ', trim(adjustl(derfile))
!         write(6,*) '-intmode        ', trim(adjustl(def_internal))
!         write(6,*) '-intfile        ', trim(adjustl(intfile))
!         write(6,*) '-rmzfile        ', trim(adjustl(rmzfile))
        write(6,*) '-[no]correct   ', do_correct_vert
        write(6,*) '-h             ',  need_help
        write(6,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program vertical2adiabatic

