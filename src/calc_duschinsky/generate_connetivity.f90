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

    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: use_symmetry=.false. ,&
               modred=.false.       ,&
               tswitch=.false.      ,&
               symaddapt=.false.    ,&
               use_hbonds=.false.   ,&
               overwrite=.false.
    character(len=4) :: def_internal='zmat'
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: state1,state2
    integer,dimension(1:NDIM) :: isym
    integer :: Nat, Nvib, Ns, N
    !====================== 

    !====================== 
    !AUXILIAR MATRICES
    real(8),dimension(NDIM,NDIM) :: Aux, Aux2
    !VECTORS
    real(8),dimension(NDIM) :: S1, S2, Vec, Vec2, Factor
    integer,dimension(NDIM) :: S_sym, bond_sym,angle_sym,dihed_sym, OrderMap
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
               i1,i2,i3,i4, i_old, i_new
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_ZMAT=11, &
               I_SYM=12,  &
               I_RMF=16,  &
               I_CNX=17,  &
               I_ORD=18,  &
               O_ICF=20
    !files
    character(len=10) :: ft ="guess", stat
    character(len=200):: inpfile  ="input.fchk",     &
                         intfile  ="int_coords.dat", &
                         rmzfile  ="none",           &
                         zmatfile ="none",           &
                         symm_file="none",           &
                         cnx_file="connect.dat",     &
                         order_file="none"
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
    call parse_input(inpfile,ft,zmatfile,rmzfile,def_internal,use_symmetry,intfile,cnx_file,order_file,use_hbonds,&
                     overwrite)
    call set_word_upper_case(def_internal)

    ! READ DATA (each element from a different file is possible)
    ! ---------------------------------
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)

    if (adjustl(ft)/="gview") then
        ! STRUCTURE FILE
        open(I_INP,file=inpfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
        call generic_strmol_reader(I_INP,ft,state1)
        close(I_INP)
        ! Shortcuts
        Nat = state1%natoms
        
        ! MANAGE INTERNAL COORDS
        ! ---------------------------------
        ! Get connectivity 
        call guess_connect(state1,use_hbonds)
    else
        ! GVIEW CONNECTIVITY
        open(I_INP,file=inpfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
        call read_connect(I_INP,state1,ft)
    endif

    ! OUTPUT CONNECTIVUTY FILE
    if (overwrite) then
        stat="unknown"
    else
        stat="new"
    endif
    open(I_CNX,file=cnx_file,status=stat,iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(cnx_file))//" for writting. Already exists?" )
    call write_connect(I_CNX,state1)
    close(I_CNX)

    call cpu_time(tf)
    write(6,'(A,F12.3)') "CPU (s) for internal vib analysis: ", tf-ti

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft,zmatfile,rmzfile,def_internal,use_symmetry,intfile,cnx_file,order_file,use_hbonds,&
                           overwrite)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft,zmatfile,&
                                          intfile,rmzfile,def_internal,cnx_file, &
                                          order_file
        logical,intent(inout)          :: use_symmetry, use_hbonds, overwrite
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

                case ("-cnx") 
                    call getarg(i+1, cnx_file)
                    argument_retrieved=.true.

                case ("-reorder")
                    call getarg(i+1, order_file)
                    argument_retrieved=.true.

                case ("-intfile") 
                    call getarg(i+1, intfile)
                    argument_retrieved=.true.

                case ("-zmatfile") 
                    call getarg(i+1, zmatfile)
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

                case ("-hbonds")
                    use_hbonds=.true.
                case ("-nohbonds")
                    use_hbonds=.false.

                case ("-ow")
                    overwrite=.true.
                case ("-noow")
                    overwrite=.false.
        
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
        write(6,'(/,A)') '                   ZMAT2INTERNALSET '    
        write(6,'(/,A)') '  Get internal set from the automatically generated Zmat  '
        write(6,'(/,A)') '           '        
        write(6,'(/,A)') '--------------------------------------------------'
        write(6,*) '-f              ', trim(adjustl(inpfile))
        write(6,*) '-ft             ', trim(adjustl(ft))
        write(6,*) '-cnx            ', trim(adjustl(cnx_file))
        write(6,*) '-reorder        ', trim(adjustl(order_file))
        write(6,*) '-intmode        ', trim(adjustl(def_internal))
        write(6,*) '-intfile        ', trim(adjustl(intfile))
        write(6,*) '-zmatfile       ', trim(adjustl(zmatfile))
        write(6,*) '-rmzfile        ', trim(adjustl(rmzfile))
        write(6,*) '-[no]hbonds    ',  use_hbonds
        write(6,*) '-[no]ow        ',  overwrite
        write(6,*) '-h             ',  need_help
        write(6,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program vertical2adiabatic

