program normal_modes_cartesian


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS 
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Generate gaussian jobs to perform derivatives in terms
    ! with respect to Cartesian coordinates. Applied to freqs
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
    use xyz_manage
    use gro_manage
    use gaussian_manage

    implicit none

    integer,parameter :: NDIM=600

    type(str_resmol) :: molec
    
    real(8),dimension(:),allocatable          :: X0,Y0,Z0
    character(len=2),dimension(:),allocatable :: AtName
    integer,dimension(:),allocatable          :: AtNum
    integer                                   :: Nat
    
    real(8) :: delta=1.d-3, disp

    integer            :: I_INP=10, &
                          O_G09=20
    integer            :: IOstatus, error
    character(len=200) :: inpfile="strucutre.log", &
                          g09file, chkfile
    character(len=50)  :: ft="guess"
    character          :: null

    integer :: i,j,k

    real(8) :: ti, tf

    call cpu_time(ti)

    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,ft)

 
    ! 1. READ DATA
    ! ---------------------------------
    !Guess filetypes
    if (ft == "guess") &
    call split_line_back(inpfile,".",null,ft)
        
    ! STRUCTURE FILE
    open(I_INP,file=inpfile,iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    call generic_strmol_reader(I_INP,ft,molec,error)
    close(I_INP)
    Nat = molec%natoms

    allocate( X0(1:Nat),Y0(1:Nat),Z0(1:Nat) )
    X0 = molec%atom(1:Nat)%X
    Y0 = molec%atom(1:Nat)%Y
    Z0 = molec%atom(1:Nat)%Z

    ! GENERATE DISPLACEMENTS
    k=0
    do i=1,Nat
        ! X
        k=k+1
        write(g09file,'(A,I0,A)') "Cart_coord_", k, ".com" 
        open(O_G09,file=g09file,iostat=IOstatus)
        ! Scan
        disp=-2.d0*delta
        molec%atom(i)%x = X0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=-delta
        molec%atom(i)%x = X0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=0
        molec%atom(i)%x = X0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=delta
        molec%atom(i)%x = X0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=2.d0*delta
        molec%atom(i)%x = X0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        ! Close file and restore input geometry
        close(O_G09)
        molec%atom(i)%x = X0(i)

        ! Y
        k=k+1
        write(g09file,'(A,I0,A)') "Cart_coord_", k, ".com"
        open(O_G09,file=g09file,iostat=IOstatus)
        ! Scan
        disp=-2.d0*delta
        molec%atom(i)%y = Y0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=-delta
        molec%atom(i)%y = Y0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=0
        molec%atom(i)%y = Y0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=delta
        molec%atom(i)%y = Y0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=2.d0*delta
        molec%atom(i)%y = Y0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        ! Close file and restore input geometry
        close(O_G09)
        molec%atom(i)%y = Y0(i)

        ! Z
        k=k+1
        write(g09file,'(A,I0,A)') "Cart_coord_", k, ".com"
        open(O_G09,file=g09file,iostat=IOstatus)
        ! Scan
        disp=-2.d0*delta
        molec%atom(i)%z = Z0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=-delta
        molec%atom(i)%z = Z0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=0
        molec%atom(i)%z = Z0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=delta
        molec%atom(i)%z = Z0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        disp=2.d0*delta
        molec%atom(i)%z = Z0(i) + disp
        write(molec%title,'(A,I0,A,2(X,F12.6))') " Coord ",k," Disp = ", disp/BOHRtoANGS
        call write_gcom(O_G09,molec,g09file,title=molec%title)
        ! Close file and restore input geometry
        close(O_G09)
        molec%atom(i)%z = Z0(i)

    enddo


    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,ft)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,ft

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

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 


       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '          CARTESIAN MODES ANIMATION '    
        write(0,'(/,A)') '         Perform vibrational analysis'
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(ft))
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input

end program normal_modes_cartesian

