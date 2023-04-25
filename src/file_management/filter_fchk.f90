program filter_fchk


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to select a subset of sections of the FCHK file
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
    !============================================
    !   File readers
    !============================================
    use gaussian_manage


    implicit none

    integer,parameter :: NDIM = 600

    !====================== 
    !Options 
    logical :: overwrite=.false.
    !======================

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=42) :: sect_label
    character(len=1) :: dtype
    integer :: error, N, lenght
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=50) :: dummy_char
    character(len=240):: line
    integer :: iat_new, iat_orig, j_orig, j_new, nswap
    !====================== 

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk, iat, nn, imin, imax, iii
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_SEC=11,  &
               O_FCHK=20
    !files
    character(len=10) :: filetype="fchk"
    character(len=200):: inpfile ="input.fchk", &
                         sectionfile='sections.dat', &
                         outfile="default"
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
    call parse_input(inpfile,outfile,overwrite,sectionfile)
 
    !================
    ! OPEN FILES
    !================
    ! INPUT FILE
    print'(X,A)', "  Input file: "//trim(adjustl(inpfile))
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )
    if (adjustl(filetype) == "guess") &
     call split_line_back(inpfile,".",null,filetype)
    ! OUTPUT FILE
    print'(X,A)', "  Output file: "//trim(adjustl(outfile))
    if (overwrite) then
        open(O_FCHK,file=outfile,status="replace",iostat=IOstatus)
    else
        open(O_FCHK,file=outfile,status="new",iostat=IOstatus)
        if (IOstatus /= 0) &
         call alert_msg("fatal","Cannot open output for writting. Use -ow to force overwrite")
    endif
    ! SECTION FILE
    print'(X,A)', "  Sections file: "//trim(adjustl(sectionfile))
    open(I_SEC,file=sectionfile,status='old',iostat=IOstatus)

    !================
    ! READ/WRITE FILES
    !================
    ! TITLE and JOB INFO
    ! Read 2 first lines, always
    read(I_INP,'(A)') line
    write(O_FCHK,'(A)') adjustl(trim(line))
    read(I_INP,'(A)') line
    write(O_FCHK,'(A)') adjustl(trim(line))
    
    ! Select sections
    do while ( len_trim(line) /= 0 )
        read(I_SEC,'(A42,X,A1)',iostat=IOstatus) sect_label, dtype
        if (IOstatus /= 0) exit
        print*, "Copying section: "//adjustl(trim(sect_label))
        ! Read
        call read_fchk(I_INP,adjustl(trim(sect_label)),dtype,N,A,IA,error,rwnd=.false.)
        if (error /= 0) then
            call alert_msg( "warning","Unable to read section "//trim(adjustl(sect_label)) )
            cycle
        endif
        ! Check if supported dtype
        if (dtype == 'C') then
            call alert_msg( "warning","Section type not yet supported: "//dtype )
            cycle
        endif
        ! Write
        call write_fchk(O_FCHK,adjustl(trim(sect_label)),dtype,N,A,IA,error)
        ! After using A or IA, we need to deallocate for next use
        if (allocated(A))  deallocate(A)
        if (allocated(IA)) deallocate(IA)
    enddo

    print'(X,A,/)', "Done"

    call cpu_time(tf)
    write(0,'(/,A,X,F12.3,/)') "CPU time (s)", tf-ti

    close(I_INP)
    close(I_SEC)
    close(O_FCHK)

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,outfile,overwrite,sectionfile)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,outfile,sectionfile
        logical,intent(inout)          :: overwrite
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
                case("-o")
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                case("-sel")
                    call getarg(i+1, sectionfile)
                    argument_retrieved=.true.
                case("-ow")
                    overwrite=.true.
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        !The default output is now xyz
        if (adjustl(outfile) == "default") then
            call split_line_back(inpfile,".",outfile,null)
            outfile=trim(adjustl(outfile))//"_F.fchk"
        endif

       !Print options (to stdx)
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '             F I L T E R   F C H K '    
        write(0,'(/,A)') '      Filter fchk picking a subset of sections'     
        call print_version()
        write(0,'(/,A)') '========================================================'
        write(0,'(/,A)') '-------------------------------------------------------------------'
        write(0,'(A)')   ' Flag         Description                      Value'
        write(0,'(A)')   '-------------------------------------------------------------------'
        write(0,*)       '-f           Input file                       ', trim(adjustl(inpfile))
        write(0,*)       '-o           Output file (fchk)               ', trim(adjustl(outfile))
        write(0,*)       '-sel         File with selected sections      ', trim(adjustl(sectionfile))
        write(0,*)       '-ow          Force overwrite output          ',  overwrite
        write(0,*)       '-h           This help                       ',  need_help
        write(0,*)       '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input
       

end program filter_fchk

