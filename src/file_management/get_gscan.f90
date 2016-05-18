program get_gscan


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Extract the steps from a scan to individual fchk
    !
    ! Change log:
    !
    ! TODO:
    ! ------
    !
    !============================================================================    

    use alerts
    use structure_types
    use line_preprocess
!     use gro_manage
!     use pdb_manage
    use gaussian_manage
    use gaussian_fchk_manage
    use xyz_manage
!     use molcas_unsym_manage
    use ff_build
    use molecular_structure
    use atomic_geom
    use symmetry_mod

    implicit none

    !====================== 
    !System variables
    type(str_resmol) :: molec, molec_aux
    real(8) :: E, geompar
    character(len=100) :: swapfile="NO"
    integer,dimension(1:1000) :: isym
    !====================== 

    !=============
    !Counters and dummies
    integer :: i,j,k,l, ii,jj,kk, iat, &
               i1,i2,i3,i4
    character(len=1) :: ctype, null
    character(len=20) :: dummy_char
    integer :: ilength
    character(len=180) :: line
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_SWP=11,  &
               O_STR=20,  &
               O_FCHK=21,  &
               O_ENR=22
    !files
    character(len=10) :: filetype_inp="guess",&
                         filetype_out="guess"
    character(len=200):: inpfile="input.log",          &
                         outfile="default"
    !status
    integer :: IOstatus
    !===================


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype_inp,outfile,filetype_out,swapfile)

 
    ! 1. READ INPUT
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

 
    !open energy file
    call split_line_back(inpfile,".",outfile,null)
    outfile=trim(adjustl(outfile))//"_scan.dat"
    open(O_ENR,file=outfile) 

    !Read logfile till ModRedundant section
    do 
       read(I_INP,'(X,A)',IOSTAT=IOstatus) line
        ! 1) End of file
        if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file")
        ! 2) ModRedundant section
        if ( INDEX(line,"The following ModRedundant input section has been read:") /= 0 ) then
            read(I_INP,'(X,A)',IOSTAT=IOstatus) line
            call split_line_back(line,"S",line,ctype)
            read(line,*) ctype
            if (ctype == "D") then
                call split_line_back(line,"D",ctype,line)
                read(line,*) i1, i2, i3, i4
            elseif (ctype == "B") then
                call split_line_back(line,"B",ctype,line)
                read(line,*) i1, i2
            endif
            exit
        endif
    enddo

    !Now continue a read every step in the scan
    i=0
    do 
       read(I_INP,'(X,A)',IOSTAT=IOstatus) line
        ! Three possible scenarios while reading:
        ! 1) End of file
        if ( IOstatus < 0 ) exit !call alert_msg("fatal","Unexpected end of file")
        ! 2) Optimized step
        if ( INDEX(line,"-- Stationary point found.") /= 0 ) then
            i=i+1
            !xyz structure
            ! Get file name from inpfile
            call split_line_back(inpfile,".",outfile,null)
            outfile=trim(adjustl(outfile))//"_step"//int20char(i,2)//".xyz"
            open(O_STR,file=outfile)           
            ! Update molec title with step and energy
            molec%title="Step "//int2char(i,2)//"  E = "//real2char(E,17,12)
            call write_xyz(O_STR,molec)
            close(O_STR)
            !FCHK file
            ! Get file name from inpfile
            call split_line_back(inpfile,".",outfile,null)
            outfile=trim(adjustl(outfile))//"_step"//int20char(i,2)//".fchk"
            open(O_FCHK,file=outfile)           
            call write_fchk_E(O_FCHK,molec,E)
            close(O_FCHK)
            !Get parameter
            if (ctype == "D") then
                geompar = calc_dihed(molec%atom(i1),&
                                     molec%atom(i2),&
                                     molec%atom(i3),&
                                     molec%atom(i4))
                geompar = geompar*180.d0/PI
            elseif (ctype == "B") then
                geompar = calc_dist(molec%atom(i1),&
                                    molec%atom(i2))
            endif
            write(O_ENR,*) geompar, E
        endif
        ! 3) Read structure and associted energy
        if ( INDEX(line,    &
             "orientation") &
              /= 0 ) then
            backspace(I_INP)
            call get_next_std_geom(I_INP,molec)
            if (swapfile/="NO") then
                molec_aux = molec
                open(I_SWP,file=swapfile)
                do ii=1,molec%natoms
                    read(I_SWP,*) jj,j
                    molec%atom(jj) = molec_aux%atom(j)
                enddo
                close(I_SWP)
            endif
        endif
        if ( INDEX(line,"SCF Done") /= 0 ) then
          call split_line_back(line,"=",null,line)
          read(line,*) E
        endif
    enddo
    close(O_ENR)

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype_inp,outfile,filetype_out,swapfile)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,outfile,&
                                          filetype_inp,filetype_out, &
                                          swapfile
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
                case ("-fti") 
                    call getarg(i+1, filetype_inp)
                    argument_retrieved=.true.

                case ("-o") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                case ("-fto") 
                    call getarg(i+1, filetype_out)
                    argument_retrieved=.true.

                case ("-swap")
                    call getarg(i+1, swapfile)
                    argument_retrieved=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        if (adjustl(outfile) == "default") then
            call split_line_back(inpfile,".",outfile,null)
            outfile=trim(adjustl(outfile))//".gro"
        endif

        ! Some checks on the input
        !----------------------------

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '        S C A N   R E A D E R '    
        write(0,'(/,A)') '        Get opt structures from a scan '  
        write(0,'(/,A)') '        Revision: get_gscan-20141022               '         
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-fti            ', trim(adjustl(filetype_inp))
        write(0,*) '-swap           ', trim(adjustl(swapfile))
        write(0,*) '-o              ', trim(adjustl(outfile))
        write(0,*) '-fto            ', trim(adjustl(filetype_out))
!         write(0,*) '-r             ',  overwrite
!         write(0,*) '-rn             ',  trim(adjustl(resname))
!         write(0,*) '-connect       ',  make_connect
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input



end program get_gscan

