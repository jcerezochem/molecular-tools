program fchk2gro


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Forms a gro structure with the one present in the fchk
    !
    ! Compilation instructions (for mymake script):
    !make$ echo "COMPILER: $FC"; sleep 1; $FC ../modules/alerts.f90 ../modules/structure_types_v2.f90 ../modules/line_preprocess.f90  ../modules/gro_manage_v2.f90 ../modules/pdb_manage_v2.f90 ../modules/constants_mod.f90 ../modules/gaussian_manage_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 ../modules/xyz_manage.f90 ../modules/ff_build_module_v3.f90 fchk2gro_v2.1.f90 -o fchk2gro_v2.1.exe -cpp 
    !
    ! Change log:
    !  July 2013: included xyz file format support
    ! -- VERSION 2 --
    !  July 2013: included std orientation from logs 
    !             included -connect to create connections (for PDB only) -- needs: ff_build_module_v3.f90
    ! -- VERSION 2.1 --
    !  July 2013: guess connect is placed out of generic_strfile_read/write
    !
    !  Version 4: addapted to v4 modules
    !  v4.1: add "swap" option (10/10/14)
    !        included fcc (structure part of state file)
    !  Jan15: Fixed the write_fchk call
    !  Jan15: Added -use-elems to force use element names and not FF names
    !  Jan15: Added -rmcog to remove the center of gravity (this is sueful 
    !         to use with avogadro, since with large distance, do not show labels)
    !
    ! TODO:
    ! ------
    !
    !============================================================================    

    use alerts
    use structure_types
    use line_preprocess
    use gro_manage
    use pdb_manage
    use gaussian_manage
    use gaussian_fchk_manage
    use xyz_manage
    use molden_manage
    use molcas_unsym_manage
    use psi4_manage
    use gamess_manage
    use ff_build
    use molecular_structure

    implicit none

    !====================== 
    !System variables
    type(str_resmol) :: molec, &
                        molec_aux
    !====================== 

    !=============
    !Counters and dummies
    integer :: i,j,k,l, jj,kk, iat
    character(len=1) :: null
    logical :: overwrite    = .false. ,&
               make_connect = .false. ,&
               use_elements = .false. ,&
               remove_com   = .false.
    !Swap related counters
    integer :: iat_new, iat_orig, nswap
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10, &
               I_SWP=11, &
               O_OUT=20  
    !files
    character(len=5) :: resname="read"
    character(len=10) :: filetype_inp="guess",&
                         filetype_out="guess"
    character(len=200):: inpfile="input.fchk",&
                         outfile="default"   ,&
                         swapfile="none"  
    !status
    integer :: IOstatus
    character(len=7) :: stat="new" !do not overwrite when writting
    !===================


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype_inp,outfile,filetype_out,overwrite,make_connect,&
                     use_elements,remove_com,resname,swapfile)

 
    ! 1. READ INPUT
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(filetype_inp) == "guess") call split_line_back(inpfile,".",null,filetype_inp)
    call generic_strfile_read(I_INP,filetype_inp,molec)
    close(I_INP)
    !Option to specify the resname from command line
    if (adjustl(resname) /= "read") molec%atom(:)%resname=resname

    ! 2. MAKE CHANGES IF REQUIRED
    ! -------------------------------
    !Swaping atoms
    if (adjustl(swapfile) /= "none") then
        open(I_SWP,file=swapfile,status='old',iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(swapfile)) )
        read(I_SWP,*) nswap
        molec_aux = molec
        do i=1,nswap
            read(I_SWP,*) iat_orig, iat_new
            molec_aux%atom(iat_new)=molec%atom(iat_orig)
        enddo
        molec = molec_aux
        close(I_SWP)
    endif
    !Forcing element names (no FF labels)
    if (use_elements) molec%atom(:)%name = molec%atom(:)%element
    !Removing center of gravity
    if (remove_com) then
        call atname2element(molec)
        call assign_masses(molec)
        call get_com(molec)
        i=molec%natoms
        molec%atom(1:i)%x = molec%atom(1:i)%x - molec%comX
        molec%atom(1:i)%y = molec%atom(1:i)%y - molec%comY
        molec%atom(1:i)%z = molec%atom(1:i)%z - molec%comZ
    endif

    ! 3. WRITE OUTPUT
    ! ---------------------------------
    if (overwrite) stat="unknown"
    open(O_OUT,file=outfile,status=stat,iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Cannot write in "//trim(adjustl(outfile))//&
                                                            ". Already exists? Use -r to overwrite" )

    if (adjustl(filetype_out) == "guess") call split_line_back(outfile,".",null,filetype_out)
    !If PDB check if connections are requested
    if (make_connect .and. adjustl(filetype_out)=="pdb") then
        call guess_connect(molec)
        filetype_out="pdb-c"
    endif
    call generic_strfile_write(O_OUT,filetype_out,molec) 
    close(O_OUT)

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype_inp,outfile,filetype_out,overwrite,make_connect,&
                           use_elements,remove_com,resname,swapfile)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,outfile,&
                                          filetype_inp,filetype_out, &
                                          resname,swapfile
        logical,intent(inout) :: overwrite, make_connect, use_elements, &
                                 remove_com
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
                case ("-fti") 
                    call getarg(i+1, filetype_inp)
                    argument_retrieved=.true.

                case ("-o") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                case ("-fto") 
                    call getarg(i+1, filetype_out)
                    argument_retrieved=.true.

                case ("-r")
                    overwrite=.true.

                case ("-rn")
                    call getarg(i+1, resname)
                    argument_retrieved=.true.

                case ("-connect")
                    make_connect=.true.

                case ("-use-elems")
                    use_elements=.true.

                case ("-rmcom")
                    remove_com=.true.

                case ("-swap")
                    call getarg(i+1, swapfile)
                    argument_retrieved=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        !The default output is now xyz
        if (adjustl(outfile) == "default") then
            call split_line_back(inpfile,".",outfile,null)
            outfile=trim(adjustl(outfile))//".xyz"
        endif

        ! Some checks on the input
        !----------------------------

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '        F O R M A T    C O N V E R T E R '    
        write(0,'(/,A)') '        Convert between geometry formats '  
        write(0,'(/,A)') '         Revision: fchk2gro-150129               '         
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-fti            ', trim(adjustl(filetype_inp))
        write(0,*) '-o              ', trim(adjustl(outfile))
        write(0,*) '-fto            ', trim(adjustl(filetype_out))
        write(0,*) '-r             ',  overwrite
        write(0,*) '-swap           ',  trim(adjustl(swapfile))
        write(0,*) '-rn             ',  trim(adjustl(resname))
        write(0,*) '-connect       ',  make_connect
        write(0,*) '-use-elems     ',  use_elements
        write(0,*) '-rmcom         ',  remove_com
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(inout) :: filetype
        type(str_resmol),intent(inout) :: molec

        !Local
        type(str_molprops) :: props

        ! Predefined filetypes
        select case (adjustl(filetype))
            case("gro")
             call read_gro(unt,molec)
            case("pdb")
             call read_pdb_new(unt,molec)
            case("g96")
             call read_g96(unt,molec)
            case("xyz")
             call read_xyz(unt,molec)
            case("log")
             call parse_summary(unt,molec,props,"struct_only")
             molec%atom(:)%resname = "UNK"
             molec%atom(:)%resseq = 1
            case("stdori")
             call get_first_std_geom(unt,molec)
             molec%atom(:)%resname = "UNK"
             molec%atom(:)%resseq = 1
            case("inpori")
             call get_first_inp_geom(unt,molec)
             molec%atom(:)%resname = "UNK"
             molec%atom(:)%resseq = 1
            case("fchk")
             call read_fchk_geom(unt,molec)
             molec%atom(:)%resname = "UNK"
             molec%atom(:)%resseq = 1
            case("UnSym")
             call read_molcas_geom(unt,molec)
            case("molden")
             call read_molden(unt,molec)
            case("psi4")
             call read_psi_geom(unt,molec)
            case("gamess")
             call read_gamess_geom(unt,molec)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

        call atname2element(molec)

        return

    end subroutine generic_strfile_read


    subroutine generic_strfile_write(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(inout) :: filetype
        type(str_resmol),intent(inout) :: molec

        ! Predefined filetypes
        select case (adjustl(filetype))
            case("gro")
             call write_gro(unt,molec)
            case("pdb")
             call write_pdb(unt,molec)
            case("pdb-c")
             call write_pdb_connect(unt,molec)
            case("g96")
             call write_g96(unt,molec)
            case("xyz")
             call write_xyz(unt,molec)
            case("fchk")
             call element2AtNum(molec)
             call write_fchk_geom(unt,molec)
            case("fcc")
             call write_fcc(unt,molec)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

        return

    end subroutine generic_strfile_write 


end program fchk2gro

