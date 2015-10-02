program rmsd_fit


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
    type(str_resmol) :: molec, molecRef, &
                        molec_filt, molecRef_filt
    !====================== 

    !=============
    !Counters and dummies
    integer :: i,j,k,l, jj,kk, iat
    character(len=1) :: null
    logical :: overwrite    = .false. ,&
               make_connect = .false. ,&
               use_elements = .false. 
    character(len=3) :: remove_mode = "COM"
    !filter related
    character(len=50) :: filter="all"
    integer,dimension(100) :: listfilter
    integer :: Nfilter
    ! Rot matri
    real(8),dimension(3,3) :: Rot
    ! COM for reference (to restore at the end)
    real(8) :: Xref, Yref, Zref, &
               Xmol, Ymol, Zmol
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
                         filetype_ref="guess",&
                         filetype_out="guess"
    character(len=200):: inpfile="input.fchk",&
                         reffile="ref.fchk"  ,&
                         outfile="default"   
    !status
    integer :: IOstatus
    character(len=7) :: stat="new" !do not overwrite when writting
    !===================


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype_inp,reffile,filetype_ref,  &
                     outfile,filetype_out,overwrite,make_connect,&
                     use_elements,remove_mode,resname,filter)
    ! Set options to upper case
    call set_word_upper_case(remove_mode)
 
    ! 1. READ INPUT
    ! ---------------------------------
    ! 1a. Rotable molecule
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(filetype_inp) == "guess") call split_line_back(inpfile,".",null,filetype_inp)
    call generic_strfile_read(I_INP,filetype_inp,molec)
    close(I_INP)
    !Option to specify the resname from command line
    if (adjustl(resname) /= "read") molec%atom(:)%resname=resname

    ! 1b. Refence molecule
    open(I_INP,file=reffile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(reffile)) )

    if (adjustl(filetype_ref) == "guess") call split_line_back(reffile,".",null,filetype_ref)
    call generic_strfile_read(I_INP,filetype_ref,molecRef)
    close(I_INP)
    !Option to specify the resname from command line
    if (adjustl(resname) /= "read") molecRef%atom(:)%resname=resname

    ! 2. MAKE CHANGES IF REQUIRED
    ! -------------------------------
    !Forcing element names (no FF labels)
    if (use_elements) molec%atom(:)%name = molec%atom(:)%element

    ! Filtering if needed
    if (adjustl(filter) == "all") then
        molec_filt = molec
        molecRef_filt = molecRef
    else
        call selection2intlist(filter,listfilter,Nfilter)
        do i=1,Nfilter
            j = listfilter(i)
            molec_filt%atom(i) = molec%atom(j)
            molecRef_filt%atom(i) = molecRef%atom(j)
        enddo
        molec_filt%natoms = Nfilter
        molecRef_filt%natoms = Nfilter
    endif

    !Removing center of gravity
    if (remove_mode == "COM") then
        call atname2element(molec_filt)
        call assign_masses(molec_filt)
        call atname2element(molecRef_filt)
        call assign_masses(molecRef_filt)
        call get_com(molec_filt)
        call get_com(molecRef_filt)
        i=molec_filt%natoms
        molec_filt%atom(1:i)%x = molec_filt%atom(1:i)%x - molec_filt%comX
        molec_filt%atom(1:i)%y = molec_filt%atom(1:i)%y - molec_filt%comY
        molec_filt%atom(1:i)%z = molec_filt%atom(1:i)%z - molec_filt%comZ
        Xmol = molec_filt%comX
        Ymol = molec_filt%comY
        Zmol = molec_filt%comZ
        molecRef_filt%atom(1:i)%x = molecRef_filt%atom(1:i)%x - molecRef_filt%comX
        molecRef_filt%atom(1:i)%y = molecRef_filt%atom(1:i)%y - molecRef_filt%comY
        molecRef_filt%atom(1:i)%z = molecRef_filt%atom(1:i)%z - molecRef_filt%comZ
        Xref = molecRef_filt%comX
        Yref = molecRef_filt%comY
        Zref = molecRef_filt%comZ
    else if (remove_mode == "COG") then
        call get_cog(molec_filt)
        call get_cog(molecRef_filt)
        i=molec_filt%natoms
        molec_filt%atom(1:i)%x = molec_filt%atom(1:i)%x - molec_filt%cogX
        molec_filt%atom(1:i)%y = molec_filt%atom(1:i)%y - molec_filt%cogY
        molec_filt%atom(1:i)%z = molec_filt%atom(1:i)%z - molec_filt%cogZ
        Xmol = molec_filt%cogX
        Ymol = molec_filt%cogY
        Zmol = molec_filt%cogZ
        molecRef_filt%atom(1:i)%x = molecRef_filt%atom(1:i)%x - molecRef_filt%cogX
        molecRef_filt%atom(1:i)%y = molecRef_filt%atom(1:i)%y - molecRef_filt%cogY
        molecRef_filt%atom(1:i)%z = molecRef_filt%atom(1:i)%z - molecRef_filt%cogZ
        Xref = molecRef_filt%cogX
        Yref = molecRef_filt%cogY
        Zref = molecRef_filt%cogZ
    else
        Xref = 0.d0
        Yref = 0.d0
        Zref = 0.d0
        Xmol = 0.d0
        Ymol = 0.d0
        Zmol = 0.d0
    endif


    ! 3. FITTING
    ! ---------------------------------
    call ROTATA1(molec_filt,molecRef_filt,Rot)
    print*, "Rotation matrix: "
    do i=1,3
        print'(3F10.3)', Rot(1:3,i)
    enddo
    print*, ""
    ! Place the whole molec according to COM mol
    do i=1,molec%natoms
        molec%atom(i)%x = molec%atom(i)%x - Xmol
        molec%atom(i)%y = molec%atom(i)%y - Ymol
        molec%atom(i)%z = molec%atom(i)%z - Zmol
    enddo
    ! Rotate
    call rotate_molec(molec,Rot)
    call rotate_molec(molec_filt,Rot)
    ! Restore to ref original position
    do i=1,molec%natoms
        molec%atom(i)%x = molec%atom(i)%x + Xref
        molec%atom(i)%y = molec%atom(i)%y + Yref
        molec%atom(i)%z = molec%atom(i)%z + Zref
    enddo

    ! 4. WRITE OUTPUT
    ! ---------------------------------
    if (overwrite) stat="unknown"
    open(O_OUT,file=outfile,status=stat,iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Cannot write in "//trim(adjustl(outfile))//&
                                                            ". Already exists? Use -ow to overwrite" )

    if (adjustl(filetype_out) == "guess") call split_line_back(outfile,".",null,filetype_out)
    !If PDB check if connections are requested
    if (make_connect .and. adjustl(filetype_out)=="pdb") then
        call guess_connect(molec)
        filetype_out="pdb-c"
    endif
    call generic_strfile_write(O_OUT,filetype_out,molec) 
    close(O_OUT)

    call generic_strfile_write(80,"xyz",molec_filt)
    call generic_strfile_write(90,"xyz",molecRef_filt)

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype_inp,reffile,filetype_ref,  &
                           outfile,filetype_out,overwrite,make_connect,&
                           use_elements,remove_mode,resname,filter)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,outfile,reffile, &
                                          filetype_inp,filetype_out, &
                                          filetype_ref,&
                                          resname,remove_mode, filter
        logical,intent(inout) :: overwrite, make_connect, use_elements
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

                case ("-r") 
                    call getarg(i+1, reffile)
                    argument_retrieved=.true.
                case ("-ftr") 
                    call getarg(i+1, filetype_ref)
                    argument_retrieved=.true.

                case ("-o") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.
                case ("-fto") 
                    call getarg(i+1, filetype_out)
                    argument_retrieved=.true.

                case ("-filter") 
                    call getarg(i+1, filter)
                    argument_retrieved=.true.

                case ("-ow")
                    overwrite=.true.

                case ("-rn")
                    call getarg(i+1, resname)
                    argument_retrieved=.true.

                case ("-connect")
                    make_connect=.true.

                case ("-use-elems")
                    use_elements=.true.

                case ("-rmode")
                    call getarg(i+1, remove_mode)
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
            outfile=trim(adjustl(outfile))//"_rot.xyz"
        endif

        ! Some checks on the input
        !----------------------------

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '         R M S D   F I T '    
        write(0,'(/,A)') '         Fit with ROTATA '  
        write(0,'(/,A)') '       Revision: rmsd_fit-150911                  '         
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-fti            ', trim(adjustl(filetype_inp))
        write(0,*) '-o              ', trim(adjustl(outfile))
        write(0,*) '-fto            ', trim(adjustl(filetype_out))
        write(0,*) '-r              ', trim(adjustl(reffile))
        write(0,*) '-ftr            ', trim(adjustl(filetype_ref))
        write(0,*) '-rn             ',  trim(adjustl(resname))
        write(0,*) '-filter         ', trim(adjustl(filter))
        write(0,*) '-connect       ',  make_connect
        write(0,*) '-use-elems     ',  use_elements
        write(0,*) '-rmode          ', remove_mode
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(in) :: filetype
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
        character(len=*),intent(in) :: filetype
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


    subroutine stringbl2vector_char(raw_vector,array_vector,n_elem,sep)

        !Description
        ! Tranforms a string of <sep> sepparated values into an
        ! array of such characters 

        character(len=*),intent(in) :: raw_vector
        character(len=*),dimension(:),intent(out) :: array_vector
        integer,intent(out) :: n_elem
        character(len=*) :: sep !separador

        !Local
        character(len=len_trim(raw_vector)) :: raw_vector_copy
        character(len=240) :: auxchar
        integer :: i
    
        
        !Copy the original vector to avoid modifying it
        raw_vector_copy = trim(adjustl(raw_vector))

        !Read unknown length vector
        i=0
        do 
            i=i+1
            if (len_trim(raw_vector_copy) == 0) then
                i=i-1
                exit
            else if ( INDEX(raw_vector_copy,sep) /= 0 ) then
                call split_line(raw_vector_copy,sep,auxchar,raw_vector_copy)
                read(auxchar,*) array_vector(i)
                ! By adjustl-ing every time we avoid double counting blank spaces
                raw_vector_copy = adjustl(raw_vector_copy)
            else 
                read(raw_vector_copy,*) array_vector(i)
                exit
            endif
        enddo  
        n_elem=i

        return

    end subroutine stringbl2vector_char


    subroutine selection2intlist(selection,list,Nlist)

        character(len=*), intent(in) :: selection
        integer, intent(out) :: Nlist
        integer,dimension(1:100) :: list
        !local 
        character(len=5),dimension(100) :: selection_split
        integer :: i, j 
        integer :: N, range_last, range_width
        logical :: is_range

        call stringbl2vector_char(selection,selection_split,N," ")

        is_range = .false.
        j = 0
        do i=1,N
            if (selection_split(i) == "to") then
                is_range =  .true.
                cycle
            endif
            ! Read number
            if (.not.is_range) then
                j = j+1
                read(selection_split(i),*) list(j)
            else
                read(selection_split(i),*) range_last
                range_width = range_last - list(j)
                do jj = 1, range_width
                    j = j + 1
                    list(j) = list(j-1) + 1
                enddo
                is_range = .false.
            endif
        enddo
        Nlist = j

        return

    end subroutine selection2intlist


end program rmsd_fit

