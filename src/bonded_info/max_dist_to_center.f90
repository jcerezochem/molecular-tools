program max_dist_to_center


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
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
    use molcas_unsym_manage
    use ff_build
    use molecular_structure

    implicit none

    !====================== 
    !System variables
    type(str_resmol) :: molec, &
                        molec_filt
    type(str_atom) :: atom
    !====================== 

    !=============
    !Counters and dummies
    integer :: i,j,k,l, jj,kk, iat
    character(len=1) :: null
    logical :: overwrite=.false.,&
               make_connect=.false.
    !Swap related counters
    integer :: iat_new, iat_orig, nswap
    !filter related
    character(len=50) :: filter="all"
    integer,dimension(100) :: listfilter
    integer :: Nfilter
    ! Dist calc
    real(8) :: d, dmax
    integer :: imax
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
    call parse_input(inpfile,filetype_inp,filter)

 
    ! 1. READ INPUT
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(filetype_inp) == "guess") call split_line_back(inpfile,".",null,filetype_inp)
    call generic_strfile_read(I_INP,filetype_inp,molec)
    close(I_INP)
    !Option to specify the resname from command line
    if (adjustl(resname) /= "read") molec%atom(:)%resname=resname
    !Get massses
    call atname2element(molec)
    call assign_masses(molec)

    ! Filtering if needed
    if (adjustl(filter) == "all") then
        molec_filt = molec
    else
        call selection2intlist(filter,listfilter,Nfilter)
        do i=1,Nfilter
            j = listfilter(i)
            molec_filt%atom(i) = molec%atom(j)
        enddo
        molec_filt%natoms = Nfilter
    endif

    ! 2. MAKE THE THING
    ! ------------------------------
    call get_cog(molec_filt)
    call get_com(molec_filt)

    dmax = 0.d0
    imax = -1
    do i=1,molec%natoms
        d = dsqrt((molec%atom(i)%x - molec_filt%cogX)**2 + &
                  (molec%atom(i)%y - molec_filt%cogY)**2 + &
                  (molec%atom(i)%z - molec_filt%cogZ)**2 )
        if ( d > dmax ) then
            atom=molec%atom(i)
            imax = i
            dmax = d
        endif
    enddo

    print*, "CENTER OF GEOMETRY"
    print*, molec_filt%cogX, molec_filt%cogY, molec_filt%cogZ
    print*, "MAX DISTANCE TO CENTER OF GEOMETRY"
    print*, dmax
    print'(X,A,I0,A,I0,A)', "(atom: ", imax, "; residue: ", atom%resseq, ")"

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype_inp,filter)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype_inp,filter
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

                case ("-filter") 
                    call getarg(i+1, filter)
                    argument_retrieved=.true.
        
                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        ! Some checks on the input
        !----------------------------

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '      D I S T   T O   C E N T E R S '    
        write(0,'(/,A)') '      Compute dist to COM and COG'  
        write(0,'(/,A)') '        Revision: max_dist_to_center            '         
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-fti            ', trim(adjustl(filetype_inp))
        write(0,*) '-filter         ', trim(adjustl(filter))
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
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

        call atname2element(molec)

        return

    end subroutine generic_strfile_read

    subroutine selection2intlist(selection,list,Nlist)

        use line_preprocess

        character(len=*), intent(in) :: selection
        integer, intent(out) :: Nlist
        integer,dimension(1:100) :: list
        !local 
        character(len=5),dimension(100) :: selection_split
        integer :: i, j 
        integer :: N, range_last, range_width
        logical :: is_range

        call string2vector_char_new(selection,selection_split,N," ")

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


end program max_dist_to_center

