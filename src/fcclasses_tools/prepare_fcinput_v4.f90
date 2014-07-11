program prepare_fcinput

    !Description
    !===========
    ! Program that prepares the input for fcclasses retrieving the
    ! some info from a log or fchk file
    !
    ! History
    ! v4: adapted to work with modules v4

    use structure_types
    use line_preprocess
    use ff_build
    use molecular_structure
    use alerts
    use constants
    use gaussian_manage
    use gaussian_fchk_manage

    implicit none

    type(str_resmol) :: molec
    integer :: N, Nlines, i, imax, imin

    real(8) :: ti, tf
    
    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               O_OUT=20
    !files
    character(len=10) :: filetype_inp="guess"
    character(len=200):: inpfile ="input.log",  &
                         outfile,               &!="fcclasses.inp", &
                         basename
!    
    !Control of stdout
    logical :: verbose=.false.
    !status
    integer :: IOstatus

    !Auxiliar
    character :: cnull
    !===================

    call cpu_time(ti)    

    call getarg(1, inpfile)
    
    ! 1. READ INPUT
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    if (adjustl(filetype_inp) == "guess") call split_line_back(inpfile,".",cnull,filetype_inp)
    call generic_strfile_read(I_INP,filetype_inp,molec)
    close(I_INP)

    !Assign masses
    call assign_masses(molec)

    ! 2. WRITE OUTPUT
    ! --------------------------------- 
    call split_line_back(inpfile,".",basename,cnull)
    outfile="fcc_"//trim(adjustl(basename))//".inp"
    open(O_OUT,file=outfile,status='unknown')
    write(O_OUT,*) molec%natoms
    write(O_OUT,*) 3*molec%natoms-6
    do i=1,molec%natoms
        write(O_OUT,*) molec%atom(i)%mass
    enddo
    write(O_OUT,*) "ENERG Adiabatic/Vertical Energy"
    write(O_OUT,*) "'abs' 'ECDNO ' 'FC'"
    write(O_OUT,*) "TEMP.0d0 1.d-1 ! Temp(K) / BoltzThr (for TI)"
    write(O_OUT,*) "'TD' 4096 1000.d0 2"
    write(O_OUT,*) "'../../../state_"//trim(adjustl(basename))//"_fchk'"
    write(O_OUT,*) "'../../../state_"//trim(adjustl(basename))//"_fchk' 'CALCTYPE'  -1.d5  1.d5"
    write(O_OUT,*) "'../../../eldip_"//trim(adjustl(basename))//"_fchk'"
    write(O_OUT,*) "'magdip_"//trim(adjustl(basename))//"_fchk'"
    write(O_OUT,*) "1 rotation to overlap S1 and S2 (1=yes)"
    write(O_OUT,*) "20"
    write(O_OUT,*) "13"
    write(O_OUT,*) "1.d8"
    write(O_OUT,*) "'Gau' 1.50d0  4.00d0 1001 0.01d0"
    write(O_OUT,*) "'NO' 'SELECT' 0.96d0"
    write(O_OUT,*) "'../../duschinsky.dat'"
    write(O_OUT,*) "'../../displacement.dat'" 

    close(O_OUT)

    stop

    contains

    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(inout) :: filetype
        type(str_resmol),intent(inout) :: molec

        !local
        type(str_molprops) :: props

        ! Predefined filetypes
        select case (adjustl(filetype))
!             case("gro")
!              call read_gro_def(unt,molec)
!             case("pdb")
!              call read_pdb_new(unt,molec)
!             case("g96")
!              call read_g96(unt,molec)
!             case("xyz")
!              call read_xyz(unt,molec)
            case("log")
             call parse_summary(unt,molec,props,"struct_only")
!             molec%atom(:)%resname = "UNK"
!             molec%atom(:)%resseq = 1
!             case("stdori")
!              call get_first_std_geom(unt,molec)
!              molec%atom(:)%resname = "UNK"
!              molec%atom(:)%resseq = 1
            case("fchk")
             call read_fchk_geom(unt,molec)
             molec%atom(:)%resname = "UNK"
             molec%atom(:)%resseq = 1
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

        call atname2element(molec)

        return

    end subroutine generic_strfile_read


end program prepare_fcinput
