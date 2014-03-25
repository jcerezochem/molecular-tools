program gen_oniom

    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    !
    ! Compilation instructions (for mymake script):
    !make$ gfortran ../modules/alerts.f90 ../modules/structure_types_v2.f90 ../modules/line_preprocess.f90 ../modules/gro_manage_v3.f90 ../modules/gaussian_manage_v2.f90 ../modules/pdb_manage_v2.f90 ../modules/molpro_manage_v2.f90 gen_oniom_v3.f90 -o gen_oniom_v3.exe -cpp
    !
    ! Change log:
    !============================================================================    
    ! Version 3 (July '13)
    !   V2 is skipped (this comes from V1)
    !   Includes support for standard GMX paths (if they are exported)
    !
    !   V4: addapted to v4 subroutines

    use structure_types
    use line_preprocess
    use gro_manage
    use gaussian_manage
    use pdb_manage
    use alerts

    implicit none

    common /ALERT_COUNT/ n_notes, n_errors

    !====================== 
    !Number of error/notes
    integer :: n_notes=0, n_errors=0
    !====================== 

    !====================== 
    !System variables
    type(str_resmol) :: molec !, molecUA
    type(str_resmol),dimension(1000) :: residue
    !
    character(len=5),dimension(1000) :: molname
    integer :: nmol
    integer,dimension(500) :: molmap, resdone=0
    integer,dimension(10000) :: frz
    !
    integer :: ires, imap
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: nelements
    character(len=1) :: null
    character(len=4) :: atname
    character(len=16) :: dummy_char
    integer :: dummy_int
    logical :: is_done
    !====================== 

    !=============
    !Counters
    integer :: i,j,k, ii
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_TOP=11,  &
               I_NDX=12,  &
               O_GJF=20
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile="input.gro",  &
                         topfile="topol.top",  &
                         ndxfile="none",       &
                         outgau="out.com" 
    !reading
    character(len=1000) :: ndxrecord    
    character(len=260)  :: line
    !status
    integer :: IOstatus
    !Control stdout
    logical :: verbose=.false.
    !===================

    !Tuning output
    character(len=50) :: resname="NKX"
    character(len=200):: jobspec="#p ONIOM(td pbe1pbe/6-31G(d):amber)=EmbedCharge Geom=connectivity"
    character(len=4)  :: nproc="8"
    character(len=50) :: mm_file="none"
    character(len=5)  :: chrgspin_h="0 1 ", &
                         chrgspin_l="0 1 "



    ! 0. GET COMMAND LINE ARGUMENTS AND OPEN FILES
    call parse_input(inpfile,filetype,topfile,ndxfile,resname,jobspec,outgau,nproc,mm_file,chrgspin_h,chrgspin_l,verbose)

    open(I_INP,file=inpfile,iostat=IOstatus,status="old")
    if (IOstatus /= 0) call alert_msg("fatal","Unable to open "//trim(adjustl(inpfile)))
    open(I_TOP,file=topfile,iostat=IOstatus,status="old")
    if (IOstatus /= 0) call alert_msg("fatal","Unable to open "//trim(adjustl(topfile)))


    ! 1. READ DATA
    call generic_strfile_read(I_INP,filetype,molec)
    call read_top(I_TOP,residue,molname,nmol)
    close(I_INP)
    close(I_TOP)
    !========================0
    !should THIS be in read_top?
    !Identify different molecules (like that if splitted in [ molecules ] will not work.
    imap=1
    k=0
    molmap=1
    do i=2,nmol
        if (molname(i) == molname(i-1)) molmap(imap)=molmap(imap)+1
        if (molname(i) /= molname(i-1)) imap=imap+1
    enddo
    !Atom names should be elemet names
    do i=1,imap
! print*, residue(i)%name
        do j=1,residue(i)%natoms
            !NOT VALID FOR TWO CHAR NAMES (NEW: introducing some exceptions)
            if (residue(i)%atom(j)%name(1:2) == "Cl") then
                residue(i)%atom(j)%name="Cl"
            elseif (residue(i)%atom(j)%name(1:2) == "NA" .or. residue(i)%atom(j)%name(1:2) == "Na") then
                residue(i)%atom(j)%name="Na"
            else
                residue(i)%atom(j)%name = residue(i)%atom(j)%name(1:1)
            endif
        enddo
    enddo
    ! Merge molec and top info (top is kept in case of non-coincident info)
    k=0
    ires=0
    do i=1,imap
    do ii=1,molmap(i)
        ires=ires+1
        do j=1,residue(i)%natoms
            k=k+1
            residue(i)%atom(j)%x = molec%atom(k)%x
            residue(i)%atom(j)%y = molec%atom(k)%y
            residue(i)%atom(j)%z = molec%atom(k)%z
            residue(i)%atom(j)%resseq=ires
            molec%atom(k) = residue(i)%atom(j)
        enddo
    enddo
    enddo
    molec%natoms = k
    !========================0

    ! 2. SELECT LAYERS (and frz index)
    frz=-1
    do i=1,molec%natoms
        if (adjustl(molec%atom(i)%resname) == adjustl(resname)) then
            molec%atom(i)%chain="H"
        else
            molec%atom(i)%chain="L"
        endif
    enddo
    !Use ndxfile
    if (adjustl(ndxfile) /= "none") then
        ndxrecord=""
        open(I_NDX,file=ndxfile,status="old",iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg("fatal","Unable to open "//trim(adjustl(ndxfile)))
        !Read file
        read(I_NDX,'(A)') line
        line=adjustl(line)
        if (line(1:1) /= '[') call alert_msg("warning","Unexpected begining of ndx file")
        do 
            read(I_NDX,'(A)',iostat=IOstatus) line
            if (IOstatus /= 0) exit
            line=adjustl(line)
            if (line(1:1) == '[') exit
            ndxrecord=trim(adjustl(ndxrecord))//" "//trim(adjustl(line))
        enddo
        call read_list_int(ndxrecord,nelements,molmap)
!         !Set the ndx-selected atoms and corresponing molec to H (ensures whole residues, even if they are not in the ndx!)
!         imap=0
!         do i=1,nelements
!             k = molec%atom(molmap(i))%resseq
!             is_done=.false.
!             do j=1,imap
!                 if (k == resdone(j)) is_done=.true.
!             enddo
!             if (is_done) cycle
!             imap = imap+1
!             resdone(imap)=k
!             if (verbose) print'(A,I5,A)', "Adding residue ",k," to the QM layer"
!             do j=1,molec%natoms
!                 if (molec%atom(j)%resseq == k) molec%atom(j)%chain='H'
!             enddo
!         enddo
        !Set the ndx-selected atoms and corresponing molec to H (assuming whole residues in ndx (more efficient...)
        imap=0
        do i=1,nelements
            j=molmap(i)
            if (verbose) print'(A,I5,A)', "Adding atom ",j," to the QM layer"
            molec%atom(j)%chain='H'
            !Also unfreeze (to be generalized)
            frz(j) = 0
        enddo
!         call sort_vec_int(selres,nelements) !it is already ordered in ndx file
     endif
        

    ! 3. WRITE OUT FILE (gaussian)
    !   (To be done on gcom_write subroutine)
    open(O_GJF,file=outgau,status="replace")
    !TO BE PLACED IN A SR
    !Link0 section
    write(O_GJF,'(A)') "%mem=2GB"
    call split_line(outgau,'.',outgau,dummy_char)
    write(O_GJF,'(A)') "%chk="//trim(adjustl(outgau))//".chk"
    write(O_GJF,'(A)') "%nproc="//trim(adjustl(nproc))
    !Route section
    write(O_GJF,'(A)') trim(adjustl(jobspec))
    write(O_GJF,'(A)') ""
    !Title
    write(O_GJF,'(A)') "File generated from "//trim(adjustl(inpfile))
    write(O_GJF,'(A)') ""
    !Charge and multiplicity
    write(O_GJF,'(X,A,X,A)') chrgspin_l,chrgspin_h
    !Geomertry
    do i=1,molec%natoms
        !First form: Atom-Type-Charge
        write(dummy_char,'(F8.4)') molec%atom(i)%q
        dummy_char = adjustl(trim(molec%atom(i)%name))//"-"//&
                     adjustl(trim(molec%atom(i)%attype))//"-"//&
                     adjustl(trim(dummy_char))
        write(O_GJF,100) dummy_char, molec%atom(i)%x,molec%atom(i)%y,molec%atom(i)%z,molec%atom(i)%chain
!         write(O_GJF,101) dummy_char, frz(i), molec%atom(i)%x,molec%atom(i)%y,molec%atom(i)%z,molec%atom(i)%chain
100 format(A20,X,3(F15.6,X),A)
101 format(A20,X,I2,X,3(F15.6,X),A)
    enddo
    write(O_GJF,*) "" 
    !Geom connectivity (faked to avoid missings) -- it can be the true one! why not? TBD
    do i=1,molec%natoms
        write(O_GJF,'(I5)') i
    enddo
    write(O_GJF,*) "" 
    if ( adjustl(mm_file) /= "none" ) &
    write(O_GJF,'(A)') "@"//trim(adjustl(mm_file))
    write(O_GJF,*) "" 

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,topfile,ndxfile,resname,jobspec,outfile,nproc,mm_file,chrgspin_h,chrgspin_l,verbose)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        common /GLOBAL_OPTS/ do_refine_charges

        character(len=*),intent(inout) :: inpfile,topfile,filetype,ndxfile,outfile,&
                                          resname,jobspec,nproc,mm_file,chrgspin_h,&
                                          chrgspin_l
        logical,intent(inout) :: verbose
        ! Local
        logical :: argument_retrieved,  &
                   need_help = .false., &
                   do_refine_charges
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
                    call getarg(i+1, filetype)
                    argument_retrieved=.true.

                case ("-o") 
                    call getarg(i+1, outfile)
                    argument_retrieved=.true.

                case ("-p") 
                    call getarg(i+1, topfile)
                    argument_retrieved=.true.

                case ("-n") 
                    call getarg(i+1, ndxfile)
                    argument_retrieved=.true.

                case ("-job") 
                    call getarg(i+1, jobspec)
                    argument_retrieved=.true.

                case ("-res")
                    call getarg(i+1, resname)
                    argument_retrieved=.true.

                case ("-nproc")
                    call getarg(i+1, nproc)
                    argument_retrieved=.true.

                case ("-mmf")
                    call getarg(i+1, mm_file)
                    argument_retrieved=.true.

                case ("-cs-l")
                    call getarg(i+1, chrgspin_l)
                    argument_retrieved=.true.

                case ("-cs-h")
                    call getarg(i+1, chrgspin_h)
                    argument_retrieved=.true.

                case ("-v") 
                    verbose=.true.

                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 
          

       !Print options (to stderr)
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,'(/,A)') '              G E N _ O N I O M '    
        write(0,'(/,A)') '         Revision: gen_oniom-140225               '           
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(filetype))
        write(0,*) '-p              ', trim(adjustl(topfile))
        write(0,*) '-n              ', trim(adjustl(ndxfile))
        write(0,*) '-o              ', trim(adjustl(outfile))
        write(0,*) '-job            ', trim(adjustl(jobspec))
        write(0,*) '-res            ', trim(adjustl(resname))
        write(0,*) '-nproc          ', trim(adjustl(nproc))
        write(0,*) '-mmf            ', trim(adjustl(mm_file))
        write(0,*) '-cs-h           ', trim(adjustl(chrgspin_h))
        write(0,*) '-cs-l           ', trim(adjustl(chrgspin_l))
        write(0,*) '-v             ',  verbose
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

        if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line(inpfile,".",null,filetype)
        select case (adjustl(filetype))
            case("gro")
             call read_gro(I_INP,molec)
            case("g96")
             call read_g96(I_INP,molec)
            case("pdb")
             call read_pdb_new(I_INP,molec)
            case("log")
             call parse_summary(I_INP,molec,props,"struct_only")
!             case("fchk")
!              call read_fchk_geom(I_INP,molec)
            case default
             call alert_msg("fatal","Trying to guess, but file type but not known: "//adjustl(trim(filetype))&
                        //". Try forcing the filetype with -ft flag (available: gro, g96, pdb)")
        end select

        else
        ! Predefined filetypes
        select case (adjustl(filetype))
            case("gro")
             call read_gro(I_INP,molec)
            case("g96")
             call read_g96(I_INP,molec)
            case("pdb")
             call read_pdb_new(I_INP,molec)
            case("log")
             call parse_summary(I_INP,molec,props,"struct_only")
!             case("fchk")
!              call read_fchk_geom(I_INP,molec)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select
        endif


        return


    end subroutine generic_strfile_read

end program gen_oniom

