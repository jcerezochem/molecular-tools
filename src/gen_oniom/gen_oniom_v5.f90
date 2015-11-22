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
    !   v4.1: Using element names where needed
    !         Solved issue: when using -n crash for large systems
    !                       with segfault. The problem is in
    !                       frz array. It is not called (should be
    !                       eliminated)
    !         Solved issue: if the ndxrecord has not atoms, exit the 
    !                       the cycle to avoid segfault
    ! (07/05/15)
    !     Add plain option, to get a "plain" input file (not oniom)
    ! Version 5
    !  Use light types
    !==============================================================

    use structure_types
    use line_preprocess
    use gro_manage_go
    use gaussian_manage, only : parse_summary
    use pdb_manage
    use alerts
    use molecular_structure

    implicit none

    !====================== 
    !System variables
    type(str_resmol_light) :: molec 
    !Residues are different molecule types (there might be several of each kind in the system).
    ! E.g. Solvent and solute are only 2 type of molecules 
    type(str_resmol),dimension(1:2) :: residue
    !
    character(len=5),dimension(10000) :: molname
    integer :: nmol
    integer,dimension(50000) :: molmap, resdone=0
    integer,dimension(10000) :: frz
    !
    integer :: ires, imap
    !
    logical :: pointcharges=.false., &
               plain=.false.
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: nelements
    character(len=1) :: null
    character(len=4) :: atname
    character(len=20) :: dummy_char
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
               I_FRZ=13,  &
               O_GJF=20
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile="input.gro",  &
                         topfile="topol.top",  &
                         ndxfile="none",       &
                         frzfile="none",       &
                         outgau="out.com" 
    !reading
    character(len=20000) :: ndxrecord    
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

    !I/O   (if not redefined here, the some labels are lost...)
    I_INP=10
    I_TOP=11
    I_NDX=12
    I_FRZ=13
    O_GJF=20

    ! 0. GET COMMAND LINE ARGUMENTS AND OPEN FILES
    call parse_input(inpfile,filetype,topfile,ndxfile,frzfile,resname,jobspec,outgau,nproc,&
                     mm_file,chrgspin_h,chrgspin_l,pointcharges,plain,verbose)
    !Check conflicts
    if (plain .and. pointcharges) then
        call alert_msg("warning","-plain and -pc cannot be call together. Setting -pc to false")
        pointcharges = .false.
    endif

    ! 1. READ DATA
    open(I_INP,file=inpfile,iostat=IOstatus,status="old")
    if (IOstatus /= 0) call alert_msg("fatal","Unable to open "//trim(adjustl(inpfile)))
    if (filetype == "guess") then
        call split_line_back(inpfile,".",null,filetype)
    endif
    call generic_strfile_read(I_INP,filetype,molec)
    close(I_INP)

    if (.not. plain) then
        if (adjustl(topfile) /= "none") then
!             allocate(residue(1:2))
!             residue(:) = molec
            open(I_TOP,file=topfile,iostat=IOstatus,status="old")
            if (IOstatus /= 0) call alert_msg("fatal","Unable to open "//trim(adjustl(topfile)))
            call read_top(I_TOP,residue,molname,nmol)
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

            ! Merge molec and top info (top is kept in case of non-coincident info)
            ! Note the coordinates from molec are kept (residue from top has not coordinates)
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
                    residue(i)%atom(j)%element=molec%atom(k)%element
                    !Get atom (need to be explicit, since in molec is light)
                    molec%atom(k)%x        = residue(i)%atom(j)%x      
                    molec%atom(k)%y        = residue(i)%atom(j)%y      
                    molec%atom(k)%z        = residue(i)%atom(j)%z      
                    molec%atom(k)%attype   = residue(i)%atom(j)%attype 
                    molec%atom(k)%fftype   = residue(i)%atom(j)%fftype 
                    molec%atom(k)%name     = residue(i)%atom(j)%name   
                    molec%atom(k)%chain    = residue(i)%atom(j)%chain  
                    molec%atom(k)%ins_code = residue(i)%atom(j)%ins_code
                    molec%atom(k)%element  = residue(i)%atom(j)%element
                    molec%atom(k)%AtNum    = residue(i)%atom(j)%AtNum  
                    molec%atom(k)%mass     = residue(i)%atom(j)%mass   
                    molec%atom(k)%q        = residue(i)%atom(j)%q      
                enddo
            enddo
            enddo
            molec%natoms = k
            !========================
        endif

        ! 2. SELECT LAYERS (and frz index)
        frz=-1
        do i=1,molec%natoms
            if (adjustl(molec%atom(i)%resname) == adjustl(resname)) then
                molec%atom(i)%chain="H"
            else
                molec%atom(i)%chain="L"
            endif
            molec%atom(i)%ins_code = "0"
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
            close(I_NDX)
            if ( trim(adjustl(ndxrecord)) /= "" ) then
                call read_list_int(ndxrecord,nelements,molmap)
!             !Set the ndx-selected atoms and corresponing molec to H (ensures whole residues, even if they are not in the ndx!)
!             imap=0
!             do i=1,nelements
!                 k = molec%atom(molmap(i))%resseq
!                 is_done=.false.
!                 do j=1,imap
!                     if (k == resdone(j)) is_done=.true.
!                 enddo
!                 if (is_done) cycle
!                 imap = imap+1
!                 resdone(imap)=k
!                 if (verbose) print'(A,I5,A)', "Adding residue ",k," to the QM layer"
!                 do j=1,molec%natoms
!                     if (molec%atom(j)%resseq == k) molec%atom(j)%chain='H'
!                 enddo
!             enddo
                !Set the ndx-selected atoms and corresponing molec to H (assuming whole residues in ndx (more efficient...)
                imap=0
                do i=1,nelements
                    j=molmap(i)
                    if (verbose) print'(A,I5,A)', "Adding atom ",j," to the QM layer"
                    molec%atom(j)%chain='H'
                    !Also unfreeze (to be generalized): why this freeze?
!       !            frz(j) = 0
                enddo
            endif
!             call sort_vec_int(selres,nelements) !it is already ordered in ndx file
        endif
    endif
    !Use freezing ndxfile
    if (adjustl(frzfile) /= "none") then
        write(0,'(/,X,A,/)') "Using freeze groups"
        ndxrecord=""
        open(I_FRZ,file=frzfile,status="old",iostat=IOstatus)
        if (IOstatus /= 0) call alert_msg("fatal","Unable to open "//trim(adjustl(frzfile)))
        !Read file
        read(I_FRZ,'(A)') line
        line=adjustl(line)
        if (line(1:1) /= '[') call alert_msg("warning","Unexpected begining of ndx file")
        do 
            read(I_FRZ,'(A)',iostat=IOstatus) line
            if (IOstatus /= 0) exit
            line=adjustl(line)
            if (line(1:1) == '[') exit
            ndxrecord=trim(adjustl(ndxrecord))//" "//trim(adjustl(line))
        enddo
        close(I_FRZ)
        if ( trim(adjustl(ndxrecord)) /= "" ) then
            call read_list_int(ndxrecord,nelements,molmap)
            !Set the ndx-selected atoms and corresponing molec to H (assuming whole residues in ndx (more efficient...)
            imap=0
            do i=1,nelements
                j=molmap(i)
                if (verbose) print'(A,I5,A)', "Adding atom ",j," to the freeze group"
                molec%atom(j)%ins_code='1'
            enddo
        endif
    else
        molec%atom(1:molec%natoms)%ins_code = ""
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
        !If pointcharges, treat L layer as point charges
        if (molec%atom(i)%chain == "L" .and. pointcharges) then
            cycle
        else
            if (pointcharges.or.plain) then
                write(O_GJF,101) molec%atom(i)%element, molec%atom(i)%x,molec%atom(i)%y,molec%atom(i)%z
            else
                !First form: Atom-Type-Charge
                write(dummy_char,'(F8.4)') molec%atom(i)%q
                dummy_char = adjustl(trim(molec%atom(i)%name))//"-"//&
                             adjustl(trim(molec%atom(i)%attype))//"-"//&
                             adjustl(trim(dummy_char))//"  "//molec%atom(i)%ins_code
                write(O_GJF,100) dummy_char, molec%atom(i)%x,molec%atom(i)%y,molec%atom(i)%z,molec%atom(i)%chain
            endif
        endif
100 format(A20,X,3(F15.6,X),A)
101 format(A5,X,3(F15.6,X))
    enddo
!   If point charges
    if (pointcharges) then
    write(O_GJF,*) "" 
    do i=1,molec%natoms
        !If pointcharges, treat L layer as point charges
        if (molec%atom(i)%chain == "L" .and. pointcharges) then
            write(O_GJF,102) molec%atom(i)%x,molec%atom(i)%y,molec%atom(i)%z,molec%atom(i)%q
        else
            cycle
        endif
102 format(4(F15.6,X))
    enddo
    endif
    write(O_GJF,*) "" 
    !Connectivity only if L layer is not pointcharges
    if (.not. pointcharges .and. .not. plain) then
        !Geom connectivity (faked to avoid missings) -- it can be the true one! why not? TBD
        do i=1,molec%natoms
            write(O_GJF,'(I5)') i
        enddo
        write(O_GJF,*) "" 
    endif
    if ( adjustl(mm_file) /= "none" ) &
    write(O_GJF,'(A)') "@"//trim(adjustl(mm_file))
    write(O_GJF,*) "" 

!     deallocate(residue)
    !Review notes and errors
    if (n_notes /= 0) &
     print'(A,I0,A,/)', "There were ", n_notes, " NOTES in this run"
    if (n_errors /= 0) &
     print'(A,I0,A,/)', "There were ", n_errors, " WARNINGS in this run"

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,topfile,ndxfile,frzfile,resname,jobspec,outfile,nproc,&
                           mm_file,chrgspin_h,chrgspin_l,pointcharges,plain,verbose)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,topfile,filetype,ndxfile,frzfile,outfile,&
                                          resname,jobspec,nproc,mm_file,chrgspin_h,&
                                          chrgspin_l
        logical,intent(inout) :: verbose, pointcharges, plain
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

                case ("-frz") 
                    call getarg(i+1, frzfile)
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

                case ("-pc")
                    pointcharges=.true.

                case ("-plain")
                    plain=.true.

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
        write(0,'(/,A)') '         Revision: gen_oniom-141105              '           
        write(0,'(/,A)') '--------------------------------------------------'
        write(0,*) '-f              ', trim(adjustl(inpfile))
        write(0,*) '-ft             ', trim(adjustl(filetype))
        write(0,*) '-p              ', trim(adjustl(topfile))
        write(0,*) '-n              ', trim(adjustl(ndxfile))
        write(0,*) '-frz            ', trim(adjustl(frzfile))
        write(0,*) '-o              ', trim(adjustl(outfile))
        write(0,*) '-job            ', trim(adjustl(jobspec))
        write(0,*) '-res            ', trim(adjustl(resname))
        write(0,*) '-nproc          ', trim(adjustl(nproc))
        write(0,*) '-mmf            ', trim(adjustl(mm_file))
        write(0,*) '-cs-h           ', trim(adjustl(chrgspin_h))
        write(0,*) '-cs-l           ', trim(adjustl(chrgspin_l))
        write(0,*) '-pc            ',  pointcharges
        write(0,*) '-v             ',  verbose
        write(0,*) '-h             ',  need_help
        write(0,*) '--------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(inout) :: filetype
        type(str_resmol_light),intent(inout) :: molec

        !Local
        type(str_molprops) :: props

        ! Guess file type
        select case (adjustl(filetype))
!             case("gro")
!              call read_gro(unt,molec)
            case("g96")
             call read_g96(unt,molec)
!             case("pdb")
!              call read_pdb_new(unt,molec)
!             case("log")
!              call parse_summary(unt,molec,props,"struct_only")
!             case("fchk")
!              call read_fchk_geom(unt,molec)
            case default
             call alert_msg("fatal","Trying to guess, but file type but not known: "//adjustl(trim(filetype))&
                        //". Try forcing the filetype with -ft flag (available: gro, g96, pdb)")
        end select

        call atname2element_light(molec)

        return


    end subroutine generic_strfile_read

end program gen_oniom

