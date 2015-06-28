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

    use structure_types
    use line_preprocess
    use gro_manage
    use gaussian_manage, only : parse_summary
    use pdb_manage
    use alerts
    use molecular_structure

    implicit none

    !====================== 
    !System variables
    type(str_resmol) :: molec, molec_L, molec_H
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
    logical :: plain=.false.
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: nelements
    character(len=1) :: null
    character(len=4) :: atname
    character(len=16) :: dummy_char
    integer :: dummy_int
    logical :: is_done
    character(len=5) :: resname
    !====================== 

    !=============
    !Counters
    integer :: i,j,k, ii
    !=============

    !================
    !I/O stuff 
    !unitsprint*, molec%atom(1)%element, molec%atom(1)%name
    integer :: I_INP=10,  &
               I_TOP=11,  &
               I_NDX=12,  &
               O_CUB=20
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile="input.gro",  &
                         topfile="topol.top",  &
                         ndxfile="none",       &
                         outcube="out.cube" 
    !reading
    character(len=20000) :: ndxrecord    
    character(len=260)  :: line
    !status
    integer :: IOstatus
    !Control stdout
    logical :: verbose=.false.


    !Elec pot things
    real(8) :: xmax=-9999.d0, ymax=-9999.d0, zmax=-9999.d0, &
               x0= 9999.d0, y0= 9999.d0, z0= 9999.d0, &
               d1x,d1y,d1z, &
               d2x,d2y,d2z, &
               d3x,d3y,d3z, &
               x,y,z, r, dd
    real(8),dimension(1:250) :: pot
    integer :: n1,n2,n3, i1, i2, i3
    !===================

    !I/O   (if not redefined here, the some labels are lost...)
    I_INP=10
    I_TOP=11
    I_NDX=12
    O_CUB=20

    !Defaults
    dd = 1.0

    ! 0. GET COMMAND LINE ARGUMENTS AND OPEN FILES
    call parse_input(inpfile,filetype,topfile,ndxfile,resname,outcube,&
                     plain,dd,verbose)

    ! 1. READ DATA
    open(I_INP,file=inpfile,iostat=IOstatus,status="old")
    if (IOstatus /= 0) call alert_msg("fatal","Unable to open "//trim(adjustl(inpfile)))
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
                    residue(i)%atom(j)%element = molec%atom(k)%element
                    molec%atom(k) = residue(i)%atom(j)
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
    
    ! 3. Take only relevant molecules (those in L layer)
    j=0
    k=0
    do i=1,molec%natoms
        !Change units
        molec%atom(i)%x = molec%atom(i)%x / BOHRtoANGS
        molec%atom(i)%y = molec%atom(i)%y / BOHRtoANGS
        molec%atom(i)%z = molec%atom(i)%z / BOHRtoANGS

        !If pointcharges, treat L layer as point charges
        if (molec%atom(i)%chain == "L") then
            j=j+1
            molec_L%atom(j) = molec%atom(i)
        else
            k=k+1
            molec_H%atom(k) = molec%atom(i)
            xmax = max(xmax,molec%atom(i)%x)
            ymax = max(ymax,molec%atom(i)%y)
            zmax = max(zmax,molec%atom(i)%z)
            x0 = min(x0,molec%atom(i)%x)
            y0 = min(y0,molec%atom(i)%y)
            z0 = min(z0,molec%atom(i)%z)
        endif
    enddo
    molec_L%natoms = j
    molec_H%natoms = k

    print*, "Will make a box with limits"
    print'(A,3F8.2)', "  X:", x0, xmax, xmax-x0
    print'(A,3F8.2)', "  Y:", y0, ymax, ymax-y0
    print'(A,3F8.2)', "  Z:", z0, zmax, zmax-z0

    !Set the limits of the box using a buffer of 3AA
    x0 = x0 - 3.d0/BOHRtoANGS
    xmax = xmax + 3.d0/BOHRtoANGS
    y0 = y0 - 3.d0/BOHRtoANGS
    ymax = ymax + 3.d0/BOHRtoANGS
    z0 = z0 - 3.d0/BOHRtoANGS
    zmax = zmax + 3.d0/BOHRtoANGS

    !We need at nums
    call atname2element(molec_H)
    call element2AtNum(molec_H)

    !For the moment we use canonical vectors
    d1x = dd
    d1y = 0.d0
    d1z = 0.d0
    d2x = 0.d0
    d2y = dd
    d2z = 0.d0
    d3x = 0.d0
    d3y = 0.d0
    d3z = dd

    n1 = int((xmax-x0)/d1x)+1
    n2 = int((ymax-y0)/d2y)+1
    n3 = int((zmax-z0)/d3z)+1

    !Print file
    open(O_CUB,file=outcube)
    write(O_CUB,*) "File generated from "//trim(adjustl(inpfile))
    write(O_CUB,*) "Electronic potential generated by low layer atoms"
    write(O_CUB,'(I5,4F12.6)') molec_H%natoms,x0,y0,z0
    write(O_CUB,'(I5,4F12.6)') n1, d1x,d1y,d1z
    write(O_CUB,'(I5,4F12.6)') n2, d2x,d2y,d2z
    write(O_CUB,'(I5,4F12.6)') n3, d3x,d3y,d3z
    do i=1,molec_H%natoms
        write(O_CUB,'(I5,4F12.6)') molec_H%atom(i)%AtNum, float(molec_H%atom(i)%AtNum), &
                                   molec_H%atom(i)%x, molec_H%atom(i)%y, molec_H%atom(i)%z
    enddo     
    do i1=1,n1
    do i2=1,n2
    do i3=1,n3
        x = x0 + (i1-1)*d1x+(i2-1)*d1y+(i3-1)*d1z
        y = y0 + (i1-1)*d2x+(i2-1)*d2y+(i3-1)*d2z
        z = z0 + (i1-1)*d3x+(i2-1)*d3y+(i3-1)*d3z
        pot(i3)=0.d0
        do i=1,molec_L%natoms
            r = dsqrt((x-molec_L%atom(i)%x)**2+&
                      (y-molec_L%atom(i)%y)**2+&
                      (z-molec_L%atom(i)%z)**2)
            pot(i3) = pot(i3)+molec_L%atom(i)%q/r
        enddo
    enddo
    write(O_CUB,'(6F12.6)') pot(1:n3)
    enddo
    enddo
    close(O_CUB)


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

    subroutine parse_input(inpfile,filetype,topfile,ndxfile,resname,outfile,&
                           plain,delta,verbose)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,topfile,filetype,ndxfile,outfile,&
                                          resname
        real(8),intent(inout) :: delta
        logical,intent(inout) :: verbose, plain
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

                case ("-res")
                    call getarg(i+1, resname)
                    argument_retrieved=.true.

                case ("-delta")
                    call getarg(i+1, arg)
                    argument_retrieved=.true.
                    read(arg,*) delta

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
        write(0,*) '-o              ', trim(adjustl(outfile))
        write(0,*) '-res            ', trim(adjustl(resname))
        
        write(0,*) '-delta          ', delta
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
             call read_gro(unt,molec)
            case("g96")
             call read_g96(unt,molec)
            case("pdb")
             call read_pdb_new(unt,molec)
            case("log")
             call parse_summary(unt,molec,props,"struct_only")
!             case("fchk")
!              call read_fchk_geom(unt,molec)
            case default
             call alert_msg("fatal","Trying to guess, but file type but not known: "//adjustl(trim(filetype))&
                        //". Try forcing the filetype with -ft flag (available: gro, g96, pdb)")
        end select

        else
        ! Predefined filetypes
        select case (adjustl(filetype))
            case("gro")
             call read_gro(unt,molec)
            case("g96")
             call read_g96(unt,molec)
            case("pdb")
             call read_pdb_new(unt,molec)
            case("log")
             call parse_summary(unt,molec,props,"struct_only")
!             case("fchk")
!              call read_fchk_geom(unt,molec)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select
        endif

        return


    end subroutine generic_strfile_read

end program gen_oniom

