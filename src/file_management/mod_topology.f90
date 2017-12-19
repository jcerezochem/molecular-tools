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

    !*****************
    !   MODULE LOAD
    !*****************
    !============================================
    !   Generic
    !============================================
    use io
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
    !============================================
    !  Structure-related modules
    !============================================
    use molecular_structure
    use atomic_geom
    use symmetry
    !============================================
    !  Internal thingies
    !============================================
    use zmat_manage 

    implicit none

    !====================== 
    !System variables
    type(str_resmol) :: molec 
    !Residues are different molecule types (there might be several of each kind in the system).
    ! E.g. Solvent and solute are only 2 type of molecules 
    type(str_resmol),dimension(:),allocatable :: residue
    !
    character(len=5),dimension(:),allocatable :: molname
    integer,dimension(:),allocatable :: molnum
    integer :: nmol
    integer,dimension(50000) :: molmap, resdone=0
    integer,dimension(10000) :: frz
    !
    integer :: ires, imap, nres
    !
    logical :: pointcharges=.false., &
               plain=.false.
    integer :: natoms = -1
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
               O_TOP=20
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
    O_TOP=20
    
    call parse_input(topfile,natoms)
    
    if (natoms<0) natoms=500

    open(I_TOP,file=topfile,iostat=IOstatus,status="old")
    if (IOstatus /= 0) call alert_msg("fatal","Unable to open "//trim(adjustl(topfile)))
    allocate(molnum(1:10))
    call read_top_nres(I_TOP,nres,molnum)
    close(I_TOP)
    ! Allocate stuff to read top
    allocate(molname(1:molec%natoms))
    allocate(residue(1:2))
    do i=1,nres
        call allocate_atoms(residue(i),molnum(i)+1000)
    enddo
    deallocate(molnum)
    open(I_TOP,file=topfile,iostat=IOstatus,status="old")
    call read_top(I_TOP,residue,molname,nmol)
!     call read_top(I_TOP,(/molec%atom(1:molnum(1)),molec%atom(molnum(1)+1:molnum(1)+molnum(2))/),molname,nmol)
    close(I_TOP)
    !========================0
    !should THIS be in read_top?
    !Identify different molecules (like that if splitted in [ molecules ] will not work.
    ! THE MANAGEMENT OF MOLECULES SHOULD BE WITH A MOLID NOT WITH MOLNAME THAT MAY BE TOO SHORT (len=5 HERE)
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

    print*, "Writing new topology for "//trim(adjustl(residue(1)%name))
    
    open(O_TOP,file="new_topol.top",status="replace")
    call write_top(O_TOP,residue(1))
    close(O_TOP)
    
    !Review notes and errors
    if (n_notes /= 0) &
     print'(/,A,I0,A,/)', "There were ", n_notes, " NOTES in this run"
    if (n_errors /= 0) &
     print'(/,A,I0,A,/)', "There were ", n_errors, " WARNINGS in this run"

    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(topfile,natoms)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: topfile
        integer,intent(inout)          :: natoms
        ! Local
        character(len=500) :: input_command
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        ! iargc type must be specified with implicit none (strict compilation)
        integer :: iargc

        call getarg(0,input_command)
        !Get input flags
        do i=1,iargc()
            call getarg(i,arg)
            input_command = trim(adjustl(input_command))//" "//trim(adjustl(arg))
        enddo

        argument_retrieved=.false.
        do i=1,iargc()
            if (argument_retrieved) then
                argument_retrieved=.false.
                cycle
            endif
            call getarg(i, arg) 
            select case (adjustl(arg))
                case ("-p") 
                    call getarg(i+1, topfile)
                    argument_retrieved=.true.

                case ("-alloc-atm")
                    call getarg(i+1,arg)
                    read(arg,*) natoms
                    argument_retrieved=.true.

                ! Control verbosity
                case ("-quiet")
                    verbose=0
                case ("-concise")
                    verbose=1
                case ("-v")
                    verbose=2
                case ("-vv")
                    verbose=3

                case ("-h")
                    need_help=.true.

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 
          

       !Print options (to stdout)    
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '              M O D _ T O P O L O G Y       '     
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag         Description                   Value'
        write(6,'(A)')   '-------------------------------------------------------------------'
        write(6,*) '-p           Topology file                 ', trim(adjustl(topfile))
        write(6,'(X,A,I0)') &
                   '-alloc-atm   Atoms to allocate(-1:default) ', natoms
        write(6,*) '-h           Show this help and quit       ',  need_help
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(A)') 'Input command:'
        write(6,'(A)') trim(adjustl(input_command))   
        write(6,'(A)') '-------------------------------------------------------------------'
        write(6,'(X,A,I0)') &
                       'Verbose level:  ', verbose        
        write(6,'(A)') '-------------------------------------------------------------------'
        if (need_help) call alert_msg("fatal", 'There is no manual (for the moment)' )

        return
    end subroutine parse_input


end program gen_oniom

