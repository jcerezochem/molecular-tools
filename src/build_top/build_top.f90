program build_top

    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to generate a topology file from gaussian log files.
    ! Both coordinates and atomic charges are read. Atom types are
    ! guessed from the (also guessed) connectivity [but it works...]
    !
    ! Compilation instructions (for mymake script):
    !make$ gfortran ../modules/alerts.f90 ../modules/structure_types_v2_ALLOC.f90 ../modules/allocation_mod.f90 ../modules/line_preprocess.f90 ../modules/ff_build_module_v2.f90 ../modules/gro_manage_v2.f90 ../modules/gaussian_manage_v2.f90 ../modules/pdb_manage_v2.f90 ../modules/molpro_manage_v2.f90 BUILD_TOP_v2.f90 -o BUILD_TOP_v2.exe -cpp
    !
    ! Change log:
    ! Feb 2012: ff_build and related subroutine now work on residue not on system (beta)
    ! May 2012: Uses version 2 of molecular tools (v2)
    !           Implement allocation_mod
    ! Jun 2012 (v2) Some features are broken (un-reverted tests?) -- SOLVED
    ! Version 4: Uses v4 modules (allocation disabled)
    !============================================================================    

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
    use ff_build
    use atomic_geom
    use symmetry
    !============================================
    !  Internal thingies
    !============================================
    use zmat_manage 

    implicit none

!     common /GLOBAL_OPTS/ do_refine_charges, I_DB, rename_atoms

    !====================== 
    !Options 
    logical :: call_vmd=.false.,          & 
               do_refine_charges=.false., & !.true.
               united_atom=.false.,       &
               rename_atoms=.false.
    !======================

    !====================== 
    !System variables
    type(str_resmol) :: molec !, molecUA
    type(str_resmol),dimension(10) :: residue, residueUA
    character(len=6) :: resname="UNK"
    character(len=8) :: q_type="mulliken"
    real :: charge
    integer :: nat_alloc=-1
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: id, ih, iat, iiat
    character(len=1) :: cnull
    character(len=4) :: atname
    character(len=16) :: dummy_char
    integer :: dummy_int
    !====================== 
    
    !=================
    !hbd stuff
    integer :: hdb_type
    integer :: numH, numNH
    integer,dimension(5) :: hhs,nhs !(hydrogens and non-hydrogens)
    logical :: reorder_hyd
    !=================

    !=============
    !Counters
    integer :: i,j,k,l, ii,jj,kk,ll
    !=============

    !================
    !I/O stuff 
    !units
    integer :: I_INP=10,  &
               I_MUL=11,  & 
               I_DB =12,  & !zero means no external DB
               S_SPDB=30, &
               O_GRO=20,  &
               O_TOP=21,  &
               O_HDB=22
    !files
    character(len=10) :: filetype="guess"
    character(len=200):: inpfile="input.log",          &
                         outgro="out.gro",             &
                         outtop="out.top",             &
                         outhdb,                       &
                         database="charmm" ! hybrid or <filename>              
    !status
    integer :: IOstatus
    !===================


    ! 0. GET COMMAND LINE ARGUMENTS
    call parse_input(inpfile,filetype,outgro,outtop,q_type,resname,database,united_atom,nat_alloc, &
                    reorder_hyd)

    ! Allocate atoms
    if (nat_alloc<0) then
        call allocate_atoms(molec)
        do i=1,10
            call allocate_atoms(residue(i))
            call allocate_atoms(residueUA(i))
        enddo
    else 
        call allocate_atoms(molec,nat_alloc)
        do i=1,10
            call allocate_atoms(residue(i),nat_alloc)
            call allocate_atoms(residueUA(i),nat_alloc)
        enddo
    endif
    
    ! 1. READ DATA
    ! ---------------------------------
    open(I_INP,file=inpfile,status='old',iostat=IOstatus)
    if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(inpfile)) )

    !Load atomic coordinates on "molecule"
    if (adjustl(filetype) == "guess") then
        ! Guess file type
        call split_line(inpfile,".",cnull,filetype)
    endif
    
    call generic_strmol_reader(I_INP,filetype,molec)

    close(I_INP)


    ! 2. ASSOCIATE ATOM TYPES AND BONDED TERMS (if not disabled)
    ! -----------------------------------------
    if ( adjustl(outtop) /= "none" )  then
        !Attypes DataBase
        if ( adjustl(database) == "charmm" ) then
            I_DB = 0
        elseif ( adjustl(database) == "hybrid" ) then
            I_DB = -1
        else
            open(I_DB, file=database, status="old", IOSTAT=IOstatus)
            if (IOstatus /= 0) call alert_msg( "fatal","Unable to open "//trim(adjustl(database)) )
        endif

        call sist2res(molec,residue)

        !Only one residue supported (for the momenfor the moment)
        if (molec%nres /= 1) call alert_msg("fatal","Only structures with one residue are currently supported")
        call guess_connect(residue(1))
        call build_ff(residue(1),I_DB)
        call gen_bonded(residue(1))
        
        if (reorder_hyd) then
            ! Rewrite structure with reordered hydrogens
            iat  = 0
            iiat = 0
            do i=1,residue(1)%natoms
                if (residue(1)%atom(i)%element == 'H') cycle
                numH=0
                do j=1,residue(1)%atom(i)%nbonds
                    jj=residue(1)%atom(i)%connect(j)
                    if (residue(1)%atom(jj)%element == 'H') then
                        numH = numH + 1
                    endif
                enddo
                ih   = 0
                iat  = iat  + 1
                iiat = iiat + 1
                molec%atom(iat) = residue(1)%atom(i)
                do j=1,residue(1)%atom(i)%nbonds
                    jj=residue(1)%atom(i)%connect(j)
                    if (residue(1)%atom(jj)%element == 'H') then
                        iat = iat + 1
                        ih = ih + 1
                        molec%atom(iat) = residue(1)%atom(jj)
                        if (numH>1) then
                            write(atname,'(A,I0,I0)') "H", iiat, ih
                        else
                            write(atname,'(A,I0,I0)') "H", iiat
                        endif
                        molec%atom(iat)%name = atname
                    endif
                enddo
            enddo
            call sist2res(molec,residue)
            call guess_connect(residue(1))
            call build_ff(residue(1),I_DB)
            call gen_bonded(residue(1))
        endif

        if (I_DB>0) close(I_DB)
    endif

    !United atoms model:
    if (united_atom) then
        call aa2ua(residue(1),residueUA(1))
        residue(1) = residueUA(1)
        !Recalculate connections 
        call guess_connect(residue(1))
        call gen_bonded(residue(1))
        call res2sist(residue,molec%nres,molec%BoxX,molec%BoxY,molec%BoxZ,molec)
    endif


    ! 3. WRITE COORDINATES (.gro) AND TOPOLOGY (.top)
    ! -------------------------------------------------  
    
    !GROFILE:
    open(O_GRO,file=outgro,status='replace')
    call write_gro(O_GRO,molec)
    close(O_GRO)

    !TOPFILE (if not disabled)
    if ( trim(adjustl(outtop)) /= "none" ) then
        !The topology builder uses residue as input:
!         call sist2res(molec,residue)
        open(O_TOP,file=outtop,status='replace')
        call split_line(outtop,".",cnull,filetype)
        select case (adjustl(filetype))
            case("itp")
             call write_top(O_TOP,residue(1))
            case("top")
             call write_top(O_TOP,residue(1))
            case("rtp")
             call write_rtp(O_TOP,residue(1))
            case default
             !By default a top format is use
             call write_top(O_TOP,residue(1))
             call alert_msg("note","Topology file written in *.top format")
        end select
        close(O_TOP)
    endif
    
    if (reorder_hyd) then
        ! Rewrite structure with reordered hydrogens
        iat  = 0
        iiat = 0
        do i=1,residue(1)%natoms
            if (residue(1)%atom(i)%element == 'H') cycle
            numH=0
            iiat = iiat + 1
            do j=1,residue(1)%atom(i)%nbonds
                jj=residue(1)%atom(i)%connect(j)
                if (residue(1)%atom(jj)%element == 'H') then
                    numH = numH + 1
                    write(atname,'(A,I0,I0)') "H", iiat
                    residue(1)%atom(jj)%name = atname
                endif
            enddo
        enddo
    endif
    
    ! Get hcb
    call split_line_back(outtop,'.',outhdb,cnull)
    outhdb=trim(adjustl(outhdb))//'.hdb'
    open(O_HDB,file=outhdb)
    kk=0
    do i=1,residue(1)%natoms
        numH=0
        do j=1,residue(1)%atom(i)%nbonds
            jj=residue(1)%atom(i)%connect(j)
            if (residue(1)%atom(jj)%element == 'H') numH=numH+1
        enddo
        if (numH==0) cycle
        kk=kk+1
    enddo
    write(O_HDB,'(A6,4X,I5)') resname, kk
    iiat = 0
    do i=1,residue(1)%natoms
        numH=0
        numNH=0
        do j=1,residue(1)%atom(i)%nbonds
            jj=residue(1)%atom(i)%connect(j)
            if (residue(1)%atom(jj)%element == 'H') then
                numH=numH+1
                hhs(numH)=jj
            else
                numNH=numNH+1
                nhs(numNH)=jj
            endif
        enddo
        if (numH==0) cycle
        call fftype_db(residue(1)%atom(i),numH, &
                       dummy_char,cnull,        &
                       -2, i)
        read(dummy_char,*) hdb_type
        
        ll=0
        select case(hdb_type)
            case(1)
                !j-i-k (i-H)
                jj=nhs(1)
                kk=nhs(2)
            case(2)
                !k-j-i (i-H)
                jj=nhs(1)
                do k=1,residue(1)%atom(jj)%nbonds
                    kk = residue(1)%atom(jj)%connect(k)
                    if (kk == i) cycle
                    if (residue(1)%atom(kk)%element /= 'H') exit 
                enddo
            case(3)
                !k-j-i (i-H1,H2)
                jj=nhs(1)
                do k=1,residue(1)%atom(jj)%nbonds
                    kk = residue(1)%atom(jj)%connect(k)
                    if (kk == i) cycle
                    if (residue(1)%atom(kk)%element /= 'H') exit 
                enddo
            case(4)
                !k-j-i (i-H1,H2,[H3])
                jj=nhs(1)
                do k=1,residue(1)%atom(jj)%nbonds
                    kk = residue(1)%atom(jj)%connect(k)
                    if (kk == i) cycle
                    if (residue(1)%atom(kk)%element /= 'H') exit 
                enddo
            case(5)
                !  l
                !  |
                !k-i-j (i-H)
                jj=nhs(1)
                kk=nhs(2)
                ll=nhs(3)
            case(6)
                !j-i-k (i-H1,H2)
                jj=nhs(1)
                kk=nhs(2)
            case default
                write(dummy_char,'(i0)') hdb_type
                call alert_msg('warning','Unkown hdb type: '//trim(adjustl(dummy_char)))
                jj=0
                kk=0
            end select
            
            if (ll==0) then
                write(O_HDB,'(10X,I5,X,I3,X,5(A5,X))') numH,hdb_type,residue(1)%atom(hhs(1))%name,&
                                                                     residue(1)%atom(i )%name,&
                                                                     residue(1)%atom(jj)%name,&
                                                                     residue(1)%atom(kk)%name
            else
                write(O_HDB,'(10X,I5,X,I3,X,5(A5,X))') numH,hdb_type,residue(1)%atom(hhs(1))%name,&
                                                                     residue(1)%atom(i )%name,&
                                                                     residue(1)%atom(jj)%name,&
                                                                     residue(1)%atom(kk)%name,&
                                                                     residue(1)%atom(ll)%name
            endif
                        
    enddo

    ! 9999. CHECK ERROR/NOTES
    ! -------------------------------------------------
    if (n_notes > 0) then 
        write(dummy_char,*) n_notes
        write(6,'(/,A,/)') "There were "//trim(adjustl(dummy_char))//" note(s) in this run"
    endif
    if (n_errors > 0) then 
        write(dummy_char,*) n_errors
        write(6,'(/,A,/)') "There were "//trim(adjustl(dummy_char))//" warning(s) in this run"
        call alert_msg("fatal", "Files generated, but with too many warnings")
    endif

    ! Summary output files
    if ( trim(adjustl(outtop)) /= "none" ) then
        write(6,'(/,A,/)') "3 output files generated: "//trim(adjustl(outgro))//",  "//trim(adjustl(outtop))//&
                                                  ",  "//trim(adjustl(outhdb))
    else
        write(6,'(/,A,/)') "1 output file generated: "//trim(adjustl(outgro))
    endif


    stop


    !==============================================
    contains
    !=============================================

    subroutine parse_input(inpfile,filetype,outfile,outfile2,q_type,resname,database,united_atom,nat_alloc,&
                           reorder_hyd)
    !==================================================
    ! My input parser (gromacs style)
    !==================================================
        implicit none

        character(len=*),intent(inout) :: inpfile,filetype,outfile,outfile2,resname,q_type,database
        logical,intent(inout) :: united_atom, reorder_hyd
        integer :: nat_alloc
        ! Local
        character(len=500) :: input_command
        logical :: argument_retrieved,  &
                   need_help = .false.
        integer:: i
        character(len=200) :: arg
        
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
                    call getarg(i+1, outfile2)
                    argument_retrieved=.true.

                case ("-chrg") 
                    call getarg(i+1, q_type)
                    argument_retrieved=.true.

                case ("-res") 
                    call getarg(i+1, resname)
                    argument_retrieved=.true.
        
                case ("-db")
                    call getarg(i+1, database)
                    argument_retrieved=.true.

                case ("-ua")
                     united_atom = .true.
                     
                case ("-reorder-hyd")
                     reorder_hyd = .true.
                     
                case ("-alloc-atm")
                    call getarg(i+1,arg)
                    read(arg,*) nat_alloc
                    argument_retrieved=.true.
                     

                case ("-h")
                    need_help=.true.
                    
                    
                ! Control verbosity
                case ("-quiet")
                    verbose=0
                case ("-concise")
                    verbose=1
                case ("-v")
                    verbose=2
                case ("-vv")
                    verbose=3
                    

                case default
                    call alert_msg("fatal","Unkown command line argument: "//adjustl(arg))
            end select
        enddo 

        ! Some checks on the input
        !----------------------------

        !-p disabling accepts different inputs
        if ( trim(adjustl(outtop)) == "none" .or.   &
             trim(adjustl(outtop)) == "NONE" .or.   &
             trim(adjustl(outtop)) == "no"   .or.   &
             trim(adjustl(outtop)) == "NO"        ) then 
            outtop = "none"
        endif

        !-chrg options
        if ( q_type /= "mulliken" .and.  &
             q_type /= "elec_pot" .and.   &
             q_type /= "all_zero" ) then
            call alert_msg("fatal", 'Options for charges (-chrg) are "mulliken", "elec_pot" and "all_zero", select one of these')
        endif
! 
!         if ( q_type == "all_zero" ) do_refine_charges=.false.
          

       !Print options (to stdout)    
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '              B U I L D _ T O P '    
        write(6,'(/,A)') '     A connectivity based topology builder '        
        write(6,'(A)')   '                for GROMACS'     
        call print_version()
        write(6,'(/,A)') '========================================================'
        write(6,'(/,A)') '-------------------------------------------------------------------'
        write(6,'(A)')   ' Flag         Description                   Value'
        write(6,'(A)')   '-------------------------------------------------------------------'


        write(0,*) '-h             ',  need_help
        
        write(6,*) '-f           Input file (structure)        ', trim(adjustl(inpfile))
        write(6,*) '-ft          \_ FileType                   ', trim(adjustl(filetype))
        write(6,*) '-o           Output (com) file             ', trim(adjustl(outgro))
        write(6,*) '-p           Topology file                 ', trim(adjustl(outtop))
        write(6,*) '-chrg        Charge type (mulliken...)     ', trim(adjustl(q_type))
        
        write(6,*) '-res         Reside name                   ', trim(adjustl(resname))
        write(6,*) '-db          Database to set atom types    ', trim(adjustl(database))
        write(6,'(X,A,I0)') &
                   '-alloc-atm   Atoms to allocate(-1:default) ', nat_alloc
        write(6,*) '-reorder-hyd Reorder hydrogens            ',  reorder_hyd
        write(6,*) '-h           Show this help and quit      ',  need_help
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

end program build_top

