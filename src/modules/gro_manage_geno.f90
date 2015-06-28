module gro_manage

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  Basic I/O subroutines to read/write GRO files.
    !    subroutine read_gro(unt,system)
    !    subroutine read_gro_def(unt,system)
    !    subroutine write_gro(unt,system)
    !    subroutine write_gro_def(unt,system)
    !    subroutine write_top_new(unt,residue)
    !
    ! Dependences
    !  Uses "structure_types" module
    !
    ! Notes
    !  Since variables molec, residue... are allocatable, they should
    !  be passed to the sr even if they are output to get the allocation.
    !_____________________________________________________________________________
    !**Default units need to be taken among all subroutines of molecular tools 
    !  unless otherwise required (e.g. atomic units in QM code). Tha will be
    !  -Length: Amstrong
    !
    !  NEW SUBROUTINES (_def) SHOULD BECOME THE DEFAULT (AND THE ONLY) CALL
    !==============================================================================
    !
    ! History
    !  2011/12/20  Added "renum" subroutine
    !  2013/01/23  Added "read_top" (still in development)
    !  2013/07/15  Improved read_top: the GMXLIB folder is now readble to look for its files (release v3)
    !  2014/02/11  Uses structure_types v4 (changed str_system by str_resmol)
    !              _def version become default (prev are removed): i.e.: \AA are used as ref. length units
    !  2014/10/23  Fixed bug in read_g96. Now it also reads properly when there is not title
    !              Fixed bug in read_g96. Now also gets the structure from POSITIONRED sections
    !-----------------------------------------------------

    !Common modules for all subroutines
    use structure_types
    use constants
    implicit none

    contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine read_gro(unt,system)

       !Input is in \AA. Conversion to nm is performed here!

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, natoms


        read(unt,'(A)') system%title
        read(unt,*) natoms

        system%natoms=natoms

        do i=1,natoms   

            read(unt,100) system%atom(i)%resseq,    &
                          system%atom(i)%resname,   &
                          system%atom(i)%name,      &
                          !serial
                          system%atom(i)%x,         &
                          system%atom(i)%y,         &
                          system%atom(i)%z
!             !Atomnames are "adjustl"ed - TODO?
!             system%atom(i)%name = adjustl(system%atom(i)%name)

        enddo
        system%atom(1:natoms)%x=system%atom(1:natoms)%x*10.d0
        system%atom(1:natoms)%y=system%atom(1:natoms)%y*10.d0
        system%atom(1:natoms)%z=system%atom(1:natoms)%z*10.d0

        read(unt,101) system%boxX, system%boxY, system%boxZ

        return

    !gro file format
100 format(i5,2a5,5X,3f8.3,3f8.4) !Serial is not read
101 format(3f10.5)

    end subroutine read_gro

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine write_gro(unt,system)

       !Output is in \AA. Conversion from nm is performed here!

        integer,intent(in)::unt
        type(str_resmol),intent(in)::system

        !local
        integer::i

        write(unt,'(A)') 'TITLE  File created with gro_manage module'
        write(unt,'(I5)') system%natoms

        do i=1,system%natoms        
        
            write(unt,200) system%atom(i)%resseq,    &
                           system%atom(i)%resname,   &
                           system%atom(i)%name,      &
                           i,                        & !This is the serial number
                           system%atom(i)%x/10.d0,   &
                           system%atom(i)%y/10.d0,   &
                           system%atom(i)%z/10.d0

        enddo

        write(unt,201) system%boxX, system%boxY, system%boxZ

        return
    
    !gro file format
200 format(i5,2a5,i5,3f8.3,3f8.4)
201 format(3f10.5)


    end subroutine write_gro

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine read_g96(unt,system)

       !Reads GROMOS-96 structure file (input/output in ANGS but the file in NM)

        integer,intent(in)::unt
        type(str_resmol_light),intent(inout)::system

        !local
        integer::i, natoms, ii, ios
        character(len=260) :: line
        !The use of get_structure=.false. fails probably due to
        ! memory issues (change without being directly accesed)
        logical :: not_get_structure=.true.

        !The file is organized in sections. Each begining with a given name
        do

            read(unt,*,iostat=ios) line
            if (ios /= 0) exit

            if (adjustl(line) == "TITLE") then
                read(unt,'(A)') system%title
                if (adjustl(system%title) == "END") then
                     system%title=""
                else
                    read(unt,*) line !END
                endif

            elseif (adjustl(line) == "POSITION" ) then
                not_get_structure=.false.
                i=0
                do 
                     read(unt,'(A)') line
                     if (adjustl(line) == "END" ) exit
                     i=i+1
                     read(line,*) system%atom(i)%resseq, &
                                  system%atom(i)%resname,&
                                  system%atom(i)%name,   &
                                  ii,                    & !serial
                                  system%atom(i)%x,      &
                                  system%atom(i)%y,      &
                                  system%atom(i)%z
                     system%atom(i)%x = system%atom(i)%x*10.d0
                     system%atom(i)%y = system%atom(i)%y*10.d0
                     system%atom(i)%z = system%atom(i)%z*10.d0
                enddo
                system%natoms = i
            elseif (adjustl(line) == "POSITIONRED" ) then
            !this section only has info about coordinates (no atom info!)
            write(0,*) "NOTE: No atomic info in g96!"
                not_get_structure=.false.
                i=0
                do 
                     read(unt,'(A)') line
                     if (adjustl(line) == "END" ) exit
                     i=i+1
                     read(line,*) system%atom(i)%x,      &
                                  system%atom(i)%y,      &
                                  system%atom(i)%z
                     system%atom(i)%x = system%atom(i)%x*10.d0
                     system%atom(i)%y = system%atom(i)%y*10.d0
                     system%atom(i)%z = system%atom(i)%z*10.d0
                     !Missing info is set to "default" values
                     system%atom(i)%resseq=1
                     system%atom(i)%resname="UNK"
                     system%atom(i)%name="X"
                enddo
                system%natoms = i
            elseif (adjustl(line) == "BOX") then
                read(unt,*) system%boxX, system%boxY, system%boxZ
                read(unt,*) line ! END
            else
                cycle
            endif

        enddo

        if (not_get_structure) then
            write(0,*) "ERROR: No structure read in g96 file"
            stop
        endif

        return

    end subroutine read_g96

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine write_g96(unt,system)

       !Writes GROMOS-96 structure file (input/output in ANGS but the file in NM)

        integer,intent(in)::unt
        type(str_resmol),intent(in)::system

        !local
        integer::i

        write(unt,'(A)') "TITLE"
        write(unt,'(A)') "Generated with gro_manage module: "//trim(adjustl(system%title))
        write(unt,'(A)') "END"

        write(unt,'(A)') "POSITION"
        do i=1,system%natoms     
        
            write(unt,300) system%atom(i)%resseq,    &
                           system%atom(i)%resname,   &
                           system%atom(i)%name,      &
                           i,                        & !This is the serial number
                           system%atom(i)%x/10.d0,   &
                           system%atom(i)%y/10.d0,   &
                           system%atom(i)%z/10.d0

        enddo
        write(unt,'(A)') "END"

        write(unt,'(A)') "BOX"
        write(unt,301) system%boxX, system%boxY, system%boxZ
        write(unt,'(A)') "END"

        return

300 format(i5,x,2(a5,x),x,i5,3f15.9,3f8.4)
301 format(3f15.9)

    end subroutine write_g96


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine write_top_new(unt,residue)
       
       !======================================================
       ! Subrutina para el la estructura de la topología (itp)
       ! de un residuo. Usa un resido (solo uno) no el sistema
       !=======================================================

        implicit none

        integer,intent(in)::unt
        type(str_resmol),intent(in)::residue
        
        !local
        integer :: nrexcl, &
                   res_number=1 !to be generalized
        real :: total_charge
        integer,dimension(3) :: date
        !Counter
        integer::i, j

        character(len=50) :: dummy_char

        !Header
!         call idate(date)
        write(dummy_char,'(i2,a1,i2,a1,i4)') date(1),"/",date(2),"/",date(3)
        write(unt,'(A,/,A)') '; Topology (itp style) generated by "build_top"', &
                             '; on '//adjustl(dummy_char)

        !1. Section: [ moleculetype ]
        write(unt,'(/,A)') '[ moleculetype ]'
        write(unt,'(A)') '; Name            nrexcl'
        if (adjustl(residue%name) == "SOL") then
            nrexcl=2
        else
            nrexcl=3
        endif
        write(unt,'(A,3X,I1)') trim(adjustl(residue%name)), nrexcl

        !2. Section: [ atoms ]
        write(unt,'(/,A)') '[ atoms ]'
        write(unt,'(A)')   ';   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB'
        total_charge=0.
        do i=1,residue%natoms
            total_charge=total_charge+residue%atom(i)%q
            write(unt,300) i,                        &
                           residue%atom(i)%fftype,   &
                           res_number            ,   &
                           residue%name          ,   &
                           residue%atom(i)%name  ,   &
                           i                     ,   & !this is CHARMM specific (no charge groups)
                           residue%atom(i)%q     ,   &
                           residue%atom(i)%mass  ,   &
                           ";"                   ,   &
                           total_charge             
        enddo

        !3. Section: [ bonds ]
        write(unt,'(/,A)') '[ bonds ]'
        write(unt,'(A)')   ';  ai    aj funct            c0            c1            c2            c3'
        do i=1,residue%natoms-1
            do j=1,residue%atom(i)%nbonds
                if ( residue%atom(i)%connect(j) > i ) &
                  write(unt,301) i, residue%atom(i)%connect(j), 1
            enddo
        enddo

        !3. Section: [ pairs ]
        write(unt,'(/,A)') '[ pairs ]'
        write(unt,'(A)')   ';  ai    aj funct            c0            c1            c2            c3'
        do i=1,residue%geom%npairs
            write(unt,301) residue%geom%pair(i,1), residue%geom%pair(i,2), 1
        enddo

        !4. Section: [ angles ]
        write(unt,'(/,A)') '[ angles ]'
        write(unt,'(A)')   ';  ai    aj    ak funct            c0            c1            c2            c3'
        do i=1,residue%geom%nangles
!                   if (adjustl(residue%atom(i)%fftype) == "CA" .and. &
!                       adjustl(residue%atom(residue%atom(i)%angle(j,1))%fftype) == "CA".and. & 
!                       adjustl(residue%atom(residue%atom(i)%angle(j,2))%fftype) == "CA") &
                write(unt,302) residue%geom%angle(i,1), & !ai
                               residue%geom%angle(i,2), & !aj
                               residue%geom%angle(i,3), & !ak
                               5                             !Funct (Charmm)
        enddo

        !5. Section: [ dihedrals ]
        write(unt,'(/,A)') '[ dihedrals ]'
        write(unt,'(A)')   ';  ai    aj    ak    al funct            c0            c1            c2            c3'
        do i=1,residue%geom%ndihed
!                   if (adjustl(residue%atom(i)%fftype) == "CA" .and. &
!                       adjustl(residue%atom(residue%atom(i)%dihed(j,1))%fftype) == "CA".and. & 
!                       adjustl(residue%atom(residue%atom(i)%dihed(j,2))%fftype) == "CA".and. &
!                       adjustl(residue%atom(residue%atom(i)%dihed(j,3))%fftype) == "CA" ) &
                !Separate rotated bonds
                if ( residue%geom%dihed(i,2) /= residue%geom%dihed(i-1,2) .or. &
                     residue%geom%dihed(i,3) /= residue%geom%dihed(i-1,3) ) then
                    write(unt,'(A8,2I3)') "; Bond: ", residue%geom%dihed(i,2),residue%geom%dihed(i,3)
                endif
                write(unt,303) residue%geom%dihed(i,1), & !ai
                               residue%geom%dihed(i,2), & !aj
                               residue%geom%dihed(i,3), & !ak
                               residue%geom%dihed(i,4), & !al
                               9                          !Funct (Charmm)
        enddo

        write(unt,'(/,A)') '; Planar/tetrahedric impropers not generated'

        !6. Position restrain option
        write(unt,'(/,A)') '; Include Position restraint file'
        write(unt,'(A)')   '#ifdef POSRES'
        write(unt,'(A)')   '#include "posre.itp"'
        write(unt,'(A,/)') '#endif'


        return
    
    !gro file format
300 format(i5,3x,a5,x,i5,x,a5,x,a5,x,i5,x,f8.3,x,f8.3,x,a1,x,f8.3)
301 format(i5,3x,i5,3x,i2)
399 format(a1,i5,3x,i5,3x,i2)
302 format(i5,3x,i5,3x,i5,3x,i2)
303 format(i5,3x,i5,3x,i5,3x,i5,3x,i2)
!201 format(3f10.5)

    end subroutine write_top_new

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine write_rtp(unt,residue)
       
       !======================================================
       ! Subrutina para generar la entrada del fichero rtp con
       ! la topología de un residuo. Usa un resido (solo uno) 
       ! no el sistema
       !=======================================================

        implicit none

        integer,intent(in)::unt
        type(str_resmol),intent(in)::residue
        
        !local
        integer :: nrexcl, &
                   res_number=1 !to be generalized
        real :: total_charge
        integer,dimension(3) :: date
        !Counter
        integer::i, j

        character(len=50) :: dummy_char

        !Header
!         call idate(date)
        write(dummy_char,'(i2,a1,i2,a1,i4)') date(1),"/",date(2),"/",date(3)
        write(unt,'(A,/,A)') '; Topology (rtp style) generated by "build_top"', &
                             '; on '//adjustl(dummy_char)

        !1. Section: [ moleculetype ]
        write(unt,'(A)') "[ "//trim(adjustl(residue%name))//" ]"

        !2. Section: [ atoms ]
        write(unt,'(A)') ' [ atoms ]'
        write(unt,'(A)')   ';   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB'
        total_charge=0.
        do i=1,residue%natoms
            total_charge=total_charge+residue%atom(i)%q
            write(unt,400) adjustl(residue%atom(i)%name)  ,   &
                           residue%atom(i)%fftype,   &
                           residue%atom(i)%q     ,   &
                           i                           !this is CHARMM specific (no charge groups)
              
        enddo                                           
                                                       
        !3. Section: [ bonds ]                          
        write(unt,'(A)') ' [ bonds ]'                
        write(unt,'(A)')   ';  ai    aj funct            c0            c1            c2            c3'
        do i=1,residue%natoms-1                        
            do j=1,residue%atom(i)%nbonds              
                if ( residue%atom(i)%connect(j) > i )  &
                  write(unt,401) residue%atom(i)%name, &
                                 residue%atom(residue%atom(i)%connect(j))%name
            enddo
        enddo

        !3. Section: [ pairs ]
        write(unt,'(A)') ' [ pairs ]'
        write(unt,'(A)')   ';  ai    aj funct            c0            c1            c2            c3'
        do i=1,residue%geom%npairs
            write(unt,401) residue%atom(residue%geom%pair(i,1))%name, &
                           residue%atom(residue%geom%pair(i,2))%name
        enddo

        !4. Section: [ angles ]
        write(unt,'(A)') ' [ angles ]'
        write(unt,'(A)')   ';  ai    aj    ak funct            c0            c1            c2            c3'
        do i=1,residue%geom%nangles
!                   if (adjustl(residue%atom(i)%fftype) == "CA" .and. &
!                       adjustl(residue%atom(residue%atom(i)%angle(j,1))%fftype) == "CA".and. & 
!                       adjustl(residue%atom(residue%atom(i)%angle(j,2))%fftype) == "CA") &
                write(unt,402) residue%atom(residue%geom%angle(i,1))%name, & !ai
                               residue%atom(residue%geom%angle(i,2))%name, & !aj
                               residue%atom(residue%geom%angle(i,3))%name    !ak
        enddo

        !5. Section: [ dihedrals ]
        write(unt,'(A)') ' [ dihedrals ]'
        write(unt,'(A)')   ';  ai    aj    ak    al funct            c0            c1            c2            c3'
        do i=1,residue%geom%ndihed
!                   if (adjustl(residue%atom(i)%fftype) == "CA" .and. &
!                       adjustl(residue%atom(residue%atom(i)%dihed(j,1))%fftype) == "CA".and. & 
!                       adjustl(residue%atom(residue%atom(i)%dihed(j,2))%fftype) == "CA".and. &
!                       adjustl(residue%atom(residue%atom(i)%dihed(j,3))%fftype) == "CA" ) &
                !Separate rotated bonds
                if ( residue%geom%dihed(i,2) /= residue%geom%dihed(i-1,2) .or. &
                     residue%geom%dihed(i,3) /= residue%geom%dihed(i-1,3) ) then
                    write(unt,'(A8,2I3)') "; Bond: ", residue%geom%dihed(i,2),residue%geom%dihed(i,3)
                endif
                write(unt,403) residue%atom(residue%geom%dihed(i,1))%name, & !ai
                               residue%atom(residue%geom%dihed(i,2))%name, & !aj
                               residue%atom(residue%geom%dihed(i,3))%name, & !ak
                               residue%atom(residue%geom%dihed(i,4))%name    !al
        enddo

!         write(unt,'(/,A)') '; Planar/tetrahedric impropers not generated'



        return
    
    !gro file format
400 format(7x,a5,3x,a5,x,f8.3,x,i5)
401 format(7x,a5,3x,a5)
402 format(7x,a5,3x,a5,3x,a5)
403 format(7x,a5,3x,a5,3x,a5,3x,a5)
!201 format(3f10.5)

    end subroutine write_rtp


    subroutine read_top(unt,residue,molmap,nmol)

       !======================================================
       ! Reads topology (top) files with the itp files therein
       ! with the include-mechanism.
       ! 
       ! LIMITATIONS
       ! Only system containing molecules with one residue are 
       ! supported. Clarifications on the molecule/residue types
       ! is needed to improve
       !
       ! CHANGES
       ! Added GMXLIB to read standard itp if not found in current dir
       ! Bug: some itp (e.g. ions.itp) contain more residues that are
       !      not in the conf file. They must be selected by name!
       !=======================================================

        use line_preprocess

        integer,intent(in)::unt
        type(str_resmol),intent(inout),dimension(:)::residue
        character(len=*),dimension(:),intent(out) :: molmap
        integer,intent(out) :: nmol
        
        !Counter
        integer::i,j,k, imap, &
                 ii, jj, kk, ll, ift

        !OTher system features
        integer :: ires, ires_new, imol

        !Reading stuff
        character(len=260) :: line, section
        logical :: is_section
        integer :: ios
        character(len=200) :: itpfile, GMXLIB
        integer :: S_TOP = 30

        !Aux
        character(len=1) :: cnull
        integer :: inull
        character(len=50) :: dummy_char 
        real(8) :: a,b,c,d,e,f 

        !Get GMX lib variable
        call GET_ENVIRONMENT_VARIABLE("GMXDATA",GMXLIB)
        if ( len_trim(GMXLIB) == 0 ) then
            print'(/,A,/)', "WARNING: GMXDATA env. variable undefined"
        else
            GMXLIB=trim(adjustl(GMXLIB))//"/gromacs/top"
        endif

        ires=-1 !Label for residues (as read in top)
        i=0     !Residue counter (no es muy buen nombre, la verdad)
        k=0     !Atom counter
        imap=0
        nmol=0
        do 

            read(unt,'(A)',iostat=ios) line
            if (ios /= 0) exit

            call split_line(line,";",line,cnull)

            !If blank line or all comment, ignore the line
            line=adjustl(line)
            if ( len_trim(line) == 0 ) cycle

            !Is it a section, which one?
            is_section=.false.
            if (line(1:1) == "[") is_section=.true.

            !Include file support
            if (line(1:8) == "#include") then
                read(line(9:260),*) itpfile
!                 print'(/,A)', "Will include: "//trim(adjustl(itpfile))
                open(S_TOP,file=itpfile,iostat=ios,status="old")
                if (ios /= 0) then
                    print*, trim(adjustl(itpfile))//" not found in curent dir"
                    print*, "Searching in "//trim(adjustl(GMXLIB))
                    itpfile=trim(adjustl(GMXLIB))//"/"//trim(adjustl(itpfile))
                    open(S_TOP,file=itpfile,iostat=ios,status="old")
                endif
                if (ios /= 0) then
                    print*, "File "//trim(adjustl(itpfile))//" should be "//&
                            "included but it was not found. Check the output!"
                    cycle
                endif
                call read_top_recurs(S_TOP,residue,molmap,nmol,i)
                ires = -1
! print*, "READ", i, residue(i)%name, residue(i)%atom(1)%name, residue(i)%natoms
                close(S_TOP)
            else if (line(1:1) == "#") then
                !Other post-processing options are ignored
                cycle
            endif


            if (is_section) then
                call split_line(line   ,"[",cnull,section)
                call split_line(section,"]",section,cnull)
! print*, "Reading section: ", trim(adjustl(section))
            else if (adjustl(section) == "moleculetype") then
                cycle
!                 read(line,*) residue(i)%name, inull !nrexcl
            else if (adjustl(section) == "atoms") then
!This will not work here!! (only in the recursive version) -- TO BE SOLVED
!   !     1   CA_CR        1 BCR   C1        1   -0.100   12.011
                read(line,*) inull,cnull,ires_new
                if (ires /= ires_new) then
                    i=i+1
                    k=0
                endif
                k=k+1
                read(line,*) inull,                      & !serial
                             residue(i)%atom(k)%attype,  &
                             ires,                       & !residue nr
                             residue(i)%atom(k)%resname, &
                             residue(i)%atom(k)%name,    &
                             inull,                      & !chrgrp
                             residue(i)%atom(k)%q,       &
                             residue(i)%atom(k)%mass
                residue(i)%natoms = k
                residue(i)%name = residue(i)%atom(k)%resname
            else if (adjustl(section) == "atomtypes") then
                write(20,*) line
            else if (adjustl(section) == "bonds") then
                read(line,*) ii, jj, ift, a, b
                write(21,*) residue(i)%atom(ii)%attype//" "//residue(i)%atom(jj)%attype//" ",ift, a, b 
            else if (adjustl(section) == "angles") then
                read(line,*) ii, jj, kk, ift, a, b
                write(22,*) residue(i)%atom(ii)%attype//" "//residue(i)%atom(jj)%attype//" "//residue(i)%atom(kk)%attype//&
                            " ",ift, a, b 
            else if (adjustl(section) == "dihedrals") then
                read(line,*) ii, jj, kk, ll, ift!, a, b, c, d, e, f
                write(23,*) residue(i)%atom(ii)%attype//" "//residue(i)%atom(jj)%attype//" "//residue(i)%atom(kk)%attype//&
                            " "//residue(i)%atom(ll)%attype//"",ift!, a, b, c, d, e, f 
            else if (adjustl(section) == "molecules") then
                read(line,*) molmap(nmol+1), imol
                do j=nmol+2,nmol+imol
                    molmap(j) = molmap(j-1)
                enddo
                nmol = nmol + imol
            else
                cycle
            endif

        enddo
        close(unt)
                 
        return

        contains
        !Recursive read (to support include files)

        subroutine read_top_recurs(unt,residue,molmap,imol,irs)

            integer,intent(in)::unt
            type(str_resmol),intent(inout),dimension(:)::residue
            character(len=*),dimension(:),intent(inout) :: molmap
            integer,intent(inout) :: imol, irs
            
            !Counter
            integer::i, k, imap

            !OTher system features
            integer :: ires, ires_new

            !Reading stuff
            character(len=260) :: line, section
            logical :: is_section
            integer :: ios

            !Aux
            character(len=1) :: cnull
            integer :: inull
            character(len=50) :: dummy_char  

            ires=-1 !Label for residues (as read in top)
            k=0     !Atom counter
            nmol = imol
            i=irs
            do 

                read(unt,'(A)',iostat=ios) line
                if (ios /= 0) exit

                call split_line(line,";",line,cnull)

                !If blank line or all comment, ignore the line
                line=adjustl(line)
                if ( len_trim(line) == 0 ) cycle
                if (line(1:8) == "#include") then
!                     print*, "Nested #include files are not yet supported. Ignoring:"
!                     print*, trim(adjustl(line))
                    cycle
                elseif (line(1:1) == "#") then
                    !Other post-processing intructions are simply ignored 
                    cycle
                endif

                !Is it a section, which one?
                is_section=.false.
                if (line(1:1) == "[") is_section=.true.


                if (is_section) then
                    call split_line(line   ,"[",cnull,section)
                    call split_line(section,"]",section,cnull)
! print*, "Reading section(itp): ", trim(adjustl(section)), i
                else if (adjustl(section) == "moleculetype") then
                    i=i+1
                    k=0
                    cycle
                    !Was trying to read resid(i_not-updated). Now read resname from [atoms]
                    !read(line,*) residue(i)%name, inull !nrexcl
                else if (adjustl(section) == "atoms") then
!   !   !     1   CA_CR        1 BCR   C1        1   -0.100   12.011
                    k=k+1
! print*, trim(adjustl(line))
                    read(line,*) inull,                      & !serial
                                 residue(i)%atom(k)%attype,  &
                                 ires,                       & !residue nr
                                 residue(i)%atom(k)%resname, &
                                 residue(i)%atom(k)%name,    &
                                 inull,                      & !chrgrp
                                 residue(i)%atom(k)%q!,       &
!                                  residue(i)%atom(k)%mass   (TO SOLVE: TIP3P in CHARMM does not have mass item!!! Temporary solution is to disable this
                    residue(i)%natoms = k
                    residue(i)%name = residue(i)%atom(k)%resname
                    !Update residue (if more than one per molecule type) -- ONGOING (need checking)
                    read(line,*) inull,dummy_char,ires_new
                    if (ires /= ires_new) then
                        i=i+1
                        read(dummy_char,*,iostat=ios) residue(i)%name
                        if (ios /= 0) then
                            print*, "ERROR: couldn't get resname", dummy_char
                            stop
                        endif
                        k=0
                    endif
                else if (adjustl(section) == "molecules") then
                    read(line,*) molmap(nmol+1), imol
                    do j=nmol+2,nmol+imol
                        molmap(j) = molmap(j-1)
                    enddo
                    nmol = nmol + imol
                else
                    cycle
                endif

            enddo
            close(unt)

           !Updato outputs
           irs = i
           imol = nmol
                     
            return

        end subroutine read_top_recurs

    end subroutine read_top

    subroutine read_gro_hess(unt,N,Hess,error)

        !Description
        ! Read Hessian from ascii file generated from an .mtx file (with gmxdump)
        ! 
        ! Error codes: 
        !  1: not square matrix

        integer,intent(in)::unt
#ifdef DOUBLE
        real(8),dimension(:,:),intent(out)::Hess
#else
        real(4),dimension(:,:),intent(out)::Hess
#endif 
        integer,intent(out)::N
        integer,intent(out)::error
        
        !local
        character :: cnull
        !Counter
        integer::i, j
         
        read(unt,'(A)') cnull
        read(unt,*) N, i
        if (N /= i) then
            error = 1
            return
        endif

        !Read in triangular form
        j=1
        do i=1,N
            read(unt,*) Hess(i,1:N)
        enddo

        ! UNIT CONVERSION                                        ! GROMACS       --> Atomic Units  
        Hess(1:N,1:N)=Hess(1:N,1:N)/CALtoJ/HtoKCALM*BOHRtoNM**2  ! KJ/mol * nm-2 --> Hartree * bohr-2

        return

    end subroutine read_gro_hess

end module gro_manage
