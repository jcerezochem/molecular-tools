module pdb_manage_molec

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information 
    !  from molpro pdb files:
    !    subroutine read_pdb(unt,natoms,system) -- deprecated, needs strip pdb
    !    subroutine read_pdb_new(unt,system)
    !    subroutine write_pdb(unt,system)
    ! Note: all subroutines rewind the file after using it
    ! Notes
    !  Since variables molec, residue... are allocatable, they should
    !  be passed to the sr even if they are output to get the allocation.
    !==============================================================


    contains

    subroutine read_pdb(unt,natoms,system)

        use structure_types

        integer,intent(in)::unt,natoms
        type(str_resmol),intent(inout)::system

        !local
        integer::i

        system%natoms=natoms

        !Rewind the file in case it was not in the gegining
        rewind(unt)
        do i=1,natoms        

            read(unt,100) system%atom(i)%attype,    &
                          !Serial
                          system%atom(i)%name,      &
                          system%atom(i)%alter_loc, &
                          system%atom(i)%resname,   &
                          system%atom(i)%chain,     &
                          system%atom(i)%resseq,    &
                          system%atom(i)%ins_code,  &
                          system%atom(i)%x,         &
                          system%atom(i)%y,         &
                          system%atom(i)%z

        enddo

        return

    !Formatos
100 format(A6,5X,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2) !Serial not read

    end subroutine read_pdb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine read_pdb_new(unt,system)

        use structure_types

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, ios=0
        character(len=100) :: line


!         system%natoms=natoms

        !Rewind the file in case it was not at the begining
        rewind(unt)

        read(unt,'(A)',iostat=ios) line
        i=0
        do while (ios==0)        

            if (line(1:6) == "ATOM  " .or. line(1:6) == "HETATM") then
                i=i+1
                read(line,100) system%atom(i)%attype,    &
                               !Serial
                               system%atom(i)%name,      &
                               system%atom(i)%alter_loc, &
                               system%atom(i)%resname,   &
                               system%atom(i)%chain,     &
                               system%atom(i)%resseq,    &
                               system%atom(i)%ins_code,  &
                               system%atom(i)%x,         &
                               system%atom(i)%y,         &
                               system%atom(i)%z
            else if (line(1:6) == "TITLE ") then
                read(line,101) system%title
            endif

            read(unt,'(A)',iostat=ios) line

        enddo
        system%natoms=i

        rewind(unt)

        return

    !Formatos
100 format(A6,5X,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2) !Serial not read
101 format(6X,A)

    end subroutine read_pdb_new

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine write_pdb(unt,system)

        use structure_types

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system !caution now is inout

        !local
        integer::i


        if (system%natoms > MAX_ATOMS) then
            print'(/,X,A26)', 'In "write_pdb" subroutine:'
            print'(X,A34,X,I10,X,A5)', 'ERROR: I cannot create a file with', system%natoms, 'atoms'
            stop 1
        endif

        write(unt,'(A)') 'TITLE  File created with pdb_manage module'

        do i=1,system%natoms        

!             !Some system info might me missing
!             if (system%atom(i)%attype /= "ATOM  " .and. &
!                 system%atom(i)%attype /= "HETATM") system%atom(i)%attype="ATOM  "
! !             if (len(system%atom(i)%alter_loc) == 0) &
!             system%atom(i)%alter_loc = "X"
!             system%atom(i)%chain = "X"
        
            write(unt,201) system%atom(i)%attype,    &
                           i,                        & !This is the serial number
                           adjustl(system%atom(i)%name),&
                           system%atom(i)%alter_loc, &
                           system%atom(i)%resname,   &
                           system%atom(i)%chain,     &
                           system%atom(i)%resseq,    &
                           system%atom(i)%ins_code,  &
                           system%atom(i)%x,         &
                           system%atom(i)%y,         &
                           system%atom(i)%z,         &
                           system%atom(i)%element      !Include the element in the B-factor column

        enddo

        return
    
    !Formatos
200 format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2) !oficial PDB
201 format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,10X,A2,A2) !modified to show element

    end subroutine write_pdb

    subroutine write_pdb_connect(unt,system)

        use structure_types

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system !caution now is inout

        !local
        integer::i,nbonds


        if (system%natoms > MAX_ATOMS) then
            print'(/,X,A26)', 'In "write_pdb" subroutine:'
            print'(X,A34,X,I10,X,A5)', 'ERROR: I cannot create a file with', system%natoms, 'atoms'
            stop 1
        endif

        write(unt,'(A)') 'TITLE  File created with pdb_manage module'

        do i=1,system%natoms        

!             !Some system info might me missing
!             if (system%atom(i)%attype /= "ATOM  " .and. &
!                 system%atom(i)%attype /= "HETATM") system%atom(i)%attype="ATOM  "
! !             if (len(system%atom(i)%alter_loc) == 0) &
!             system%atom(i)%alter_loc = "X"
!             system%atom(i)%chain = "X"
        
            write(unt,300) system%atom(i)%attype,    &
                           i,                        & !This is the serial number
                           system%atom(i)%name,      &
                           system%atom(i)%alter_loc, &
                           system%atom(i)%resname,   &
                           system%atom(i)%chain,     &
                           system%atom(i)%resseq,    &
                           system%atom(i)%ins_code,  &
                           system%atom(i)%x,         &
                           system%atom(i)%y,         &
                           system%atom(i)%z

        enddo

        do i=1,system%natoms 
            nbonds=system%atom(i)%nbonds
            write(unt,301) "CONECT", i,system%atom(i)%connect(1:nbonds)
!             do j=1,system%atom(i)%nbonds
!                 k=system%atom(i)%connect(j)
!                 if (k<i) cycle
!                 write(unt,301) "CONECT", i,k
!             enddo
        enddo

        return
    
    !Formatos
300 format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2)
301 format(A6,10I5)

    end subroutine write_pdb_connect



!     subroutine pdb_connect(unt,system)
!     ! 1.1. Coordinates (pdb file)
!     open(I_PDB,file=pdbfile,status='old')
!     open(S_SPDB,status='scratch')
!     ! First, strip the pdb file so that we can use the module
!     i=0
!     do
!         read(10,'(A)',IOSTAT=IOstatus) line
!         ! End of file
!         if ( IOstatus < 0 ) exit
! 
!         read(line(1:6),'(A)') label
! 
!         !Get the stripped line
!         if ( label == "ATOM  " .or. label == "HETATM" ) then
!             i=i+1
!             write(S_SPDB,'(A)') line
! 
!        elseif ( label == "CONECT" ) then
!             call parse_line(line,narg,arg)
!             !The first arg is the CONECT label, then the id followed by the connects
!             read(arg(2),*) id
!             do j=3,narg 
!                 read(arg(j),*) connect(id,j-2)
!             enddo
!             nbonds(id) = narg-2
! 
!        endif
!     enddo
!     close(I_PDB)
!     molec%natoms = i
! 
! 
!     end subroutine pdb_connect 


end module pdb_manage_molec
