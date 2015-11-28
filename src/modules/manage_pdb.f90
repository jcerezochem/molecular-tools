module pdb_manage

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  Version not based on str_resmol types
    !==============================================================


    contains

    subroutine read_pdb_geom(unt,Nat,AtName,X,Y,Z, &
                             !optional
                             title,attype,resname,alter_loc,chain,ins_code,resseq)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from pdb. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vertor    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        !
        !Notes
        ! Natoms is guessed from the number of entries starting with 
        ! ATOM or HETATM
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out) :: X,Y,Z
        character(len=*),intent(out),optional :: title
        character(len=*),dimension(:),intent(out),optional :: attype,    &
                                                 resname,   &
                                                 alter_loc, &
                                                 chain,     &
                                                 ins_code
        integer,dimension(:),intent(out),optional :: resseq
                                                
        !local
        integer::i, natoms, ii, ios
        character(len=260) :: line
        character :: dummy_char
        !pdb attributes 
        character(len=99):: title_local
        character(len=6),dimension(1:size(X)) :: attype_local,    &
                            resname_local
        character(len=1),dimension(1:size(X)) :: alter_loc_local, &
                                                 chain_local,     &
                                                 ins_code_local
        integer,dimension(1:size(X))          :: resseq_local
        logical :: have_read_structure=.false.

        ! Compared with g96, we don't know the end of the section. To animations we
        ! would need to implement TER or MODEL entries. For the moment, we stop
        ! reading when the second TITLE section is found (not that good, though, we
        ! don't know how generic PDB files are writting. AFAIK The standard do not 
        ! prevent from reusing TITLE section 

        read(unt,'(A)',iostat=ios) line
        i=0
        do while (ios==0)        

            if (line(1:6) == "ATOM  " .or. line(1:6) == "HETATM") then
                i=i+1
                read(line,100) attype_local(i),    &
                               !Serial
                               AtName(i),          &
                               alter_loc_local(i), &
                               resname_local(i),   &
                               chain_local(i),     &
                               resseq_local(i),    &
                               ins_code_local(i),  &
                               X(i),         &
                               Y(i),         &
                               Z(i)
                have_read_structure=.true.
            else if (line(1:6) == "TITLE ") then
                if (have_read_structure) exit
                read(line,101) title_local
            endif

            read(unt,'(A)',iostat=ios) line

        enddo
        Nat=i

        ! Take requested optional data
        if (present(title))     title = title_local
        if (present(attype))    attype = attype_local
        if (present(resname))   resname = resname_local
        if (present(alter_loc)) alter_loc = alter_loc_local
        if (present(chain))     chain = chain_local
        if (present(ins_code))  ins_code = ins_code_local
        if (present(resseq))    resseq = resseq_local

        return

    !Formatos
100 format(A6,5X,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2) !Serial not read
101 format(6X,A)

    end subroutine read_pdb_geom


    subroutine write_pdb_geom(unt,Nat,AtName,X,Y,Z, &
                             !optional
                             title,attype,resname,alter_loc,chain,ins_code,resseq)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Get geometry and atom names from pdb. The number of atoms
        ! is also taken
        !
        !Arguments
        ! unt     (inp) int /scalar    unit for the file 
        ! Nat     (out) int /scalar    Number of atoms
        ! AtName  (out) char/vertor    Atom names
        ! X,Y,Z   (out) real/vectors   Coordinate vectors (ANGSTRONG)
        !
        !Notes
        ! Natoms is guessed from the number of entries starting with 
        ! ATOM or HETATM
        !==============================================================

        integer,intent(in)  :: unt
        integer,intent(in)  :: Nat
        character(len=*), dimension(:), intent(in) :: AtName
        real(kind=8),dimension(:), intent(in) :: X,Y,Z
        character(len=*),intent(in),optional :: title
        character(len=*),dimension(:),intent(in),optional :: attype,    &
                                                 resname,   &
                                                 alter_loc, &
                                                 chain,     &
                                                 ins_code
        integer,dimension(:),intent(in),optional :: resseq
                                                
        !local
        integer::i, natoms, ii, ios
        character(len=260) :: line
        character :: dummy_char
        !pdb attributes 
        character(len=99):: title_local
        character(len=6),dimension(1:Nat) :: attype_local,    &
                                             resname_local
        character(len=1),dimension(1:Nat) :: alter_loc_local, &
                                             chain_local,     &
                                             ins_code_local
        integer,dimension(1:Nat)          :: resseq_local

        ! Take defaults for non-defined optional args
        title_local = "Generated with write_pdb_geom"
        attype_local(1:Nat) = "ATOM"
        resname_local(1:Nat) = "MOL"
        alter_loc_local(1:Nat) = ""
        chain_local(1:Nat) = ""
        ins_code_local(1:Nat) = ""
        resseq_local(1:Nat) = 1
        ! Update with defined ones
        if (present(title))     title_local            = title
        if (present(attype))    attype_local(1:Nat)    = attype(1:Nat)
        if (present(resname))   resname_local(1:Nat)   = resname(1:Nat)
        if (present(alter_loc)) alter_loc_local(1:Nat) = alter_loc(1:Nat)
        if (present(chain))     chain_local(1:Nat)     = chain(1:Nat)
        if (present(ins_code))  ins_code_local(1:Nat)  = ins_code(1:Nat)
        if (present(resseq))    resseq_local(1:Nat)    = resseq(1:Nat)

        write(unt,'(A)') "TITLE "//trim(adjustl(title_local))
        i=0
        do i=1,Nat
                write(unt,200) attype_local(i),    &
                               i,                  &
                               AtName(i),          &
                               alter_loc_local(i), &
                               resname_local(i),   &
                               chain_local(i),     &
                               resseq_local(i),    &
                               ins_code_local(i),  &
                               X(i),               &
                               Y(i),               &
                               Z(i)
        enddo

        return

    !Formatos
200 format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,10X,A2,A2) !oficial PDB
201 format(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,10X,A2,A2) !modified to show element

    end subroutine write_pdb_geom

end module pdb_manage
