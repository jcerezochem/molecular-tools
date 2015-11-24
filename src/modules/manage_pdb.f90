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
        character(len=*),intent(out),optional :: title,     &
                                                 attype,    &
                                                 resname,   &
                                                 alter_loc, &
                                                 chain,     &
                                                 ins_code
        integer,intent(out),optional :: resseq
                                                
        !local
        integer::i, natoms, ii, ios
        character(len=260) :: line
        character :: dummy_char
        !pdb attributes 
        character(len=99):: title_local
        character(len=6) :: attype_local,    &
                            resname_local
        character(len=1) :: alter_loc_local, &
                            chain_local,     &
                            ins_code_local
        integer          :: resseq_local

        read(unt,'(A)',iostat=ios) line
        i=0
        do while (ios==0)        

            if (line(1:6) == "ATOM  " .or. line(1:6) == "HETATM") then
                i=i+1
                read(line,100) attype_local,    &
                               !Serial
                               AtName,          &
                               alter_loc_local, &
                               resname_local,   &
                               chain_local,     &
                               resseq_local,    &
                               ins_code_local,  &
                               X,         &
                               Y,         &
                               Z
            else if (line(1:6) == "TITLE ") then
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

end module pdb_manage
