module xyz_manage

    contains

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  Basic I/O subroutines to read/write XYZ files.
    !  WITHOUT structure_types
    !
    ! Notes

    !-----------------------------------------------------


    ! ====================
    ! GENERAL INTERFACE WITHOUT structure_types
    !======================
    subroutine read_xyz_natoms(unt,Nat)

        !Read standard xyz files
        use constants

        integer,intent(in)::unt
        integer,intent(out)::Nat

        read(unt,*) Nat

        return

    end subroutine read_xyz_natoms

    subroutine read_xyz_geom(unt,Nat,AtName,X,Y,Z,title)

        !Read standard xyz files

        use constants

        integer,intent(in)  :: unt
        integer,intent(out) :: Nat
        character(len=*), dimension(:), intent(out) :: AtName
        real(kind=8), dimension(:), intent(out) :: X,Y,Z
        character(len=*),intent(out),optional :: title

        !local
        integer::i, AtNum
        character(len=1) :: test_name


        read(unt,*) Nat
        if (present(title)) then
            read(unt,'(A)') title
        else
            read(unt,'(A)') test_name
        endif

        do i=1,Nat   

            read(unt,*) AtName(i),    &
                        X(i),         &
                        Y(i),         &
                        Z(i)

        enddo

        !If atomic numbers are used instead of atom name, make the transformation
        test_name = adjustl(AtName(1))
        if (test_name == "1" .or.&
            test_name == "2" .or.&
            test_name == "3" .or.&
            test_name == "4" .or.&
            test_name == "5" .or.&
            test_name == "6" .or.&
            test_name == "7" .or.&
            test_name == "8" .or.&
            test_name == "9") then
            
            do i=1,Nat
                read(AtName(i),*) AtNum
                AtName(i) = atname_from_atnum(AtNum)
                !-1 are dummies
                if (AtNum == -1) AtName(i)="X"
            enddo

        endif

        return

    end subroutine read_xyz_geom



end module xyz_manage
