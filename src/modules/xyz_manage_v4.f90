module xyz_manage

    contains

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  Basic I/O subroutines to read/write XYZ files.
    !    subroutine read_xyz(unt,system)
    !    subroutine write_xyz(unt,system)
    !
    ! Dependences
    !  Uses "structure_types" module
    !
    ! Notes
    !  Since variables molec, residue... are allocatable, they should
    !  be passed to the sr even if they are output to get the allocation.
    !
    ! History
    !  2013/03/08  Born (el adjetivo, no la persona)
    !  2013/06/14  Added gaussian com writer
    !-----------------------------------------------------


    subroutine read_xyz(unt,system)

        !Read standard xyz files

        use structure_types

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, natoms
        character(len=1) :: test_name

        character(len=5),dimension(103) :: atom_names_from_atnum
 
        !This should be elsewhere (constants_mod?)
        data atom_names_from_atnum(1:103) &
         /'H' ,                                                                                'He',&
          'Li','Be',                                                  'B' ,'C' ,'N' ,'O' ,'F' ,'Ne',&
          'Na','Mg',                                                  'Al','Si','P' ,'S' ,'Cl','Ar',&
          'K' ,'Ca','Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',&
          'Rb','Sr','Y' ,'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I' ,'Xe',&
          'Cs','Ba','La',& !Lantanides:  
!                  ---------------------------------------------------
                    'Ce','Pr','Nd','Pm','Sm','Eu','Gd',&
                    'Tb','Dy','Ho','Er','Tm','Yb','Lu',&
!                  ---------------------------------------------------
                         'Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',&
         'Fr','Ra','Ac',& !Actinides:
!                  ---------------------------------------------------
                   'Th','Pa','U' ,'Np','Pu','Am','Cm',&
                   'Bk','Cf','Es','Fm','Md','No','Lr'&
!                  ---------------------------------------------------
         /


        read(unt,*) natoms
        system%natoms=natoms
        read(unt,*) system%title

        do i=1,natoms   

            read(unt,*) system%atom(i)%name,      &
                        system%atom(i)%x,         &
                        system%atom(i)%y,         &
                        system%atom(i)%z

        enddo

        !If atomic numbers are used instead of atom name, make the transformation
        test_name = adjustl(system%atom(1)%name)
        if (test_name == "1" .or.&
            test_name == "2" .or.&
            test_name == "3" .or.&
            test_name == "4" .or.&
            test_name == "5" .or.&
            test_name == "6" .or.&
            test_name == "7" .or.&
            test_name == "8" .or.&
            test_name == "9") then
            
            do i=1,natoms
                read(system%atom(i)%name,*) system%atom(i)%AtNum
                system%atom(i)%name = atom_names_from_atnum(system%atom(i)%AtNum)
                !-1 are dummies
                if (system%atom(i)%AtNum == -1) system%atom(i)%name="X"
            enddo

        endif

        return

    end subroutine read_xyz


    subroutine write_xyz(unt,system)

        !Write standard xyz files

        use structure_types

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i


        write(unt,*) system%natoms
        if (len_trim(system%title) == 0) then
            write(unt,*) "File generated with write_xyz subroutine"
        else
            write(unt,*) system%title
        endif

        do i=1,system%natoms   

            write(unt,*) system%atom(i)%name,      &
                         system%atom(i)%x,         &
                         system%atom(i)%y,         &
                         system%atom(i)%z

        enddo

        return

    end subroutine write_xyz

    subroutine write_gcom(unt,system)

        !Write gaussian com file (in cartesian coord)

        use structure_types
        use line_preprocess

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, natoms
        character(len=100) :: chkname
        character(len=1) :: null

        write(unt,'(A)') "--link1--"
        write(unt,'(A)') "%mem=2GB"
        !Supposing that file name is passed through title
        call split_line(system%title,".",chkname,null)
        write(unt,'(A)') "%chk="//trim(adjustl(chkname))//".chk"
        write(unt,'(A)') "%nproc=8"
        !Indicate a SP with "standard" model chem
        if (adjustl(system%job%type) == "XX") then
            !Sensible default
            write(unt,'(A)') "#p B3LYP/6-31G(d)"
        else 
            write(unt,'(A)') "#p "//trim(adjustl(system%job%type))//" "&
                           //trim(adjustl(system%job%method))//"/"&
                           //trim(adjustl(system%job%basis))
        endif
        write(unt,'(A)') ""
        write(unt,'(A)') trim(adjustl(system%job%title))
        write(unt,'(A)') ""
        !Charge and multiplicity (TODO: read from system attributes)
        write(unt,'(A)') "0 1"
        !Geomertry
        do i=1,system%natoms
            write(unt,'(A10,X,3(F15.6,X))') system%atom(i)%name, system%atom(i)%x,system%atom(i)%y,system%atom(i)%z
        enddo
        write(unt,'(A)') ""

        return

    end subroutine write_gcom


end module xyz_manage
