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
        read(unt,'(A)') system%title

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
            write(unt,'(A)') system%title
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


    subroutine write_fcc(unt,system)

        !Writes the structure related part of a fcc input (state)
        !
        !

        use structure_types
        use line_preprocess

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, natoms, j
        character(len=180) :: line
        character(len=1) :: null
        character(len=9) :: mark
        real(8) :: x,y,z

        !
        do i=1,system%natoms   

            write(unt,*) system%atom(i)%x
            write(unt,*) system%atom(i)%y
            write(unt,*) system%atom(i)%z

        enddo

        return

    end subroutine write_fcc


    subroutine read_gcom_multi2fccinp(unt,system)

        !Read gaussian com file (in cartesian coord, not Zmat..)
        !CAUTION: this need exahustive testing to cope with different 
        ! input styles. 
        !
        !Notes
        ! Also reads nested (--link1--) files. The first one is taken?
        !

        use structure_types
        use line_preprocess

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, natoms, j
        character(len=180) :: line, output
        character(len=1) :: null
        character(len=9) :: mark
        real(8) :: x,y,z

        mark="--link1--"
        i=0
        do while (mark=="--link1--")
            i=i+1
            write(output,*) i
            output="fcc_"//trim(adjustl(output))//".dat"
            open(20,file=output)
            read(unt,'(A)') line
            do while (index(line,'#')==0)        
                read(unt,'(A)') line
            enddo
            read(unt,'(A)') line
            !Read till title section
            do while (trim(line)=="")        
                read(unt,'(A)') line
            enddo
            !Read till title section
            read(unt,'(A)') line
            do while (trim(line)=="")        
                read(unt,'(A)') line
            enddo

            !From now, read the strucuture
            do        
                read(unt,'(A)') line
                if (trim(line) == "") exit
                read(line,*) null, x, y, z
                write(20,*) x
                write(20,*) y
                write(20,*) z
            enddo
            close(20)
            !Look for another mark
            read(unt,*) mark
            print*, mark
        enddo

        return

    end subroutine read_gcom_multi2fccinp

    subroutine read_gcom_multi2xyz(unt,system)

        !Read gaussian com file (in cartesian coord, not Zmat..)
        !CAUTION: this need exahustive testing to cope with different 
        ! input styles. 
        !
        !Notes
        ! Also reads nested (--link1--) files. The first one is taken?
        !

        use structure_types
        use line_preprocess

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, natoms, j
        character(len=180) :: line, output, title
        character(len=1) :: null
        character(len=9) :: mark, atname
        real(8) :: x,y,z

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

        mark="--link1--"
        i=0
        do while (mark=="--link1--")
            i=i+1
            write(output,*) i
            output="fcc_"//trim(adjustl(output))//".xyz"
            open(20,file=output)
            read(unt,'(A)') line
            do while (index(line,'#')==0)        
                read(unt,'(A)') line
            enddo
            read(unt,'(A)') line
            !Read till title section
            do while (trim(line)=="")        
                read(unt,'(A)') line
            enddo
            title=line
            !Read till title section
            read(unt,'(A)') line
            do while (trim(line)=="")        
                read(unt,'(A)') line
            enddo

            !From now, read the strucuture
            write(20,*) 21
            write(20,*) trim(adjustl(title))
            do        
                read(unt,'(A)') line
                if (trim(line) == "") exit
                read(line,*) j, x, y, z
                atname =  atom_names_from_atnum(j)
                write(20,*) trim(adjustl(atname)), x, y, z
            enddo
            close(20)
            !Look for another mark
            read(unt,*) mark
        enddo

        return

    end subroutine read_gcom_multi2xyz

end module xyz_manage
