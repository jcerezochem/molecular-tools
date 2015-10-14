module molden_manage

    use constants

    contains

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  Basic I/O subroutines to read/write XYZ files.
    !    subroutine read_molden(unt,system)
    !    subroutine write_molden(unt,system)
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


    subroutine read_molden(unt,system)

        use structure_types
        use alerts

       !Reads MOLDEN structure file from [ATOMS] (input/output in ANGS)

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, natoms, ii, ios
        character(len=260) :: line
        character(len=10) :: units
        character :: cnull
        !The use of get_structure=.false. fails probably due to
        ! memory issues (change without being directly accesed)
        logical :: not_get_structure=.true.

        !The file is organized in sections. Each begining with a given name
        ! This SR scans all the file till finding [ATOMS]. Some files  
        ! have [N_ATOMS], but others doesn't. So we set natoms from the 
        ! number of lines read in [ATOMS] section
        system%natoms = 0
        do

            read(unt,'(A)',iostat=ios) line
            if (ios /= 0) exit

            if (index(line,'[N_ATOMS]') /= 0) then
                read(unt,*) system%natoms
            elseif ( (index(line,'[ATOMS]') /= 0) .or. &
                     (index(line,'[Atoms]') /= 0) ) then
                not_get_structure=.false.
                read(line,*) cnull, units
                i=0
                do 
                     read(unt,'(A)') line
                     if ( (adjustl(line) == "END") .or. &
                          (index(line,'[') /= 0  ) ) exit
                     i=i+1
                     read(line,*) system%atom(i)%name,   &
                                  ii,                    & !serial
                                  system%atom(i)%AtNum,  &
                                  system%atom(i)%x,      &
                                  system%atom(i)%y,      &
                                  system%atom(i)%z
                     !Get elements from AtNum
                     system%atom(i)%element = atname_from_atnum(system%atom(i)%AtNum)
                     
                enddo
                if (i /= system%natoms) call alert_msg("note","[N_ATOMS] not found of conflictive in molden file")
                system%natoms = i
                !Change units if needed
                if (adjustl(units) == "(AU)") then
                    system%atom(1:i)%x = system%atom(1:i)%x*BOHRtoANGS
                    system%atom(1:i)%y = system%atom(1:i)%y*BOHRtoANGS
                    system%atom(1:i)%z = system%atom(1:i)%z*BOHRtoANGS
                endif
      
                !We are done with the file
                exit

            endif

        enddo

        rewind(unt)

        return 

    end subroutine read_molden


end module molden_manage
