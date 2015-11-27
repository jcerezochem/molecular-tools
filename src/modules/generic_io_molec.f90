module generic_io_molec

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to manage output files 
    !  from different QM codes. Relies in specific modules for
    !  each package.
    !    
    !==============================================================

    !Common declarations:
    !===================
    use structure_types
    use generic_io
    implicit none

    contains


    subroutine generic_strmol_reader(unt,filetype,molec,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic geometry reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! molec   (io)   str_resmol    Molecule
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        use molecular_structure

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        type(str_resmol),intent(inout)  :: molec
        integer,intent(out),optional    :: error_flag

        !Local
        integer   :: error_local
        character(len=200) :: msg

        call generic_structure_reader(unt,filetype,molec%natoms,      &
                                                   molec%atom(:)%x,   &
                                                   molec%atom(:)%y,   &
                                                   molec%atom(:)%z,   &
                                                   molec%atom(:)%mass,&
                                                   molec%atom(:)%name,&
                                      error_local)

        call atname2element(molec)

        ! Reader provides coordiantes in Angstrong
        molec%units = "Angs"

        ! Error handling
        if (error_local /= 0) then
            write(msg,'(A,I0)') "ERROR readig structure from file. Error code: ", error_local
            call alert_msg("fatal",msg)
        endif

        if (present(error_flag)) error_flag=error_local
        

        return

    end subroutine generic_strmol_reader


    subroutine generic_strmol_writer(unt,filetype,molec,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic geometry reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! molec   (io)   str_resmol    Molecule
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        use molecular_structure

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        type(str_resmol),intent(inout)  :: molec
        integer,intent(out),optional    :: error_flag

        !Local
        integer   :: error_local
        character(len=200) :: msg
        character(len=4) ::  current_units

        ! Save current_units and ensure Angs
        current_units=molec%units
        call set_geom_units(molec,"Angs")

        call generic_structure_writer(unt,filetype,molec%natoms,      &
                                                   molec%atom(:)%x,   &
                                                   molec%atom(:)%y,   &
                                                   molec%atom(:)%z,   &
                                                   molec%atom(:)%mass,&
                                                   molec%atom(:)%name,&
                                      error_local)

        ! Error handling
        if (error_local /= 0) then
            write(msg,'(A,I0)') "ERROR writting structure to file. Error code: ", error_local
            call alert_msg("fatal",msg)
        endif

        if (present(error_flag)) error_flag=error_local
        
        ! Reset input units before leaving
        call set_geom_units(molec,adjustl(current_units))

        return

    end subroutine generic_strmol_writer

end module generic_io_molec

