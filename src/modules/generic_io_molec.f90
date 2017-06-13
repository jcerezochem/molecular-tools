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
    ! The following need to be Supported by generic_io module
    use gro_manage
    use xyz_manage_molec
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

        ! Readers still unsupported by generic_io module
        if (adjustl(filetype)=="xxx") then
            ! all parser are now in generic_structure_reader
            ! but keep this if structure for the moment

        ! any other format is assumed to be supported
        else 
            call generic_structure_reader(unt,filetype,molec%natoms,      &
                                                       molec%atom(:)%x,   &
                                                       molec%atom(:)%y,   &
                                                       molec%atom(:)%z,   &
                                                       molec%atom(:)%mass,&
                                                       molec%atom(:)%name,&
                                                      ! Additional info (not available for all filetypes)
                                                      molec%atom(:)%resname,&
                                          error_flag=error_local)
        endif

        call atname2element(molec)
        call assign_atnum_molec(molec)
        ! Assign masses for formats where the name can be different from the element (so the assigment
        ! is not possible at generic_structure_reader)
        if (adjustl(filetype) == 'gro' .or. adjustl(filetype) == 'pdb' .or. adjustl(filetype) == 'g96') then
            call assign_masses_molec(molec)
        endif

        ! Reader provides coordiantes in Angstrong
        molec%units = "Angs"

        ! Error handling
        ! Catch it if error_flag is present. Otherwise stop if error/=0
        if (present(error_flag)) then
            error_flag=error_local
        else if (error_local /= 0) then
            write(msg,'(A,I0)') "reading structure from file. Error code: ", error_local
            call alert_msg("fatal",msg)
        endif
        

        return

    end subroutine generic_strmol_reader


    subroutine generic_strmol_writer(unt,filetype,molec,error_flag,title)

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
        character(len=*),intent(in),optional :: title

        !Local
        integer   :: error_local
        character(len=200) :: msg, title_local

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with Angstrong
        current_units=molec%units
        call set_geom_units(molec,"Angs")
        !----------------------------------------------

        ! Initialize error local (old readers don't use it)
        error_local=0

        if (present(title)) then
            title_local = adjustl(title)
        else
            title_local = "Structure written with generic_io_molec"
        endif

        ! Readers still unsupported by generic_io module (more that in the case of readers)
        if (adjustl(filetype)=="gro") then
            molec%title=title_local
            call write_gro(unt,molec)
        elseif (adjustl(filetype)=="g96") then
            molec%title=title_local
            call write_g96(unt,molec)
        elseif (adjustl(filetype)=="xyz") then
            molec%title=title_local
            call write_xyz(unt,molec)

        ! any other format is assumed to be supported
        else 
            call generic_structure_writer(unt,filetype,molec%natoms,      &
                                                       molec%atom(:)%x,   &
                                                       molec%atom(:)%y,   &
                                                       molec%atom(:)%z,   &
                                                       molec%atom(:)%mass,&
                                                       molec%atom(:)%name,&
                                          error_flag=error_local,         &
                                          title=title_local)
        endif

        ! Error handling
        ! Catch it if error_flag is present. Otherwise stop if error/=0
        if (present(error_flag)) then
            error_flag=error_local
        else if (error_local /= 0) then
            write(msg,'(A,I0)') "writting structure to file. Error code: ", error_local
            call alert_msg("fatal",msg)
        endif
        
        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! Revert original units
        call set_geom_units(molec,adjustl(current_units))
        !----------------------------------------------

        return

    end subroutine generic_strmol_writer

end module generic_io_molec

