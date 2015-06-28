module chem_io

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS
    !==============================================================
    ! Description
    !  Subroutine to interface all read/write structure subroutines
    !
    !-----------------------------------------------------

    contains

    subroutine generic_strfile_read(unt,filetype,molec)

        integer, intent(in) :: unt
        character(len=*),intent(inout) :: filetype
        type(str_resmol),intent(inout) :: molec

        !local
        type(str_molprops) :: props

        select case (adjustl(filetype))
            case("gro")
             call read_gro(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("g96")
             call read_g96(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("pdb")
             call read_pdb_new(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("log")
             call parse_summary(I_INP,molec,props,"struct_only")
             call atname2element(molec)
             call assign_masses(molec)
            case("fchk")
             call read_fchk_geom(I_INP,molec)
             call atname2element(molec)
!              call assign_masses(molec) !read_fchk_geom includes the fchk masses
            case("UnSym")
             call read_molcas_geom(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case("psi4")
             call read_psi_geom(I_INP,molec)
             call atname2element(molec)
             call assign_masses(molec)
            case default
             call alert_msg("fatal","File type not supported: "//filetype)
        end select

        return

    end subroutine generic_strfile_read

end module chem_io
