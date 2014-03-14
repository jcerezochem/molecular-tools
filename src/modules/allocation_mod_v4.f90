module allocation

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to allocate system and residue
    !  arrays:
    !    sys_alloc(molecule,context)
    !    sys_dealloc(molecule,context)
    !    res_alloc(molecule,context)
    !    res_dealloc(molecule,context)
    !==============================================================

    !Common declarations:
    !===================
    use structure_types
    use alerts
    implicit none

    contains

    subroutine mol_alloc(molecule,context)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! This subroutine is called before using molecular (system) arrays to  
        ! allocate them. Calls specify the context where they are called in 
        ! order to control which arrays are allocated (not all are needed) 
        !Arguments
          ! molec (str_system,inout): molecule 
          ! context (character(len=10))
        !==============================================================

        type(str_resmol),intent(inout) :: molecule
        character(len=*),intent(in) :: context
        !Local
        integer :: i

        select case (adjustl(context))
            case ("atoms") 
                allocate(molecule%atom(1:MAX_ATOMS))

            case ("conections")
                if (.not.allocated(molecule%atom)) allocate(molecule%atom(1:MAX_ATOMS))
                do i=1,MAX_ATOMS
                    allocate(molecule%atom(i)%connect(1:MAX_CONNEXIONS))
                    allocate(molecule%atom(i)%env(1:MAX_CONNEXIONS+1))
                enddo

!             THIS IS A RESIDUE ATTRIBUTE!
!             case ("bond_param") 
!                 allocate(residue%geom%pair(1:MAX_DERIVED_CNX,2))
!                 allocate(residue%geom%angle(1:MAX_DERIVED_CNX,3))
!                 allocate(residue%geom%dihed(1:MAX_DERIVED_CNX,4))

            case default
                call alert_msg("fatal","Wrong use of mol_alloc subroutine")
        end select

        return

    end subroutine mol_alloc


    subroutine prop_alloc(props,context)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! This subroutine is called before using molecular (system) arrays to  
        ! allocate them. Calls specify the context where they are called in 
        ! order to control which arrays are allocated (not all are needed) 
        !Arguments
          ! molec (str_system,inout): molecule 
          ! context (character(len=10))
        !==============================================================

        type(str_molprops),intent(inout) :: props
        character(len=*),intent(in) :: context
        !Local
        integer :: i

        select case (adjustl(context))

            case ("td_properties") 
                allocate(props%ddip(1:3*MAX_ATOMS,1:3))
                allocate(props%trans_ddip(1:3*MAX_ATOMS,1:3))

            case ("frequencies") 
                allocate(props%freq(1:1000))
                allocate(props%L(1:1000,1:1000))
                allocate(props%H(1:100000))

            case ("bond_param") 
                allocate(molecule%geom%bond(1:MAX_DERIVED_CNX,2))
                allocate(molecule%geom%pair(1:MAX_DERIVED_CNX,2))
                allocate(molecule%geom%angle(1:MAX_DERIVED_CNX,3))
                allocate(molecule%geom%dihed(1:MAX_DERIVED_CNX,4))

            case default
                call alert_msg("fatal","Wrong use of prop_alloc subroutine")
        end select

        return

    end subroutine prop_alloc


    subroutine mol_dealloc(molecule,context)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! This subroutine is called after using molecular (system) arrays to  
        ! deallocate them. Calls specify the context where they are called in 
        ! order to control which arrays are allocated (if not all).
        !Arguments
          ! molec (str_system,inout): molecule 
          ! context (character(len=10))
        !==============================================================

        use structure_types
!        use gaussian_manage
        use alerts

        implicit none

        type(str_resmol),intent(inout) :: molecule
        character(len=*),intent(in) :: context
        !Local
        integer :: i

        select case (adjustl(context))
            case ("atoms") 
                deallocate(molecule%atom)

            case ("conections")
                do i=1,MAX_ATOMS
                    if (allocated(molecule%atom(i)%connect)) deallocate(molecule%atom(i)%connect)
                    if (allocated(molecule%atom(i)%env))     deallocate(molecule%atom(i)%env)
                enddo

            case ("bond_param") 
                if (allocated(molecule%geom%bond))  deallocate(molecule%geom%bond)
                if (allocated(molecule%geom%pair))  deallocate(molecule%geom%pair)
                if (allocated(molecule%geom%angle)) deallocate(molecule%geom%angle)
                if (allocated(molecule%geom%dihed)) deallocate(molecule%geom%dihed)

            case ("all")
                if (allocated(molecule%atom))       deallocate(molecule%atom)
                if (allocated(molecule%geom%bond))  deallocate(molecule%geom%bond)
                if (allocated(molecule%geom%pair))  deallocate(molecule%geom%pair)
                if (allocated(molecule%geom%angle)) deallocate(molecule%geom%angle)
                if (allocated(molecule%geom%dihed)) deallocate(molecule%geom%dihed)


            case default
                call alert_msg("fatal","Wrong use of sys_dealloc subroutine")
        end select

        return

    end subroutine mol_dealloc

    subroutine prop_dealloc(props,context)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! This subroutine is called after using molecular (system) arrays to  
        ! deallocate them. Calls specify the context where they are called in 
        ! order to control which arrays are allocated (if not all).
        !Arguments
          ! molec (str_system,inout): props 
          ! context (character(len=10))
        !==============================================================

        use structure_types
!        use gaussian_manage
        use alerts

        implicit none

        type(str_molprops),intent(inout) :: props
        character(len=*),intent(in) :: context
        !Local
        integer :: i

        select case (adjustl(context))

            case ("td_properties") 
                deallocate(props%ddip)
                deallocate(props%trans_ddip)

            case ("frequencies") 
                deallocate(props%freq)
                deallocate(props%L)
                deallocate(props%H)

            case ("all")
                if (allocated(props%ddip))       deallocate(props%ddip)
                if (allocated(props%trans_ddip)) deallocate(props%trans_ddip)
                if (allocated(props%freq))       deallocate(props%freq)
                if (allocated(props%L))          deallocate(props%L)
                if (allocated(props%H))          deallocate(props%H)

            case default
                call alert_msg("fatal","Wrong use of prop_dealloc subroutine")
        end select

        return

    end subroutine prop_dealloc


end module allocation
