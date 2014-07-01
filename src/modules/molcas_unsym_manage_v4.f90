module molcas_UnSym_manage

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information 
    !  from MOLCAS UnSym files:
    !    
    ! Note: all subroutines rewind the file after using it
    ! Notes
    !  Since variables molec, residue... are allocatable, they should
    !  be passed to the sr even if they are output to get the allocation.
    !==============================================================
    !History
    ! V2.1: add get_jobtype_fchk(I_FCHK,job,method,basis) SR
    !Version 4: uses structure_types v4
    !   note that only one subroutine uses structure_types: read_fchk_geom
    !
    ! History
    ! 27/02/14: get_jobtype_fchk now takes molec(str_resmol) as input
    !           using the new molec%job attribute

    !Common declarations:
    !===================
!     use structure_types
    use alerts
!     use line_preprocess
    implicit none
    !Practically zero:
        !Practically zero charge
#ifdef DOUBLE
        double precision,parameter :: ZERO = 0.0d0, &
                                      PREC=1.d-11
#else
        real,parameter :: ZERO = 0.0e0, &
                          PREC=1.e-11
#endif



    contains

    subroutine read_molcas_hess(unt,Npert,hess,error_flag)


        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        !
        !Arguments
        !==============================================================

        use constants

        integer,intent(in) :: unt
        integer,intent(out) :: Npert
#ifdef DOUBLE
        double precision, dimension(:,:), intent(out) :: hess
#else
        real, dimension(:,:), intent(out) :: hess
#endif
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i
        
        
        ! Search section
        error_flag = 0
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    call alert_msg("note","End of file but Hessian not read")
                    error_flag=1
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( adjustl(line) == "*BEGIN HESSIAN" ) then
                    read(unt,'(A16,I20)') line, Npert
                    exit
                endif
        enddo

        do i=1,Npert
            read(unt,*) line
            read(unt,*) hess(1:Npert,i)
        enddo
!         hess(1:Npert,1:Npert) = hess(1:Npert,1:Npert)*BOHRtoANGS**2

        rewind(unt)
        return

    end subroutine read_molcas_hess

    subroutine read_molcas_geom(unt,molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        !
        !Arguments
        !==============================================================

        use constants
        use structure_types

        integer,intent(in) :: unt
        type(str_resmol),intent(inout) :: molec

        !Local stuff
        !=============
        character(len=240) :: line=""
        !I/O
        integer :: IOstatus
        !Counters
        integer :: i
        
        
        ! Search section
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    call alert_msg("fatal","End of file but Coordinates not read")

                endif
                ! 2) Found what looked for!      
                if ( adjustl(line) == "*BEGIN COORDINATES" ) then
                    read(unt,'(A16,I20)') line !Encabezado (¿Cambia la info segun el cálculo? Si es así el encabezado es importante
                    exit
                endif
        enddo
       
        i=0
        do
            read(unt,'(A)',IOSTAT=IOstatus) line
            if ( IOstatus /= 0 ) call alert_msg("fatal","while reading coordinates from UnSym")
            if (adjustl(line) == "*END COORDINATES") exit
            i=i+1
            read(line,*) molec%atom(i)%name, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z, molec%atom(i)%atnum
        enddo
        molec%natoms=i
        molec%atom(1:i)%x = molec%atom(1:i)%x*BOHRtoANGS
        molec%atom(1:i)%y = molec%atom(1:i)%y*BOHRtoANGS
        molec%atom(1:i)%z = molec%atom(1:i)%z*BOHRtoANGS

        rewind(unt)
        return

    end subroutine read_molcas_geom


end module molcas_UnSym_manage

