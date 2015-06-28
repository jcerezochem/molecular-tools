module gaussian_fchk_manage_core

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information 
    !  from gaussian fchk files:
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
    ! Jan2015 : added write_fchk subroutine (a generic fchk section writer)
    ! ******************************************************************************
    ! TODO: Build a module only with read_fchk and a generic (section oriented)
    !       write_fchk subroutine. Then, the oder "relevant" subroutines should
    !       be included in a "extended" fchk module if necesary
    ! ******************************************************************************
    !
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

    subroutine read_fchk(unt,section,data_type,N,A,I,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Generic SR to read any section of the checkpoint
        ! Enter deallocated arrays
        !Arguments
        ! unt (int;in): unit number of the log file
        ! section(char,in): name of the section to be read
        ! data_type(char,out): Integer (I) or Real (R) data read
        ! N(int,in): Number of elements to be read
        ! A(real,dimension(:)): Real array to store real data
        ! I(integer,dimension(:)): Int array to store int data
        ! error_flag(integer,out): 0: success
        !                          1: section not found
        !==============================================================

        integer,intent(in) :: unt
        character(len=*),intent(in) :: section
        character(len=1),intent(out) :: data_type
        integer,intent(out) :: N
#ifdef DOUBLE
        double precision, dimension(:), allocatable, intent(out) :: A
#else
        real, dimension(:), allocatable, intent(out) :: A
#endif
        integer,dimension(:), allocatable, intent(out) :: I
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=42) :: section_full
        character(len=1) :: is_array
        character(len=40) :: cdata
        !I/O
        integer :: IOstatus
        
        
        ! Search section
        error_flag = 0
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    call alert_msg("note","Section '"//trim(adjustl(section))//"' not present in the FCHK file.")
                    error_flag=1
                    rewind(unt)
                    return
                endif
                ! 2) Found what looked for!      
                if ( INDEX(line,trim(adjustl(section))) /= 0 ) then
                    read(line,'(A42)') section_full
                    if (adjustl(section_full) == adjustl(section)) exit
                endif
        enddo

        !Get info about section from line just read
        read(line,'(A42,X,A1,3X,A1,X,A)') section_full, data_type, is_array, cdata
        if (is_array /= "N") then
            !Is not an array
            N=1
            if ( data_type == "R" ) then 
                allocate( A(1:1) )
                read(cdata,*) A(1)
            elseif ( data_type == "I" ) then
                allocate( I(1:1) )
                read(cdata,*) I(1) 
            endif
        else
            read(cdata,*) N
            if ( data_type == "R" ) then
                allocate( A(1:N) )
                read(unt,*) A(1:N)
            elseif ( data_type == "I" ) then
                allocate( I(1:N) )
                read(unt,*) I(1:N)
            endif
        endif 

        rewind(unt)
        return

    end subroutine read_fchk


end module gaussian_fchk_manage_core

