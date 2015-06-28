module gaussian_manage_lowlevel
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012!

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get information from 
    !  Gaussian log files. This is a med-low level module that
    !  only requires low-level subroutines (i.e. do not use 
    !  structure_types module) 
    ! 
    !==============================================================

    !Common declarations:
    !===================
    use line_preprocess
    implicit none

    contains

    subroutine summary_parser(unt,isect,section,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! New summary parser, that simply reads in the summary section 
        ! and passes a piece of it back to the calling program as a string.
        ! So, there is no need of any special types in the subroutine.
        ! Only line_preprocess moudule is used.
        ! The length of the summary lines is assumed to be 240
        !
        !Arguments
        ! unt     (inp) int /scalar   unit for the file 
        ! isect   (inp) int /scalar   number of the section
        !                             1: General info (host,user,date and 
        !                                job info: job type, method, basis set
        !                                molecular formula)
        !                             2: Route section
        !                             3: Job title
        !                             4: Molecular geometry (and atnames)
        !                             5: Electronic properties
        !                             6: Hessian (lower triangular form)
        !                             7: Gradient
        ! section (out) char/scalar   string with the selected section
        ! error_flag (out) flag  0 : success
        !                        1 : Requested secction not present in the summary
        !                            Returns the last section read
        !                     10+i : read error in line i
        !                       -i : the section size is i characters larger then
        !                            length of the character. Additional i character 
        !                            are required. The output includes the last
        !                            part of the section that fits the length
        !
        !==============================================================


        integer,intent(in) :: unt
        integer,intent(in) :: isect
        character(len=*),intent(out) :: section
        integer,intent(out) :: error_flag

        !Local
        ! Counters
        integer :: i, isection
        ! I/O
        integer :: IOstatus
        character(len=240) :: line

        !Locate summary
        error_flag = 0
        i=0
        do
            i=i+1
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                error_flag = i+10
                rewind(unt)
                return
            endif
            if ( INDEX(line,"GINC") /= 0 ) exit
        enddo

        !Read sections
        isection = 1
        do isection=1,7
            error_flag = 0
            section    = ""
            do
                i=i+1
                !Check if there is space for a line. Otherwise, make 
                ! space and continue
                if (len(section)-len_trim(section) <= 240) then
                    error_flag = error_flag - 240
                    section=section(241:len(section))
                endif
                !Append last line to section
                section=adjustl(trim(section))//line
                !If we reach a section end "\\", split and finish section
                if ( INDEX(section,'\\') /= 0 ) then
                    call split_line(section,'\\',section,line)
                    exit
                endif
                !Read in new line
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                if ( IOstatus < 0 ) then
                    error_flag = i+10
                    rewind(unt)
                    return
                endif
            enddo
            if (isection == isect) then
                rewind(unt)
                return
            endif
            if (line == "@") then
                error_flag = 1
                rewind(unt)
                return
            endif

        enddo

        return

    end subroutine summary_parser


    subroutine estimate_section_length(unt,isect,length,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! New summary parser, that simply reads in the summary section 
        ! and passes a piece of it back to the calling program as a string.
        ! So, there is no need of any special types in the subroutine.
        ! Only line_preprocess moudule is used.
        ! The length of the summary lines is assumed to be 240
        !
        !Arguments
        ! unt     (inp) int /scalar   unit for the file 
        ! isect   (inp) int /scalar   number of the section
        !                             1: General info (host,user,date and 
        !                                job info: job type, method, basis set
        !                                molecular formula)
        !                             2: Route section
        !                             3: Job title
        !                             4: Molecule specifications (atname,geom)
        !                             5: Electronic properties
        !                             6: Hessian (lower triangular form)
        !                             7: Gradient
        ! section (out) int /scalar   Length of the string to store the section isect
        ! error_flag (out) flag  0 : success
        !                        1 : Requested secction not present in the summary
        !                            Returns the last section read
        !                     10+i : read error in line i
        !                       -i : the section size is i characters larger then
        !                            length of the character. Additional i character 
        !                            are required. The output includes the last
        !                            part of the section that fits the length
        !
        !==============================================================


        integer,intent(in) :: unt
        integer,intent(in) :: isect
        integer,intent(out) :: length
        integer,intent(out) :: error_flag

        !Local
        ! Counters
        integer :: i, isection
        ! I/O
        integer :: IOstatus
        character          :: cnull
        character(len=240) :: line

        !Locate summary
        error_flag = 0
        i=0
        do
            i=i+1
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                error_flag = i+10
                rewind(unt)
                return
            endif
            if ( INDEX(line,"GINC") /= 0 ) exit
        enddo

        !Read sections
        isection = 1
        do isection=1,7
            length     = 0
            error_flag = 0
            do
                !Update length
                length=length+len_trim(line)

                !If we reach a section end "\\", pass to the next section
                if ( INDEX(line,'\\') /= 0 ) then
                    call split_line(line,'\\',cnull,line)
                    exit
                endif

                !Read in new line
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                if ( IOstatus < 0 ) then
                    error_flag = i+10
                    rewind(unt)
                    return
                endif
            enddo
            if (isection == isect) then
                rewind(unt)
                return
            endif
            if (line == "@") then
                error_flag = 1
                rewind(unt)
                return
            endif

        enddo

        return

    end subroutine estimate_section_length

    subroutine read_gausslog_property(unt,property,value,io_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Get electrinic property from the summary (section=5). Returs the
        ! the value in a string (so as to avoid specifying the type)
        ! USES ROUTINES IN THIS MODULE
        !
        !Arguments
        ! unt      (inp) int /scalar   unit for the file 
        ! property (inp) char/scalar   Name of the property to get
        ! value    (out) char/scalar   Value for the property
        ! io_flag  (io ) flag          On input : length of properties section to allocate
        !                                   an auxiliar character (prop_section)
        !                              On output: error_flag
        !                                   1 : Secction 5 not present in the summary
        !                                   2 : Requested property not found
        !                                   3 : Requested property not found, and 
        !                                       'prop_section' was not allocatated properly
        !                                   4 : Warning: 'prop_section' was not allocatated 
        !                                       properly and the value might be incorrect
        !                                10+i : read error in line i
        !                                  -i : 'value' length is not enough to store the value
        !                                        It should be at least len=i
        !
        !Notes:
        ! -Allocatable characters:
        !  Fortran2003 supports allocatable strings (implemented in gfortran >=4.8)
        !  This would overcome the need of explicetily allocating the
        !  auxiliar string with io_flag (but would impose a requirement to the compiler version)
        !
        !==============================================================

        integer,intent(in)           :: unt
        character(len=*),intent(in)  :: property
        character(len=*),intent(out) :: value
        integer,intent(inout)        :: io_flag
        
        !local
        character(len=io_flag)       :: prop_section
        character(len=len(property)+len(value))    :: current_property
        integer :: i,j


        call summary_parser(unt,5,prop_section,io_flag)
        if (io_flag > 0) return
        !If the error_flag<0, prop_section was not allocated properly
        !(not large enough). Why try to read the property anyway
        !since it might be present

        !Run over all properties
        do while ( len_trim(prop_section) /= 0 )
            i=len_trim(prop_section)
            call split_line(prop_section,'\',current_property,prop_section)
            !Get the value
            call split_line(current_property,'=',current_property,value)

            !Check if it is the requested property          
            if ( adjustl(current_property) == adjustl(property) ) then
                !Check if 'value' was large enough
                j = i-len_trim(prop_section)-len_trim(property)-1
                if ( j > len(value) ) then
                    io_flag=-j
                    return
                endif
                !Warning about summary_parser errors
                if ( (io_flag < 0) .and. len_trim(prop_section) == 0) io_flag=4
                return
            endif
        enddo

        if (io_flag < 0) then
            io_flag = 3
        else
            io_flag = 2
        endif

        return
        
    end subroutine read_gausslog_property


    subroutine read_gausslog_natoms(unt,natoms,io_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Get natoms from molecular formula. It is taken from section 1, 
        ! which contains
        !  1: flag? =1?
        !  2: flag? =1?
        !  3: GINC-$HOST
        !  4: Job Type (Opt, Freq...)
        !  5: Method (RB3LYP...)
        !  6: Basis set
        !  7: Molecular formula   <===
        !  8: $USER
        !  9: $date
        ! 10: flag? =0?
        ! USES ROUTINES IN THIS MODULE
        !
        !Arguments
        ! unt      (inp) int /scalar   unit for the file 
        ! natoms   (out) char/scalar   Number of atoms
        ! io_flag  (io ) flag          On input : length of properties section to allocate
        !                                   an auxiliar character (prop_section)
        !                              On output: error_flag
        !                                   1 : Secction 5 not present in the summary
        !                                   2 : Requested property not found
        !                                   3 : Requested property not found, and 
        !                                       'prop_section' was not allocatated properly
        !                                   4 : Warning: 'prop_section' was not allocatated 
        !                                       properly and the value might be incorrect
        !                                10+i : read error in line i
        !                                  -i : 'value' length is not enough to store the value
        !                                        It should be at least len=i
        !
        !Notes:
        ! -Allocatable characters:
        !  Fortran2003 supports allocatable strings (implemented in gfortran >=4.8)
        !  This would overcome the need of explicetily allocating the
        !  auxiliar strings
        !
        !==============================================================
        
        integer,intent(in)           :: unt
        integer,intent(out)          :: natoms
        integer,intent(inout)        :: io_flag
        
        !local
        character(len=500)           :: info_section
        character(len=100)           :: formula
        integer                      :: i,j,n
        ! array with number of atoms per element in the formula
        !  len set the maximum number (len=3 => 999)
        !  dimension set the maximum different elements (50 seems enough)
        character(len=3),dimension(50) :: nel_char

        call summary_parser(unt,1,info_section,io_flag)

        !Get formula (7th entry)
        do i=1,7
            call split_line(info_section,'\',formula,info_section)
        enddo

        do i=1,len_trim(formula)
            !Convert to blank if it is not a number (48==0 to 57==9)
            if (ichar(formula(i:i)) < 48 .or. &
                ichar(formula(i:i)) > 57) formula(i:i)=" "
        enddo
 
        !Use line parser to get the numbers reparatelly
        call parse_line(formula,n,nel_char)

        !Now read each one and sum them
        natoms = 0
        do i=1,n
            read(nel_char(i),*) j
            natoms = natoms + j
        enddo

        return

    end subroutine read_gausslog_natoms


end module gaussian_manage_lowlevel