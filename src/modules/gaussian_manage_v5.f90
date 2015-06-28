module gaussian_manage
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012!

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information 
    !  from gaussian log files:
    !    get_std_geom(unt,molec)
    !    read_charge(unt,molec,q_type)
    !    parse_summary(unt,molec,tune_access)
    !    read_freq(unt,molec)
    ! Note: all subroutines rewind the file after using it
    ! Notes
    !  Since variables molec, residue... are allocatable, they should
    !  be passed to the sr even if they are output to get the allocation.
    !  Updated parse_summary to include Z-matrix reading
    !
    ! History
    ! Version 4: uses structure_types v4
    ! Version 5: Simplyfies subroutines relying on gaussian_manage_lowlevel
    !            Namely, removed parse_summary subroutine
    ! 
    !==============================================================

    !Common declarations:
    !===================
    use structure_types
    use alerts
    use line_preprocess
    implicit none
    !Practically zero:
#ifdef DOUBLE
    double precision,parameter :: ZERO = 1.D-10
#else
    real,parameter :: ZERO = 1.E-10
#endif



    contains


    subroutine get_std_geom(unt,molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the gaussian standard geometry from the last record (optimized)
        ! from the log file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! molec (str_resmol,inout): molecule 
        !==============================================================

        !Practically zero charge
#ifdef DOUBLE
        double precision,parameter :: ZERO_q = 2.d-6
#else
        real,parameter :: ZERO_q = 2.d-6
#endif

        !System variables
        type(str_resmol),intent(inout) :: molec

        !Reading stuff
        character(len=50) :: header_of_section, &
                             end_of_section
        integer :: n_useles_lines
        character(len=240) :: line=""
        logical :: final_section, found_standard_orientation

        !Auxiliar variables and Dummies
        character(len=30) :: dummy_char
        integer :: dummy_int

        !=============
        !Counters
        integer :: i,j
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================


        !Set variables
        header_of_section="Standard orientation"
        end_of_section=   "----------"
        n_useles_lines=4

        ! Take the last Standard orientation in file
        final_section=.false.
        found_standard_orientation=.false.
        do
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching Standard Orientation")
                ! 2) Final section reached
                if ( INDEX(line,"GINC") /= 0 ) then
                    final_section=.true.
                    exit
                endif
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    found_standard_orientation=.true.
                    exit
                endif
            enddo
            if (final_section) exit

            ! Overpass lines until reaching the target table
            do j=1,n_useles_lines
                read(unt,'(X,A)',IOSTAT=IOstatus) line 
            enddo

            ! Read Standard orientation to a predefined end of table (no need to know the number of atoms)
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while scanning Standard orientation")
                ! 2) End of table
                if ( INDEX(line,trim(adjustl(end_of_section))) /= 0 ) exit    
                ! 3) Table entry
                read(line,*) i, molec%atom(i)%AtNum, dummy_int, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z
!                 read(line,*) i, molec%atom(i)%AtNum, dummy_int, molec%atom(i)%R(1), molec%atom(i)%R(2), molec%atom(i)%R(3)
            enddo

        enddo

        if ( .not.found_standard_orientation ) call alert_msg("warning","No Standard orientation retrieved")

        rewind(unt)
        return

    end subroutine get_std_geom


    subroutine get_first_std_geom(unt,molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the gaussian standard geometry from the last record (optimized)
        ! from the log file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! molec (str_resmol,inout): molecule 
        !==============================================================

        !Practically zero charge
#ifdef DOUBLE
        double precision,parameter :: ZERO_q = 2.d-6
#else
        real,parameter :: ZERO_q = 2.d-6
#endif

        !System variables
        type(str_resmol),intent(inout) :: molec

        !Reading stuff
        character(len=50) :: header_of_section, &
                             end_of_section
        integer :: n_useles_lines
        character(len=240) :: line=""
        logical :: final_section, found_standard_orientation

        !Auxiliar variables and Dummies
        character(len=30) :: dummy_char
        integer :: dummy_int

        !=============
        !Counters
        integer :: i,j
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================

        character(len=5),dimension(100) :: atom_names_from_atnum
 
        !This should be elsewhere (constants_mod?)
        data atom_names_from_atnum(1:18) &
         /'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar'/


        !Set variables
        header_of_section="Standard orientation"
        end_of_section=   "----------"
        n_useles_lines=4

        ! Take the last Standard orientation in file
        final_section=.false.
        found_standard_orientation=.false.
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Three possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching Standard Orientation")
!             ! 2) Final section reached
!             if ( INDEX(line,"GINC") /= 0 ) then
!                 final_section=.true.
!                 exit
!             endif
            ! 3) Found what looked for!      
            if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                found_standard_orientation=.true.
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        do j=1,n_useles_lines
            read(unt,'(X,A)',IOSTAT=IOstatus) line 
        enddo

        ! Read Standard orientation to a predefined end of table (no need to know the number of atoms)
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Three possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while scanning Standard orientation")
            ! 2) End of table
            if ( INDEX(line,trim(adjustl(end_of_section))) /= 0 ) exit    
            ! 3) Table entry
            read(line,*) i, molec%atom(i)%AtNum, dummy_int, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z
!             read(line,*) i, molec%atom(i)%AtNum, dummy_int, molec%atom(i)%R(1), molec%atom(i)%R(2), molec%atom(i)%R(3)
        enddo
        molec%natoms = i

        do i=1,molec%natoms
            molec%atom(i)%name = atom_names_from_atnum(molec%atom(i)%AtNum)
        enddo

        if ( .not.found_standard_orientation ) call alert_msg("warning","No Standard orientation retrieved")

        rewind(unt)
        return

    end subroutine get_first_std_geom

    subroutine get_first_inp_geom(unt,molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the gaussian standard geometry from the last record (optimized)
        ! from the log file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! molec (str_resmol,inout): molecule 
        !==============================================================

        !Practically zero charge
#ifdef DOUBLE
        double precision,parameter :: ZERO_q = 2.d-6
#else
        real,parameter :: ZERO_q = 2.d-6
#endif

        !System variables
        type(str_resmol),intent(inout) :: molec

        !Reading stuff
        character(len=50) :: header_of_section, &
                             end_of_section
        integer :: n_useles_lines
        character(len=240) :: line=""
        logical :: final_section, found_standard_orientation

        !Auxiliar variables and Dummies
        character(len=30) :: dummy_char
        integer :: dummy_int

        !=============
        !Counters
        integer :: i,j
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================

        character(len=5),dimension(100) :: atom_names_from_atnum
 
        !This should be elsewhere (constants_mod?)
        data atom_names_from_atnum(1:18) &
         /'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar'/


        !Set variables
        header_of_section="Input orientation"
        end_of_section=   "----------"
        n_useles_lines=4

        ! Take the last Standard orientation in file
        final_section=.false.
        found_standard_orientation=.false.
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Three possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching Standard Orientation")
!             ! 2) Final section reached
!             if ( INDEX(line,"GINC") /= 0 ) then
!                 final_section=.true.
!                 exit
!             endif
            ! 3) Found what looked for!      
            if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                found_standard_orientation=.true.
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        do j=1,n_useles_lines
            read(unt,'(X,A)',IOSTAT=IOstatus) line 
        enddo

        ! Read Standard orientation to a predefined end of table (no need to know the number of atoms)
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Three possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while scanning Standard orientation")
            ! 2) End of table
            if ( INDEX(line,trim(adjustl(end_of_section))) /= 0 ) exit    
            ! 3) Table entry
            read(line,*) i, molec%atom(i)%AtNum, dummy_int, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z
!             read(line,*) i, molec%atom(i)%AtNum, dummy_int, molec%atom(i)%R(1), molec%atom(i)%R(2), molec%atom(i)%R(3)
        enddo
        molec%natoms = i

        do i=1,molec%natoms
            molec%atom(i)%name = atom_names_from_atnum(molec%atom(i)%AtNum)
        enddo

        if ( .not.found_standard_orientation ) call alert_msg("warning","No Standard orientation retrieved")

        rewind(unt)
        return

    end subroutine get_first_inp_geom

    subroutine get_next_std_geom(unt,molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the NEXT gaussian standard geometry from the line where
        ! is called, while reading a log file. DO NOT REWIND
        !Arguments
        ! unt (int;in): unit number of the log file
        ! molec (str_resmol,inout): molecule 
        !==============================================================

        !Practically zero charge
#ifdef DOUBLE
        double precision,parameter :: ZERO_q = 2.d-6
#else
        real,parameter :: ZERO_q = 2.d-6
#endif

        !System variables
        type(str_resmol),intent(inout) :: molec

        !Reading stuff
        character(len=50) :: header_of_section, &
                             header_of_section2,&
                             end_of_section
        integer :: n_useles_lines
        character(len=240) :: line=""
        logical :: final_section, found_standard_orientation

        !Auxiliar variables and Dummies
        character(len=30) :: dummy_char
        integer :: dummy_int

        !=============
        !Counters
        integer :: i,j
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================

        character(len=5),dimension(104) :: atom_names_from_atnum
 
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


        !Set variables
        header_of_section="Standard orientation"
        header_of_section2="Input orientation"
        end_of_section=   "----------"
        n_useles_lines=4

        ! Take the last Standard orientation in file
        final_section=.false.
        found_standard_orientation=.false.
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Three possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching the next Standard Orientation")
!             ! 2) Final section reached
!             if ( INDEX(line,"GINC") /= 0 ) then
!                 final_section=.true.
!                 exit
!             endif
            ! 3) Found what looked for!      
            if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 .or. &
                 INDEX(line,trim(adjustl(header_of_section2))) /= 0) then
                found_standard_orientation=.true.
                exit
            endif
        enddo

        ! Overpass lines until reaching the target table
        do j=1,n_useles_lines
            read(unt,'(X,A)',IOSTAT=IOstatus) line 
        enddo

        ! Read Standard orientation to a predefined end of table (no need to know the number of atoms)
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Three possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while scanning Standard orientation")
            ! 2) End of table
            if ( INDEX(line,trim(adjustl(end_of_section))) /= 0 ) exit    
            ! 3) Table entry
            read(line,*) i, molec%atom(i)%AtNum, dummy_int, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z
!             read(line,*) i, molec%atom(i)%AtNum, dummy_int, molec%atom(i)%R(1), molec%atom(i)%R(2), molec%atom(i)%R(3)
        enddo
        molec%natoms = i

        do i=1,molec%natoms
            molec%atom(i)%name = atom_names_from_atnum(molec%atom(i)%AtNum)
        enddo

        if ( .not.found_standard_orientation ) call alert_msg("warning","No Standard orientation retrieved")

!       This sr does not rewid the file
!        rewind(unt)
        return

    end subroutine get_next_std_geom


    subroutine read_charge(unt,molec,q_type)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read charge from the last population analysis in a log file
        ! Charges can be selected between the electronic potential fit
        ! (e.g. Merz-Kollman) or Mulliken population analysis. This 
        ! subroutine can also be used to set all charges to zero
        !Arguments
        ! unt (int;in): unit number of the log file
        ! molec (str_resmol,inout): molecule 
        ! q_type (char,in): indicate charges to be read
        !==============================================================

        integer,intent(in) :: unt
        type(str_resmol),intent(inout) :: molec
        character(len=*),intent(in) :: q_type        

        !Lookup auxiliar variables
        character(len=50) :: header_of_section, &
                             end_of_section, &
                             end_of_section2
        character(len=240) :: line
        integer :: n_useles_lines
#ifdef DOUBLE
        double precision :: q_tot
#else
        real :: q_tot
#endif

        !Logical switchs
        logical :: found_charge_table, final_section

        !Counters and dummies (TODO: dummies module)
        character(len=50) :: dummy_char
        integer :: i,j, IOstatus


        if (q_type == "mulliken") then
            ! MULLIKEN
            header_of_section="Mulliken atomic charges"
            end_of_section="Sum of Mulliken"
            end_of_section2="Sum of Mulliken"
            n_useles_lines = 1
        elseif (q_type == "elec_pot") then
            ! ELECTRIC POTENTIAL DERIVED CHARGES
            header_of_section="Charges from ESP fit"
            end_of_section="---------------------"
            end_of_section2="Sum of ESP charges"
            n_useles_lines = 2
        elseif (q_type == "all_zero") then
            molec%atom(:)%q = 0
            rewind(unt)
            return
        else
            call alert_msg("warning","Unrecognize q_type: "//q_type//" All charges set to zero")
            molec%atom(:)%q = 0
            rewind(unt)
            return
        endif
         
        q_tot=0.
        final_section=.false.
        found_charge_table = .false.
        do

            ! Locate charge section identifier
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file")
                ! 2) Final section reached
                if ( INDEX(line,"GINC") /= 0 ) then
                    final_section=.true.
                    exit
                endif
                ! 3) Found a charge section
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    found_charge_table = .true.
                    exit
                endif

            enddo
            !In case 2) leave the main loop
            if ( final_section ) exit

            ! Overpass lines until reaching the target table
            do j=1,n_useles_lines
                read(unt,'(X,A)',IOSTAT=IOstatus) line 
            enddo

            ! Read charge table to a predefined end of table
            q_tot=0.
            do
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file:
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while reading charges")
                ! 2) End of table:
                if ( INDEX(line,trim(adjustl(end_of_section))) /= 0 .or. &
                     INDEX(line,trim(adjustl(end_of_section2))) /= 0) exit    
                ! 3) Table entry:
                read(line,*) i, dummy_char, molec%atom(i)%q
                ! Monitor total charge:
                q_tot=q_tot+molec%atom(i)%q
            enddo

        enddo

        if (.not.found_charge_table) then
            call alert_msg( "error", 'Charge table not found for q_type "'//trim(adjustl(q_type))//&
                                     '", all charges set to zero' )  
            molec%atom(:)%q = 0
            rewind(unt)
            return 
        endif

        if ( q_tot > ZERO ) then
             write(dummy_char,*) q_tot
             call alert_msg("note", "Non-zero total charge ("//trim(adjustl(dummy_char))//")" )
        endif

    end subroutine read_charge

    subroutine read_charge_new(unt,molec,dens_type,q_type)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read charge from the last population analysis in a log file
        ! Charges can be selected between the electronic potential fit
        ! (e.g. Merz-Kollman) or Mulliken population analysis. This 
        ! subroutine can also be used to set all charges to zero
        !Arguments
        ! unt (int;in): unit number of the log file
        ! molec (str_resmol,inout): molecule 
        ! dens_type (char,in): type of density (SCF, CI)
        ! q_type (char,in): indicate charges to be read
        !==============================================================

        integer,intent(in) :: unt
        type(str_resmol),intent(inout) :: molec
        character(len=*),intent(in) :: q_type, dens_type      

        !Lookup auxiliar variables
        character(len=50) :: header_of_section, &
                             header, &
                             end_of_section, &
                             end_of_section2
        character(len=240) :: line, pop_header
        integer :: n_useles_lines
#ifdef DOUBLE
        double precision :: q_tot
#else
        real :: q_tot
#endif

        !Logical switchs
        logical :: found_charge_table, final_section

        !Counters and dummies (TODO: dummies module)
        character(len=50) :: dummy_char
        integer :: i,j, IOstatus


        if (adjustl(dens_type) == "SCF") then
            pop_header = "Population analysis using the SCF density."
        elseif (adjustl(dens_type) == "CI") then
            pop_header = "Population analysis using the CI density."
        elseif (adjustl(dens_type) == "SAC") then
            pop_header = "Population analysis using the QCI/CC density."
        else   
            call alert_msg("fatal","Unsupported density"//trim(adjustl(dens_type)))
        endif

        if (q_type == "mulliken") then
            ! MULLIKEN
            header_of_section="Mulliken charges"
            end_of_section="Sum of Mulliken"
            end_of_section2="Sum of Mulliken"
            n_useles_lines = 1
        elseif (q_type == "elec_pot") then
            ! ELECTRIC POTENTIAL DERIVED CHARGES
            header_of_section="Charges from ESP fit"
            end_of_section="---------------------"
            end_of_section2="Sum of ESP charges"
            n_useles_lines = 2
        elseif (q_type == "all_zero") then
            molec%atom(:)%q = 0
            rewind(unt)
            return
        else
            call alert_msg("warning","Unrecognize q_type: "//q_type//" All charges set to zero")
            molec%atom(:)%q = 0
            rewind(unt)
            return
        endif
         
        q_tot=0.
        final_section=.false.
        found_charge_table = .false.
        do

            ! Locate charge section identifier
            header="X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X_X"
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file")
                ! 2) Final section reached
                if ( INDEX(line,"GINC") /= 0 ) then
                    final_section=.true.
                    exit
                endif
                ! 3a) Found a charge section
                if ( INDEX(line,trim(adjustl(pop_header))) /= 0 ) then
                    header=header_of_section
print*, "Found", pop_header
                endif
                ! 3b) Found a charge section
                if ( INDEX(line,trim(adjustl(header))) /= 0 ) then
                    found_charge_table = .true.
                    exit
                endif

            enddo
            !In case 2) leave the main loop
            if ( final_section ) exit

            ! Overpass lines until reaching the target table
            do j=1,n_useles_lines
                read(unt,'(X,A)',IOSTAT=IOstatus) line 
            enddo

            ! Read charge table to a predefined end of table
            q_tot=0.
            do
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file:
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while reading charges")
                ! 2) End of table:
                if ( INDEX(line,trim(adjustl(end_of_section))) /= 0 .or. &
                     INDEX(line,trim(adjustl(end_of_section2))) /= 0) exit    
                ! 3) Table entry:
                read(line,*) i, dummy_char, molec%atom(i)%q
                ! Monitor total charge:
                q_tot=q_tot+molec%atom(i)%q
            enddo

        enddo

        if (.not.found_charge_table) then
            call alert_msg( "error", 'Charge table not found for q_type "'//trim(adjustl(q_type))//&
                                     '", all charges set to zero' )  
            molec%atom(:)%q = 0
            rewind(unt)
            return 
        endif

        if ( q_tot > ZERO ) then
             write(dummy_char,*) q_tot
             call alert_msg("note", "Non-zero total charge ("//trim(adjustl(dummy_char))//")" )
        endif

    end subroutine read_charge_new


    subroutine read_gauslog_geom(unt,molec)

        !TODO: Add Z-matrix to be read

        use gaussian_manage_lowlevel

        integer,intent(in)               :: unt
        type(str_molprops),intent(inout) :: molec

        !Local
        integer :: error, nitems, i
        character(len=100000) :: string 
        character(len=500)    :: geom_char 
#ifdef DOUBLE
        double precision,dimension(4) :: dummy_real
#else
        real,dimension(5)  :: dummy_real
#endif

        ! Natoms
        call read_gausslog_natoms(10,molec%natoms,error)
        if (error /= 0) call alert_msg("fatal","Reading natoms in read_gauslog_geom")

        !Geom
        call summary_parser(10,4,string,error)
        if (error /= 0) call alert_msg("fatal","Reading geom in read_gauslog_geom")
        !Throw "charge mult" away
        call split_line(string,'\',geom_char,string)

        !Read the first geometry entry
        call split_line(string,'\',geom_char,string)
        !Get the number of elements
        call string2vector(geom_char,dummy_real,nitems)
        !The decide the format
        if (nitems == 4) then
            !Read first atom
            read(geom_char,*) molec%atom(1)%name, &
                              molec%atom(1)%x,    &
                              molec%atom(1)%y,    &
                              molec%atom(1)%z
            do i=2,molec%natoms
                call split_line(string,'\',geom_char,string)
                read(geom_char,*) molec%atom(i)%name, &
                                  molec%atom(i)%x,    &
                                  molec%atom(i)%y,    &
                                  molec%atom(i)%z
            enddo
        elseif (nitems == 5) then
            !Read first atom
            read(geom_char,*) molec%atom(1)%name, &
                              i,                  &
                              molec%atom(1)%x,    &
                              molec%atom(1)%y,    &
                              molec%atom(1)%z
            do i=2,molec%natoms
                call split_line(string,'\',geom_char,string)
                read(geom_char,*) molec%atom(i)%name, &
                                  i,                  &
                                  molec%atom(i)%x,    &
                                  molec%atom(i)%y,    &
                                  molec%atom(i)%z
            enddo
        elseif (nitems == 1) then
            !Read Z-mat...
            call alert_msg("fatal","Reading Z-mat still not supported in read_gauslog_geom")
        else 
            call alert_msg("fatal","Confused about the format when reading geom in read_gauslog_geom")
        endif

        return

    end subroutine read_gauslog_geom


    subroutine read_freq(unt,molec,props)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read normal mode information from Gaussian log file. It reads 
        ! frequencies and normal mode description through the cartesian
        ! L matrix (might be in low --2 demimals-- or high precission)
        !Arguments
        ! unt (int;in): unit number of the log file
        ! molec (str_resmol,inout): molecule 
        ! props (str_molprops,inout) :: props -- NEW
        !==============================================================

        integer,intent(in) :: unt
        type(str_resmol),intent(inout) :: molec
        type(str_molprops),intent(inout) :: props

        !Lookup auxiliar variables
        character(len=240) :: line, subline
        integer :: Nnm, ncol, extra_col
        character(len=4) splitter
        character(len=20) :: property, value

        ! Switches
        integer :: get_frq

        !Counters and dummies (TODO: dummies module)
        character(len=50) :: dummy_char
        integer :: dummy_int
        integer :: i,j, IOstatus


        !Supose non-linear molecule
        props%Nnm=3*molec%natoms-5
        !Check if we have HPmodes
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file")
            ! 2) Frequencies section reached
            if ( INDEX(line,"and normal coordinates") /= 0 ) exit
        enddo
        read(unt,'(X,A)',IOSTAT=IOstatus) line  ! NM counter
        read(unt,'(X,A)',IOSTAT=IOstatus) line  ! Symmetry
        read(unt,'(X,A)',IOSTAT=IOstatus) line  ! Frequency line
        if ( INDEX(line,'---') /= 0 ) then
            !High precision modes
            get_frq = 1
            ncol=5
            splitter="---"
        else
            !Low precision modes
            get_frq = 2
            ncol=3
            splitter="--"
            call alert_msg("note","No high precision modes available You may want to get them from fchk file.")
        endif
        !Place the scan back to the begining of frequency section (3 lines back)
        backspace(unt)
        backspace(unt)
        backspace(unt)

        !Read complete frequency blocks
        !------------------------------
        do i=1,Nnm/ncol
            read(unt,'(X,A)',IOSTAT=IOstatus) line  ! NM counter
            read(unt,'(X,A)',IOSTAT=IOstatus) line  ! Symmetry
            ! Following information is variable (whether Raman is activated or not)
            ! The format is:  "Property --(-) Value" 
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                if ( INDEX(line,trim(adjustl(splitter))) == 0 ) exit
                call split_line(line,trim(adjustl(splitter)),property,value)
                select case ( adjustl(property) )
                    case ("Frequencies")
                        read(value,*) props%freq(ncol*i-ncol+1:ncol*i)
                    !Open to include more data to retrieve...
                end select
            enddo
            !Now read normal modes
            if ( get_frq == 1 ) then
                do j=1,3*molec%natoms
                    read(unt,*,IOSTAT=IOstatus) dummy_int, dummy_int, dummy_int, props%L(j,ncol*i-ncol+1:ncol*i) !L
                enddo
            elseif (get_frq == 2 ) then
                do j=1,molec%natoms
                    read(unt,*,IOSTAT=IOstatus) dummy_int, dummy_int, props%L(3*j-2:3*j,ncol*i-ncol+1:ncol*i) !L
                enddo
            endif

        enddo
        i=min(i,Nnm/ncol)

        !Read incomplete frequency blocks
        !---------------------------------
        extra_col = Nnm - Nnm/ncol*ncol
        if ( extra_col /= 0 ) then
            read(unt,'(X,A)',IOSTAT=IOstatus) line      ! NM counter
            read(unt,'(X,A)',IOSTAT=IOstatus) line  ! Symmetry
            ! Folloing information is variable (whether Raman is activated or not)
            ! The format is:  "Property --(-) Value" 
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                if ( INDEX(line,trim(adjustl(splitter))) == 0 ) exit
                call split_line(line,trim(adjustl(splitter)),property,value)
                select case ( adjustl(property) )
                    case ("Frequencies")
                        read(value,*) props%freq(ncol*i+1:ncol*i+extra_col)
                    !Open to include more data to retrieve...
                end select
            enddo
            !Now read normal modes
            if ( get_frq == 1 ) then
                do j=1,3*molec%natoms
                    read(unt,*,IOSTAT=IOstatus) dummy_int, dummy_int, dummy_int, props%L(j,ncol*i+1:ncol*i+extra_col) !L
                enddo
            elseif (get_frq == 2 ) then
                do j=1,molec%natoms
                    read(unt,*,IOSTAT=IOstatus) dummy_int, dummy_int, props%L(3*j-2:3*j,ncol*i+1:ncol*i+extra_col) !L
                enddo
            endif
        endif

        rewind(unt)
        return

    end subroutine read_freq



end module gaussian_manage

