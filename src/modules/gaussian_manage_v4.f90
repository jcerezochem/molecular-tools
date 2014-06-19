module gaussian_manage
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012!

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
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


    subroutine parse_summary(unt,molec,props,tune_access)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Parse summary section from a gaussia log file. Data are shown to 
        ! stdout and are stored in molecule structure. The summary is structured
        ! in sections. The data to read can be selected by tunning options.
        !Arguments
        ! unt (int;in): unit number of the log file
        ! molec (str_resmol,inout): molecule 
        ! props (str_molprops,inout) :: props -- NEW in v4
        ! tune_access (char,in): indicate if only some parts are needed
        !                        ("struct_only", "read_hess")
        !
        !==============================================================


        integer,intent(in) :: unt
        type(str_resmol),intent(inout) :: molec
        type(str_molprops),intent(inout) :: props
        character(len=*),intent(in) :: tune_access

        integer,parameter :: SECTION_LENGTH=6000, &
                             LONG_SECTION_LENGTH=600000

        !Calculation info (not stored in "system")
        character(len=20) :: jobtype, method, basis, formula
        integer :: scan_steps, q_tot, mult, geom_style

        !Lookup auxiliar variables
        character(len=240) :: line, subline
        character(len=20) :: property
        character(len=50000) :: value !Can be as long as 3Natoms (dipole derivatives)
        !Log sections can be huge if the hessian is read
        !Thus, we allocate a (huge) character 
        character(len=LONG_SECTION_LENGTH) :: log_section

        !Counters and dummies (TODO: dummies module)
        character(len=50) :: dummy_char
        integer :: dummy_int
        character(len=1) :: sep="-"
        integer :: n_elem
        character(len=20),dimension(1:5) :: dummy_char_array
        integer :: i,j, IOstatus


        !The summary is read section by section
        line=""
        log_section=""

        !Locate summary
        i=0
        do
            i=i+1
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) then
                write(dummy_char,*) i
                call alert_msg("fatal","Unexpected end before reaching summary section. Line "//dummy_char)
            endif
            if ( INDEX(line,"GINC") /= 0 ) exit
        enddo

        !---------------------------------------
        !SECTION 1: job type
        !---------------------------------------
        do
            log_section=adjustl(trim(log_section))//line
            if ( INDEX(log_section,'\\') /= 0 ) then
                call split_line(log_section,'\\',log_section,line)
                exit
            endif
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file in log section 1")
        enddo

        !Parse section
        call split_line(log_section,'\',dummy_char,log_section) ! 1
        call split_line(log_section,'\',dummy_char,log_section) ! 1
        call split_line(log_section,'\',dummy_char,log_section) ! GINC_[SERVER]
        !Get the important:
        call split_line(log_section,'\',molec%job%type,log_section)    ! JOBTYPE
        jobtype=adjustl(molec%job%type)
        call split_line(log_section,'\',molec%job%method,log_section)     ! METHOD
        call split_line(log_section,'\',molec%job%basis,log_section)      ! BASIS
        call split_line(log_section,'\',formula,log_section)    ! FORMULA
        !---
        call split_line(log_section,'\',dummy_char,log_section) ! USER
        call split_line(log_section,'\',dummy_char,log_section) ! DATE
        call split_line(log_section,'\',dummy_char,log_section) ! GEOM-STYLE
        read(dummy_char,*) geom_style
        !geom_style might be 0(cartesian) or 1(Z-matrix)

        !.....................
        !Refine method info (this shuld be a separate subroutine)
        !.....................
        ! It can be a dash separated list: TD-B3LYP-FC...
        call string2vector_char(molec%job%method,dummy_char_array,n_elem,sep)

        if (n_elem > 1) then
            molec%job%method = ""
            sep=""
            do i=1,n_elem
                !Remove the FC card (frozen core) and separate TD
                if (adjustl(dummy_char_array(i))=="FC") then
                    dummy_char_array(i) = ""
                    sep=""
                endif

                !Construct method
                molec%job%method = trim(adjustl(molec%job%method))//sep//&
                             trim(adjustl(dummy_char_array(i)))

                !Set the next separator properly
                if (adjustl(dummy_char_array(i))=="TD" .or.&
                    adjustl(dummy_char_array(i))=="RTD".or.&
                    adjustl(dummy_char_array(i))=="UTD") then
                    !space-separate TD instruction (could also be a comma)
                    sep=" "
                else
                    !separate by dash (e.g. CAM-B3LYP)
                    sep="-"
                endif            
            enddo
        endif
        !.....................
        
        log_section=""


        !---------------------------------------
        !SECTION2: commands card
        !---------------------------------------
        do
            log_section=adjustl(trim(log_section))//line
            if ( INDEX(log_section,'\\') /= 0 ) then
                call split_line(log_section,'\\',log_section,line)
                exit
            endif
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file in log section 2")
        enddo
        !We don't care about the command line
        log_section=""


        !---------------------------------------
        !SECTION3: title card
        !---------------------------------------
        do
            log_section=adjustl(trim(log_section))//line
            if ( INDEX(log_section,'\\') /= 0 ) then
                call split_line(log_section,'\\',log_section,line)
                exit
            endif
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file in log section 3")
        enddo
        molec%title=trim(adjustl(log_section))
        log_section=""


        !---------------------------------------
        !SECTION4: molecule specifications
        !---------------------------------------
        do
            log_section=adjustl(trim(log_section))//line
            if ( INDEX(log_section,'\\') /= 0 ) then
                call split_line(log_section,'\\',log_section,line)
                exit
            endif
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file in log section 4")
        enddo

        !Get data from the section:
        ! charge and multiplicity
        call split_line(log_section,'\',subline,log_section)
        read(subline,*) q_tot, mult
        ! molecular geometry: note Z-matix need to be supported
        if (geom_style == 0) then
            i=0
            do while ( len_trim(log_section) /= 0 )
                call split_line(log_section,'\',subline,log_section)
                i=i+1
                !The format depends on the type of job (apparently): With SP includes a "layer=0"
                if ( adjustl(jobtype) == "SP" ) then
                    read(subline,*) molec%atom(i)%name, dummy_int, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z
!   !                 read(subline,*) molec%atom(i)%name, dummy_int, molec%atom(i)%R(1), molec%atom(i)%R(2), molec%atom(i)%R(3)
                else
                    read(subline,*) molec%atom(i)%name, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z
!   !                 read(subline,*) molec%atom(i)%name, molec%atom(i)%R(1), molec%atom(i)%R(2), molec%atom(i)%R(3)
                endif
            enddo
            molec%natoms=i
        else
                print*, "CAUTION: Only cartesian geometry currently supported for summary reading"
                print*, "         Coordinates will not be read"
            i=0
            do while ( len_trim(log_section) /= 0 )
                call split_line(log_section,'\',subline,log_section)
                i=i+1
                !Only taking the atomnames, for the moment...
                if ( adjustl(jobtype) == "SP" ) then
                    read(subline,*) molec%atom(i)%name
!   !                 read(subline,*) molec%atom(i)%name, dummy_int, molec%atom(i)%R(1), molec%atom(i)%R(2), molec%atom(i)%R(3)
                else
                    read(subline,*) molec%atom(i)%name
!   !                 read(subline,*) molec%atom(i)%name, molec%atom(i)%R(1), molec%atom(i)%R(2), molec%atom(i)%R(3)
                endif
            enddo
            molec%natoms=i
        endif
        log_section=""

        if ( trim(adjustl(tune_access)) == "struct_only" ) then
            rewind(unt)
            return
        endif


        !---------------------------------------
        !SECTION5: electronic properties
        !---------------------------------------
        do
            log_section=adjustl(trim(log_section))//line
            if ( INDEX(log_section,'\\') /= 0 ) then
                call split_line(log_section,'\\',log_section,line)
                exit
            endif
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file reading lof section 5")
        enddo

        !Section 5 contains assimentem with "=" as Property=Value
        !The number of items depends on the type of calculation. Let's
        !do it as general as possible
        !
        do while ( len_trim(log_section) /= 0 )
            call split_line(log_section,'\',subline,log_section)
            call split_line(subline,'=',property,value)
            !TODO: cover all possible cases
            select case ( adjustl(property) )
                case ( "Version" )
                    write(0,*)  "Program version is "//trim(adjustl(value))
                case ( "State" )
                    write(0,*)  "State of the system is "//trim(adjustl(value))
                case ( "HF" )
                    write(0,*)  "HF energy of the system is "//trim(adjustl(value))
                    if ( trim(adjustl(jobtype)) == "Scan" ) then
                        call string2vector_gm(value,props%energy_scan,scan_steps)
                        props%energy=props%energy_scan(1)
                    else 
                        read(value,*) props%energy
                    endif
                case ( "MP2" )
                    write(0,*)  "Energy with MP2 correlation of the system is "//trim(adjustl(value))
                    if ( trim(adjustl(jobtype)) == "Scan" ) then
                        call string2vector_gm(value,props%energy_scan,scan_steps)
                        props%energy=props%energy_scan(1)
                    else 
                        read(value,*) props%energy
                    endif
                case ( "MP3" )
                    write(0,*)  "Energy with MP3 correlation of the system is "//trim(adjustl(value))
                    if ( trim(adjustl(jobtype)) == "Scan" ) then
                        call string2vector_gm(value,props%energy_scan,scan_steps)
                        props%energy=props%energy_scan(1)
                    else 
                        read(value,*) props%energy
                    endif
                case ( "MP4" )
                    if ( trim(adjustl(jobtype)) == "Scan" ) then
                        call string2vector_gm(value,props%energy_scan,scan_steps)
                        props%energy=props%energy_scan(1)
                    else 
                        read(value,*) props%energy
                    endif
                case ( "RMSD" )
                    write(0,*)  "RMSD of the system is "//trim(adjustl(value))
                case ( "RMSF" )
                    write(0,*)  "RMSF of the system is "//trim(adjustl(value))
                case ( "Dipole" )
                    write(0,*)  "Dipole of the system is "//trim(adjustl(value))
                    read(value,*) props%dip(1:3)
                case ( "Quadrupole" )
                    write(0,*)  "Quadrupole of the system is "//trim(adjustl(value))
                case ("PG" )
                    write(0,*)  "Point Group of the system is "//trim(adjustl(value))
                !The following come with TD calculations
                case ("QuadrupoleDeriv" )
                    write(0,*)  "Quadrupole Derivatives present but not read"!//trim(adjustl(value))
                !The following ones correspond to a frequency job:
                case ("ZeroPoint" )
                    write(0,*)  "ZeroPoint of the system is "//trim(adjustl(value))
                case ("Thermal" )
                    write(0,*)  "Thermal energy of the system is "//trim(adjustl(value))
                case ("DipoleDeriv" )
                    write(0,*)  "Dipole Derivatives (not read)"! read in but not shown"!//trim(adjustl(value))
!                     read(value,*) props%ddip(1:3*molec%natoms,1:3)
                case ("Polar" )
                    write(0,*)  "Polarizability of the system is "//trim(adjustl(value))
                case ("PolarDeriv" )
                    write(0,*)  "Polarizability Derivatives present but not read"!//trim(adjustl(value))
                case ("HyperPolar" )
                    write(0,*)  "Hyperpolarizability of the system is "//trim(adjustl(value))
                case ("NImag" )
                    write(0,*)  "The number of imaginary frequencies is "//trim(adjustl(value))
                case default
                    call alert_msg( "fatal","Unknown item in Gaussian log section 5: "//trim(adjustl(property)) )
            end select
        enddo
        log_section=""

        ! J O B S   W I T H O U T  F R E Q U E N C I E S: leave now
        if (trim(adjustl(line)) == "@" ) then
            rewind(unt)
            return
        endif
        !Hessian is not read by default (very ineficient algorithm makes the IO proccess too costly)
        if ( trim(adjustl(tune_access)) /= "read_hess" ) then
            rewind(unt)
            return
        endif


        ! J O B S   W I T H   F R E Q U E N C I E S: continue
        !---------------------------------------
        !SECTION6: Hessian (upper triangular form)
        !---------------------------------------
        do
            log_section=adjustl(trim(log_section))//line
            if ( INDEX(log_section,'\\') /= 0 ) then
                call split_line(log_section,'\\',log_section,line)
                exit
            endif
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file reading log section 6")
        enddo

        !Read the hessian, here:
        i=0
        do while ( len_trim(log_section) /= 0 )
            i=i+1
            call split_line(log_section,',',value,log_section)
            read(value,*) props%H(i)
        enddo
        write(0,*) i, "hessian elements read in lower triangular form"
        log_section=""


        !---------------------------------------
        !SECTION7: Gradient
        !---------------------------------------
        do
            log_section=adjustl(trim(log_section))//line
            if ( INDEX(log_section,'\\') /= 0 ) then
                call split_line(log_section,'\\',log_section,line)
                exit
            endif
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file")
        enddo
        !Read-in the gradient (needs structure_types version 3)
        !
        i=0
        do while ( len_trim(log_section) /= 0 )
            i=i+1
            call split_line(log_section,',',value,log_section)
            read(value,*) props%grad(i)
        enddo
        write(0,*) i, "Gradient elements read"
        rewind(unt)

        return

        contains

        subroutine string2vector_gm(raw_vector,array_vector,n_elem)

            !Description
            ! Tranforms a string of comma sepparated values into an
            ! array of such real vaues 

            character(len=*),intent(inout) :: raw_vector
#ifdef DOUBLE
            double precision,dimension(:),intent(out) :: array_vector
#else
            real,dimension(:),intent(out) :: array_vector
#endif
            integer,intent(out) :: n_elem

            !Local
            character(len=240) :: auxchar
            integer :: i
        
            
            !Read unknown length vector (comma sepparated)
            i=0
            do 
                i=i+1
                if ( INDEX(raw_vector,',') /= 0 ) then
                    call split_line(raw_vector,',',auxchar,raw_vector)
                    read(auxchar,*) array_vector(i)
                else 
                    read(raw_vector,*) array_vector(i)
                    exit
                endif
            enddo  
            n_elem=i

            return

        end subroutine string2vector_gm

    end subroutine parse_summary


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

