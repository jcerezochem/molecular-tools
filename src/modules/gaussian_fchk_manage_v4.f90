module gaussian_fchk_manage

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

    subroutine write_fchk(unt,section,data_type,N,A,I,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Generic SR to read any section of the checkpoint
        ! Enter allocated arrays
        !Arguments
        ! unt (int;in): unit number of the log file
        ! section(char,in): name of the section to be written
        ! data_type(char,in): Integer (I) or Real (R) data write
        ! N(int,in): Number of elements to be written
        ! A(real,dimension(:)): Real array to store real data
        ! I(integer,dimension(:)): Int array to store int data
        ! error_flag(integer,out): 0: success
        !                          1: write failure
        !
        ! NOTES
        ! Should is_array be an input?
        ! Not really, we can specify that is a scalar, e.g. by setting 
        ! N=0 (so we let N=1 fro vectors with length equal 1, if this 
        ! is ever the case)
        !==============================================================

        integer,intent(in) :: unt
        character(len=*),intent(in) :: section
        character(len=1),intent(in) :: data_type
        integer,intent(in) :: N
#ifdef DOUBLE
        double precision, dimension(:), intent(in) :: A
#else
        real, dimension(:), intent(in) :: A
#endif
        integer,dimension(:), allocatable, intent(in) :: I
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=43) :: section_full
        !I/O
        integer :: IOstatus
        
        error_flag = 0
        section_full = adjustl(section)

        !If N=0, it is an scalar
        if (N == 0) then
            if ( data_type == "I" ) then 
                write(unt,'(A43,A,I17)') section_full, data_type, I(1)
            elseif ( data_type == "R" ) then
                write(unt,'(A43,A,ES27.15)') section_full, data_type, A(1)
            endif
        else
            write(unt,'(A43,A,A5,I12)') section_full, data_type,"   N=",N
            if ( data_type == "I" ) then 
                write(unt,'(6I12)') I(1:N)
            elseif ( data_type == "R" ) then
                write(unt,'(5ES16.8)') A(1:N)
            endif
        endif 

        return

    end subroutine write_fchk


    subroutine read_fchk_geom(unt,molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Kind of wrapper of read_fchk addapted to read geommetry with
        ! derived types and atom names. Coordinates are given in ANGS
        !==============================================================

        use structure_types

#ifdef DOUBLE
        double precision, parameter :: BOHRtoANGS = 5.2917720859D-1
#else
        real, parameter :: BOHRtoANGS = 5.2917720859D-1
#endif

        integer,intent(in) :: unt
        type(str_resmol),intent(inout) :: molec

        !read_fchk stuff
        character(len=1) :: dtype
        integer :: N
#ifdef DOUBLE
        double precision, dimension(:), allocatable :: A
#else
        real, dimension(:), allocatable :: A
#endif
        integer,dimension(:), allocatable :: IA
        integer :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=42) :: section_full
        character(len=1) :: is_array
        character(len=40) :: cdata
        !I/O
        integer :: i,j

        character(len=5),dimension(103) :: atom_names_from_atnum
 
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
        
        ! Natoms
        call read_fchk(unt,"Number of atoms",dtype,N,A,IA,error_flag)
        if (error_flag == 0) then
            molec%natoms = IA(1)
            deallocate(IA)
        endif

        ! Geometry
        call read_fchk(unt,"Current cartesian coordinates",dtype,N,A,IA,error_flag)
        if (error_flag == 0) then
            do i=1,N/3
                j=3*i
                molec%atom(i)%x = A(j-2)*BOHRtoANGS
                molec%atom(i)%y = A(j-1)*BOHRtoANGS
                molec%atom(i)%z = A(j)  *BOHRtoANGS
            enddo
            deallocate(A)
        endif

        ! Atom names
        call read_fchk(unt,"Atomic numbers",dtype,N,A,IA,error_flag)
        if (error_flag == 0) then
            do j=1,molec%natoms
                molec%atom(j)%atnum = IA(j)
                molec%atom(j)%name = atom_names_from_atnum(IA(j))
            enddo
            deallocate(IA)
        endif

        ! Atom masses
        call read_fchk(unt,"Real atomic weights",dtype,N,A,IA,error_flag)
        if (error_flag == 0) then
            do j=1,molec%natoms
                molec%atom(j)%mass = A(j)
            enddo
            deallocate(A)
        endif

        !Management of non-existing information
        !No symmetry available in fchk. It's indicated with XX
        molec%PG = "XX"
        !Resname, resseq
        molec%atom(:)%resseq = 1
        molec%atom(:)%resname = "UNK"

        return

    end subroutine read_fchk_geom

!===========================================================================

    subroutine sort_jobtype_fchk(unt,jobname,jobtype,method,basis,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read second line fo fchk files (job type info) and shows the info
        !Arguments
        ! unt (int;in): unit number of the log file
        ! jobname(char,out): name of the job
        ! jobtype(char,out): jobtype (Freq, SP...)
        ! method(char,out): method (B3LYP, HF...)
        ! basis(char,out): basis set
        ! error_flag(integer,out): 0: success
        !                          1: section not found
        !==============================================================

        integer,intent(in) :: unt
        character(len=*),intent(out) :: jobname, &
                                        jobtype, &
                                        method,  &
                                        basis
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        !I/O
        integer :: IOstatus

        error_flag = 0
        read(unt,'(A)',IOSTAT=IOstatus) jobname
        if ( IOstatus < 0 ) then
            call alert_msg("note","Cannot read FCHK file (line 1). Exiting subroutine")
            error_flag=1
            rewind(unt)
            return
        endif
        read(unt,*,IOSTAT=IOstatus) jobtype,method,basis
        if ( IOstatus < 0 ) then
            call alert_msg("note","Cannot read FCHK file (line 2). Exiting subroutine")
            error_flag=1
            rewind(unt)
            return
        endif

        return

    end subroutine sort_jobtype_fchk

!===========================================================================

    subroutine get_jobtype_fchk(unt,job,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read second line fo fchk files (job type info) and get the info
        ! to job structure_type
        !Arguments
        ! unt (int;in): unit number of the log file
        ! job including..
        !       job%name(char,out): name of the job
        !       job%type(char,out): jobtype (Freq, SP...)
        !       job%method(char,out): method (B3LYP, HF...)
        !       job%basis(char,out): basis set
        ! error_flag(integer,out): 0: success
        !                          1: section not found
        !==============================================================

        use structure_types
        use line_preprocess

        integer,intent(in) :: unt
        type(str_job),intent(inout) :: job
        integer,intent(out) :: error_flag

        !Local stuff
        !=============
        character(len=240) :: line=""
        character(len=1) :: sep
        integer :: n_elem, i
        character(len=100),dimension(1:5) :: dummy_char_array
        !I/O
        integer :: IOstatus

        error_flag = 0
        read(unt,'(A)',IOSTAT=IOstatus) job%title
        if ( IOstatus < 0 ) then
            call alert_msg("note","Cannot read FCHK file (line 1). Exiting subroutine")
            error_flag=1
            rewind(unt)
            return
        endif
!         read(unt,*,IOSTAT=IOstatus) job%type,job%method,job%basis
        ! We cannot make a standard read since basis includes ","
        ! and fortran would think they separate entries. So read and 
        ! split the line
        read(unt,'(A)',IOSTAT=IOstatus) line
        sep=" "
        call string2vector_char_new(line,dummy_char_array,n_elem,sep)
        if (n_elem /= 3) then
            call alert_msg("note","Cannot read job info from fchk file")
            error_flag=-1
            return
        endif
        job%type   = trim(adjustl(dummy_char_array(1)))
        job%method = trim(adjustl(dummy_char_array(2)))
        job%basis  = trim(adjustl(dummy_char_array(3)))
        if ( IOstatus < 0 ) then
            call alert_msg("note","Cannot read FCHK file (line 2). Exiting subroutine")
            error_flag=1
            rewind(unt)
            return
        endif

        !Refine method info
        ! It can be a dash separated list: TD-B3LYP-FC...
        sep="-"
        call string2vector_char(job%method,dummy_char_array,n_elem,sep)

        if (n_elem > 1) then
            job%method = ""
            sep=""
            do i=1,n_elem
                !Remove the FC card (frozen core) and separate TD
                if (adjustl(dummy_char_array(i))=="FC") then
                    dummy_char_array(i) = ""
                    sep=""
                endif

                !Construct method
                job%method = trim(adjustl(job%method))//sep//&
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

        return

    end subroutine get_jobtype_fchk

!===========================================================================

    subroutine AO_atom_map(unt,AOmap)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the mapping between atomic orbitals and atoms from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! AOmap (int,dim(:);out): mapping matrix
        !==============================================================

        integer,dimension(:),intent(out) :: AOmap

        integer,dimension(5000) :: ShellType,AO_shells_map

        !Reading stuff
        character(len=50) :: header_of_section
        character(len=240) :: line=""

        !Auxiliar variables and Dummies
!         character(len=30) :: dummy_char
!         integer :: dummy_int

        !=============
        !Counters
        integer :: i,k
        integer :: Nshells, Shell_size
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================



        !Read shell types
        header_of_section="Shell types"
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching shell types")
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') Nshells
                    exit
                endif
        enddo

        read(unt,*) ShellType(1:Nshells)

        !Read mapping
        header_of_section="Shell to atom map"
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching AO mappings")
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') Nshells
                    exit
                endif
        enddo

        read(unt,*) AO_shells_map(1:Nshells)

        k=0
        do i=1,Nshells
            !Translate shell type to shell size
            select case (ShellType(i))
             case (0)
                 ! S
                 Shell_size = 1
             case (1)
                 ! P
                 Shell_size = 3
             case (-1)
                 ! SP
                 Shell_size = 4
             case (2)
                 ! D (cartesian)
                 Shell_size = 6
             case (-2)
                 ! D (pure)
                 Shell_size = 5
             case (3)
                 ! F (cartesian)
                 Shell_size = 10
             case (-3)
                 ! F (pure)
                 Shell_size = 7
             case default
                 Shell_size = 2*iabs(ShellType(i)) + 1
            end select

            AOmap(k+1:k+Shell_size) = AO_shells_map(i)
            k = k + Shell_size

        enddo

        rewind(unt)

        return   

    end subroutine AO_atom_map

!===========================================================================

    subroutine write_fchk_geom(unt,system)

        !Write fchk file

        use structure_types
        use constants

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system

        !local
        integer::i, j
!         integer,dimension(:), allocatable :: aux_int
#ifdef DOUBLE
        double precision, dimension(:), allocatable :: aux 
#else
        real, dimension(:), allocatable :: aux 
#endif
        character :: vartype
        character(len=43) :: section
        
        !Title
         write(unt,'(A)') "File converted to fchk. Title:"//trim(adjustl(system%title))//&
                                                 " Name:"//trim(adjustl(system%name))

        !Number of atoms
        section="Number of atoms"
        vartype="I"
        write(unt,'(A43,A,I17)') section, vartype, system%natoms

        !Atomic numbers
        section="Atomic numbers"
        vartype="I"
        j=system%natoms
        write(unt,'(A43,A,A5,I12)') section, vartype,"   N=",j
        write(unt,'(6I12)') system%atom(1:j)%AtNum

        !Nuclear charges
        section="Nuclear charges"
        vartype="R"
        j=system%natoms
        write(unt,'(A43,A,A5,I12)') section, vartype,"   N=",j
        write(unt,'(5ES16.8)') float(system%atom(1:j)%AtNum)

        !Cartesian coordinates
        allocate( aux(1:3*system%natoms) )
        section="Current cartesian coordinates"
        vartype="R"
        j=0
        do i=1,system%natoms
            j=j+1
            aux(j)=system%atom(i)%x/BOHRtoANGS
            j=j+1
            aux(j)=system%atom(i)%y/BOHRtoANGS
            j=j+1
            aux(j)=system%atom(i)%z/BOHRtoANGS
        enddo       
        write(unt,'(A43,A,A5,I12)') section, vartype,"   N=",j
        write(unt,'(5ES16.8)') aux(1:j)

        deallocate(aux)

        return

    end subroutine write_fchk_geom

!===========================================================================

    subroutine write_fchk_E(unt,system,Energ)

        !Write fchk file (including energy)

        use structure_types
        use constants

        integer,intent(in)::unt
        type(str_resmol),intent(inout)::system
#ifdef DOUBLE
        double precision,intent(in) :: Energ
#else
        real,intent(in) :: Energ
#endif

        !local
        integer::i, j
!         integer,dimension(:), allocatable :: aux_int
#ifdef DOUBLE
        double precision, dimension(:), allocatable :: aux 
#else
        real, dimension(:), allocatable :: aux 
#endif
        character :: vartype
        character(len=43) :: section
        character(len=90) :: job_sect
        
        !Title
         write(unt,'(A)') "File converted to fchk. Title:"//trim(adjustl(system%title))//&
                                                 " Name:"//trim(adjustl(system%name))
        !jog type and method
         job_sect=adjustl(system%job%type)//&
                  adjustl(system%job%method)//&
                  adjustl(system%job%basis)
         write(unt,'(A)') job_sect

        !Number of atoms
        section="Number of atoms"
        vartype="I"
        write(unt,'(A43,A,I17)') section, vartype, system%natoms

        !Atomic numbers
        section="Atomic numbers"
        vartype="I"
        j=system%natoms
        write(unt,'(A43,A,A5,I12)') section, vartype,"   N=",j
        write(unt,'(6I12)') system%atom(1:j)%AtNum

        !Nuclear charges
        section="Nuclear charges"
        vartype="R"
        j=system%natoms
        write(unt,'(A43,A,A5,I12)') section, vartype,"   N=",j
        write(unt,'(5ES16.8)') float(system%atom(1:j)%AtNum)

        !Cartesian coordinates
        allocate( aux(1:3*system%natoms) )
        section="Current cartesian coordinates"
        vartype="R"
        j=0
        do i=1,system%natoms
            j=j+1
            aux(j)=system%atom(i)%x/BOHRtoANGS
            j=j+1
            aux(j)=system%atom(i)%y/BOHRtoANGS
            j=j+1
            aux(j)=system%atom(i)%z/BOHRtoANGS
        enddo       
        write(unt,'(A43,A,A5,I12)') section, vartype,"   N=",j
        write(unt,'(5ES16.8)') aux(1:j)

        !SCF Energy
        section="SCF Energy"
        vartype="R"
        write(unt,'(A43,A,ES27.15)') section, vartype, Energ

        !Total Energy
        section="Total Energy"
        vartype="R"
        write(unt,'(A43,A,ES27.15)') section, vartype, Energ

        deallocate(aux)

        return

    end subroutine write_fchk_E


!===========================================================================
!   The following are deprecated (the general read_fchk SR should be used)
!===========================================================================

    subroutine read_dens(unt,P,N,dtype)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read density matrix (in a.o. basis) from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! P(real,dim(:);out): vector corresponding to density matrix in lower 
        !             triangular form
        ! N(int;out): number of elements in the vector [=Nbasis(Nbasis+1)/2]
        ! dtype(char;in): type of density (SCF, spin, alpha, beta, CI)
        !==============================================================


        !Arguments
#ifdef DOUBLE
        double precision, dimension(:), intent(out) :: P
#else
        real, dimension(:), intent(out) :: P
#endif
        integer, intent(out) :: N
        character(len=*), intent(in) :: dtype

        !Reading stuff
        character(len=50) :: header_of_section
        character(len=240) :: line=""

        !Auxiliar variables and Dummies
!         character(len=30) :: dummy_char
!         integer :: dummy_int

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
        select case ( adjustl(dtype) )
            case ( "SCF" )
             header_of_section="Total SCF Density"
            case ( "spinSCF" )
             header_of_section="Spin SCF Density"
            case ( "CI" )
             header_of_section="Total CI Density"
            case ( "spinCI" )
             header_of_section="Spin CI Density"
            case default
             call alert_msg( "fatal","Incorrect usage of read_dens subroutine. dtype cannot be: "//trim(adjustl(dtype)) )
        end select

        ! Search density
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching density: "//trim(adjustl(dtype)))
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') N
                    exit
                endif
        enddo

        ! Read density
         read(unt,*) P(1:N)
 

        rewind(unt)
        return

    end subroutine read_dens

    subroutine read_MO(unt,MO,N,dtype,isel)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read density matrix (in a.o. basis) from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! MO(real,dim(:,:);out): array which columns correspond to MO coefficients
        !                        in rows: MO(1:M,i) ar MO_i coefficients
        ! N(int;out): number of elements read
        ! dtype(char;out) : type of orbitals (alpha/beta)
        ! isel(char;in): target a given MO. If -1, retrieve all
        !==============================================================


        !Arguments
#ifdef DOUBLE
        double precision, dimension(:,:), intent(out) :: MO
#else
        real, dimension(:,:), intent(out) :: MO
#endif
        integer, intent(in)  :: isel
        integer, intent(out) :: N
        character(len=*), intent(in) :: dtype

        !Local auxiliar array (will store ALL coefficients in the SR)
#ifdef DOUBLE
        double precision, dimension(:,:),allocatable :: MO_L
#else
        real, dimension(:,:),allocatable :: MO_L
#endif
        integer :: M

        !Reading stuff
        character(len=50) :: header_of_section
        character(len=240) :: line=""

        !Auxiliar variables and Dummies
!         character(len=30) :: dummy_char
!         integer :: dummy_int

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
        select case ( adjustl(dtype) )
            case ( "alpha" )
             header_of_section="Alpha MO coefficients"
            case ( "beta" )
             header_of_section="Beta MO coefficients"
            case default
             call alert_msg( "fatal","Incorrect usage of read_MO subroutine. dtype cannot be: "//trim(adjustl(dtype)) )
        end select

        ! Search MO type
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching density: "//trim(adjustl(dtype)))
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') N
                    exit
                endif
        enddo

        !Calculate the number of states as the sqrt of N
        M = int ( sqrt(float(N)) )

        !Allocate MO_L
        allocate(MO_L(1:M,1:M))

        ! Read density
         read(unt,*) MO_L(1:M,1:M)

        if (isel == -1) then
            MO(1:M,1:M) = MO_L(1:M,1:M)
        else
            MO(1,1:M)   = MO_L(isel,1:M)
        endif
 

        rewind(unt)
        return

    end subroutine read_MO


    subroutine read_MOenergy(unt,E_MO,N,dtype,isel)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read density matrix (in a.o. basis) from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! E_MO(real,dim(:);out): vector with MO energies
        ! N(int;out): number of elements read
        ! dtype(char;out) : type of orbitals (alpha/beta)
        ! isel(char;in): target a given MO. If -1, retrieve all
        !==============================================================


        !Arguments
#ifdef DOUBLE
        double precision, dimension(:), intent(out) :: E_MO
#else
        real, dimension(:), intent(out) :: E_MO
#endif
        integer, intent(in)  :: isel
        integer, intent(out) :: N
        character(len=*), intent(in) :: dtype

        !Local auxiliar array (will store ALL coefficients in the SR)
#ifdef DOUBLE
        double precision, dimension(:),allocatable :: E_MO_L
#else
        real, dimension(:),allocatable :: E_MO_L
#endif
        integer :: M

        !Reading stuff
        character(len=50) :: header_of_section
        character(len=240) :: line=""

        !Auxiliar variables and Dummies
!         character(len=30) :: dummy_char
!         integer :: dummy_int

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
        select case ( adjustl(dtype) )
            case ( "alpha" )
             header_of_section="Alpha Orbital Energies"
            case ( "beta" )
             header_of_section="Beta Orbital Energies"
            case default
             call alert_msg( "fatal","Incorrect usage of read_MOenergy subroutine. dtype cannot be: "//trim(adjustl(dtype)) )
        end select

        ! Search MO type
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while searching density: "//trim(adjustl(dtype)))
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') N
                    exit
                endif
        enddo

        !Allocate E_MO_L
        allocate(E_MO_L(1:N))

        ! Read density
         read(unt,*) E_MO_L(1:N)

        if (isel == -1) then
            E_MO(1:N) = E_MO_L(1:N)
        else
            E_MO(1)   = E_MO_L(isel)
        endif
 

        rewind(unt)
        return

    end subroutine read_MOenergy


    subroutine read_NrEl(unt,N,dtype)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Get nr of electrons from fchk file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! Nel(int;out): number of electrons
        ! dtype(char;out) : type of orbitals (alpha/beta)
        !==============================================================


        !Arguments
        integer, intent(out) :: N
        character(len=*), intent(in) :: dtype

        !Reading stuff
        character(len=50) :: header_of_section
        character(len=240) :: line=""

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================


        !Set variables
        select case ( adjustl(dtype) )
            case ( "total" )
             header_of_section="Number of electrons"
            case ( "alpha" )
             header_of_section="Number of alpha electrons"
            case ( "beta" )
             header_of_section="Number of beta electrons"
            case default
             call alert_msg( "fatal","Incorrect usage of read_NrEl subroutine. dtype cannot be: "//trim(adjustl(dtype)) )
        end select

        ! Search MO type
        do
                read(unt,'(A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while getting NrEl: "//trim(adjustl(dtype)))
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    read(line,'(50X,I12)') N
                    exit
                endif
        enddo

        rewind(unt)
        return

    end subroutine read_NrEl



!     subroutine read_freq_fckh
! 
!         !==============================================================
!         ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
!         !==============================================================
!         !Description
!         ! Read normal mode information from Gaussian fchk file. It reads 
!         ! frequencies and normal mode description through the cartesian
!         ! L matrix 
!         !Arguments
!         ! unt (int;in): unit number of the log file
!         ! molec (str_system,inout): molecule 
!         !==============================================================
! 
!     end subroutine read_freq_fckh


end module gaussian_fchk_manage

