module gaussian_manage_notypes

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
    !==============================================================

    !Common declarations:
    !===================
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


    subroutine read_freq_only(unt,N,freq)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read normal mode information from Gaussian log file. It reads 
        ! frequencies and normal mode description through the cartesian
        !Arguments
        !==============================================================

        integer,intent(in) :: unt
        integer, intent(in) :: N
#ifdef DOUBLE
        double precision, dimension(:), intent(out) :: freq
#else
        real, dimension(:), intent(out) :: freq
#endif

        !Lookup auxiliar variables
        character(len=240) :: line, subline
        character(len=4) splitter
        character(len=10000) :: cfreq

        !Counters and dummies (TODO: dummies module)
        integer :: i,j, IOstatus

        cfreq = ""
        do
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus /= 0) exit

            if ( INDEX(line,"Frequencies") /= 0 ) then
                if ( INDEX(line,'---') /= 0 ) then
                    !High precision modes
                    splitter="---"
                else
                    !Low precision modes
                    splitter="--"
                endif

                call split_line(line,trim(adjustl(splitter)),line,subline)
                cfreq = trim(adjustl(cfreq))//" "//trim(adjustl(subline))
            endif
        enddo

        read(cfreq,*) freq(1:N)

        rewind(unt)
        return

    end subroutine read_freq_only


    subroutine read_freq_NT(unt,Nvib,Nat,freq,Lvector,err_label)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read normal mode information from Gaussian log file. It reads 
        ! frequencies and normal mode description through the cartesian
        !Arguments
        !==============================================================

        integer,intent(in) :: unt
        integer, intent(in) :: Nvib, Nat
#ifdef DOUBLE
        double precision, dimension(:), intent(out) :: freq, Lvector
#else
        real, dimension(:), intent(out) :: freq, Lvector
#endif
        integer,intent(out) :: err_label

        !Lookup auxiliar variables
        character(len=240) :: line, subline, cnull
        character(len=2) :: modes, modes_prev
        character(len=4) splitter
        character(len=10000) :: cfreq
        character(len=10000),dimension(1:1000) :: cL
#ifdef DOUBLE
        double precision, dimension(1:500,1:500) :: L
#else
        real, dimension(1:500,1:500) :: L
#endif
        integer :: nlines, N


        !Counters and dummies (TODO: dummies module)
        integer :: i,j, IOstatus, k

        cfreq = ""
        cL = ""
        modes=""
        modes_prev=""
        err_label = -1
        do
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus /= 0) exit

            if ( INDEX(line,"Frequencies") /= 0 ) then
                err_label = 0

                if ( INDEX(line,'---') /= 0 ) then
                    !High precision modes
                    splitter="---"
                    nlines=3*Nat
                    modes="HP"
                else
                    !Low precision modes
                    splitter="--"
                    nlines=Nat
                    modes="LP"
                endif

                ! When HP modes are available:
                ! LP modes come after HP, but they are not interesting
                if ( modes == "LP" .and. modes_prev == "HP" ) then
                    modes = "HP"
                    nlines=3*Nat
                    exit
                endif
                modes_prev = modes

                call split_line(line,trim(adjustl(splitter)),line,subline)
                cfreq = trim(adjustl(cfreq))//" "//trim(adjustl(subline))

                !Look for the Lcart matrix
                do 
                    read(unt,'(X,A)',IOSTAT=IOstatus) line
                    if ( IOstatus /= 0) stop
                    if ( INDEX(line,"Atom") /= 0 ) exit
                enddo
                do i=1,nlines
                    if (modes == "HP") then
                        read(unt,'(A23,A)') cnull, line
                    elseif (modes == "LP") then
                        read(unt,'(A13,A)') cnull, line
                    endif
                    cL(i) = trim(adjustl(cL(i)))//" "//trim(adjustl(line))
                enddo

            endif
        enddo

        read(cfreq,*) freq(1:Nvib)
        N = Nvib * 3*Nat
        
        do i=1,nlines
            if ( modes == "LP" ) then
                ! In LP modes coordinates for a given atom are grouped therefore
                ! it's read using the fast index (row) running i-2 --> i
                j = i*3
                read(cL(i),*) L(j-2:j,1:Nvib)
            elseif ( modes == "HP" ) then
                ! HP modes are read line by line through the 3Nat coordinates
                read(cL(i),*) L(i,1:Nvib)
            endif
        enddo

        !The subroutine returns a vector (in the same order as stored in fchk)
        !(so we forget about the two indexed...)
        k=0
        do i=1,Nvib
            do j=1,3*Nat
                k=k+1
                Lvector(k) = L(j,i)
            enddo
        enddo

        if ( modes == "LP" ) err_label = 1

        rewind(unt)
        return

    end subroutine read_freq_NT

    subroutine get_ori_geom(unt,geom,orientation,err_label)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the gaussian standard geometry from the last record (optimized)
        ! from the log file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! geom  
        !==============================================================

        !Arguments
#ifdef DOUBLE
        double precision,dimension(:),intent(inout) :: geom
#else
        real,dimension(:),intent(inout) :: geom
#endif
        character(len=*),intent(in) :: orientation
        integer,intent(out) :: err_label


        !Reading stuff
        character(len=50) :: header_of_section, &
                             end_of_section
        integer :: n_useles_lines
        character(len=240) :: line=""

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
        header_of_section=trim(adjustl(orientation))
        end_of_section="----------"
        n_useles_lines=4

        ! Take the last Standard orientation in file (of the last job --NEW)
        err_label = 1
        do
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    if (err_label == 1) call alert_msg("note","Exit without finding "//trim(adjustl(orientation)))
                    rewind(unt)
                    return
                endif 
                ! 2) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) exit
            enddo
            err_label = 0
            ! Overpass lines until reaching the target table
            do j=1,n_useles_lines
                read(unt,'(X,A)',IOSTAT=IOstatus) line 
            enddo

            ! Read Standard orientation to a predefined end of table (no need to know the number of atoms)
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while scanning "//trim(adjustl(orientation)))
                ! 2) End of table
                if ( INDEX(line,trim(adjustl(end_of_section))) /= 0 ) exit    
                ! 3) Table entry
                read(line,*) i, dummy_int, dummy_int, geom(3*i-2:3*i)
            enddo

        enddo

        rewind(unt)
        return

    end subroutine get_ori_geom

    subroutine get_ori_geom2(unt,geom,orientation,isel,err_label)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the gaussian standard geometry from the last record (optimized)
        ! from the log file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! geom  
        !==============================================================

        !Arguments
#ifdef DOUBLE
        double precision,dimension(:),intent(inout) :: geom
#else
        real,dimension(:),intent(inout) :: geom
#endif
        character(len=*),intent(in) :: orientation
        integer,intent(in) :: isel
        integer,intent(out) :: err_label


        !Reading stuff
        character(len=50) :: header_of_section, &
                             end_of_section
        integer :: n_useles_lines
        character(len=240) :: line=""

        !Auxiliar variables and Dummies
        character(len=30) :: dummy_char
        integer :: dummy_int

        !=============
        !Counters
        integer :: i,j, ifound
        !=============

        !================
        !I/O stuff
        !units
        integer,intent(in) :: unt
        !status
        integer :: IOstatus
        !===================


        !Set variables
        header_of_section=trim(adjustl(orientation))
        end_of_section="----------"
        n_useles_lines=4
        ifound = 0

        ! Take the last Standard orientation in file (of the last job --NEW)
        err_label = 1
        do
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Two possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    if (err_label == 1) call alert_msg("note","Exit without finding "//trim(adjustl(orientation)))
                    rewind(unt)
                    return
                endif 
                ! 2) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) exit
            enddo
            ifound = ifound + 1
            err_label = 0
            ! Overpass lines until reaching the target table
            do j=1,n_useles_lines
                read(unt,'(X,A)',IOSTAT=IOstatus) line 
            enddo

            ! Read Standard orientation to a predefined end of table (no need to know the number of atoms)
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while scanning "//trim(adjustl(orientation)))
                ! 2) End of table
                if ( INDEX(line,trim(adjustl(end_of_section))) /= 0 ) exit    
                ! 3) Table entry
                read(line,*) i, dummy_int, dummy_int, geom(3*i-2:3*i)
            enddo
            if (ifound == isel) then
                print*, "Geom block retrieved: ", ifound
                rewind(unt)
                return
            endif

        enddo

        rewind(unt)
        return

    end subroutine get_ori_geom2


    subroutine get_Natoms(unt,Nat)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the gaussian standard geometry from the last record (optimized)
        ! from the log file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! geom  
        !==============================================================

        !Arguments
        integer,intent(out) :: Nat

        integer :: err_label

        !Reading stuff
        character(len=240) :: line=""

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

        err_label = 1
        do 
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            ! 1) End of file
            if ( IOstatus < 0 ) then
                if (err_label == 1) call alert_msg("note","Exit without retrieving Nat")
                rewind(unt)
                return
            endif 
            ! 2) Found what looked for!      
            if ( INDEX(line,"NAtoms=") /= 0 ) exit
        enddo
        err_label = 0

        read(line,*) dummy_char, Nat, dummy_char

        rewind(unt)
        return

    end subroutine get_Natoms

    subroutine get_PG(unt,PG)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read the gaussian standard geometry from the last record (optimized)
        ! from the log file
        !Arguments
        ! unt (int;in): unit number of the log file
        ! geom  
        !==============================================================

        !Arguments
        character(len=*),intent(out) :: PG
        
        integer :: err_label

        !Reading stuff
        character(len=240) :: line=""

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

        err_label = 1
        do 
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    if (err_label == 1) call alert_msg("note","Exit without retrieving PG")
                    rewind(unt)
                    return
                endif 
                ! 2) Found what looked for!      
                if ( INDEX(line,"Full point group") /= 0 ) exit
            enddo
            err_label = 0

            read(line,*) dummy_char,dummy_char,dummy_char,PG, dummy_char

        enddo

        rewind(unt) 
        return

    end subroutine get_PG



end module gaussian_manage_notypes

