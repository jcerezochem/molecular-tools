module molpro_manage

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information 
    !  from molpro log/out files:
    !    read_log_xyz(unt,molec)
    !  Note: all subroutines rewind the file after using it
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


    subroutine get_log_xyz(unt,molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read xyz geometry from the last record (optimized) on molpro
        ! log file (note this is the long output file, not the main one)
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
        logical :: final_section, found_geometry

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
        header_of_section="Current geometry"
!         end_of_section="***************"
        n_useles_lines=1

        ! Take the last Standard orientation in file
        final_section=.false.
        found_geometry=.false.
        do
            do 
                read(unt,'(X,A)',IOSTAT=IOstatus) line
                ! Three possible scenarios while reading:
                ! 1) End of file
                if ( IOstatus < 0 ) then
                    final_section=.true.
                    exit
                endif
                ! 2) Final section reached
!                 if ( INDEX(line,"GINC") /= 0 ) then
!                     final_section=.true.
!                     exit
!                 endif
                ! 3) Found what looked for!      
                if ( INDEX(line,trim(adjustl(header_of_section))) /= 0 ) then
                    found_geometry=.true.
                    exit
                endif
            enddo
            if (final_section) exit

            ! Overpass lines until reaching the target table
            do j=1,n_useles_lines
                read(unt,'(X,A)',IOSTAT=IOstatus) line 
            enddo

            read(unt,*) molec%natoms
            read(unt,*) molec%title

            ! Read Standard orientation to a predefined end of table (no need to know the number of atoms)
            do i=1,molec%natoms
                read(unt,*,IOSTAT=IOstatus) molec%atom(i)%name, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z
                if ( IOstatus < 0 ) call alert_msg("fatal","Unexpected end of file while scanning for xyz geometry")
            enddo

        enddo

        if ( .not.found_geometry ) call alert_msg("fatal","Unexpected end of file while searching xyz geometry")

        rewind(unt)
        return

    end subroutine get_log_xyz


    subroutine read_freq_molpro(unt,Nvib,Nat,freq,L,err_label)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read frequencies. The format is like HP modes in G09. But 
        ! Lcart matrix vectors are not normalized
        !Arguments
        !==============================================================

        integer,intent(in) :: unt
        integer, intent(in) :: Nvib, Nat
#ifdef DOUBLE
        double precision, dimension(:), intent(out) :: freq
        double precision, dimension(:,:), intent(out) :: L
#else
        real, dimension(:), intent(out) :: freq
        real, dimension(:,:), intent(out) :: L
#endif
        integer,intent(out) :: err_label

        !Lookup auxiliar variables
        character(len=240) :: line, subline, cnull
        character(len=2) :: modes, modes_prev
        character(len=4) splitter
        character(len=10000) :: cfreq
        character(len=10000),dimension(1:1000) :: cL
! #ifdef DOUBLE
!         double precision, dimension(1:500,1:500) :: L
! #else
!         real, dimension(1:500,1:500) :: L
! #endif
        integer :: nlines, ngroups, N, Nremain


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

            if ( INDEX(line,"Normal Modes") /= 0 ) then
                err_label = 0
                exit
            endif
        enddo
        
        if (err_label == -1) then
            call alert_msg("fatal","No normal modes read in molpro logfile")
        endif

        !Modes are written five-by-five
        ngroups = Nvib/5
        if (ngroups*5 < Nvib ) ngroups=ngroups+1
        do i=1,ngroups
            read(unt,'(X,A)',IOSTAT=IOstatus) cnull ! blank line
            read(unt,'(X,A)',IOSTAT=IOstatus) cnull ! Vibratio symmetry
            read(unt,'(A27,A)',IOSTAT=IOstatus) cnull, line ! Wavenumbers
            Nremain = Nvib - (i-1)*5
            N = min(Nremain,5)
            read(line,*) freq((i-1)*5+1:(i-1)*5+N)
            read(unt,'(X,A)',IOSTAT=IOstatus) cnull ! Intensities
            read(unt,'(X,A)',IOSTAT=IOstatus) cnull ! Intensities (relative)
            !Reading L matrix (Nvib x 3Nat)
            do j=1,3*Nat
                read(unt,'(A27,A)',IOSTAT=IOstatus) cnull, line
                read(line,*)  L((i-1)*5+1:(i-1)*5+N,j)
            enddo

        enddo
! print*, Freq(1:Nvib)
! do i=1,Nvib
!     print'(100F10.2)', L(i,1:3*Nat)
! enddo
! stop
!         !The subroutine returns a vector (in the same order as stored in fchk)
!         !(so we forget about the two indexed...)
!         k=0
!         do i=1,Nvib
!             do j=1,3*Nat
!                 k=k+1
!                 Lvector(k) = L(j,i)
!             enddo
!         enddo
! 
!         if ( modes == "LP" ) err_label = 1

        rewind(unt)
        return

    end subroutine read_freq_molpro


    subroutine read_freq_only_molpro(unt,Nvib,freq,err_label)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Read frequencies. The format is like HP modes in G09. But 
        ! Lcart matrix vectors are not normalized
        !Arguments
        !==============================================================

        integer,intent(in) :: unt
        integer, intent(in) :: Nvib
#ifdef DOUBLE
        double precision, dimension(:), intent(out) :: freq
#else
        real, dimension(:), intent(out) :: freq
#endif
        integer,intent(out) :: err_label

        !Lookup auxiliar variables
        character(len=240) :: line, subline, cnull
        character(len=4) splitter
! #ifdef DOUBLE
!         double precision, dimension(1:500,1:500) :: L
! #else
!         real, dimension(1:500,1:500) :: L
! #endif
        integer :: nlines, ngroups, N, Nremain


        !Counters and dummies (TODO: dummies module)
        integer :: i,j, IOstatus, k

        err_label = -1
        do
            read(unt,'(X,A)',IOSTAT=IOstatus) line
            if ( IOstatus /= 0) exit

            if ( INDEX(line,"  Vibration") /= 0 ) then
                err_label = 0
                exit
            endif
        enddo
        
        if (err_label == -1) then
            call alert_msg("fatal","No normal modes read in molpro logfile")
        endif

        !Read column titles
        read(unt,*) cnull

        !Read normal frequencies
        do i=1,Nvib
            write(cnull,*) i
            read(unt,*,IOSTAT=IOstatus) j, freq(i)
            if ( IOstatus /= 0) call alert_msg("fatal","while reading frequency: "//trim(adjustl(cnull))//" from Molpro log")
        enddo

        rewind(unt)
        return

    end subroutine read_freq_only_molpro



end module molpro_manage

