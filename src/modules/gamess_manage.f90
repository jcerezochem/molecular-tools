module gamess_manage

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to get molecular information
    !   from GAMESS output
    !    
    ! Notes  
    !  All subroutines rewind the file after using it
    !==============================================================

    !Common declarations:
    !===================
    use alerts
    use constants
    use structure_types
    implicit none


    contains

    subroutine read_gamess_geom(unt,molec)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        !
        !Arguments
        !==============================================================

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
                if ( adjustl(line) == "ATOM      ATOMIC                      COORDINATES (BOHR)" ) then
                    read(unt,'(A)') line !part of the header
                    exit
                endif
        enddo
        
        i=0
        do
            read(unt,'(A)',IOSTAT=IOstatus) line
            if ( IOstatus /= 0 ) call alert_msg("fatal","while reading coordinates from gamess out")
            if (len_trim(line) == 0 ) exit
            i=i+1
            read(line,*) molec%atom(i)%name, molec%atom(i)%q, molec%atom(i)%x, molec%atom(i)%y, molec%atom(i)%z
            molec%atom(i)%AtNum = int(molec%atom(i)%q)
        enddo
        molec%natoms=i
        molec%atom(1:i)%x = molec%atom(1:i)%x*BOHRtoANGS
        molec%atom(1:i)%y = molec%atom(1:i)%y*BOHRtoANGS
        molec%atom(1:i)%z = molec%atom(1:i)%z*BOHRtoANGS

        rewind(unt)
        return

    end subroutine read_gamess_geom


    subroutine read_gamess_hess(unt,N,hess,error_flag)


        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        !
        !Arguments
        !==============================================================

        integer,intent(in) :: unt
        integer,intent(in) :: N
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
        integer :: i, j, imax, imin, jmin, &
                   ib, nblocks
        
        
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
                if ( adjustl(line) == "CARTESIAN FORCE CONSTANT MATRIX" ) then
                    read(unt,*) line
                    exit
                endif
        enddo

        !Organized in blocks, displaying 6 columns each (2atoms)
        nblocks = N/6
        if (6*nblocks /= N) nblocks=nblocks+1
        do ib=1,nblocks
            !Pass headers
            read(unt,'(A)') line !blank
            read(unt,'(A)') line !index
            read(unt,'(A)') line !atom names
            read(unt,'(A)') line !axis labels
            !Parse hessian elements
            imin = (ib-1)*6 + 1
            imax = ib    *6
            imax=min(imax,N)
            read(unt,'(20X,6F9.6)') hess(imin:imax,imin:N)
        enddo

        !Complete symmetric blocks
        do i=1,N
            do j=1,i
                hess(i,j) = hess(j,i)
            enddo
        enddo

        rewind(unt)
        return

    end subroutine read_gamess_hess


!     subroutine read_gamess_nm(unt,Nv,N,Freq,mu,L,error_flag)
! 
! 
!         !==============================================================
!         ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
!         !==============================================================
!         !Description
!         ! Read normal modes from GAMESS output
!         ! Returs: 
!         !   Freq(cm-1)
!         !   mu (AMU)
!         !   Lcar' (normalized, dimesionelss): the same as G09 outputs
!         ! 
!         !Arguments
!         ! unt   (inp) scalar   unit for the file
!         ! N     (inp) scalar   degrees of freedom (3*Nat)
!         ! Nv    (inp) scalar   vibrational degrees of freedom
!         ! Freq  (out) vector   Frequencies (cm-1)
!         ! mu    (out) vector   Red. masses (AMU)
!         ! L     (out) matrix   normal modes arranged as L(NvxN)
!         ! error_flag (out) scalar  error_flag 
!         !
!         ! Notes
!         !  Using variable format as suggested in:
!         !  http://stackoverflow.com/questions/9881186/fortran-output-format-dependent-on-a-variable
!         !  (Answer by eriktous)
!         !==============================================================
! 
!         integer,intent(in) :: unt
!         integer,intent(in) :: Nv, N
! #ifdef DOUBLE
!         double precision, dimension(:,:), intent(out) :: L
!         double precision, dimension(:), intent(out) :: Freq, mu
! #else
!         real, dimension(:,:), intent(out) :: L
!         real, dimension(:), intent(out) :: Freq
! #endif
!         integer,intent(out) :: error_flag
! 
!         !Local stuff
!         !=============
!         character(len=240) :: line=""
!         double precision, dimension(1:N) :: F_local, mu_local
!         double precision, dimension(1:N,1:N) :: L_local
!         !I/O
!         integer :: IOstatus
!         character(len=100) :: fmt
!         !Counters
!         integer :: i, j, imax, imin, jmin, &
!                    ib, nblocks, icols, &
!                    irt
!         
!         
!         ! Search section
!         error_flag = 0
!         do
!                 read(unt,'(A)',IOSTAT=IOstatus) line
!                 ! Two possible scenarios while reading:
!                 ! 1) End of file
!                 if ( IOstatus < 0 ) then
!                     call alert_msg("note","End of file but normal modes not read")
!                     error_flag=1
!                     rewind(unt)
!                     return
!                 endif
!                 ! 2) Found what looked for!      
!                 if ( index(line,"FREQUENCIES IN CM**-1, IR INTENSITIES IN DEBYE**2/AMU-ANGSTROM**2") /= 0 ) then
!                     read(unt,*) line
!                     exit
!                 endif
!         enddo
! 
!         !Organized in blocks, displaying 5 columns each 
!         nblocks = N/5
!         if (5*nblocks /= N) nblocks=nblocks+1
!         do ib=1,nblocks
!             imin = (ib-1)*5 + 1 
!             imax = ib    *5 
!             imax=min(imax,N)
!             icols = 1 + (imax-imin)
!             write(fmt,'(a,i0,a)') '(20X,',icols,'F12.8)'
!             !Pass headers
!             read(unt,'(A)') line !blank
!             read(unt,'(A)') line !index
!             !Parse frequencies
!             read(unt,'(18X,5F12.2)') F_local(imin:imax)
!             read(unt,'(A)') line !symmetry
!             read(unt,'(18X,5F12.5)') mu_local(imin:imax)
!             read(unt,'(A)') line !IR int 
!             read(unt,'(A)') line !blank 
!             !Parse L matrix elements
!             read(unt,fmt) L_local(imin:imax,1:N)
!             read(unt,'(A)') line !blank
!             read(unt,'(A)') line !Trans X
!             read(unt,'(A)') line !Trans Y
!             read(unt,'(A)') line !Trans Z
!             read(unt,'(A)') line !Trans tot
!             read(unt,'(A)') line !blank
!             read(unt,'(A)') line !Rot X
!             read(unt,'(A)') line !Rot Y
!             read(unt,'(A)') line !Rot Z
!             read(unt,'(A)') line !Rot tot
!         enddo
! 
! do i=1,N
! print'(20X,5F12.8)', L_local(1,i), L_local(N-3:N,i)
! enddo
! 
!         !Discard rotationa and translation
!         ! NOTE: GAMESS prints Lcar(in AMU^-1/2)
!         !       so we multiply by mu^1/2 to get the 
!         !       normalized vectors (as in G09)
!         irt = N-Nv
!         do i=1,Nv
!             Freq(i) = F_local(i+irt)
!             mu(i)   = mu_local(i+irt) 
!         do j=1,N
!             L(i,j) = L_local(i+irt,j)*sqrt(mu(i))
!         enddo
!         enddo
! 
!         rewind(unt)
!         return
! 
!     end subroutine read_gamess_nm

end module gamess_manage

