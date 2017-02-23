module gaussian_fchk_analisys

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to perform wavefunction analysis
    !    
    ! Notes
    !  Since variables molec, residue... are allocatable, they should
    !  be passed to the sr even if they are output to get the allocation.
    !==============================================================

    !Common declarations:
    !===================
    use alerts
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

    subroutine DA_analysis(Pdif,N,D,A)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! 
        !Arguments
        ! 
        !==============================================================

        use MatrixMod

        !Input
        real,dimension(:),intent(in) :: Pdif
        integer, intent(in) :: N

        !Output
        real,dimension(:),intent(out) :: D,A

        !Local
        integer :: NBasis
        real :: suma
        !Arrays Nbasis x Nbasis
        real,dimension(:,:),allocatable :: Pdif_F, A_F, D_F, U, Ut
        real,dimension(:),allocatable :: delta, dp, ap
 
        !Counters
        integer :: i,j,k
        

        !Get NBasis from N
        NBasis=INT((-1.+sqrt(1.+8.*N))/2)

        !Then allocate arrays
        allocate( &!Pdif_F(NBasis,NBasis), &
                  A_F(NBasis,NBasis),    &
                  D_F(NBasis,NBasis),    &
                  U(NBasis,NBasis),      &
                  Ut(NBasis,NBasis),     &
                  delta(NBasis),         &
                  ap(NBasis),            &
                  dp(NBasis)             )

        ! Reconstruct the full simmetric Pdif
!         k=0
!         do i=1,NBasis
!             do j=1,i
!                 k=k+1
!                 Pdif_F(i,j) = Pdif(k) 
!                 Pdif_F(j,i) = Pdif_F(i,j)
!             enddo
!         enddo

        !Diagonalize matrix
!         print*, "Diagonalizing density matrix"
! print*, "Pdif"
! print'(1000F8.3)', Pdif(1:N) 
! print*, "Pdif_F"
! do i=1,NBasis
!     print'(1000F8.3)', Pdif_F(i,1:NBasis)
! enddo
! print*, ""
!         call JACOBI_SIM(NBasis,PREC,Pdif_F(1:NBasis,1:NBasis),U(1:NBasis,1:NBasis),delta(1:NBasis))
        call diagonalize_lt(Pdif(1:N),NBasis,U(1:NBasis,1:NBasis),delta(1:NBasis),"emtccm") ! eigrs emtccm
! print*, "Pdfif"
! print*, Pdif(1:N)
! print*, "delta"
! print*, delta(1:NBasis)
        !Get the transpose (=inverse) of U to rotate the density
        do i=1,NBasis
            do j=1,NBasis
                Ut(i,j) = U(j,i)
            enddo
        enddo

!     call EIGRS(Pdif(N),              &
!                NBasis,               &
!                1,                    &
!                delta(1:NBasis),      &
!                U(1:NBasis,1:NBasis), &
!                NBasis,               &
!                dummy_array,          &
!                IOstatus)  

        !Detachement
        dp = -min(delta,0.)
        print*, "Detached electrons"
        suma=0.
        do i=1,NBasis
            suma=suma+dp(i)
            if (dp(i) > 0.1) print*, "AO", i, "Eigenvalue", dp(i)
        enddo
        print'(X,A6,3X,F15.8,/)', "Total: ",suma
        D_F=0
        do i=1, NBasis
            D_F(i,i) = dp(i)
        enddo
        D_F(1:NBasis,1:NBasis)=matmul(U(1:NBasis,1:NBasis),D_F(1:NBasis,1:NBasis))
        D_F(1:NBasis,1:NBasis)=matmul(D_F(1:NBasis,1:NBasis),Ut(1:NBasis,1:NBasis))

        !Attachement
        ap =  max(delta,0.)
        print*, "Attached electrons"
        suma=0.
        do i=1,NBasis
            suma=suma+ap(i)
            if (ap(i) > 0.1) print*, "AO", i, "Eigenvalue", ap(i)
        enddo
        print'(X,A6,3X,F15.8,/)', "Total: ",suma
        A_F=0
        do i=1, NBasis
            A_F(i,i) = ap(i)
        enddo
        A_F(1:NBasis,1:NBasis)=matmul(U(1:NBasis,1:NBasis),A_F(1:NBasis,1:NBasis))
        A_F(1:NBasis,1:NBasis)=matmul(A_F(1:NBasis,1:NBasis),Ut(1:NBasis,1:NBasis))


! print*, "Pdif_F"
! do i=1,NBasis
!     print'(1000F8.3)', Pdif_F(i,1:NBasis)
! enddo
! print*, ""
! print*, "A-D"
! do i=1,NBasis
!     print'(1000F8.3)', A_F(i,1:NBasis)-D_F(i,1:NBasis)
! enddo
! print*, ""

        !Lower triangular form
        k=0
        do i=1,NBasis
            do j=1,i
                k=k+1
                D(k) = D_F(i,j)
                A(k) = A_F(i,j)
            enddo
        enddo    

!         deallocate(Pdif_F,A_F,D_F,U,Ut,delta,dp,ap)

        return

    end subroutine DA_analysis


    subroutine replace_DA(unt,unt2,D,A,N,tune)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
        !==============================================================
        !Description
        ! Generates a fchk where the density specified with "tune" is replaced by the DA density
        ! 
        !Arguments
        ! 
        !==============================================================

        !Input
        integer,intent(in) :: unt, unt2, N
        character(len=*),intent(in) :: tune
        real, dimension(:), intent(in) :: D,A

        !Local
        real,dimension(1:N) :: Aux
        character(len=17) :: SCF_header, CI_header
        integer :: IOstatus
        character(len=80) :: line

        !Where to replace
        if ( adjustl(tune) == "Total" ) then
            SCF_header="Total SCF Density"
            CI_header= "Total CI Density "
        elseif ( adjustl(tune) == "Spin" ) then
            SCF_header="Spin SCF Density"
            CI_header= "Spin CI Density"
        endif


        ! Search density
        do
            read(unt2,'(A)',IOSTAT=IOstatus) line
            ! Two possible scenarios while reading:
            ! 1) End of file
            if ( IOstatus < 0 ) then
                print'(X,A,/)', "DONE"
                exit
            else
                write(unt,'(A)') line
                if ( INDEX(line,trim(adjustl(SCF_header))) /= 0 ) then
                    read(unt2,*) Aux(1:N)
                    write(unt,'(5(X,E15.8))') D(1:N)
                elseif ( INDEX(line,trim(adjustl(CI_header))) /= 0 ) then
                    read(unt2,*) Aux(1:N)
                    write(unt,'(5(X,E15.8))') A(1:N)
                endif
            endif
        enddo

        rewind(unt)
        rewind(unt2)
        return

    end subroutine replace_DA



end module gaussian_fchk_analisys

