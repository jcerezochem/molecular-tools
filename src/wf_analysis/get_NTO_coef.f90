program get_NTO_coef

    !==============================================================
    ! This code uses MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    ! Description
    !  Utility to retrieve NTO coeficients
    ! Compilation intructions:
    !  gfortran ../modules/gaussian_fchk_manage_v2.f90 get_NTO_coef.f90 -o get_NTO_coef.exe -cpp 
    !
    !==============================================================

    use gaussian_fchk_manage


    implicit none

    

    !====================== 
    !Read fchk auxiliars
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: IA
    character(len=1) :: dtype
    integer :: error, N
    !====================== 

    character(len=150) :: fileinp 
    integer :: I_FCHK=10

    real(8),dimension(1000) :: C_NTO
    integer :: i

    call getarg(1, fileinp)
    open(I_FCHK,file=fileinp,status="old")

    !Alpha
    N=0
    call read_fchk(I_FCHK,"Alpha Orbital Energies",dtype,N,A,IA,error)
    C_NTO(1:N) = A(1:N)
    deallocate(A)

    if (N/=0) then
        print*, "Alpha"
        call sort_vec_max(C_NTO,N)
        do i=1,N,2
            if (C_NTO(i) < 0.1d0) exit
            print*,  C_NTO(i)
        enddo

    endif

    !Beta
    N=0
    call read_fchk(I_FCHK,"Beta Orbital Energies",dtype,N,A,IA,error)
    C_NTO(1:N) = A(1:N)
    deallocate(A)

    if (N/=0) then
        print*, "Beta"
        call sort_vec_max(C_NTO,N)
        do i=1,N,2
            if (C_NTO(i) < 0.1d0) exit
            print*,  C_NTO(i)
        enddo

    endif
    

    stop

    contains

    subroutine sort_vec_max(V,N)

        implicit none

#ifdef DOUBLE
        double precision,dimension(:),intent(inout) :: V
        double precision :: aux
#else
        real,dimension(:),intent(inout) :: V
        real :: aux
#endif  
!         integer,dimension(:),intent(inout) :: IORD
        integer,intent(in) :: N
        integer :: i,j, iaux

!         !Intialize IORD
!         do i=1,N
!             IORD(i) = i
!         enddo

        do i=1,N-1
            do j=i+1,N
                if (V(j)>V(i)) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                    !Track the index permutations in IORD
!                     iaux = IORD(i)
!                     IORD(i) = IORD(j)
!                     IORD(j) = iaux
                endif
            enddo
        enddo

        return

    end subroutine sort_vec_max      

end program get_NTO_coef

    

