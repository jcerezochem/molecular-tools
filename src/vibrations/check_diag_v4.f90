program check_diag

    !Compilation 
    ! $FC ../modules/alerts.f90 ../modules/MatrixMod.f90 ../modules/constants_mod.f90 ../modules/structure_types_v2.f90 ../modules/gaussian_fchk_manage_v2.f90 check_diag.f90 -llapack -o check_diag.exe -cpp -DDOUBLE

    use alerts
    use structure_types
    use molecular_structure
    use line_preprocess
    use constants
    use gro_manage
    use gaussian_fchk_manage
    use molcas_UnSym_manage
    use MatrixMod

    integer,parameter :: NDIM=5000, NDIM2=600

    type(str_resmol) ::  molec

    real(8),dimension(NDIM) :: H_lt, freq, amass
    real(8),dimension(NDIM2,NDIM2) :: H, L
    integer :: Nat, Nvib

    !Stuf for fchk reading
    integer,dimension(:),allocatable :: IA
    real(8),dimension(:),allocatable :: A
    integer :: N, error, N_lt
    character(len=1) :: dtype

    logical :: debug = .false.

    !Commandline input stuff
     character(len=100) :: input
     character(len=6) :: ext

    call getarg(1, input) 
    open(10,file=input,status="old")

    call split_line(input,".",input,ext)

    if (iargc() > 1) call getarg(2, input)
    if (adjustl(input) == "dbg") debug=.true.

    if (adjustl(ext) == "UnSym") then

        call read_molcas_geom(10,molec)
        call atname2element(molec)
        call assign_masses(molec)
        Nat=molec%natoms

        call read_molcas_hess(10,N,H,error)

        k=0
        do i=1,Nat
        do j=1,3
            k=k+1
            amass(k) = molec%atom(i)%mass
        enddo
        enddo

        k=0
        do i=1,3*Nat
            do j=1,i
                k=k+1
                H(i,j) = H(i,j)/dsqrt(amass(i))/dsqrt(amass(j))
                H(j,i) = H(i,j)
                H_lt(k) = H(i,j)
            enddo
        enddo

    elseif (adjustl(ext) == "fchk") then

        ! Natoms
        call read_fchk(10,"Number of atoms",dtype,N,A,IA,error)
        if (error == 0) then
            Nat = IA(1)
            deallocate(IA)
        endif

        ! Atom masses
        call read_fchk(10,"Real atomic weights",dtype,N,A,IA,error)
        if (error == 0) then
            k=0
            do i=1,N
            do j=1,3
                k=k+1
                amass(k) = A(i)
            enddo
            enddo
            deallocate(A)
        endif


        call read_fchk(10,"Cartesian Force Constants",dtype,N,A,IA,error)
        if (error == 0) then
            k=0
            do i=1,3*Nat
                do j=1,i
                    k=k+1
                    H(i,j) = A(k)/dsqrt(amass(i))/dsqrt(amass(j))
                    H(j,i) = H(i,j)
                    H_lt(k) = H(i,j)
                enddo
            enddo
            N_lt = N
            deallocate(A)
        endif

        call read_fchk(10,"Vib-E2",dtype,N,A,IA,error)
        if (error == 0) then
            print*, ""
            print*, "Freq (Gaussian)="
            do i=1,3*Nat-6
                print'(100(F10.3,2X))', A(i)
            enddo
            deallocate(A)
        endif

    else

        call alert_msg("fatal","Qué carajo me das para leer")

    endif

    if (debug) then
    print*, ""
    print*, "H="
    do i=1,3*Nat
        print'(100(F10.3,2X))', H(i,1:3*Nat)
    enddo 
    endif

    print*, ""
    print*, "==================="
    print*, " FULL - lapack"
    print*, "==================="
    call diagonalize_full(H(1:3*Nat,1:3*Nat),3*Nat,L(1:3*Nat,1:3*Nat),freq(1:3*Nat),"lapack")

    if (debug) then
    print*, ""
    print*, "L="
    do i=1,3*Nat
        print'(100(F10.3,2X))', L(i,1:3*Nat)
    enddo 
    endif

    call sort_vec(Freq,3*Nat)
    print*, ""
    print*, "Freq="
    do i=1,3*Nat
        print'(100(F10.3,2X))', dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
    enddo 


!     print*, ""
!     print*, "==================="
!     print*, " FULL - emtccm"
!     print*, "==================="
!     call diagonalize_full(H(1:3*Nat,1:3*Nat),3*Nat,L(1:3*Nat,1:3*Nat),freq(1:3*Nat),"emtccm")
! 
!     print*, ""
!     print*, "L="
!     do i=1,3*Nat
!         print'(100(F10.3,2X))', L(i,1:3*Nat)
!     enddo 
! 
!     call sort_vec(Freq,3*Nat)
!     print*, ""
!     print*, "Freq="
!     do i=1,3*Nat
!         print'(100(F10.3,2X))', dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
!     enddo 
! 
!     print*, ""
!     print*, "==================="
!     print*, " LT - emtccm"
!     print*, "==================="
!     call diagonalize_lt(H_lt(1:N_lt),3*Nat,L(1:3*Nat,1:3*Nat),freq(1:3*Nat),"emtccm")
! 
!     print*, ""
!     print*, "L="
!     do i=1,3*Nat
!         print'(100(F10.3,2X))', L(i,1:3*Nat)
!     enddo 
! 
!     call sort_vec(Freq,3*Nat)
!     print*, ""
!     print*, "Freq="
!     do i=1,3*Nat
!         print'(100(F10.3,2X))', dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
!     enddo 


!     print*, ""
!     print*, "==================="
!     print*, " FULL - eigrs"
!     print*, "==================="
!     call diagonalize_full(H(1:3*Nat,1:3*Nat),3*Nat,L(1:3*Nat,1:3*Nat),freq(1:3*Nat),"eigrs")
! 
!     print*, ""
!     print*, "L="
!     do i=1,3*Nat
!         print'(100(F10.3,2X))', L(i,1:3*Nat)
!     enddo 
! 
!     call sort_vec(Freq,3*Nat)
!     print*, ""
!     print*, "Freq="
!     do i=1,3*Nat
!         print'(100(F10.3,2X))', dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
!     enddo 
! 
! 
!     print*, ""
!     print*, "==================="
!     print*, " LT - eigrs"
!     print*, "==================="
!     call diagonalize_lt(H_lt(1:N_lt),3*Nat,L(1:3*Nat,1:3*Nat),freq(1:3*Nat),"eigrs")
! 
!     print*, ""
!     print*, "L="
!     do i=1,3*Nat
!         print'(100(F10.3,2X))', L(i,1:3*Nat)
!     enddo 
! 
!     call sort_vec(Freq,3*Nat)
!     print*, ""
!     print*, "Freq="
!     do i=1,3*Nat
!         print'(100(F10.3,2X))', dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
!     enddo 

    stop

end program check_diag