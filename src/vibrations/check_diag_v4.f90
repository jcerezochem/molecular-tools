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
    use gaussian_manage_lowlevel
    use molcas_UnSym_manage
    use psi4_manage
    use gamess_manage
    use MatrixMod

    integer,parameter :: NDIM=1000, NDIM2=300

    type(str_resmol) ::  molec

    real(8),dimension(NDIM) :: freq, amass
    real(8),dimension(320000) :: H_lt
    real(8),dimension(NDIM2,NDIM2) :: H, L, D
    real(8),dimension(3,3) :: MI, X
    real(8),dimension(3) :: R
    real(8) :: Px, Py, Pz, det, pes
    integer :: Nat, Nvib

    !Stuf for fchk reading
    integer,dimension(:),allocatable :: IA
    real(8),dimension(:),allocatable :: A
    integer :: N, error, N_lt
    character(len=1) :: dtype
    character(len=700000) :: string
    character(len=50) :: geom_char

    logical :: debug = .false.

    !Commandline input stuff
     character(len=100) :: input, input_add
     character(len=6) :: ext, ext_add
     character :: cnull

    call getarg(1, input) 
    open(10,file=input,status="old")

    call split_line_back(input,".",input,ext)

    if (iargc() > 1) then
        call getarg(2, input)
        if (adjustl(input) == "dbg") then
            debug=.true.
        else
            call getarg(2, input_add)
            open(11,file=input_add,status="old")
            call split_line_back(input_add,".",input_add,ext_add)
            if (iargc() > 2) call getarg(3, input)
            if (adjustl(input) == "dbg") debug=.true.
        endif
    endif


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


    elseif (adjustl(ext) == "log") then

        ! Natoms
        call read_gausslog_natoms(10,Nat,error)
        molec%natoms = Nat

        !Geom
        call summary_parser(10,4,string,error)
        !Throw "charge mult" away
        call split_line(string,'\',geom_char,string)
        do i=1,Nat
            call split_line(string,'\',geom_char,string)
            read(geom_char,*) molec%atom(i)%name, &
                              molec%atom(i)%x,    &
                              molec%atom(i)%y,    &
                              molec%atom(i)%z
        enddo
        !Get mass...
        call atname2element(molec)
        call assign_masses(molec)
        k=0
        do i=1,Nat
        do j=1,3
            k=k+1
            amass(k) = molec%atom(i)%mass
        enddo
        enddo

        !Place at COM 
        call get_com(molec)
        do i=1,molec%natoms
            molec%atom(i)%x = molec%atom(i)%x - molec%comX
            molec%atom(i)%y = molec%atom(i)%y - molec%comY
            molec%atom(i)%z = molec%atom(i)%z - molec%comZ
        enddo

        !Get moment of intertia
        call inertia(molec,MI)

        call diagonalize_full(MI(1:3,1:3),3,X(1:3,1:3),freq(1:3),"lapack")

        X=transpose(X)

        !Following G09 white paper
        ! Note that there is a typo:
        !  * Rotational coordinates should have m^1/2 factor multiplied, not divided
        ! Furthermore, there additional issues are  
        !  * Confusing  matrix indices. Note that X should be transposed first to use the
        !    orther they use
        D(1:3*Nat,1:3*Nat) = 0.d0
        !Traslation
        do i=1,3*Nat,3
            D(i  ,1) = dsqrt(amass(i)) 
            D(i+1,2) = dsqrt(amass(i+1)) 
            D(i+2,3) = dsqrt(amass(i+2)) 
        enddo
        !Rotation
        do i=1,3*Nat,3
            j=(i-1)/3+1
            R=(/molec%atom(j)%x,molec%atom(j)%y,molec%atom(j)%z/)
            Px=dot_product(R(1:3),X(1,1:3))
            Py=dot_product(R(1:3),X(2,1:3))
            Pz=dot_product(R(1:3),X(3,1:3))

            D(i  ,4) = (Py*X(3,1) - Pz*X(2,1))*dsqrt(amass(i))
            D(i+1,4) = (Py*X(3,2) - Pz*X(2,2))*dsqrt(amass(i))
            D(i+2,4) = (Py*X(3,3) - Pz*X(2,3))*dsqrt(amass(i))

            D(i  ,5) = (Pz*X(1,1) - Px*X(3,1))*dsqrt(amass(i))
            D(i+1,5) = (Pz*X(1,2) - Px*X(3,2))*dsqrt(amass(i))
            D(i+2,5) = (Pz*X(1,3) - Px*X(3,3))*dsqrt(amass(i))

            D(i  ,6) = (Px*X(2,1) - Py*X(1,1))*dsqrt(amass(i))
            D(i+1,6) = (Px*X(2,2) - Py*X(1,2))*dsqrt(amass(i))
            D(i+2,6) = (Px*X(2,3) - Py*X(1,3))*dsqrt(amass(i))
        enddo

! print*, ""
!         print'(F9.5)', dot_product(D(1:3*Nat,1),D(1:3*Nat,4))
!         print'(F9.5)', dot_product(D(1:3*Nat,1),D(1:3*Nat,5))
!         print'(F9.5)', dot_product(D(1:3*Nat,1),D(1:3*Nat,6))
! print*, ""
!         print'(F9.5)', dot_product(D(1:3*Nat,2),D(1:3*Nat,4))
!         print'(F9.5)', dot_product(D(1:3*Nat,2),D(1:3*Nat,5))
!         print'(F9.5)', dot_product(D(1:3*Nat,2),D(1:3*Nat,6))
! print*, ""
!         print'(F9.5)', dot_product(D(1:3*Nat,3),D(1:3*Nat,4))
!         print'(F9.5)', dot_product(D(1:3*Nat,3),D(1:3*Nat,5))
!         print'(F9.5)', dot_product(D(1:3*Nat,3),D(1:3*Nat,6))
! print*, ""
!         print'(F15.5)', dot_product(D(1:3*Nat,4),D(1:3*Nat,5))
!         print'(F15.5)', dot_product(D(1:3*Nat,4),D(1:3*Nat,6))
!         print'(F15.5)', dot_product(D(1:3*Nat,5),D(1:3*Nat,6))
! print*, ""
! print*, "Rotations:"
!         print'(F15.5)', dot_product(D(1:3*Nat,4),D(1:3*Nat,4))
!         print'(F15.5)', dot_product(D(1:3*Nat,5),D(1:3*Nat,5))
!         print'(F15.5)', dot_product(D(1:3*Nat,6),D(1:3*Nat,6))

        !Normalize
        D(1:3*Nat,1) = D(1:3*Nat,1)/sqrt(dot_product(D(1:3*Nat,1),D(1:3*Nat,1)))
        D(1:3*Nat,2) = D(1:3*Nat,2)/sqrt(dot_product(D(1:3*Nat,2),D(1:3*Nat,2)))
        D(1:3*Nat,3) = D(1:3*Nat,3)/sqrt(dot_product(D(1:3*Nat,3),D(1:3*Nat,3)))
        D(1:3*Nat,4) = D(1:3*Nat,4)/sqrt(dot_product(D(1:3*Nat,4),D(1:3*Nat,4)))
        D(1:3*Nat,5) = D(1:3*Nat,5)/sqrt(dot_product(D(1:3*Nat,5),D(1:3*Nat,5)))
        D(1:3*Nat,6) = D(1:3*Nat,6)/sqrt(dot_product(D(1:3*Nat,6),D(1:3*Nat,6)))


        !Get remaining vectors by G-S orthogonalization
        ! We add a complete canonical basis in R^N.
        ! The additional 6 (5) vectors are set to zero
        ! by the G-S procedure
        nd = 0
        do i=7,3*Nat+6
            ii = i - nd
            D(i-6,ii) = 1.d0
            Freq(1:3*Nat) = D(1:3*Nat,ii)
            do j=1,ii-1
                Freq(1:3*Nat) = Freq(1:3*Nat) - dot_product(D(1:3*Nat,ii),D(1:3*Nat,j))*D(1:3*Nat,j)
            enddo
            D(1:3*Nat,ii) = Freq(1:3*Nat)
            pes = dot_product(D(1:3*Nat,ii),D(1:3*Nat,ii))
            if (pes < 1.d-5) then
                print*, "Discarded",i,  pes
                nd = nd + 1
            else 
                D(1:3*Nat,ii) = D(1:3*Nat,ii)/sqrt(pes)
            endif
        enddo
        print*, ""
        print*, "Discarded vectors:", nd
        print*, ""

        !Hessian
        call summary_parser(10,6,string,error)
        if (error /= 0) call alert_msg("fatal","Error on summary_parser")
        allocate(A(1:3*Nat*(3*Nat+1)/2))
        read(string,*) A(1:3*Nat*(3*Nat+1)/2)

        k=0
        do i=1,3*Nat
            do j=1,i
                k=k+1
                H(i,j) = A(k)/dsqrt(amass(i))/dsqrt(amass(j))
                H(j,i) = H(i,j)
                H_lt(k) = H(i,j)
            enddo
        enddo
        deallocate(A)

        !Rotate Hessian
        H(1:3*Nat,1:3*Nat-6) = matmul(H(1:3*Nat,1:3*Nat),D(1:3*Nat,7:3*Nat))
        D = transpose(D)
        H(1:3*Nat-6,1:3*Nat-6) = matmul(D(7:3*Nat,1:3*Nat),H(1:3*Nat,1:3*Nat-6))


    print*, ""
    print*, "==================="
    print*, " FULL - lapack"
    print*, "==================="
    call diagonalize_full(H(1:3*Nat-6,1:3*Nat-6),3*Nat-6,L(1:3*Nat-6,1:3*Nat-6),freq(1:3*Nat-6),"lapack")

    if (debug) then
    print*, ""
    print*, "L="
    do i=1,3*Nat
        print'(100(F10.3,2X))', L(i,1:3*Nat)
    enddo 
    endif

    call sort_vec(Freq,3*Nat-6)
    print*, ""
    print*, "Freq="
    do i=1,3*Nat-6
        print'(100(F10.3,2X))', dsign(dsqrt(dabs(Freq(i))*HARTtoJ/BOHRtoM**2/UMAtoKG)/2.d0/pi/clight/1.d2,Freq(i))
    enddo 


    elseif (adjustl(ext) == "ghess") then
        !you should add a g96 file with the structure
        if (adjustl(ext_add) /= "g96") &
         call alert_msg("fatal","With ghess you should provide a g96 structure file as second input")

        call read_g96(11,molec)

        call atname2element(molec)
        call assign_masses(molec)
        Nat=molec%natoms

        call read_gro_hess(10,N,H,error)

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

    elseif (adjustl(ext) == "psi4") then

        call read_psi_geom(10,molec)
        !Get mass...
        call atname2element(molec)
        call assign_masses(molec)
        Nat=molec%natoms
        N=3*Nat

        call read_psi_hess(10,N,H,error)

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

    elseif (adjustl(ext) == "gamess") then

        call read_gamess_geom(10,molec)
        !Get mass...
        call atname2element(molec)
        call assign_masses(molec)
        Nat=molec%natoms
        N=3*Nat

        call read_gamess_hess(10,N,H,error)

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

    else

        call alert_msg("fatal","Qué carajo me das para leer")

    endif

    if (debug) then
    print*, ""
    print*, "H="
    do i=1,3*Nat
        print'(100(F10.3,2X))', H(i,1:3*Nat)
    enddo 
    print*, "mass="
    do i=1,Nat
        print'(100(F10.3,2X))', amass(i)
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