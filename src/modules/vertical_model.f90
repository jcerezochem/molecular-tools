module vertical_model
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012!

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines with macros useful for
    !  frequency analysis on non-stationary points
    ! 
    !==============================================================

    !Common declarations:
    !===================
    use matrix
    use matrix_print
    use verbosity
    use constants
    use alerts
    use structure_types
    use symmetry

    implicit none

    contains

    subroutine check_symm_gsBder(molecule,gsBder)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Diagonalizes a Hessian after mass-weighing and translation to 
        ! the internal frame defined by satifying the Eckart-Saytvez conditions.
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom
        ! X,Y,Z   (inp) real/vectors  Coordinate vectors (ANGSTRONG)
        ! Mass    (inp) real/vector   Atomic masses (AMU)
        ! Freq    (inp) real/vector   Frequencies (cm-1)
        !
        !
        !==============================================================

        type(str_resmol),intent(inout)       :: molecule
        real(8),dimension(:,:),intent(in)    :: gsBder

        !Local
        integer,parameter :: NDIM = 600

        character(len=5) :: PG_current
        integer,dimension(1:NDIM) :: isym
        integer,dimension(1:4,1:NDIM,1:NDIM) :: Osym
        integer :: Nsym, Nat
        !Counters
        integer :: iop, i, j
        ! Aux scalar and matrix
        real(8) :: Theta, Theta2
        real(8),dimension(NDIM,NDIM) :: Aux

        Nat = molecule%natoms

        print'(/,X,A)', "---------------------------------------"
        print'(X,A  )', " Check effect of symmetry operations"
        print'(X,A  )', " on the correction term gs^t\beta"
        print'(X,A  )', "---------------------------------------"
        PG_current = molecule%PG
        molecule%PG="XX"
        call symm_atoms(molecule,isym,Osym,rotate=.false.,nsym_ops=nsym)
        ! Check the symmetry of the correction term
        ! Check all detected symmetry ops
        do iop=1,Nsym
            Aux(1:3*Nat,1:3*Nat) = dfloat(Osym(iop,1:3*Nat,1:3*Nat))
            Aux(1:3*Nat,1:3*Nat) = matrix_basisrot(3*Nat,3*Nat,Aux,gsBder,counter=.true.)
            Theta=0.d0
            do i=1,3*Nat 
            do j=1,3*Nat 
                if (Theta < abs(Aux(i,j)-gsBder(i,j))) then
                    Theta = abs(Aux(i,j)-gsBder(i,j))
                    Theta2=gsBder(i,j)
                endif
            enddo
            enddo
            print'(X,A,I0)', "Symmetry operation :   ", iop
            print'(X,A,F10.6)',   " Max abs difference : ", Theta
            print'(X,A,F10.6,/)', " Value before sym op: ", Theta2
        enddo
        print'(X,A,/)', "---------------------------------------"

        molecule%PG = PG_current

        return

    end subroutine check_symm_gsBder

    function projection_matrix(Nat,Ns,B) result(P)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Diagonalizes a Hessian after mass-weighing and translation to 
        ! the internal frame defined by satifying the Eckart-Saytvez conditions.
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Ns      (inp) int /scalar   Number of internal coordinates
        !
        !
        !==============================================================

        integer,intent(in)                 :: Nat, Ns
        real(8),dimension(:,:),intent(in)  :: B
        real(8),dimension(3*Nat,3*Nat)     :: P

        !Local
        integer,parameter :: NDIM = 600

        real(8),dimension(NDIM,NDIM) :: Aux

        Aux(1:Ns,1:Ns)     = matrix_product(Ns,Ns,3*Nat,B,B,tB=.true.)
        Aux(1:Ns,1:Ns)     = inverse_realgen(Ns,Aux)
        Aux(1:3*Nat,1:Ns)  = matrix_product(3*Nat,Ns,Ns,B,Aux,tA=.true.)
        P(1:3*Nat,1:3*Nat) = matrix_product(3*Nat,3*Nat,Ns,Aux,B)

        return

    end function projection_matrix

end module vertical_model