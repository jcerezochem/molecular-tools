program compute_overlap

    !=================================================================
    ! Description
    ! -----------
    ! Program to compute the overlap matrix from the 
    ! coefficients:
    !   C' = X C # orthogonalizaton of the basis (X can be S^1/2)
    !   C'^t C' = I = C^t X^tX C = C^t S C
    !   S = (C^t)^-1 C^-1
    ! If the basis set if not linearly independent, the number 
    ! of orbitals is not the same as the number of basis functions
    ! and C is rectangular (Nb x No). We need to use the generalized
    ! inverse:
    !   C^+ = (C^tC)^-1 C^t
    !   (C^t)^+ = (C^+)^t = C (C^tC)^-1
    ! And the overlap matrix is given by
    !   
    !   S = C (C^tC)^-1 (C^tC)^-1 C^t
    !
    !=================================================================


    subroutine overlap_matrix(Nb,No,MO,S)

        integer,intent(in)                 :: Nb, No
        real(8),dimension(:,:),intent(in)  :: MO
        real(8),dimension(:,:),intent(out) :: S

        ! S = (C^t)^-1 C^-1 = (C^-1)^t C^-1
        ! In the case of non-independent basis, Norb<Nbasis
        ! we use the generalized inverse
        ! S = C (C^tC)^-1 (C^tC)^-1 C^t
        S(1:No,1:No) = matrix_product(No,No,Nb,MO,MO,tA=.true.)
        S(1:No,1:No) = inverse_realgen(No,S)
        S(1:No,1:No) = matrix_product(No,No,No,S,S)
        S(1:Nb,1:No) = matrix_product(Nb,No,No,MO,S)
        S(1:Nb,1:Nb) = matrix_product(Nb,Nb,No,S,MO,tB=.true.)
        
        return

    end subroutine overlap_matrix

end program compute_overlap