module vibrational_analysis
!23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012!

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to perform the vibrational
    !  analysis, by properly Diagonalizing the Hessian
    ! 
    !==============================================================

    !Common declarations:
    !===================
    use matrix
    use matrix_print
    use verbosity
    use constants

    implicit none

    contains

    subroutine vibrations_Cart(Nat,X,Y,Z,Mass,Hlt,Nvib,L,Freq,error_flag,Din,Dout)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Diagonalizes a Hessian after mass-weighing and translation to 
        ! the internal frame defined by satifying the Eckart-Saytvez conditions.
        !
        !Arguments
        ! Nat     (int) int /scalar   Number of atoms
        ! X,Y,Z   (inp) real/vectors  Coordinate vectors (ANGSTRONG)
        ! Mass    (inp) real/vector   Atomic masses (AMU)
        ! Hlt     (inp) real/vector   Lower triangular part of the Hessian in Cartesian coordinates (AU)
        ! Nvib    (out) int /scalar   Number of vibrational degrees of freedom
        ! L       (out) real/matrix   Normal modes (vibrations only) in Cartesian coordinates (Dimless)
        ! Freq    (out) real/vector   Frequencies (cm-1)
        ! Optional kwargs:
        ! error_flag (out) flag  0 : success
        !                        1 : Wrong number of Trans+Rot (<5)
        !                        2 : Wrong number of independent vectors
        ! Dout    (out) real/matrix   Rotation matrix to Eckart frame (get as output)
        ! Din     (out) real/matrix   Rotation matrix to Eckart frame (provide as input)
        !
        !          
        !Notes
        ! Uses built-in matrix manipulation functions in fortran (dot_product,matmul,transpose)
        ! We assume the compiler would recognise them (true for gfortran, ifort)  
        !
        ! General note: better not to "reuse" output arrays as work arrays if the output
        ! size is not as large as the work to do. E.g., Freq can simply be (1:Nvib), so 
        ! do not use is as a work array for (1:3*Nat)     
        !
        !==============================================================

        !Approximate zero
        real(kind=8),parameter :: ZERO=1.d-10

        integer,intent(in)                      :: Nat
        real(kind=8),dimension(:),intent(in)    :: X,Y,Z
        real(kind=8),dimension(:),intent(in)    :: Mass
        real(kind=8),dimension(:),intent(in)    :: Hlt
        integer,intent(out)                     :: Nvib
        real(kind=8),dimension(:,:),intent(out) :: L
        real(kind=8),dimension(:),intent(out)   :: Freq
        integer,intent(out),optional :: error_flag
        ! Rotation to Ekart frame
        real(kind=8),dimension(:,:),intent(out),optional  :: Dout
        real(kind=8),dimension(:,:),intent(in),optional   :: Din

        !Local
        ! Counters
        integer :: i, j, k, ii, jj
        integer :: Nrt
        integer :: error_local
        ! Auxiliar scalars and arrays
        real(kind=8)                :: pes, MassTot
        real(kind=8),dimension(3)   :: R, RCOM
        real(kind=8),dimension(3,3) :: MI, Xrot
        real(kind=8),dimension(1:3*Nat,1:3*Nat+6) :: D
        real(kind=8),dimension(3*Nat,3*Nat)       :: H
        real(kind=8),dimension(Nat)               :: Mass_local


        !Working in AU: transform mass to AU
        Mass_local(1:Nat) = Mass(1:Nat) * AMUtoAU

        ! Compute rotation matrix to the Eckart frame
        if (.not.present(Din)) then
            
            !Get COM 
            RCOM(1:3) = 0.d0
            MassTot   = 0.d0
            do i=1,Nat
                RCOM(1) = RCOM(1) + X(i)*Mass_local(i)
                RCOM(2) = RCOM(2) + Y(i)*Mass_local(i)
                RCOM(3) = RCOM(3) + Z(i)*Mass_local(i)
                MassTot = MassTot + Mass_local(i)
            enddo
            RCOM(1:3) = RCOM(1:3)/MassTot
            
            !Get moment of intertia
            MI=0.d0
            do i=1,Nat
                R=(/X(i)-RCOM(1),Y(i)-RCOM(2),Z(i)-RCOM(3)/)
                !diag
                MI(1,1)=MI(1,1)+Mass_local(i)*(R(2)**2+R(3)**2)
                MI(2,2)=MI(2,2)+Mass_local(i)*(R(1)**2+R(3)**2)
                MI(3,3)=MI(3,3)+Mass_local(i)*(R(1)**2+R(2)**2)
                !off-diag
                MI(2,1)=MI(2,1)-Mass_local(i)*(R(2)*R(1))
                MI(3,1)=MI(3,1)-Mass_local(i)*(R(3)*R(1))
                MI(3,2)=MI(3,2)-Mass_local(i)*(R(3)*R(2))
            enddo
            do i=1,3
              do j=1,i-1
                  MI(j,i) = MI(i,j)
              enddo
            enddo
            
            !Diagonalize to get the rotation to the principal axes
            call diagonalize_full(MI(1:3,1:3),3,Xrot(1:3,1:3),Freq(1:3),"lapack")
            ! MI = Xrot^t * Freq * Xrot
            !Note we need to transpose to follow G09 white paper
            Xrot=transpose(Xrot)
            
            !Get the orthogonal transformation to the internal frame
            ! we follow G09 white paper
            ! Note that there is a typo:
            !  * Rotational coordinates should have m^1/2 factor multiplied, not divided
            ! Furthermore, there additional issues are  
            !  * Confusing  matrix indices. Note that X should be transposed first to use the
            !    order they use
            D(1:3*Nat,1:3*Nat+6) = 0.d0
            !Traslation
            do i=1,3*Nat,3
                j=(i-1)/3+1
                !D(1)
                D(i  ,1) = dsqrt(Mass_local(j)) 
                !D(2)
                D(i+1,2) = dsqrt(Mass_local(j)) 
                !D(3)
                D(i+2,3) = dsqrt(Mass_local(j)) 
            enddo
            !Normalize
            D(1:3*Nat,1) = D(1:3*Nat,1)/sqrt(dot_product(D(1:3*Nat,1),D(1:3*Nat,1)))
            D(1:3*Nat,2) = D(1:3*Nat,2)/sqrt(dot_product(D(1:3*Nat,2),D(1:3*Nat,2)))
            D(1:3*Nat,3) = D(1:3*Nat,3)/sqrt(dot_product(D(1:3*Nat,3),D(1:3*Nat,3)))
            
            !Rotation
            do i=1,3*Nat,3
                j=(i-1)/3+1
                !Get Equil. coordinates in the principal axis frame 
                R=(/X(j)-RCOM(1),Y(j)-RCOM(2),Z(j)-RCOM(3)/)
                R(1:3) = matmul(Xrot(1:3,1:3),R(1:3))
                !D(4)
                D(i  ,4) = (R(2)*Xrot(3,1) - R(3)*Xrot(2,1))*dsqrt(Mass_local(j)) 
                D(i+1,4) = (R(2)*Xrot(3,2) - R(3)*Xrot(2,2))*dsqrt(Mass_local(j)) 
                D(i+2,4) = (R(2)*Xrot(3,3) - R(3)*Xrot(2,3))*dsqrt(Mass_local(j)) 
                !D(5)
                D(i  ,5) = (R(3)*Xrot(1,1) - R(1)*Xrot(3,1))*dsqrt(Mass_local(j)) 
                D(i+1,5) = (R(3)*Xrot(1,2) - R(1)*Xrot(3,2))*dsqrt(Mass_local(j)) 
                D(i+2,5) = (R(3)*Xrot(1,3) - R(1)*Xrot(3,3))*dsqrt(Mass_local(j)) 
                !D(5)
                D(i  ,6) = (R(1)*Xrot(2,1) - R(2)*Xrot(1,1))*dsqrt(Mass_local(j)) 
                D(i+1,6) = (R(1)*Xrot(2,2) - R(2)*Xrot(1,2))*dsqrt(Mass_local(j)) 
                D(i+2,6) = (R(1)*Xrot(2,3) - R(2)*Xrot(1,3))*dsqrt(Mass_local(j)) 
            enddo
            !Normalize (and determine if there is one equal to zero: linear molecules)
            ii = 3
            do i=4,6 
                ii = ii + 1
                pes=dot_product(D(1:3*Nat,ii),D(1:3*Nat,ii))
                if (abs(pes) < ZERO) then
                    print*, "NOTE: linear molecule detected"
                    !Shift comlumns
                    do j=ii,5
                        D(1:3*Nat,j) = D(1:3*Nat,j+1)
                    enddo
                    ii = ii - 1
                else
                    D(1:3*Nat,ii) = D(1:3*Nat,ii)/sqrt(pes)
                endif
            enddo
            Nrt = ii
            if (Nrt < 5) then
                error_local = 1
                if (present(error_flag)) error_flag = error_local
                return
            endif
            
            !Get remaining vectors (vibration) by G-S orthogonalization
            ! Adding a complete canonical basis in R^N. The additional 6 (5)
            ! vectors are set to zero by the G-S procedure
            do i=Nrt+1,3*Nat+Nrt
                ii = ii + 1
                !Initial guess: canonical vector (0,0,...0,1,0,...,0,0)
                D(i-Nrt,ii) = 1.d0
                !Grand-Schmdit orthogonalization step. (we use H(:,1) as work vector)
                H(1:3*Nat,1) = D(1:3*Nat,ii)
                do j=1,ii-1
                    H(1:3*Nat,1) = H(1:3*Nat,1) - dot_product(D(1:3*Nat,ii),D(1:3*Nat,j))*D(1:3*Nat,j)
                enddo
                D(1:3*Nat,ii) = H(1:3*Nat,1)
                !Check if the new vector was linearly independent, if so, normalize it
                pes = dot_product(D(1:3*Nat,ii),D(1:3*Nat,ii))
                if (abs(pes) < ZERO) then
                    !If not linear independent, shift columns
                    ii = ii - 1
                else 
                    D(1:3*Nat,ii) = D(1:3*Nat,ii)/sqrt(pes)
                endif
            enddo
            !Compute the number of vibrational degrees of freedom
            Nvib = ii-Nrt
            !Check if the G-S discarded the right number of vectors
            if ( 3*Nat+Nrt - ii /= Nrt ) then
                error_local = 2
                if (present(error_flag)) error_flag = error_local
                return
            endif

        else
            Nrt = 6
            D(1:3*Nat,1:3*Nat+6) = Din(1:3*Nat,1:3*Nat+6)
        endif
        ! At this point, also get it if requested
        if (present(Dout)) then
            Dout(1:3*Nat,1:3*Nat+6) = D(1:3*Nat,1:3*Nat+6)
        endif

        !Massweight the Hessian
        k=0
        do i=1,3*Nat
        do j=1,i
            k=k+1
            ii = (i-1)/3+1
            jj = (j-1)/3+1
            H(i,j) = Hlt(k)/sqrt(Mass_local(ii)*Mass_local(jj)) 
            H(j,i) = H(i,j)
        enddo 
        enddo

        !Rotate Hessian to the internal frame (3Nat-6/5 coordinates)
        H(1:3*Nat,1:Nvib) = matmul(H(1:3*Nat,1:3*Nat),D(1:3*Nat,Nrt+1:3*Nat))
        D = transpose(D)
        H(1:Nvib,1:Nvib)  = matmul(D(Nrt+1:3*Nat,1:3*Nat),H(1:3*Nat,1:Nvib))
        !Restore D
        D = transpose(D)

        !Diagonalize
        call diagonalize_full(H(1:Nvib,1:Nvib),Nvib,L(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
        !Transform L from internal frame (Nvib x Nvib) into Cartesian (3Nat x Nvib) using Mass and D(3Nat x Nvib)
        ! Lcart = m^1/2 D L
        L(1:3*Nat,1:Nvib) = matmul(D(1:3*Nat,Nrt+1:3*Nat),L(1:Nvib,1:Nvib))

        !Check FC
        if (verbose>1) &
            call print_vector(6,Freq*1.d6,Nvib,"FORCE CONSTANTS x 10^6 (A.U.)")

        !Transform to FC to Freq
        do i=1,Nvib
              Freq(i) = sign(dsqrt(abs(Freq(i))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
                             Freq(i))
        enddo
        if (verbose>0) &
            call print_vector(6,Freq,Nvib,"Frequencies (cm-1)")

        if (verbose>2) &
            call MAT0(6,L,3*Nat,Nvib,"Lmwc matrix")

        return

    end subroutine vibrations_Cart


    !=====================================================
    ! Conversion tools between different L definitions
    !=====================================================

    subroutine LcartNrm_to_Lmwc(Nat,Nvib,Mass,LcartNrm,Lmwc)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom
        ! Mass    (inp) real/vector   Atomic masses (AMU)
        ! LcartNrm(inp) real/matrix   Normalized Normal modes in Cartesian (adim) (3Nat x Nvib)
        ! Lmwc    (out) real/matrix   Normal modes in mxc          (mass_AU^-1/2) (3Nat x Nvib)
        !          
        !Notes
        ! We can use the same matrix as input and as output as there is
        ! no conflict
        !
        !==============================================================

        integer,intent(in)                      :: Nat, Nvib
        real(kind=8),dimension(:),intent(in)    :: Mass
        real(kind=8),dimension(:,:),intent(in)  :: LcartNrm
        real(kind=8),dimension(:,:),intent(out) :: Lmwc

        !Local
        integer :: i, j, k, kk
        real(8) :: p1
        real(kind=8),dimension(Nat) :: Mass_local

        !Working in AU: transform mass to AU
        Mass_local(1:Nat) = Mass(1:Nat) * AMUtoAU

        do k=1,Nvib
            kk=0
            p1=0.d0
            do i=1,Nat
            do j=1,3
                kk=kk+1
                Lmwc(kk,k)=LcartNrm(kk,k)*dsqrt(Mass_local(i))
                p1=p1+Lmwc(kk,k)**2
            enddo
            enddo

            ! Redefine T=T*MQ^(-0.5) i.e T=lmwc
            do I=1,3*Nat
                Lmwc(I,k)=Lmwc(I,k)/dsqrt(p1)
            enddo
        enddo

        return

    end subroutine LcartNrm_to_Lmwc


    subroutine Lmwc_to_Lcart(Nat,Nvib,Mass,Lmwc,Lcart,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Tranform the Lmwc into the Lcart.
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom
        ! Mass    (inp) real/vector   Atomic masses (AMU)
        ! Lmwc    (inp) real/matrix   Normal modes in mxc (adim)           (3Nat x Nvib)
        ! Lcart   (inp) real/matrix   Normal modes in Cartesian (mass_AU^-1/2) (3Nat x Nvib)
        ! error_flag (out) flag  0 : success
        !
        !Notes
        ! We can use the same matrix as input and as output as there is
        ! no conflict
        !
        !==============================================================

        integer,intent(in)                      :: Nat, Nvib
        real(kind=8),dimension(:),intent(in)    :: Mass
        real(kind=8),dimension(:,:),intent(in)  :: Lmwc
        real(kind=8),dimension(:,:),intent(out) :: Lcart
        integer,intent(out),optional            :: error_flag

        !Local
        integer :: i, ii
        integer :: error_local
        real(kind=8),dimension(Nat) :: Mass_local

        !Working in AU: transform mass to AU
        Mass_local(1:Nat) = Mass(1:Nat) * AMUtoAU

        error_local = 0
        do i=1,3*Nat
            ii = (i-1)/3+1
            Lcart(i,1:Nvib) = Lmwc(i,1:Nvib)/dsqrt(Mass_local(ii))
        enddo
        if (present(error_flag)) error_flag=error_local

        return

    end subroutine Lmwc_to_Lcart

    subroutine Lcart_to_Lmwc(Nat,Nvib,Mass,Lcart,Lmwc,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Tranform the Lcart into the Lmwc.
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom
        ! Mass    (inp) real/vector   Atomic masses (AMU)
        ! Lcart   (inp) real/matrix   Normal modes in Cartesian (mass_AU^-1/2) (3Nat x Nvib)
        ! Lmwc    (inp) real/matrix   Normal modes in mxc (adim)           (3Nat x Nvib)
        ! error_flag (out) flag  0 : success
        !
        !Notes
        ! We can use the same matrix as input and as output as there is
        ! no conflict
        !
        !==============================================================

        integer,intent(in)                      :: Nat, Nvib
        real(kind=8),dimension(:),intent(in)    :: Mass
        real(kind=8),dimension(:,:),intent(in)  :: Lcart
        real(kind=8),dimension(:,:),intent(out) :: Lmwc
        integer,intent(out),optional            :: error_flag

        !Local
        integer :: i, ii
        integer :: error_local
        real(kind=8),dimension(Nat) :: Mass_local

        !Working in AU: transform mass to AU
        Mass_local(1:Nat) = Mass(1:Nat) * AMUtoAU

        error_local = 0
        do i=1,3*Nat
            ii = (i-1)/3+1
            Lmwc(i,1:Nvib) = Lcart(i,1:Nvib)*dsqrt(Mass_local(ii))
        enddo
        if (present(error_flag)) error_flag=error_local

        return

    end subroutine Lcart_to_Lmwc


    subroutine Lcart_to_LcartNrm(Nat,Nvib,Lcart,LcartNrm,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Tranform the Lcart into the normalized Lcart as in G09. We could output the reduced mas from here
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom
        ! Lcart   (inp) real/matrix   Normal modes in mxc              (AU^-1/2) (3Nat x Nvib)
        ! LcartNrm(out) real/matrix   Normalized Normal modes in Cartesian (adim) (3Nat x Nvib)
        ! error_flag (out) flag  0 : success
        !                        1 : 
        !          
        !Notes
        ! We can use the same matrix as input and as output as there is
        ! no conflict
        !
        !==============================================================

        integer,intent(in)                      :: Nat, Nvib
        real(kind=8),dimension(:,:),intent(in)  :: Lcart
        real(kind=8),dimension(:,:),intent(out) :: LcartNrm
        integer,intent(out),optional            :: error_flag

        !Local
        integer :: i, j
        real(kind=8) :: Factor
        integer :: error_local

        error_local = 0
        do i=1,Nvib
            Factor=0.d0
            do j=1,3*Nat
                Factor = Factor + Lcart(j,i)**2
            enddo
            !Reduced_mass(i) = 1.d0/Factor
            Factor = 1.d0/dsqrt(Factor)
            LcartNrm(1:3*Nat,i) = Lcart(1:3*Nat,i)*Factor
        enddo
        if (present(error_flag)) error_flag=error_local

        return

    end subroutine Lcart_to_LcartNrm


    subroutine Ls_to_Lcart(Nat,Nvib,Mass,B,G,Ls,Lcart,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Tranform the Ls into the Lcart:
        !    Lcart = m^-1 B^t Ginv Ls
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom (Ns should be = Nvib!!)
        ! Mass    (inp) real/vector   Atomic masses (AMU)
        ! B       (inp) real/matrix   B Wilsn matrix (AU)
        ! G       (inp) real/matrix   G metric matrix (AU)
        ! Ls      (inp) real/matrix   Normal modes in internal coords (Nvib x Nvib)
        ! Lcart   (inp) real/matrix   Normal modes in Cartesian (mass_AU^-1/2) (3Nat x Nvib)
        ! error_flag (out) flag  0 : success
        !
        !Notes
        ! We can use the same matrix as input and as output as there is
        ! no conflict
        ! Only valid for Ns=Nvib
        !
        !==============================================================

        integer,intent(in)                      :: Nat, Nvib
        real(kind=8),dimension(:),intent(in)    :: Mass
        real(kind=8),dimension(:,:),intent(in)  :: B
        real(kind=8),dimension(:,:),intent(in)  :: G
        real(kind=8),dimension(:,:),intent(in)  :: Ls
        real(kind=8),dimension(:,:),intent(out) :: Lcart
        integer,intent(out),optional            :: error_flag

        !Local
        integer :: i, ii
        integer :: error_local
        real(kind=8),dimension(Nat) :: Mass_local
        real(kind=8),dimension(Nvib,Nvib) :: Ginv
        
        !Working in AU: transform mass to AU
        Mass_local(1:Nat) = Mass(1:Nat) * AMUtoAU

        ! Compute G^-1 
        Ginv(1:Nvib,1:Nvib) = inverse_realgen(Nvib,G)

        ! Compute B^t G^-1
        Lcart(1:3*Nat,1:Nvib)= matrix_product(3*Nat,Nvib,Nvib,B,Ginv,tA=.true.)
        ! Compute [B^t G^-1] L
        Lcart(1:3*Nat,1:Nvib) = matrix_product(3*Nat,Nvib,Nvib,Lcart,Ls)
        ! Compute m^-1 [B^t G^-1 L] 
        do i=1,3*Nat
            ii = (i-1)/3+1
            Lcart(i,1:Nvib) = Lcart(i,1:Nvib)/Mass_local(ii)
        enddo

        return 

    end subroutine Ls_to_Lcart

    subroutine Lcart_to_Ls(Nat,Nvib,B,Lcart,Ls,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Tranform the Lcart into the Ls:
        !    Ls = B Lcart 
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom (Ns should be = Nvib!!)
        ! B       (inp) real/matrix   B Wilsn matrix (AU)
        ! Lcart   (inp) real/matrix   Normal modes in internal coords (Nvib x Nvib)
        ! Ls      (inp) real/matrix   Normal modes in Cartesian (mass_AU^-1/2) (3Nat x Nvib)
        ! error_flag (out) flag  0 : success
        !
        !Notes
        ! We can use the same matrix as input and as output as there is
        ! no conflict
        ! Only valid for Ns=Nvib
        !
        !==============================================================

        integer,intent(in)                      :: Nat, Nvib
        real(kind=8),dimension(:,:),intent(in)  :: B
        real(kind=8),dimension(:,:),intent(in)  :: Lcart
        real(kind=8),dimension(:,:),intent(out) :: Ls
        integer,intent(out),optional            :: error_flag

        !Local
        integer :: i, ii
        integer :: error_local

        ! Compute B Lcart
        Ls(1:Nvib,1:Nvib)= matrix_product(Nvib,Nvib,3*Nat,B,Lcart)

        return 

    end subroutine Lcart_to_Ls

    subroutine mu_from_Lcart(Nat,Nvib,Lcart,mu,error_flag)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS
        !==============================================================
        !Description
        ! Tranform the Lcart into the normalized Lcart as in G09. We could output the reduced mas from here
        !
        !Arguments
        ! Nat     (inp) int /scalar   Number of atoms
        ! Nvib    (inp) int /scalar   Number of vibrational degrees of freedom
        ! Lcart   (inp) real/matrix   Normal modes in mxc              (mass_AU^-1/2) (3Nat x Nvib)
        ! mu      (out) real/vector   reduced mass array (AMU) (Nvib)
        ! error_flag (out) flag  0 : success
        !                        1 : 
        !          
        !Notes
        ! We can use the same matrix as input and as output as there is
        ! no conflict
        !
        !==============================================================

        integer,intent(in)                      :: Nat, Nvib
        real(kind=8),dimension(:,:),intent(in)  :: Lcart
        real(kind=8),dimension(:),intent(out)   :: mu
        integer,intent(out),optional            :: error_flag

        !Local
        integer :: i, j
        integer :: error_local

        error_local = 0
        do i=1,Nvib
            mu(i)=0.d0
            do j=1,3*Nat
                mu(i) = mu(i) + Lcart(j,i)**2
            enddo
            !Reduced_mass(i) = 1.d0/Factor
            mu(i) = 1.d0/dsqrt(mu(i)) / AMUtoAU
        enddo
        if (present(error_flag)) error_flag=error_local

        return

    end subroutine mu_from_Lcart

    function FC2Freq(Nvib,FC) result(Freq)

        integer,intent(in) :: Nvib
        real(8),dimension(:),intent(in) :: FC
        real(8),dimension(Nvib) :: Freq

        !Local
        integer :: i

        !Transform FC to Freq
        do i=1,Nvib
            Freq(i) = sign(dsqrt(abs(FC(i))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
                                 FC(i))
!             if (FC(i)<0) then
!                 print*, i, FC(i)
!                 call alert_msg("warning","A negative FC found")
!             endif
        enddo

        return

    end function FC2Freq


    function Freq2FC(Nvib,Freq) result(FC)

        integer,intent(in) :: Nvib
        real(8),dimension(:),intent(in) :: Freq
        real(8),dimension(Nvib) :: FC

        !Local
        integer :: i

        !Transform FC to Freq
        do i=1,Nvib
            FC(i) = sign((Freq(i)*2.d0*pi*clight*1.d2)**2/HARTtoJ*BOHRtoM**2*AUtoKG,Freq(i))
!             if (FC(i)<0) then
!                 print*, i, FC(i)
!                 call alert_msg("warning","A negative FC found")
!             endif
        enddo

        return

    end function Freq2FC

    function Hess_to_Hlt(N,Hess) result(Hlt)

        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: Hess
        real(8),dimension(1:(N*(N+1))/2) :: Hlt 

        !Local
        integer :: i,j,k


        k=0
        do i=1,N 
        do j=1,i
            k=k+1
            Hlt(k) = Hess(i,j)
        enddo
        enddo

        return

    end function Hess_to_Hlt

    function Hlt_to_Hess(N,Hlt) result(Hess)

        ! Where 
        ! N=3*Nat 

        integer,intent(in) :: N
        real(8),dimension(:),intent(in) :: Hlt
        real(8),dimension(N,N) :: Hess 

        !Local
        integer :: i,j,k


        k=0
        do i=1,N 
        do j=1,i
            k=k+1
            Hess(i,j) = Hlt(k)
            Hess(j,i) = Hess(i,j)
        enddo
        enddo

        return

    end function Hlt_to_Hess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ANALYSIS


    subroutine analyze_duschinsky(unt,Nvib,G,K,Freq1,Freq2)

        use constants
        use matrix

        integer,intent(in)                :: unt
        integer,intent(in)                :: Nvib
        real(8),dimension(:,:),intent(in) :: G
        real(8),dimension(:)  ,intent(in) :: K
        real(8),dimension(:)  ,intent(in) :: Freq1,Freq2
        
        !Local
        integer :: k90, k95, k99
        integer :: i, j, jj
        real(8),dimension(size(K)) :: AuxVec
        real(8),dimension(size(K)) :: Factor
        integer,dimension(size(K)) :: ipiv
        real(8),dimension(size(K)) ::         Kinv
        real(8),dimension(size(K),size(K)) :: Ginv
        real(8) :: Theta, Qindex

        !Define the Factor to convert into shift into addimensional displacements
        ! for the shift in SI units:
        Factor(1:Nvib) = dsqrt(dabs(Freq1(1:Nvib))*1.d2*SL*2.d0*PI/plankbar)
        ! but we have it in au
        Factor(1:Nvib)=Factor(1:Nvib)*BOHRtoM*dsqrt(AUtoKG)

        write(unt,*) ""
        write(unt,*) "==============================================================================================="
        write(unt,*) " DUSCHINSKI MATRIX (STATE1 WITH RESPECT TO STATE2) (SUMMARY)"
        write(unt,*) "==============================================================================================="
        write(unt,*) "   NM     FREQ1      I2      C2^2       I90     I95     I99       K     K-Dimless   Q-index "
        write(unt,*) "-----------------------------------------------------------------------------------------------"
        do i=1,Nvib
            k90=0
            k95=0
            k99=0
            !Copy the row and reorder 
            AuxVec(1:Nvib) = abs(G(i,1:Nvib))
            call sort_vec_max(AuxVec(1:Nvib),ipiv(1:Nvib),Nvib)
            !Compute quality index
            Qindex = 0.d0
            do j=1,Nvib
                Qindex = Qindex + G(i,j)**2
            enddo
            Theta = 0.d0
            do j=1,Nvib
                if (Theta > 0.9d0) exit
                jj = ipiv(j)
                Theta = Theta + G(i,jj)**2
                k90=k90+1
            enddo 
            Theta = 0.d0
            do j=1,Nvib
                if (Theta > 0.95d0) exit
                jj = ipiv(j)
                Theta = Theta + G(i,jj)**2
                k95=k95+1
            enddo 
            Theta = 0.d0
            do j=1,Nvib
                if (Theta > 0.99d0) exit
                jj = ipiv(j)
                Theta = Theta + G(i,jj)**2
                k99=k99+1
            enddo 
            write(unt,'(X,I5,3X,F8.2,2X,I5,3X,F7.2,5X,3(I5,3X),F7.2,X,F9.4,2X,F9.4)') &
                i, Freq1(i), ipiv(1),G(i,ipiv(1))**2, k90, k95, k99, K(i), K(i)*Factor(i),Qindex
        enddo
        write(unt,*) "==============================================================================================="
        write(unt,*) ""


        if ( verbose==0 ) return

        !============================================
        ! Now repeat for the iniverse transition (if verbose>0)
        !============================================
        Ginv(1:Nvib,1:Nvib) = inverse_realgen(Nvib,G)
        ! Kinv = -Ginv * K
        do i=1,Nvib
            Kinv(i) = 0.d0
            do j=1,Nvib
                Kinv(i) = Kinv(i) - Ginv(i,j)*K(j)
            enddo
        enddo

       !Define the Factor to convert into shift into addimensional displacements
        ! for the shift in SI units:
        Factor(1:Nvib) = dsqrt(dabs(Freq2(1:Nvib))*1.d2*SL*2.d0*PI/plankbar)
        ! but we have it in au
        Factor(1:Nvib)=Factor(1:Nvib)*BOHRtoM*dsqrt(AUtoKG)

        write(unt,*) ""
        write(unt,*) "==============================================================================================="
        write(unt,*) " DUSCHINSKI MATRIX (STATE2 WITH RESPECT TO STATE1) (SUMMARY)"
        write(unt,*) "==============================================================================================="
        write(unt,*) "   NM     FREQ1      I2      C2^2       I90     I95     I99       K     K-Dimless   Q-index "
        write(unt,*) "-----------------------------------------------------------------------------------------------"
        do i=1,Nvib
            k90=0
            k95=0
            k99=0
            !Copy the row and reorder 
            AuxVec(1:Nvib) = abs(Ginv(i,1:Nvib))
            call sort_vec_max(AuxVec(1:Nvib),ipiv(1:Nvib),Nvib)
            !Compute quality index
            Qindex = 0.d0
            do j=1,Nvib
                Qindex = Qindex + Ginv(i,j)**2
            enddo
            Theta = 0.d0
            do j=1,Nvib
                if (Theta > 0.9d0) exit
                jj = ipiv(j)
                Theta = Theta + Ginv(i,jj)**2
                k90=k90+1
            enddo 
            Theta = 0.d0
            do j=1,Nvib
                if (Theta > 0.95d0) exit
                jj = ipiv(j)
                Theta = Theta + Ginv(i,jj)**2
                k95=k95+1
            enddo 
            Theta = 0.d0
            do j=1,Nvib
                if (Theta > 0.99d0) exit
                jj = ipiv(j)
                Theta = Theta + Ginv(i,jj)**2
                k99=k99+1
            enddo 
            write(unt,'(X,I5,3X,F8.2,2X,I5,3X,F7.2,5X,3(I5,3X),F7.2,X,F9.4,2X,F9.4)') &
                i, Freq2(i), ipiv(1),Ginv(i,ipiv(1))**2, k90, k95, k99, Kinv(i), Kinv(i)*Factor(i),Qindex
        enddo
        write(unt,*) "==============================================================================================="
        write(unt,*) ""


    end subroutine analyze_duschinsky


end module vibrational_analysis