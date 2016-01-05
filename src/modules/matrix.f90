module matrix

    !==============================================================
    ! This code is part of FCclasses2
    !==============================================================
    !
    ! Despription
    ! -----------
    ! Subroutines to perform matrix manipulations:
    ! *diagonalization of symmetric matrices
    ! *Lowding orthogonalization
    ! *Determinant
    ! *Matrix printing
    !
    ! Notes
    ! LAPACK subroutines are used
    !==============================================================

    implicit none

    CONTAINS

    function determinant_realsym(N,A) result(det)
    
        !=============================================================
        ! Description
        ! -------------
        ! Computes the determinant by first applying a PLU-like 
        ! decomposition. Actually, the LAPACK symmetric routine.
        ! dsytrf is used performing Bunch-Kaufman diagonal pivoting 
        ! method. It performs the transformation of input A as (*):
        !  A = L * D * L^t
        ! where L is a unit lower triangular matrix multiplied by a
        ! permutation and D is a diagonal matrix with only 1x1 or 
        ! 2x2 non-zero diagonal blocks. Since L is a unit triangular
        ! matrix multiplied by a permutation:
        ! det(L)*det(L^t) = 1
        ! So det(A) = det(D)
        ! The only complications may arise with 2x2 blocks, which are 
        ! marked by negative values of ipiv(i) = ipiv(i+1). Since the 
        ! block is symmetric, only one of the non-diagonal elements is
        ! printed: the D(k+1,k).
        ! Note that it might arrise that ipiv(i+2)=ipiv(i+1), but it
        ! conforms another block. This will confuse an if statement based
        ! only on ipiv(i) and ipiv(i+1), so we explicitely indicate a 
        ! block with a logical variable
        !
        ! (*) this is valid if 'L' is used. For 'U' the actual details change
        !==============================================================

        real(8) :: det 
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8),dimension(N,N) :: AuxArray
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Other
        integer :: alloc_status
        logical :: cycle_next=.false.
        ! Counters
        integer :: i

        ! Copy input matrix
        AuxArray(1:N,1:N) = A(1:N,1:N)

        ! Initialize LAPACK (calculating optimal size first)
        call dsytrf('L', N, AuxArray, N, ipiv, AuxScalar, -1, info)
        lwork = int(AuxScalar)
        allocate(work(1:lwork),stat=alloc_status) 
        if (alloc_status /= 0) then
            write(0,*) "ERROR: Memory cannot be allocated (DETERMINANT ROUTINE)"
            stop
        endif

        ! Perform factorization
        call dsytrf('L', N, AuxArray, N, ipiv, work, lwork, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dsytrf (DETERMINANT ROUTINE)"
            stop
        endif

        ! And compute the determinant
        det = 1.d0
        do i=1,N
            if (cycle_next) then
                cycle_next=.false.
                cycle
            endif
            if (ipiv(i) < 0 .and. ipiv(i) == ipiv(i+1)) then
                ! Compute the 2x2 block
                cycle_next=.true.
                det = det * &
                      ( AuxArray(i,i)*AuxArray(i+1,i+1) - AuxArray(i+1,i)**2 )
             else if (ipiv(i) < 0 .and. ipiv(i) == ipiv(i-1)) then
                ! Second part of a 2x2 block
                cycle
             else
                det = det * AuxArray(i,i)
            endif
        enddo
        
        return

    end function determinant_realsym

    subroutine Log_determinant_realsym(N,A,Log_det,sign_det)
    
        !=============================================================
        ! Description
        ! -------------
        ! Logaritmic version of determinant_realsym()
        ! Returns the log of the absulute value and the sign
        !==============================================================

        real(8),intent(out) :: Log_det, sign_det
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8),dimension(N,N) :: AuxArray
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Other
        integer :: alloc_status
        logical :: cycle_next=.false.
        ! Counters
        integer :: i

        ! Copy input matrix
        AuxArray(1:N,1:N) = A(1:N,1:N)

        ! Initialize LAPACK (calculating optimal size first)
        call dsytrf('L', N, AuxArray, N, ipiv, AuxScalar, -1, info)
        lwork = int(AuxScalar)
        allocate(work(1:lwork),stat=alloc_status) 
        if (alloc_status /= 0) then
            write(0,*) "ERROR: Memory cannot be allocated (DETERMINANT ROUTINE)"
            stop
        endif

        ! Perform factorization
        call dsytrf('L', N, AuxArray, N, ipiv, work, lwork, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dsytrf (DETERMINANT ROUTINE)"
            stop
        endif

        ! And compute the determinant
        Log_det  = 0.d0
        sign_det = 1.d0
        do i=1,N
            if (cycle_next) then
                cycle_next=.false.
                cycle
            endif
            if (ipiv(i) < 0 .and. ipiv(i) == ipiv(i+1)) then
                ! Compute the 2x2 block
                cycle_next=.true.
                AuxScalar = AuxArray(i,i)*AuxArray(i+1,i+1) - AuxArray(i+1,i)**2
                Log_det  = Log_det + Log(abs((AuxScalar)))
                sign_det = sign_det*sign(1.d0,AuxScalar)
             else if (ipiv(i) < 0 .and. ipiv(i) == ipiv(i-1)) then
                ! Second part of a 2x2 block
                cycle
             else
                AuxScalar = AuxArray(i,i)
                Log_det  = Log_det + Log(abs((AuxScalar)))
                sign_det = sign_det*sign(1.d0,AuxScalar)
            endif
        enddo
        
        return

    end subroutine Log_determinant_realsym

    function determinant_realgen(N,A) result(det)
    
        !=============================================================
        ! Description
        ! -------------
        ! Computes the determinant by first applying a PLU-like 
        ! decomposition, using the LAPACK for generic matrix
        !==============================================================

        real(8) :: det 
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8),dimension(N,N) :: AuxArray
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Counters
        integer :: i

        ! Copy input matrix
        AuxArray(1:N,1:N) = A(1:N,1:N)

        ! Perform factorization
        call dgetrf(N, N, AuxArray, N, ipiv, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dgetrf (GEN DETERMINANT ROUTINE)"
            stop
        endif

        ! And compute the determinant
        det = 1.d0
        do i=1,N
            det = det*AuxArray(i,i)
            if (ipiv(i) /= i) det=-det
        enddo
        
        return

    end function determinant_realgen

    function inverse_realsym(N,A) result(Ainv)
    
        !========================================
        ! Description
        ! Computes the determinant by first 
        ! applying a PLU decomposition.
        ! Using LAPACK
        !=======================================

        real(8),dimension(N,N) :: Ainv 
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Other
        integer :: alloc_status
        ! Counters
        integer :: i,j

        ! Copy input matrix
        Ainv(1:N,1:N) = A(1:N,1:N)

        ! Initialize LAPACK (calculating optimal size first)
        call dsytrf('L', N, Ainv, N, ipiv, AuxScalar, -1, info)
        lwork = int(AuxScalar)
        allocate(work(1:lwork),stat=alloc_status) 
        if (alloc_status /= 0) then
            print*, "ERROR: Memory cannot be allocated (LAPACK)"
            stop
        endif

        ! Inverse
        call dsytrf('L', N, Ainv, N, ipiv, work, lwork, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dsytrf (INVERSION ROUTINE)"
            stop
        endif
        call dsytri('L', N, Ainv, N, ipiv, work, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dsytri (INVERSION ROUTINE)"
            stop
        endif
        do i=1,N
        do j=i+1,N
            Ainv(i,j) = Ainv(j,i)
        enddo
        enddo

        return

    end function inverse_realsym

    function inverse_realgen(N,A) result(Ainv)
    
        !========================================
        ! Description
        ! Computes the determinant by first 
        ! applying a PLU decomposition.
        ! Using LAPACK
        !=======================================

        real(8),dimension(N,N) :: Ainv 
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: A
        ! Local
        real(8) :: AuxScalar
        ! LAPACK things
        double complex,dimension(:),allocatable :: work
        integer,dimension(1:N) :: ipiv
        integer :: info, lwork
        ! Other
        integer :: alloc_status

        ! Copy input matrix
        Ainv(1:N,1:N) = A(1:N,1:N)

        ! Perform factorization
        call dgetrf(N, N, Ainv, N, ipiv, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dsytrf (INVERSION ROUTINE)"
            stop
        endif
        ! Initialize LAPACK (calculating optimal size first)
        call dgetri(N, Ainv, N, ipiv, AuxScalar, -1, info)
        lwork = int(AuxScalar)
        allocate(work(1:lwork),stat=alloc_status) 
        if (alloc_status /= 0) then
            print*, "ERROR: Memory cannot be allocated (LAPACK)"
            stop
        endif
        call dgetri(N, Ainv, N, ipiv, work, lwork, info)
        if (info /= 0) then
            write(0,*) "ERROR IN dsytri (INVERSION ROUTINE)"
            stop
        endif

        return

    end function inverse_realgen

    subroutine diagonalize_full(A,N,U,d,alg)
        !
        ! =========================================================
        ! Wrapper code: it calls a given diagonalization SR,
        !               given a lower triangular matrix as input
        !               For symmetric matices only
        ! =========================================================
        !  Input:
        !    A(:,:)  - full matrix to diagonalize. It is not destroyed
        !            2D array dimension: (N,N)
        !    N     - matrix dimension
        !    alg   - string to select the algorith
        !  Output:
        !    U(:,:)- orthogonal matrix with eigenvectors as columns
        !    d(:)  - eigenvalues
        !

        implicit none

        !Input
        character(len=*),intent(in) :: alg
        integer, intent(in) :: N
        real(8),dimension(:,:), intent(in) :: A
        !Output
        real(8),dimension(:,:),intent(out) :: U
        real(8),dimension(:),intent(out) :: d 

        !Local
        real(8), parameter :: PREC=1.d-11
        integer :: IER

        !Auxiliar
        character(len=50) :: dummy_char
        integer :: i,j,k


        select case (trim(adjustl(alg)))
            !The one from the master (M1 year). This is left for testing ONLY
            case("emtccm")
              !We need to feed the whole matrix, as we have
              call JACOBI(N,PREC,A,U,d)

! #ifdef USE_LAPACK
            !Only available if compiled in duble precision and if lapack libs are required
            case("lapack")
              !We need to feed the whole matrix, as we have
              !We don't destroy the original
              U(1:N,1:N)=A(1:N,1:N)
              call diasym(U(1:N,1:N),d(1:N),N)
! #endif

            case default
              print*, "Unsupported diagonalization algorith:"//alg

        end select   

        return      

    end subroutine diagonalize_full

! #ifdef USE_LAPACK
    subroutine diasym(a,eig,n)
    !Taken from http://physics.bu.edu/py502/lectures4/examples/diatest.f90
    ! needs lapack!
    !---------------------------------------------------------!
    !Calls the LAPACK diagonalization subroutine DSYEV        !
    !input:  a(n,n) = real symmetric matrix to be diagonalized!
    !            n  = size of a                               !
    !output: a(n,n) = orthonormal eigenvectors of a           !
    !        eig(n) = eigenvalues of a in ascending order     !
    !---------------------------------------------------------!
        implicit none 

        integer,intent(in) :: n
        real(8),dimension(N,N),intent(inout) :: a
        real(8),dimension(N),intent(out) :: eig
        !Local
        integer :: l,inf
        real(8),dimension(1:1) :: work_estimate
        real(8),dimension(:),allocatable :: work 
        integer :: i

        !Needs LAPACK
        external dsyev

!          print'(/,X,A,/)', "Entering lapack diagonalization subroutine for symmetric matrices..."

        ! First estimate of the WORK size
        call dsyev('V','U',n,a,n,eig,work_estimate,-1,inf)
        ! Now we use the value estimated by LAPACK
        ! which is stored in work_estimate(1)
        l=work_estimate(1)
        allocate(work(1:l))
        call dsyev('V','U',n,a,n,eig,work,l,inf)
        if (inf /= 0) then
            print'(2X,A,I0,A,/)', "ERROR in diagonalization. (LAPACK error code: ", inf, ")"
            stop
!          else  
!              print'(2X,A,/)', "Successfull diagonalization"
        endif

    end subroutine diasym
! #endif


    subroutine JACOBI(N,RHO,FF,V,d)
    
    !=====================================================================
    ! Description
    ! Program to diagonalize a real symmetric matrix
    ! Source: EMTCCM - M1. Pais Vasco, 2010
    !
    ! Arguaments
    ! FF  is the matrix to be diagonalized (it is not modified)
    ! V   is the matrix of eigenvectors (by column)
    ! d   is the vector of eigenvalues
    ! RHO is the convergence criterium. The process ends successfully
    !     if TE<RHO, where TE is the norm of non-diagonal elements
    !
    ! Notes
    ! There is not maximum number of iterations (infinite)
    !=====================================================================
      
        implicit none
        
        integer,intent(in)::N
        real(8),dimension(:,:),intent(in)::FF
        real(8),dimension(:,:),intent(out)::V
        real(8),dimension(:),intent(out)::d
        real(8),intent(in)::RHO
        !Local
        real(8),dimension(1:N,1:N)::F !this is an auxiliar now
        real(8) :: A, COST, TE, TEN, omega, sint, U, V1, V2, V3, Z, TEM
        integer :: iter
        integer :: i, ii, ij, j, jj
        
        logical::reduced,converged
    
        print'(/,X,A)', "Entering a not so efficient diagonalization subroutine for symmetric matrices"
        print'(2X,A,E13.6,/)', "Threshol in norm:", RHO
        print'(2X,A)', "ITERATIONS..."
    
        
        ! Copy matrix to local auxiliar
        F(1:N,1:N)=FF(1:N,1:N)
        
        !C...INICIALIZACION DE LA MATRIZ DE VECTORES PROPIOS
        do I=1,N
          do J=1,N
            V(I,J)=0.0D0
          enddo
          V(I,I)=1.0D0
        enddo
        
        !C...INICIALIZACION DE LAS VARIABLES DE LA ITERACION
        A=dfloat(N)
        ITER=0
        call RMS(TE,F,N)
        TEN=TE/A
        
            
        !C...PROCESO ITERATIVO
        !Cada ciclo de iteración reduce la norma de los elementos no diagonales de forma "equilibrada"
        write(6,'(3X,A,I0)')     ' ITERATION: ',ITER
        converged=.false. 
        do while (.not.converged)
        
            ITER=ITER+1
            write(6,'(4X,A,E13.6,/)')' ERROR:     ',TE
            write(6,'(3X,A,I0)')     ' ITERATION: ',ITER
        
            reduced=.false.
            do while (.not.reduced)
            !Cada ciclo reduce el valor de los elementos que se pasen del valor medio de la norma de los elementos fuera de la diagonal
            !el ciclo de iteración se acaba cuando hayamos reducido convenientemente el valor de todos los elementos diagonales, tendremos que hacer varias pasadas
            !la variable MA controlaba si se ha producido-->cambiada por la variable lógica reduced.
        
            reduced=.true.
        
            !Buscamos en todos los elementos II,JJ no diagonales (solo la mitad, ya que es simétrica) para reducirlos 
            !hasta que todos son menores que la media de la norma (criterio adoptado)
            do II=2,N !14
                IJ=II-1
                do JJ=1,IJ !14-2
        
                !Comprobamos si el elemento no diagonal (II,JJ) cumple el criterio adoptado: si lo cumple no se aplica el algoritmo
                !si alguno no lo cumple, se repite el ciclo hasta que todos lo cumplan.
                if (DABS(F(II,JJ)).gt.TEN) then
    
                  reduced=.false.  !obliga a repetir el ciclo de reducción
            
                  !Algoritmo de Jacobi (reducción de los elementos no diagonales):
                  V1=F(JJ,JJ)
                  V2=F(II,JJ)
                  V3=F(II,II)
                  U=.5D0*(F(JJ,JJ)-F(II,II))
      
                  if (DABS(U).lt.1.d-10) then
                        OMEGA=-1.0D0
                  else
    
                    OMEGA=-F(II,JJ)/DSQRT(F(II,JJ)*F(II,JJ)+U*U)
                    Z=1.D0
                    if(U.LT.0.D0) Z=-Z
                      OMEGA=OMEGA*Z
                    endif
        
                  SINT=OMEGA/DSQRT(2.D0*(1.D0+DSQRT(1.D0-OMEGA*OMEGA)))
                  COST=DSQRT(1.D0-SINT*SINT)
        
                  !Alteramos todos los elementos de la matriz en las filas y columnas II y JJ
                  do I=1,N !13
                      !Operamos sobre la matriz F...
                      if(I.ge.II) then
                        TEM=F(I,JJ)*COST-F(I,II)*SINT
                        F(I,II)=F(I,JJ)*SINT+F(I,II)*COST
                        F(I,JJ)=TEM
                      else
                          if(I.lt.JJ) then
                              TEM=F(JJ,I)*COST-F(II,I)*SINT
                              F(II,I)=F(JJ,I)*SINT+F(II,I)*COST
                              F(JJ,I)=TEM
                          else
                              TEM=F(I,JJ)*COST-F(II,I)*SINT
                              F(II,I)=F(I,JJ)*SINT+F(II,I)*COST
                              F(I,JJ)=TEM
                          endif
                      endif
                      !Y actualizamos los vecotores propios V...
                      TEM=V(I,JJ)*COST-V(I,II)*SINT
                      V(I,II)=V(I,JJ)*SINT+V(I,II)*COST
                      V(I,JJ)=TEM
                  enddo !13
        
                  !Valor de los elementos invloucrados en el paso (II&JJ)
                  F(JJ,JJ)=V1*COST*COST+V3*SINT*SINT-2.D0*V2*SINT*COST
                  F(II,II)=V1*SINT*SINT+V3*COST*COST+2.D0*V2*SINT*COST
                  F(II,JJ)=2.D0*U*SINT*COST+V2*(COST*COST-SINT*SINT) 
    
                endif
            
              enddo !14-2
            enddo !14
        
          enddo !do while interno (ciclos de reducción)
        
        !Actualiza los valores de la norma y la norma media
        call RMS(TE,F,N)
        TEN=TE/A
        !Comprueba si ha convergido
        IF(TE.LT.RHO) then
            write(6,'(4X,A,E13.6,/)')' ERROR:     ',TE
            print'(2X,A,/)', "Successfull diagonalization"
            converged=.true.
        endif
        
        enddo !do while externo (convergencia global)
    
        d(1:N) = (/ (F(i,i), i=1,N ) /)
        
        return
    
    end subroutine JACOBI
    
    subroutine RMS(TE,F,N)

        ! This is used by JACOBI 
        
        integer,intent(in) :: N
        real(8),dimension(:,:),intent(in) :: F
        real(8),intent(out) :: TE
        !Local
        integer :: i, j, k
        
        TE=0.0D0
        do I=2,N
          K=I-1
          do J=1,K
            TE=TE+2.D0*F(I,J)*F(I,J)
          enddo
        enddo
        
        TE=DSQRT(TE)
        
        return
    
    end subroutine RMS

    subroutine orth_Lowdin(Nr,Nc,T)

        !======================================
        ! Description
        ! ------------
        ! Performs the orthogonalization of the
        ! matrix T(NrxNc). The orthogonal eigenvectors
        ! correspond to the columns: T(:,Nc).
        !
        ! Uses some Subroutines from this module
        !=======================================

        integer,intent(in) :: Nr, Nc
        real(8),dimension(:,:),intent(inout) :: T
        !Local
        real(8),dimension(1:Nr,1:Nr) :: AA
        real(8),dimension(1:Nr,1:Nr) :: BB
        real(8),dimension(1:Nr)      :: V
        !
        integer :: i, j, k, ii, jj, kk
        real(8),dimension(1:Nr,1:Nr) :: p


        ! Compute the metric matrix of T
        do i=1,Nc 
        do j=1,Nc 
            AA(i,j) = 0.d0
            do k=1,Nr
                AA(i,j)=AA(i,j)+T(k,i)*T(k,j)
            enddo
        enddo
        enddo
        ! Diagonalize metric matrix
        call diagonalize_full(AA,Nc,BB,V,"lapack")
! c     compute s^-0.5
        do ii=1,Nc
        do jj=1,Nc
        p(ii,jj)=0.d0
        do kk=1,Nc
        p(ii,jj)=p(ii,jj)+BB(jj,kk)/dsqrt(V(kk))*BB(ii,kk)
        enddo
        enddo
        enddo
! c     compute t1*s^(-1/2)
        do ii=1,Nr
        do jj=1,Nc
        BB(ii,jj)=0.d0
        do kk=1,Nc
        BB(ii,jj)=BB(ii,jj)+T(ii,kk)*p(kk,jj)
        enddo
        enddo
        enddo 
        do ii=1,Nr
        do jj=1,Nc
        T(ii,jj)=BB(ii,jj)
        enddo
        enddo

        return

    end subroutine orth_Lowdin

    !-----------------------------------------
    ! WRAPPER FUNCTION TO BLAS MATRIX PRODUCTS
    !-----------------------------------------

    function matrix_product(NA,NB,NK,A,B,tA,tB) result(P)

        ! A Wrapper to dgemm to multiply 
        ! matrices
        ! P(NA,NB) = A(NA,K) * B(K,NB)
        ! P(NA,NB) = A^t(K,NA) * B(K,NB)
        ! and modificaitions alike

        integer,intent(in)                   :: NA,NB,NK
        real(8),dimension(:,:),intent(in)    :: A,B
        logical,intent(in),optional          :: tA,tB
        real(8),dimension(NA,NB)             :: P
        !Local
        character :: opA,opB
        integer :: LDA, LDB

        !Needs BLAS
        external dgemm

        !First dimension as specified in the calling program
        ! equal to NA and NB if opA,opB = 'N'
        LDA = size(A,1)
        LDB = size(B,1)
        opA = 'N'
        opB = 'N'

        if (present(tA).and.tA) opA='T'
        if (present(tB).and.tB) opB='T'

        call dgemm(opA,opB,NA,NB,NK,1.d0,A,LDA,B,LDB,0.d0,P,NA)  
        return

!         if (opA == 'N'.and. opB == 'N') then    
!             do i=1,NA
!             do j=1,NB
!                 P(i,j) = 0.d0
!                 do k=1,NK
!                     P(i,j) = P(i,j) + A(i,k) * B(k,j)
!                 enddo
!             enddo
!             enddo
!         else if (opA == 'T'.and. opB == 'N') then  
!             do i=1,NA
!             do j=1,NB
!                 P(i,j) = 0.d0
!                 do k=1,NK
!                     P(i,j) = P(i,j) + A(k,i) * B(k,j)
!                 enddo
!             enddo
!             enddo
!         else if (opA == 'N'.and. opB == 'T') then  
!             do i=1,NA
!             do j=1,NB
!                 P(i,j) = 0.d0
!                 do k=1,NK
!                     P(i,j) = P(i,j) + A(i,k) * B(j,k)
!                 enddo
!             enddo
!             enddo
!         else if (opA == 'T'.and. opB == 'T') then  
!             do i=1,NA
!             do j=1,NB
!                 P(i,j) = 0.d0
!                 do k=1,NK
!                     P(i,j) = P(i,j) + A(k,i) * B(j,k)
!                 enddo
!             enddo
!             enddo
!         endif

        return

    end function matrix_product

    function matrix_basisrot(M,N,X,A,counter) result(P)

        !=============================================
        ! Rotation by X of basis set A: 
        ! NORMAL (counter=.false.)
        ! P = X A X^t, where X(M,N) and A(N,N)
        ! INVERSE (counter=.true.)
        ! P = X^t A X, where X(N,M) and A(N,N)
        !=============================================

        integer,intent(in)                   :: M,N
        real(8),dimension(:,:),intent(in)    :: X,A
        logical,intent(in),optional          :: counter
        real(8),dimension(M,M)               :: P
        !Local
        real(8),dimension(M,N)               :: Aux

        if (present(counter) .and. counter) then
            Aux = matrix_product(M,N,N,X,A,tA=.true.)
            P   = matrix_product(M,M,N,Aux,X)
        else
            Aux = matrix_product(M,N,N,X,A)
            P   = matrix_product(M,M,N,Aux,X,tB=.true.)
        endif

        return

    end function matrix_basisrot

    function diag_basisrot(M,N,X,a,counter) result(P)

        !=============================================
        ! Rotation by X of basis set T: 
        ! NORMAL (counter=.false.)
        ! P = X A X^t, where X(M,N) and T(N,N)
        ! INVERSE (counter=.true.)
        ! P = X^t A X, where X(N,M) and T(N,N)
        !
        ! A is diagonal, and the subroutine uses 
        ! the vector a=diag(A) as input
        !=============================================

        integer,intent(in)                   :: M,N
        real(8),dimension(:,:),intent(in)    :: X
        real(8),dimension(:),intent(in)      :: a
        logical,intent(in),optional          :: counter
        real(8),dimension(M,M)               :: P
        !Local
        real(8),dimension(M,N)               :: Aux
        integer :: j

        if (present(counter) .and. counter) then
        ! NOT TESTED
            do j=1,N
                Aux(:,j) = X(j,1:M)*a(j)
            enddo
            P   = matrix_product(M,M,N,Aux,X)
        else
            do j=1,N
                Aux(:,j) = X(1:M,j)*a(j)
            enddo
            P   = matrix_product(M,M,N,Aux,X,tB=.true.)
        endif


        return

    end function diag_basisrot


     ! SORTING ROUTINES

    subroutine sort_diag(N,D,V)
    !=========================================================
    ! Sorts eigenvalues (from Lower to Higher) in diag(D) and
    ! applies the same permutations to the eigenvector matrix (V)
    !=========================================================
    
            implicit none

            INTEGER, intent(in):: N
            real(8), DIMENSION(:,:),intent(inout) :: D,V

            real(8) :: Amin, aux
            real(8),DIMENSION(1:N) :: aux2
                    
            integer :: i,j,imin
                    
            do  i=1,N-1
              Amin=D(i,i)
              imin=i
              
              do j=i+1,N
                if (D(j,j).lt.Amin) then
                 Amin=D(j,j)
                 imin=j
                endif
              enddo

              aux=D(i,i)
              aux2=V(1:N,i)
              D(i,i)=D(imin,imin)
              V(1:N,i)=V(1:N,imin)
              D(imin,imin)=aux
              V(1:N,imin)=aux2
            enddo

            RETURN
    END SUBROUTINE sort_diag


    SUBROUTINE sort_diag_rev(N,D,V)

    !=========================================================
    ! Sorts eigenvalues in reverse order (from Higher to Lower) 
    ! in diag(D) and applies the same permutations to the 
    ! eigenvector matrix (V)
    !=========================================================
    
            implicit none

            INTEGER, intent(in):: N
            real(8), DIMENSION(:,:),intent(inout) :: D,V

            real(8) :: Amax, aux
            real(8),DIMENSION(1:N) :: aux2
                    
            integer :: i,j,imax
                    
            do  i=1,N-1
              Amax=D(i,i)
              imax=i
              
              do j=i+1,N
                if (D(j,j).gt.Amax) then
                 Amax=D(j,j)
                 imax=j
                endif
              enddo

              aux=D(i,i)
              aux2=V(1:N,i)
              D(i,i)=D(imax,imax)
              V(1:N,i)=V(1:N,imax)
              D(imax,imax)=aux
              V(1:N,imax)=aux2
            enddo

            RETURN

    END SUBROUTINE sort_diag_rev

    !Simple routines
    subroutine sort_vec(V,N,IORD,reverse)

        ! From min to max

        implicit none

        real(8),dimension(:),intent(inout) :: V
        integer,intent(in) :: N
        integer,dimension(:),intent(inout),optional :: IORD
        logical,intent(in),optional :: reverse

        real(8) :: aux
        integer :: i,j
        logical :: rev=.false.

        if (present(reverse)) then
            rev=reverse
        endif

        do i=1,N-1
            do j=i+1,N
                ! Normal ordering (Low to high)
                if (V(j)<V(i) .and. .not.rev) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                ! Reverse ordering (High to low)
                elseif (V(j)>V(i) .and. rev) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                endif
            enddo
        enddo

        return

    end subroutine sort_vec         

    subroutine sort_vec_max(V,IORD,N)

        !Order from max to min and track indices in IORD matrix

        implicit none

        real(8),dimension(:),intent(inout) :: V
        real(8) :: aux
        integer,dimension(:),intent(inout) :: IORD
        integer,intent(in) :: N
        integer :: i,j, iaux

        !Intialize IORD
        do i=1,N
            IORD(i) = i
        enddo

        do i=1,N-1
            do j=i+1,N
                if (V(j)>V(i)) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                    !Track the index permutations in IORD
                    iaux = IORD(i)
                    IORD(i) = IORD(j)
                    IORD(j) = iaux
                endif
            enddo
        enddo

        return

    end subroutine sort_vec_max        


    subroutine sort_vec_int(V,N,reverse)

        !Integer version of sort_vec subroutine

        implicit none

        integer,dimension(:),intent(inout) :: V
        integer,intent(in) :: N
        logical,intent(in),optional :: reverse

        integer :: aux
        integer :: i,j

        logical :: rev=.false.

        if (present(reverse)) then
            rev=reverse
        endif

        do i=1,N-1
            do j=i+1,N
                ! Normal ordering (Low to high)
                if (V(j)<V(i) .and. .not.rev) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                ! Reverse ordering (High to low)
                elseif (V(j)>V(i) .and. rev) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                endif
            enddo
        enddo

    end subroutine sort_vec_int    

    
    subroutine rotation_3D(vx,vy,vz,tx,ty,tz,Theta) 
    
        !Description:
        ! Subroutine to rotate the vector (vx,vy,vz) around the axis
        ! defined by (tx,ty,tz) an angle Theta (rad).
    
        real(8), intent(inout) :: vx,vy,vz
        real(8), intent(in)    :: tx,ty,tz, Theta
    
        !Local
        real(8),dimension(1:3,1:3) :: R
        real(8) :: vx_tmp, vy_tmp, tmod
        real(8) :: ux, uy, uz 
    
        ! Vector u must be unitary
        tmod = sqrt(tx**2 + ty**2 + tz**2)
        ux = tx/tmod
        uy = ty/tmod
        uz = tz/tmod
    
        ! Form 3D-rotation matrix (from Wikipedia)
        R(1,1) = cos(Theta) + ux**2*(1.0-cos(Theta))
        R(1,2) = ux*uy*(1.0-cos(Theta)) - uz*sin(Theta)
        R(1,3) = ux*uz*(1.0-cos(Theta)) + uy*sin(Theta)
        R(2,1) = ux*uy*(1.0-cos(Theta)) + uz*sin(Theta)
        R(2,2) = cos(Theta) + uy**2*(1.0-cos(Theta))
        R(2,3) = uy*uz*(1.0-cos(Theta)) - ux*sin(Theta)
        R(3,1) = ux*uz*(1.0-cos(Theta)) - uy*sin(Theta)
        R(3,2) = uy*uz*(1.0-cos(Theta)) + ux*sin(Theta)
        R(3,3) = cos(Theta) + uz**2*(1.0-cos(Theta))
    
        ! Apply rotaion
        vx_tmp = vx*R(1,1) + vy*R(1,2) + vz*R(1,3)
        vy_tmp = vx*R(2,1) + vy*R(2,2) + vz*R(2,3)
        vz =     vx*R(3,1) + vy*R(3,2) + vz*R(3,3)
        vx = vx_tmp 
        vy = vy_tmp 
    
       return
    
    end subroutine rotation_3D

end module matrix
