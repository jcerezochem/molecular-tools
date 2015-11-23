module MatrixMod

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS
    !==============================================================
    !
    ! Despription
    ! -----------
    ! Subroutines to perform matrix manipulations
    !==============================================================


!MODULO CON LA SUBRUTINA DE DIAGONALIZACIÓN DE JACOBI Y RELACIONADAS
!Todas las funciones y subrutinas están definidas en doble precisión
!Fuente: material del curso intensivo del Master en Química Teórica y Modelización Computacional (Semana 2). San Sebastián, Enero-Febrero 2010 
! Original en F77 --> Reescrito en F90. Se reescribe en f90 para poder compilarlo junto al programa principal sin tener problemas (gfortran). 
! Otra opción con este compilador hubiera sido compilar la subrutina por separado (crear un objeto). El objeto se enlaza (link) al compilar el 
! programa principal.

    CONTAINS

    subroutine diagonalize_lt(A,N,U,d,alg)
        !
        ! =========================================================
        ! Wrapper code: it calls a given diagonalization SR,
        !               given a lower triangular matrix as input
        !               For symmetric matices only
        ! =========================================================
        !  Input:
        !    A(:)  - lt matrix to diagonalize. It is not destroyed
        !            1D array dimension: (N*(N+1)/2)
        !    N     - matrix dimension
        !    alg   - string to select the algorith
        !  Output:
        !    U(:,:)- orthogonal matrix with eigenvectors as columns
        !    d(:)  - eigenvalues
        !

        use alerts

        implicit none

        !Input
        character(len=*),intent(in) :: alg
        integer, intent(in) :: N
#ifdef DOUBLE
        double precision,dimension(:), intent(in) :: A
        !Output
        double precision,dimension(:,:),intent(out) :: U
        double precision,dimension(:),intent(out) :: d
#else
        real,dimension(:), intent(in) :: A
        !Output
        real,dimension(:,:),intent(out) :: U
        real,dimension(:),intent(out) :: d
#endif   

        !Local
#ifdef DOUBLE
        double precision, dimension(N,N) :: A_full
        double precision, dimension( (N*(N+1))/2 ) :: WK
        double precision, dimension( (N*(N+1))/2 ) :: AA
        double precision, parameter :: PREC=1.d-11
#else
        real, dimension(N,N) :: A_full
        real, dimension( (N*(N+1))/2 ) :: WK
        real, dimension( (N*(N+1))/2 ) :: AA
        real, parameter :: PREC=1.d-11
#endif
        integer :: IER

        !Auxiliar
        character(len=50) :: dummy_char
        integer :: i,j,k

        

        select case (trim(adjustl(alg)))
            !The one from the master (M1 year)
            case("emtccm")
              !We need to feed the whole matrix
              k=0
              do i=1,N
                  do j=1,i
                      k=k+1
                      A_full(i,j) = A(k)
                      A_full(j,i) = A(k)
                  enddo
              enddo
              call JACOBI_SIM(N,PREC,A_full,U,d)

            !The IMSL eigrs rewritten to F90
            case("eigrs")
              !We don't want to destroy our valuable initial matrix
              AA=A
              call eigrs(AA,N,2,d,U,N,WK,IER)
              write(dummy_char,*) IER
              if (IER > 128) call alert_msg("fatal","Diagonalization with eigrs not converged: "//trim(adjustl(dummy_char)))
              write(dummy_char,*) WK(1)
              if (WK(1) > 1.d0 .and. WK(1) < 100.d0 ) call alert_msg("note","Diagonalization is just satisfactory: PI = "&
                                                                      //trim(adjustl(dummy_char)))
              if (WK(1) > 100.d0) call alert_msg("fatal","Poor diagonalization: PI = "//trim(adjustl(dummy_char)))

            case default
              call alert_msg("fatal","Error using the diagonalize_lt subroutine.")

        end select   

        return      

    end subroutine diagonalize_lt


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

        use alerts

        implicit none

        !Input
        character(len=*),intent(in) :: alg
        integer, intent(in) :: N
#ifdef DOUBLE
        double precision,dimension(:,:), intent(in) :: A
        !Output
        double precision,dimension(:,:),intent(out) :: U
        double precision,dimension(:),intent(out) :: d
#else
        real,dimension(:,:), intent(in) :: A
        !Output
        real,dimension(:,:),intent(out) :: U
        real,dimension(:),intent(out) :: d
#endif   

        !Local
#ifdef DOUBLE
        double precision, dimension( (N*(N+1))/2 ) :: A_lt
        double precision, dimension( (N*(N+1))/2 ) :: WK
        double precision, parameter :: PREC=1.d-11
#else
        real, dimension( (N*(N+1))/2 ) :: A_lt
        real, dimension( (N*(N+1))/2 ) :: WK
        real, parameter :: PREC=1.d-11
#endif
        integer :: IER

        !Auxiliar
        character(len=50) :: dummy_char
        integer :: i,j,k

        

        select case (trim(adjustl(alg)))
            !The one from the master (M1 year)
            case("emtccm")
              !We need to feed the whole matrix, as we have
              call JACOBI_SIM(N,PREC,A,U,d)

#ifdef DOUBLE
#ifdef USE_LAPACK
            !Only available if compiled in duble precision and if lapack libs are required
            case("lapack")
              !We need to feed the whole matrix, as we have
              !We don't destroy the original
              U(1:N,1:N)=A(1:N,1:N)
              call diasym(U,d,N)
#endif
#endif

            !The IMSL eigrs rewritten to F90
            case("eigrs")
              !We don't want to destroy our valuable initial matrix
              !and need a lt matrix
              k=0
              do i=1,N
              do j=1,i
                  k=k+1
                  A_lt(k) = A(i,j)
              enddo
              enddo
              call eigrs(A_lt,N,2,d,U,N,WK,IER)
              write(dummy_char,*) IER
              if (IER > 128) call alert_msg("fatal","Diagonalization with eigrs not converged: "//trim(adjustl(dummy_char)))
              write(dummy_char,*) WK(1)
              if (WK(1) > 1.d0 .and. WK(1) < 100.d0 ) call alert_msg("note","Diagonalization is just satisfactory: PI = "&
                                                                      //trim(adjustl(dummy_char)))
              if (WK(1) > 100.d0) call alert_msg("fatal","Poor diagonalization: PI = "//trim(adjustl(dummy_char)))

            case default
              call alert_msg("fatal","Diagonalize_full subroutine called with an unsupported option: "//trim(adjustl(alg)))

        end select   

        return      

    end subroutine diagonalize_full

!===============================================
    SUBROUTINE JACOBI(A,N,NP,D,V,NROT) 
! From Numerical Recipes (only changed CONTINUEs for ENDDOs, and DIMENSION..)
!
! Description
! Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
! by n, stored in a physical np by np (this featured changed by the use of implicit dimension)
!  array. On output, elements of a above the diagonal are
! destroyed. d returns the eigenvalues of a in its first n elements. v is a matrix with the same
! logical and physical dimensions as a, whose columns contain, on output, the normalized
! eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
#ifdef DOUBLE
    implicit double precision (A-H,O-Z)
#else
    implicit real (A-H,O-Z)
#endif

    !In NR the limit was 100. Maybe indicates that this algorith is not so good..
    ! (eigrs is not using it, so its faster...)
    integer,PARAMETER :: NMAX=500

    integer,intent(in) :: N
    integer,intent(out) :: NROT
#ifdef DOUBLE
    double precision,dimension(NP,NP),intent(inout) :: A
    double precision,dimension(NP,NP),intent(out) :: V
    double precision,dimension(NP),intent(out) :: D
    !Local
    double precision,dimension(NMAX) :: B,Z
#else
    real,dimension(NP,NP),intent(inout) :: A
    real,dimension(NP,NP),intent(out) :: V
    real,dimension(NP),intent(out) :: D
    !Local
    real,dimension(NMAX) :: B,Z
#endif

 
    DO IP=1,N 
      DO IQ=1,N 
        V(IP,IQ)=0. 
      ENDDO 
      V(IP,IP)=1. 
    ENDDO 
    DO IP=1,N 
      B(IP)=A(IP,IP) 
      D(IP)=B(IP) 
      Z(IP)=0. 
    ENDDO 
    NROT=0 
    DO I=1,50 
      SM=0. 
      DO IP=1,N-1 
        DO IQ=IP+1,N 
#ifdef DOUBLE
          SM=SM+DABS(A(IP,IQ)) 
#else
          SM=SM+ABS(A(IP,IQ)) 
#endif
        ENDDO 
      ENDDO 
      IF(SM.EQ.0.d0)RETURN 
      IF(I.LT.4)THEN 
        TRESH=0.2*SM/N**2 
      ELSE 
        TRESH=0. 
      ENDIF 
      DO IP=1,N-1 
        DO IQ=IP+1,N 
#ifdef DOUBLE
          G=100.*DABS(A(IP,IQ)) 
          IF((I.GT.4).AND.(DABS(D(IP))+G.EQ.DABS(D(IP))) &
             .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ))))THEN 
            A(IP,IQ)=0. 
          ELSE IF(DABS(A(IP,IQ)).GT.TRESH)THEN 
            H=D(IQ)-D(IP) 
            IF(DABS(H)+G.EQ.DABS(H))THEN 
              T=A(IP,IQ)/H 
            ELSE 
              THETA=0.5*H/A(IP,IQ) 
              T=1./(DABS(THETA)+DSQRT(1.+THETA**2)) 
              IF(THETA.LT.0.)T=-T 
            ENDIF 
            C=1./DSQRT(1+T**2) 
#else
          G=100.*ABS(A(IP,IQ)) 
          IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP))) &
             .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN 
            A(IP,IQ)=0. 
          ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN 
            H=D(IQ)-D(IP) 
            IF(ABS(H)+G.EQ.ABS(H))THEN 
              T=A(IP,IQ)/H 
            ELSE 
              THETA=0.5*H/A(IP,IQ) 
              T=1./(ABS(THETA)+SQRT(1.+THETA**2)) 
              IF(THETA.LT.0.)T=-T 
            ENDIF 
            C=1./SQRT(1+T**2) 
#endif
            S=T*C 
            TAU=S/(1.+C) 
            H=T*A(IP,IQ) 
            Z(IP)=Z(IP)-H 
            Z(IQ)=Z(IQ)+H 
            D(IP)=D(IP)-H 
            D(IQ)=D(IQ)+H 
            A(IP,IQ)=0. 
            DO J=1,IP-1 
              G=A(J,IP) 
              H=A(J,IQ) 
              A(J,IP)=G-S*(H+G*TAU) 
              A(J,IQ)=H+S*(G-H*TAU) 
            ENDDO 
            DO J=IP+1,IQ-1 
              G=A(IP,J) 
              H=A(J,IQ) 
              A(IP,J)=G-S*(H+G*TAU) 
              A(J,IQ)=H+S*(G-H*TAU) 
            ENDDO 
            DO J=IQ+1,N 
              G=A(IP,J) 
              H=A(IQ,J) 
              A(IP,J)=G-S*(H+G*TAU) 
              A(IQ,J)=H+S*(G-H*TAU) 
            ENDDO 
            DO J=1,N 
              G=V(J,IP) 
              H=V(J,IQ) 
              V(J,IP)=G-S*(H+G*TAU) 
              V(J,IQ)=H+S*(G-H*TAU) 
            ENDDO 
            NROT=NROT+1 
          ENDIF 
        ENDDO 
      ENDDO 
      DO IP=1,N 
        B(IP)=B(IP)+Z(IP) 
        D(IP)=B(IP) 
        Z(IP)=0. 
      ENDDO 
    ENDDO 
    STOP '50 iterations should never happen' 
    RETURN 
    END SUBROUTINE JACOBI



#ifdef USE_LAPACK
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
         integer l,inf
         real*8  work(n*(3+n/2)) 

!          print'(/,X,A,/)', "Entering lapack diagonalization subroutine for symmetric matrices"

         l=n*(3+n/2)
         call dsyev('V','U',n,a,n,eig,work,l,inf)

     end subroutine diasym
#endif


!C-----------------------------------------------------------------------
SUBROUTINE JACOBI_SIM(N,RHO,FF,V,d)
!C-----------------------------------------------------------------------
!CPROGRAMA PARA DIGONALIZAR UNA MATRIZ REAL Y SIMETRICA DE ORDEN N
!C
!C  F ES LA MATRIZ A DIAGONALIZAR. (NO SE MODIFICA -- nueva versión)
!C  V ES LA MATRIZ DE VECTORES PROPIOS (CADA COLUMNA UN VECTOR)
!C  d  es el vector de valores propios
!C  RHO ES EL CRITERIO DE CONVERGENCIA. EL PROCESO FINALIZA SI TE<RHO
!C    DONDE TE ES LA NORMA DE LOS ELEMENTOS NO DIAGONALES.
!C 

    !use matmod !Para la escritura de las matrices con WriteMatrix
      
#ifdef DOUBLE
    implicit double precision (A-H,O-Z)
#else
    implicit real (A-H,O-Z)
#endif
    implicit integer (i-n)
    
    integer,intent(in)::N
#ifdef DOUBLE
    double precision,dimension(:,:),intent(inout)::V
    double precision,dimension(:,:),intent(in)::FF
    double precision,dimension(:),intent(out)::d
    double precision,intent(in)::RHO
    double precision,dimension(N,N)::F !this is an auxiliar now
#else
    real,dimension(:,:),intent(inout)::V
    real,dimension(:,:),intent(in)::FF
    real,dimension(:),intent(out)::d
    real,dimension(N,N)::F !this is an auxiliar now
    real,intent(in)::RHO
#endif
    
    logical::reduced,converged

    print'(/,X,A,/)', "Entering a not so efficient diagonalization subroutine for symmetric matrices"
    
    ! Copy matrix to local auxiliar
    F=FF
    
    !C...INICIALIZACION DE LA MATRIZ DE VECTORES PROPIOS
    do I=1,N
      do J=1,N
        V(I,J)=0.0D0
      enddo
      V(I,I)=1.0D0
    enddo
    
    !C...INICIALIZACION DE LAS VARIABLES DE LA ITERACION
    A=N
    ITER=0
    call RMS(TE,F,N)
    TEN=TE/A
    
        
    !C...PROCESO ITERATIVO
    !Cada ciclo de iteración reduce la norma de los elementos no diagonales de forma "equilibrada"
    converged=.false. 
    do while (.not.converged)
    
        ITER=ITER+1
!         write(6,*)' ITERACION ..',ITER
!         write(6,*)' ERROR ......',TE
    
        reduced=.false.
        do while (.not.reduced)
        !Cada ciclo reduce el valor de los elementos que se pasen del valor medio de la norma de los elementos fuera de la diagonal
        !el ciclo de iteración se acaba cuando hayamos reducido convenientemente el valor de todos los elementos diagonales, tendremos que hacer varias pasadas
        !la variable MA controlaba si se ha producido-->cambiada por la variable lógica reduced.
    
        reduced=.true.
    
        !Buscamos en todos los elementos II,JJ no diagonales (solo la mitad, ya que es simétrica) para hacer reducirlos 
        !hasta que todos son menores que la media de la norma (criterio adoptado)
        do II=2,N !14
            IJ=II-1
            do JJ=1,IJ !14-2
    
            !Comprobamos si el elemento no diagonal (II,JJ) cumple el criterio adoptado: si lo cumple no se aplica el algoritmo
            !si alguno no lo cumple, se repite el ciclo hasta que todos lo cumplan.
#ifdef DOUBLE
            if (DABS(F(II,JJ)).gt.TEN) then
#else
            if (ABS(F(II,JJ)).gt.TEN) then
#endif

              reduced=.false.  !obliga a repetir el ciclo de reducción
        
              !Algoritmo de Jacobi (reducción de los elementos no diagonales):
              V1=F(JJ,JJ)
              V2=F(II,JJ)
              V3=F(II,II)
              U=.5D0*(F(JJ,JJ)-F(II,II))
  
#ifdef DOUBLE
              if (DABS(U).lt.1.d-10) then
#else
              if (ABS(U).lt.1.d-10) then
#endif  
                    OMEGA=-1.0D0
              else
#ifdef DOUBLE
                OMEGA=-F(II,JJ)/DSQRT(F(II,JJ)*F(II,JJ)+U*U)
#else
                OMEGA=-F(II,JJ)/SQRT(F(II,JJ)*F(II,JJ)+U*U)
#endif  
                Z=1.D0
                if(U.LT.0.D0) Z=-Z
                  OMEGA=OMEGA*Z
                endif
    
#ifdef DOUBLE
              SINT=OMEGA/DSQRT(2.D0*(1.D0+DSQRT(1.D0-OMEGA*OMEGA)))
              COST=DSQRT(1.D0-SINT*SINT)
#else
              SINT=OMEGA/SQRT(2.D0*(1.D0+SQRT(1.D0-OMEGA*OMEGA)))
              COST=SQRT(1.D0-SINT*SINT)
#endif 
    
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
    IF(TE.LT.RHO) converged=.true.
    
    enddo !do while externo (convergencia global)
    
    
!     do I=2,N
!       II=I-1
!       do J=1,II
!         F(J,I)=F(I,J)
!       enddo
!     enddo

! print*, "Original"
! do i=1,N
!     print'(1000F8.3)', FF(i,1:N)
! enddo

! print*, "Diagonal"
! do i=1,N
!     print'(1000F8.3)', F(i,1:N)
! enddo

    d(1:N) = (/ (F(i,i), i=1,N ) /)
    
! print*, "Autovalores"
! print'(1000F8.3)', d(1:N)

    !write(6,*)' MATRIZ DIAGONAL'
      !call WriteMatrix(N,N,F,6)
        !write(6,*)' MATRIZ DE VECTORES PROPIOS (NORMALIZADOS)'
    !  call WriteMatrix(N,N,V,6)
    
    !WRITE(6,*)' *** FIN DE LA DIAGONALIZACION ***'
    
    RETURN

END SUBROUTINE JACOBI_SIM


!C-----------------------------------------------------------------------
SUBROUTINE RMS(TE,F,N)
!C-----------------------------------------------------------------------
#ifdef DOUBLE
    implicit double precision (A-H,O-Z)
#else
    implicit real (A-H,O-Z)
#endif
    implicit integer (i-n)
    
    integer,intent(in)::N
#ifdef DOUBLE
    double precision,dimension(:,:),intent(in)::F
    double precision,intent(out)::TE
#else
    real,dimension(:,:),intent(in)::F
    real,intent(out)::TE
#endif
    
    TE=0.0D0
    do I=2,N
      K=I-1
      do J=1,K
        TE=TE+2.D0*F(I,J)*F(I,J)
      enddo
    enddo
    
#ifdef DOUBLE
    TE=DSQRT(TE)
#else
    TE=SQRT(TE)
#endif
    
    RETURN

END SUBROUTINE RMS

!C-----------------------------------------------------------------------
        SUBROUTINE ORDEN_DIAG(N,D,V)
!C-----------------------------------------------------------------------
        !Subroutina para ordenar los elemntos diagonales aplicar las mismas permutaciones a las columnas de la matriz de vecotores propios
        
                implicit none

                INTEGER, intent(in):: N
#ifdef DOUBLE
                double precision, DIMENSION(:,:),intent(inout) :: D,V

                double precision :: Amin, aux
                double precision,DIMENSION(1:N) :: aux2
#else
                real, DIMENSION(:,:),intent(inout) :: D,V

                real :: Amin, aux
                real,DIMENSION(1:N) :: aux2
#endif

                        
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
        END SUBROUTINE ORDEN_DIAG

!C-----------------------------------------------------------------------
        SUBROUTINE REV_ORDEN_DIAG(N,D,V)
!C-----------------------------------------------------------------------
        !Subroutina para ordenar los elemntos diagonales aplicar las mismas permutaciones a las columnas de la matriz de vecotores propios
        !Reverse order: from Larger to Lower
        
                implicit none

                INTEGER, intent(in):: N

#ifdef DOUBLE
                double precision, DIMENSION(:,:),intent(inout) :: D,V

                double precision :: Amax, aux
                double precision,DIMENSION(1:N) :: aux2
#else
                real, DIMENSION(:,:),intent(inout) :: D,V

                real :: Amax, aux
                real,DIMENSION(1:N) :: aux2
#endif

                        
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
        END SUBROUTINE REV_ORDEN_DIAG

!C-----------------------------------------------------------------------  
#ifdef DOUBLE
        double precision function PRODES(A,B) RESULT(C)
#else
        real function PRODES(A,B) RESULT(C)
#endif
!C-----------------------------------------------------------------------  
        !Función producto escalar (A · B = C), A, B --vectores  //  C--escalar
        !
                implicit none

#ifdef DOUBLE
                double precision, DIMENSION(:),intent(in) :: A,B
#else
                real, DIMENSION(:),intent(in) :: A,B
#endif
                        
                integer :: i,N

                if ( size(A) /= size(B) ) then
                    print*, "Error in PRODES subroutine: size(A) /= size(B)"
                    stop
                endif
 
                N=size(A)

                C=0.d0
                do  i=1,N
                  C = C + A(i)*B(i)
                enddo

                RETURN
        END function PRODES

  
!C-----------------------------------------------------------------------        
        SUBROUTINE GRAMSTH_COL(NA,Am,Vm)
        !Subroutina de ortogonaización de matrices con el método Gram-Schmidt

!       Descripción
!       Ortoganalización de la matriz Am[IN] de dimension NA * NA[IN]. La matriz ortogonalizada es Vm[OUT]

                implicit double precision (a-h,o-z)
                implicit integer (i-n)

#ifdef DOUBLE
                double precision,dimension(:,:),intent(in)::Am
                double precision,dimension(:,:),intent(out)::Vm
#else
                real,dimension(:,:),intent(in)::Am
                real,dimension(:,:),intent(out)::Vm
#endif

                integer,intent(in)::NA

                integer::i,j,k

                !Primer vector (igual al original)
                Vm(1:NA,1)=Am(1:NA,1)

                !Ciclos de ortogonalizacion de los siguientes vectores
                do i=2,NA
                  !Implementación de la fórmula: 
                  !" v_i = a_i – SUMA_k=1-->i-1{ [(a_i · v_k) / (v_k · v_k)] v_k } ", donde v y a son vectores 
                  Vm(1:NA,i)=Am(1:NA,i)
                  do k=1,i-1
                    rnumerador=PRODES(Am(1:NA,i),Vm(1:NA,k))
                    rdenominador=PRODES(Vm(1:NA,k),Vm(1:NA,k))
                    Vm(1:NA,i)=Vm(1:NA,i)-rnumerador/rdenominador*Vm(1:NA,k)
                  enddo
                enddo
                
                RETURN

        END SUBROUTINE GRAMSTH_COL



        subroutine eigrs(A,N,JOBN,D,Z,IZ,WK,IER)
            !
            ! This subroutine is a f90 translation of the eigrs IMSL subroutine.
            !
            ! Original Documentation follows:
!C Short description:
!C                                  SPECIFICATIONS FOR ARGUMENTS         EIRS0880
!C      INTEGER            N,JOBN,IZ,IER                                  EIRS0890
!C      DOUBLE PRECISION   A(1),D(1),WK(1),Z(IZ,1)                        EIRS0900
!C                                  SPECIFICATIONS FOR LOCAL VARIABLES   EIRS0910
!C      INTEGER            IJOB,IR,JR,IJ,JI,NP1                           EIRS0920
!C      INTEGER            JER,NA,ND,IIZ,IBEG,IL,KK,LK,I,J,K,L            EIRS0930
!C      DOUBLE PRECISION   ANORM,ASUM,PI,SUMZ,SUMR,AN,S,TEN,RDELP,ZERO,   EIRS0940
!C     1                   ONE,THOUS                                      EIRS0950
!C The full IMSL documentation for EIGRS follows.
!
!C   IMSL ROUTINE NAME   - EIGRS
!C
!C-----------------------------------------------------------------------
!C
!C   COMPUTER            - IBM/DOUBLE
!C
!C   LATEST REVISION     - JUNE 1, 1980
!C
!C   PURPOSE             - EIGENVALUES AND (OPTIONALLY) EIGENVECTORS OF
!C                           A REAL SYMMETRIC MATRIX
!C
!C   USAGE               - CALL EIGRS (A,N,JOBN,D,Z,IZ,WK,IER)
!C
!C   ARGUMENTS    A      - INPUT REAL SYMMETRIC MATRIX OF ORDER N,
!C                           WHOSE EIGENVALUES AND EIGENVECTORS
!C                           ARE TO BE COMPUTED. INPUT A IS
!C                           DESTROYED IF IJOB IS EQUAL TO 0 OR 1.
!C                N      - INPUT ORDER OF THE MATRIX A.
!C                JOBN   - INPUT OPTION PARAMETER.  IF JOBN.GE.10
!C                         A IS ASSUMED TO BE IN FULL STORAGE MODE
!C                         (IN THIS CASE, A MUST BE DIMENSIONED EXACTLY
!C                         N BY N IN THE CALLING PROGRAM).
!C                         IF JOBN.LT.10 THEN A IS ASSUMED TO BE IN
!C                         SYMMETRIC STORAGE MODE.  DEFINE
!C                         IJOB=MOD(JOBN,10).  THEN WHEN
!C                           IJOB = 0, COMPUTE EIGENVALUES ONLY
!C                           IJOB = 1, COMPUTE EIGENVALUES AND EIGEN-
!C                             VECTORS.
!C                           IJOB = 2, COMPUTE EIGENVALUES, EIGENVECTORS
!C                             AND PERFORMANCE INDEX.
!C                           IJOB = 3, COMPUTE PERFORMANCE INDEX ONLY.
!C                           IF THE PERFORMANCE INDEX IS COMPUTED, IT IS
!C                           RETURNED IN WK(1). THE ROUTINES HAVE
!C                           PERFORMED (WELL, SATISFACTORILY, POORLY) IF
!C                           WK(1) IS (LESS THAN 1, BETWEEN 1 AND 100,
!C                           GREATER THAN 100).
!C                D      - OUTPUT VECTOR OF LENGTH N,
!C                           CONTAINING THE EIGENVALUES OF A.
!C                Z      - OUTPUT N BY N MATRIX CONTAINING
!C                           THE EIGENVECTORS OF A.
!C                           THE EIGENVECTOR IN COLUMN J OF Z CORRES-
!C                           PONDS TO THE EIGENVALUE D(J).
!C                           IF IJOB = 0, Z IS NOT USED.
!C                IZ     - INPUT ROW DIMENSION OF MATRIX Z EXACTLY AS
!C                           SPECIFIED IN THE DIMENSION STATEMENT IN THE
!C                           CALLING PROGRAM.
!C                WK     - WORK AREA, THE LENGTH OF WK DEPENDS
!C                           ON THE VALUE OF IJOB, WHEN
!C                           IJOB = 0, THE LENGTH OF WK IS AT LEAST N.
!C                           IJOB = 1, THE LENGTH OF WK IS AT LEAST N.
!C                           IJOB = 2, THE LENGTH OF WK IS AT LEAST
!C                             N(N+1)/2+N.
!C                           IJOB = 3, THE LENGTH OF WK IS AT LEAST 1.
!C                IER    - ERROR PARAMETER (OUTPUT)
!C                         TERMINAL ERROR
!C                           IER = 128+J, INDICATES THAT EQRT2S FAILED
!C                             TO CONVERGE ON EIGENVALUE J. EIGENVALUES
!C                             AND EIGENVECTORS 1,...,J-1 HAVE BEEN
!C                             COMPUTED CORRECTLY, BUT THE EIGENVALUES
!C                             ARE UNORDERED. THE PERFORMANCE INDEX
!C                             IS SET TO 1000.0
!C                         WARNING ERROR (WITH FIX)
!C                           IN THE FOLLOWING, IJOB = MOD(JOBN,10).
!C                           IER = 66, INDICATES IJOB IS LESS THAN 0 OR
!C                             IJOB IS GREATER THAN 3. IJOB SET TO 1.
!C                           IER = 67, INDICATES IJOB IS NOT EQUAL TO
!C                             ZERO, AND IZ IS LESS THAN THE ORDER OF
!C                             MATRIX A. IJOB IS SET TO ZERO.
!C
!C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32
!C                       - SINGLE/H36,H48,H60
!C
!C   REQD. IMSL ROUTINES - EHOBKS,EHOUSS,EQRT2S,UERTST,UGETIO
!C
!C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND
!C                           CONVENTIONS IS AVAILABLE IN THE MANUAL
!C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP
!C
!C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.
!C
!C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN
!C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,
!C                           EXPRESSED OR IMPLIED, IS APPLICABLE.


#ifdef DOUBLE
        implicit double precision (A-H,O-Z)  
#else
        implicit real (A-H,O-Z)  
#endif      
        implicit integer (I-N)        

        !Input
        integer, intent(in) :: N, JOBN, IZ
#ifdef DOUBLE
        double precision, dimension(:),intent(inout) :: A
#else
        real, dimension(:),intent(inout) :: A
#endif  
        
        !Output
        integer, intent(out) :: IER
#ifdef DOUBLE
        double precision, dimension(:),intent(out) :: D, WK
        double precision, dimension(:,:),intent(out) :: Z
#else
        real, dimension(:),intent(out) :: D, WK
        real, dimension(:,:),intent(out) :: Z
#endif  

        !Parameters in the subroutine
#ifdef DOUBLE
        double precision :: rdelp, ZERO, ONE, TEN, THOUS
#else
        real :: rdelp, ZERO, ONE, TEN, THOUS
#endif 

        data rdelp/2.775557562d-17/ 
        data ZERO,ONE/0.0D0,1.0D0/,TEN/10.0D0/,THOUS/1000.0d0/

        print'(/,X,A,/)', "Entering eigrs (addapted to F90)"

        ier = 0
        jer = 0
       

        !Convert to symmetric storage mode (only if jobn<10)
        if (JOBN>=10) then
            k = 1        
            ji= N-1
            ij= 1
            do j=1,N
                do i=1,j
                    A(k) = A(ij)
                    ij   = ij + 1
                    k    = k + 1
                enddo
                ij = ij + ji
                ji = ji - 1
            enddo
        endif

        ijob=mod(JOBN,10)
        !IJOB shold be 0<=IJOB<=3
        if (ijob<0 .or. ijob>3) then
            !IJOB out of range (reset to 1)
            IER  = 66
            IJOB = 1
        elseif (ijob /= 0 .and. iz < N) then
            !IZ IS LESS THAN N.EIGENVECTORS CAN NOT BE COMPUTED
            ! IJOB SET TO ZERO
            IER  = 67
            IJOB = 0
        endif
                
        if (IJOB /= 3) then
            !Then do the work (3 only checks performance)
            NA = (N*(N+1))/2

            if (IJOB == 2) then
                !SAVE INPUT A IF IJOB=2 
                do i=1,NA
                    WK(i) = A(i)
                enddo
            endif

            ND=1
            if (IJOB == 2) ND = NA + 1

            !Reduce to symmetric triagonal form
            call EHOUSS(A,N,D,WK(1:ND),WK(1:ND))

            IIZ = 1
            if ( IJOB /= 0) then
                !Compute eigenvectors too
                IIZ = IZ
                !Identity matrix
                do i=1,N
                    do j=1,N
                        Z(1:N,1:N) = ZERO
                     enddo
                     Z(i,i) = ONE
                enddo
            endif

           !==========================================
            !Compute eigenvalues and eigenvectors
            call EQRT2S(D,WK(1:ND),N,Z,IIZ,JER)
           !==========================================

            !Exit of only insterested eigenvalues
            if (IJOB == 0) then
                if (JER /= 0) IER=JER 
                return
            endif  
            !Back tranform eigenvectors if JER<128 ¿?
            if (JER < 128) call EHOBKS(A,N,1,N,Z,IZ)
            if (IJOB >= 1 ) then
                if (JER /= 0) IER=JER 
                return
            endif 

            do i=1,NA
                A(i) = WK(i)
            enddo
            WK(1)=THOUS
            if (JER /= 0) then
                if (JER /= 0) IER=JER 
                return
            endif 
    
        endif

        !-------------------------------------------------------
        !From now on, just things related to Performance Index
        
        !Compute 1-Norm of A
        ANORM = ZERO
        IBEG  = 1
        do i=1,N
            ASUM = ZERO
            IL   = IBEG
            KK   = 1
            do L=1,N
#ifdef DOUBLE
                ASUM =  ASUM + DABS(A(IL))
#else
                ASUM =  ASUM + ABS(A(IL))
#endif
                if (L>=i) kk=L
                IL = IL + kk
            enddo
#ifdef DOUBLE
            ANORM = dmax1(ANORM,ASUM)
#else
            ANORM = max1(ANORM,ASUM)
#endif
            IBEG = IBEG + 1
        enddo
        if (ANORM == ZERO) ANORM = ONE

        !Actually compute performance index
        PI = ZERO                                                  
        do I=1,N                                                      
            IBEG = 1                                                       
            S = ZERO                                                       
            SUMZ = ZERO                                                    
                DO L=1,N                                                    
                    LK = IBEG                                                   
                    KK = 1                                                      
#ifdef DOUBLE
                    SUMZ = SUMZ+DABS(Z(L,I))
#else
                    SUMZ = SUMZ+ABS(Z(L,I))
#endif                                    
                    SUMR = -D(I)*Z(L,I)                                         
                    DO K=1,N                                                 
                        SUMR = SUMR+A(LK)*Z(K,I)                                 
                        IF (K.GE.L) KK = K                                       
                        LK = LK+KK                                               
                    enddo                                                    
#ifdef DOUBLE
                    S = S+DABS(SUMR)   
#else
                    S = S+ABS(SUMR)
#endif   
                    IBEG = IBEG+L                                               
                enddo
#ifdef DOUBLE
                IF (SUMZ /= ZERO) PI = DMAX1(PI,S/SUMZ)   
#else
                IF (SUMZ /= ZERO) PI = MAX1(PI,S/SUMZ)   
#endif                                       
        enddo                                                          
        AN = N                                                            
        PI = PI/(ANORM*TEN*AN*RDELP)                                      
        WK(1) = PI                                                        
        IF (JOBN.LT.10) then
            if (JER /= 0) IER=JER 
            return
        endif                                        

        !CONVERT BACK TO FULL STORAGE MODE    
        NP1 = N+1                                                         
        IJ = (N-1)*NP1 + 2                                                
        K = (N*(NP1))/2                                                   
        DO JR=1,N                                                     
            J = NP1-JR                                                     
            DO IR=1,J                                                  
                IJ = IJ-1                                                   
                A(IJ) = A(K)                                                
                K = K-1   
            enddo
            IJ = IJ-JR                                                     
        enddo                                                          
        JI = 0                                                            
        K = N-1                                                           
        DO I=1,N                                                      
            IJ = I-N                                                       
            DO J=1,I                                                   
                IJ = IJ+N                                                   
                JI = JI+1                                                   
                A(IJ) = A(JI)
            enddo
            JI = JI + K                                                    
            K = K-1
        enddo 

        if (JER /= 0) IER=JER                                                       

        return

    end subroutine eigrs


    SUBROUTINE EHOUSS(A,N,D,E,E2)            
!C   IMSL ROUTINE NAME   - EHOUSS                                        
!C                                                                       
!C-----------------------------------------------------------------------
!C                                                                       
!C   COMPUTER            - VAX/DOUBLE                                    
!C                                                                       
!C   LATEST REVISION     - OCTOBER 1, 1980                               
!C                                                                       
!C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRS     
!C                                                                       
!C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
!C                       - SINGLE/H36,H48,H60                            
!C                                                                       
!C   REQD. IMSL ROUTINES - NONE REQUIRED                                 
!C                                                                       
!C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
!C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
!C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
!C                                                                       
!C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.       
!C                                                                       
!C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
!C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
!C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
!C                                                                       
!C-----------------------------------------------------------------------
!C        

#ifdef DOUBLE
        implicit double precision (A-H,O-Z)  
#else
        implicit real (A-H,O-Z)  
#endif      
        implicit integer (I-N) 
                                                               
        !Input
        integer, intent(in) :: N
#ifdef DOUBLE
        double precision, dimension(:), intent(inout) :: A
        double precision, dimension(:), intent(out) :: D, E, E2
#else
        real, dimension(:), intent(inout) :: A
        real, dimension(:), intent(out) :: D, E, E2
#endif   
            
#ifdef DOUBLE
        double precision :: ZERO,H,SCALE,ONE,SCALE1,F,G,HH  
#else
        real :: ZERO,H,SCALE,ONE,SCALE1,F,G,HH  
#endif                                  
 
        DATA               ZERO/0.0D0/,ONE/1.0D0/                         

!C      FIRST EXECUTABLE STATEMENT         
        NP1 = N+1                                                         
        NN = (N*NP1)/2-1                                                  
        NBEG = NN+1-N                                                     
        DO II = 1,N                                                    
            I = NP1-II                                                     
            L = I-1                                                        
            H = ZERO                                                       
            SCALE = ZERO                                                   
!             IF (L .LT. 1) GO TO 10   
            IF (L > 1) then                                         
!C              SCALE ROW (ALGOL TOL THEN NOT NEEDED)
                NK = NN                                                        
                DO K = 1,L   
#ifdef DOUBLE
                    SCALE = SCALE+DABS(A(NK))  
#else
                    SCALE = SCALE+ABS(A(NK))  
#endif                 
                    NK = NK-1                                                   
                enddo
            endif                                                   
            IF (SCALE == ZERO .and. L > 1) then                                  
                E(I) = ZERO                                                    
                E2(I) = ZERO  
                !And exit                                                 
                D(I) = A(NBEG+I)                                               
                A(NBEG+I) = H*SCALE*SCALE                                      
                NBEG = NBEG-I+1                                                
                NN = NN-I
                cycle                             
            endif
                                                       
            NK = NN                                                        
            SCALE1 = ONE/SCALE                                             
            DO K = 1,L                                                  
                A(NK) = A(NK)*SCALE1                                        
                H = H+A(NK)*A(NK)                                           
                NK = NK-1
            enddo
            E2(I) = SCALE*SCALE*H                                          
            F = A(NN)   
#ifdef DOUBLE
            G = -DSIGN(DSQRT(H),F) 
#else
            G = -SIGN(SQRT(H),F) 
#endif                                          
            E(I) = SCALE*G                                                 
            H = H-F*G                                                      
            A(NN) = F-G                                                    
            IF (L .EQ. 1) then
                DO K = 1,L                                                  
                    A(NBEG+K) = SCALE*A(NBEG+K)                                 
                enddo                                                         
                D(I) = A(NBEG+I)                                               
                A(NBEG+I) = H*SCALE*SCALE                                      
                NBEG = NBEG-I+1                                                
                NN = NN-I   
                cycle 
            endif
                                       
            F = ZERO                                                       
            JK1 = 1                                                        
            DO J = 1,L                                                  
                G = ZERO                                                    
                IK = NBEG+1                                                 
                JK = JK1                                                    
!C              FORM ELEMENT OF A*U                  
                DO K = 1,J                                               
                    G = G+A(JK)*A(IK)                                        
                    JK = JK+1                                                
                    IK = IK+1 
                enddo                                               
                JP1 = J+1                                                   
!                 IF (L .LT. JP1) GO TO 35                                    
                IF (L > JP1) then       
                    JK = JK+J-1                                                 
                    DO K = JP1,L                                             
                        G = G+A(JK)*A(IK)                                        
                        JK = JK+K                                                
                        IK = IK+1
                    enddo            
                endif                                                  
!C              FORM ELEMENT OF P                    
                E(J) = G/H                                                  
                F = F+E(J)*A(NBEG+J)                                        
                JK1 = JK1+J
            enddo                                                 
                                                 
            HH = F/(H+H)                                                   
!C          FORM REDUCED A                       
            JK = 1                                                         
            DO J = 1,L                                                  
                F = A(NBEG+J)                                               
                G = E(J)-HH*F                                               
                E(J) = G                                                    
                DO K = 1,J                                               
                    A(JK) = A(JK)-F*E(K)-G*A(NBEG+K)                         
                    JK = JK+1    
                enddo
            enddo                                            
                                                   
            do K = 1,L                                                  
                A(NBEG+K) = SCALE*A(NBEG+K)
            enddo                                                                                    
            D(I) = A(NBEG+I)                                               
            A(NBEG+I) = H*SCALE*SCALE                                      
            NBEG = NBEG-I+1                                                
            NN = NN-I

        enddo                                                      
                                                         
        RETURN                                                            
    end subroutine EHOUSS


    SUBROUTINE EQRT2S (D,E,N,Z,IZ,IER)                                                       
!C   IMSL ROUTINE NAME   - EQRT2S                                        
!C                                                                       
!C-----------------------------------------------------------------------
!C                                                                       
!C   COMPUTER            - VAX/DOUBLE                                    
!C                                                                       
!C   LATEST REVISION     - OCTOBER 15,1980                               
!C                                                                       
!C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGCH     
!C                                                                       
!C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
!C                       - SINGLE/H36,H48,H60                            
!C                                                                       
!C   REQD. IMSL ROUTINES - UERTST,UGETIO                                 
!C                                                                       
!C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
!C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
!C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
!C                                                                       
!C   COPYRIGHT           - 1978 BY IMSL, INC. ALL RIGHTS RESERVED.       
!C                                                                       
!C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
!C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
!C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
!C                                                                       
!C-----------------------------------------------------------------------
!C        

#ifdef DOUBLE
        implicit double precision (A-H,O-Z)  
#else
        implicit real (A-H,O-Z)  
#endif      
        implicit integer (I-N) 

        !Input
        integer,intent(in) :: N, IZ
#ifdef DOUBLE
        double precision, dimension(:), intent(inout) :: E
#else
        real, dimension(:), intent(inout) :: E
#endif

        !Output
        integer,intent(out) :: IER
#ifdef DOUBLE
        double precision, dimension(:), intent(out) :: D
        double precision, dimension(:,:), intent(out) :: Z
#else
        real, dimension(:), intent(out) :: D
        real, dimension(:,:), intent(out) :: Z
#endif  
                                                               
                           
#ifdef DOUBLE
        double precision :: B,C,F,G,H,P,R,S,RDELP,ONE,ZERO     
#else
        real :: B,C,F,G,H,P,R,S,RDELP,ONE,ZERO     
#endif  
         
        DATA               RDELP/2.775557562D-17/                         
        DATA               ZERO,ONE/0.0D0,1.0D0/    
                      
!C      MOVE THE LAST N-1 ELEMENTS           
!C      OF E INTO THE FIRST N-1 LOCATIONS    
!C      FIRST EXECUTABLE STATEMENT           
        IER  = 0                                                          
        IF (N .EQ. 1) return                                          
        DO I=2,N                                                       
            E(I-1) = E(I)
        enddo                                                   
        E(N) = ZERO                                                       
        B = ZERO                                                          
        F = ZERO                                                          
        DO  L=1,N                                                     
            J = 0
#ifdef DOUBLE
            H = RDELP*(DABS(D(L))+DABS(E(L)))     
#else
            H = RDELP*(ABS(D(L))+ABS(E(L)))     
#endif                                                                         
            IF (B.LT.H) B = H                                              
!C          LOOK FOR SMALL SUB-DIAGONAL ELEMENT  
            DO M=L,N                                                   
                K=M
#ifdef DOUBLE
                IF (DABS(E(K)) <= B) exit
#else
                IF (ABS(E(K)) <= B) exit
#endif               
            enddo
            M = K                                                          
            IF (M.EQ.L) then
                D(L) = D(L) + F
                cycle
            endif       
            !==============           
            ! ITERATIONS
            !==============
#ifdef DOUBLE
            do while (DABS(E(L)) .GT. B)
#else
            do while (ABS(E(L)) .GT. B)
#endif

                IF (J == 30) then
                    !Max number of iterations exeeded
                    IER = 128 + L
                    return
                endif
                                     
                J = J+1                                                        
                L1 = L+1                                                       
                G = D(L)                                                       
                P = (D(L1)-G)/(E(L)+E(L))
#ifdef DOUBLE
                R = DSQRT(P*P+ONE)                                             
                D(L) = E(L)/(P+DSIGN(R,P))   
#else
                R = SQRT(P*P+ONE)                                             
                D(L) = E(L)/(P+SIGN(R,P))   
#endif                                        
                                  
                H = G-D(L)                                                     
                DO I = L1,N                                                 
                    D(I) = D(I)-H  
                enddo                                             
                F = F+H                                                        
!C                                  QL TRANSFORMATION                    
                P = D(M)                                                       
                C = ONE                                                        
                S = ZERO                                                       
                MM1 = M-1                                                      
                MM1PL = MM1+L                                                  
                IF (L < MM1) then
                    DO II=L,MM1
                        I = MM1PL-II                                                
                        G = C*E(I)                                                  
                        H = C*P
#ifdef DOUBLE
                        IF (DABS(P) > DABS(E(I))) then                         
                            C = E(I)/P                                                  
                            R = DSQRT(C*C+ONE)   
#else
                        IF (ABS(P) > ABS(E(I))) then                           
                            C = E(I)/P                                                  
                            R = SQRT(C*C+ONE)   
#endif                                        
                            E(I+1) = S*P*R                                              
                            S = C/R                                                     
                            C = ONE/R
                        else                                                                                                   
                            C = P/E(I)                                                  
#ifdef DOUBLE
                            R = DSQRT(C*C+ONE)     
#else
                            R = SQRT(C*C+ONE)     
#endif                                      
                            E(I+1) = S*E(I)*R                                           
                            S = ONE/R                                                   
                            C = C*S
                        endif                                                     
                        P = C*D(I)-S*G                                              
                        D(I+1) = H+S*(C*G+S*D(I))                                   
                        IF (IZ .LT. N) cycle                                     
!C                                  FORM VECTOR                          
                        DO K=1,N                                                 
                            H = Z(K,I+1)                                             
                            Z(K,I+1) = S*Z(K,I)+C*H                                  
                            Z(K,I) = C*Z(K,I)-S*H
                        enddo
                    enddo
                endif                               
                E(L) = S*P                                                     
                D(L) = C*P 

            enddo !ITERATIONS TO CONVERGE EVERY ELEMENT (MAX=30, CRITERIO: E(L)<B)                                                
                            
            D(L) = D(L) + F

        enddo !Over all (J) elements

! TODO: Add a option to order if needed                                                    
! !C      ORDER EIGENVALUES AND EIGENVECTORS   
!         DO I=1,N                                                     
!             K = I                                                          
!             P = D(I)                                                       
!             IP1 = I+1                                                      
! !             IF (IP1.GT.N) GO TO 70  
!             IF (IP1 < N) then                                       
!                 DO J=IP1,N                                                 
!                     IF (D(J) .GE. P) cycle                                   
!                     K = J                                                       
!                     P = D(J)              
!                 enddo
!             endif
!             IF (K.EQ.I) cycle                                           
!             D(K) = D(I)                                                    
!             D(I) = P                                                       
!             IF (IZ .LT. N) cycle                                        
!             DO J = 1,N                                                  
!                 P = Z(J,I)                                                  
!                 Z(J,I) = Z(J,K)                                             
!                 Z(J,K) = P 
!             enddo                                                 
!         enddo                                                    
                                                       
        return                                                         

    end subroutine EQRT2S  



    SUBROUTINE EHOBKS (A,N,M1,M2,Z,IZ)                                                            
!C   IMSL ROUTINE NAME   - EHOBKS                                        
!C                                                                       
!C-----------------------------------------------------------------------
!C                                                                       
!C   COMPUTER            - VAX/DOUBLE                                    
!C                                                                       
!C   LATEST REVISION     - OCTOBER 1, 1980                               
!C                                                                       
!C   PURPOSE             - NUCLEUS CALLED ONLY BY IMSL ROUTINE EIGRS     
!C                                                                       
!C   PRECISION/HARDWARE  - SINGLE AND DOUBLE/H32                         
!C                       - SINGLE/H36,H48,H60                            
!C                                                                       
!C   REQD. IMSL ROUTINES - NONE REQUIRED                                 
!C                                                                       
!C   NOTATION            - INFORMATION ON SPECIAL NOTATION AND           
!C                           CONVENTIONS IS AVAILABLE IN THE MANUAL      
!C                           INTRODUCTION OR THROUGH IMSL ROUTINE UHELP  
!C                                                                       
!C   COPYRIGHT           - 1980 BY IMSL, INC. ALL RIGHTS RESERVED.       
!C                                                                       
!C   WARRANTY            - IMSL WARRANTS ONLY THAT IMSL TESTING HAS BEEN 
!C                           APPLIED TO THIS CODE. NO OTHER WARRANTY,    
!C                           EXPRESSED OR IMPLIED, IS APPLICABLE.        
!C                                                                       
!C-----------------------------------------------------------------------
!C                                                                       
          
#ifdef DOUBLE
        implicit double precision (A-H,O-Z)  
#else
        implicit real (A-H,O-Z)  
#endif      
        implicit integer (I-N) 

        !Input
        integer, intent(in) :: N,IZ,M1,M2
#ifdef DOUBLE
        double precision, dimension(:), intent(in) :: A
        !Output
        double precision, dimension(:,:),intent(out) :: Z
#else
        real, dimension(:), intent(in) :: A
        !Output
        real, dimension(:,:),intent(out) :: Z
#endif 

#ifdef DOUBLE
        double precision :: H,S 
#else
        real :: H,S 
#endif                                              
!C                                  SPECIFICATIONS FOR LOCAL VARIABLES   
!C                                  FIRST EXECUTABLE STATEMENT           
        IF (N==1) return
                                              
        DO I=2,N
                                                       
            L = I-1                                                        
            IA = (I*L)/2                                                   
            H = A(IA+I)                                                    
            IF (H.EQ.0.D0) cycle                                        
!C                                  DERIVES EIGENVECTORS M1 TO M2 OF     
!C                                  THE ORIGINAL MATRIX FROM EIGENVECTORS
!C                                  M1 TO M2 OF THE SYMMETRIC            
!C                                  TRIDIAGONAL MATRIX                   
            DO J=M1,M2                                                  
                S = 0.0D0                                                   
                DO K=1,L                                                  
                    S = S+A(IA+K)*Z(K,J)  
                enddo                                                                                     
                S = S/H                                                     
                DO K=1,L                                                 
                    Z(K,J) = Z(K,J)-S*A(IA+K)    
                enddo  
            enddo                          

        enddo                                                                                                       

        RETURN    
                                                        
    end subroutine EHOBKS  

!=======================================================

    !Simple routines
    subroutine sort_vec(V,N)

        ! From min to max

        implicit none

#ifdef DOUBLE
        double precision,dimension(:),intent(inout) :: V
        double precision :: aux
#else
        real,dimension(:),intent(inout) :: V
        real :: aux
#endif  
        integer,intent(in) :: N
        integer :: i,j

        do i=1,N-1
            do j=i+1,N
                if (V(j)<V(i)) then
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

#ifdef DOUBLE
        double precision,dimension(:),intent(inout) :: V
        double precision :: aux
#else
        real,dimension(:),intent(inout) :: V
        real :: aux
#endif  
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


    subroutine sort_ivec(V,N)

        !Integer version of sort_vec subroutine

        implicit none

        integer,dimension(:),intent(inout) :: V
        integer :: aux
        integer,intent(in) :: N
        integer :: i,j

        do i=1,N-1
            do j=i+1,N
                if (V(j)<V(i)) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                endif
            enddo
        enddo

        return

    end subroutine sort_ivec

    subroutine sort_ivec_max(V,N)

        !Integer version of sort_vec_max subroutine

        implicit none

        integer,dimension(:),intent(inout) :: V
        integer :: aux
        integer,intent(in) :: N
        integer :: i,j

        do i=1,N-1
            do j=i+1,N
                if (V(j)>V(i)) then
                    aux=V(i)
                    V(i) = V(j)
                    V(j) = aux
                endif
            enddo
        enddo

        return

    end subroutine sort_ivec_max                                
                                   


end module MatrixMod
