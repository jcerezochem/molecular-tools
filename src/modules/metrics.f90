module metrics

    !This module is a  structure_types-free version of atomic_geom
    ! Uses only double precision

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS 
    !==============================================================
    ! MODULE with functions to meassure molecular geometric parameters  
    ! (bonds angles...) using str_atom as input
    !     function calc_dist(atom1,atom2) result(b) 
    !     function calc_angle(atom1,atom2,atom3) result(a) 
    !     function calc_dihed(atom1,atom2,atom3,atom4) result(d) 
    !     subroutine rotation_3D(vx,vy,vz,tx,ty,tz,Theta) 
    ! NOTE: this code is re-collected from ancient (beguiner) trials
    ! thus more error prone than usual
    ! Changelog:
    !   Angles are returned in radians. Conversion to deg in main
    !    (somentimes are intermendiate values, which should be in rad)
    !==============================================================

    !Comon stuff
    use constants

    implicit none

    contains

    function calc_dist(X1,Y1,Z1, &
                       X2,Y2,Z2) result(b)

        real(8),intent(in) :: X1,Y1,Z1,&
                              X2,Y2,Z2
        real(8) :: b
        
        !Local
        real(8),dimension(3)::r1,r2, raux

        !Interface: (X,Y,Z) to vectors
        r1(1:3)=(/X1,Y1,Z1/)
        r2(1:3)=(/X2,Y2,Z2/)

        raux=r1-r2
        b=dsqrt(dot_product(raux,raux))

        return

    end function calc_dist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function calc_angle(X1,Y1,Z1, &
                        X2,Y2,Z2, &
                        X3,Y3,Z3) result(a)

        real(8),intent(in) :: X1,Y1,Z1,&
                              X2,Y2,Z2,&
                              X3,Y3,Z3
        real(8) :: a

        !Local
        double precision b1,b2
        double precision,dimension(1:3)::r1,r2,r3
        double precision,dimension(1:3)::raux1,raux2
        integer :: i

        !Interface: (X,Y,Z) to vectors
        r1(1:3)=(/X1,Y1,Z1/)
        r2(1:3)=(/X2,Y2,Z2/)
        r3(1:3)=(/X3,Y3,Z3/)

        do i=1,3
            raux1(i)=r1(i)-r2(i)
            raux2(i)=r3(i)-r2(i)
        enddo
        raux1 = r1-r2
        raux2 = r3-r2

        b1=dsqrt(dot_product(raux1,raux1))
        b2=dsqrt(dot_product(raux2,raux2))
        a=dacos(dot_product(raux1,raux2)/b1/b2)

        return

    end function calc_angle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function calc_cosangle(X1,Y1,Z1, &
                           X2,Y2,Z2, &
                           X3,Y3,Z3) result(a)

        real(8),intent(in) :: X1,Y1,Z1,&
                              X2,Y2,Z2,&
                              X3,Y3,Z3
        real(8) :: a

        !Local
        double precision b1,b2
        double precision,dimension(1:3)::r1,r2,r3
        double precision,dimension(1:3)::raux1,raux2
        integer :: i

        !Interface: (X,Y,Z) to vectors
        r1(1:3)=(/X1,Y1,Z1/)
        r2(1:3)=(/X2,Y2,Z2/)
        r3(1:3)=(/X3,Y3,Z3/)

        do i=1,3
            raux1(i)=r1(i)-r2(i)
            raux2(i)=r3(i)-r2(i)
        enddo

        b1 = dsqrt(raux1(1)**2+raux1(2)**2+raux1(3)**2)
        b2 = dsqrt(raux2(1)**2+raux2(2)**2+raux2(3)**2)
        a  = raux2(1)*raux1(1)+raux2(2)*raux1(2)+raux2(3)*raux1(3)


        return

    end function calc_cosangle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function calc_dihed(X1,Y1,Z1, &
                        X2,Y2,Z2, &
                        X3,Y3,Z3, &
                        X4,Y4,Z4) result(dh)

        !----------------------------------------------------------------
        ! Esta función coloca el dihedro en una orientación estandard
        ! y mide el dihdro. Se obtienen valores con signo de acuerod
        ! con el criterio establecido.
        ! NOTE:
        ! The SR works with individual x(), y(), z() vectors, not atomic
        ! 3D vectors. Thus the SR is first interfaced to inpunt 3D atomic
        ! vectors as default in this module.
        !----------------------------------------------------------------

        real(8),intent(in) :: X1,Y1,Z1,&
                              X2,Y2,Z2,&
                              X3,Y3,Z3,&
                              X4,Y4,Z4
        real(8) :: dh

        !Local
        double precision,dimension(1:4)::x,y,z
        double precision::dx,dy,dz,z_aux,y_aux,theta
        integer::i

        !Interface: (X,Y,Z) to vectors (in this case x,y,z)
        x(1:4)=(/X1,X2,X3,X4/)
        y(1:4)=(/Y1,Y2,Y3,Y4/)
        z(1:4)=(/Z1,Z2,Z3,Z4/)


        !1) Center the system at Atom3
        dx=x(3)
        dy=y(3)
        dz=z(3)
        do i=1,4
            x(i)=x(i)-dx
            y(i)=y(i)-dy
            z(i)=z(i)-dz
        enddo


        !2) Lie Atom2--Atom3 on the z-axis.
        !   a. Rotation around X
        theta=datan(y(2)/z(2))
        z_aux=y(2)*dsin(theta)+z(2)*dcos(theta)
        if (z_aux>0) theta=theta+pi ! Note: 2 is behind 3
        do i=1,4
            z_aux=z(i)
            z(i)=y(i)*dsin(theta)+z_aux*dcos(theta)
            y(i)=y(i)*dcos(theta)-z_aux*dsin(theta)
        enddo
        !   b. Rotation around Y
        theta=datan(x(2)/z(2))
        z_aux=x(2)*dsin(theta)+z(2)*dcos(theta)
        if (z_aux>0) theta=theta+pi
        do i=1,4
            z_aux=z(i)
            z(i)=x(i)*dsin(theta)+z_aux*dcos(theta)
            x(i)=x(i)*dcos(theta)-z_aux*dsin(theta)
        enddo

        !3) Put Atom1 on the YZ plane (y>0) 
        !   a. Rotation around Z
        theta=datan(x(1)/y(1))
        y_aux=x(1)*dsin(theta)+y(1)*dcos(theta)
        if (y_aux<0) theta=theta+pi ! Note: 1 has positive y
        do i=1,4
            y_aux=y(i)
            y(i)=x(i)*dsin(theta)+y_aux*dcos(theta)
            x(i)=x(i)*dcos(theta)-y_aux*dsin(theta)
        enddo
        theta=datan(x(1)/y(1))
        y_aux=x(1)*dsin(theta)+y(1)*dcos(theta)
        if (y_aux<0) theta=theta+pi ! Note: 1 has positive y
        do i=1,4
            y_aux=y(i)
            y(i)=x(i)*dsin(theta)+y_aux*dcos(theta)
            x(i)=x(i)*dcos(theta)-y_aux*dsin(theta)
        enddo

        dh=datan(x(4)/y(4))
        ! Comprobamos que de los dos ángulos posibles (en +/- 180)
        ! obtenemos el que nos interesa comparando con el criterio
        ! según el cual signo(dh)=signo(x(4))
        if (dh>0 .and. x(4)<0) dh=dh-pi
        if (dh<0 .and. x(4)>0) dh=dh+pi

        ! WARNING!!
        ! Sale con el signo cambiado! (si no se sabe por qué, esto no vale..)
        dh=-dh
        
        return

     end function calc_dihed

    function calc_improper(X1,Y1,Z1, &
                           X2,Y2,Z2, &
                           X3,Y3,Z3, &
                           X4,Y4,Z4) result(im)

        !----------------------------------------------------------------
        ! Ańgulo formado por el vector 4-1 con el plano 2-3-4
        !        2
        !       /
        !   1--4
        !       \
        !       3
        !----------------------------------------------------------------

        use constants

        implicit none

        real(8),intent(in) :: X1,Y1,Z1,&
                              X2,Y2,Z2,&
                              X3,Y3,Z3,&
                              X4,Y4,Z4
        real(8) :: im

        !Local
        real(8) :: v1,v1x,v1y,v1z,&
                      v2x,v2y,v2z,&
                   vn,vnx,vny,vnz
        integer::i

        !Vectors 4-2 y 4-3
        v1x = X2 - X4
        v1y = Y2 - Y4
        v1z = Z2 - Z4
        v2x = X3 - X4
        v2y = Y3 - Y4
        v2z = Z3 - Z4
        !Vector normal al plano
        vnx = v1y*v2z - v1z*v2y
        vny = v1z*v2x - v1x*v2z
        vnz = v1x*v2y - v1y*v2x
        !Se normaliza
        vn = dsqrt(vnx**2+vny**2+vnz**2)
        vnx=vnx/vn
        vny=vny/vn
        vnz=vnz/vn

        !Se calcula y normaliza 4-1 (en v1)
        v1x = X1 - X4
        v1y = Y1 - Y4
        v1z = Z1 - Z4
        v1 = dsqrt(v1x**2+v1y**2+v1z**2)
        v1x=v1x/v1
        v1y=v1y/v1
        v1z=v1z/v1

        !Se calcula el producto escalar de vn y v1
        im = v1x*vnx + v1y*vny + v1z*vnz
        im = acos(im)

        !El ángulo es el complementario del anterior
        im = PI/2.d0 - im 

        return

    end function calc_improper

    function calc_dihed_new(X1,Y1,Z1, &
                            X2,Y2,Z2, &
                            X3,Y3,Z3, &
                            X4,Y4,Z4) result(dh)

        !----------------------------------------------------------------
        ! Esta función coloca el dihedro en una orientación estandard
        ! con el enlace 2-3 sobre el eje X con 3 en la parte negativa
        ! el signo se obtiene con el criterio de que la rotación de la
        ! proyección del eje 1-2 al eje 3-4 en el sentido de las agujas
        ! del reloj indica el diedro positivo
        !----------------------------------------------------------------

        use constants

        implicit none

        real(8),intent(in) :: X1,Y1,Z1,&
                              X2,Y2,Z2,&
                              X3,Y3,Z3,&
                              X4,Y4,Z4
        real(8) :: dh

        !Local parameters
        real(8),parameter :: ZERO = 1.d-10
        !Local variables
        real(8) :: aux_X1,aux_Y1,aux_Z1,&
                   aux_X2,aux_Y2,aux_Z2,&
                   aux_X3,aux_Y3,aux_Z3,&
                   aux_X4,aux_Y4,aux_Z4

        real(8) :: v1,v2,v3,Theta,sinTheta,baux
        integer::i

        !New frame 
        !1. Traslate to have atom 3 at origin
        aux_X1 = X1 - X3
        aux_Y1 = Y1 - Y3
        aux_Z1 = Z1 - Z3
        aux_X2 = X2 - X3
        aux_Y2 = Y2 - Y3
        aux_Z2 = Z2 - Z3
        aux_X4 = X4 - X3
        aux_Y4 = Y4 - Y3
        aux_Z4 = Z4 - Z3
        !2. Rotate to place atom 2 on -x axis
        ! Find rotation axis and angle with vector product of position(atom i_a) X (-1,0,0)
        v1 =  0.d0
        v2 = -aux_Z2
        v3 =  aux_Y2
        baux = dsqrt(aux_X2**2+aux_Y2**2+aux_Z2**2)
        sinTheta = dsqrt(aux_Y2**2+aux_Z2**2)/baux
        Theta =  dasin(sinTheta)
        
        !Fix rotation (no need to chech y and z??? see below
        if (aux_X2 > 0) Theta = PI-Theta
        
        !Only rotate if needed
        if (Theta == PI) then
            aux_X1 = -aux_X1
            aux_X1 = -aux_X1
            aux_X1 = -aux_X1
            aux_X4 = -aux_X4
            aux_X4 = -aux_X4
            aux_X4 = -aux_X4
            aux_X2 = -aux_X2
            aux_X2 = -aux_X2
            aux_X2 = -aux_X2
        else if (Theta /= 0.d0) then
            call rotation_3D(aux_X1,&
                             aux_Y1,&
                             aux_Z1,&
                             v1,v2,v3,Theta)
            call rotation_3D(aux_X4,&
                             aux_Y4,&
                             aux_Z4,&
                             v1,v2,v3,Theta)
            call rotation_3D(aux_X2,&
                             aux_Y2,&
                             aux_Z2,&
                             v1,v2,v3,Theta)
        endif
        if (aux_X2>0.d0) then
            print*, "ERROR IN calc_dihed_new - 1"
            stop
        endif
        
        !Normalize projections
        baux = dsqrt(aux_Y1**2+aux_Z1**2)
        
        if (baux == 0.d0) then
            print*, "Singularity when computing dihedral"
            stop
        endif
        aux_Y1 = aux_Y1/baux
        aux_Z1 = aux_Z1/baux
        baux = dsqrt(aux_Y4**2+aux_Z4**2)
        
        aux_Y4 = aux_Y4/baux
        aux_Z4 = aux_Z4/baux
        !Scalar product
        baux = aux_Y1*aux_Y4+aux_Z1*aux_Z4
        if (dabs(dabs(baux) - 1.d0) < ZERO) baux=sign(1.d0,baux)
        dh = dacos(baux)
        
        
        !Check sign
        !Rotation around X to place 1-2 projection on +y
        ! Find rotation axis and angle with vector product of position(atom i_a) Y (0,+1,0) => (1,0,0)
        baux = dsqrt(aux_Y1**2+aux_Z1**2)
        if (baux == 0.d0) then
            print*, "Singularity when computing dihedral"
            stop
        endif
        sinTheta = dabs(aux_Z1)/baux
        Theta =  dasin(sinTheta)
        !Fix rotation (it's clockwise if not inverted)
        if (aux_Y1 < 0) Theta = PI-Theta
        if (aux_Z1 > 0) Theta = -Theta
        call rotation_3D(aux_X1,&
                         aux_Y1,&
                         aux_Z1,&
                         1.d0,0.d0,0.d0,Theta)
        call rotation_3D(aux_X4,&
                         aux_Y4,&
                         aux_Z4,&
                         1.d0,0.d0,0.d0,Theta)
        if (aux_Y1<0.d0) then
            print*, "ERROR IN calc_dihed_new - 2"
            stop
        endif
        
        !Ajustamos el signo (esta vez es sistmático)
        dh = sign(dabs(dh),aux_Z4)
        
        return

    end function calc_dihed_new


    subroutine rotation_3D(vx,vy,vz,tx,ty,tz,Theta,RR) 
    
            !Description:
            ! Subroutine to rotate the vector (vx,vy,vz) around the axis
            ! defined by (tx,ty,tz) an angle Theta (rad).
            real(8), intent(inout) :: vx,vy,vz
            real(8), intent(in)    :: tx,ty,tz, Theta
            real(8), dimension(1:3,1:3), intent(out), optional :: RR 
    
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
    
            if (present(RR)) RR = R
    
            return
    
    end subroutine rotation_3D


end module metrics
