module atomic_geom

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
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
    use structure_types
    use constants

    implicit none
!#ifdef DOUBLE
!    double precision,parameter:: pi=3.14159265358979323846d0
!#else
!    real,parameter:: pi=3.14159265358979323846
!#endif

    contains

    function calc_dist(atom1, atom2) result(b)

        type(str_atom),intent(in) :: atom1, atom2
#ifdef DOUBLE
        double precision b
        double precision,dimension(3)::r1,r2
        double precision,dimension(3)::raux
#else
        real b
        real,dimension(3)::r1,r2
        real,dimension(3)::raux
#endif

        !Interface: atom to vectors
        r1(1:3)=(/atom1%x,atom1%y,atom1%z/)
        r2(1:3)=(/atom2%x,atom2%y,atom2%z/)

        raux=r1-r2
#ifdef DOUBLE
        b=dsqrt(dot_product(raux,raux))
!         b = dsqrt(raux(1)**2+raux(2)**2+raux(3)**2)
#else
        b=sqrt(dot_product(raux,raux))
!         b = sqrt(raux(1)**2+raux(2)**2+raux(3)**2)
#endif

        return

    end function calc_dist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function calc_angle(atom1,atom2,atom3) result(a)

        type(str_atom),intent(in) :: atom1, atom2, atom3
#ifdef DOUBLE
        double precision a,b1,b2
        double precision,dimension(1:3)::r1,r2,r3
        double precision,dimension(1:3)::raux1,raux2

#else
        real a,b1,b2
        real,dimension(1:3)::r1,r2,r3
        real,dimension(1:3)::raux1,raux2
#endif
        integer :: i

        !Interface: atom to vectors
        r1(1:3)=(/atom1%x,atom1%y,atom1%z/)
        r2(1:3)=(/atom2%x,atom2%y,atom2%z/)
        r3(1:3)=(/atom3%x,atom3%y,atom3%z/)

        do i=1,3
            raux1(i)=r1(i)-r2(i)
            raux2(i)=r3(i)-r2(i)
        enddo
        raux1 = r1-r2
        raux2 = r3-r2

#ifdef DOUBLE
!         b1 = dsqrt(raux1(1)**2+raux1(2)**2+raux1(3)**2)
!         b2 = dsqrt(raux2(1)**2+raux2(2)**2+raux2(3)**2)
!         a  = raux2(1)*raux1(1)+raux2(2)*raux1(2)+raux2(3)*raux1(3)
!         a = dacos(a/b1/b2)
        b1=dsqrt(dot_product(raux1,raux1))
        b2=dsqrt(dot_product(raux2,raux2))
        a=dacos(dot_product(raux1,raux2)/b1/b2)
#else
        b1=sqrt(dot_product(raux1,raux1))
        b2=sqrt(dot_product(raux2,raux2))
        a=acos(dot_product(raux1,raux2)/b1/b2)
#endif

        return

    end function calc_angle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function calc_cosangle(atom1,atom2,atom3) result(a)

        type(str_atom),intent(in) :: atom1, atom2, atom3
#ifdef DOUBLE
        double precision a,b1,b2
        double precision,dimension(1:3)::r1,r2,r3
        double precision,dimension(1:3)::raux1,raux2

#else
        real a,b1,b2
        real,dimension(1:3)::r1,r2,r3
        real,dimension(1:3)::raux1,raux2
#endif
        integer :: i

        !Interface: atom to vectors
        r1(1:3)=(/atom1%x,atom1%y,atom1%z/)
        r2(1:3)=(/atom2%x,atom2%y,atom2%z/)
        r3(1:3)=(/atom3%x,atom3%y,atom3%z/)

        do i=1,3
            raux1(i)=r1(i)-r2(i)
            raux2(i)=r3(i)-r2(i)
        enddo

#ifdef DOUBLE
        b1 = dsqrt(raux1(1)**2+raux1(2)**2+raux1(3)**2)
        b2 = dsqrt(raux2(1)**2+raux2(2)**2+raux2(3)**2)
        a  = raux2(1)*raux1(1)+raux2(2)*raux1(2)+raux2(3)*raux1(3)
!         b1=dsqrt(dot_product(raux1,raux1))
!         b2=dsqrt(dot_product(raux2,raux2))
!         a=dacos(dot_product(raux1,raux2)/b1/b2)
#else
        b1=sqrt(dot_product(raux1,raux1))
        b2=sqrt(dot_product(raux2,raux2))
        a=dot_product(raux1,raux2)/b1/b2
#endif

        return

    end function calc_cosangle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function calc_dihed(atom1,atom2,atom3,atom4) result(dh)

        !----------------------------------------------------------------
        ! Esta función coloca el dihedro en una orientación estandard
        ! y mide el dihdro. Se obtienen valores con signo de acuerod
        ! con el criterio establecido.
        ! NOTE:
        ! The SR works with individual x(), y(), z() vectors, not atomic
        ! 3D vectors. Thus the SR is first interfaced to inpunt 3D atomic
        ! vectors as default in this module.
        !----------------------------------------------------------------

        type(str_atom),intent(in) :: atom1, atom2, atom3, atom4
#ifdef DOUBLE
        double precision::dh
        double precision,dimension(1:4)::x,y,z
        double precision::dx,dy,dz,z_aux,y_aux,theta
#else
        real::dh
        real,dimension(1:4)::x,y,z
        real::dx,dy,dz,z_aux,y_aux,theta
#endif
        integer::i

        !Interface: atom to vectors (in this case x,y,z)
        x(1:4)=(/atom1%x,atom2%x,atom3%x,atom4%x/)
        y(1:4)=(/atom1%y,atom2%y,atom3%y,atom4%y/)
        z(1:4)=(/atom1%z,atom2%z,atom3%z,atom4%z/)


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
#ifdef DOUBLE
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
#else
        theta=atan(y(2)/z(2))
        z_aux=y(2)*sin(theta)+z(2)*cos(theta)
        if (z_aux>0) theta=theta+pi ! Note: 2 is behind 3
        do i=1,4
            z_aux=z(i)
            z(i)=y(i)*sin(theta)+z_aux*cos(theta)
            y(i)=y(i)*cos(theta)-z_aux*sin(theta)
        enddo
        !   b. Rotation around Y
        theta=atan(x(2)/z(2))
        z_aux=x(2)*sin(theta)+z(2)*cos(theta)
        if (z_aux>0) theta=theta+pi
        do i=1,4
            z_aux=z(i)
            z(i)=x(i)*sin(theta)+z_aux*cos(theta)
            x(i)=x(i)*cos(theta)-z_aux*sin(theta)
        enddo

        !3) Put Atom1 on the YZ plane (y>0) 
        !   a. Rotation around Z
        theta=atan(x(1)/y(1))
        y_aux=x(1)*sin(theta)+y(1)*cos(theta)
        if (y_aux<0) theta=theta+pi ! Note: 1 has positive y
        do i=1,4
            y_aux=y(i)
            y(i)=x(i)*sin(theta)+y_aux*cos(theta)
            x(i)=x(i)*cos(theta)-y_aux*sin(theta)
        enddo
        theta=atan(x(1)/y(1))
        y_aux=x(1)*sin(theta)+y(1)*cos(theta)
        if (y_aux<0) theta=theta+pi ! Note: 1 has positive y
        do i=1,4
            y_aux=y(i)
            y(i)=x(i)*sin(theta)+y_aux*cos(theta)
            x(i)=x(i)*cos(theta)-y_aux*sin(theta)
        enddo

        dh=atan(x(4)/y(4))
        ! Comprobamos que de los dos ángulos posibles (en +/- 180)
        ! obtenemos el que nos interesa comparando con el criterio
        ! según el cual signo(dh)=signo(x(4))
        if (dh>0 .and. x(4)<0) dh=dh-pi
        if (dh<0 .and. x(4)>0) dh=dh+pi
#endif

        ! Sale con el signo cambiado!
        dh=-dh
        
        return

     end function calc_dihed

    function calc_improper(atom_1,atom_2,atom_3,atom_4) result(im)

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

        type(str_atom),intent(in) :: atom_1, atom_2, atom_3, atom_4
#ifdef DOUBLE
        double precision::im
        real(8) :: v1,v1x,v1y,v1z,&
                      v2x,v2y,v2z,&
                   vn,vnx,vny,vnz
#else
        real::im
        real :: v1,v1x,v1y,v1z,&
                   v2x,v2y,v2z,&
                vn,vnx,vny,vnz
#endif
        integer::i

        !Vectors 4-2 y 4-3
        v1x = atom_2%x - atom_4%x
        v1y = atom_2%y - atom_4%y
        v1z = atom_2%z - atom_4%z
        v2x = atom_3%x - atom_4%x
        v2y = atom_3%y - atom_4%y
        v2z = atom_3%z - atom_4%z
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
        v1x = atom_1%x - atom_4%x
        v1y = atom_1%y - atom_4%y
        v1z = atom_1%z - atom_4%z
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

    function calc_dihed_new(atom_1,atom_2,atom_3,atom_4) result(dh)

        !----------------------------------------------------------------
        ! Esta función coloca el dihedro en una orientación estandard
        ! con el enlace 2-3 sobre el eje X con 3 en la parte negativa
        ! el signo se obtiene con el criterio de que la rotación de la
        ! proyección del eje 1-2 al eje 3-4 en el sentido de las agujas
        ! del reloj indica el diedro positivo
        !----------------------------------------------------------------

        use constants

        implicit none

        type(str_atom),intent(in) :: atom_1, atom_2, atom_3, atom_4
        type(str_atom) :: aux_atom_1,aux_atom_2,aux_atom_3,aux_atom_4
#ifdef DOUBLE
        real(8),parameter :: ZERO = 1.d-8
        double precision::dh
        real(8) :: v1,v2,v3,Theta,sinTheta,baux
#else
        real::dh
        real :: v1,v2,v3,Theta,sinTheta,baux
#endif
        integer::i

    !New frame 
    !1. Traslate to have atom 3 at origin
    aux_atom_1%x = atom_1%x - atom_3%x
    aux_atom_1%y = atom_1%y - atom_3%y
    aux_atom_1%z = atom_1%z - atom_3%z
    aux_atom_2%x = atom_2%x - atom_3%x
    aux_atom_2%y = atom_2%y - atom_3%y
    aux_atom_2%z = atom_2%z - atom_3%z
    aux_atom_4%x = atom_4%x - atom_3%x
    aux_atom_4%y = atom_4%y - atom_3%y
    aux_atom_4%z = atom_4%z - atom_3%z
    !2. Rotate to place atom 2 on -x axis
    ! Find rotation axis and angle with vector product of position(atom i_a) X (-1,0,0)
    v1 =  0.d0
    v2 = -aux_atom_2%z
    v3 =  aux_atom_2%y
#ifdef DOUBLE
    baux = dsqrt(aux_atom_2%x**2+aux_atom_2%y**2+aux_atom_2%z**2)
    sinTheta = dsqrt(aux_atom_2%y**2+aux_atom_2%z**2)/baux
    Theta =  dasin(sinTheta)
#else
    baux = sqrt(aux_atom_2%x**2+aux_atom_2%y**2+aux_atom_2%z**2)
    sinTheta = sqrt(aux_atom_2%y**2+aux_atom_2%z**2)/baux
    Theta =  asin(sinTheta)
#endif
    !Fix rotation (no need to chech y and z??? see below
    if (aux_atom_2%x > 0) Theta = PI-Theta


    !Only rotate if needed
    if (Theta == PI) then
        aux_atom_1%x = -aux_atom_1%x
        aux_atom_1%x = -aux_atom_1%x
        aux_atom_1%x = -aux_atom_1%x
        aux_atom_4%x = -aux_atom_4%x
        aux_atom_4%x = -aux_atom_4%x
        aux_atom_4%x = -aux_atom_4%x
        aux_atom_2%x = -aux_atom_2%x
        aux_atom_2%x = -aux_atom_2%x
        aux_atom_2%x = -aux_atom_2%x
    else if (Theta /= 0.d0) then
        call rotation_3D_(aux_atom_1%x,&
                     aux_atom_1%y,&
                     aux_atom_1%z,&
                     v1,v2,v3,Theta)
        call rotation_3D_(aux_atom_4%x,&
                     aux_atom_4%y,&
                     aux_atom_4%z,&
                     v1,v2,v3,Theta)
        call rotation_3D_(aux_atom_2%x,&
                     aux_atom_2%y,&
                     aux_atom_2%z,&
                     v1,v2,v3,Theta)
    endif
    if (aux_atom_2%x>0.d0) then
        print*, "ERROR IN calc_dihed_new - 1"
        stop
    endif

    !Normalize projections
#ifdef DOUBLE
    baux = dsqrt(aux_atom_1%y**2+aux_atom_1%z**2)
#else
    baux = sqrt(aux_atom_1%y**2+aux_atom_1%z**2)
#endif
    if (baux == 0.d0) then
        print*, "Singularity when computing dihedral"
        stop
    endif
    aux_atom_1%y = aux_atom_1%y/baux
    aux_atom_1%z = aux_atom_1%z/baux
#ifdef DOUBLE
    baux = dsqrt(aux_atom_4%y**2+aux_atom_4%z**2)
#else
    baux = sqrt(aux_atom_4%y**2+aux_atom_4%z**2)
#endif
    aux_atom_4%y = aux_atom_4%y/baux
    aux_atom_4%z = aux_atom_4%z/baux
    !Scalar product
    baux = aux_atom_1%y*aux_atom_4%y+aux_atom_1%z*aux_atom_4%z
#ifdef DOUBLE
    if (dabs(dabs(baux) - 1.d0) < 1.d-10) baux=sign(1.d0,baux)
    dh = dacos(baux)
#else
    if (abs(abs(baux) - 1.) < 1.e-10) baux=sign(1.,baux)
    dh = acos(baux)
#endif

    !Check sign
    !Rotation around X to place 1-2 projection on +y
    ! Find rotation axis and angle with vector product of position(atom i_a) Y (0,+1,0) => (1,0,0)
#ifdef DOUBLE
    baux = dsqrt(aux_atom_1%y**2+aux_atom_1%z**2)
#else
    baux = sqrt(aux_atom_1%y**2+aux_atom_1%z**2)
#endif
    if (baux == 0.d0) then
        print*, "Singularity when computing dihedral"
        stop
    endif
#ifdef DOUBLE
    sinTheta = dabs(aux_atom_1%z)/baux
    Theta =  dasin(sinTheta)
#else
    sinTheta = abs(aux_atom_1%z)/baux
    Theta =  asin(sinTheta)
#endif
    !Fix rotation (it's clockwise if not inverted)
    if (aux_atom_1%y < 0) Theta = PI-Theta
    if (aux_atom_1%z > 0) Theta = -Theta
    call rotation_3D_(aux_atom_1%x,&
                     aux_atom_1%y,&
                     aux_atom_1%z,&
                     1.d0,0.d0,0.d0,Theta)
    call rotation_3D_(aux_atom_4%x,&
                     aux_atom_4%y,&
                     aux_atom_4%z,&
                     1.d0,0.d0,0.d0,Theta)
    if (aux_atom_1%y<0.d0) then
        print*, "ERROR IN calc_dihed_new - 2"
        stop
    endif
 
! print*, ""
!     print*, aux_atom_4%x,aux_atom_4%y,aux_atom_4%z
! print*, "Dihed", dh
#ifdef DOUBLE
    dh = sign(dabs(dh),aux_atom_4%z)
#else
    dh = sign(abs(dh),aux_atom_4%z)
#endif
! print*, "Dihed", dh

    return

    contains

    subroutine rotation_3D_(vx,vy,vz,tx,ty,tz,Theta) 

        !Description:
        ! Subroutine to rotate the vector (vx,vy,vz) around the axis
        ! defined by (tx,ty,tz) an angle Theta (rad).

#ifdef DOUBLE
        real(8), intent(inout) :: vx,vy,vz
        real(8), intent(in)    :: tx,ty,tz, Theta

        !Local
        real(8),dimension(1:3,1:3) :: R
        real(8) :: vx_tmp, vy_tmp, tmod
        real(8) :: ux, uy, uz 
#else
        real, intent(inout) :: vx,vy,vz
        real, intent(in)    :: tx,ty,tz, Theta

        !Local
        real,dimension(1:3,1:3) :: R
        real :: vx_tmp, vy_tmp, tmod
        real :: ux, uy, uz 
#endif

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

    end subroutine rotation_3D_

end function calc_dihed_new

subroutine rotation_3D(vx,vy,vz,tx,ty,tz,Theta,RR) 

        !Description:
        ! Subroutine to rotate the vector (vx,vy,vz) around the axis
        ! defined by (tx,ty,tz) an angle Theta (rad).

#ifdef DOUBLE
        real(8), intent(inout) :: vx,vy,vz
        real(8), intent(in)    :: tx,ty,tz, Theta
        real(8), dimension(1:3,1:3), intent(out), optional :: RR 

        !Local
        real(8),dimension(1:3,1:3) :: R
        real(8) :: vx_tmp, vy_tmp, tmod
        real(8) :: ux, uy, uz 
#else
        real, intent(inout) :: vx,vy,vz
        real, intent(in)    :: tx,ty,tz, Theta
        real, dimension(1:3,1:3), intent(out), optional :: RR 

        !Local
        real,dimension(1:3,1:3) :: R
        real :: vx_tmp, vy_tmp, tmod
        real :: ux, uy, uz 
#endif

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

       RR = R

       return

end subroutine rotation_3D


end module atomic_geom
