module geom_meter

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    ! MODULE with functions to meassure molecular geometric parameters  
    ! (bonds angles...) using 3D vector positions (not independet x,y,z).
    ! So, it is not compatible with deprecated "structure_types".
    !     function dist(r1,r2) result(b) 
    !     function angle(r1,r2,r3) result(a) 
    !     function dihed(r1,r2,r3,r4) result(d) 
    ! NOTE: this code is re-collected from ancient (beguiner) trials
    ! thus more error prone than usual
    !==============================================================

    use constants

!    implicit none
!#ifdef DOUBLE
!    double precision,parameter:: pi=3.14159265358979323846d0
!#else
!    real,parameter:: pi=3.14159265358979323846
!#endif

    contains

    function dist(r1,r2) result(b)

#ifdef DOUBLE
        double precision b
        double precision,intent(in),dimension(3)::r1,r2
        double precision,dimension(3)::raux
#else
        real b
        real,intent(in),dimension(3)::r1,r2
        real,dimension(3)::raux
#endif

        raux=r1-r2
#ifdef DOUBLE
        b=dsqrt(dot_product(raux,raux))
#else
        b=sqrt(dot_product(raux,raux))
#endif

        return

    end function dist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function angle(r1,r2,r3) result(a)

#ifdef DOUBLE
        double precision a,b1,b2
        double precision,intent(in),dimension(3)::r1,r2,r3
        double precision,dimension(3)::raux1,raux2

#else
        real a,b1,b2
        real,intent(in),dimension(3)::r1,r2,r3
        real,dimension(3)::raux1,raux2
#endif

        raux1=r1-r2
        raux2=r3-r2

#ifdef DOUBLE
        b1=dsqrt(dot_product(raux1,raux1))
        b2=dsqrt(dot_product(raux2,raux2))
        a=dacos(dot_product(raux1,raux2)/b1/b2)*360.d0/2.d0/pi
#else
        b1=sqrt(dot_product(raux1,raux1))
        b2=sqrt(dot_product(raux2,raux2))
        a=acos(dot_product(raux1,raux2)/b1/b2)*360./2./pi
#endif

        return

    end function angle

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dihed(r1,r2,r3,r4) result(d)

       !This function does not work properly. Use new_dihed instead

#ifdef DOUBLE
        double precision d,b1,b2
        double precision,intent(in),dimension(3)::r1,r2,r3,r4
        double precision,dimension(3)::raux1,raux2,raux3,raux4
#else
        real d,b1,b2
        real,intent(in),dimension(3)::r1,r2,r3,r4
        real,dimension(3)::raux1,raux2,raux3,raux4
#endif

        raux1=r1-r2
        raux2=r3-r2

        raux3(1)=raux1(2)*raux2(3)-raux1(3)*raux2(2)
        raux3(2)=raux1(3)*raux2(1)-raux1(1)*raux2(3)
        raux3(3)=raux1(1)*raux2(2)-raux1(2)*raux2(1)

        raux1=r2-r3
        raux2=r4-r3

        raux4(1)=raux1(2)*raux2(3)-raux1(3)*raux2(2)
        raux4(2)=raux1(3)*raux2(1)-raux1(1)*raux2(3)
        raux4(3)=raux1(1)*raux2(2)-raux1(2)*raux2(1)

#ifdef DOUBLE
        b1=dsqrt(dot_product(raux3,raux3))
        b2=dsqrt(dot_product(raux4,raux4))
        d=dacos(dot_product(raux3,raux4)/b1/b2)*360.d0/2.d0/pi
#else
        b1=sqrt(dot_product(raux3,raux3))
        b2=sqrt(dot_product(raux4,raux4))
        d=acos(dot_product(raux3,raux4)/b1/b2)*360./2./pi
#endif

        return

    end function dihed

    function new_dihed(r1,r2,r3,r4) result(dh)

        !----------------------------------------------------------------
        ! Esta función coloca el dihedro en una orientación estandard
        ! y mide el dihdro. Se obtienen valores con signo de acuerod
        ! con el criterio establecido.
        ! NOTE:
        ! The SR works with individual x(), y(), z() vectors, not atomic
        ! 3D vectors. Thus the SR is first interfaced to inpunt 3D atomic
        ! vectors as default in this module.
        !----------------------------------------------------------------

#ifdef DOUBLE
        double precision::dh
        double precision,intent(in),dimension(1:3)::r1,r2,r3,r4
        double precision,dimension(1:4)::x,y,z
        double precision::dx,dy,dz,z_aux,y_aux,theta
#else
        real::dh
        real,intent(in),dimension(1:3)::r1,r2,r3,r4
        real,dimension(1:4)::x,y,z
        real::dx,dy,dz,z_aux,y_aux,theta
#endif
        integer::i

        !Interface 3D vectors with x,y,z arrays
        x(1)=r1(1)
        y(1)=r1(2)
        z(1)=r1(3)
        x(2)=r2(1)
        y(2)=r2(2)
        z(2)=r2(3)
        x(3)=r3(1)
        y(3)=r3(2)
        z(3)=r3(3)
        x(4)=r4(1)
        y(4)=r4(2)
        z(4)=r4(3)

! !Bug: If any point equals 0,0,0 this is because it was not read (dangerous)!
! !     If so, we force NaN to be returned
!         do i=1,4
!             if (x(i)==0. .and. y(i)==0. .and. z(i)==0.) then
!                 dh=sin(i/0.)
!                 return
!             endif
!         enddo

        !1) Center the str_system at A3
        dx=x(3)
        dy=y(3)
        dz=z(3)
        do i=1,4
            x(i)=x(i)-dx
            y(i)=y(i)-dy
            z(i)=z(i)-dz
        enddo


        !2) Lie A2--A3 on the z-axis.
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

        !3) Put A1 on the YZ plane (y>0) 
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

        dh=datan(x(4)/y(4))*360.d0/2.d0/pi
        ! Comprobamos que de los dos ángulos posibles (en +/- 180)
        ! obtenemos el que nos interesa comparando con el criterio
        ! según el cual signo(dh)=signo(x(4))
        if (dh>0 .and. x(4)<0) dh=dh-180.d0
        if (dh<0 .and. x(4)>0) dh=dh+180.d0
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

        !3) Put A1 on the YZ plane (y>0) 
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

        dh=atan(x(4)/y(4))*360./2./pi
        ! Comprobamos que de los dos ángulos posibles (en +/- 180)
        ! obtenemos el que nos interesa comparando con el criterio
        ! según el cual signo(dh)=signo(x(4))
        if (dh>0 .and. x(4)<0) dh=dh-180.
        if (dh<0 .and. x(4)>0) dh=dh+180.
#endif

     end function new_dihed

!     function new_dihed(r1,r2,r3,r4) result(dh)
!! TRANSLATION OF THE WHOLE SR TO 3D VECTOR FASHION
! 
!         !----------------------------------------------------------------
!         ! Esta función coloca el dihedro en una orientación estandard
!         ! y mide el dihdro. Se obtienen valores con signo de acuerdo
!         ! con el criterio establecido.
!         !----------------------------------------------------------------
! 
!         implicit none
!         real::dh
!         real,intent(in),dimension(1:4)::x0,y0,z0
!         real,dimension(1:4)::x,y,z
!         real,dimension(1:3) :: d
!         real::dx,dy,dz,z_aux,y_aux,theta
!         real,parameter:: pi=3.14159265358979323846
! 
!         integer::i
! 
!         x=x0
!         y=y0
!         z=z0
! 
! ! !Bug: If any point equals 0,0,0 this is because it was not read (dangerous)!
! ! !     If so, we force NaN to be returned
! !         do i=1,4
! !             if (x(i)==0. .and. y(i)==0. .and. z(i)==0.) then
! !                 dh=sin(i/0.)
! !                 return
! !             endif
! !         enddo
! 
!         !1) Center the str_system at A3
!         d(1)=r3(1)
!         d(2)=r3(2)
!         d(3)=r3(3)
!         do i=1,3
!             r1(i)=r1(i)-d(i)
!             r2(i)=r2(i)-d(i)
!             r3(i)=r3(i)-d(i)
!             r4(i)=r4(i)-d(i)
!         enddo
! 
! 
!         !2) Lie A2--A3 on the z-axis.
!         !   a. Rotation around X
!         theta=atan(r2(2)/r2(3))
!         z_aux=r2(2)*sin(theta)+r2(3)*cos(theta)
!         if (z_aux>0) theta=theta+pi ! Note: 2 is behind 3
!             z_aux=z(i)
!         z_aux=r1(3)
!         r1(3) = r1(2)*sin(theta)+z_aux*cos(theta)
!             z(i)=y(i)*sin(theta)+z_aux*cos(theta)
!             y(i)=y(i)*cos(theta)-z_aux*sin(theta)
! 
!         !   b. Rotation around Y
!         theta=atan(r2(1)/r2(3))
!         z_aux=r2(1)*sin(theta)+r2(3)*cos(theta)
!         if (z_aux>0) theta=theta+pi
!         do i=1,4
!             z_aux=z(i)
!             z(i)=x(i)*sin(theta)+z_aux*cos(theta)
!             x(i)=x(i)*cos(theta)-z_aux*sin(theta)
!         enddo
! 
!         !3) Put A1 on the YZ plane (y>0) 
!         !   a. Rotation around Z
!         theta=atan(r1(1)/r1(2))
!         y_aux=r1(1)*sin(theta)+r1(2)*cos(theta)
!         if (y_aux<0) theta=theta+pi ! Note: 1 has positive y
!         do i=1,4
!             y_aux=y(i)
!             y(i)=x(i)*sin(theta)+y_aux*cos(theta)
!             x(i)=x(i)*cos(theta)-y_aux*sin(theta)
!         enddo
! 
!         dh=atan(r4(1)/r4(2))*360./2./pi
!         ! Comprobamos que de los dos ángulos posibles (en +/- 180)
!         ! obtenemos el que nos interesa comparando con el criterio
!         ! según el cual signo(dh)=signo(r4(1))
!         if (dh>0 .and. r4(1)<0) dh=dh-180.
!         if (dh<0 .and. r4(1)>0) dh=dh+180.
! 
!      end function new_dihed


end module geom_meter
