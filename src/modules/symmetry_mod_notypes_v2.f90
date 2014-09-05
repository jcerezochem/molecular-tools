module symmetry_mod_notypes

    !Uses geom as a 3Nat vector with (x1,y1,z1,x2,y2,z2 ..., xN,yN,zN)
    !Only accepts symmetry guessing (sym group is not an input/output)
    !
    !History
    ! Version 2: allows multiple symmetries to be detected (Ci and C2),
    !            thus allowing several equivalent atoms. See symm_atoms_nt

    implicit none

    contains

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine symm_atoms_nt(geom,natoms,isym)

        !This routine tries to guess the symmetry (only Ci and C2 for the moment)
        ! Multiple symmetries can be simultaneously detected, thus identifying 
        ! several equivalent atoms. The output (isym) only reports the atom with
        ! the lower serial number.

        use alerts
        use MatrixMod

        real(8),parameter :: THRS=1.d-1 !loose

        !System variables
        real(8),dimension(:),intent(in)  :: geom
        integer,intent(in)               :: natoms
        integer,dimension(:),intent(out) :: isym
        integer,dimension(200,10)        :: isym_v2

        double precision :: X_COG, Y_COG, Z_COG, dist, x_sym, y_sym, z_sym
        logical :: located

        !local
        integer :: i,j, ii,jj, nsym=0

        !====================
        !Center of Geometry  (?? shouldt be COM ??)
        !====================
        X_COG=0.d0
        Y_COG=0.d0
        Z_COG=0.d0
        do i=1,natoms
            ii = i*3-2
            X_COG=X_COG+geom(ii)
            Y_COG=Y_COG+geom(ii+1)
            Z_COG=Z_COG+geom(ii+2)
        enddo
        X_COG=X_COG/dfloat(natoms)
        Y_COG=Y_COG/dfloat(natoms)
        Z_COG=Z_COG/dfloat(natoms)

        ! Axes...


!         !==================
!         !Use symmetry 
!         !==================
!         if ( molec%PG(1:2) == "CI" ) then
!         
!             do i=1,molec%natoms
!                 ii = i*3-2
!                 x_sym=2.d0*X_COG-geom(ii)
!                 y_sym=2.d0*Y_COG-geom(ii+1)
!                 z_sym=2.d0*Z_COG-geom(ii+2)
! 
!                 !Search the atom located at the symmetric position (with a given threshold)
!                 do j=1,natoms
!                     dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
! !                     if ( dist < 1d-3 ) then
!                     if ( dist < THRS ) then
!                         isym(i)=j
!                         exit
!                     endif
!                 enddo
! 
!             enddo
! 
!         elseif ( molec%PG(1:3) == "C02" ) then
!         
!             do i=1,molec%natoms
!                 !Only works if C2 is the Z axis (default for standard orientation)
!                 x_sym=2.d0*X_COG-molec%atom(i)%x
!                 y_sym=2.d0*Y_COG-molec%atom(i)%y
!                 z_sym=molec%atom(i)%z
! 
!                 !Search the atom located at the symmetric position (with a given threshold)
!                 do j=1,molec%natoms
!                     dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
! !                     if ( dist < 1d-3 ) then
!                     if ( dist < THRS ) then
!                         isym(i)=j
!                         exit
!                     endif
!                 enddo
! 
!             enddo
! 
!         elseif ( molec%PG(1:2) == "XX" ) then

            ! CI
            !If there is an equivalent atom at -r, is CI symmetry (at least)
            nsym = nsym+1
            do i=1,natoms
                ii = 3*i-2
                located=.false.
                x_sym=2.d0*X_COG-geom(ii)
                y_sym=2.d0*Y_COG-geom(ii+1)
                z_sym=2.d0*Z_COG-geom(ii+2)

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,natoms
                    jj = 3*j-2
                    dist=dsqrt((geom(jj)-x_sym)**2+(geom(jj+1)-y_sym)**2+(geom(jj+2)-z_sym)**2)
                    if ( dist < THRS ) then
                        located=.true.
                        isym(i)=j
                        isym_v2(i,nsym)=j
                        exit
                    endif
                enddo
                if (.not.located) exit

            enddo
            if (located) then
                 print*, ""
                 print*, "CI symmetry element found"
!                  molec%PG="CI"   
!                  return
            else
                 nsym = nsym-1
            endif

            ! C2
            !Only works if C2 is the Z axis (default for standard orientation)
            nsym = nsym+1
            do i=1,natoms
                ii = 3*i-2
                located=.false.
                x_sym=2.d0*X_COG-geom(ii)
                y_sym=2.d0*Y_COG-geom(ii+1)
                z_sym=geom(ii+2)

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,natoms
                    jj = 3*j-2
                    dist=dsqrt((geom(jj)-x_sym)**2+(geom(jj+1)-y_sym)**2+(geom(jj+2)-z_sym)**2)
                    if ( dist < THRS ) then
                        located=.true.
                        isym(i)=j
                        isym_v2(i,nsym)=j
                        exit
                    endif
                enddo
                if (.not.located) exit

            enddo
            if (located) then
!                  nsym = nsym+1
                 print*, ""
                 print*, "C02(z) symmetry element found"
!                  molec%PG="C02"   
!                  return
            else
                 nsym = nsym-1
            endif

            !Trying C2 along Y axis 
            nsym = nsym+1
            do i=1,natoms
                ii = 3*i-2
                located=.false.
                x_sym=2.d0*X_COG-geom(ii)
                y_sym=geom(ii+1)
                z_sym=2.d0*Z_COG-geom(ii+2)

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,natoms
                    jj = 3*j-2
                    dist=dsqrt((geom(jj)-x_sym)**2+(geom(jj+1)-y_sym)**2+(geom(jj+2)-z_sym)**2)
                    if ( dist < THRS ) then
                        located=.true.
                        isym(i)=j
                        isym_v2(i,nsym)=j
                        exit
                    endif
                enddo
                if (.not.located) exit

            enddo
            if (located) then
!                  nsym = nsym+1
                 print*, ""
                 print*, "C02(y) symmetry element found"
!                  molec%PG="C02"   
!                  return
            else
                 nsym = nsym-1
            endif

            !Trying C2 along Y axis 
            nsym = nsym+1
            do i=1,natoms
                ii = 3*i-2
                located=.false.
                x_sym=geom(ii)
                y_sym=2.d0*Y_COG-geom(ii+1)
                z_sym=2.d0*Z_COG-geom(ii+2)

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,natoms
                    jj = 3*j-2
                    dist=dsqrt((geom(jj)-x_sym)**2+(geom(jj+1)-y_sym)**2+(geom(jj+2)-z_sym)**2)
                    if ( dist < THRS ) then
                        located=.true.
                        isym(i)=j
                        isym_v2(i,nsym)=j
                        exit
                    endif
                enddo
                if (.not.located) exit

            enddo
            if (located) then
!                  nsym = nsym+1
                 print*, ""
                 print*, "C02(x) symmetry element found"
!                  molec%PG="C02"   
!                  return
            else
                 nsym = nsym-1
            endif
 
            !No sym found
            if (nsym == 0) then
                call alert_msg("note","Cannot determine symmetry. Set to C1 (no symmetric atoms)")
                do i=1,natoms
                    isym(i) = i
                enddo
!                 molec%PG="C1"
            else
                do i=1,natoms
                    call sort_vec_int(isym_v2(i,1:nsym),nsym)
                    isym(i) = isym_v2(i,1)
                enddo
            endif
            print*, ""

!             do i=1,natoms
!                 print*, i, isym(i)
!             enddo

            !Higher symmetry groups shuld also be tested... (TODO)
            

!         else
!             call alert_msg("note","This kind of symmetry in not supported (yet): "//molec%PG)
! 
!         endif


        return

    end subroutine symm_atoms_nt



end module symmetry_mod_notypes