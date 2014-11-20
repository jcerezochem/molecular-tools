module symmetry_mod

    !CAUTION: UNCHECKED MODULE (14/11/2012)

    implicit none

    contains

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine symm_atoms_old(molec,isym)

        use structure_types
        use alerts

        !System variables
        type(str_resmol),intent(inout) :: molec
        integer,dimension(:),intent(out) :: isym

        double precision :: X_COG, Y_COG, Z_COG, dist, x_sym, y_sym, z_sym

        !local
        integer :: i,j

        !====================
        !Center of Geometry  (?? shouldt be COM ??)
        !====================
        X_COG=0.d0
        Y_COG=0.d0
        Z_COG=0.d0
        do i=1,molec%natoms
            X_COG=X_COG+molec%atom(i)%x
            Y_COG=Y_COG+molec%atom(i)%y
            Z_COG=Z_COG+molec%atom(i)%z
        enddo
        X_COG=X_COG/dfloat(molec%natoms)
        Y_COG=Y_COG/dfloat(molec%natoms)
        Z_COG=Z_COG/dfloat(molec%natoms)

        ! Axes...


        !==================
        !Use symmetry 
        !==================
        if ( molec%PG(1:2) == "CI" ) then
        
            do i=1,molec%natoms
                x_sym=2.d0*X_COG-molec%atom(i)%x
                y_sym=2.d0*Y_COG-molec%atom(i)%y
                z_sym=2.d0*Z_COG-molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
                    if ( dist < 1d-3 ) then
                        isym(i)=j
                        exit
                    endif
                enddo

            enddo

        else
            call alert_msg("note","This kind of symmetry in not supported (yet): "//molec%PG)

        endif


        return

    end subroutine symm_atoms_old


    subroutine symm_atoms(molec,isym)

        !This routine also try to guess the symmetry (only Ci and C2 along Cartesian axis, for the moment)
        ! (in compatible with the old routine calls, so this is taken as default)
        ! Multiple symmetries can be simultaneously detected, thus identifying 
        ! several equivalent atoms. The output (isym) only reports the atom with
        ! the lower serial number.

        use structure_types
        use alerts
        use MatrixMod, only:sort_vec_int

        real(8),parameter :: THRS=1.d-1 !loose

        !System variables
        type(str_resmol),intent(inout) :: molec
        integer,dimension(:),intent(out) :: isym
        integer,dimension(200,10)        :: isym_v2

        double precision :: X_COG, Y_COG, Z_COG, dist, x_sym, y_sym, z_sym
        logical :: located

        !local
        integer :: i,j, nsym

        !====================
        !Center of Geometry  (?? shouldt be COM ??)
        !====================
        X_COG=0.d0
        Y_COG=0.d0
        Z_COG=0.d0
        do i=1,molec%natoms
            X_COG=X_COG+molec%atom(i)%x
            Y_COG=Y_COG+molec%atom(i)%y
            Z_COG=Z_COG+molec%atom(i)%z
        enddo
        X_COG=X_COG/dfloat(molec%natoms)
        Y_COG=Y_COG/dfloat(molec%natoms)
        Z_COG=Z_COG/dfloat(molec%natoms)

        ! Axes...


        !==================
        !Use symmetry (note: only the main symmetry element is used!)
        !==================
        if ( molec%PG(1:2) == "CI" ) then
        
            do i=1,molec%natoms
                x_sym=2.d0*X_COG-molec%atom(i)%x
                y_sym=2.d0*Y_COG-molec%atom(i)%y
                z_sym=2.d0*Z_COG-molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
!                     if ( dist < 1d-3 ) then
                    if ( dist < THRS ) then
                        isym(i)=j
                        exit
                    endif
                enddo

            enddo

        elseif ( molec%PG(1:3) == "C02" ) then
        
            do i=1,molec%natoms
                !Only works if C2 is the Z axis (default for standard orientation)
                x_sym=2.d0*X_COG-molec%atom(i)%x
                y_sym=2.d0*Y_COG-molec%atom(i)%y
                z_sym=molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
!                     if ( dist < 1d-3 ) then
                    if ( dist < THRS ) then
                        isym(i)=j
                        exit
                    endif
                enddo

            enddo

        elseif ( molec%PG(1:2) == "XX" ) then
        !All possible symmetry elements are used
            nsym = 0

            ! CI
            !If there is an equivalent atom at -r, is CI symmetry (at least)
            nsym = nsym+1
            do i=1,molec%natoms
                located=.false.
                x_sym=2.d0*X_COG-molec%atom(i)%x
                y_sym=2.d0*Y_COG-molec%atom(i)%y
                z_sym=2.d0*Z_COG-molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
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
                 molec%PG="CI"   
!                  return
            else
                 nsym = nsym-1
            endif

            ! C2
            !Only works if C2 is the Z axis (default for standard orientation)
            nsym = nsym+1
            do i=1,molec%natoms
                located=.false.
                x_sym=2.d0*X_COG-molec%atom(i)%x
                y_sym=2.d0*Y_COG-molec%atom(i)%y
                z_sym=molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
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
                 molec%PG="C02"   
!                  return
            else
                 nsym = nsym-1
            endif

            ! C2 along Y axis
            nsym = nsym+1
            do i=1,molec%natoms
                located=.false.
                x_sym=2.d0*X_COG-molec%atom(i)%x
                y_sym=molec%atom(i)%y
                z_sym=2.d0*Z_COG-molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
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
                 print*, "C02(y) symmetry element found"
                 molec%PG="C02"   
!                  return
            else
                 nsym = nsym-1
            endif

            ! C2 along X axis
            nsym = nsym+1
            do i=1,molec%natoms
                located=.false.
                x_sym=molec%atom(i)%x
                y_sym=2.d0*Y_COG-molec%atom(i)%y
                z_sym=2.d0*Z_COG-molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
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
                 print*, "C02(x) symmetry element found"
                 molec%PG="C02"   
!                  return
            else
                 nsym = nsym-1
            endif
 
            !No sym found
            if (nsym == 0) then
                call alert_msg("note","Cannot determine symmetry. Set to C1 (no symmetric atoms)")
                do i=1,molec%natoms
                    isym(i) = i
                enddo
                molec%PG="C1"
            else
                do i=1,molec%natoms
                    call sort_vec_int(isym_v2(i,1:nsym),nsym)
                    isym(i) = isym_v2(i,1)
                enddo
            endif
            print*, ""


            !Higher symmetry groups shuld also be test... (TODO)
            

        else
            call alert_msg("note","This kind of symmetry in not supported (yet): "//molec%PG)

        endif


        return

    end subroutine symm_atoms



end module symmetry_mod