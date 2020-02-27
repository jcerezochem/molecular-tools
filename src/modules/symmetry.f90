module symmetry

    use verbosity
    use line_preprocess
    use alerts

    implicit none

    character(len=6),save :: sym_thr="normal"  !"loose", "normal",  "tight", "vtight"

    contains

    subroutine symm_atoms(molec,isym,Osym,rotate,nsym_ops)

        !This routine also try to guess the symmetry (only Ci and C2 along Cartesian axis, for the moment)
        ! (in compatible with the old routine calls, so this is taken as default)
        ! Multiple symmetries can be simultaneously detected, thus identifying 
        ! several equivalent atoms. The output (isym) only reports the atom with
        ! the lower serial number.

        use structure_types
        use matrix, only:sort_vec_int, diagonalize_full
        use matrix_print
        use molecular_structure


        integer,parameter :: NDIM = 600

        !System variables
        type(str_resmol),intent(inout)   :: molec
        integer,dimension(:),intent(out) :: isym
        integer,dimension(NDIM,10)        :: isym_v2
        integer,dimension(4,NDIM,NDIM),intent(out),optional :: Osym
        logical,intent(in),optional :: rotate ! default is true
        integer,intent(out),optional :: nsym_ops

        double precision :: X_COM, Y_COM, Z_COM, dist, x_sym, y_sym, z_sym, Mass
        logical :: located
        ! threshold to identify symmetric atoms (now can be changed)
        real(8) :: THRS

        !local
        integer :: i,j, ii,jj, nsym=0
!         type(str_resmol) :: molec_local
        real(8),dimension(:),allocatable :: X0,Y0,Z0
        real(8),dimension(3,3) :: Rot, MI
        real(8),dimension(3) :: R
        integer :: Nat
        character(len=6) :: sym_thr_local
        logical          :: rotate_local

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with Angstrong 
        current_units=molec%units
        call set_geom_units(molec,"Angs")
        !----------------------------------------------

        rotate_local = .true.
        if (present(rotate)) then
            rotate_local = rotate
        endif

        ! Set THRS according to the sym_thr variable
        sym_thr_local = adjustl(sym_thr)
        call set_word_upper_case(sym_thr_local)
        select case (sym_thr_local)
            case("LOOSE" ) 
             THRS=1.d-1
            case("NORMAL") 
             THRS=1.d-2
            case("TIGHT" ) 
             THRS=1.d-3
            case("VTIGHT") 
             THRS=1.d-5
            case default
             call alert_msg("fatal","Thr for symmetry detection cannot be set")
        end select

        !Save original coordinates
        Nat=molec%natoms
        allocate(X0(1:Nat),Y0(1:Nat),Z0(1:Nat))
        X0=molec%atom(1:Nat)%x
        Y0=molec%atom(1:Nat)%y
        Z0=molec%atom(1:Nat)%z

        ! FIRST PLACE THE MOLECULE ON THE COM
        if (.not.rotate_local) then
            print'(/,X,A,/)', "Using input geometry instead of rotating to principal axis of intertia"
        else
            !====================
            !Center of Masses
            !====================
            call get_com(molec)
            ! Translate
            do i=1,molec%natoms
                molec%atom(i)%x = molec%atom(i)%x - molec%comX
                molec%atom(i)%y = molec%atom(i)%y - molec%comY
                molec%atom(i)%z = molec%atom(i)%z - molec%comZ
            enddo
            
            !====================
            !Axes of intertia
            !====================
            !Get moment of intertia
            call inertia(molec,MI)
            !Diagonalize to get the rotation to the principal axes
            call diagonalize_full(MI(1:3,1:3),3,Rot(1:3,1:3),R(1:3),"lapack")
            ! And rotate
            Rot=transpose(Rot)
            if (verbose>2) &
             call MAT0(6,Rot,3,3,"Rotation matrix (symm_atoms)")
            call rotate_molec(molec,Rot)
        endif


        !==================
        !Use symmetry (note: only the main symmetry element is used!)
        !==================
        if ( molec%PG(1:2) == "CI" ) then
        
            do i=1,molec%natoms
                x_sym=-molec%atom(i)%x
                y_sym=-molec%atom(i)%y
                z_sym=-molec%atom(i)%z

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
                x_sym=-molec%atom(i)%x
                y_sym=-molec%atom(i)%y
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
                x_sym=-molec%atom(i)%x
                y_sym=-molec%atom(i)%y
                z_sym=-molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
                    if ( dist < THRS ) then
                        located=.true.
                        isym(i)=j
                        isym_v2(i,nsym)=j
                        if (present(Osym)) then
                            ! Form the corresponding element of the matrix representation
                            ii = 3*i-2
                            jj = 3*j-2
                            Osym(nsym,ii+0,jj+0) = -1
                            Osym(nsym,ii+1,jj+1) = -1
                            Osym(nsym,ii+2,jj+2) = -1
                        endif
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
                 if (verbose>2.and.present(Osym)) &
                  call MAT0(6,dfloat(Osym(nsym,1:3*Nat,1:3*Nat)),3*Nat,3*Nat,"Matrix Representation (Cartesian basis)")
            else
                 nsym = nsym-1
            endif

            ! C2
            !Only works if C2 is the Z axis (default for standard orientation)
            nsym = nsym+1
            do i=1,molec%natoms
                located=.false.
                x_sym=-molec%atom(i)%x
                y_sym=-molec%atom(i)%y
                z_sym=molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
                    if ( dist < THRS ) then
                        located=.true.
                        isym(i)=j
                        isym_v2(i,nsym)=j
                        if (present(Osym)) then
                            ! Form the corresponding element of the matrix representation
                            ii = 3*i-2
                            jj = 3*j-2
                            Osym(nsym,ii+0,jj+0) = -1
                            Osym(nsym,ii+1,jj+1) = -1
                            Osym(nsym,ii+2,jj+2) =  1
                        endif
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
                 if (verbose>2.and.present(Osym)) &
                  call MAT0(6,dfloat(Osym(nsym,1:3*Nat,1:3*Nat)),3*Nat,3*Nat,"Matrix Representation (Cartesian basis)")
            else
                 nsym = nsym-1
            endif

            ! C2 along Y axis
            nsym = nsym+1
            do i=1,molec%natoms
                located=.false.
                x_sym=-molec%atom(i)%x
                y_sym=molec%atom(i)%y
                z_sym=-molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
                    if ( dist < THRS ) then
                        located=.true.
                        isym(i)=j
                        isym_v2(i,nsym)=j
                        if (present(Osym)) then
                            ! Form the corresponding element of the matrix representation
                            ii = 3*i-2
                            jj = 3*j-2
                            Osym(nsym,ii+0,jj+0) = -1
                            Osym(nsym,ii+1,jj+1) =  1
                            Osym(nsym,ii+2,jj+2) = -1
                        endif
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
                 if (verbose>2.and.present(Osym)) &
                  call MAT0(6,dfloat(Osym(nsym,1:3*Nat,1:3*Nat)),3*Nat,3*Nat,"Matrix Representation (Cartesian basis)")
            else
                 nsym = nsym-1
            endif

            ! C2 along X axis
            nsym = nsym+1
            do i=1,molec%natoms
                located=.false.
                x_sym=molec%atom(i)%x
                y_sym=-molec%atom(i)%y
                z_sym=-molec%atom(i)%z

                !Search the atom located at the symmetric position (with a given threshold)
                do j=1,molec%natoms
                    dist=dsqrt((molec%atom(j)%x-x_sym)**2+(molec%atom(j)%y-y_sym)**2+(molec%atom(j)%z-z_sym)**2)
                    if ( dist < THRS ) then
                        located=.true.
                        isym(i)=j
                        isym_v2(i,nsym)=j
                        if (present(Osym)) then
                            ! Form the corresponding element of the matrix representation
                            ii = 3*i-2
                            jj = 3*j-2
                            Osym(nsym,ii+0,jj+0) =  1
                            Osym(nsym,ii+1,jj+1) = -1
                            Osym(nsym,ii+2,jj+2) = -1
                        endif
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
                 if (verbose>2.and.present(Osym)) &
                  call MAT0(6,dfloat(Osym(nsym,1:3*Nat,1:3*Nat)),3*Nat,3*Nat,"Matrix Representation (Cartesian basis)")
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

        ! Restore original coordinates
        molec%atom(1:Nat)%x=X0
        molec%atom(1:Nat)%y=Y0
        molec%atom(1:Nat)%z=Z0
        deallocate(X0,Y0,Z0)

        if (present(nsym_ops)) nsym_ops = nsym

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! Revert original units
        call set_geom_units(molec,adjustl(current_units))
        !----------------------------------------------

        return

    end subroutine symm_atoms



end module symmetry
