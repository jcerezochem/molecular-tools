module zmat_manage

    use verbosity
    use matrix
    use matrix_print

    contains

    subroutine modredundant(I_RED,molec)
    
        use structure_types
        use alerts
        use constants
        use atomic_geom
        use symmetry
    
        type(str_resmol),intent(inout) :: molec
        integer,intent(in) :: I_RED
        !Local
        character(len=1) :: int_type
        character(len=100) :: line

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with no Units
        !------------------------------------------------
    
        nbonds  = 0
        nangles = 0
        ndihed  = 0
        nimprop = 0    

        do
            read(I_RED,'(A,A)',iostat=ios) int_type, line
            if (ios /=0) exit
            if (int_type == "B") then
                read(line,*) i1, i2
                nbonds = nbonds + 1
                molec%geom%bond(nbonds,1:2) = (/i1,i2/)
            elseif (int_type == "A") then
                read(line,*) i1, i2, i3
                nangles = nangles + 1
                molec%geom%angle(nangles,1:3) = (/i1,i2,i3/)
            elseif (int_type == "D") then
                read(line,*) i1, i2, i3, i4
                ndihed = ndihed + 1
                molec%geom%dihed(ndihed,1:4) = (/i1,i2,i3,i4/)
            elseif (int_type == "I") then
                read(line,*) i1, i2, i3, i4
                nimprop = nimprop + 1
                molec%geom%improp(nimprop,1:4) = (/i1,i2,i3,i4/)
            elseif (int_type == "!") then
                !this is a comment
                cycle
            elseif (int_type == "") then
            !Blank line is the end
                exit
            else
                call alert_msg("fatal","Unknows internal type: "//int_type)
            endif
        enddo
    
        molec%geom%nbonds  = nbonds
        molec%geom%nangles = nangles
        molec%geom%ndihed  = ndihed
        molec%geom%nimprop = nimprop    

        return
    
    end subroutine modredundant
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    subroutine check_colinear(molec,bond_s,angle_s,natoms0,is_colinear)
    
        !not tested in molecular_tool (a working code is in internal project)
    
        use structure_types
        use alerts
        use constants
        use atomic_geom
        use symmetry
    
        type(str_resmol),intent(inout) :: molec
        integer,dimension(:,:),intent(inout) :: bond_s,angle_s
        integer,intent(in) :: natoms0
        logical,intent(out) :: is_colinear
    
        real(8) :: ang, dist, dx, dy, dz, vx, vy, vz
        integer :: i_1, i_2, i_3
        integer :: i,j,ii, natoms, iat_1, iat_2

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with no Units
        !------------------------------------------------
    
        is_colinear = .false.
        !Run over angle_s triads
        natoms=molec%natoms
        do i=3,natoms
            i_1 = angle_s(i,1)
            i_3 = angle_s(i,2)
            i_2 = angle_s(i,3)
            if (i_1.gt.natoms0 .or. i_2.gt.natoms0 .or. i_3.gt.natoms0) cycle
            ang = calc_atm_angle(molec%atom(i_1),molec%atom(i_3),molec%atom(i_2))
            if (ang.gt.0.97*PI) then
                print*, "Caution: colinear atoms"
                is_colinear=.true.
                print'(3(A,I2,X))', trim(molec%atom(i_1)%name),i_1,&
                                    trim(molec%atom(i_3)%name),i_3,&
                                    trim(molec%atom(i_2)%name),i_2
                print*, "Angle:", ang*180.d0/PI
                !Identify the bonded term in atom%geom
                print*, "==================="
                print*, angle_s(i,1:3)
                print*, "I will cut the connection between:", angle_s(i,2), angle_s(i,3)
                print*, "==================="
                !Add dummy atom
                molec%natoms=molec%natoms+1
                molec%atom(molec%natoms)%name="XX"
                do j=1,molec%geom%nbonds
    !                 print*, molec%geom%bond(j,1:2)
                    if((bond_s(i,1) == molec%geom%bond(j,1) .and. &
                        bond_s(i,2) == molec%geom%bond(j,2)) .or. &
                       (bond_s(i,2) == molec%geom%bond(j,1) .and. &
                        bond_s(i,1) == molec%geom%bond(j,2))) then
                        print*, "bond identified (", j,")", molec%geom%bond(j,1:2)
                        iat_1=molec%geom%bond(j,1)
                        iat_2=molec%geom%bond(j,2)
                        exit
                     endif
                enddo
                !Erase conection iat_1 -- iat_2, and redo as iat_1 -- XX --iat_2
                print*, "Cutting..."
                do ii=1,molec%atom(iat_1)%nbonds
                    if (molec%atom(iat_1)%connect(ii) == iat_2) molec%atom(iat_1)%connect(ii) = molec%natoms
                enddo
                do ii=1,molec%atom(iat_2)%nbonds
                    if (molec%atom(iat_2)%connect(ii) == iat_1) molec%atom(iat_2)%connect(ii) = molec%natoms
                enddo
                molec%atom(molec%natoms)%nbonds = 2
                molec%atom(molec%natoms)%connect(1) = iat_1
                molec%atom(molec%natoms)%connect(2) = iat_2
                dist=calc_atm_dist(molec%atom(iat_1),molec%atom(iat_2))
                vx=(molec%atom(iat_1)%x - molec%atom(iat_2)%x)/2.d0
                vy=(molec%atom(iat_1)%y - molec%atom(iat_2)%y)/2.d0
                vz=(molec%atom(iat_1)%z - molec%atom(iat_2)%z)/2.d0
                if (abs(vx) .lt. 1.d-2) then
                    dy = -vz
                    dz = vy
                elseif (abs(vy) .lt. 1.d-2) then
                    dx = -vz
                    dz = vx
                else
                    dx = -vy
                    dy = vx
                endif
                molec%atom(molec%natoms)%x = (molec%atom(iat_1)%x + molec%atom(iat_2)%x)/2.d0 + dx
                molec%atom(molec%natoms)%y = (molec%atom(iat_1)%y + molec%atom(iat_2)%y)/2.d0 + dy
                molec%atom(molec%natoms)%z = (molec%atom(iat_1)%z + molec%atom(iat_2)%z)/2.d0 + dz
                molec%atom(molec%natoms)%mass = 10000000.d0
                exit
            endif
        enddo
        return
    end subroutine check_colinear
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    subroutine build_Z(molec,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
    
        ! Construct the Z-matrix information from the connectivity
        !
        !Note
        ! Z-matrix is stored inversely as printed. I.e. for atom X
        ! X - i - j - k
        ! bond = /i,X/
        ! angle = /j,i,X/ ...
    
        use structure_types
        use molecular_structure
        use alerts
        use constants
        use atomic_geom
        use symmetry
    
        implicit none
    
        type(str_resmol),intent(inout) :: molec
        integer,dimension(:,:),intent(out) :: bond_s, angle_s, dihed_s
        character(len=*),intent(inout) :: PG    
        integer,dimension(:),intent(in) :: isym
        integer,dimension(:),intent(out) :: bond_sym,angle_sym,dihed_sym
        integer :: Nat
        integer :: i,j,k,m, ii,kk, ibonds, &
                   i_1,i_2,i_3,i_4, ii_1,ii_2,ii_3,ii_4
        logical :: done

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with Atomic Units 
        current_units=molec%units
        call set_geom_units(molec,"Bohr")
        !----------------------------------------------
    
        ! BUILD Z-MATRIX on molecule
        ! New numbering is stored in resseq
        !Otra posibilidad es correr por las conexiones "cerrando" los pasos con las un indentificador
        ! cuando no queden nbonds o cuando nbonds=0 (hay que hacer una cuenta atras de estos...). Pero
        ! hay que repensar como hacer lo simétricos entonces...
        molec%atom(:)%chain = "O"
        !Fisrt atom
        !For symmetric a bond along centre of symmetry should be located (already done for beta-car example fchk)
        ! (30/10/2013: add support for custom symmetry)
        i_1 = -1
        i_2 = -1
        i_3 = -1
        if (adjustl(PG) == "CI" .or. adjustl(PG) == "C02" .or. adjustl(PG) == "CUS") then
            !If natoms is even, a bond must exist between two equiv atoms
            if ( mod(molec%natoms,2) == 0 ) then
                do i=1,molec%natoms
                    do j=1,molec%atom(i)%nbonds
                        k=molec%atom(i)%connect(j)
                        if (k == isym(i)) then
                            if (verbose>0) &
                             print'(X,A,I0,A,I0,/)', "Fisrt bond will be ", i, "--", isym(i)
                            i_1 = i 
                            ! And we know that the first bond have not symmetric
                            bond_sym(2) = 2
                            exit
                        endif
                        if (i_1>0) exit
                    enddo
                    if (i_1>0) exit
                enddo
                if (i_1<0) call alert_msg("fatal","Failure while trying to use symmetry")
            else
                call alert_msg("warning","This combination of geom/symm is under contruction. Take care of your results.")
    
                !Find three atoms that make the brige between symmetric parts
                do i=1,molec%natoms
                    if (i == isym(i)) cycle
                    do j=1,molec%atom(i)%nbonds
                        k=molec%atom(i)%connect(j)
                        do kk=1,molec%atom(isym(i))%nbonds
                            m = molec%atom(isym(i))%connect(kk)
                            if (k == m) then
                                if (verbose>0) then
                                    print*, "Fisrt bond will be ", i, "--", m
                                    print*, "Second bond will be", m, "--", isym(i)
                                    print*, ""
                                endif
                                i_1 = i
                                i_2 = m
                                i_3 = isym(i)
                                ! And we know that the bonds are symmetric
                                bond_sym(2) = 3
                                bond_sym(3) = 2
                                exit
                            endif
                            if (i_1>0) exit
                        enddo
                        if (i_1>0) exit
                    enddo
                    if (i_1>0) exit
                enddo
                if (i_1<0) call alert_msg("fatal","Failure while trying to use symmetry")
            endif
        else
            i_1 = 1
        endif
    
        if (verbose>0) then
            print*, ""
            print*, "PSEUDO Z-MATRIX"
            print*, "================"
        endif

        !k is the counter for the atoms that are being included in Z-matrix
        k=1
        if (verbose>0) &
         print*, k,molec%atom(i_1)%name, i_1
        molec%atom(i_1)%chain = "P"
        molec%atom(i_1)%resseq = k
    
        ! Reset input units before leaving
        call set_geom_units(molec,adjustl(current_units))
        if (molec%natoms == 1) return
    
        !Second atom 
        k=k+1
        if (adjustl(PG) == "CI" .or. adjustl(PG) == "C02" .or. adjustl(PG) == "CUS") then
        !For symmetric is the symmetric... (only if it is not at the centre of symmetry!!!) -- to be fixed
            if (i_2<0) i_2 = isym(i_1)
        else
            i_2 = molec%atom(i_1)%connect(1)
        endif
    
        if (verbose>0) &
         print*, k,molec%atom(i_2)%name, i_2, i_1
        !Save bond
        bond_s(k,1:2) = (/i_1,i_2/)
        !===========================
        molec%atom(i_2)%chain = "P"
        molec%atom(i_2)%resseq = k 
    
        ! Reset input units before leaving
        call set_geom_units(molec,adjustl(current_units))
        if (molec%natoms == 2) return
    
        !Third atom
        ! We need to avoid a blockin situation here, from whihc only an inproper dihedral
        ! would be valid (which are not scanned here). So the third atom should preferably
        ! not be terminal  
        k=k+1
        !Select candidates with the largest number of bonds
        ibonds = 0
        if (i_3<0) then
    
            do i=1,molec%atom(i_1)%nbonds
                ii = molec%atom(i_1)%connect(i)
                if (ii == i_2) cycle
                if (molec%atom(ii)%nbonds > ibonds) then
                    !Pre-selection (double index)
                    ii_1 = i_2
                    ii_2 = i_1
                    i_3 = ii
                    ibonds = molec%atom(ii)%nbonds
                endif
            enddo
            do i=1,molec%atom(i_2)%nbonds
                ii = molec%atom(i_2)%connect(i)
                if (ii == i_1) cycle
                if (molec%atom(ii)%nbonds > ibonds) then
                    ii_1 = i_1
                    ii_2 = i_2
                    i_3 = ii
                    ibonds = molec%atom(ii)%nbonds
                endif
            enddo
            i_1 = ii_1
            i_2 = ii_2
        endif
    ! 
    !     if (molec%atom(i_1)%nbonds > 1) i_3 = molec%atom(i_1)%connect(2)
    !     if (
    ! 
    !     if (molec%atom(i_1)%nbonds > 1) thenS2
    !         ii_1 = i_1
    !         i_1 = i_2
    !         i_2 = ii_1
    !         i_3 = molec%atom(i_1)%connect(2)
    !     else
    !         i_1 = i_1
    !         i_2 = i_2
    !         i_3 = molec%atom(i_2)%connect(1)
    !         if (i_3 == i_1) i_3 = molec%atom(i_2)%connect(2)
    !     endif
        if (verbose>0) &
         print*, k,molec%atom(i_3)%name, i_3, i_2, i_1 
        !Save bonded
        bond_s(k,1:2) = (/i_2,i_3/) 
        angle_s(k,1:3) = (/i_1,i_2,i_3/) 
        !=======================
        molec%atom(i_3)%chain = "P"
        molec%atom(i_3)%resseq = k
    
        if (molec%natoms == 3) then
            angle_sym(3) = 3 
            ! Reset input units before leaving
             call set_geom_units(molec,adjustl(current_units))
            return
        endif
    
        !Only for symmetric if natoms is even!
        if ((adjustl(PG) == "CI" .or. adjustl(PG) == "C02" .or. adjustl(PG) == "CUS") &
            .and. mod(molec%natoms,2) == 0 ) then
            k=k+1
            ii_4 = isym(i_3)
            ii_3 = isym(i_2)
            ii_2 = isym(i_1)
            ii_1 = i_3
            if (verbose>0) &
             print*, k,molec%atom(ii_4)%name, ii_4, ii_3, ii_2, ii_1, "!sym"
            molec%atom(ii_4)%chain = "P"
            molec%atom(ii_4)%resseq = k
            !Save bonded
            bond_s(k,1:2)  = (/ii_3,ii_4/) 
            angle_s(k,1:3) = (/ii_2,ii_3,ii_4/) 
            dihed_s(k,1:4) = (/ii_1,ii_2,ii_3,ii_4/)
            !=======================
            bond_sym(k-1) = k
            bond_sym(k) = k-1
            angle_sym(k-1) = k
            angle_sym(k) = k-1
            !This has not symmetric... (only if it is not at the centre of symmetry!!!) -- to be fixed
            dihed_sym(k) = k
        endif
    
        kk=0
        done=.false.
        do while (.not.done)
            kk=kk+1
            done=.true.
            do i=2,molec%natoms
                if (molec%atom(i)%chain == "P") cycle
                !Look for the dihedrals
                do j=1,molec%geom%ndihed
                    i_1 = molec%geom%dihed(j,1)
                    i_2 = molec%geom%dihed(j,2)
                    i_3 = molec%geom%dihed(j,3)
                    i_4 = molec%geom%dihed(j,4)
                    if((molec%atom(i_1)%chain == "P" .and.&
                        molec%atom(i_2)%chain == "P" .and.& 
                        molec%atom(i_3)%chain == "P" .and.&
                        i_4 == i )) then
                        !Add line to Z-mat (and symmetric if flag active)
                        k=k+1
                        if (verbose>0) &
                         print*, k,molec%atom(i)%name, molec%geom%dihed(j,4:1:-1)
                        !Save bonded
                        bond_s(k,1:2)  = (/i_3,i_4/) 
                        angle_s(k,1:3) = (/i_2,i_3,i_4/) 
                        dihed_s(k,1:4) = (/i_1,i_2,i_3,i_4/)
                        !=======================
                        molec%atom(i)%chain = "P"
                        molec%atom(i)%resseq = k
                      if (adjustl(PG) == "CI" .or. adjustl(PG) == "C02" .or. adjustl(PG) == "CUS") then
                       !Only for symmetric...
                        k=k+1
                        ii_1 = isym(i_1)
                        ii_2 = isym(i_2)
                        ii_3 = isym(i_3)
                        ii_4 = isym(i_4)
                        if (verbose>0) &
                         print*, k,molec%atom(i_4)%name, ii_4, ii_3, ii_2, ii_1, "!sym"
                        molec%atom(ii_4)%chain = "P"
                        molec%atom(ii_4)%resseq = k
                        !Save bonded
                        bond_s(k,1:2)  = (/ii_3,ii_4/) 
                        angle_s(k,1:3) = (/ii_2,ii_3,ii_4/) 
                        dihed_s(k,1:4) = (/ii_1,ii_2,ii_3,ii_4/)
                        !=======================
                        bond_sym(k-1) = k
                        bond_sym(k) = k-1
                        angle_sym(k-1) = k
                        angle_sym(k) = k-1
                        dihed_sym(k-1) = k
                        dihed_sym(k) = k-1
                      endif
                        !This exit ensures a .not.done if any happen to not be done
                        exit
                    !if clause splitted to ensure correct order, thus correct selection of bonds and angles
                    else if((molec%atom(i_4)%chain == "P" .and.&
                             molec%atom(i_3)%chain == "P" .and.& 
                             molec%atom(i_2)%chain == "P" .and.&
                             i_1 == i )) then
                        !Add line to Z-mat (and symmetric)
                        k=k+1
                        if (verbose>0) &
                         print*, k,molec%atom(i)%name, molec%geom%dihed(j,1:4) !, "!reversed"
                        !Save bonded
                        bond_s(k,1:2)  = (/i_2,i_1/) 
                        angle_s(k,1:3) = (/i_3,i_2,i_1/) 
                        dihed_s(k,1:4) = (/i_4,i_3,i_2,i_1/)
                        !=======================
                        molec%atom(i)%chain = "P"
                        molec%atom(i)%resseq = k
                      if (adjustl(PG) == "CI" .or. adjustl(PG) == "C02" .or. adjustl(PG) == "CUS") then
                       !Only for symmetric...
                        k=k+1
                        ii_1 = isym(i_1)
                        ii_2 = isym(i_2)
                        ii_3 = isym(i_3)
                        ii_4 = isym(i_4)
                        if (verbose>0) &
                         print*, k,molec%atom(i_4)%name, ii_1, ii_2, ii_3, ii_4, "!sym"
                        !Save bonded
                        bond_s(k,1:2)  = (/ii_2,ii_1/) 
                        angle_s(k,1:3) = (/ii_3,ii_2,ii_1/) 
                        dihed_s(k,1:4) = (/ii_4,ii_3,ii_2,ii_1/)
                        !=======================
                        bond_sym(k-1) = k
                        bond_sym(k) = k-1
                        angle_sym(k-1) = k
                        angle_sym(k) = k-1
                        dihed_sym(k-1) = k
                        dihed_sym(k) = k-1
                        !=======================
                        molec%atom(ii_1)%chain = "P"
                        molec%atom(ii_1)%resseq = k
                      endif
                        !This exist ensures a .not.done if any happen to not be done
                        exit
                    endif   
                    
                enddo !Cycle over stored dihedrals
                if (j>molec%geom%ndihed) done = .false.
            enddo
            if (kk>50) call alert_msg("fatal","Z-matrix could not be formed")
        enddo
        if (verbose>0) &
         print'(X,A,I0,A,/)', "Success in ", kk, " cycles"

        !Update molec%geom with Zmat
        Nat = molec%natoms
        molec%geom%bond(1:Nat-1,1:2)  = bond_s(2:Nat,1:2)
        molec%geom%angle(1:Nat-2,1:3) = angle_s(3:Nat,1:3)
        molec%geom%dihed(1:Nat-3,1:4) = dihed_s(4:Nat,1:4)
        molec%geom%nbonds  = Nat-1
        molec%geom%nangles = Nat-2
        molec%geom%ndihed  = Nat-3

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! Revert original units
        call set_geom_units(molec,adjustl(current_units))
        !----------------------------------------------
    
        return
    
    end subroutine build_Z
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    subroutine read_Z(unitZ,molec,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
    
        ! Read the Z-matrix information from the pseudo-Zmat constructed by build_Z
        !
        !Note
        ! Z-matrix is stored inversely as printed. I.e. for atom X
        ! X - i - j - k
        ! bond = /i,X/
        ! angle = /j,i,X/ ...
        !
        ! Some features need to be included (seems easy): 
        !   - Inherit sym internals (given by !sym comment)
        !
        ! TODO: bond_s,angle_s,dihed_s can be local
    
        use structure_types
        use alerts
        use constants
        use atomic_geom
        use symmetry
    
        implicit none
    
        type(str_resmol),intent(inout) :: molec
        integer,dimension(:,:),intent(out) :: bond_s, angle_s, dihed_s
        character(len=*),intent(inout) :: PG    
        integer,dimension(:),intent(in) :: isym
        integer,dimension(:),intent(out) :: bond_sym,angle_sym,dihed_sym
        integer,intent(in) :: unitZ
        integer :: Nat
        integer :: i,j,k, ii,kk, ibonds, &
                   i_1,i_2,i_3,i_4, ii_1,ii_2,ii_3,ii_4
        logical :: done
        character(len=10) :: dummy_char

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with no Units
        !------------------------------------------------
    
        !First element
        read(unitZ,*) k, dummy_char, i_1
        !Second element: start bond_s counting
        read(unitZ,*) k, dummy_char, i_2, i_1
        bond_s(k,1:2) = (/i_1,i_2/) 
        !Third element: start angle_s counting
        read(unitZ,*) k, dummy_char, i_3, i_2, i_1
        bond_s(k,1:2) = (/i_2,i_3/) 
        angle_s(k,1:3) = (/i_1,i_2,i_3/)
        !Rest of element: include dihed_s 
        do kk=k+1,molec%natoms
            read(unitZ,*) k, dummy_char, i_4, i_3, i_2, i_1
            bond_s(k,1:2)  = (/i_3,i_4/) 
            angle_s(k,1:3) = (/i_2,i_3,i_4/) 
            dihed_s(k,1:4) = (/i_1,i_2,i_3,i_4/)
        enddo
    
        !Update molec%geom with Zmat
        Nat = molec%natoms
        molec%geom%bond(1:Nat-1,1:2)  = bond_s(2:Nat,1:2)
        molec%geom%angle(1:Nat-2,1:3) = angle_s(3:Nat,1:3)
        molec%geom%dihed(1:Nat-3,1:4) = dihed_s(4:Nat,1:4)
        molec%geom%nbonds  = Nat-1
        molec%geom%nangles = Nat-2
        molec%geom%ndihed  = Nat-3

        !To be done: inherit symmetry. For the moment:
        bond_sym=0
        angle_sym=0
        dihed_sym=0
     
        return
    
    end subroutine read_Z
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    subroutine zmat2cart(molec,S)
    
        ! Transform z-matrix to cartesian coordinates
        ! Run with atomic coordinates (for S)
        !  * input : bohr (length) and radians (angles)
        !  * outout: bohr (cartesian coods)
        ! NOTES
        ! Molec has to include the ic information (atom%geom)
    
        use structure_types
        use constants
        use matrix
        use verbosity
    
        integer,parameter :: NDIM = 600
    
        !Arguments
        type(str_resmol),intent(inout) :: molec
        real(8),dimension(NDIM),intent(in) :: S
       
        !Local
        integer,dimension(1:NDIM,1:4) :: bond_s,angle_s,dihed_s
        real(8) :: det

        real(8) :: R,angle,dihed, xcom,ycom,zcom, MASS,xaux,yaux,zaux
    
        real(8),dimension(3,3) :: PIner,Aux, Aux3
        real(8),dimension(3) :: Vec
    
        real(8),dimension(NDIM,NDIM) :: work
        integer,dimension(NDIM) :: ipiv
        integer :: info
    
        integer :: Nat
    
        integer :: i, i_at, i_b,i_a,i_d

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine returns data in Atomic Units
        !------------------------------------------------
    
        ! Setup things 
        Nat = molec%natoms
        bond_s(2:Nat,1:2)  = molec%geom%bond(1:Nat-1,1:2)
        angle_s(3:Nat,1:3) = molec%geom%angle(1:Nat-2,1:3)
        dihed_s(4:Nat,1:4) = molec%geom%dihed(1:Nat-3,1:4)
    
        !Atom1
    !     i_at = atom_order(1)
        i_at = bond_s(2,1)
        i1 = bond_s(2,1)
        molec%atom(i_at)%x = 0.d0
        molec%atom(i_at)%y = 0.d0
        molec%atom(i_at)%z = 0.d0
        !Atom2
        R0 = S(1)*BOHRtoAMS
    !     i_at = atom_order(2)
        i_at = bond_s(2,2)
        molec%atom(i_at)%x = -R0
        molec%atom(i_at)%y = 0.d0
        molec%atom(i_at)%z = 0.d0
        !Atom3
        R = S(2)*BOHRtoAMS
        angle = S(Nat)
        i_a = angle_s(3,1)
        i_b = angle_s(3,2)
    !     i_at = atom_order(3)
        i_at = angle_s(3,3) !<-- this info is here!
        if (i_b == i1 ) then !atom_order(1)) then
    ! print*, "connected to i1", R
            molec%atom(i_at)%x = R*dcos(PI-angle)
            molec%atom(i_at)%y = R*dsin(PI-angle)
            molec%atom(i_at)%z = 0.d0
        else
    ! print*, i_at, "connected to other than i1",i_b, R
            molec%atom(i_at)%x = -R*dcos(PI-angle)-R0
            molec%atom(i_at)%y = R*dsin(PI-angle)
            molec%atom(i_at)%z = 0.d0
        endif
    
        !COLOCACIÓN DEL RESTO DE ÁTOMOS
        k=0
        do i=4,molec%natoms
            k=k+1
            i_d = dihed_s(i,1)
            i_a = dihed_s(i,2)
            i_b = dihed_s(i,3)
            i_at = dihed_s(i,4)
            R = S(2+k)*BOHRtoAMS
            angle = S(Nat+k)
            dihed = S(2*Nat-3+k)
    ! print'(A,4I3,A,F8.2)', "Setting", i_at, i_b, i_a, i_d, " with dihed =", dihed*180./PI
    ! print*, "Atom", i_at, "conected to", i_b, i_a, i_d, R, angle*180./PI, dihed*180./PI
            call addcart(molec%atom(i_b),molec%atom(i_a),molec%atom(i_d),&
                         R,angle,dihed, molec%atom(i_at))
        enddo
    
        !THE FOLLOWING PLACES ALL IN  STD ORI (IS NOT FIX!
        !Place on COM
        xcom = 0.d0
        ycom = 0.d0
        zcom = 0.d0
        MASS = 0.d0
        do i=1,molec%natoms
            xcom = xcom + molec%atom(i)%x*molec%atom(i)%mass
            ycom = ycom + molec%atom(i)%y*molec%atom(i)%mass
            zcom = zcom + molec%atom(i)%z*molec%atom(i)%mass
            MASS = MASS + molec%atom(i)%mass
        enddo
        xcom = xcom/MASS
        ycom = ycom/MASS
        zcom = zcom/MASS
        do i=1,molec%natoms
            molec%atom(i)%x = molec%atom(i)%x - xcom
            molec%atom(i)%y = molec%atom(i)%y - ycom
            molec%atom(i)%z = molec%atom(i)%z - zcom
        enddo
        !Calculate the intertia tensor
        PIner(1,1) = 0.d0
        do i=1,molec%natoms
            PIner(1,1) = PIner(1,1) + molec%atom(i)%mass*&
                         (molec%atom(i)%y**2+molec%atom(i)%z**2)
        enddo
        PIner(2,2) = 0.d0
        do i=1,molec%natoms
            PIner(2,2) = PIner(2,2) + molec%atom(i)%mass*&
                         (molec%atom(i)%x**2+molec%atom(i)%z**2)
        enddo
        PIner(3,3) = 0.d0
        do i=1,molec%natoms
            PIner(3,3) = PIner(3,3) + molec%atom(i)%mass*&
                         (molec%atom(i)%y**2+molec%atom(i)%x**2)
        enddo
        PIner(1,2) = 0.d0
        do i=1,molec%natoms
            PIner(1,2) = PIner(1,2) - molec%atom(i)%mass*&
                         (molec%atom(i)%x*molec%atom(i)%y)
        enddo
        PIner(2,1) = PIner(1,2)
        PIner(1,3) = 0.d0
        do i=1,molec%natoms
            PIner(1,3) = PIner(1,3) - molec%atom(i)%mass*&
                         (molec%atom(i)%x*molec%atom(i)%z)
        enddo
        PIner(3,1) = PIner(1,3)
        PIner(2,3) = 0.d0
        do i=1,molec%natoms
            PIner(2,3) = PIner(2,3) - molec%atom(i)%mass*&
                         (molec%atom(i)%y*molec%atom(i)%z)
        enddo
        PIner(3,2) = PIner(2,3)
        !Obtein the trasformation to principal axis of inertia
        call diagonalize_full(PIner(1:3,1:3),3,Aux(1:3,1:3),Vec(1:3),"lapack")
        if (verbose>1) then
        write(6,*) ""
        write(6,*) "I"
        do i=1,3
            write(6,'(100(F9.5,X))') PIner(i,1:3)
        enddo
        write(6,*) ""
        write(6,*) "Xrot"
        do i=1,3
            write(6,'(100(F9.5,X))') Aux(i,1:3)
        enddo
        endif
        !Check Determinant of the rotation
        det = Aux(1,1)*Aux(2,2)*Aux(3,3)
        det = det + Aux(1,2)*Aux(2,3)*Aux(3,1)
        det = det + Aux(2,1)*Aux(3,2)*Aux(1,3)
        det = det - Aux(3,1)*Aux(2,2)*Aux(1,3)
        det = det - Aux(2,1)*Aux(1,2)*Aux(3,3)
        det = det - Aux(3,2)*Aux(2,3)*Aux(1,1)
    !     print*, "DETERMINANT=", det
        
        !Rotate the structure_types to the principal axis of inertia frame
        !Using the traspose of the eigenvalue matrix. Right?!
        do i=1,molec%natoms
            xaux = molec%atom(i)%x
            yaux = molec%atom(i)%y
            zaux = molec%atom(i)%z
            molec%atom(i)%x = Aux(1,1)*xaux+Aux(2,1)*yaux+zaux*Aux(3,1)
            molec%atom(i)%y = Aux(1,2)*xaux+Aux(2,2)*yaux+zaux*Aux(3,2)
            molec%atom(i)%z = Aux(1,3)*xaux+Aux(2,3)*yaux+zaux*Aux(3,3)
    !         molec%atom(i)%x = Aux(1,1)*xaux+Aux(1,2)*yaux+zaux*Aux(1,3)
    !         molec%atom(i)%y = Aux(2,1)*xaux+Aux(2,2)*yaux+zaux*Aux(2,3)
    !         molec%atom(i)%z = Aux(3,1)*xaux+Aux(3,2)*yaux+zaux*Aux(3,3)
        enddo
    
        if (det < 0) then
        do i=1,molec%natoms
            molec%atom(i)%x = -molec%atom(i)%x 
            molec%atom(i)%y = -molec%atom(i)%y 
            molec%atom(i)%z = -molec%atom(i)%z 
        enddo
        endif
    
        do i=1,molec%natoms
            !Convert back to a.u.
            molec%atom(i)%x = molec%atom(i)%x/BOHRtoAMS
            molec%atom(i)%y = molec%atom(i)%y/BOHRtoAMS
            molec%atom(i)%z = molec%atom(i)%z/BOHRtoAMS
        enddo
    
        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine get data in Atomin Units
        molec%units="Bohr"
        !----------------------------------------------
    
        return
    
    
    end subroutine zmat2cart
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    subroutine addcart(atom_1,atom_2,atom_3,&
                       bond,angle,dihed,atom_New)
    
        !Description: add the cartesian coordinates of atom_New given the position of atoms 1,2&3, so that
        ! the connectivity is:
        !   atom_New -- atom_1 -- atom_2 -- atom_3
    ! Angles in radians
    
        use structure_types
        use constants
        use matrix
    
        implicit none
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: ZERO = 1.d-10
    
        type(str_atom),intent(in) :: atom_1,atom_2,atom_3
        type(str_atom),intent(inout) :: atom_New 
        real(8),intent(in) :: bond,angle,dihed
    
        type(str_atom) :: aux_atom_1,aux_atom_2,aux_atom_3
        real(8) :: v1,v2,v3,sinTheta,Theta,baux,baux2,dihed_aux
    
        integer :: i
    
        !New frame (only for the atoms involved)
        !1. Traslate to have atom i_b (atom1) at origin
    !     aux_atom_1%x = atom_1%x - atom_1%x
    !     aux_atom_1%y = atom_1%y - atom_1%y
    !     aux_atom_1%z = atom_1%z - atom_1%z
    ! !======
        aux_atom_2%x = atom_2%x - atom_1%x
        aux_atom_2%y = atom_2%y - atom_1%y
        aux_atom_2%z = atom_2%z - atom_1%z
        aux_atom_3%x = atom_3%x - atom_1%x
        aux_atom_3%y = atom_3%y - atom_1%y
        aux_atom_3%z = atom_3%z - atom_1%z
        !2. Rotate to place atom i_a on -x axis
        ! Find rotation axis and angle with vector product of position(atom i_a) X (-1,0,0)
        v1 =  0.d0
        v2 = -aux_atom_2%z
        v3 =  aux_atom_2%y
        !If v2=v3=0, any arbitray vector perpedicular to xaxis is valid
        if (v3==0.d0 .and. v2==0.d0) v3=1.d0
        baux     = dsqrt(aux_atom_2%x**2+aux_atom_2%y**2+aux_atom_2%z**2)
        sinTheta = dsqrt(aux_atom_2%y**2+aux_atom_2%z**2)/baux
        Theta    = dasin(sinTheta)
        !Fix rotation
    
        if (aux_atom_2%x > 0) Theta = PI-Theta
    !     if (aux_atom_2%x > 0) Theta = Theta+PI
    !     if (aux_atom_2%x > 0) then
    !         Theta = Theta+PI
    !         print*, "In int2cart. WARNING: rotation angle changed to:", theta*180./PI 
    !     endif
    
    
        !Rotate
        call rotation_3D(aux_atom_3%x,&
                         aux_atom_3%y,&
                         aux_atom_3%z,&
                         v1,v2,v3,Theta)
        call rotation_3D(aux_atom_2%x,&
                         aux_atom_2%y,&
                         aux_atom_2%z,&
                         v1,v2,v3,Theta)
        if (aux_atom_2%x>0.d0) then
            print*, "ERROR IN int2cart"
            stop
        endif
    
        atom_New%x = bond * dcos(PI-angle)
        baux  = bond * dsin(PI-angle)
        baux2 = dsqrt(aux_atom_3%y**2+aux_atom_3%z**2)
        !This will be unstable for linear molecules
        atom_New%y = aux_atom_3%y/baux2*baux
        atom_New%z = aux_atom_3%z/baux2*baux
        !Rotation around X axis 
        ! (the sign of the rotation is adjusted empirically! Caution!!!)
    !     dihed_aux = dsign(dihed,aux_atom_3%z)
        call rotation_3D(atom_New%x,&
                         atom_New%y,&
                         atom_New%z,&
                         1.d0,0.d0,0.d0,dihed)
    
        !Rotate and traslate back
        !Only rotate if needed
    !     if (Theta == PI) then
    !         atom_New%x = -atom_New%x
    !         atom_New%y = -atom_New%y
    !         atom_New%z = -atom_New%z
    !     else if (Theta /= 0.d0) then
        call rotation_3D(atom_New%x,&
                         atom_New%y,&
                         atom_New%z,&
                         v1,v2,v3,-Theta)
    !     endif
        atom_New%x = atom_New%x + atom_1%x
        atom_New%y = atom_New%y + atom_1%y
        atom_New%z = atom_New%z + atom_1%z
    
    
        return
    
    end subroutine addcart

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

    subroutine rmsd_fit_frame(molec,molec2,info)
    
        ! THIS SR ARE HERE SINCE THEY ARE USED ALONG WITH zmat2cart
        ! Check if two molecules share the same orientation and, if not, reorinent molec to fit molec2
        ! Thia version (4) uses ROTATA to get a first stimation of the rotation
        ! and then converts it to a rotation keeping the same axis
    
        use structure_types
        use molecular_structure
    
        implicit none
    
        type(str_resmol),intent(inout) :: molec
        type(str_resmol),intent(in)    :: molec2
        integer,intent(out) :: info
    
        real(8),dimension(3,3) :: T0, R
        real(8),dimension(3) :: Vec2
    
        integer :: i,j, Nat
        real(8) :: xaux, yaux, zaux, x2aux, y2aux, z2aux,&
                   ux,uy,uz, umod, baux, baux2, Theta, costheta, dist

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with any Units
        !------------------------------------------------
    
        Nat = molec%natoms
        info=0
    
        ! Get rotation to minimize RMSD
        call ROTATA1(molec,molec2,R)
    
        ! Analyze the rotation
        if (verbose>1) &
         call MAT0(6,R,3,3,"R (to minimize RMSD)")
        do i=1,3
            Theta=0.d0
            do j=1,3
                if (abs(R(i,j)) > 0.8) then
                    T0(i,j) = sign(1.d0,R(i,j))
                    Theta = Theta + 1.d0
                else 
                    T0(i,j) = 0.d0
                endif
            enddo
            if (Theta /= 1.d0) info=-1
        enddo
        if (verbose>1) &
         call MAT0(6,T0,3,3,"T0 (matching stdori frame)")
    
        ! If successful, rotate the molecule
        if (info /= 0) then
            call alert_msg("warning","New fast method failed.")
        else
            ! Rotate
            call rotate_molec(molec,T0)
        endif
        
    
        return
    
    end subroutine rmsd_fit_frame


    subroutine rmsd_fit_frame_brute(molec,molec2,dist)
    
        ! THIS SR ARE HERE SINCE THEY ARE USED ALONG WITH zmat2cart
        ! Check if two molecules share the same orientation and, if not, reorinent molec to fit molec2
        ! This version "ensures" proper relative orientation (but is very time consuming!!)
        ! New: version b includes a distance threshold to avoid checking all the orientations
    
        use structure_types
    
        implicit none
    
        type(str_resmol),intent(inout) :: molec
        type(str_resmol),intent(in)    :: molec2
        real(8),intent(inout)          :: dist

        integer                        :: info
    
!         type(str_resmol)    :: molec_aux
        real(8),dimension(3,3) :: T
        real(8),dimension(3,3) :: T0
        real(8),dimension(3) :: Vec2
    
        integer :: iT, iTT, i, Nat
        real(8) :: xaux, yaux, zaux, det, rsum, rcheck, &
                   xrot, yrot, zrot

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with any Units
        !------------------------------------------------
    
        if (verbose>0) &
         print*, "Using brute force RMSD fit subroutine"

        ! Set initially to failure
        info=1

        Nat = molec%natoms
        rcheck=1.d10
        do iTT=1,7
            !Check all 8 possible solutions for T0 (see Sando, 2001) 
            do iT=1,8
            
                Vec2=1.d0
                if (iT==2) then
                    Vec2(1) = -1.d0
                else if (iT==3) then
                    Vec2(2) = -1.d0
                else if (iT==4) then
                    Vec2(3) = -1.d0
                else if (iT==5) then
                    Vec2(1) = -1.d0
                    Vec2(2) = -1.d0
                else if (iT==6) then
                    Vec2(1) = -1.d0
                    Vec2(3) = -1.d0
                else if (iT==7) then
                    Vec2(2) = -1.d0
                    Vec2(3) = -1.d0
                else if (iT==8) then
                    Vec2(1) = -1.d0
                    Vec2(2) = -1.d0
                    Vec2(3) = -1.d0
                endif
            
                !T is diagonal 
                T=0
                if (iTT ==1 ) then
                T(1,1) = Vec2(1)
                T(2,2) = Vec2(2)
                T(3,3) = Vec2(3)
                elseif (iTT ==2 ) then
                T(1,2) = Vec2(1)
                T(2,1) = Vec2(2)
                T(3,3) = Vec2(3)
                elseif (iTT ==3 ) then
                T(2,1) = Vec2(1)
                T(2,1) = Vec2(2)
                T(3,3) = Vec2(3)
                elseif (iTT ==4 ) then
                T(1,3) = Vec2(1)
                T(2,2) = Vec2(2)
                T(3,1) = Vec2(3)
                elseif (iTT ==5 ) then
                T(3,1) = Vec2(1)
                T(2,2) = Vec2(2)
                T(1,3) = Vec2(3)
                elseif (iTT ==6 ) then
                T(1,1) = Vec2(1)
                T(2,3) = Vec2(2)
                T(3,2) = Vec2(3)
                elseif (iTT ==7 ) then
                T(1,1) = Vec2(1)
                T(3,2) = Vec2(2)
                T(2,3) = Vec2(3)
                endif
                
                ! Compute Determinant
                det =       T(1,1)*T(2,2)*T(3,3)
                det = det + T(2,1)*T(3,2)*T(1,3)
                det = det + T(1,2)*T(2,3)*T(3,1)
                det = det - T(3,1)*T(2,2)*T(1,3)
                det = det - T(2,1)*T(1,2)*T(3,3)
                det = det - T(3,2)*T(2,3)*T(1,1)
!                 print*, det
            
                if (det<0.d0) cycle
                
                !Rotate and check sum of distances
                rsum=0.d0
                do i=1,Nat
                    xaux = molec%atom(i)%x
                    yaux = molec%atom(i)%y
                    zaux = molec%atom(i)%z
                    xrot = T(1,1)*xaux + T(1,2)*yaux + T(1,3)*zaux - molec2%atom(i)%x
                    yrot = T(2,1)*xaux + T(2,2)*yaux + T(2,3)*zaux - molec2%atom(i)%y
                    zrot = T(3,1)*xaux + T(3,2)*yaux + T(3,3)*zaux - molec2%atom(i)%z
                    rsum=rsum + dsqrt( xrot**2+ yrot**2+ zrot**2)
                enddo
    !             print*, "Sum of distance ", rsum
                
                if (rsum < rcheck) then
                    T0=T
                    rcheck=rsum
                endif
                
                if (rsum < dist) exit
            
            enddo !
            
            !If exited due to criterion met, exit here
            if (rsum < dist) then
                info=0
                exit
            endif
    
        enddo ! Possible T values

        !Update the distance criterion
        if (verbose>1.and.info/=0) then
            print*, "Criterion not met"
            print*, "  Distance criterion was: ", dist
            print*, "  Updated to: ", rcheck
        endif
        dist=rcheck
        if (info/=0) return
     
        if (verbose>1) &
         call MAT0(6,T0,3,3,"T0")

        do i=1,Nat
            xaux = molec%atom(i)%x
            yaux = molec%atom(i)%y
            zaux = molec%atom(i)%z
            molec%atom(i)%x = T0(1,1)*xaux + T0(1,2)*yaux + T0(1,3)*zaux 
            molec%atom(i)%y = T0(2,1)*xaux + T0(2,2)*yaux + T0(2,3)*zaux
            molec%atom(i)%z = T0(3,1)*xaux + T0(3,2)*yaux + T0(3,3)*zaux 
        enddo
    
        return
    
    end subroutine rmsd_fit_frame_brute


end module zmat_manage