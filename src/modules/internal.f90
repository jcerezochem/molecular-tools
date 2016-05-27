module internal_module

    use molecular_structure
    use matrix
    use matrix_print
    use verbosity

    implicit none

    contains


    subroutine define_internal_set(molec,def_internal,intfile,rmzfile,use_symmetry,isym,S_sym,Ns,Nf,Fltr)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  
        !--------------------------------------------------------------

        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use matrix
        use zmat_manage
    
        implicit none
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp = 1.d-10 !practically zero
    
        !====================== 
        !ARGUMENTS
        type(str_resmol),intent(inout)     :: molec             ! Input molecule (but only use geom...) 
        character(len=*),intent(in)        :: def_internal      ! switch with the way internal are selected
        character(len=*),intent(in)        :: intfile           ! file with the def of internal coords
        character(len=*),intent(in)        :: rmzfile           ! additional file with coordinates to remove from Zmat
        logical,intent(in)                 :: use_symmetry      ! logical to use or not symmetry
        integer,dimension(:),intent(in)    :: isym              ! array with the simmetric atom (in)
        integer,dimension(:),intent(out)   :: S_sym             ! array with symmetric internals (out)
        integer,intent(out)                :: Ns                ! Total number of internal coordiantes to use
        integer,intent(out)       :: Nf
        real(8),dimension(:,:)    :: Fltr
        !====================== 
    
        !======================
        !LOCAL 
        !System info
        integer :: Nat
        character(len=len(def_internal)) :: def_internal_local
        !Counters
        integer :: i,j,k, ii
        !Matrices to store/manage internal coordianates
        integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
        integer,dimension(NDIM) :: bond_sym,angle_sym,dihed_sym
        !Symm
        character(len=5) :: PG
        ! Zmat elements removal
        character :: rm_type
        integer :: rm_zline, Nrm, nbonds_rm, nangles_rm, ndiheds_rm
        integer,dimension(100) :: bond_rm, angle_rm, dihed_rm
        ! File IO
        integer :: I_FILE=70
        integer :: IOstatus
        !=============

        ! Preliminar things
        def_internal_local=def_internal
        call set_word_upper_case(def_internal_local)
        Nat = molec%natoms
        PG = molec%PG

        Nf = 0

        !GEN BONDED SET FOR INTERNAL COORD
        if (adjustl(def_internal_local) == "SEL") then
            if (verbose>0) &
             print*, "Reading internal coordianates from: "//trim(adjustl(intfile))
            open(I_FILE,file=intfile,iostat=IOstatus,status='old') 
            if (IOstatus /= 0) call alert_msg("fatal","Cannot open file: "//trim(adjustl(intfile)))
            ! Preprocess the file to get combination labels
            call pp_modredundant(I_FILE,molec,Nf,Fltr)
            ! Get internal coords (using modredundant sr)
!             call modredundant(I_FILE,molec)
            close(I_FILE)

        elseif (adjustl(def_internal_local) == "ZMAT") then
            if (adjustl(intfile) == "none") then
                call build_Z(molec,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
            else
                open(I_FILE,file=intfile,status="old")
                print*, "Z-matrix read from "//trim(adjustl(intfile))
                call read_Z(I_FILE,molec,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
                close(I_FILE)
            endif

        else if (adjustl(def_internal_local) == "ALL") then!otherwise all parameters are used
            if (verbose>0) &
             print*, "Using all internal coordinates", molec%geom%nbonds+molec%geom%nangles+molec%geom%ndihed
        else
            call alert_msg("fatal","Unkownn option for internal definition. Valid options are 'zmat', 'sel', 'all'")
        endif 

        ! Compute the number of internal coords
        Ns = molec%geom%nbonds  + &
             molec%geom%nangles + &
             molec%geom%ndihed  + &
             molec%geom%nimprop


        ! Remove some Zmat elements if required
        if (def_internal_local=="ZMAT" .and. rmzfile /= "none") then
            open(I_FILE,file=rmzfile,status='old')
            read(I_FILE,*) Nrm
            ! First, read lines to remove
            nbonds_rm  = 0
            nangles_rm = 0
            ndiheds_rm = 0
            do i=1,Nrm
                read(I_FILE,*) rm_zline, rm_type
                if (rm_type == "B") then
                    nbonds_rm = nbonds_rm + 1
                    bond_rm(nbonds_rm) = rm_zline
                elseif (rm_type == "A") then
                    nangles_rm = nangles_rm + 1
                    angle_rm(nangles_rm) = rm_zline
                elseif (rm_type == "D") then
                    ndiheds_rm = ndiheds_rm + 1
                    dihed_rm(ndiheds_rm) = rm_zline
                endif
            enddo
            close(I_FILE)
            ! Short remove lists
            call sort_vec_int(bond_rm(1:nbonds_rm),nbonds_rm,reverse=.true.)
            call sort_vec_int(angle_rm(1:nangles_rm),nangles_rm,reverse=.true.)
            call sort_vec_int(dihed_rm(1:ndiheds_rm),ndiheds_rm,reverse=.true.)
            ! The remove
            if (verbose>0) &
             print*, "Removing bonds from lines:"
            do i=1,nbonds_rm
                if (verbose>0) &
                 print*, bond_rm(i)
                rm_zline = bond_rm(i) - 1
                do ii = rm_zline, molec%geom%nbonds-1
                    molec%geom%bond(ii,1:2) = molec%geom%bond(ii+1,1:2)
                enddo
                molec%geom%nbonds  = molec%geom%nbonds  - 1
                Ns = Ns - 1
            enddo
            if (verbose>0) &
             print*, "Removing angles from lines:"
            do i=1,nangles_rm
                if (verbose>0) &
                 print*, angle_rm(i)
                rm_zline = angle_rm(i) - 2
                do ii = rm_zline, molec%geom%nangles-1
                    molec%geom%angle(ii,1:3) = molec%geom%angle(ii+1,1:3)
                enddo
                molec%geom%nangles  = molec%geom%nangles  - 1
                Ns = Ns - 1
            enddo
            if (verbose>0) &
             print*, "Removing dihedrals from lines:"
            do i=1,ndiheds_rm
                if (verbose>0) & 
                 print*, dihed_rm(i)
                rm_zline = dihed_rm(i) - 3
                do ii = rm_zline, molec%geom%ndihed-1
                    molec%geom%dihed(ii,1:4) = molec%geom%dihed(ii+1,1:4)
                enddo
                molec%geom%ndihed  = molec%geom%ndihed  - 1
                Ns = Ns - 1
            enddo
        endif

        !Set symmetry of internal (only if symmetry is detected)
        if (use_symmetry) then
            do i=1,Nat-1
                S_sym(i) = bond_sym(i+1)-1
            enddo
            do i=1,Nat-2
                S_sym(i+Nat-1) = angle_sym(i+2)+Nat-3
            enddo
            do i=1,Nat-3
                S_sym(i+2*Nat-3) = dihed_sym(i+3)+2*Nat-6
            enddo
        endif

        if (verbose > 0) then
            print'(/,X,A)', "----------------------------------------------------------"
            print*, "Constructed set of valence internal coordianated"
            print'(X,A,A )', "  Type    = ", def_internal_local
            print'(X,A,I0)', "  Ns      = ", Ns
            print'(X,A,I0)', "  Nbonds  = ", molec%geom%nbonds
            print'(X,A,I0)', "  Nangles = ", molec%geom%nangles
            print'(X,A,I0)', "  Ndiheds = ", molec%geom%ndihed
            print'(X,A,I0)', "  Nimprop = ", molec%geom%nimprop
            print'(X,A,/)', "----------------------------------------------------------"
        endif
            

        return

    end subroutine define_internal_set


    subroutine compute_internal(molec,Ns,S, &
!                                           Optional:
                                            ICDef)
    
        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  Computes internal coordinates. Similar to internal_Wilson
        !  but without computing the B matrix
        !
        !--------------------------------------------------------------

        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp = 1.d-10 !practically zero
    
        !====================== 
        !ARGUMENTS
        type(str_resmol),intent(in)        :: molec    ! Input molecule (but only use geom...) - 
                                                       !    maybe ic structure used in fcclasses2 might be useful, setting
                                                       !    the values of the IC out of this SR
        integer,intent(in)                 :: Ns       ! Total number of internal coordiantes to use
        real(8),dimension(NDIM),intent(out)      :: S       ! Vector of internal coordinates
        character(len=100),dimension(NDIM),intent(out),optional :: ICDef !Definition of ICs
        !====================== 
    
        !======================
        !LOCAL 
        !System info
        integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s, improp_s
        integer :: nbonds, ndihed, nangles, nimprop
        integer :: Nat
        !Intenal parameters ant unitary vectors
        real(8) :: ang1, r21
        !Counters
        integer :: i,j,k
        integer :: i_1, i_2, i_3, i_4
        !=============
    
        if (verbose>0) &
            print'(/,X,A,/)', "Computing internal coordianates"
    
        !Set bonded
        nbonds  = molec%geom%nbonds
        nangles = molec%geom%nangles
        ndihed  = molec%geom%ndihed
        nimprop = molec%geom%nimprop
        bond_s(1:nbonds,1:2)  =  molec%geom%bond(1:nbonds,1:2)
        angle_s(1:nangles,1:3) =  molec%geom%angle(1:nangles,1:3)
        dihed_s(1:ndihed,1:4)  =  molec%geom%dihed(1:ndihed,1:4)
        improp_s(1:nimprop,1:4)  =  molec%geom%improp(1:nimprop,1:4)    

        !Initialize matrices
        Nat = molec%natoms
    
        !k-index runs over internal coordinates
        k=0
        if (verbose>0) &
        write(6,'(/,A,I3)') "BONDS", nbonds 
        do i=1,nbonds
            k=k+1

            i_1 = bond_s(i,1)
            i_2 = bond_s(i,2)
            r21 = calc_atm_dist(molec%atom(i_1),molec%atom(i_2))
            S(k) = r21
            if (verbose>0) &
            write(6,'(I5,X,A,2(I3,A),2F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                  trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                                               r21*BOHRtoAMS, r21
            if (present(ICDef)) &
            write(ICDef(k),'(A,2(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")"
        enddo
    
        if (verbose>0) &
        write(6,'(/,A,I3)') "ANGLES", nangles
        do i=1,nangles
            k=k+1
    
            i_1 = angle_s(i,1)
            i_3 = angle_s(i,2)
            i_2 = angle_s(i,3)
            ang1 = calc_atm_angle(molec%atom(i_1),molec%atom(i_3),molec%atom(i_2))
            S(k) = ang1
            if (verbose>0) &
            write(6,'(I5,X,A,3(I3,A),F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                 trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                                 trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                                              ang1*360.d0/2.d0/pi
            if (present(ICDef)) &
            write(ICDef(k),'(A,3(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")"
        enddo
    
        if (verbose>0) &
        write(6,'(/,A,I3)') "DIHEDRALS", ndihed
        do i=1,ndihed
            k=k+1
    
            i_1 = dihed_s(i,1)
            i_2 = dihed_s(i,2)
            i_3 = dihed_s(i,3)
            i_4 = dihed_s(i,4)
            ang1 = calc_atm_dihed_new(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
            S(k) = ang1
            if (verbose>0) &
            write(6,'(I5,X,A,4(I3,A),F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                 trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                                                 trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                                 trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")", &
                                              ang1*360.d0/2.d0/pi
            if (present(ICDef)) &
            write(ICDef(k),'(A,4(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                                          trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                          trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")"
        enddo
    
        if (verbose>0) &
        write(6,'(/,A,I3)') "IMPROPERS", nimprop
        do i=1,nimprop
            k=k+1
    
            i_1 = improp_s(i,1)
            i_2 = improp_s(i,2)
            i_3 = improp_s(i,3)
            i_4 = improp_s(i,4)
            ang1 = calc_atm_improper(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
            S(k) = ang1
            if (verbose>0) &
            write(6,'(I5,X,A,4(I3,A),F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                 trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                                                 trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                                 trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")", &
                                              ang1*360.d0/2.d0/pi
            if (present(ICDef)) &
            write(ICDef(k),'(A,4(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                                          trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                          trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")"
    
        enddo

        return
    
    end subroutine compute_internal


    subroutine internals_mapping(geom1,geom2,Map)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        ! Map geom2 elements into geom1 elements
        ! For zmat mapping: 
        !    zmatgeom = geom2
        !    allgeom  = geom1
        !
        ! Normally, geom2 should be a subgroup of geom1
        ! Otherwise, zeroes will arise in the Map 
        ! for non-coincident atoms
        !  
        !--------------------------------------------------------------

        use structure_types

        integer,parameter :: NDIM = 600

        type(str_bonded),intent(in)         :: geom1, geom2
        integer,dimension(NDIM),intent(out) :: Map

        ! Local
        integer :: i,j, k1, k2, Ns2

        print*, "Mapping internal sets..."

        !Initialize Zmap
        Ns2=geom2%nbonds+geom2%nangles+geom2%ndihed
        Map(1:Ns2) = 0

        ! Map(geom2) = geom1
        do j=1,geom2%nbonds  !Nat-1  ! Zmat loop 
        do i=1,geom1%nbonds ! Redundant loop
            if (geom2%bond(j,1)==geom1%bond(i,1).and.&
                geom2%bond(j,2)==geom1%bond(i,2)) then
                Map(j) = i
            endif
            if (geom2%bond(j,2)==geom1%bond(i,1).and.&
                geom2%bond(j,1)==geom1%bond(i,2)) then
                Map(j) = i
            endif
        enddo
        enddo
        k1=geom1%nbonds
        k2=geom2%nbonds
        do j=1,geom2%nangles !Nat-2   ! Zmat loop 
        do i=1,geom1%nangles ! Redundant loop
            if (geom2%angle(j,1)==geom1%angle(i,1).and.&
                geom2%angle(j,2)==geom1%angle(i,2).and.&
                geom2%angle(j,3)==geom1%angle(i,3)) then
                Map(k2+j) = k1+i
            endif
            if (geom2%angle(j,3)==geom1%angle(i,1).and.&
                geom2%angle(j,2)==geom1%angle(i,2).and.&
                geom2%angle(j,1)==geom1%angle(i,3)) then
                Map(k2+j) = k1+i
            endif
        enddo
        enddo
        k1=k1+geom1%nangles
        k2=k2+geom2%nangles
        do j=1,geom2%ndihed !Nat-3  ! Zmat loop 
        do i=1,geom1%ndihed ! Redundant loop
            if (geom2%dihed(j,1)==geom1%dihed(i,1).and.&
                geom2%dihed(j,2)==geom1%dihed(i,2).and.&
                geom2%dihed(j,3)==geom1%dihed(i,3).and.&
                geom2%dihed(j,4)==geom1%dihed(i,4)) then
                Map(k2+j) = k1+i
            endif
            if (geom2%dihed(j,4)==geom1%dihed(i,1).and.&
                geom2%dihed(j,3)==geom1%dihed(i,2).and.&
                geom2%dihed(j,2)==geom1%dihed(i,3).and.&
                geom2%dihed(j,1)==geom1%dihed(i,4)) then
                Map(k2+j) = k1+i
            endif
        enddo
        enddo

        return

    end subroutine internals_mapping

    function map_Zmatrix(Nvib,S,Zmap,S0) result(Sz)

        ! Get S(Zmatrix) from S(redudant) using the Zmap associations
        ! and S0 to get frozen values (from rmzfile cases)
        ! S0 has the shape of the output Sz

        integer,parameter :: NDIM = 600

        integer,intent(in)              :: Nvib
        real(8),dimension(:),intent(in) :: S
        integer,dimension(:),intent(in) :: Zmap
        real(8),dimension(:),intent(in),optional :: S0
        ! Output:                       
        real(8),dimension(Nvib)         :: Sz

        ! Local
        integer :: i

        do i=1,Nvib
            if (Zmap(i) == 0) then
                if (present(S0)) then
                    Sz(i) = S0(i)
                else
                    call alert_msg("fatal","Wrong Zmap index")
                endif
            else
                Sz(i) = S(Zmap(i))
            endif
        enddo

        return

    end function map_Zmatrix


!=========================================================
! Wilson B matrix elements and their derivatives
!=========================================================

    function Bstre(X1,Y1,Z1,X2,Y2,Z2) result(Bi)
    
        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  STRETCHING:
        !    (1)-(2)
        !
        !   Ref: Decius, Cross and Wilson (Section 4.2) --
        !    The same nomenclature is used (index refer to the same
        !    atoms, even if reversed)
        !--------------------------------------------------------------

        use line_preprocess
        use alerts
        use constants
        use metrics
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp = 1.d-10 !practically zero

        ! ARGUMENTS
        real(8),intent(in)   :: X1,Y1,Z1, &
                                X2,Y2,Z2
        real(8),dimension(6) :: Bi
        
    
        !======================
        !LOCAL 
        !AUXILIAR MATRICES
        !Intenal parameters ant unitary vectors
        real(8) :: ang1, ang2, ang3, r21, r31, r32, r43
        real(8) :: e21x, e21y, e21z,&
                   e31x, e31y, e31z,&
                   e32x, e32y, e32z,&
                   e43x, e43y, e43z
        real(8) :: e21Pe32x,e21Pe32y,e21Pe32z,&
                   e43Pe32x,e43Pe32y,e43Pe32z,&
                   e32Pe21Pe32x,e32Pe21Pe32y,e32Pe21Pe32z,&
                   e32Pe43Pe32x,e32Pe43Pe32y,e32Pe43Pe32z
        !Counters
        integer :: i,j,k
        !=============
    

        r21 = calc_dist(X1,Y1,Z1,&
                        X2,Y2,Z2)
    
        !Two cart displacements different from zero.
        e21x = (X1-X2)/r21
        e21y = (Y1-Y2)/r21
        e21z = (Z1-Z2)/r21
        ! s1 = e21 = -e12
        Bi(1) = e21x
        Bi(2) = e21y
        Bi(3) = e21z
        ! s2 = e12 = -e21
        Bi(4) =-e21x
        Bi(5) =-e21y
        Bi(6) =-e21z

        return

    end function Bstre

    subroutine dernumBstre(derB,X1,Y1,Z1,X2,Y2,Z2)

        ! ARGUMENTS
        real(8),dimension(6,6),intent(out) :: derB
        real(8),intent(in)                 :: X1,Y1,Z1,X2,Y2,Z2

        !Local
        real(8),parameter :: delta=1.889726133d-3 !for numerical ders, in bohr(=10^-3 \AA, as Num freq in G09)


        ! X1
        derB(1,1:6) = (Bstre(X1+delta,Y1,Z1,X2,Y2,Z2) -  &
                       Bstre(X1-delta,Y1,Z1,X2,Y2,Z2) )/ & 
                      (2.d0*delta)
        ! Y1
        derB(2,1:6) = (Bstre(X1,Y1+delta,Z1,X2,Y2,Z2) -  &
                       Bstre(X1,Y1-delta,Z1,X2,Y2,Z2) )/ & 
                      (2.d0*delta)
        ! Z1
        derB(3,1:6) = (Bstre(X1,Y1,Z1+delta,X2,Y2,Z2) -  &
                       Bstre(X1,Y1,Z1-delta,X2,Y2,Z2) )/ & 
                      (2.d0*delta)
        ! X2
        derB(4,1:6) = (Bstre(X1,Y1,Z1,X2+delta,Y2,Z2) -  &
                       Bstre(X1,Y1,Z1,X2-delta,Y2,Z2) )/ & 
                      (2.d0*delta)
        ! Y2
        derB(5,1:6) = (Bstre(X1,Y1,Z1,X2,Y2+delta,Z2) -  &
                       Bstre(X1,Y1,Z1,X2,Y2-delta,Z2) )/ & 
                      (2.d0*delta)
        ! Z2
        derB(6,1:6) = (Bstre(X1,Y1,Z1,X2,Y2,Z2+delta) -  &
                       Bstre(X1,Y1,Z1,X2,Y2,Z2-delta) )/ & 
                      (2.d0*delta)

        return

    end subroutine dernumBstre

    function Bbend(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3) result(Bi)
    
        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  ANGLE:
        !    (1)    (2)
        !       \   /
        !        (3)
        !
        !   Ref: Decius, Cross and Wilson (Section 4.2) --
        !    The same nomenclature is used (index refer to the same
        !    atoms, even if reversed)
        !
        ! NOTE: on input, the arguments are ordered as: 
        !       (XYZ)_1, (XYZ)_3, (XYZ)_2
        !       i.e., following a linear order in terms of connectivity
        !       But in the subroutine, we use the same order as in the book
        !--------------------------------------------------------------

        use line_preprocess
        use alerts
        use constants
        use metrics
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp = 1.d-10 !practically zero

        ! ARGUMENTS
        real(8),intent(in)   :: X1,Y1,Z1, &
                                X2,Y2,Z2, &
                                X3,Y3,Z3
        real(8),dimension(9) :: Bi
        
    
        !======================
        !LOCAL 
        !AUXILIAR MATRICES
        !Intenal parameters ant unitary vectors
        real(8) :: ang1, r21, r31, r32
        real(8) :: e21x, e21y, e21z,&
                   e31x, e31y, e31z,&
                   e32x, e32y, e32z
        !Counters
        integer :: i,j,k
        !=============


        ang1 = calc_angle(X1,Y1,Z1,X3,Y3,Z3,X2,Y2,Z2)
! print*, X2, ang1*180.d0/PI
    
        !Three cart displacements different from zero.
        r31=calc_dist(X1,Y1,Z1,X3,Y3,Z3)
        e31x = (X1-X3)/r31
        e31y = (Y1-Y3)/r31
        e31z = (Z1-Z3)/r31
        r32=calc_dist(X2,Y2,Z2,X3,Y3,Z3)
        e32x = (X2-X3)/r32
        e32y = (Y2-Y3)/r32
        e32z = (Z2-Z3)/r32
        ! s1 = [ cos(ang1)*e31 - e32 ] / [ r31 sin(ang1)]
        Bi(1) = (dcos(ang1)*e31x-e32x)/(r31*dsin(ang1))
! print*, "B1", Bi(1)
        Bi(2) = (dcos(ang1)*e31y-e32y)/(r31*dsin(ang1))
        Bi(3) = (dcos(ang1)*e31z-e32z)/(r31*dsin(ang1))
        ! s2 = [ cos(ang1)*e32 - e31 ] / [ r32 sin(ang1)]
        Bi(4) = (dcos(ang1)*e32x-e31x)/(r32*dsin(ang1))
        Bi(5) = (dcos(ang1)*e32y-e31y)/(r32*dsin(ang1))
        Bi(6) = (dcos(ang1)*e32z-e31z)/(r32*dsin(ang1))
        ! s3 = [(r31-r32 cos(ang1))e31 + (r32-r31 cos(ang1))e32 / [ rr3132 sin(ang1)]
        Bi(7) = ( (r31-r32*dcos(ang1))*e31x + (r32-r31*dcos(ang1))*e32x )&
                 / ( r31*r32*dsin(ang1))
        Bi(8) = ( (r31-r32*dcos(ang1))*e31y + (r32-r31*dcos(ang1))*e32y )&
                 / ( r31*r32*dsin(ang1))
        Bi(9) = ( (r31-r32*dcos(ang1))*e31z + (r32-r31*dcos(ang1))*e32z )&
                 / ( r31*r32*dsin(ang1))

        return
    
    end function Bbend

    subroutine dernumBbend(derB,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)

        ! Note that here, the arguments (XYZ)_i are in order (i=1,2,3)

        ! ARGUMENTS
        real(8),dimension(9,9),intent(out) :: derB
        real(8),intent(in)                 :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3

        !Local
        real(8),parameter :: delta=1.889726133d-3 !for numerical ders, in bohr(=10^-3 \AA, as Num freq in G09)


        ! X1
        derB(1,1:9) = (Bbend(X1+delta,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3) -  &
                       Bbend(X1-delta,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3) )/ & 
                      (2.d0*delta)
        ! Y1
        derB(2,1:9) = (Bbend(X1,Y1+delta,Z1,X2,Y2,Z2,X3,Y3,Z3) -  &
                       Bbend(X1,Y1-delta,Z1,X2,Y2,Z2,X3,Y3,Z3) )/ & 
                      (2.d0*delta)
        ! Z1
        derB(3,1:9) = (Bbend(X1,Y1,Z1+delta,X2,Y2,Z2,X3,Y3,Z3) -  &
                       Bbend(X1,Y1,Z1-delta,X2,Y2,Z2,X3,Y3,Z3) )/ & 
                      (2.d0*delta)
        ! X2
        derB(4,1:9) = (Bbend(X1,Y1,Z1,X2+delta,Y2,Z2,X3,Y3,Z3) -  &
                       Bbend(X1,Y1,Z1,X2-delta,Y2,Z2,X3,Y3,Z3) )/ & 
                      (2.d0*delta)
        ! Y2
        derB(5,1:9) = (Bbend(X1,Y1,Z1,X2,Y2+delta,Z2,X3,Y3,Z3) -  &
                       Bbend(X1,Y1,Z1,X2,Y2-delta,Z2,X3,Y3,Z3) )/ & 
                      (2.d0*delta)
        ! Z2
        derB(6,1:9) = (Bbend(X1,Y1,Z1,X2,Y2,Z2+delta,X3,Y3,Z3) -  &
                       Bbend(X1,Y1,Z1,X2,Y2,Z2-delta,X3,Y3,Z3) )/ & 
                      (2.d0*delta)
        ! X3
        derB(7,1:9) = (Bbend(X1,Y1,Z1,X2,Y2,Z2,X3+delta,Y3,Z3) -  &
                       Bbend(X1,Y1,Z1,X2,Y2,Z2,X3-delta,Y3,Z3) )/ & 
                      (2.d0*delta)
        ! Y3
        derB(8,1:9) = (Bbend(X1,Y1,Z1,X2,Y2,Z2,X3,Y3+delta,Z3) -  &
                       Bbend(X1,Y1,Z1,X2,Y2,Z2,X3,Y3-delta,Z3) )/ & 
                      (2.d0*delta)
        ! Z3
        derB(9,1:9) = (Bbend(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3+delta) -  &
                       Bbend(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3-delta) )/ & 
                      (2.d0*delta)

        return

    end subroutine dernumBbend

    function Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) result(Bi)
    
        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  DIHEDRAL:
        !    (1)       (4)
        !      \       /
        !       (2)-(3)
        !
        !   Ref: Decius, Cross and Wilson (Section 4.2) --
        !    The same nomenclature is used (index refer to the same
        !    atoms, even if reversed)
        !--------------------------------------------------------------

        use line_preprocess
        use alerts
        use constants
        use metrics
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp = 1.d-10 !practically zero

        ! ARGUMENTS
        real(8),intent(in)   :: X1,Y1,Z1, &
                                X2,Y2,Z2, &
                                X3,Y3,Z3, &
                                X4,Y4,Z4
        real(8),dimension(12):: Bi
        
    
        !======================
        !LOCAL 
        !AUXILIAR MATRICES
        !Intenal parameters ant unitary vectors
        real(8) :: ang1, ang2, ang3, r21, r31, r32, r43
        real(8) :: e21x, e21y, e21z,&
                   e31x, e31y, e31z,&
                   e32x, e32y, e32z,&
                   e43x, e43y, e43z
        real(8) :: e21Pe32x,e21Pe32y,e21Pe32z,&
                   e43Pe32x,e43Pe32y,e43Pe32z,&
                   e32Pe21Pe32x,e32Pe21Pe32y,e32Pe21Pe32z,&
                   e32Pe43Pe32x,e32Pe43Pe32y,e32Pe43Pe32z
        !Counters
        integer :: i,j,k
        !=============
   

            ang1 = calc_dihed(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4)
    
            !Four cart displacements different from zero (some index intercheged with Decius..)
            ang2 = calc_angle(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)
            ang3 = calc_angle(X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4)
            r21  = calc_dist(X1,Y1,Z1,X2,Y2,Z2)

            e21x = (X1-X2)/r21
            e21y = (Y1-Y2)/r21
            e21z = (Z1-Z2)/r21
            r32  = calc_dist(X2,Y2,Z2,X3,Y3,Z3)
            e32x = (X2-X3)/r32
            e32y = (Y2-Y3)/r32
            e32z = (Z2-Z3)/r32
            r43  = calc_dist(X3,Y3,Z3,X4,Y4,Z4)
            e43x = (X3-X4)/r43
            e43y = (Y3-Y4)/r43
            e43z = (Z3-Z4)/r43
            e21Pe32x=e21y*e32z-e21z*e32y
            e21Pe32y=e21z*e32x-e21x*e32z
            e21Pe32z=e21x*e32y-e21y*e32x
            e32Pe21Pe32x=e32y*e21Pe32z-e32z*e21Pe32y
            e32Pe21Pe32y=e32z*e21Pe32x-e32x*e21Pe32z
            e32Pe21Pe32z=e32x*e21Pe32y-e32y*e21Pe32x
            e43Pe32x=e43y*e32z-e43z*e32y
            e43Pe32y=e43z*e32x-e43x*e32z
            e43Pe32z=e43x*e32y-e43y*e32x
            e32Pe43Pe32x=e32y*e43Pe32z-e32z*e43Pe32y
            e32Pe43Pe32y=e32z*e43Pe32x-e32x*e43Pe32z
            e32Pe43Pe32z=e32x*e43Pe32y-e32y*e43Pe32x
    
            !s1
            Bi(1)  =  -e21Pe32x/(r21*dsin(ang2)**2)
            Bi(2)  =  -e21Pe32y/(r21*dsin(ang2)**2)
            Bi(3)  =  -e21Pe32z/(r21*dsin(ang2)**2)
            !s2   
            Bi(4)  = ((r32-r21*dcos(ang2))*e21Pe32x/(r32*r21*dsin(ang2)**2) &
                   +  dcos(ang3)*e43Pe32x/(r32*dsin(ang3)**2))
            Bi(5)  = ((r32-r21*dcos(ang2))*e21Pe32y/(r32*r21*dsin(ang2)**2) &
                   +  dcos(ang3)*e43Pe32y/(r32*dsin(ang3)**2))
            Bi(6)  = ((r32-r21*dcos(ang2))*e21Pe32z/(r32*r21*dsin(ang2)**2) &
                   +  dcos(ang3)*e43Pe32z/(r32*dsin(ang3)**2))
            !s3   
            Bi(7)  = ((r32-r43*dcos(ang3))*e43Pe32x/(r32*r43*dsin(ang3)**2) &
                   +  dcos(ang2)*e21Pe32x/(r32*dsin(ang2)**2))
            Bi(8)  = ((r32-r43*dcos(ang3))*e43Pe32y/(r32*r43*dsin(ang3)**2) &
                   +  dcos(ang2)*e21Pe32y/(r32*dsin(ang2)**2))
            Bi(9)  = ((r32-r43*dcos(ang3))*e43Pe32z/(r32*r43*dsin(ang3)**2) &
                   +  dcos(ang2)*e21Pe32z/(r32*dsin(ang2)**2))
            !s4
            Bi(10) =  -e43Pe32x/(r43*dsin(ang3)**2)
            Bi(11) =  -e43Pe32y/(r43*dsin(ang3)**2)
            Bi(12) =  -e43Pe32z/(r43*dsin(ang3)**2)
    
        return
    
    end function Bdihe

    subroutine dernumBdihe(derB,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4)

        ! ARGUMENTS
        real(8),dimension(12,12),intent(out) :: derB
        real(8),intent(in)                 :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4

        !Local
        real(8),parameter :: delta=1.889726133d-3 !for numerical ders, in bohr(=10^-3 \AA, as Num freq in G09)


        ! X1
        derB(1, 1:12)= (Bdihe(X1+delta,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bdihe(X1-delta,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Y1
        derB(2, 1:12)= (Bdihe(X1,Y1+delta,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bdihe(X1,Y1-delta,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Z1
        derB(3, 1:12)= (Bdihe(X1,Y1,Z1+delta,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bdihe(X1,Y1,Z1-delta,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! X2
        derB(4, 1:12)= (Bdihe(X1,Y1,Z1,X2+delta,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bdihe(X1,Y1,Z1,X2-delta,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Y2
        derB(5, 1:12)= (Bdihe(X1,Y1,Z1,X2,Y2+delta,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bdihe(X1,Y1,Z1,X2,Y2-delta,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Z2
        derB(6, 1:12)= (Bdihe(X1,Y1,Z1,X2,Y2,Z2+delta,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bdihe(X1,Y1,Z1,X2,Y2,Z2-delta,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! X3
        derB(7, 1:12)= (Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3+delta,Y3,Z3,X4,Y4,Z4) -  &
                        Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3-delta,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Y3
        derB(8, 1:12)= (Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3+delta,Z3,X4,Y4,Z4) -  &
                        Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3-delta,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Z3
        derB(9, 1:12)= (Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3+delta,X4,Y4,Z4) -  &
                        Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3-delta,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! X4
        derB(10,1:12)= (Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4+delta,Y4,Z4) -  &
                        Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4-delta,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Y4
        derB(11,1:12)= (Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4+delta,Z4) -  &
                        Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4-delta,Z4) )/ & 
                       (2.d0*delta)
        ! Z4
        derB(12,1:12)= (Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4+delta) -  &
                        Bdihe(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4-delta) )/ & 
                       (2.d0*delta)

        return

    end subroutine dernumBdihe


    function Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) result(Bi)
    
        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  IMPROPER (bond 1-4 with plane 2-4-3:
        !           (2)
        !           /
        !  (1) - (4) 
        !           \
        !           (3)
        !
        !   Ref: Decius, Cross and Wilson (Section 4.2) --
        !    The same nomenclature is used (index refer to the same
        !    atoms, even if reversed)
        !--------------------------------------------------------------

        use line_preprocess
        use alerts
        use constants
        use metrics
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp = 1.d-10 !practically zero

        ! ARGUMENTS
        real(8),intent(in)   :: X1,Y1,Z1, &
                                X2,Y2,Z2, &
                                X3,Y3,Z3, &
                                X4,Y4,Z4
        real(8),dimension(12):: Bi
        
    
        !======================
        !LOCAL 
        !AUXILIAR MATRICES
        !Intenal parameters ant unitary vectors
        real(8) :: theta, phi1, r41, r42, r43
        real(8) :: e41x, e41y, e41z,&
                   e42x, e42y, e42z,&
                   e43x, e43y, e43z
        real(8) :: e42Pe43x,e42Pe43y,e42Pe43z,&
                   e43Pe41x,e43Pe41y,e43Pe41z,&
                   e41Pe42x,e41Pe42y,e41Pe42z
        !Counters
        integer :: i,j,k
        !=============
   

            theta = calc_improper(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4)
    
            !Four cart displacements different from zero (some index intercheged with Decius..)
            phi1 = calc_angle(X2,Y2,Z2,X4,Y4,Z4,X3,Y3,Z3)

            r41  = calc_dist(X1,Y1,Z1,X4,Y4,Z4)
            e41x = (X1-X4)/r41
            e41y = (Y1-Y4)/r41
            e41z = (Z1-Z4)/r41
            r42  = calc_dist(X2,Y2,Z2,X4,Y4,Z4)
            e42x = (X2-X4)/r42
            e42y = (Y2-Y4)/r42
            e42z = (Z2-Z4)/r42
            r43  = calc_dist(X3,Y3,Z3,X4,Y4,Z4)
            e43x = (X3-X4)/r43
            e43y = (Y3-Y4)/r43
            e43z = (Z3-Z4)/r43
            e42Pe43x=e42y*e43z-e42z*e43y
            e42Pe43y=e42z*e43x-e42x*e43z
            e42Pe43z=e42x*e43y-e42y*e43x
            e43Pe41x=e43y*e41z-e43z*e41y
            e43Pe41y=e43z*e41x-e43x*e41z
            e43Pe41z=e43x*e41y-e43y*e41x
            e41Pe42x=e41y*e42z-e41z*e42y
            e41Pe42y=e41z*e42x-e41x*e42z
            e41Pe42z=e41x*e42y-e41y*e42x    

            !s1
            Bi(1)  = (e42Pe43x/dcos(theta)/dsin(phi1) - dtan(theta)*e41x)/r41
            Bi(2)  = (e42Pe43y/dcos(theta)/dsin(phi1) - dtan(theta)*e41y)/r41
            Bi(3)  = (e42Pe43z/dcos(theta)/dsin(phi1) - dtan(theta)*e41z)/r41
            !s2   
            Bi(4)  = (e43Pe41x/dcos(theta)/dsin(phi1) - dtan(theta)/dsin(phi1)**2*(e42x-dcos(phi1)*e43x))/r42
            Bi(5)  = (e43Pe41y/dcos(theta)/dsin(phi1) - dtan(theta)/dsin(phi1)**2*(e42y-dcos(phi1)*e43y))/r42
            Bi(6)  = (e43Pe41z/dcos(theta)/dsin(phi1) - dtan(theta)/dsin(phi1)**2*(e42z-dcos(phi1)*e43z))/r42
            !s3   
            Bi(7)  = (e41Pe42x/dcos(theta)/dsin(phi1) - dtan(theta)/dsin(phi1)**2*(e43x-dcos(phi1)*e42x))/r43
            Bi(8)  = (e41Pe42y/dcos(theta)/dsin(phi1) - dtan(theta)/dsin(phi1)**2*(e43y-dcos(phi1)*e42y))/r43
            Bi(9)  = (e41Pe42z/dcos(theta)/dsin(phi1) - dtan(theta)/dsin(phi1)**2*(e43z-dcos(phi1)*e42z))/r43
            !s4
            Bi(10) =  -Bi(1) - Bi(4) - Bi(7)
            Bi(11) =  -Bi(2) - Bi(5) - Bi(8)
            Bi(12) =  -Bi(3) - Bi(6) - Bi(9)

        return
    
    end function Bimpr

    subroutine dernumBimpr(derB,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4)

        ! ARGUMENTS
        real(8),dimension(12,12),intent(out) :: derB
        real(8),intent(in)                 :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4

        !Local
        real(8),parameter :: delta=1.889726133d-3 !for numerical ders, in bohr(=10^-3 \AA, as Num freq in G09)


        ! X1
        derB(1, 1:12)= (Bimpr(X1+delta,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bimpr(X1-delta,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Y1
        derB(2, 1:12)= (Bimpr(X1,Y1+delta,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bimpr(X1,Y1-delta,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Z1
        derB(3, 1:12)= (Bimpr(X1,Y1,Z1+delta,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bimpr(X1,Y1,Z1-delta,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! X2
        derB(4, 1:12)= (Bimpr(X1,Y1,Z1,X2+delta,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bimpr(X1,Y1,Z1,X2-delta,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Y2
        derB(5, 1:12)= (Bimpr(X1,Y1,Z1,X2,Y2+delta,Z2,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bimpr(X1,Y1,Z1,X2,Y2-delta,Z2,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Z2
        derB(6, 1:12)= (Bimpr(X1,Y1,Z1,X2,Y2,Z2+delta,X3,Y3,Z3,X4,Y4,Z4) -  &
                        Bimpr(X1,Y1,Z1,X2,Y2,Z2-delta,X3,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! X3
        derB(7, 1:12)= (Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3+delta,Y3,Z3,X4,Y4,Z4) -  &
                        Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3-delta,Y3,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Y3
        derB(8, 1:12)= (Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3+delta,Z3,X4,Y4,Z4) -  &
                        Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3-delta,Z3,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Z3
        derB(9, 1:12)= (Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3+delta,X4,Y4,Z4) -  &
                        Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3-delta,X4,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! X4
        derB(10,1:12)= (Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4+delta,Y4,Z4) -  &
                        Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4-delta,Y4,Z4) )/ & 
                       (2.d0*delta)
        ! Y4
        derB(11,1:12)= (Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4+delta,Z4) -  &
                        Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4-delta,Z4) )/ & 
                       (2.d0*delta)
        ! Z4
        derB(12,1:12)= (Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4+delta) -  &
                        Bimpr(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4-delta) )/ & 
                       (2.d0*delta)

        return

    end subroutine dernumBimpr


    subroutine internal_Wilson(molec,Ns,S,B, &
!                              Optional:
                               ICDef)
    
        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  Computes Wilson B matrix. Also get computed the values of the
        !  internal coordiantes (ICs) and, optionally, gets a written
        !  description of the ICs
        !
        !   Ref: Decius, Cross and Wilson (Section 4.2) --
        !    The same nomenclature is used (index refer to the same
        !    atoms, even if reversed)
        !--------------------------------------------------------------

        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp = 1.d-10 !practically zero
    
        !====================== 
        !ARGUMENTS
        type(str_resmol),intent(inout)     :: molec    ! Input molecule (but only use geom...) - 
                                                       !    maybe ic structure used in fcclasses2 might be useful, setting
                                                       !    the values of the IC out of this SR
        integer,intent(in)                 :: Ns       ! Total number of internal coordiantes to use
        real(8),dimension(NDIM,NDIM),intent(out) :: B  ! B Wilson matrix
        real(8),dimension(NDIM),intent(out)      :: S       ! Vector of internal coordinates
        character(len=100),dimension(NDIM),intent(out),optional :: ICDef !Definition of ICs
        !====================== 
    
        !======================
        !LOCAL 
        !System info
        integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s, improp_s
        integer :: nbonds, ndihed, nangles, nimprop
        integer :: Nat
        !AUXILIAR MATRICES
        real(8),dimension(1:12) :: Baux
        !Intenal parameters ant unitary vectors
        real(8) :: ang1, r21
        !Counters
        integer :: i,j,k
        integer :: i_1, i_2, i_3, i_4
        !=============

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with Atomic Units 
        current_units=molec%units
        call set_geom_units(molec,"Bohr")
        !----------------------------------------------
    
        if (verbose>0) then
            write(6,'(/,2X,A)') "Calculating internal coordianates and Wilson B matrix..."
            write(6,'(X,A)')    "---------------------------------------------------------"
        endif
    
        !Set bonded
        nbonds  = molec%geom%nbonds
        nangles = molec%geom%nangles
        ndihed  = molec%geom%ndihed
        nimprop = molec%geom%nimprop
        bond_s(1:nbonds,1:2)  =  molec%geom%bond(1:nbonds,1:2)
        angle_s(1:nangles,1:3) =  molec%geom%angle(1:nangles,1:3)
        dihed_s(1:ndihed,1:4)  =  molec%geom%dihed(1:ndihed,1:4)
        improp_s(1:nimprop,1:4)  =  molec%geom%improp(1:nimprop,1:4)    

        !Initialize matrices
        Nat = molec%natoms
        B(1:Ns,1:3*Nat) = 0.d0
    
        !k-index runs over internal coordinates
        k=0
        if (verbose>0) &
        write(6,'(/,A,I3)') "BONDS", nbonds 
        do i=1,nbonds
            k=k+1

            i_1 = bond_s(i,1)
            i_2 = bond_s(i,2)
            r21 = calc_atm_dist(molec%atom(i_1),molec%atom(i_2))
            S(k) = r21
            if (verbose>0) &
            write(6,'(I5,X,A,2(I3,A),2F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                  trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                                               r21*BOHRtoAMS, r21
            if (present(ICDef)) &
            write(ICDef(k),'(A,2(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")"

            ! Compute with external functions
            Baux(1:6) = Bstre(molec%atom(i_1)%x,molec%atom(i_1)%y,molec%atom(i_1)%z,&
                              molec%atom(i_2)%x,molec%atom(i_2)%y,molec%atom(i_2)%z)
            B(k,3*i_1-2:3*i_1) = Baux(1:3)
            B(k,3*i_2-2:3*i_2) = Baux(4:6)
        enddo
    
        if (verbose>0) &
        write(6,'(/,A,I3)') "ANGLES", nangles
        do i=1,nangles
            k=k+1
    
            i_1 = angle_s(i,1)
            i_3 = angle_s(i,2)
            i_2 = angle_s(i,3)
            ang1 = calc_atm_angle(molec%atom(i_1),molec%atom(i_3),molec%atom(i_2))
            S(k) = ang1
            if (verbose>0) &
            write(6,'(I5,X,A,3(I3,A),F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                 trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                                 trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                                              ang1*360.d0/2.d0/pi
            if (present(ICDef)) &
            write(ICDef(k),'(A,3(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")"
    
            ! Compute with external functions
            Baux(1:9) = Bbend(molec%atom(i_1)%x,molec%atom(i_1)%y,molec%atom(i_1)%z,&
                              molec%atom(i_2)%x,molec%atom(i_2)%y,molec%atom(i_2)%z,&
                              molec%atom(i_3)%x,molec%atom(i_3)%y,molec%atom(i_3)%z)
            B(k,3*i_1-2:3*i_1) = Baux(1:3)
            B(k,3*i_2-2:3*i_2) = Baux(4:6)
            B(k,3*i_3-2:3*i_3) = Baux(7:9)
        enddo
    
        if (verbose>0) &
        write(6,'(/,A,I3)') "DIHEDRALS", ndihed
        do i=1,ndihed
            k=k+1
    
            i_1 = dihed_s(i,1)
            i_2 = dihed_s(i,2)
            i_3 = dihed_s(i,3)
            i_4 = dihed_s(i,4)
            ang1 = calc_atm_dihed_new(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
            S(k) = ang1
            if (verbose>0) &
            write(6,'(I5,X,A,4(I3,A),F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                 trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                                                 trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                                 trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")", &
                                              ang1*360.d0/2.d0/pi
            if (present(ICDef)) &
            write(ICDef(k),'(A,4(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                                          trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                          trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")"
    
            ! Compute with external functions
            Baux(1:12)= Bdihe(molec%atom(i_1)%x,molec%atom(i_1)%y,molec%atom(i_1)%z,&
                              molec%atom(i_2)%x,molec%atom(i_2)%y,molec%atom(i_2)%z,&
                              molec%atom(i_3)%x,molec%atom(i_3)%y,molec%atom(i_3)%z,&
                              molec%atom(i_4)%x,molec%atom(i_4)%y,molec%atom(i_4)%z)
            B(k,3*i_1-2:3*i_1) = Baux(1:3)
            B(k,3*i_2-2:3*i_2) = Baux(4:6)
            B(k,3*i_3-2:3*i_3) = Baux(7:9)
            B(k,3*i_4-2:3*i_4) = Baux(10:12)
    
        enddo
    
        if (verbose>0) &
        write(6,'(/,A,I3)') "IMPROPERS", nimprop
        do i=1,nimprop
            k=k+1
    
            i_1 = improp_s(i,1)
            i_2 = improp_s(i,2)
            i_3 = improp_s(i,3)
            i_4 = improp_s(i,4)
            ang1 = calc_atm_improper(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
            S(k) = ang1
            if (verbose>0) &
            write(6,'(I5,X,A,4(I3,A),F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                 trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                                                 trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                                 trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")", &
                                              ang1*360.d0/2.d0/pi
            if (present(ICDef)) &
            write(ICDef(k),'(A,4(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                                          trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                                          trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")"
    
            ! Compute with external functions
            Baux(1:12)= Bdihe(molec%atom(i_1)%x,molec%atom(i_1)%y,molec%atom(i_1)%z,&
                              molec%atom(i_2)%x,molec%atom(i_2)%y,molec%atom(i_2)%z,&
                              molec%atom(i_3)%x,molec%atom(i_3)%y,molec%atom(i_3)%z,&
                              molec%atom(i_4)%x,molec%atom(i_4)%y,molec%atom(i_4)%z)
            B(k,3*i_1-2:3*i_1) = Baux(1:3)
            B(k,3*i_2-2:3*i_2) = Baux(4:6)
            B(k,3*i_3-2:3*i_3) = Baux(7:9)
            B(k,3*i_4-2:3*i_4) = Baux(10:12)
    
        enddo

        if (verbose>0) &
            print*, ""
        if (verbose>1) &
            call MAT0(6,B,Ns,3*Nat,"B MATRIX")

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! Revert original units
        call set_geom_units(molec,adjustl(current_units))
        !----------------------------------------------

        return
    
    end subroutine internal_Wilson


    subroutine internal_Gmetric(Nat,Ns,Mass,B,G)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  Compute G metric matrix
        !------------------------------------------------------------------

        use structure_types
        use line_preprocess
        use alerts
        use constants
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
    
        !====================== 
        !ARGUMENTS
        integer,intent(in)                          :: Nat    ! Number of atoms (in)
        integer,intent(in)                          :: Ns     ! Number of internal coordinates (in)
        real(8),dimension(1:NDIM),intent(in)        :: Mass   ! Wilson matrices (in AMU)
        real(8),dimension(1:NDIM,1:NDIM),intent(in) :: B
        real(8),dimension(1:NDIM,1:NDIM),intent(out):: G
        !====================== 

        !======================
        !LOCAL 
        !Counters
        integer :: i,j,k, kk, iat
        real(8) :: mu
        !=============

        !CONCTRUCT G
        do i=1,Ns
            do j=1,Ns
                G(i,j) = 0.d0
                k=0
                do kk=1,Nat
                do iat=1,3
                    k=k+1
                    mu = 1.d0/Mass(kk)/UMAtoAU
                    G(i,j) = G(i,j) + mu*B(i,k)*B(j,k)
                enddo
                enddo
            enddo
        enddo

        if (verbose>1) &
            call MAT0(6,G*1.d5,Ns,Ns,"G MATRIX x1e5")

        return

    end subroutine internal_Gmetric



    subroutine redundant2nonredundant(Nred,Nvib,G,Asel)
    
        !IF WE USE ALL BONDED PARAMETERS,WE HAVE REDUNDANCY. WE SELECT A NON-REDUNDANT                                                                  
        !COMBINATION FROM THE NON-ZERO EIGENVALUES OF G (Reimers 2001, JCP)   


        use matrix
        use matrix_print

        implicit none

        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp=1.d-10

        integer,intent(inout) :: Nred, Nvib
        real(8),dimension(NDIM,NDIM),intent(in) :: G
        real(8),dimension(NDIM,NDIM),intent(out):: Asel

        ! Local
        real(8),dimension(NDIM,NDIM) :: Aux
        real(8),dimension(NDIM)      :: Vec, Vec2
        integer :: i, k,kk,kkk

        !Get a non-redundant set from the non-zero eigenvalues of G
        call diagonalize_full(G(1:Nred,1:Nred),Nred,Aux(1:Nred,1:Nred),Vec(1:Nred),"lapack")

        if (verbose>2) &
         call MAT0(6,Aux,Nred,Nred,"A MATRIX (before reordering)")
        if (verbose>1) &
         call print_vector(6,Vec,Nred,"A MATRIX Eigenvalues (before reordering)")
 
        ! The non-redundant set is formed by eigenvectors with non-zero eigenvalue
        kk=0 
        kkk=0
        do k=1,Nred
            if (dabs(Vec(k)) > ZEROp) then
                kk=kk+1
                !The eigenvectors are stored in columns (since D = P^-1 A P)
                Asel(1:Nred,kk) = Aux(1:Nred,k)
                Vec2(kk)        = Vec(k)
                ! .. for symmetric A matrices, either rows or columns can be selected as eigenvectors anyway
            else
                !Redudant eigenvectors
                kkk=kkk+1
                Asel(1:Nred,Nvib+kkk) = Aux(1:Nred,k)
                Vec2(Nvib+kkk)        = Vec(k)
            endif
        enddo

        if (Nred /= Nvib+kkk) then
            call sort_vec(Vec,Nred)
            call print_vector(6,Vec*1e5,Nred,"A MATRIX Eigenvalues (x10^5)")
            if (kkk > Nred-Nvib) then
                print*, "Internal-space dimension is reduced"
                print*, " Deleted modes:    ", kkk - Nred + Nvib
                print*, " Initial:", Nvib
                print*, " Reduced:", Nred-kkk
                call alert_msg("note","Redundant to non-redundant trasformation"//&
                                     " resulted in a reduced space")
                Nvib=Nred-kkk
            else 
                call alert_msg("fatal","Redundant to non-redundant trasformation failed")
            endif
        else
            if (verbose>0) then
                print'(X,A,I0,/)', "Internal-space dimension:", Nvib
            endif
        endif

        if (verbose>2) &
         call MAT0(6,Asel,Nred,Nred,"A MATRIX")
        if (verbose>1) then
            call print_vector(6,Vec2(Nvib+1:Nred)*1.d15,Nred-Nvib,"A MATRIX Zero-Eigenvalues (x10^15)")
            call print_vector(6,Vec2(1:Nvib)*1.d5,Nvib,"A MATRIX NonZero-Eigenvalues (x10^5)")
        endif

        return

    end subroutine redundant2nonredundant


    subroutine generalized_inv(Nred,Nnonred,G,Ginv)
    
        !IF WE USE ALL BONDED PARAMETERS,WE HAVE REDUNDANCY. WE SELECT A NON-REDUNDANT                                                                  
        !COMBINATION FROM THE NON-ZERO EIGENVALUES OF G (Reimers 2001, JCP) 
        ! *** Needs more checks ***
        ! Alternative: use redundant2nonredundant + matrix rotations  


        use matrix
        use matrix_print

        implicit none

        integer,parameter :: NDIM = 600
        real(8),parameter :: ZEROp=1.d-10

        integer,intent(inout) :: Nred
        integer,intent(out)   :: Nnonred
        real(8),dimension(NDIM,NDIM),intent(in) :: G
        real(8),dimension(NDIM,NDIM),intent(out):: Ginv

        ! Local
        real(8),dimension(NDIM,NDIM) :: Aux, Aux2
        real(8),dimension(NDIM)      :: Vec, Vec2, Vec3
        integer :: i, k,kk,kkk

        if (verbose>0) then
           print'(/,X,A,/)', "Generalized Inverse subroutine"
        endif

        !Get a non-redundant set from the non-zero eigenvalues of G
        call diagonalize_full(G(1:Nred,1:Nred),Nred,Aux(1:Nred,1:Nred),Vec(1:Nred),"lapack")
 
        ! Inverse the non-zero eigenvalues
        kk=0 
        kkk=0
        do k=1,Nred
            if (dabs(Vec(k)) > ZEROp) then
                kk=kk+1
                Vec2(k)        = 1.d0/Vec(k)
            else
                !Redudant eigenvectors
                kkk=kkk+1
                Vec2(k)        = 0.d0
            endif
        enddo
        Nnonred = kk
        kk=0 
        kkk=0
        do k=1,Nred
            if (dabs(Vec(k)) > ZEROp) then
                kk=kk+1
                Vec3(kk)          = Vec(k)
            else
                !Redudant eigenvectors
                kkk=kkk+1
                Vec3(Nnonred+kkk) = Vec(k)
            endif
        enddo
        
        if (verbose>2) then
            call MAT0(6,Aux,Nred,Nred,"A MATRIX")
            call print_vector(6,Vec,Nred,"A MATRIX Eigenvalues")
        endif
        if (verbose>1) then
            call print_vector(6,Vec3(Nnonred+1:Nred)*1e15,Nred-Nnonred,"Zero-eigenvalues(x1e15)")
            call print_vector(6,Vec3(1:Nnonred)*1e5,Nnonred,"NonZero-eigenvalues(x1e5)")
        endif

        if (verbose>0) then
            print*, "Original row-dimension: ", Nred
            print*, "Non-zero eigenvalues  : ", Nnonred
        endif 

        ! The diag inverse is stored (as vector of the diagonal elements) in Vec2
        ! Now rotate back to the original frame
        Ginv(1:Nred,1:Nred) = diag_basisrot(Nred,Nred,Aux,Vec2,counter=.false.)

        return

    end subroutine generalized_inv


    subroutine HessianCart2int(Nat,Ns,Hess,Mass,B,G,method)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  HESSIAN IN INTERNAL COORDINATES (JCC, 17, 49-56, by Frisch et al)
        !    Hint = G^- Bu(Hx+B'^Tg_q)u^TB^T G^-
        !   g_q is the gradient, so g_q=0 in a minimum
        !   G^- is the generalized inverse (for redundant internal) or simply the
        !   inverse for nonredundant
        !
        ! Arguments
        !  Nat    Int /Scalar    Number of atoms
        !  Ns     Int /Scalar    Number of internal coordianates
        !  Hess   Real/Matrix    Hessian in Cartesian (corrected or not)
        !  Mass   Real/Vector    Mass vector (Nat)
        !  B      Real/Matrix    B matrix
        !  G      Real/Matrix    Metric matrix
        !------------------------------------------------------------------

        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
    
        !====================== 
        !ARGUMENTS
        integer,intent(in)                          :: Nat    ! Number of atoms (in)
        integer,intent(in)                          :: Ns     ! Number of internal coordinates (in)
        real(8),dimension(1:NDIM),intent(in)        :: Mass   ! Wilson matrices (in)
        real(8),dimension(1:NDIM,1:NDIM),intent(in) :: G,B    ! Wilson matrices (in)
        real(8),dimension(1:NDIM,1:NDIM),intent(inout) :: Hess   !Hessian: cart(in)-intern(out)
        character(len=*),intent(in),optional        :: method
        !====================== 
    
        !====================== 
        !LOCAL
        integer                             :: Nvib
        character(len=5)                    :: method_local
        !Internal analysis 
        real(8),dimension(1:NDIM,1:NDIM)    :: Ginv
        !Auxiliar arrays
        real(8),dimension(1:NDIM,1:NDIM)    :: AuxT,Aux
        real(8),dimension(NDIM)             :: Vec
        !Counters
        integer :: i,j,k, ii
        !====================== 

        if (present(method)) then
            method_local = adjustl(method)
            call set_word_upper_case(method_local)
        else
            method_local="MASS"
        endif
    
        ! If Ns > Nvib, then we need to extract the non-zero eigen values 
        ! from G in order to compute the generalized G inverse
        ! (TBD) 
        ! 
        Nvib = Ns

        !Two choice with u=M^-1 or u=I
        if (adjustl(method_local) == "MASS") then
            !Inverse of G
            Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,G)
            !Compute G^-1Bu  (where u is the inverse mass matrix)
            Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
            do i=1,3*Nat
                ii = (i-1)/3+1
                Aux(1:Ns,i) = Aux(1:Ns,i)/Mass(ii)/UMAtoAU
            enddo           
        else
            !Inverse of (BB^t)
            Ginv(1:Ns,1:Ns) = matrix_product(Ns,Ns,3*Nat,B,B,tB=.true.)
            Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,Ginv)
            !Compute G^-1Bu  (where u is the inverse mass matrix)
            Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
        endif

        
        ! Hint = Aux ([~Hx]) Aux^T (this is "matrix_basisrot")
        Hess(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,3*Nat,Aux,Hess)
        Hess(1:Ns,1:Ns)    = matrix_product(Ns,3*Nat,3*Nat,Hess,Aux,tB=.true.)
    
        if (verbose>1) then
            Vec(1:Ns) = (/(Hess(i,i), i=1,Ns)/)
            call print_vector(6,Vec,Ns,"F MATRIX (diagonal)")
        endif
        if (verbose>2) &
         call MAT0(6,Hess,Ns,Ns,"F MATRIX")

        return

    end subroutine HessianCart2int


    subroutine HessianCart2int_red(Nat,Ns,Hess,Mass,B,G,method)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  HESSIAN IN INTERNAL COORDINATES (JCC, 17, 49-56, by Frisch et al)
        !    Hint = G^- Bu(Hx+B'^Tg_q)u^TB^T G^-
        !   g_q is the gradient, so g_q=0 in a minimum
        !   G^- is the generalized inverse (for redundant internal) or simply the
        !   inverse for nonredundant
        !
        ! Arguments
        !  Nat    Int /Scalar    Number of atoms
        !  Ns     Int /Scalar    Number of internal coordianates
        !  Hess   Real/Matrix    Hessian in Cartesian (corrected or not)
        !  Mass   Real/Vector    Mass vector (Nat)
        !  B      Real/Matrix    B matrix
        !  G      Real/Matrix    Metric matrix
        !
        ! _red for redundant
        !------------------------------------------------------------------

        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
    
        !====================== 
        !ARGUMENTS
        integer,intent(in)                          :: Nat    ! Number of atoms (in)
        integer,intent(in)                          :: Ns     ! Number of internal coordinates (in)
        real(8),dimension(1:NDIM),intent(in)        :: Mass   ! Wilson matrices (in)
        real(8),dimension(1:NDIM,1:NDIM),intent(in) :: G,B    ! Wilson matrices (in)
        real(8),dimension(1:NDIM,1:NDIM),intent(inout) :: Hess   !Hessian: cart(in)-intern(out)
        character(len=*),intent(in),optional        :: method
        !====================== 
    
        !====================== 
        !LOCAL
        integer                             :: Nvib
        character(len=5)                    :: method_local
        !Internal analysis 
        real(8),dimension(1:NDIM,1:NDIM)    :: Ginv
        !Auxiliar arrays
        real(8),dimension(1:NDIM,1:NDIM)    :: AuxT,Aux
        real(8),dimension(NDIM)             :: Vec
        !Counters
        integer :: i,j,k, ii, NN
        !====================== 

        if (present(method)) then
            method_local = adjustl(method)
            call set_word_upper_case(method_local)
        else
            method_local="MASS"
        endif
    
        ! If Ns > Nvib, then we need to extract the non-zero eigen values 
        ! from G in order to compute the generalized G inverse
        ! (TBD) 
        ! 
!         Nvib = Ns

        !Two choice with u=M^-1 or u=I
        if (adjustl(method_local) == "MASS") then
            !Inverse of G
!             Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,G)
            NN=Ns
            call generalized_inv(NN,Nvib,G,Ginv)
            !Compute G^-1Bu  (where u is the inverse mass matrix)
            Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
            do i=1,3*Nat
                ii = (i-1)/3+1
                Aux(1:Ns,i) = Aux(1:Ns,i)/Mass(ii)/UMAtoAU
            enddo           
        else
            !Inverse of (BB^t)
            Ginv(1:Ns,1:Ns) = matrix_product(Ns,Ns,3*Nat,B,B,tB=.true.)
!             Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,Ginv)
            NN=Ns
            call generalized_inv(NN,Nvib,Ginv,Ginv)
            !Compute G^-1Bu  (where u is the inverse mass matrix)
            Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
        endif

        
        ! Hint = Aux ([~Hx]) Aux^T (this is "matrix_basisrot")
        Hess(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,3*Nat,Aux,Hess)
        Hess(1:Ns,1:Ns)    = matrix_product(Ns,3*Nat,3*Nat,Hess,Aux,tB=.true.)
    
        if (verbose>1) then
            Vec(1:Ns) = (/(Hess(i,i), i=1,Ns)/)
            call print_vector(6,Vec,Ns,"F MATRIX (diagonal)")
        endif
        if (verbose>2) &
         call MAT0(6,Hess,Ns,Ns,"F MATRIX")

        return

    end subroutine HessianCart2int_red


    subroutine Gradcart2int(Nat,Ns,Grad,Mass,B,G,method)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        ! Convert the Gradient from cartesian to internal 
        ! useful when the Gradient is not used in HessianCart2int
        ! (if used there, it is converted there!)
        !------------------------------------------------------------------

        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
    
        !====================== 
        !ARGUMENTS
        integer,intent(in)                          :: Nat    ! Number of atoms (in)
        integer,intent(in)                          :: Ns     ! Number of internal coordinates (in)
        real(8),dimension(1:NDIM),intent(in)        :: Mass   ! Wilson matrices (in)
        real(8),dimension(1:NDIM,1:NDIM),intent(in) :: G,B    ! Wilson matrices (in)
        real(8),dimension(1:NDIM),intent(inout)     :: Grad   ! Gradient 
        character(len=*),intent(in),optional        :: method
        !====================== 
    
        !====================== 
        !LOCAL
        integer                             :: Nvib
        character(len=5)                    :: method_local
        !Internal analysis 
        real(8),dimension(1:NDIM,1:NDIM)    :: Ginv
        !Auxiliar arrays
        real(8),dimension(1:NDIM,1:NDIM)    :: Aux
        real(8),dimension(NDIM)             :: Vec
        !Counters
        integer :: i,j,k, ii
        !====================== 

        if (present(method)) then
            method_local = adjustl(method)
            call set_word_upper_case(method_local)
        else
            method_local="MASS"
        endif

        ! If Ns > Nvib, then we need to extract the non-zero eigen values 
        ! from G in order to compute the generalized G inverse
        ! (TBD) 
        ! 
        Nvib = Ns

        !Two choice with u=M^-1 or u=I
        if (adjustl(method_local) == "MASS") then
            !Inverse of G
            Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,G)
            !Compute G^-1Bu  (where u is the inverse mass matrix)
            Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
            do i=1,3*Nat
                ii = (i-1)/3+1
                Aux(1:Ns,i) = Aux(1:Ns,i)/Mass(ii)/UMAtoAU
            enddo           
        else
            !Inverse of (BB^t)
            Ginv(1:Ns,1:Ns) = matrix_product(Ns,Ns,3*Nat,B,B,tB=.true.)
            Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,Ginv)
            !Compute G^-1Bu  (where u is the inverse mass matrix)
            Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
        endif

        ! Get the gradient in internal coords: gq = G^-1Bu(gx)
        do i=1,Ns
            Vec(i) = 0.d0
            do j=1,3*Nat
                Vec(i) = Vec(i) + Aux(i,j) * Grad(j)
            enddo
        enddo
        ! Update the gradient on output
        Grad(1:3*Nat) = 0.d0
        Grad(1:Ns) = Vec(1:Ns)

        return

    end subroutine Gradcart2int

    subroutine Gradcart2int_red(Nat,Ns,Grad,Mass,B,G,method)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        ! Convert the Gradient from cartesian to internal 
        ! useful when the Gradient is not used in HessianCart2int
        ! (if used there, it is converted there!)
        !
        ! _red for redundant
        !------------------------------------------------------------------

        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
    
        !====================== 
        !ARGUMENTS
        integer,intent(in)                          :: Nat    ! Number of atoms (in)
        integer,intent(in)                          :: Ns     ! Number of internal coordinates (in)
        real(8),dimension(1:NDIM),intent(in)        :: Mass   ! Wilson matrices (in)
        real(8),dimension(1:NDIM,1:NDIM),intent(in) :: G,B    ! Wilson matrices (in)
        real(8),dimension(1:NDIM),intent(inout)     :: Grad   ! Gradient 
        character(len=*),intent(in),optional        :: method
        !====================== 
    
        !====================== 
        !LOCAL
        integer                             :: Nvib
        character(len=5)                    :: method_local
        !Internal analysis 
        real(8),dimension(1:NDIM,1:NDIM)    :: Ginv
        !Auxiliar arrays
        real(8),dimension(1:NDIM,1:NDIM)    :: Aux
        real(8),dimension(NDIM)             :: Vec
        !Counters
        integer :: i,j,k, ii, NN
        !====================== 

        if (present(method)) then
            method_local = adjustl(method)
            call set_word_upper_case(method_local)
        else
            method_local="MASS"
        endif

        !Two choice with u=M^-1 or u=I
        if (adjustl(method_local) == "MASS") then
            !Inverse of G
!             Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,G)
            NN=Ns
            call generalized_inv(NN,Nvib,G,Ginv)
            !Compute G^-1Bu  (where u is the inverse mass matrix)
            Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
            do i=1,3*Nat
                ii = (i-1)/3+1
                Aux(1:Ns,i) = Aux(1:Ns,i)/Mass(ii)/UMAtoAU
            enddo           
        else
            !Inverse of (BB^t)
            Ginv(1:Ns,1:Ns) = matrix_product(Ns,Ns,3*Nat,B,B,tB=.true.)
!             Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,Ginv)
            NN=Ns
            call generalized_inv(NN,Nvib,Ginv,Ginv)
            !Compute G^-1Bu  (where u is the inverse mass matrix)
            Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
        endif

        ! Get the gradient in internal coords: gq = G^-1Bu(gx)
        do i=1,Ns
            Vec(i) = 0.d0
            do j=1,3*Nat
                Vec(i) = Vec(i) + Aux(i,j) * Grad(j)
            enddo
        enddo
        ! Update the gradient on output
        Grad(1:3*Nat) = 0.d0
        Grad(1:Ns) = Vec(1:Ns)

        return

    end subroutine Gradcart2int_red

    
    subroutine gf_method(Nvib,G,Hess,L,Freq,X,Xinv)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        !FG PROCEDURE (Wilson, Decius and Cross, Section 4.3)
        ! |GF - \lamda| = 0  <==> GFL = L\lamda   means the diagonalize of GF, but is not symmetric
        ! Solutions
        ! *Orthogonalization of FL = G^-1L\lamda similar to orthogonalization fo S in Roothan. 
        ! *Symmetrize GF and diagonalize following Mizayawa (Sect. 4.6, Califano, Vibrational States)
        !
        ! Here: FOLLOWING FIRST ALTERNATIVE (similar to diagonalization of basis set)
        !------------------------------------------------------------------
    
        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
    
        !====================== 
        !ARGUMENTS
        !Input
        integer,intent(in)                          :: Nvib      !Number of vib nm
        real(8),dimension(1:NDIM,1:NDIM),intent(in) :: G         !Wilson metric matrix
        real(8),dimension(1:NDIM,1:NDIM),intent(in) :: Hess      !Hessian: internal
        !Output
        real(8),dimension(NDIM,NDIM),intent(out)    :: L,X,Xinv  !Normal modes, G^1/2, G^-1/2 (internal)
        real(8),dimension(NDIM),intent(out)         :: Freq      !Frequencies
        !====================== 
    
        !====================== 
        !LOCAL
        !Auxiliar matrices
        real(8),dimension(1:NDIM,1:NDIM) :: Aux,Aux3
        !Counters
        integer :: i,j,k, ii,jj,kk, iat
        !=============
    
    
        call diagonalize_full(G(1:Nvib,1:Nvib),Nvib,Aux(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
    
        !The non-unitary rotation which orthogonalize G^-1 could be X=G^(1/2)  (inverse compared with the case of Roothan eqs)
        ! where G^{1/2} = Ug^{1/2}U^T
        !Store the rotation (X) in AuxT    !G. Note X=X^T
        do i=1,Nvib
            do j=1,Nvib
                X(i,j) = 0.d0
                do k=1,Nvib
                    X(i,j) = X(i,j) + Aux(i,k)*dsqrt(Freq(k))*Aux(j,k)
                enddo
            enddo
        enddo
        !And the inverse
        do i=1,Nvib
            do j=1,Nvib
                Xinv(i,j) = 0.d0
                do k=1,Nvib
                    Xinv(i,j) = Xinv(i,j) + Aux(i,k)/dsqrt(Freq(k))*Aux(j,k)
                enddo
            enddo
        enddo
    
        if (verbose>1) &
            call MAT0(6,X,Nvib,Nvib,"X (G^1/2)")
        if (verbose>1) &
            call MAT0(6,Xinv,Nvib,Nvib,"Xinv (G^-1/2)")

    
        !Now rotate F, F'=X^TFX. Store the rotated matrix in Aux3 (temporary array)
        Aux3(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Nvib,X,Hess,counter=.true.)
    
        !We can now diagonalize F'
        call diagonalize_full(Aux3(1:Nvib,1:Nvib),Nvib,L(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
    
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

        !Normal mode matrix
        if (verbose>1) &
            call MAT0(6,L,Nvib,Nvib,"L' (ORTH NORMAL MODES)")
    
        if (verbose>2) then
            Aux(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,L,L,tB=.true.)
            call MAT0(6,L,Nvib,Nvib,"L'L'^t")
        endif
    
        ! The L' matrix is given in the an orthogonalized basis. 
        ! To have a description based in the initial internal, restore in the original basis:
        ! L = XL'   (with X the orhogonalizing matrix, stored in X).
        L(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,X,L)
    
        if (verbose>1) &
            call MAT0(6,L,Nvib,Nvib,"L (NON-ORTH NORMAL MODES)")
    
        ! Basic checks on the normal mode matrix
        if (verbose>2) then
            !L L^t
            Aux3(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,L,L,tB=.true.)
            call MAT0(6,Aux3,Nvib,Nvib,"L L^t")
            !L^t L
            Aux3(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,L,L,tA=.true.)
            call MAT0(6,Aux3,Nvib,Nvib,"L^t L")
            !L^tFL (should provide the FC, in this case Freq)
            Aux3(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Nvib,L,Hess,counter=.true.)
            Aux3(1:Nvib,1:Nvib) = dsqrt(dabs(Aux3(1:Nvib,1:Nvib))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2
            call MAT0(6,Aux3,Nvib,Nvib,"L^t L (cm-1)")
            !L^tG^-1L (this test only works if G has dimension Nvib x Nvib
            Aux3(1:Nvib,1:Nvib) = inverse_realsym(Nvib,G)
            Aux3(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Nvib,L,Aux3,counter=.true.)
            call MAT0(6,Aux3,Nvib,Nvib,"L^tG^-1L")
        endif
    
        return
    
    end subroutine gf_method


    subroutine analyze_internal(Nvib,Ns,L,Freq,ICDef,Ssym)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  Analysis of normal modes in terms of the contributing internal
        !  coordinates
        !------------------------------------------------------------------
    
        use structure_types
        use line_preprocess
        use alerts
        use constants
        use atomic_geom
        use matrix
        use verbosity
    
        implicit none
    
        integer,parameter :: NDIM = 600
    
        !====================== 
        !ARGUMENTS
        !Input
        integer,intent(in)                            :: Nvib, Ns
        real(8),dimension(NDIM,NDIM),intent(in)       :: L       !Normal modes
        real(8),dimension(NDIM),intent(in)            :: Freq    !Frequencies
        character(len=100),dimension(NDIM),intent(in) :: ICDef   !Definition of ICs
        integer,dimension(NDIM),intent(in),optional   :: Ssym    !array with symmetric ICs
        !====================== 
    
        !====================== 
        !LOCAL
        !System variables
        integer :: Nat
        character(len=10) :: ModeSymm
        !Auxiliar matrices
        real(8),dimension(1:NDIM,1:NDIM) :: Aux,Aux3,AuxT
        integer,dimension(1:NDIM)        :: ipiv, ipiv2
        !Auxiliar variables
        real(8) :: Theta, Theta2, mu
        !Counters
        integer :: i,j,k, ii,jj,kk
        !=============

        ! Guess Nat from Nvib (will not always work...)
        Nat = (Nvib+6)/3

        print*, ""
        print*, "L MATRIX (LARGER ELEMENTS) - by columns"
        do i=1,Nvib
            !Copy the row, normalize and reorder 
            Aux3(1,1:Ns) = abs(L(1:Ns,i))
            Theta  = 0.d0
            Theta2 = 0.d0
            do ii=1,Ns
                !Displacement-weighted elements
                if (ii <= Nat-1) then
                    Aux(1,ii) = Aux3(1,ii)
                elseif (ii <= 2*Nat-3) then
                    Aux(1,ii) = Aux3(1,ii)*2.3d0
                elseif (ii <= 3*Nat-6) then
                    Aux(1,ii) = Aux3(1,ii)*0.8d0
                endif 
                !Standard contribution calc
                Theta = Theta + Aux3(1,ii)!**2
                !Displacement-weighted calc (only valid for Zmat)
                if (Nvib==Ns) &
                 Theta2 = Theta2 + Aux(1,ii)
            enddo
            Aux3(1,1:Ns) = Aux3(1,1:Ns)/Theta
            call sort_vec_max(Aux3(1,1:Ns),ipiv(1:Ns),Ns)
            if (Nvib==Ns) then
                Aux(1,1:Ns)  = Aux(1,1:Ns)/Theta2
                call sort_vec_max(Aux(1,1:Ns),ipiv2(1:Ns),Ns)
            else
                Aux(1,1:Ns) = 0.d0
            endif
            !
            if (present(Ssym)) then
                !Determine symmetry
                do j=1,Ns
                    jj=ipiv(j)
                    if (Ssym(jj) == jj) cycle
                    ii = Ssym(jj)
                    Theta = L(jj,i)/L(ii,i)
                    if (Theta > 0.d0) then 
                        ModeSymm="Symm: A' "
                    else
                        ModeSymm="Symm: A''"
                    endif
                    exit
                    !Use stretching coordinates to set symmetry
                    if (ii<Nat) exit
                enddo
            else
                ModeSymm=""
            endif
    
            print*, ""
            print'(A,I4,6X,A,F8.3,4X,A)', "Mode ", i, " Freq(cm^-1) = ", Freq(i), ModeSymm
            print*, "      S        Coef.     Contrib.(%)  ContribCorr(%)          Description"
            print*, " ======================================================================="
            Theta = 0.d0
            kk=0
            do j=1,Ns
                if (Theta > 0.9d0) exit
                jj = ipiv(j)
                Theta = Theta + Aux3(1,j)!**2
                print'(5x,i4,4x,g10.3,4x,2(f8.3,4x),a)', jj, L(jj,i), Aux3(1,j)*100, Aux(1,j)*100, trim(adjustl(ICDef(jj)))
                kk=kk+1
            enddo 
            print*, " ========================================================================"
            write(6,'(A,I3)') "     Total Number of internal to describe >90% of the mode: ", kk
        enddo

        ! Only give analysis by rows if verbosity level >= 1
        if (verbose<=2) return
    
        if (Nvib /= Ns) then
            call alert_msg("warning","Analysis by rows not implemented for redundant sets")
            return
        endif

        !The former trasformation related S in terms of Q's. To obtain the inverse relation:
        !Inverse of L (needs Asel for redundant sets (TBD))
        Aux(1:Ns,1:Nvib)=inverse_realgen(Nvib,L)

        print*, ""
        print*, "L^-1 MATRIX (LARGER ELEMENTS) -by rows"
        do i=1,Nvib
            !Copy the row, normalize and reorder 
            Aux3(1,1:Nvib) = abs(Aux(i,1:Nvib))
            Theta = 0.d0
            do ii=1,Nvib
                Theta = Theta + Aux3(1,ii)!**2
            enddo
            Aux3(1,1:Nvib) = Aux3(1,1:Nvib)/Theta
    !         Aux3(1,1:Nvib) = Aux3(1,1:Nvib)/dsqrt(Theta)
            call sort_vec_max(Aux3(1,1:Nvib),ipiv(1:Nvib),Nvib)
    

            if (present(Ssym)) then
                !Determine symmetry
                do j=1,Nvib
                    jj=ipiv(j)
                    if (Ssym(jj) == jj) cycle
                    ii = Ssym(jj)
                    Theta = Aux(i,jj)/Aux(i,ii)
                    if (Theta > 0.d0) then 
                        ModeSymm="A' "
                    else
                        ModeSymm="A''"
                    endif
                    exit
                enddo
            else
                ModeSymm=""
            endif

            print*, ""
            print'(A,I4,6X,A,F8.3,4X,A)', "Mode ", i, " Freq(cm^-1) = ", Freq(i), ModeSymm
            print*, "      S        Coef.     Contrib.(%)          Description"
            print*, " ======================================================================="
            Theta = 0.d0
            kk=0
            do j=1,Nvib
                if (Theta > 0.9d0) exit
                jj = ipiv(j)
                Theta = Theta + Aux3(1,j)!**2
                print'(5x,i4,4x,f9.3,4x,f8.3,4x,a)', jj, Aux(i,jj), Aux3(1,j)*100, trim(adjustl(ICDef(jj)))
                kk=kk+1
            enddo 
            print*, " ========================================================================"
            write(6,'(A,I3)') "     Total Number of internal to describe >90% of the mode: ", kk
        enddo
    
        return

    end subroutine analyze_internal


    subroutine calc_BDer(molec,Ns,Bder,analytical)
    
        use structure_types
        use verbosity
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: delta = 1.889726133d-3 !for numerical ders, in bohr(=10^-3 \AA, as Num freq in G09)
    
        !====================== 
        !ARGUMENTS
        type(str_resmol),intent(inout)       :: molec
        integer,intent(in)                   :: Ns
        real(8),dimension(:,:,:),intent(out) :: Bder    ! (Ns x 3Nat x 3Nat) the last index is the second der
        logical,intent(in),optional          :: analytical
        !======================
        ! Parts of Bder
        real(8),dimension(6,6)   :: BderStre
        real(8),dimension(9,9)   :: BderBend
        real(8),dimension(12,12) :: BderDihe, BderImpr
        character(len=50)        :: title

        !======================  
        !LOCAL
        integer :: Nat
        real(8),dimension(NDIM) :: S
        real(8),dimension(NDIM,NDIM) :: Bplus, Bmin
        real(8) :: X1,Y1,Z1,X2,Y2,Z2
        !Counters
        integer :: i,j,k, ii,jj,kk, irow,icol, is, i1,i2,i3,i4
        !swith mode
        logical :: do_analytical
        !======================  

        if (present(analytical) .and. analytical) then
            do_analytical=.true.
        else
            do_analytical=.false.
        endif

#ifndef WITH_LIBBDERS
        if (do_analytical) then
            call alert_msg("warning","Analytical derivatives of Bder not available "//&
                           "when not compiled --with-libbders. Using numerical code")
            do_analytical=.false.
        endif
#endif    

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with Atomic Units 
        current_units=molec%units
        call set_geom_units(molec,"Bohr")
        !----------------------------------------------

        !shortcuts
        Nat = molec%natoms
    
        if (verbose>0) then
            write(6,'(/,2X,A)') "Computing derivatives of Wilson B matrix..."
            write(6,'(X,A)')    "---------------------------------------------------------"
            if (do_analytical) then
                write(6,'(3X,A)') "with analytical derivatives"
            else   
                write(6,'(3X,A)') "with numerical derivatives"
            endif
        endif
    
        ! Computing only non-zero elements
        Bder(1:Ns,1:3*Nat,1:3*Nat) = 0.d0

        i=0
        ! STRETCHING
        do is=1,molec%geom%nbonds
            i=i+1
            i1 = molec%geom%bond(is,1)
            i2 = molec%geom%bond(is,2)
            if (do_analytical) then
#ifdef WITH_LIBBDERS
                call derBstre(BderStre,                                         &
                             molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
                             molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z)
#endif
!                 call DERSTRE(BderStre,                                      &
!                              molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
!                              molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z)
            else
                call dernumBstre(BderStre,                                      &
                             molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
                             molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z)
            endif
            !Retrieve actual Bder elements
            do irow=1,6
            do icol=1,6
                ! Select block
                if (irow<=3) then
                    jj = 3*i1-3  
                else
                    jj = 3*i2-6 
                endif
                if (icol<=3) then
                    kk = 3*i1-3  
                else
                    kk = 3*i2-6 
                endif
                Bder(i,jj+irow,kk+icol) = BderStre(irow,icol)
            enddo
            enddo

            if (verbose>2) then
                write(title,'(A,I0)') "BderStre - ",i 
                call MAT0(6,BderStre,6,6,title)
            endif  
        enddo 

        ! BENDING
        do is=1,molec%geom%nangles
            i=i+1
            i1 = molec%geom%angle(is,1)
            i3 = molec%geom%angle(is,2)
            i2 = molec%geom%angle(is,3)
            if (do_analytical) then
#ifdef WITH_LIBBDERS
                call derBbend(BderBend,                                         &
                             molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
                             molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z,&
                             molec%atom(i3)%x,molec%atom(i3)%y,molec%atom(i3)%z)
#endif
!                 call DERBEND(BderBend,                                          &
!                              molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
!                              molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z,&
!                              molec%atom(i3)%x,molec%atom(i3)%y,molec%atom(i3)%z)
            else
                call dernumBbend(BderBend,                                      &
                             molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
                             molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z,&
                             molec%atom(i3)%x,molec%atom(i3)%y,molec%atom(i3)%z)
            endif
            !Retrieve actual Bder elements
            do irow=1,9
            do icol=1,9
                ! Select block
                if (irow<=3) then
                    jj = 3*i1-3  
                elseif (irow<=6) then
                    jj = 3*i2-6 
                else
                    jj = 3*i3-9
                endif
                if (icol<=3) then
                    kk = 3*i1-3  
                elseif (icol<=6) then
                    kk = 3*i2-6 
                else
                    kk = 3*i3-9
                endif
                Bder(i,jj+irow,kk+icol) = BderBend(irow,icol)
            enddo
            enddo

            if (verbose>2) then
                write(title,'(A,I0)') "BderBend - ",i 
                call MAT0(6,BderBend,9,9,title)
            endif  
        enddo 

        ! DIHEDRAL
        do is=1,molec%geom%ndihed
            i=i+1
            i1 = molec%geom%dihed(is,1)
            i2 = molec%geom%dihed(is,2)
            i3 = molec%geom%dihed(is,3)
            i4 = molec%geom%dihed(is,4)
            if (do_analytical) then
#ifdef WITH_LIBBDERS
                call derBdihe(BderDihe,                                         &
                             molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
                             molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z,&
                             molec%atom(i3)%x,molec%atom(i3)%y,molec%atom(i3)%z,&
                             molec%atom(i4)%x,molec%atom(i4)%y,molec%atom(i4)%z)
#endif
!                 call DERTORS(BderDihe,                                          &
!                              molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
!                              molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z,&
!                              molec%atom(i3)%x,molec%atom(i3)%y,molec%atom(i3)%z,&
!                              molec%atom(i4)%x,molec%atom(i4)%y,molec%atom(i4)%z)
            else
                call dernumBdihe(BderDihe,                                      &
                             molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
                             molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z,&
                             molec%atom(i3)%x,molec%atom(i3)%y,molec%atom(i3)%z,&
                             molec%atom(i4)%x,molec%atom(i4)%y,molec%atom(i4)%z)
            endif
            !Retrieve actual Bder elements
            do irow=1,12
            do icol=1,12
                ! Select block
                if (irow<=3) then
                    jj = 3*i1-3  
                elseif (irow<=6) then
                    jj = 3*i2-6 
                elseif (irow<=9) then
                    jj = 3*i3-9 
                else
                    jj = 3*i4-12
                endif
                if (icol<=3) then
                    kk = 3*i1-3  
                elseif (icol<=6) then
                    kk = 3*i2-6 
                elseif (icol<=9) then
                    kk = 3*i3-9 
                else
                    kk = 3*i4-12
                endif
                Bder(i,jj+irow,kk+icol) = BderDihe(irow,icol)
            enddo
            enddo

            if (verbose>2) then
                write(title,'(A,I0)') "BderDihe - ",i 
                call MAT0(6,BderDihe,12,12,adjustl(title))
            endif  
        enddo  

        ! IMPROPERS
        do is=1,molec%geom%nimprop
            i=i+1
            i1 = molec%geom%improp(is,1)
            i2 = molec%geom%improp(is,2)
            i3 = molec%geom%improp(is,3)
            i4 = molec%geom%improp(is,4)
            ! No analytical derivatives (for the moment)
            if (do_analytical) & !then
                call alert_msg("note","No alasytical ders for impropers. Using numerical")
! #ifdef WITH_LIBBDERS
!                 call derBimpr(BderImpr,                                         &
!                              molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
!                              molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z,&
!                              molec%atom(i3)%x,molec%atom(i3)%y,molec%atom(i3)%z,&
!                              molec%atom(i4)%x,molec%atom(i4)%y,molec%atom(i4)%z)
! #endif
!             else
                call dernumBimpr(BderImpr,                                      &
                             molec%atom(i1)%x,molec%atom(i1)%y,molec%atom(i1)%z,&
                             molec%atom(i2)%x,molec%atom(i2)%y,molec%atom(i2)%z,&
                             molec%atom(i3)%x,molec%atom(i3)%y,molec%atom(i3)%z,&
                             molec%atom(i4)%x,molec%atom(i4)%y,molec%atom(i4)%z)
!             endif
            !Retrieve actual Bder elements
            do irow=1,12
            do icol=1,12
                ! Select block
                if (irow<=3) then
                    jj = 3*i1-3  
                elseif (irow<=6) then
                    jj = 3*i2-6 
                elseif (irow<=9) then
                    jj = 3*i3-9 
                else
                    jj = 3*i4-12
                endif
                if (icol<=3) then
                    kk = 3*i1-3  
                elseif (icol<=6) then
                    kk = 3*i2-6 
                elseif (icol<=9) then
                    kk = 3*i3-9 
                else
                    kk = 3*i4-12
                endif
                Bder(i,jj+irow,kk+icol) = BderImpr(irow,icol)
            enddo
            enddo

            if (verbose>2) then
                write(title,'(A,I0)') "BderImpr - ",i 
                call MAT0(6,BderImpr,12,12,adjustl(title))
            endif  
        enddo  

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! Revert original units
        call set_geom_units(molec,adjustl(current_units))
        !----------------------------------------------

        return
    
    end subroutine calc_BDer

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

    subroutine intshif2cart(molecule,DeltaS,&
                                            thr_set,maxiter_set,converged)

        !===========================================================
        ! Description
        ! -----------
        ! Subroutine to transform a shift in internal coordinates
        ! to Cartesian coordinates. The shift in applied onto the
        ! input molecule, so changing its coordianates
        !
        ! Notes
        ! -----
        ! * DeltaS need to be in AU/RAD
        ! * The basic input can replace zmat2cart. Additional options can be
        !   used to tune the behaviour
        ! * IMPORTANT: the molecule%geom must be consistent with DeltaS
        ! * It uses subroutines in this module to build B and Ginv
        !================================================================

        use constants
        use matrix

        integer,parameter :: NDIM = 600

        type(str_resmol),intent(inout)    :: molecule
        real(8),dimension(:),intent(in)   :: DeltaS
        real(8),intent(in),optional       :: thr_set
        integer,intent(in),optional       :: maxiter_set
        logical,intent(out),optional      :: converged

        ! Local
        real(8),dimension(NDIM,NDIM) :: A, B, G, Ginv
        real(8),dimension(NDIM) :: S0, S1, DSx, DStarget, DX, Snull
        real(8)                 :: Theta, Theta2, Theta3, thr, rmsd,rmsd0, Sbond,Sangle,Sdihed
        integer                 :: Nat, Ns, Nvib
        integer                 :: i,j,k, ii,jj,kk, iter, maxiter

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! This subroutine works with Atomic Units 
        current_units=molecule%units
        call set_geom_units(molecule,"Bohr")
        !----------------------------------------------

        ! Set optional arguments
        thr = 1.d-10
        maxiter = 100
        if (present(thr_set))     thr     = thr_set
        if (present(maxiter_set)) maxiter = maxiter_set

        ! shortcuts
        Ns  = molecule%geom%nbonds +&
              molecule%geom%nangles+&
              molecule%geom%ndihed
        Nat = molecule%natoms
      
        if (verbose>0) then  
            print'(/,A)', "================================================="
            print'(A)', " Iterative process to get the DeltaX from DeltaS"
            print'(A,/)', "================================================="
        endif
        if (verbose>1) then
            print'(/,A)', " Iter     RMSD-str(Bohr)   RMSD-bond(Bohr)  RMSD-angl(rad)   RMSD-dihed(rad)"
        endif

        !------------------------------------------------------------
        ! ** Update transformation matrix ** 
        ! Build B and G (silently)
        call verbose_mute()
        call internal_Wilson(molecule,Ns,S1,B)
        call internal_Gmetric(Nat,Ns,molecule%atom(:)%mass,B,G)
        call generalized_inv(Ns,Nvib,G,Ginv)
        call verbose_continue()
        
        ! A = M^-1 B^t G^-
        do i=1,Ns
        do j=1,3*Nat
            jj=(j-1)/3+1
            A(i,j) = B(i,j) / molecule%atom(jj)%mass / AMUtoAU
        enddo
        enddo
        A(1:3*Nat,1:Ns) = matrix_product(3*Nat,Ns,Ns,A,Ginv,tA=.true.)
        !------------------------------------------------------------
        
        ! Initialization
        iter=0
        rmsd=1.d10
        DStarget(1:Ns) = DeltaS(1:Ns)

        do while (rmsd>thr .and. iter<maxiter)
            S0(1:Ns) = S1(1:Ns)
            iter = iter+1
            rmsd = 0.d0

            ! And get DeltaX as DX
            do i=1,3*Nat
                DX(i) = 0.d0
                do k=1,Ns
                    DX(i) = DX(i) + A(i,k) * DStarget(k) 
                enddo
                rmsd = rmsd + DX(i)**2
            enddo
            rmsd = dsqrt(rmsd/dfloat(Nat))
            if (iter==1) rmsd0=rmsd 
            if (rmsd>rmsd0) then
                exit
            endif
            do i=1,Nat
                ii = 3*(i-1)
                molecule%atom(i)%x = molecule%atom(i)%x + DX(ii+1) 
                molecule%atom(i)%y = molecule%atom(i)%y + DX(ii+2)
                molecule%atom(i)%z = molecule%atom(i)%z + DX(ii+3)
            enddo
        
            !------------------------------------------------------------
            ! ** Update transformation matrix ** 
            ! Build B and G (silently)
            call verbose_mute()
            call internal_Wilson(molecule,Ns,S1,B)
            call internal_Gmetric(Nat,Ns,molecule%atom(:)%mass,B,G)
            call generalized_inv(Ns,Nvib,G,Ginv)
            call verbose_continue()
            
            ! A = M^-1 B^t G^-
            do i=1,Ns
            do j=1,3*Nat
                jj=(j-1)/3+1
                A(i,j) = B(i,j) / molecule%atom(jj)%mass / AMUtoAU
            enddo
            enddo
            A(1:3*Nat,1:Ns) = matrix_product(3*Nat,Ns,Ns,A,Ginv,tA=.true.)
            !------------------------------------------------------------

            ! DDs = Ds(curv) - Ds(x)
            Theta = 0.d0
            Sbond = 0.d0
            k=0
            do i=1,molecule%geom%nbonds
                k=k+1
                DSx(k) = S1(k)-S0(k)
                DStarget(k) = DStarget(k) - DSx(k)
                Sbond = max(Sbond,dabs(DStarget(k)))
                Theta = Theta+DStarget(k)**2
            enddo
            Theta = dsqrt(Theta/Nat)
            Theta2 = 0.d0
            Sangle = 0.d0
            do i=1, molecule%geom%nangles
                k=k+1
                DSx(k) = S1(k)-S0(k)
                DStarget(k) = DStarget(k) - DSx(k)
                Sangle = max(Sangle,dabs(DStarget(k)))
                Theta2 = Theta2+DStarget(k)**2
            enddo
            Theta2 = dsqrt(Theta2/Nat)
            Theta3 = 0.d0
            Sdihed = 0.d0
            do i=1,molecule%geom%ndihed
                k=k+1
                DSx(k) = S1(k)-S0(k)
                if (abs(DSx(k)) > abs(DSx(k)-2*PI) ) then
                    DSx(k) = DSx(k)-2*PI
                else if (abs(DSx(k)) > abs(DSx(k)+2*PI) ) then
                    DSx(k) = DSx(k)+2*PI
                endif
                DStarget(k) = DStarget(k) - DSx(k)
                Sdihed = max(Sdihed,dabs(DStarget(k)))
                Theta3 = Theta3+DStarget(k)**2
            enddo
            Theta3 = dsqrt(Theta3/Nat)
        
            if (verbose>1) then
                print'(X,I4,2X,4(ES16.6,X))', iter, rmsd, Theta, Theta2, Theta3
            endif
        
        enddo

        if (iter>=maxiter .and. rmsd>=thr) then
            print*, "Not converged!"
            print*, "Iterations = ", iter
            print'(X,A,ES16.6)', "RMSD(Bohr) = ", rmsd 
            print'(X,A,F12.6)',  "MaxDev-bond(AA) = ", Sbond * BOHRtoANGS
            print'(X,A,F12.6)',  "MaxDev-angl(deg)= ", Sangle* 180.d0/PI
            print'(X,A,F12.6)',  "MaxDev-dihe(deg)= ", Sdihed* 180.d0/PI
            print*, ""
            call alert_msg("warning","DeltaS to DeltaX not converged")
            if (present(converged)) converged=.false.
        elseif (rmsd0<rmsd .and. iter==2) then
            print*, "Diverged!"
            print*, "Iterations = ", iter
            print'(X,A,ES16.6)', "RMSD(Bohr) = ", rmsd 
            print'(X,A,F12.6)',  "MaxDev-bond(AA) = ", Sbond * BOHRtoANGS
            print'(X,A,F12.6)',  "MaxDev-angl(deg)= ", Sangle* 180.d0/PI
            print'(X,A,F12.6)',  "MaxDev-dihe(deg)= ", Sdihed* 180.d0/PI
            print*, ""
            call alert_msg("warning","DeltaS to DeltaX diverged. Iter 1 taken")
        elseif (rmsd0<rmsd) then
            print*, "Diverged!"
            print*, "Iterations = ", iter
            print'(X,A,ES16.6)', "RMSD(Bohr) = ", rmsd 
            print'(X,A,F12.6)',  "MaxDev-bond(AA) = ", Sbond * BOHRtoANGS
            print'(X,A,F12.6)',  "MaxDev-angl(deg)= ", Sangle* 180.d0/PI
            print'(X,A,F12.6)',  "MaxDev-dihe(deg)= ", Sdihed* 180.d0/PI
            print*, ""
            call alert_msg("fatal","DeltaS to DeltaX diverged after second iter")
        else
            if (verbose>0) then
                print*, "Converged!"
                print*, "Iterations = ", iter
                print'(X,A,ES16.6)', "RMSD(Bohr) = ", rmsd
                print'(X,A,F12.6)',  "MaxDev-bond(AA) = ", Sbond * BOHRtoANGS
                print'(X,A,F12.6)',  "MaxDev-angl(deg)= ", Sangle* 180.d0/PI
                print'(X,A,F12.6)',  "MaxDev-dihe(deg)= ", Sdihed* 180.d0/PI
                print*, ""
            endif
            if (present(converged)) converged=.true.
        endif

        !----------------------------------------------
        ! UNITS MANAGEMENT
        ! Revert original units
        call set_geom_units(molecule,adjustl(current_units))
        !----------------------------------------------

        return

    end subroutine intshif2cart

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
!     
!     subroutine axis_swithching(molec,molec2,T0)
!     
!         ! Subroutine to account for the axis swithching effect. molec2 is rotated to
!         ! fulfil the Eckart frame simulataneously with molec
!     
!         use structure_types
!         use matrix
!     
!     
!         type(str_resmol),intent(in)    :: molec
!         type(str_resmol),intent(inout) :: molec2 
!         real(8),dimension(3,3)          :: T0
!     
!         !Local
!         type(str_resmol)   :: molec_aux
!         real(8),dimension(1:3,1:3)     :: C,Cinv
!         real(8) :: xaux,yaux,zaux, xrot,yrot,zrot, rsum, rcheck=99999.d0
!         integer :: Nat
!         !Other auxiliar
!         real(8),dimension(3,3) :: Aux,AuxT,T
!         real(8),dimension(3)   :: Vec,Vec2
!     
!         !Stuff to use lapack inversion (should be in matrix)
!         integer,parameter :: NDIM=6
!         real(8),dimension(NDIM,NDIM) :: work
!         integer,dimension(NDIM,NDIM) :: ipiv
!         integer :: info
!     
!         !Load shortcuts
!         Nat = molec%natoms
!     
!         !Following the expresion given by Borrelli 2006
!     
!     !=========================
!         C(1,1) = 0.d0
!         C(2,2) = 0.d0
!         C(3,3) = 0.d0
!         C(1,2) = 0.d0
!         C(2,1) = 0.d0
!         C(1,3) = 0.d0
!         C(3,1) = 0.d0
!         C(2,3) = 0.d0
!         C(3,2) = 0.d0
!         do i=1,Nat
!             C(1,1) = C(1,1) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%x
!             C(2,2) = C(2,2) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%y
!             C(3,3) = C(3,3) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%z
!     
!             C(1,2) = C(1,2) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%y
!             C(1,3) = C(1,3) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%z
!             C(2,3) = C(2,3) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%z
!     
!             C(2,1) = C(2,1) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%x
!             C(3,1) = C(3,1) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%x
!             C(3,2) = C(3,2) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%y
!         enddo
!     !     print*, "C_ini"
!     !     do i=1,3
!     !         print'(100(F8.3,2X))', C(i,1:3)
!     !     enddo
!         !Inverse of C
!     !     nn=3
!         Cinv=C
!         call dgetrf(3,3, Cinv, 3, ipiv, info)
!         call dgetri(3,   Cinv, 3, ipiv, work, NDIM, info)
!     
!     
!     !     Aux=0.d0
!     !     Aux(1:3,1:3) = matmul(C(1:3,1:3),Cinv(1:3,1:3))
!     !     print*, "CinvC", info
!     !     do i=1,3
!     !         print'(100(F8.3,2X))', Aux(i,1:3)
!     !     enddo
!     
!     
!         do i=1,3
!             do j=1,3
!                 Aux(i,j) = 0.d0
!                 do  k=1,3
!                     Aux(i,j) = Aux(i,j) + C(k,i)*C(k,j)
!                 enddo
!            enddo
!         enddo
!         call diagonalize_full(Aux(1:3,1:3),3,AuxT(1:3,1:3),Vec(1:3),"lapack")
!     
!     
!     !     do i=1,8
!     ! 
!     !     Vec2=1.d0
!     ! 
!     !         iT = 2*i/3
!     !         jT = 2*(i+1)/3
!     !         kT = 2*(i+2)/3
!     ! 
!     !         Vec2(1) = (-1.d0)**iT
!     !         Vec2(2) = (-1.d0)**jT
!     !         Vec2(3) = (-1.d0)**kT
!     ! 
!     !         print'(3F8.3,X,3I3)', Vec2(1:3), iT,jT,kT
!     ! 
!     !     enddo
!     ! stop
!     
!         !Check all 8 possible solutions for T0 (see Sando, 2001) 
!         do iT=1,8
!     
!             Vec2=1.d0
!             if (iT==2) then
!                 Vec2(1) = -1.d0
!             else if (iT==3) then
!                 Vec2(2) = -1.d0
!             else if (iT==4) then
!                 Vec2(3) = -1.d0
!             else if (iT==5) then
!                 Vec2(1) = -1.d0
!                 Vec2(2) = -1.d0
!             else if (iT==6) then
!                 Vec2(1) = -1.d0
!                 Vec2(3) = -1.d0
!             else if (iT==7) then
!                 Vec2(2) = -1.d0
!                 Vec2(3) = -1.d0
!             else if (iT==8) then
!                 Vec2(1) = -1.d0
!                 Vec2(2) = -1.d0
!                 Vec2(3) = -1.d0
!             endif
!         
!     
!         !(C^t C)^1/2
!         do i=1,3
!             do j=1,3
!                 C(i,j) = 0.d0
!                 do k=1,3
!                     if (dabs(Vec(k)) < 1d-10) Vec(k) = 0.d0
!                     C(i,j) = C(i,j) + AuxT(i,k)*dsqrt(Vec(k))*AuxT(j,k)*Vec2(k)
!                 enddo
!             enddo
!         enddo
!     
!         !T = C^t C)^1/2 C^-1
!         do i=1,3
!             do j=1,3
!                 T(i,j) = 0.d0
!                 do k=1,3
!                     T(i,j) = T(i,j) + C(i,k)*Cinv(k,j)
!                 enddo
!             enddo
!         enddo
!         !Manage 0-eigenvalues in diagonal matrix
!         Vec2(1:3)=0.d0
!         do i=1,3
!             do j=1,3
!                Vec2(i) = Vec2(i) + T(i,j)**2
!             enddo
!         enddo
!         !To do: add support for one-off diagonal eigenvalue in C
!         ! If one eigenvalue is 0
!         if (abs(Vec2(1)) < 1.d-10) then
!             print*, "Warning: eingen 1 value is zero"
!             T(1,1) = T(2,2)*T(3,3)-T(2,3)*T(3,2)       
!             T(1,2) = T(2,3)*T(3,1)-T(2,1)*T(3,3)       
!             T(1,3) = T(2,1)*T(3,2)-T(2,2)*T(3,1)       
!         else if (abs(Vec2(2)) < 1.d-10) then
!             print*, "Warning: eingen 2 value is zero"
!             T(2,1) = T(1,2)*T(3,3)-T(1,3)*T(3,2)        
!             T(2,2) = T(1,3)*T(3,1)-T(1,1)*T(3,3)        
!             T(2,3) = T(1,1)*T(3,2)-T(1,2)*T(3,1)       
!         else if (abs(Vec2(3)) < 1.d-10) then
!             print*, "Warning: eingen 3 value is zero"
!             T(3,1) = T(1,2)*T(2,3)-T(1,3)*T(2,2)        
!             T(3,2) = T(1,3)*T(2,1)-T(1,1)*T(2,3)        
!             T(3,3) = T(1,1)*T(2,2)-T(1,2)*T(2,1)   
!         endif 
!             
!         print*, "T",iT
!         do i=1,3
!             print'(100(F8.3,2X))', T(i,1:3)
!         enddo
!     
!         det =       T(1,1)*T(2,2)*T(3,3)
!         det = det + T(2,1)*T(3,2)*T(1,3)
!         det = det + T(1,2)*T(2,3)*T(3,1)
!         det = det - T(3,1)*T(2,2)*T(1,3)
!         det = det - T(2,1)*T(1,2)*T(3,3)
!         det = det - T(3,2)*T(2,3)*T(1,1)
!     
!         print*, det
!     
!         if (det<0.d0) cycle
!     
!         !Rotate and check sum of distances
!         molec_aux=molec2
!         rsum=0.d0
!         do i=1,Nat
!             xaux = molec2%atom(i)%x
!             yaux = molec2%atom(i)%y
!             zaux = molec2%atom(i)%z
!             xrot = T(1,1)*xaux + T(2,1)*yaux + T(3,1)*zaux - molec%atom(i)%x
!             yrot = T(1,2)*xaux + T(2,2)*yaux + T(3,2)*zaux - molec%atom(i)%y
!             zrot = T(1,3)*xaux + T(2,3)*yaux + T(3,3)*zaux - molec%atom(i)%z
!             rsum=rsum + dsqrt(xrot**2+yrot**2+zrot**2)
!         enddo
!         print*, "Sum of distance ", rsum
!     
!         if (rsum < rcheck) then
!             T0=T
!             rcheck=rsum
!         endif
!     
!         enddo ! Possible T values
!     
!         !Rotate State 2 
!         !---------------
!         print*, "T0"
!         do i=1,3
!             print'(100(F8.3,2X))', T0(i,1:3)
!         enddo
!         det =       T0(1,1)*T0(2,2)*T0(3,3)
!         det = det + T0(2,1)*T0(3,2)*T0(1,3)
!         det = det + T0(1,2)*T0(2,3)*T0(3,1)
!         det = det - T0(3,1)*T0(2,2)*T0(1,3)
!         det = det - T0(2,1)*T0(1,2)*T0(3,3)
!         det = det - T0(3,2)*T0(2,3)*T0(1,1)
!         print*, "Det", det
!         rsum=0.d0
!         do i=1,Nat
!             xaux = molec2%atom(i)%x
!             yaux = molec2%atom(i)%y
!             zaux = molec2%atom(i)%z
!             molec2%atom(i)%x = T0(1,1)*xaux + T0(2,1)*yaux + T0(3,1)*zaux
!             molec2%atom(i)%y = T0(1,2)*xaux + T0(2,2)*yaux + T0(3,2)*zaux
!             molec2%atom(i)%z = T0(1,3)*xaux + T0(2,3)*yaux + T0(3,3)*zaux
!             rsum=rsum + dsqrt((molec2%atom(i)%x-molec%atom(i)%x)**2+&
!                               (molec2%atom(i)%y-molec%atom(i)%y)**2+&
!                               (molec2%atom(i)%z-molec%atom(i)%z)**2)
!         enddo
!         print*, "Sum of distance ", rsum
!         !Check new C matrix
!         !C(1,1)
!         C(1,1) = 0.d0
!         C(2,2) = 0.d0
!         C(3,3) = 0.d0
!         C(1,2) = 0.d0
!         C(2,1) = 0.d0
!         C(1,3) = 0.d0
!         C(3,1) = 0.d0
!         C(2,3) = 0.d0
!         C(3,2) = 0.d0
!         do i=1,Nat
!             C(1,1) = C(1,1) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%x
!             C(2,2) = C(2,2) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%y
!             C(3,3) = C(3,3) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%z
!     
!             C(1,2) = C(1,2) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%y
!             C(1,3) = C(1,3) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%z
!             C(2,3) = C(2,3) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%z
!     
!             C(2,1) = C(2,1) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%x
!             C(3,1) = C(3,1) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%x
!             C(3,2) = C(3,2) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%y
!         enddo
!         print*, "C_new"
!         do i=1,3
!             print'(100(F8.3,2X))', C(i,1:3)
!         enddo
!     !=========================
!     
!         return
!     end subroutine axis_swithching
!     


end module internal_module
