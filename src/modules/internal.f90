module internal_module

    use matrix
    use matrix_print
    use verbosity

    implicit none

    contains


    subroutine define_internal_set(molec,def_internal,intfile,rmzfile,use_symmetry,isym,S_sym,Ns)

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
        character(len=*),intent(inout)     :: def_internal      ! switch with the way internal are selected
        character(len=*),intent(in)        :: intfile           ! file with the def of internal coords
        character(len=*),intent(in)        :: rmzfile           ! additional file with coordinates to remove from Zmat
        logical,intent(in)                 :: use_symmetry      ! logical to use or not symmetry
        integer,dimension(:),intent(in)    :: isym              ! array with the simmetric atom (in)
        integer,dimension(:),intent(out)   :: S_sym             ! array with symmetric internals (out)
        integer,intent(out)                :: Ns                ! Total number of internal coordiantes to use
        !====================== 
    
        !======================
        !LOCAL 
        !System info
        integer :: Nat
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
        call set_word_upper_case(def_internal)
        Nat = molec%natoms
        PG = molec%PG

        !GEN BONDED SET FOR INTERNAL COORD
        if (adjustl(def_internal) == "SEL") then
            print*, "Reading internal coordianates from: "//trim(adjustl(intfile))
            if (verbose>0) &
             open(I_FILE,file=intfile,iostat=IOstatus) 
            if (IOstatus /= 0) call alert_msg("fatal","Cannot open file: "//trim(adjustl(intfile)))
            ! Get internal coords (using modredundant sr)
            call modredundant(I_FILE,molec)
            close(I_FILE)

        elseif (adjustl(def_internal) == "ZMAT") then
            if (adjustl(intfile) == "none") then
                call build_Z(molec,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
            else
                open(I_FILE,file=intfile,status="old")
                print*, "Z-matrix read from "//trim(adjustl(intfile))
                call read_Z(I_FILE,molec,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym)
                close(I_FILE)
            endif

        else if (adjustl(def_internal) == "ALL") then!otherwise all parameters are used
            if (verbose>0) &
             print*, "Using all internal coordinates", molec%geom%nangles
        else
            call alert_msg("fatal","Unkownn option for internal definition. Valid options are 'zmat', 'sel', 'all'")
        endif 

        ! Compute the number of internal coords
        Ns = molec%geom%nbonds  + &
             molec%geom%nangles + &
             molec%geom%ndihed

        ! Remove some Zmat elements if required
        if (def_internal=="ZMAT" .and. rmzfile /= "none") then
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
                Ns = NS - 1
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


        return

    end subroutine define_internal_set



    subroutine internal_Wilson(molec,Ns,S,B, &
!                                           Optional:
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
        type(str_resmol),intent(in)        :: molec    ! Input molecule (but only use geom...) - 
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
        integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
        integer :: nbonds, ndihed, nangles
        integer :: Nat
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
        integer :: i_1, i_2, i_3, i_4
        !=============
    
    
        !Set bonded
        nbonds  = molec%geom%nbonds
        nangles = molec%geom%nangles
        ndihed  = molec%geom%ndihed
        bond_s(1:nbonds,1:2)  =  molec%geom%bond(1:nbonds,1:2)
        angle_s(1:nangles,1:3) =  molec%geom%angle(1:nangles,1:3)
        dihed_s(1:ndihed,1:4)  =  molec%geom%dihed(1:ndihed,1:4)
    
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
            r21 = calc_dist(molec%atom(i_1),molec%atom(i_2))
            S(k) = r21
            if (verbose>0) &
            write(6,'(I5,X,A,2(I3,A),2F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                                  trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                                               r21*BOHRtoAMS, r21
            if (present(ICDef)) &
            write(ICDef(k),'(A,2(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                                          trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")"
    
            !Two cart displacements different from zero.
             e21x = (molec%atom(i_1)%x-molec%atom(i_2)%x)/r21
             e21y = (molec%atom(i_1)%y-molec%atom(i_2)%y)/r21
             e21z = (molec%atom(i_1)%z-molec%atom(i_2)%z)/r21
            ! s1 = e21 = -e12
             !B, index: internal,cartesian -> xn = 3n-2; yn=3n-1; zn=3n
             B(k,3*i_1-2) = e21x
             B(k,3*i_1-1) = e21y
             B(k,3*i_1  ) = e21z
            ! s2 = e12 = -e21
             B(k,3*i_2-2) = -e21x
             B(k,3*i_2-1) = -e21y
             B(k,3*i_2  ) = -e21z
        enddo
    
        if (verbose>0) &
        write(6,'(/,A,I3)') "ANGLES", nangles
        do i=1,nangles
            k=k+1
    
            i_1 = angle_s(i,1)
            i_3 = angle_s(i,2)
            i_2 = angle_s(i,3)
            ang1 = calc_angle(molec%atom(i_1),molec%atom(i_3),molec%atom(i_2))
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
    
            !Three cart displacements different from zero.
             r31=calc_dist(molec%atom(i_1),molec%atom(i_3))
             e31x = (molec%atom(i_1)%x-molec%atom(i_3)%x)/r31
             e31y = (molec%atom(i_1)%y-molec%atom(i_3)%y)/r31
             e31z = (molec%atom(i_1)%z-molec%atom(i_3)%z)/r31
             r32=calc_dist(molec%atom(i_2),molec%atom(i_3))
             e32x = (molec%atom(i_2)%x-molec%atom(i_3)%x)/r32
             e32y = (molec%atom(i_2)%y-molec%atom(i_3)%y)/r32
             e32z = (molec%atom(i_2)%z-molec%atom(i_3)%z)/r32
             ! s1 = [ cos(ang1)*e31 - e32 ] / [ r31 sin(ang1)]
             B(k,3*i_1-2) = (dcos(ang1)*e31x-e32x)/(r31*dsin(ang1))
             B(k,3*i_1-1) = (dcos(ang1)*e31y-e32y)/(r31*dsin(ang1))
             B(k,3*i_1  ) = (dcos(ang1)*e31z-e32z)/(r31*dsin(ang1))
             ! s2 = [ cos(ang1)*e32 - e31 ] / [ r32 sin(ang1)]
             B(k,3*i_2-2) = (dcos(ang1)*e32x-e31x)/(r32*dsin(ang1))
             B(k,3*i_2-1) = (dcos(ang1)*e32y-e31y)/(r32*dsin(ang1))
             B(k,3*i_2  ) = (dcos(ang1)*e32z-e31z)/(r32*dsin(ang1))
             ! s3 = [(r31-r32 cos(ang1))e31 + (r32-r31 cos(ang1))e32 / [ rr3132 sin(ang1)]
             B(k,3*i_3-2) = ( (r31-r32*dcos(ang1))*e31x + (r32-r31*dcos(ang1))*e32x )&
                          / ( r31*r32*dsin(ang1))
             B(k,3*i_3-1) = ( (r31-r32*dcos(ang1))*e31y + (r32-r31*dcos(ang1))*e32y )&
                          / ( r31*r32*dsin(ang1))
             B(k,3*i_3  ) = ( (r31-r32*dcos(ang1))*e31z + (r32-r31*dcos(ang1))*e32z )&
                          / ( r31*r32*dsin(ang1))
        enddo
    
        if (verbose>0) &
        write(6,'(/,A,I3)') "DIHEDRALS", ndihed
        do i=1,ndihed
            k=k+1
    
            i_1 = dihed_s(i,1)
            i_2 = dihed_s(i,2)
            i_3 = dihed_s(i,3)
            i_4 = dihed_s(i,4)
            ang1 = calc_dihed_new(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
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
    
            !Four cart displacements different from zero (some index intercheged with Decius..)
            ! this was buggy. Changed on 10/06/2014
            ang2 = calc_angle(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3))
            ang3 = calc_angle(molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
            r21  = calc_dist(molec%atom(i_1),molec%atom(i_2))
            e21x = (molec%atom(i_1)%x-molec%atom(i_2)%x)/r21
            e21y = (molec%atom(i_1)%y-molec%atom(i_2)%y)/r21
            e21z = (molec%atom(i_1)%z-molec%atom(i_2)%z)/r21
            r32  = calc_dist(molec%atom(i_2),molec%atom(i_3))
            e32x = (molec%atom(i_2)%x-molec%atom(i_3)%x)/r32
            e32y = (molec%atom(i_2)%y-molec%atom(i_3)%y)/r32
            e32z = (molec%atom(i_2)%z-molec%atom(i_3)%z)/r32
            r43  = calc_dist(molec%atom(i_3),molec%atom(i_4))
            e43x = (molec%atom(i_3)%x-molec%atom(i_4)%x)/r43
            e43y = (molec%atom(i_3)%y-molec%atom(i_4)%y)/r43
            e43z = (molec%atom(i_3)%z-molec%atom(i_4)%z)/r43
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
            B(k,3*i_1-2) =  -e21Pe32x/(r21*dsin(ang2)**2)
            B(k,3*i_1-1) =  -e21Pe32y/(r21*dsin(ang2)**2)
            B(k,3*i_1  ) =  -e21Pe32z/(r21*dsin(ang2)**2)
            !s2 
            B(k,3*i_2-2) = ((r32-r21*dcos(ang2))*e21Pe32x/(r32*r21*dsin(ang2)**2) &
                         +  dcos(ang3)*e43Pe32x/(r32*dsin(ang3)**2))
            B(k,3*i_2-1) = ((r32-r21*dcos(ang2))*e21Pe32y/(r32*r21*dsin(ang2)**2) &
                         +  dcos(ang3)*e43Pe32y/(r32*dsin(ang3)**2))
            B(k,3*i_2  ) = ((r32-r21*dcos(ang2))*e21Pe32z/(r32*r21*dsin(ang2)**2) &
                         +  dcos(ang3)*e43Pe32z/(r32*dsin(ang3)**2))
            !s3
            B(k,3*i_3-2) = ((r32-r43*dcos(ang3))*e43Pe32x/(r32*r43*dsin(ang3)**2) &
                         +  dcos(ang2)*e21Pe32x/(r32*dsin(ang2)**2))
            B(k,3*i_3-1) = ((r32-r43*dcos(ang3))*e43Pe32y/(r32*r43*dsin(ang3)**2) &
                         +  dcos(ang2)*e21Pe32y/(r32*dsin(ang2)**2))
            B(k,3*i_3  ) = ((r32-r43*dcos(ang3))*e43Pe32z/(r32*r43*dsin(ang3)**2) &
                         +  dcos(ang2)*e21Pe32z/(r32*dsin(ang2)**2))
            !s4
            B(k,3*i_4-2) =  -e43Pe32x/(r43*dsin(ang3)**2)
            B(k,3*i_4-1) =  -e43Pe32y/(r43*dsin(ang3)**2)
            B(k,3*i_4  ) =  -e43Pe32z/(r43*dsin(ang3)**2)
    
        enddo
    
!         write(6,'(/,A,I3)') "IMPROPERS", nimprop
!         do i=1,nimprop
!             k=k+1
!     
!             i_1 = dihed_s(i,1)
!             i_2 = dihed_s(i,2)
!             i_3 = dihed_s(i,3)
!             i_4 = dihed_s(i,4)
!             ang1 = calc_dihed_new(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
!             S(k) = ang1
!             if (verbose>0) &
!             write(6,'(I5,X,A,4(I3,A),F15.8)') k, trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
!                                                  trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
!                                                  trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
!                                                  trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")", &
!                                               ang1*360.d0/2.d0/pi
!             if (present(ICDef)) &
!             write(ICDef(k),'(A,4(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
!                                           trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
!                                           trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
!                                           trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")"
!     
!             !Four cart displacements different from zero (some index intercheged with Decius..)
!             ang2 = calc_angle(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3))
!             ang3 = calc_angle(molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
!             r21  = calc_dist(molec%atom(i_1),molec%atom(i_2))
!             e21x = (molec%atom(i_1)%x-molec%atom(i_2)%x)/r21
!             e21y = (molec%atom(i_1)%y-molec%atom(i_2)%y)/r21
!             e21z = (molec%atom(i_1)%z-molec%atom(i_2)%z)/r21
!             r32  = calc_dist(molec%atom(i_2),molec%atom(i_3))
!             e32x = (molec%atom(i_2)%x-molec%atom(i_3)%x)/r32
!             e32y = (molec%atom(i_2)%y-molec%atom(i_3)%y)/r32
!             e32z = (molec%atom(i_2)%z-molec%atom(i_3)%z)/r32
!             r43  = calc_dist(molec%atom(i_3),molec%atom(i_4))
!             e43x = (molec%atom(i_3)%x-molec%atom(i_4)%x)/r43
!             e43y = (molec%atom(i_3)%y-molec%atom(i_4)%y)/r43
!             e43z = (molec%atom(i_3)%z-molec%atom(i_4)%z)/r43
!             e21Pe32x=e21y*e32z-e21z*e32y
!             e21Pe32y=e21z*e32x-e21x*e32z
!             e21Pe32z=e21x*e32y-e21y*e32x
!             e32Pe21Pe32x=e32y*e21Pe32z-e32z*e21Pe32y
!             e32Pe21Pe32y=e32z*e21Pe32x-e32x*e21Pe32z
!             e32Pe21Pe32z=e32x*e21Pe32y-e32y*e21Pe32x
!             e43Pe32x=e43y*e32z-e43z*e32y
!             e43Pe32y=e43z*e32x-e43x*e32z
!             e43Pe32z=e43x*e32y-e43y*e32x
!             e32Pe43Pe32x=e32y*e43Pe32z-e32z*e43Pe32y
!             e32Pe43Pe32y=e32z*e43Pe32x-e32x*e43Pe32z
!             e32Pe43Pe32z=e32x*e43Pe32y-e32y*e43Pe32x
!     
!             !s1...
!             
!         enddo
    
        if (verbose>1) &
            call MAT0(6,B,Ns,3*Nat,"B MATRIX")

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
            call MAT0(6,G,Ns,Ns,"G MATRIX")

        return

    end subroutine internal_Gmetric


    subroutine HessianCart2int(Nat,Ns,Hess,Mass,B,G, &
                                                 !Optional:
                                                 Grad,Bder)

        !==============================================================
        ! This code is part of MOLECULAR_TOOLS 
        !==============================================================
        ! Description
        !  HESSIAN IN INTERNAL COORDINATES (JCC, 17, 49-56, by Frisch et al)
        !    Hint = G^- Bu(Hx+B'^Tg_q)u^TB^T G^-
        !   g_q is the gradient, so g_q=0 in a minimum
        !   G^- is the generalized inverse (for redundant internal) or simply the
        !   inverse for nonredundant
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
        real(8),dimension(1:NDIM,1:NDIM,1:NDIM),intent(in),optional :: Bder ! Bmatrxi derivatives
        real(8),dimension(1:NDIM,1:NDIM),intent(inout) :: Hess   !Hessian: cart(in)-intern(out)
        real(8),dimension(1:NDIM),intent(inout),optional :: Grad !Gradient 
        !====================== 
    
        !====================== 
        !LOCAL
        integer                             :: Nvib
        !Internal analysis 
        real(8),dimension(1:NDIM,1:NDIM)    :: Ginv
        !Auxiliar arrays
        real(8),dimension(1:NDIM,1:NDIM)    :: AuxT,Aux
        real(8),dimension(NDIM)             :: Vec
        !Counters
        integer :: i,j,k, ii
        !====================== 

    
        ! If Ns > Nvib, then we need to extract the non-zero eigen values 
        ! from G in order to compute the generalized G inverse
        ! (TBD) 
        ! 
        Nvib = Ns

        !Inverse of G
        Ginv(1:Ns,1:Ns) = inverse_realsym(Ns,G)
    
        !Compute G^-1Bu  (where u is the inverse mass matrix)
        Aux(1:Ns,1:3*Nat) = matrix_product(Ns,3*Nat,Ns,Ginv,B)
        do i=1,3*Nat
            ii = (i-1)/3+1
            Aux(1:Ns,i) = Aux(1:Ns,i)/Mass(ii)/UMAtoAU
        enddo
    
        if (present(Grad)) then
            if (.not.present(Bder)) call alert_msg("fatal","API Error: Bder needed with Grad present")
            if (verbose>0) &
            print*, "Using Gradient..."
            ! Get the gradient in internal coords first: gq = G^-1Bu(gx)
            do i=1,Nvib
                Vec(i) = 0.d0
                do j=1,3*Nat
                    Vec(i) = Vec(i) + Aux(i,j) * Grad(j)
                enddo
            enddo
            ! Update the gradient on output
            Grad(1:3*Nat) = 0.d0
            Grad(1:Nvib) = Vec(1:Nvib)
            ! .. and multiply: Bder(i,j,K)^t * gq(K)
            do i=1,3*Nat
            do j=1,3*Nat
                AuxT(i,j) = 0.d0
                do k=1,Nvib
                    AuxT(i,j) = AuxT(i,j) + Bder(k,i,j)*Vec(k)
                enddo
                ! Apply correction to the Hessian term
                Hess(i,j) = Hess(i,j) - AuxT(i,j)
            enddo
            enddo
        endif ! gradient correction
    
        ! Hint = Aux (Hcart-Bder*gq) Aux^T
        Hess(1:3*Ns,1:3*Nat) = matrix_product(Ns,3*Nat,3*Nat,Aux,Hess)
        Hess(1:3*Ns,1:3*Ns)  = matrix_product(Ns,3*Nat,3*Nat,Hess,Aux,tB=.true.)
    
        if (verbose>1) &
            call MAT0(6,Hess,Ns,Ns,"F MATRIX")

        return

    end subroutine HessianCart2int

    
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

    
        !Now rotate F, F'=X^TFX. Store the rotated matrix in Aux3 (temporary array)
        Aux3(1:Nvib,1:Nvib) = matrix_basisrot(Nvib,Nvib,X,Hess,counter=.true.)
    
        !We can now diagonalize F'
        call diagonalize_full(Aux3(1:Nvib,1:Nvib),Nvib,L(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
    
        !Check FC
        if (verbose>0) &
            call print_vector(6,Freq,Nvib,"FORCE CONSTANTS (A.U.)")

        !Check freqcuencies
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


    subroutine analyze_internal(Nvib,L,Freq,ICDef,Ssym)

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
        integer,intent(in)                            :: Nvib
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
            Aux3(1,1:Nvib) = abs(L(1:Nvib,i))
            Theta  = 0.d0
            Theta2 = 0.d0
            do ii=1,Nvib
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
                !Displacement-weighted calc
                Theta2 = Theta2 + Aux(1,ii)
            enddo
            Aux3(1,1:Nvib) = Aux3(1,1:Nvib)/Theta
            Aux(1,1:Nvib)  = Aux(1,1:Nvib)/Theta2
    !         Aux3(1,1:Nvib) = Aux3(1,1:Nvib)/dsqrt(Theta)
            call sort_vec_max(Aux3(1,1:Nvib),ipiv(1:Nvib),Nvib)
            call sort_vec_max(Aux(1,1:Nvib),ipiv2(1:Nvib),Nvib)
    
            if (present(Ssym)) then
                !Determine symmetry
                do j=1,Nvib
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
            print'(A,I4,6X,A,F8.3,4X,A,A,F8.3)', "Mode ", i, " Freq(cm^-1) = ", Freq(i), "Symm ", ModeSymm
            print*, "      S        Coef.     Contrib.(%)  ContribCorr(%)          Description"
            print*, " ======================================================================="
            Theta = 0.d0
            kk=0
            do j=1,Nvib
                if (Theta > 0.9d0) exit
                jj = ipiv(j)
                Theta = Theta + Aux3(1,j)!**2
                print'(5x,i4,4x,g10.3,4x,2(f8.3,4x),a)', jj, L(jj,i), Aux3(1,j)*100, Aux(1,j)*100, trim(adjustl(ICDef(jj)))
                kk=kk+1
            enddo 
            print*, " ========================================================================"
            write(6,'(A,I3)') "     Total Number of internal to describe >90% of the mode: ", kk
        enddo
    
    
        !The former trasformation related S in terms of Q's. To obtain the inverse relation:
        !Inverse of L 
        Aux(1:Nvib,1:Nvib)=inverse_realgen(Nvib,L)
        ! Check inversion
        Aux3(1:Nvib,1:Nvib) = matrix_product(Nvib,Nvib,Nvib,Aux,L)
        if (verbose>2) & 
            call MAT0(6,Aux3,Nvib,Nvib,"L L^-1")

    
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
            print'(A,I4,6X,A,F8.3,4X,A,A,F8.3)', "Mode ", i, " Freq(cm^-1) = ", Freq(i), "Symm ", ModeSymm
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
    


    subroutine NumBDer(molec,Ns,Bder)
    
        use structure_types
        use verbosity
    
        integer,parameter :: NDIM = 600
        real(8),parameter :: delta = 1.889726133d-3 !for numerical ders, in bohr(=10^-3 \AA, as Num freq in G09)
    
        !====================== 
        !ARGUMENTS
        type(str_resmol),intent(in)          :: molec
        integer,intent(in)                   :: Ns
        real(8),dimension(:,:,:),intent(out) :: Bder    ! (Ns x 3Nat x 3Nat) the last index is the second der
        !======================

        !======================  
        !LOCAL
        integer :: Nat
        real(8),dimension(NDIM) :: S
        type(str_resmol) :: molecB
        real(8),dimension(NDIM,NDIM) :: Bplus, Bmin
        integer :: verbose_current
        !Counters
        integer :: i,j,k, ii
        !======================  

        !Get current verbose level
        verbose_current = verbose
        !And set to quiet
        verbose = 0
    
        !shortcuts
        Nat = molec%natoms
    
        if (verbose_current>0) then
            print*, ""
            print*, "COMPUTING NUMERICAL DERIVATIVES FOR B..."
            print*, ""
        endif
    
        do i=1,Nat
    
            !Displace X
            molecB = molec
            ii = 3*i-2
            molecB%atom(i)%x = molec%atom(i)%x + delta
            !Call B matrix at this geometry   
            call internal_Wilson(molecB,Ns,S,Bplus)
            molecB%atom(i)%x = molec%atom(i)%x - delta
            !Call B matrix at this geometry   
            call internal_Wilson(molecB,Ns,S,Bmin)
    
            do j=1,Ns
            do k=1,3*Nat
               Bder(j,k,ii) = ( Bplus(j,k) - Bmin(j,k) ) / (2.d0*delta)
               if (verbose_current>1) &
                 print*, j,k,ii, Bder(j,k,ii), Bplus(j,k), Bmin(j,k)
            enddo
            enddo
    
            !Displace Y
            molecB = molec
            ii = 3*i-1
            molecB%atom(i)%y = molec%atom(i)%y + delta
            !Call B matrix at this geometry   
            call internal_Wilson(molecB,Ns,S,Bplus)
            molecB%atom(i)%y = molec%atom(i)%y - delta
            !Call B matrix at this geometry   
            call internal_Wilson(molecB,Ns,S,Bmin)
    
            do j=1,Ns
            do k=1,3*Nat
               Bder(j,k,ii) = ( Bplus(j,k) - Bmin(j,k) ) / (2.d0*delta)
               if (verbose_current>1) &
                 print*, j,k,ii, Bder(j,k,ii), Bplus(j,k), Bmin(j,k)
            enddo
            enddo
    
            !Displace Z
            molecB = molec
            ii = 3*i
            molecB%atom(i)%z = molec%atom(i)%z + delta
            !Call B matrix at this geometry   
            call internal_Wilson(molecB,Ns,S,Bplus)
            molecB%atom(i)%z = molec%atom(i)%z - delta
            !Call B matrix at this geometry   
            call internal_Wilson(molecB,Ns,S,Bmin)
    
            do j=1,Ns
            do k=1,3*Nat
               Bder(j,k,ii) = ( Bplus(j,k) - Bmin(j,k) ) / (2.d0*delta)
               if (verbose_current>1) &
                 print*, j,k,ii, Bder(j,k,ii), Bplus(j,k), Bmin(j,k)
            enddo
            enddo
    
        enddo
    
        !Restore verbose level
        verbose =  verbose_current

        return
    
    end subroutine NumBDer


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
