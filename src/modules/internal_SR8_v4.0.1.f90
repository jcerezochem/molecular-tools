module internal_module

 !The hessian is passed to gf_method instead of the input label
 ! Version 3: rearragement to tidy up things and solve compatibility issues
 ! Version 4: L(non orth) is passed back in gf_method instead of L'(orth)
 ! Version 5: str_residue replaced by str_resmol. V4 with general clean up
 ! Version 6: Added determination of the symmetric internal coordinates
 ! Version 7: Added %contrib to Normal mode analysis in terms of internal (gf_method)
 !            NDIM increased to 600
 ! Version 8: Added proper redundant treatment (testing)
 !            Include readZ SR

 !Contains
! subroutine build_Z(molec,bond_s,angle_s,dihed_s,PG,isym)
! 
! subroutine read_Z(molec,bond_s,angle_s,dihed_s,PG,isym,bond_sym,angle_sym,dihed_sym,unitZ)
!
! subroutine internal_Wilson(molec,S,ModeDef,B,G,Asel,verbose)
! 
! subroutine gf_method(Hess,molec,ModeDef,L,B,G,Freq,Asel,X,Xinv,verbose)
! 
! subroutine axis_swithching(molec,molec2,T0)
! 
! subroutine zmat2cart(molec,bond_s,angle_s,dihed_s,S,verbose)
! 
! subroutine addcart(atom_1,atom_2,atom_3,&
!     subroutine rotation_3D(vx,vy,vz,tx,ty,tz,Theta) 
! 
! subroutine addcart(atom_1,atom_2,atom_3,bond,angle,dihed,atom_New)
!
! subroutine check_ori2(molec,molec2)
!
! subroutine check_Ci(molec,molecP,molec2)
!
! History
! v0.4.1: Add support to non-redundat GF method

 contains

subroutine internal_Wilson(molec,S,S_sym,ModeDef,B,G,Asel,verbose)

    !This routine needs to be put in a real SR shape!
    !Performs vibrational analysis using ghe Wilson B matrix and the GF method.
    !Returns L matrix (in Aux) and B matrix (in B)

    use structure_types
    use line_preprocess
    use alerts
    use constants
    use atomic_geom
    use MatrixMod

    implicit none

    integer,parameter :: NDIM = 800
    real(8),parameter :: ZEROp = 1.d-10 !practically zero

    !====================== 
    !ARGUMENTS
    !Input
    type(str_resmol),intent(in) :: molec            ! Input molecule to which the vib analysis is performed
    logical,intent(in) :: verbose                   ! Control Output details
    integer,dimension(NDIM),intent(in) :: S_sym     ! Manages symmetric internals
    !Input-Output
    real(8),dimension(NDIM,NDIM),intent(inout) :: Asel !Used to manage nonredundant internal variables (disabled option)
    !Output
    real(8),dimension(NDIM,NDIM),intent(out) :: B, G          !B -Wilson- and G matrices
    real(8),dimension(NDIM),intent(out) :: S                  !Vector of internal coordinates
    character(len=100),dimension(NDIM),intent(out) :: ModeDef !Definition of normal modes
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=100) :: dummy_char
    !====================== 

    !====================== 
    !System info
    integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
    integer :: nbonds, ndihed, nangles
    integer :: Nat, Nvib, Nred
    character(len=5) :: PG
    !====================== 

    !====================== 
    !INTERNAL ANALYSIS
    !B and G matrices
!     real(8),dimension(NDIM,NDIM),intent(out) :: B, G
    !AUXILIAR MATRICES
    real(8),dimension(NDIM,NDIM) :: AuxT,Aux, Aux3
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
    real(8) :: mu
    real(8),dimension(NDIM) :: Vec, Vec2
    !====================== 

    !====================== 
    !Auxiliars for LAPACK matrix nversion
    integer :: info
    integer,dimension(NDIM) :: ipiv
    real(8),dimension(NDIM,NDIM) :: work
    !====================== 

    !=============
    !Counters
    integer :: i,j,k, ii,jj,kk,kkk, iat
    integer :: i_1, i_2, i_3, i_4, imax, jmax
    !=============

! (End of variables declaration) 
!==================================================================================


    !Set bonded
    nbonds  = molec%geom%nbonds
    nangles = molec%geom%nangles
    ndihed  = molec%geom%ndihed
    bond_s(1:nbonds,1:2)  =  molec%geom%bond(1:nbonds,1:2)
    angle_s(1:nangles,1:3) =  molec%geom%angle(1:nangles,1:3)
    dihed_s(1:ndihed,1:4)  =  molec%geom%dihed(1:ndihed,1:4)


    !Initialize matrices
    Nat = molec%natoms
    Nvib = 3*Nat-6
    B(1:Nvib,1:3*Nat) = 0.d0
    G(1:Nvib,1:Nvib)  = 0.d0

!=============================================================
!   Ref: Decius, Cross and Wilson (Section 4.2) --
!    The same nomenclature is used (index refer to the same
!    atoms, even if reversed)
!==============================================================

    !k-index runs over internal coordinates
    k=0
    write(6,'(/,A,I3)') "BONDS", nbonds !molec%natoms-1
!     do i=2,molec%natoms
    do i=1,nbonds
        k=k+1

        i_1 = bond_s(i,1)
        i_2 = bond_s(i,2)
        r21 = calc_dist(molec%atom(i_1),molec%atom(i_2))
        S(k) = r21
        write(6,'(I3,A,X,A,2(I3,A),2F15.8)') k, ".-",&
                            trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                            trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                            r21*BOHRtoAMS, r21
        write(ModeDef(k),'(A,2(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
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

    write(6,'(/,A,I3)') "ANGLES", nangles! molec%natoms-2
!     do i=3,molec%natoms
    do i=1,nangles
        k=k+1

        i_1 = angle_s(i,1)
        i_3 = angle_s(i,2)
        i_2 = angle_s(i,3)
        ang1 = calc_angle(molec%atom(i_1),molec%atom(i_3),molec%atom(i_2))
        S(k) = ang1
        write(6,'(I3,A,X,A,3(I3,A),F15.8)') k, ".-",& 
                            trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                            trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                            trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                            ang1*360.d0/2.d0/pi
        write(ModeDef(k),'(A,3(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
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

    write(6,'(/,A,I3)') "DIHEDRALS", ndihed!molec%natoms-3
!     do i=4,molec%natoms
    do i=1,ndihed
        k=k+1

        i_1 = dihed_s(i,1)
        i_2 = dihed_s(i,2)
        i_3 = dihed_s(i,3)
        i_4 = dihed_s(i,4)
        ang1 = calc_dihed_new(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
        S(k) = ang1
        write(6,'(I3,A,X,A,4(I3,A),F15.8)') k, ".- ",&
                            trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                            trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                            trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                            trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")", &
                            ang1*360.d0/2.d0/pi
        write(ModeDef(k),'(A,4(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                            trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                            trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                            trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")"

        !Four cart displacements different from zero (some index intercheged with Decius..)
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
        B(k,3*i_1-2) =  e21Pe32x/(r21*dsin(ang2)**2)
        B(k,3*i_1-1) =  e21Pe32y/(r21*dsin(ang2)**2)
        B(k,3*i_1  ) =  e21Pe32z/(r21*dsin(ang2)**2)
        !s2 
        B(k,3*i_2-2) = -((r32-r21*dcos(ang2))*e21Pe32x/(r32*r21*dsin(ang2)**2) &
                     +  dcos(ang3)*e43Pe32x/(r32*dsin(ang3)**2))
        B(k,3*i_2-1) = -((r32-r21*dcos(ang2))*e21Pe32y/(r32*r21*dsin(ang2)**2) &
                     +  dcos(ang3)*e43Pe32y/(r32*dsin(ang3)**2))
        B(k,3*i_2  ) = -((r32-r21*dcos(ang2))*e21Pe32z/(r32*r21*dsin(ang2)**2) &
                     +  dcos(ang3)*e43Pe32z/(r32*dsin(ang3)**2))
        !Following Decius paper (instead of W-D-C book formulation)
!         B(k,3*i_2-2) = ( 1.d0/(r21*dsin(ang2)) - (dcos(ang1)/dtan(ang3) + &
!                          1.d0/dtan(ang2))/r32 )*e21Pe32x/dsin(ang2) +     &
!                          dsin(ang1)/dtan(ang3)/dsin(ang2)/r32*e32Pe21Pe32x
!         B(k,3*i_2-1) = ( 1.d0/(r21*dsin(ang2)) - (dcos(ang1)/dtan(ang3) + &
!                          1.d0/dtan(ang2))/r32 )*e21Pe32y/dsin(ang2) +     &
!                          dsin(ang1)/dtan(ang3)/dsin(ang2)/r32*e32Pe21Pe32y
!         B(k,3*i_2  ) = ( 1.d0/(r21*dsin(ang2)) - (dcos(ang1)/dtan(ang3) + &
!                          1.d0/dtan(ang2))/r32 )*e21Pe32z/dsin(ang2) +     &
!                          dsin(ang1)/dtan(ang3)/dsin(ang2)/r32*e32Pe21Pe32z
        !s3
        B(k,3*i_3-2) = -((r32-r43*dcos(ang3))*e43Pe32x/(r32*r43*dsin(ang3)**2) &
                     +  dcos(ang2)*e21Pe32x/(r32*dsin(ang2)**2))
        B(k,3*i_3-1) = -((r32-r43*dcos(ang3))*e43Pe32y/(r32*r43*dsin(ang3)**2) &
                     +  dcos(ang2)*e21Pe32y/(r32*dsin(ang2)**2))
        B(k,3*i_3  ) = -((r32-r43*dcos(ang3))*e43Pe32z/(r32*r43*dsin(ang3)**2) &
                     +  dcos(ang2)*e21Pe32z/(r32*dsin(ang2)**2))
        !Following Decius paper (instead of W-D-C book formulation)
!         B(k,3*i_3-2) = ( 1.d0/(r43*dsin(ang3)) - (dcos(ang1)/dtan(ang2) + &
!                          1.d0/dtan(ang3))/r32 )*e43Pe32x/dsin(ang3) -     &
!                          dsin(ang1)/dtan(ang2)/dsin(ang3)/r32*e32Pe43Pe32x
!         B(k,3*i_3-1) = ( 1.d0/(r43*dsin(ang3)) - (dcos(ang1)/dtan(ang2) + &
!                          1.d0/dtan(ang3))/r32 )*e43Pe32y/dsin(ang3) -     &
!                          dsin(ang1)/dtan(ang2)/dsin(ang3)/r32*e32Pe43Pe32y
!         B(k,3*i_3  ) = ( 1.d0/(r43*dsin(ang3)) - (dcos(ang1)/dtan(ang2) + &
!                          1.d0/dtan(ang3))/r32 )*e43Pe32z/dsin(ang3) -     &
!                          dsin(ang1)/dtan(ang2)/dsin(ang3)/r32*e32Pe43Pe32z
        !s4
        B(k,3*i_4-2) =  e43Pe32x/(r43*dsin(ang3)**2)
        B(k,3*i_4-1) =  e43Pe32y/(r43*dsin(ang3)**2)
        B(k,3*i_4  ) =  e43Pe32z/(r43*dsin(ang3)**2)

    enddo

!     write(6,'(/,A,I3)') "IMPROPERS", ndihed!molec%natoms-3
! !     do i=4,molec%natoms
!     do i=1,ndihed
!         k=k+1
! 
!         i_1 = dihed_s(i,1)
!         i_2 = dihed_s(i,2)
!         i_3 = dihed_s(i,3)
!         i_4 = dihed_s(i,4)
!         ang1 = calc_dihed_new(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
!         S(k) = ang1
!         write(6,'(A,4(I3,A),F15.8)') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
!                             trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
!                             trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
!                             trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")", &
!                             ang1*360.d0/2.d0/pi
!         write(ModeDef(k),'(A,4(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
!                             trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
!                             trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
!                             trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")"
! 
!         !Four cart displacements different from zero (some index intercheged with Decius..)
!         ang2 = calc_angle(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3))
!         ang3 = calc_angle(molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
!         r21  = calc_dist(molec%atom(i_1),molec%atom(i_2))
!         e21x = (molec%atom(i_1)%x-molec%atom(i_2)%x)/r21
!         e21y = (molec%atom(i_1)%y-molec%atom(i_2)%y)/r21
!         e21z = (molec%atom(i_1)%z-molec%atom(i_2)%z)/r21
!         r32  = calc_dist(molec%atom(i_2),molec%atom(i_3))
!         e32x = (molec%atom(i_2)%x-molec%atom(i_3)%x)/r32
!         e32y = (molec%atom(i_2)%y-molec%atom(i_3)%y)/r32
!         e32z = (molec%atom(i_2)%z-molec%atom(i_3)%z)/r32
!         r43  = calc_dist(molec%atom(i_3),molec%atom(i_4))
!         e43x = (molec%atom(i_3)%x-molec%atom(i_4)%x)/r43
!         e43y = (molec%atom(i_3)%y-molec%atom(i_4)%y)/r43
!         e43z = (molec%atom(i_3)%z-molec%atom(i_4)%z)/r43
!         e21Pe32x=e21y*e32z-e21z*e32y
!         e21Pe32y=e21z*e32x-e21x*e32z
!         e21Pe32z=e21x*e32y-e21y*e32x
!         e32Pe21Pe32x=e32y*e21Pe32z-e32z*e21Pe32y
!         e32Pe21Pe32y=e32z*e21Pe32x-e32x*e21Pe32z
!         e32Pe21Pe32z=e32x*e21Pe32y-e32y*e21Pe32x
!         e43Pe32x=e43y*e32z-e43z*e32y
!         e43Pe32y=e43z*e32x-e43x*e32z
!         e43Pe32z=e43x*e32y-e43y*e32x
!         e32Pe43Pe32x=e32y*e43Pe32z-e32z*e43Pe32y
!         e32Pe43Pe32y=e32z*e43Pe32x-e32x*e43Pe32z
!         e32Pe43Pe32z=e32x*e43Pe32y-e32y*e43Pe32x
! 
!         !s1
!         
! 
!     enddo

    Nred = k
print*, "Actual number of vibrational coordinates", Nred

    if (verbose) then
    write(6,*) ""
    write(6,*) "B MATRIX-"
    do i=1,Nred
        write(6,'(100(F9.5,X))') B(i,1:3*Nat)
    enddo
    print*, ""
    print*, "B MATRIX (all elements)"
    do i=1,Nred
        do j=1,3*Nat
            print'(G15.8)', B(i,j)
        enddo
    enddo
    print*, "--END OF B MATRIX (all elements)"
    endif

    !Symmetry adapted coordinates
!     Aux(1,1:3*Nat) = (B(1,1:3*Nat)+B(6,1:3*Nat))
!     Aux(2,1:3*Nat) = (B(1,1:3*Nat)-B(6,1:3*Nat))
!     B(1,1:3*Nat) = Aux(1,1:3*Nat)
!     B(6,1:3*Nat) = Aux(2,1:3*Nat)
!     Aux(1,1:3*Nat) = (B(2,1:3*Nat)+B(3,1:3*Nat))
!     Aux(2,1:3*Nat) = (B(2,1:3*Nat)-B(3,1:3*Nat))
!     B(2,1:3*Nat) = Aux(1,1:3*Nat)
!     B(3,1:3*Nat) = Aux(2,1:3*Nat)
!     Aux(1,1:3*Nat) = (B(4,1:3*Nat)+B(5,1:3*Nat))
!     Aux(2,1:3*Nat) = (B(4,1:3*Nat)-B(5,1:3*Nat))
!     B(4,1:3*Nat) = Aux(1,1:3*Nat)
!     B(5,1:3*Nat) = Aux(2,1:3*Nat)
    if (S_sym(3*Nat) /= 0) then
         print*, ""
         print*, "Symmetry addapted coordinates will be used"
         do i=1,Nvib
             if (S_sym(i) <= i) cycle
             j=S_sym(i)
             Aux(1,1:3*Nat) = (B(i,1:3*Nat)+B(j,1:3*Nat))
             Aux(2,1:3*Nat) = (B(i,1:3*Nat)-B(j,1:3*Nat))
             B(i,1:3*Nat) = Aux(1,1:3*Nat)
             B(j,1:3*Nat) = Aux(2,1:3*Nat)
             dummy_char = trim(adjustl(ModeDef(i)))//"+"//trim(adjustl(ModeDef(j)))
             ModeDef(j) = trim(adjustl(ModeDef(i)))//"-"//trim(adjustl(ModeDef(j)))
             ModeDef(i) = trim(adjustl(dummy_char))
         enddo
    endif


    !CONCTRUCT G
    do i=1,Nred
        do j=1,Nred
            G(i,j) = 0.d0
            k=0
            do kk=1,Nat
            do iat=1,3
                k=k+1
                mu = 1.d0/molec%atom(kk)%mass/UMAtoAU
                G(i,j) = G(i,j) + mu*B(i,k)*B(j,k)
            enddo 
            enddo
        enddo
    enddo

    ! Check symmetry
    mu = 0.d0
    jmax=1
    imax=1
    do i=1,Nred
        do j=i,Nred
            Aux(i,j) = (G(i,j)-G(j,i))/G(i,j) *100.d0
            mu = max(mu,dabs(Aux(i,j)))
            if (mu == Aux(i,j)) then
                imax = i
                jmax = j
            endif
        enddo
    enddo
    if (verbose) then
    write(6,*) ""
    print*, "Maximum difference in symmetric matrix (%):", mu, "for", imax,jmax
    print*, G(imax,jmax), G(jmax,imax)
    write(6,*) ""

    write(6,*) ""
    write(6,*) "G MATRIX"
    do i=1,Nred
        write(6,'(100(E10.2,3X))') G(i,1:Nred)
    enddo
    print*, ""
    print*, "G MATRIX (lower triangular)"
    do i=1,Nred
        do j=1,i
            print'(G15.8)', G(i,j)
        enddo
    enddo
    print*, "--END OF G MATRIX (lower triangular)"
    endif

    print*, "Redundant coordinates", Nred-Nvib

    !IF WE USE ALL BONDED PARAMETERS,WE HAVE REDUNDANCY. WE SELECT A NON-REDUNDANT
    !COMBINATION FROM THE NON-ZERO EIGENVALUES OF G (Reimers 2001, JCP)
    !--- this is not tested for new versions and must not be used ---
    if (Nred==-9999) then !if (Nred-Nvib /= 0) then    
!         if (Asel(1,1) == 99.d0) then
            Asel = 0.d0
            !Get a non-redundant set from the non-zero eigenvalues of G
            call diagonalize_full(G(1:Nred,1:Nred),Nred,Aux(1:Nred,1:Nred),Vec(1:Nred),"lapack")
            if (verbose) then
            write(6,*) ""
            write(6,*) "A MATRIX before reordering"
            do i=1,Nred
                write(6,'(100(F9.5,X))') Aux(i,1:Nred)
            enddo
            write(6,*) "Eigenvalues"
            do i=1,Nred
                write(6,'(100(F9.5,X))') Vec(i)
            enddo
            endif

            kk=0 
            kkk=0
            do k=1,Nred
                if (dabs(Vec(k)) > ZEROp) then
                    kk=kk+1
                    !Does the diagonalization return eigenvectos in columns of rows????
                    ! should be columns as: D = P^-1 A P
                    Asel(1:Nred,kk) = Aux(1:Nred,k)
                    Vec2(kk) = Vec(k)
                    ! .. for symmetric A matrices, either rows or columns can be selected as eigenvectors anyway
                else
                    !Redudant eigenvectors
                    kkk=kkk+1
                    Asel(1:Nred,Nvib+kkk) = Aux(1:Nred,k)
                    Vec2(kk) = Vec(k)
                endif
            enddo
            print*, "Check Nred Nvib+kkk", Nred, Nvib+kkk
            if (kk /= Nvib) then
                call sort_vec(Vec,Nred)
                do i=1,Nred
                    print*, Vec(i)
                enddo
                print*, "Error", kk, Nvib
                stop
            endif
            print*, "New set of non-redundant internal coordinates", kk
!         else
!             print'(/,X,A,/)', "Using A determined for other state."
! ! print*, "Asel(1,1)", Asel(1,1)
!         endif


        !Definitions of the modes are no longer valid
        ModeDef(:) = ""

        if (verbose) then
        write(6,*) ""
        write(6,*) "A MATRIX (redundant eigenvectors)"
        do i=1,Nred
            write(6,'(100(F9.5,X))') Asel(i,1:Nred)
        enddo
        write(6,*) "Eigenvalues"
        do i=1,Nred
            write(6,'(100(F9.5,X))') Vec2(i)
        enddo
        endif

        !Check diagonalization
        do i=1,Nred
            Aux(i,1:Nred) = Vec2(i)*Asel(1:Nred,i)
        enddo
        Aux(1:Nred,1:Nred) = matmul(Asel(1:Nred,1:Nred),Aux(1:Nred,1:Nred))
        
!         write(6,*) ""
!         write(6,*) "G MATRIX"
!         do i=1,Nred
!             write(6,'(100(F9.5,X))') G(1:Nred,i)
!         enddo
!         write(6,*) ""
!         write(6,*) "G MATRIX from diag"
!         do i=1,Nred
!             write(6,'(100(F9.5,X))') Aux(1:Nred,i)
!         enddo

        !B' = A^T B
        do i=1,Nvib
            do j=1,3*Nat
                Aux(i,j) = 0.d0
                do k=1,Nred
                    Aux(i,j) = Aux(i,j) + Asel(k,i)*B(k,j)
                enddo
            enddo
        enddo
        B(1:Nvib,1:3*Nat)=Aux(1:Nvib,1:3*Nat)
!         !Trasform coordinates -- need validation!!
!         do i=1,Nvib
!             Vec(i) = 0.d0
!             do k=1,Nred
!                 Vec(i) = Vec(i) + Asel(i,k)*S(k)
!             enddo
!         enddo
!         S(1:Nvib) = Vec(1:Nvib)

        if (verbose) then
        write(6,*) ""
        write(6,*) "B MATRIX FOR NON REDUNDANT SET"
        do i=1,Nvib
            write(6,'(100(F9.5,X))') Aux(i,1:3*Nat)
        enddo
        endif


        !CONCTRUCT G
        do i=1,Nvib
            do j=1,Nvib
                G(i,j) = 0.d0
                k=0
                do kk=1,Nat
                do iat=1,3
                    k=k+1
                    mu = 1.d0/molec%atom(kk)%mass/UMAtoAU
                    G(i,j) = G(i,j) + mu*B(i,k)*B(j,k)
                enddo 
                enddo
            enddo
        enddo

        if (verbose) then
        write(6,*) ""
        write(6,*) "G MATRIX"
        do i=1,Nvib
            write(6,'(100(E10.3,X))') G(i,1:Nvib)
        enddo
        endif

    endif

    return

end subroutine internal_Wilson

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine gf_method(Hess,molec,S_sym,ModeDef,L,B,G,Freq,Asel,X,Xinv,verbose)

    !Performs vibrational analysis using ghe Wilson B matrix and the GF method.
    !Returns L matrix (in Aux) and B matrix (in B)

    use structure_types
    use line_preprocess
    use alerts
    use constants
    use atomic_geom
    use MatrixMod

    implicit none

    integer,parameter :: NDIM = 800

    !====================== 
    !ARGUMENTS
    !Inpu
    type(str_resmol),intent(in) :: molec                     !Input molecule
    real(8),dimension(1:NDIM,1:NDIM),intent(in) ::    G,B    !Wilson matrices
    character(len=100),dimension(NDIM),intent(in) :: ModeDef !Definition of modes
    integer,dimension(NDIM),intent(in) :: S_sym              !Manages symmetric internals
    logical,intent(in) :: verbose                            !Verbose level
    !Input-Output
    real(8),dimension(1:NDIM,1:NDIM),intent(inout) :: Hess   !Hessian: cart(in)-intern(out)
    real(8),dimension(NDIM,NDIM),intent(inout) :: Asel       !Manage nonredundant (disabled)
    !Output
    real(8),dimension(NDIM,NDIM),intent(out) :: L,X,Xinv     !Normal modes, G^1/2, G^-1/2 (internal)
    !====================== 

    !====================== 
    !System variables
    integer :: Nat, Nvib
    character(len=3) :: ModeSymm
    !====================== 

    !====================== 
    !Internal analysis 
    real(8),dimension(1:NDIM,1:NDIM) :: Ginv
    !Auxiliar matrices
    real(8),dimension(1:NDIM,1:NDIM) :: AuxT,Aux,Aux3
    !Freq
    real(8),dimension(NDIM),intent(out) :: Freq
    !====================== 

    !====================== 
    !Auxiliars for LAPACK matrix nversion
    real(8),dimension(NDIM,NDIM) :: work
    integer,dimension(NDIM) :: ipiv, ipiv2
    integer :: info
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: ierr
    character(1) :: null
    character(len=16) :: dummy_char
    real(8) :: Theta, Theta2, mu
    !TESTING ARRAYS (auxiliar)
    real(8),dimension(NDIM,NDIM) :: test1,test2,test3
    !====================== 

    !=============
    !Counters
    integer :: i,j,k, ii,jj,kk, iat
    !=============

! (End of variables declaration) 
!==================================================================================

    Nat = molec%natoms
    Nvib = molec%geom%nbonds + &!3*Nat-6
           molec%geom%nangles+ &
           molec%geom%ndihed
print*, "Number of vibrational coordinates", Nvib

    !=====================================================================
    ! HESSIAN IN INTERNAL COORDINATES (JCC, 17, 49-56, by Frisch et al)
    !======================================================================
    !From JCC, 17, 49-56,:
     ! Hint = G^- Bu(Hx+B'^Tg_q)u^TB^T G^-
     !g_q is the gradient, so g_q=0
     !G^- is the generalized inverse (for redundant internal) or simply the
     !inverse for nonredundant

    !Inverse of G (Generalized inverse)
    call diagonalize_full(G(1:Nvib,1:Nvib),Nvib,Aux(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
    if (verbose) then
        write(6,*) ""
        write(6,*) "G eigenvalues"
        do i=1,Nvib
            write(6,'(100(E10.2,3X))') Freq(i)
        enddo
    endif

    !Set low values to zero
    ii=0
    do i=1,Nvib
        if (dabs(Freq(i)) < 1.d-15) then
            Freq(i)=0.d0
        else
            ii=ii+1            
        endif
    enddo
    print*, "Nvib:", 3*Nat-6
    print*, "Non-zero G eigenvalues", ii
    if (3*Nat-6 /= ii) then
        print*, "Error: number of G eigenvalues does not match vibrational degrees of freedom"
        stop
    endif

!     do i=1,Nvib
!         do j=1,Nvib
!             Ginv(i,j) = 0.d0
!             do k=1,Nvib
!                 if (Freq(k)/=0.d0) Ginv(i,j) = Ginv(i,j) + Aux(i,k)/dsqrt(Freq(k))*Aux(j,k)
!             enddo
!         enddo
!     enddo

    Ginv=0.d0
    do i=1,Nvib
        if (Freq(i) /= 0.d0) Ginv(i,i) = 1.d0/Freq(i)
    enddo
    do i=1,Nvib
        do j=1,Nvib
            AuxT(i,j) = 0.d0
            do k=1,Nvib
                AuxT(i,j) = AuxT(i,j) + Aux(i,k)*Ginv(k,j)                
            enddo
        enddo
    enddo
    do i=1,Nvib
        do j=1,Nvib
            Ginv(i,j) = 0.d0
            do k=1,Nvib
                Ginv(i,j) = Ginv(i,j) + AuxT(i,k)*Aux(j,k)                
            enddo
        enddo
    enddo

    if (verbose) then
        write(6,*) ""
        write(6,*) "Ginv MATRIX"
        do i=1,Nvib
            write(6,'(100(E10.2,3X))') Ginv(i,1:Nvib)
        enddo
    endif


!     Ginv(1:Nvib,1:Nvib)=G(1:Nvib,1:Nvib)
!     call dsytrf('L', Nvib, Ginv, NDIM, ipiv, work, NDIM, info)
!     call dsytri('L', Nvib, Ginv, NDIM, ipiv, work, info)
! !   RECONSTRUCT UPPER PART
!     do i=1,Nvib
!     do j=i+1,Nvib
!         Ginv(i,j) = Ginv(j,i)
!     enddo
!     enddo

! stop

    !Compute G^-1Bu  (where u is the inverse mass matrix)
    do i=1,Nvib
        j=0
        do jj=1,Nat
        do iat=1,3
            j=j+1
            AuxT(i,j) = 0.d0
            mu=1.d0/molec%atom(jj)%mass/UMAtoAU
            do k=1,Nvib
                AuxT(i,j) = AuxT(i,j) + Ginv(i,k)*B(k,j)                
            enddo
            AuxT(i,j) = AuxT(i,j)*mu
        enddo
        enddo
    enddo

    !Now: Hint = Aux Hcart Aux^T
    do i=1,Nvib
        do j=1,3*Nat
            work(i,j) = 0.d0
            do k=1,3*Nat
                work(i,j) = work(i,j) + AuxT(i,k)*Hess(k,j)                
            enddo
        enddo
    enddo
    do i=1,Nvib
        do j=1,3*Nat
            Hess(i,j) = 0.d0
            do k=1,3*Nat
                Hess(i,j) = Hess(i,j) + work(i,k)*AuxT(j,k)                
            enddo
        enddo
    enddo

    if (verbose) then
    print*, ""
    print*, "F MATRIX"
    do i=1,Nvib
        print'(100(F8.3,X))', Hess(i,1:Nvib)
    enddo 
    print*, ""
    print*, "F MATRIX (lower triangular)"
    do i=1,Nvib
        do j=1,i
            print'(G15.8)', Hess(i,j)
        enddo
    enddo
    print*, "--END OF F MATRIX (lower triangular)"
    endif
!     test1(1:Nvib,1:Nvib) = Hess(1:Nvib,1:Nvib)


    !FG PROCEDURE (Wilson, Decius and Cross, Section 4.3)
    ! |GF - \lamda| = 0  <==> GFL = L\lamda   means the diagonalize of GF, but is not symmetric
    ! Solutions
    ! *Orthogonalization of FL = G^-1L\lamda similar to orthogonalization fo S in Roothan. 
    ! *Symmetrize GF and diagonalize following Mizayawa (Sect. 4.6, Califano, Vibrational States)

    ! FOLLOWING FIRST ALTERNATIVE (similar to diagonalization of basis set)
    !Note this diagonalization was performed to get the generalized inverse: G^- (NOT REPEATED, THEN)
!     call diagonalize_full(G(1:Nvib,1:Nvib),Nvib,Aux(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")

    !----------------------------------------------------------------
    !REFLEXIONES DEL PROGRAMADOR
    !The internal coordinates also need to be rotated. But using G^-1/2? Para que??
    ! No, they should not be changed since the ones finally used are the "original" ones, not this orthogonalized set
    !----------------------------------------------------------------

    if (verbose) then
    print*, ""
    print*, "EigenVaules="
    print*, Freq(1:Nvib)
    endif

    !The non-unitary rotation which orthogonalize G^-1 could be X=G^(1/2)  (inverse compared with the case of Roothan eqs)
    ! where G^{1/2} = Ug^{1/2}U^T
    !Store the rotation (X) in AuxT    !G. Note X=X^T
!     X=0.d0
!     do i=1,Nvib
!         if (Freq(i) /= 0.d0) X(i,i) = Freq(i)
!     enddo
!     do i=1,Nvib
!         do j=1,Nvib
!             AuxT(i,j) = 0.d0
!             do k=1,Nvib
!                 AuxT(i,j) = AuxT(i,j) + Aux(i,k)*X(k,j)                
!             enddo
!         enddo
!     enddo
!     do i=1,Nvib
!         do j=1,Nvib
!             X(i,j) = 0.d0
!             do k=1,Nvib
!                 X(i,j) = X(i,j) + AuxT(i,k)*Aux(j,k)                
!             enddo
!         enddo
!     enddo

    do i=1,Nvib
        do j=1,Nvib
            X(i,j) = 0.d0
            do k=1,Nvib
                X(i,j) = X(i,j) + Aux(i,k)*dsqrt(Freq(k))*Aux(j,k)
            enddo
        enddo
    enddo
    !And the inverse
!     Xinv=0.d0
!     do i=1,Nvib
!         if (Freq(i) /= 0.d0) Xinv(i,i) = 1.d0/Freq(i)
!     enddo
!     do i=1,Nvib
!         do j=1,Nvib
!             AuxT(i,j) = 0.d0
!             do k=1,Nvib
!                 AuxT(i,j) = AuxT(i,j) + Aux(i,k)*Xinv(k,j)                
!             enddo
!         enddo
!     enddo
!     do i=1,Nvib
!         do j=1,Nvib
!             Xinv(i,j) = 0.d0
!             do k=1,Nvib
!                 Xinv(i,j) = Xinv(i,j) + AuxT(i,k)*Aux(j,k)                
!             enddo
!         enddo
!     enddo
    do i=1,Nvib
        do j=1,Nvib
            Xinv(i,j) = 0.d0
            do k=1,Nvib
                if (Freq(k)/=0.d0) Xinv(i,j) = Xinv(i,j) + Aux(i,k)/dsqrt(Freq(k))*Aux(j,k)
            enddo
        enddo
    enddo

    if (verbose) then
    print*, ""
    print*, "X="
    do i=1,Nvib
        print'(100(F10.3,2X))', X(i,1:Nvib)
    enddo
    endif

    !Now rotate F, F'=X^TFX. Store the rotated matrix in F
    do i=1,Nvib
        do j=1,Nvib
            Aux(i,j) = 0.d0
            do k=1,Nvib
                Aux(i,j) = Aux(i,j) + X(k,i)*Hess(k,j)
            enddo
        enddo
    enddo
    do i=1,Nvib
        do j=1,Nvib
            Hess(i,j) = 0.d0
            do k=1,Nvib
                Hess(i,j) = Hess(i,j) + Aux(i,k)*X(k,j)
            enddo
        enddo
    enddo

    !We can now diagonalize F
    call diagonalize_full(Hess(1:Nvib,1:Nvib),Nvib,L(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
    !WE NEED TO PASS "L" AS IT IS NOW

    !Check freqcuencies
    print*, ""
    print*, "FORCE CONSTANTS (A.U.)"
!     call sort_vec(Freq,Nvib)
    do i=1,Nvib
          write(6,*) Freq(i)
    enddo

    !Check freqcuencies
    print*, ""
    print*, "FREQUENCIES (cm^-1)"
!     call sort_vec(Freq,Nvib)
    do i=1,Nvib
          Freq(i) = sign(dsqrt(abs(Freq(i))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
                         Freq(i))
          write(6,*) Freq(i)
    enddo
!     print*, "L'="
!     do i=1,Nvib
!         print'(100(F10.3,2X))', Aux(i,1:Nvib)
!     enddo
    if (verbose) then
    print*, "L'="
    do i=1,Nvib
        print'(100(F8.3,2X))', L(i,1:Nvib)
    enddo
    endif

    if (verbose) then
    do i=1,Nvib
        do j=1,Nvib
            Aux(i,j) = 0.d0
            do k=1,Nvib 
                Aux(i,j) = Aux(i,j) + L(i,k)*L(j,k)
            enddo
        enddo
    enddo
    print*, "L'L'^t="
    do i=1,Nvib
        print'(100(F8.3,2X))', Aux(i,1:Nvib)
    enddo
    endif

    ! The L' matrix is given in the an orthogonalized basis. 
    ! To have a description based in the initial internal, restore in the original basis:
    ! L = XL'   (with X the orhogonalizing matrix, stored in X).
      do i=1,Nvib
          do j=1,Nvib
              Aux3(i,j) = 0.d0
              do k=1,Nvib 
                  Aux3(i,j) = Aux3(i,j) + X(i,k)*L(k,j) !Que lio de exponenetes!!
              enddo
          enddo
      enddo

    if (verbose) then
    print*, "L="
    do i=1,Nvib
        print'(100(F8.3,2X))', Aux3(i,1:Nvib)
    enddo
    endif

    !We are going to pass the L matrix in the original internal frame
    L(1:Nvib,1:Nvib) = Aux3(1:Nvib,1:Nvib)

    if (verbose) then
        !L^t
        do i=1,Nvib
        do j=1,Nvib
            test3(i,j) = Aux3(j,i)
        enddo
        enddo
        test2(1:Nvib,1:Nvib) = matmul(Aux3(1:Nvib,1:Nvib),test3(1:Nvib,1:Nvib))
        print*, "L L^t="
        do i=1,Nvib
            print'(100(F8.3,2X))', test2(i,1:Nvib)
        enddo
        test2(1:Nvib,1:Nvib) = matmul(test3(1:Nvib,1:Nvib),Aux3(1:Nvib,1:Nvib))
        print*, "L^t L="
        do i=1,Nvib
            print'(100(F8.3,2X))', test2(i,1:Nvib)
        enddo
        !L^tFL
        test2(1:Nvib,1:Nvib) = matmul(test3(1:Nvib,1:Nvib),test1(1:Nvib,1:Nvib))
        test2(1:Nvib,1:Nvib) = matmul(test2(1:Nvib,1:Nvib),Aux3(1:Nvib,1:Nvib))
        test2(1:Nvib,1:Nvib) = dsqrt(dabs(test2(1:Nvib,1:Nvib))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2
        print*, "L^tFL="
        do i=1,Nvib
            print'(100(F8.3,2X))', test2(i,1:Nvib)
        enddo
        !L^G^-1L
        test2(1:Nvib,1:Nvib) = matmul(test3(1:Nvib,1:Nvib),Ginv(1:Nvib,1:Nvib))
        test2(1:Nvib,1:Nvib) = matmul(test2(1:Nvib,1:Nvib),Aux3(1:Nvib,1:Nvib))
        print*, "L^tG^-1L="
        do i=1,Nvib
            print'(100(F8.3,2X))', test2(i,1:Nvib)
        enddo
    endif

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
!           Aux3(1,1:Nvib) = Aux3(1,1:Nvib)/dsqrt(Theta)
          call sort_vec_max(Aux3(1,1:Nvib),ipiv(1:Nvib),Nvib)
          call sort_vec_max(Aux(1,1:Nvib),ipiv2(1:Nvib),Nvib)

          !Determine symmetry
          do j=1,Nvib
              jj=ipiv(j)
              if (S_sym(jj) == jj) cycle
              ii = S_sym(jj)
              Theta = L(jj,i)/L(ii,i)
              if (Theta > 0.d0) then 
                  ModeSymm="A' "
              else
                  ModeSymm="A''"
              endif
              exit
              !Use stretching coordinates to set symmetry
              if (ii<Nat) exit
          enddo

          print*, ""
          print'(A,I4,6X,A,F8.3,4X,A,A,F8.3)', "Mode ", i, " Freq(cm^-1) = ", Freq(i), "Symm ", ModeSymm, Theta
          print*, "      S        Coef.     Contrib.(%)  ContribCorr(%)          Description"
          print*, " ======================================================================="
          Theta = 0.d0
          kk=0
          do j=1,Nvib
              if (Theta > 0.9d0) exit
              jj = ipiv(j)
              Theta = Theta + Aux3(1,j)!**2
              print'(5x,i4,4x,g10.3,4x,2(f8.3,4x),a)', jj, L(jj,i), Aux3(1,j)*100, Aux(1,j)*100, trim(adjustl(ModeDef(jj)))
              kk=kk+1
          enddo 
          print*, " ========================================================================"
          write(6,'(A,I3)') "     Total Number of internal to describe >90% of the mode: ", kk
      enddo


    !The former trasformation related S in terms of Q's. To obtain the inverse relation:
    !Inverse of L 
    Aux(1:Nvib,1:Nvib)=L(1:Nvib,1:Nvib)
    call dgetrf(Nvib,Nvib, Aux, NDIM, ipiv, info)
    call dgetri(Nvib, Aux, NDIM, ipiv, work, NDIM, info)
!           call zgetrf(N36, N36, Aux, mad, ipiv, info)
!           call zgetri(N36, Aux, mad, ipiv, work, mad, info)
!     Aux(1:Nvib,1:Nvib)=L(1:Nvib,1:Nvib)

      !Compute LL^-1
      do i=1,Nvib
          do j=1,Nvib
              AuxT(i,j) = 0.d0
              do k=1,Nvib 
                  AuxT(i,j) = AuxT(i,j) + Aux3(i,k)*Aux(k,j)
              enddo
          enddo
      enddo

    if (verbose) then
        print*, "L^-1="
        do i=1,Nvib
            print'(100(F8.3,2X))', Aux(i,1:Nvib)
        enddo

        print*, "L L^-1="
        do i=1,Nvib
            print'(100(F8.3,2X))', AuxT(i,1:Nvib)
        enddo
    endif


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
!           Aux3(1,1:Nvib) = Aux3(1,1:Nvib)/dsqrt(Theta)
          call sort_vec_max(Aux3(1,1:Nvib),ipiv(1:Nvib),Nvib)

          !Determine symmetry
          do j=1,Nvib
              jj=ipiv(j)
              if (S_sym(jj) == jj) cycle
              ii = S_sym(jj)
              Theta = Aux(i,jj)/Aux(i,ii)
              if (Theta > 0.d0) then 
                  ModeSymm="A' "
              else
                  ModeSymm="A''"
              endif
              exit
          enddo

          print*, ""
          print'(A,I4,6X,A,F8.3,4X,A,A,F8.3)', "Mode ", i, " Freq(cm^-1) = ", Freq(i), "Symm ", ModeSymm, Theta
          print*, "      S        Coef.     Contrib.(%)          Description"
          print*, " ======================================================================="
          Theta = 0.d0
          kk=0
          do j=1,Nvib
              if (Theta > 0.9d0) exit
              jj = ipiv(j)
              Theta = Theta + Aux3(1,j)!**2
              print'(5x,i4,4x,f9.3,4x,f8.3,4x,a)', jj, Aux(i,jj), Aux3(1,j)*100, trim(adjustl(ModeDef(jj)))
              kk=kk+1
          enddo 
          print*, " ========================================================================"
          write(6,'(A,I3)') "     Total Number of internal to describe >90% of the mode: ", kk
      enddo


      return

end subroutine gf_method


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


subroutine axis_swithching(molec,molec2,T0)

    ! Subroutine to account for the axis swithching effect. molec2 is rotated to
    ! fulfil the Eckart frame simulataneously with molec

    use structure_types
    use MatrixMod


    type(str_resmol),intent(in)    :: molec
    type(str_resmol),intent(inout) :: molec2 
    real(8),dimension(3,3)          :: T0

    !Local
    type(str_resmol)   :: molec_aux
    real(8),dimension(1:3,1:3)     :: C,Cinv
    real(8) :: xaux,yaux,zaux, xrot,yrot,zrot, rsum, rcheck=99999.d0
    integer :: Nat
    !Other auxiliar
    real(8),dimension(3,3) :: Aux,AuxT,T
    real(8),dimension(3)   :: Vec,Vec2

    !Stuff to use lapack inversion (should be in MatrixMod)
    integer,parameter :: NDIM=6
    real(8),dimension(NDIM,NDIM) :: work
    integer,dimension(NDIM,NDIM) :: ipiv
    integer :: info

    !Load shortcuts
    Nat = molec%natoms

    !Following the expresion given by Borrelli 2006

!=========================
    C(1,1) = 0.d0
    C(2,2) = 0.d0
    C(3,3) = 0.d0
    C(1,2) = 0.d0
    C(2,1) = 0.d0
    C(1,3) = 0.d0
    C(3,1) = 0.d0
    C(2,3) = 0.d0
    C(3,2) = 0.d0
    do i=1,Nat
        C(1,1) = C(1,1) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%x
        C(2,2) = C(2,2) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%y
        C(3,3) = C(3,3) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%z

        C(1,2) = C(1,2) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%y
        C(1,3) = C(1,3) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%z
        C(2,3) = C(2,3) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%z

        C(2,1) = C(2,1) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%x
        C(3,1) = C(3,1) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%x
        C(3,2) = C(3,2) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%y
    enddo
!     print*, "C_ini"
!     do i=1,3
!         print'(100(F8.3,2X))', C(i,1:3)
!     enddo
    !Inverse of C
!     nn=3
    Cinv=C
    call dgetrf(3,3, Cinv, 3, ipiv, info)
    call dgetri(3,   Cinv, 3, ipiv, work, NDIM, info)


!     Aux=0.d0
!     Aux(1:3,1:3) = matmul(C(1:3,1:3),Cinv(1:3,1:3))
!     print*, "CinvC", info
!     do i=1,3
!         print'(100(F8.3,2X))', Aux(i,1:3)
!     enddo


    do i=1,3
        do j=1,3
            Aux(i,j) = 0.d0
            do  k=1,3
                Aux(i,j) = Aux(i,j) + C(k,i)*C(k,j)
            enddo
       enddo
    enddo
    call diagonalize_full(Aux(1:3,1:3),3,AuxT(1:3,1:3),Vec(1:3),"lapack")


!     do i=1,8
! 
!     Vec2=1.d0
! 
!         iT = 2*i/3
!         jT = 2*(i+1)/3
!         kT = 2*(i+2)/3
! 
!         Vec2(1) = (-1.d0)**iT
!         Vec2(2) = (-1.d0)**jT
!         Vec2(3) = (-1.d0)**kT
! 
!         print'(3F8.3,X,3I3)', Vec2(1:3), iT,jT,kT
! 
!     enddo
! stop

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
    

    !(C^t C)^1/2
    do i=1,3
        do j=1,3
            C(i,j) = 0.d0
            do k=1,3
                if (dabs(Vec(k)) < 1d-10) Vec(k) = 0.d0
                C(i,j) = C(i,j) + AuxT(i,k)*dsqrt(Vec(k))*AuxT(j,k)*Vec2(k)
            enddo
        enddo
    enddo

    !T = C^t C)^1/2 C^-1
    do i=1,3
        do j=1,3
            T(i,j) = 0.d0
            do k=1,3
                T(i,j) = T(i,j) + C(i,k)*Cinv(k,j)
            enddo
        enddo
    enddo
    !Manage 0-eigenvalues in diagonal matrix
    Vec2(1:3)=0.d0
    do i=1,3
        do j=1,3
           Vec2(i) = Vec2(i) + T(i,j)**2
        enddo
    enddo
    !To do: add support for one-off diagonal eigenvalue in C
    ! If one eigenvalue is 0
    if (abs(Vec2(1)) < 1.d-10) then
        print*, "Warning: eingen 1 value is zero"
        T(1,1) = T(2,2)*T(3,3)-T(2,3)*T(3,2)       
        T(1,2) = T(2,3)*T(3,1)-T(2,1)*T(3,3)       
        T(1,3) = T(2,1)*T(3,2)-T(2,2)*T(3,1)       
    else if (abs(Vec2(2)) < 1.d-10) then
        print*, "Warning: eingen 2 value is zero"
        T(2,1) = T(1,2)*T(3,3)-T(1,3)*T(3,2)        
        T(2,2) = T(1,3)*T(3,1)-T(1,1)*T(3,3)        
        T(2,3) = T(1,1)*T(3,2)-T(1,2)*T(3,1)       
    else if (abs(Vec2(3)) < 1.d-10) then
        print*, "Warning: eingen 3 value is zero"
        T(3,1) = T(1,2)*T(2,3)-T(1,3)*T(2,2)        
        T(3,2) = T(1,3)*T(2,1)-T(1,1)*T(2,3)        
        T(3,3) = T(1,1)*T(2,2)-T(1,2)*T(2,1)   
    endif 
        
    print*, "T",iT
    do i=1,3
        print'(100(F8.3,2X))', T(i,1:3)
    enddo

    det =       T(1,1)*T(2,2)*T(3,3)
    det = det + T(2,1)*T(3,2)*T(1,3)
    det = det + T(1,2)*T(2,3)*T(3,1)
    det = det - T(3,1)*T(2,2)*T(1,3)
    det = det - T(2,1)*T(1,2)*T(3,3)
    det = det - T(3,2)*T(2,3)*T(1,1)

    print*, det

    if (det<0.d0) cycle

    !Rotate and check sum of distances
    molec_aux=molec2
    rsum=0.d0
    do i=1,Nat
        xaux = molec2%atom(i)%x
        yaux = molec2%atom(i)%y
        zaux = molec2%atom(i)%z
        xrot = T(1,1)*xaux + T(2,1)*yaux + T(3,1)*zaux - molec%atom(i)%x
        yrot = T(1,2)*xaux + T(2,2)*yaux + T(3,2)*zaux - molec%atom(i)%y
        zrot = T(1,3)*xaux + T(2,3)*yaux + T(3,3)*zaux - molec%atom(i)%z
        rsum=rsum + dsqrt(xrot**2+yrot**2+zrot**2)
    enddo
    print*, "Sum of distance ", rsum

    if (rsum < rcheck) then
        T0=T
        rcheck=rsum
    endif

    enddo ! Possible T values

    !Rotate State 2 
    !---------------
    print*, "T0"
    do i=1,3
        print'(100(F8.3,2X))', T0(i,1:3)
    enddo
    det =       T0(1,1)*T0(2,2)*T0(3,3)
    det = det + T0(2,1)*T0(3,2)*T0(1,3)
    det = det + T0(1,2)*T0(2,3)*T0(3,1)
    det = det - T0(3,1)*T0(2,2)*T0(1,3)
    det = det - T0(2,1)*T0(1,2)*T0(3,3)
    det = det - T0(3,2)*T0(2,3)*T0(1,1)
    print*, "Det", det
    rsum=0.d0
    do i=1,Nat
        xaux = molec2%atom(i)%x
        yaux = molec2%atom(i)%y
        zaux = molec2%atom(i)%z
        molec2%atom(i)%x = T0(1,1)*xaux + T0(2,1)*yaux + T0(3,1)*zaux
        molec2%atom(i)%y = T0(1,2)*xaux + T0(2,2)*yaux + T0(3,2)*zaux
        molec2%atom(i)%z = T0(1,3)*xaux + T0(2,3)*yaux + T0(3,3)*zaux
        rsum=rsum + dsqrt((molec2%atom(i)%x-molec%atom(i)%x)**2+&
                          (molec2%atom(i)%y-molec%atom(i)%y)**2+&
                          (molec2%atom(i)%z-molec%atom(i)%z)**2)
    enddo
    print*, "Sum of distance ", rsum
    !Check new C matrix
    !C(1,1)
    C(1,1) = 0.d0
    C(2,2) = 0.d0
    C(3,3) = 0.d0
    C(1,2) = 0.d0
    C(2,1) = 0.d0
    C(1,3) = 0.d0
    C(3,1) = 0.d0
    C(2,3) = 0.d0
    C(3,2) = 0.d0
    do i=1,Nat
        C(1,1) = C(1,1) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%x
        C(2,2) = C(2,2) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%y
        C(3,3) = C(3,3) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%z

        C(1,2) = C(1,2) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%y
        C(1,3) = C(1,3) + molec%atom(i)%mass*molec%atom(i)%x*molec2%atom(i)%z
        C(2,3) = C(2,3) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%z

        C(2,1) = C(2,1) + molec%atom(i)%mass*molec%atom(i)%y*molec2%atom(i)%x
        C(3,1) = C(3,1) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%x
        C(3,2) = C(3,2) + molec%atom(i)%mass*molec%atom(i)%z*molec2%atom(i)%y
    enddo
    print*, "C_new"
    do i=1,3
        print'(100(F8.3,2X))', C(i,1:3)
    enddo
!=========================

    return
end subroutine axis_swithching

!=========================================================================
! The following are somehow redundant and should be merged in the 
!  corresponding "leading" subroutine
!=========================================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine internal_fc(Hess,molec,S_sym,ModeDef,B,G,verbose)

    ! Compute bond and angle force constants from the Cartesian Hessian using the B matrix
    ! Internal coordinates are those used to construct the B matrix (dihedrals are not considered)
    ! Based on gf_method SR (till comput. of Fint)

    use structure_types
    use line_preprocess
    use alerts
    use constants
    use atomic_geom
    use MatrixMod

    implicit none

    integer,parameter :: NDIM = 800

    !====================== 
    !ARGUMENTS
    !Inpu
    type(str_resmol),intent(in) :: molec                     !Input molecule
    real(8),dimension(1:NDIM,1:NDIM),intent(in) ::    G,B    !Wilson matrices
    character(len=100),dimension(NDIM),intent(in) :: ModeDef !Definition of modes
    integer,dimension(NDIM),intent(in) :: S_sym              !Manages symmetric internals
    logical,intent(in) :: verbose                            !Verbose level
    !Input-Output
    real(8),dimension(1:NDIM,1:NDIM),intent(inout) :: Hess   !Hessian: cart(in)-intern(out)
    !====================== 

    !====================== 
    !System variables
    integer :: Nat, Nvib, Nbonds, Nangles
    character(len=3) :: ModeSymm
    !====================== 

    !====================== 
    !Internal analysis 
    real(8),dimension(1:NDIM,1:NDIM) :: Ginv, X, Xinv, L, HessP
    !Auxiliar matrices
    real(8),dimension(NDIM,NDIM) :: AuxT,Aux, Aux3
    real(8),dimension(1:NDIM) :: Freq
    !====================== 

    !====================== 
    !Auxiliars for LAPACK matrix nversion
    real(8),dimension(NDIM,NDIM) :: work
    integer,dimension(NDIM) :: ipiv
    integer :: info
    !====================== 

    !====================== 
    !Auxiliar variables
    integer :: ierr
    character(1) :: null
    character(len=16) :: dummy_char
    real(8) :: Theta, mu
    !TESTING ARRAYS (auxiliar)
    real(8),dimension(NDIM,NDIM) :: test1,test2,test3
    !====================== 

    !=============
    !Counters
    integer :: i,j,k, ii,jj,kk, iat
    !=============

! (End of variables declaration) 
!==================================================================================

    Nat = molec%natoms
    Nvib = 3*Nat-6
    Nangles = S_sym(3*Nat+2) ! 
    Nbonds = S_sym(3*Nat+1)

    !=====================================================================
    ! HESSIAN IN INTERNAL COORDINATES (JCC, 17, 49-56, by Frisch et al)
    !======================================================================
    !From JCC, 17, 49-56,:
     ! Hint = G^- Bu(Hx+B'^Tg_q)u^TB^T G^-
     !g_q is the gradient, so g_q=0
     !G^- is the generalized inverse (for redundant internal) or simply the
     !inverse for nonredundant

    !Inverse of G
    Ginv(1:Nvib,1:Nvib)=G(1:Nvib,1:Nvib)
    call dsytrf('L', Nvib, Ginv, NDIM, ipiv, work, NDIM, info)
    call dsytri('L', Nvib, Ginv, NDIM, ipiv, work, info)
!   RECONSTRUCT UPPER PART
    do i=1,Nvib
    do j=i+1,Nvib
        Ginv(i,j) = Ginv(j,i)
    enddo
    enddo

    !Compute G^-1Bu  (where u is the inverse mass matrix)
    do i=1,Nvib
        j=0
        do jj=1,Nat
        do iat=1,3
            j=j+1
            Aux(i,j) = 0.d0
            mu=1.d0/molec%atom(jj)%mass/UMAtoAU
            do k=1,Nvib
                Aux(i,j) = Aux(i,j) + Ginv(i,k)*B(k,j)                
            enddo
            Aux(i,j) = Aux(i,j)*mu
        enddo
        enddo
    enddo

    !Now: Hint = Aux Hcart Aux^T
    do i=1,Nvib
        do j=1,3*Nat
            work(i,j) = 0.d0
            do k=1,3*Nat
                work(i,j) = work(i,j) + Aux(i,k)*Hess(k,j)                
            enddo
        enddo
    enddo
    do i=1,Nvib
        do j=1,3*Nat
            Hess(i,j) = 0.d0
            do k=1,3*Nat
                Hess(i,j) = Hess(i,j) + work(i,k)*Aux(j,k)                
            enddo
        enddo
    enddo

    if (verbose) then
    print*, ""
    print*, "F MATRIX"
    do i=1,Nvib
        print'(100(F10.3,2X))', Hess(i,1:Nvib)
    enddo 
    endif
    test1(1:Nvib,1:Nvib) = Hess(1:Nvib,1:Nvib)


    !FG PROCEDURE (Wilson, Decius and Cross, Section 4.3)
    ! |GF - \lamda| = 0  <==> GFL = L\lamda   means the diagonalize of GF, but is not symmetric
    ! Solutions
    ! *Orthogonalization of FL = G^-1L\lamda similar to orthogonalization fo S in Roothan. 
    ! *Symmetrize GF and diagonalize following Mizayawa (Sect. 4.6, Califano, Vibrational States)

    ! FOLLOWING FIRST ALTERNATIVE (similar to diagonalization of basis set)
    call diagonalize_full(G(1:Nvib,1:Nvib),Nvib,Aux(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
    !The internal coordinates also need to be rotated. But using G^-1/2? Para que??
    ! No, they should not be changed since the ones finally used are the "original" ones, not this orthogonalized set 

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

    if (verbose) then
    print*, ""
    print*, "X="
    do i=1,Nvib
        print'(100(F10.3,2X))', X(i,1:Nvib)
    enddo
    endif

    !Now rotate F, F'=X^TFX. Store the rotated matrix in F
    do i=1,Nvib
        do j=1,Nvib
            Aux(i,j) = 0.d0
            do k=1,Nvib
                Aux(i,j) = Aux(i,j) + X(k,i)*Hess(k,j)
            enddo
        enddo
    enddo
    do i=1,Nvib
        do j=1,Nvib
            HessP(i,j) = 0.d0
            do k=1,Nvib
                HessP(i,j) = HessP(i,j) + Aux(i,k)*X(k,j)
            enddo
        enddo
    enddo

    !We can now diagonalize F
    call diagonalize_full(HessP(1:Nvib,1:Nvib),Nvib,L(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
    !WE NEED TO PASS "L" AS IT IS NOW

    !Check freqcuencies
    print*, ""
    print*, "NORMAL FORCE CONSTANTS (A.U.)"
!     call sort_vec(Freq,Nvib)
    do i=1,Nvib
          write(6,*) Freq(i)
    enddo

    !Check freqcuencies
    print*, ""
    print*, "FREQUENCIES (cm^-1)"
!     call sort_vec(Freq,Nvib)
    do i=1,Nvib
          Freq(i) = sign(dsqrt(abs(Freq(i))*HARTtoJ/BOHRtoM**2/AUtoKG)/2.d0/pi/clight/1.d2,&
                         Freq(i))
          write(6,*) Freq(i)
    enddo
    
!     print*, ""
!     print*, "FORCE CONSTANTS (NO DIHEDRALS):"
!     do i=1,Nbonds
!         print'(I4,X,A10,F10.2,A,F10.2)', i, adjustl(ModeDef(i)), Hess(i,i)*938211.3414d0, " (KJ/mol/nm^2)  Scaled: ", Hess(i,i)*938211.3414d0*0.89
!     enddo
!     do i=Nbonds+1,Nangles
!         print'(I4,X,A10,F10.2,A,F10.2)', i, adjustl(ModeDef(i)), Hess(i,i)*2625.5d0,      " (KJ/mol/rad^2) Scaled: ", Hess(i,i)*2625.5d0*0.89
!     enddo

!     stop
    return

end subroutine internal_fc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine internal_Wilson_fc(molec,S,S_sym,ModeDef,B,G,Asel,verbose)

    !Adapted from internal_Wilson. This version is addapted for its used with internal_fc

    use structure_types
    use line_preprocess
    use alerts
    use constants
    use atomic_geom
    use MatrixMod

    implicit none

    integer,parameter :: NDIM = 800
    real(8),parameter :: ZEROp = 1.d-10 !practically zero

    !====================== 
    !ARGUMENTS
    !Input
    type(str_resmol),intent(in) :: molec            ! Input molecule to which the vib analysis is performed
    logical,intent(in) :: verbose                   ! Control Output details
    integer,dimension(NDIM),intent(in) :: S_sym     ! Manages symmetric internals
    !Input-Output
    real(8),dimension(NDIM,NDIM),intent(inout) :: Asel !Used to manage nonredundant internal variables (disabled option)
    !Output
    real(8),dimension(NDIM,NDIM),intent(out) :: B, G          !B -Wilson- and G matrices
    real(8),dimension(NDIM),intent(out) :: S                  !Vector of internal coordinates
    character(len=100),dimension(NDIM),intent(out) :: ModeDef !Definition of normal modes
    !====================== 

    !====================== 
    !Auxiliar variables
    character(1) :: null
    character(len=100) :: dummy_char
    !====================== 

    !====================== 
    !System info
    integer,dimension(1:NDIM,1:4) :: bond_s, angle_s, dihed_s
    integer :: nbonds, ndihed, nangles
    integer :: Nat, Nvib, Nred
    character(len=5) :: PG
    character(len=1) :: connector1, connector2
    !====================== 

    !====================== 
    !INTERNAL ANALYSIS
    !B and G matrices
!     real(8),dimension(NDIM,NDIM),intent(out) :: B, G
    !AUXILIAR MATRICES
    real(8),dimension(NDIM,NDIM) :: AuxT,Aux, Aux3
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
    real(8) :: mu
    real(8),dimension(NDIM) :: Vec, Vec2
    !====================== 

    !====================== 
    !Auxiliars for LAPACK matrix nversion
    integer :: info
    integer,dimension(NDIM) :: ipiv
    real(8),dimension(NDIM,NDIM) :: work
    !====================== 

    !=============
    !Counters
    integer :: i,j,k, ii,jj,kk,kkk, iat
    integer :: i_1, i_2, i_3, i_4, imax, jmax
    !=============

! (End of variables declaration) 
!==================================================================================


    !Set bonded
    nbonds  = molec%geom%nbonds
    nangles = molec%geom%nangles
    ndihed  = molec%geom%ndihed
    bond_s(1:nbonds,1:2)  =  molec%geom%bond(1:nbonds,1:2)
    angle_s(1:nangles,1:3) =  molec%geom%angle(1:nangles,1:3)
    dihed_s(1:ndihed,1:4)  =  molec%geom%dihed(1:ndihed,1:4)


    !Initialize matrices
    Nat = molec%natoms
    Nvib = 3*Nat-6
    B(1:Nvib,1:3*Nat) = 0.d0
    G(1:Nvib,1:Nvib)  = 0.d0

!=============================================================
!   Ref: Decius, Cross and Wilson (Section 4.2) --
!    The same nomenclature is used (index refer to the same
!    atoms, even if reversed)
!==============================================================

    !k-index runs over internal coordinates
    k=0
    write(6,'(/,A,I3)') "BONDS", nbonds !molec%natoms-1
!     do i=2,molec%natoms
    do i=1,nbonds
        k=k+1

        i_1 = bond_s(i,1)
        i_2 = bond_s(i,2)
        r21 = calc_dist(molec%atom(i_1),molec%atom(i_2))
        !Set bond order
        if (r21*BOHRtoAMS < 1.42d0                .and. &
            adjustl(molec%atom(i_1)%name) == "C"   .and. &
            adjustl(molec%atom(i_2)%name) == "C") then
            connector1 = "="
        else if (r21*BOHRtoAMS < 1.5d0                .and. &
            adjustl(molec%atom(i_1)%name) == "C"   .and. &
            adjustl(molec%atom(i_2)%name) == "C") then
            connector1 = "#"
        else
            connector1 = "-"
        endif
           
        S(k) = r21
        write(6,'(I4,X,A,2(I3,A),2F15.8)') k,trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                            trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                            r21*BOHRtoAMS, r21
        write(ModeDef(k),'(A)') trim(adjustl(molec%atom(i_1)%name))//connector1//&
                            trim(adjustl(molec%atom(i_2)%name))

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

    write(6,'(/,A,I3)') "ANGLES", nangles! molec%natoms-2
!     do i=3,molec%natoms
    do i=1,nangles
        k=k+1

        i_1 = angle_s(i,1)
        i_3 = angle_s(i,2)
        i_2 = angle_s(i,3)

        !Set bond order
        r31=calc_dist(molec%atom(i_1),molec%atom(i_3))
        if (r31*BOHRtoAMS < 1.42d0                .and.&
            adjustl(molec%atom(i_1)%name) == "C"   .and. &
            adjustl(molec%atom(i_2)%name) == "C") then
            connector1 = "="
        else if (r31*BOHRtoAMS < 1.5d0                .and. &
            adjustl(molec%atom(i_1)%name) == "C"   .and. &
            adjustl(molec%atom(i_2)%name) == "C") then
            connector1 = "#"
        else
            connector1 = "-"
        endif
        r32=calc_dist(molec%atom(i_2),molec%atom(i_3))
        if (r32*BOHRtoAMS < 1.42d0                .and.&
            adjustl(molec%atom(i_1)%name) == "C"   .and. &
            adjustl(molec%atom(i_2)%name) == "C") then
            connector2 = "="
        else if (r32*BOHRtoAMS < 1.5d0                .and. &
            adjustl(molec%atom(i_1)%name) == "C"   .and. &
            adjustl(molec%atom(i_2)%name) == "C") then
            connector2 = "#"
        else
            connector2 = "-"
        endif

        ang1 = calc_angle(molec%atom(i_1),molec%atom(i_3),molec%atom(i_2))
        S(k) = ang1
        write(6,'(I4,X,A,3(I3,A),F15.8)') k,trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                            trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                            trim(adjustl(molec%atom(i_2)%name))//"(",i_2,")", &
                            ang1*360.d0/2.d0/pi
        write(ModeDef(k),'(A)') trim(adjustl(molec%atom(i_1)%name))//connector1//&
                            trim(adjustl(molec%atom(i_3)%name))//connector2//&
                            trim(adjustl(molec%atom(i_2)%name))

        !Three cart displacements different from zero.
!          r31=calc_dist(molec%atom(i_1),molec%atom(i_3))
         e31x = (molec%atom(i_1)%x-molec%atom(i_3)%x)/r31
         e31y = (molec%atom(i_1)%y-molec%atom(i_3)%y)/r31
         e31z = (molec%atom(i_1)%z-molec%atom(i_3)%z)/r31
!          r32=calc_dist(molec%atom(i_2),molec%atom(i_3))
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

!     DIHEDRALS ARE NOT CONSIDERED (BUT COMPUTED!!!)
    write(6,'(/,A,I3)') "DIHEDRALS", ndihed!molec%natoms-3
!     do i=4,molec%natoms
    do i=1,ndihed
        k=k+1

        i_1 = dihed_s(i,1)
        i_2 = dihed_s(i,2)
        i_3 = dihed_s(i,3)
        i_4 = dihed_s(i,4)
        ang1 = calc_dihed_new(molec%atom(i_1),molec%atom(i_2),molec%atom(i_3),molec%atom(i_4))
        S(k) = ang1
        write(6,'(A,4(I3,A),F15.8)') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                            trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                            trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                            trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")", &
                            ang1*360.d0/2.d0/pi
        write(ModeDef(k),'(A,4(I3,A))') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
                            trim(adjustl(molec%atom(i_2)%name))//"(",i_2,") -- "//&
                            trim(adjustl(molec%atom(i_3)%name))//"(",i_3,") -- "//&
                            trim(adjustl(molec%atom(i_4)%name))//"(",i_4,")"

        !Four cart displacements different from zero (some index intercheged with Decius..)
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
        B(k,3*i_1-2) =  e21Pe32x/(r21*dsin(ang2)**2)
        B(k,3*i_1-1) =  e21Pe32y/(r21*dsin(ang2)**2)
        B(k,3*i_1  ) =  e21Pe32z/(r21*dsin(ang2)**2)
        !s2 
        B(k,3*i_2-2) = -((r32-r21*dcos(ang2))*e21Pe32x/(r32*r21*dsin(ang2)**2) &
                     +  dcos(ang3)*e43Pe32x/(r32*dsin(ang3)**2))
        B(k,3*i_2-1) = -((r32-r21*dcos(ang2))*e21Pe32y/(r32*r21*dsin(ang2)**2) &
                     +  dcos(ang3)*e43Pe32y/(r32*dsin(ang3)**2))
        B(k,3*i_2  ) = -((r32-r21*dcos(ang2))*e21Pe32z/(r32*r21*dsin(ang2)**2) &
                     +  dcos(ang3)*e43Pe32z/(r32*dsin(ang3)**2))
        !Following Decius paper (instead of W-D-C book formulation)
!         B(k,3*i_2-2) = ( 1.d0/(r21*dsin(ang2)) - (dcos(ang1)/dtan(ang3) + &
!                          1.d0/dtan(ang2))/r32 )*e21Pe32x/dsin(ang2) +     &
!                          dsin(ang1)/dtan(ang3)/dsin(ang2)/r32*e32Pe21Pe32x
!         B(k,3*i_2-1) = ( 1.d0/(r21*dsin(ang2)) - (dcos(ang1)/dtan(ang3) + &
!                          1.d0/dtan(ang2))/r32 )*e21Pe32y/dsin(ang2) +     &
!                          dsin(ang1)/dtan(ang3)/dsin(ang2)/r32*e32Pe21Pe32y
!         B(k,3*i_2  ) = ( 1.d0/(r21*dsin(ang2)) - (dcos(ang1)/dtan(ang3) + &
!                          1.d0/dtan(ang2))/r32 )*e21Pe32z/dsin(ang2) +     &
!                          dsin(ang1)/dtan(ang3)/dsin(ang2)/r32*e32Pe21Pe32z
        !s3
        B(k,3*i_3-2) = -((r32-r43*dcos(ang3))*e43Pe32x/(r32*r43*dsin(ang3)**2) &
                     +  dcos(ang2)*e21Pe32x/(r32*dsin(ang2)**2))
        B(k,3*i_3-1) = -((r32-r43*dcos(ang3))*e43Pe32y/(r32*r43*dsin(ang3)**2) &
                     +  dcos(ang2)*e21Pe32y/(r32*dsin(ang2)**2))
        B(k,3*i_3  ) = -((r32-r43*dcos(ang3))*e43Pe32z/(r32*r43*dsin(ang3)**2) &
                     +  dcos(ang2)*e21Pe32z/(r32*dsin(ang2)**2))
        !Following Decius paper (instead of W-D-C book formulation)
!         B(k,3*i_3-2) = ( 1.d0/(r43*dsin(ang3)) - (dcos(ang1)/dtan(ang2) + &
!                          1.d0/dtan(ang3))/r32 )*e43Pe32x/dsin(ang3) -     &
!                          dsin(ang1)/dtan(ang2)/dsin(ang3)/r32*e32Pe43Pe32x
!         B(k,3*i_3-1) = ( 1.d0/(r43*dsin(ang3)) - (dcos(ang1)/dtan(ang2) + &
!                          1.d0/dtan(ang3))/r32 )*e43Pe32y/dsin(ang3) -     &
!                          dsin(ang1)/dtan(ang2)/dsin(ang3)/r32*e32Pe43Pe32y
!         B(k,3*i_3  ) = ( 1.d0/(r43*dsin(ang3)) - (dcos(ang1)/dtan(ang2) + &
!                          1.d0/dtan(ang3))/r32 )*e43Pe32z/dsin(ang3) -     &
!                          dsin(ang1)/dtan(ang2)/dsin(ang3)/r32*e32Pe43Pe32z
        !s4
        B(k,3*i_4-2) =  e43Pe32x/(r43*dsin(ang3)**2)
        B(k,3*i_4-1) =  e43Pe32y/(r43*dsin(ang3)**2)
        B(k,3*i_4  ) =  e43Pe32z/(r43*dsin(ang3)**2)

    enddo
    Nred = k

    if (verbose) then
    write(6,*) ""
    write(6,*) "B MATRIX-"
    do i=1,Nred
        write(6,'(100(F9.5,X))') B(i,1:3*Nat)
    enddo
    endif

    !Symmetry adapted coordinates
!     Aux(1,1:3*Nat) = (B(1,1:3*Nat)+B(6,1:3*Nat))
!     Aux(2,1:3*Nat) = (B(1,1:3*Nat)-B(6,1:3*Nat))
!     B(1,1:3*Nat) = Aux(1,1:3*Nat)
!     B(6,1:3*Nat) = Aux(2,1:3*Nat)
!     Aux(1,1:3*Nat) = (B(2,1:3*Nat)+B(3,1:3*Nat))
!     Aux(2,1:3*Nat) = (B(2,1:3*Nat)-B(3,1:3*Nat))
!     B(2,1:3*Nat) = Aux(1,1:3*Nat)
!     B(3,1:3*Nat) = Aux(2,1:3*Nat)
!     Aux(1,1:3*Nat) = (B(4,1:3*Nat)+B(5,1:3*Nat))
!     Aux(2,1:3*Nat) = (B(4,1:3*Nat)-B(5,1:3*Nat))
!     B(4,1:3*Nat) = Aux(1,1:3*Nat)
!     B(5,1:3*Nat) = Aux(2,1:3*Nat)
    if (S_sym(3*Nat) /= 0) then
         print*, ""
         print*, "Symmetry addapted coordinates will be used"
         do i=1,Nvib
             if (S_sym(i) <= i) cycle
             j=S_sym(i)
             Aux(1,1:3*Nat) = (B(i,1:3*Nat)+B(j,1:3*Nat))
             Aux(2,1:3*Nat) = (B(i,1:3*Nat)-B(j,1:3*Nat))
             B(i,1:3*Nat) = Aux(1,1:3*Nat)
             B(j,1:3*Nat) = Aux(2,1:3*Nat)
             dummy_char = trim(adjustl(ModeDef(i)))//"+"//trim(adjustl(ModeDef(j)))
             ModeDef(j) = trim(adjustl(ModeDef(i)))//"-"//trim(adjustl(ModeDef(j)))
             ModeDef(i) = trim(adjustl(dummy_char))
         enddo
    endif


    !CONCTRUCT G
    do i=1,Nred
        do j=1,Nred
            G(i,j) = 0.d0
            k=0
            do kk=1,Nat
            do iat=1,3
                k=k+1
                mu = 1.d0/molec%atom(kk)%mass/UMAtoAU
                G(i,j) = G(i,j) + mu*B(i,k)*B(j,k)
            enddo 
            enddo
        enddo
    enddo

    ! Check symmetry
    mu = 0.d0
    do i=1,Nred
        do j=i,Nred
            Aux(i,j) = (G(i,j)-G(j,i))/G(i,j) *100.d0
            mu = max(mu,dabs(Aux(i,j)))
            if (mu == Aux(i,j)) then
                imax = i
                jmax = j
            endif
        enddo
    enddo
    write(6,*) ""
    print*, "Maximum difference in symmetric matrix (%):", mu, "for", imax,jmax
    print*, G(imax,jmax), G(jmax,imax)
    write(6,*) ""

    if (verbose) then
    write(6,*) ""
    write(6,*) "G MATRIX"
    do i=1,Nred
        write(6,'(100(E10.2,3X))') G(i,1:Nvib)
    enddo
    print*, ""
    print*, "G MATRIX (lower triangular)"
    do i=1,Nvib
        do j=1,i
            print'(G15.8)', G(i,j)
        enddo
    enddo
    print*, "--END OF G MATRIX (lower triangular)"
    endif

    return

end subroutine internal_Wilson_fc



end module internal_module
