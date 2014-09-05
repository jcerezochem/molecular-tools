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
 ! Version 9: Include Bders (this branch has B elements corrected)
 !**********************
 ! SR9_v4: addpated to molecular tools


 contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

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

    integer,parameter :: NDIM = 400
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
        write(6,'(A,2(I3,A),2F15.8)') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
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
        write(6,'(A,3(I3,A),F15.8)') trim(adjustl(molec%atom(i_1)%name))//"(",i_1,") -- "//&
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
        B(k,3*i_3-2) = ((r32-r43*dcos(ang3))*e43Pe32x/(r32*r43*dsin(ang3)**2) &
                     +  dcos(ang2)*e21Pe32x/(r32*dsin(ang2)**2))
        B(k,3*i_3-1) = ((r32-r43*dcos(ang3))*e43Pe32y/(r32*r43*dsin(ang3)**2) &
                     +  dcos(ang2)*e21Pe32y/(r32*dsin(ang2)**2))
        B(k,3*i_3  ) = ((r32-r43*dcos(ang3))*e43Pe32z/(r32*r43*dsin(ang3)**2) &
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
        B(k,3*i_4-2) =  -e43Pe32x/(r43*dsin(ang3)**2)
        B(k,3*i_4-1) =  -e43Pe32y/(r43*dsin(ang3)**2)
        B(k,3*i_4  ) =  -e43Pe32z/(r43*dsin(ang3)**2)

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
    endif

    print*, "Redundant coordinates", Nred-Nvib

    !IF WE USE ALL BONDED PARAMETERS,WE HAVE REDUNDANCY. WE SELECT A NON-REDUNDANT
    !COMBINATION FROM THE NON-ZERO EIGENVALUES OF G (Reimers 2001, JCP)
    !--- this is not tested for new versions and must not be used ---
    if (Nred-Nvib /= 0) then    
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

subroutine internal_Wilson_fc(molec,S,S_sym,ModeDef,B,G,Asel,verbose)

    !Adapted from internal_Wilson. This version is addapted for its used with internal_fc

    use structure_types
    use line_preprocess
    use alerts
    use constants
    use atomic_geom
    use MatrixMod

    implicit none

    integer,parameter :: NDIM = 400
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
    endif

    return

end subroutine internal_Wilson_fc

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

    integer,parameter :: NDIM = 400

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

    integer,parameter :: NDIM = 400

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
    Nvib = 3*Nat-6

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

    !Now rotate F, F'=X^TFX. Store the rotated matrix in Aux3 (temporary array)
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
            Aux3(i,j) = 0.d0
            do k=1,Nvib
                Aux3(i,j) = Aux3(i,j) + Aux(i,k)*X(k,j)
            enddo
        enddo
    enddo

    !We can now diagonalize F
    call diagonalize_full(Aux3(1:Nvib,1:Nvib),Nvib,L(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine gf_method_V(Hess,Grad,molec,S_sym,ModeDef,L,B,Bder,G,Freq,Asel,X,Xinv,verbose)

    !Performs vibrational analysis using ghe Wilson B matrix and the GF method.
    !Returns L matrix (in Aux) and B matrix (in B)
    ! This version includes the use of Vertical models!

    use structure_types
    use line_preprocess
    use alerts
    use constants
    use atomic_geom
    use MatrixMod

    implicit none

    integer,parameter :: NDIM = 400

    !====================== 
    !ARGUMENTS
    !Inpu
    type(str_resmol),intent(in) :: molec                     !Input molecule
    real(8),dimension(1:NDIM,1:NDIM),intent(in) ::    G,B    !Wilson matrices
    real(8),dimension(1:NDIM,1:NDIM,1:NDIM),intent(in) ::    Bder
    character(len=100),dimension(NDIM),intent(in) :: ModeDef !Definition of modes
    integer,dimension(NDIM),intent(in) :: S_sym              !Manages symmetric internals
    logical,intent(in) :: verbose                            !Verbose level
    !Input-Output
    real(8),dimension(1:NDIM,1:NDIM),intent(inout) :: Hess   !Hessian: cart(in)-intern(out)
    real(8),dimension(1:NDIM),intent(inout) :: Grad
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
    Nvib = 3*Nat-6

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

    ! Get the gradient in internal coords first: gq = G^-1Bu(gx)
    ! use Freq as axiliar vector
    do i=1,Nvib
        Freq(i) = 0.d0
        do j=1,3*Nat
            Freq(i) = Freq(i) + Aux(i,j) * Grad(j)
        enddo
    enddo
    ! .. and multiply: Bder(i,j,K)^t * gq(K)
    do i=1,3*Nat
    do j=1,3*Nat
        AuxT(i,j) = 0.d0
        do k=1,Nvib
            AuxT(i,j) = AuxT(i,j) + Bder(k,i,j)*Freq(k)
        enddo
    enddo
    enddo


    !Now (modified for Vertical): Hint = Aux (Hcart-Bder*gq) Aux^T
    do i=1,Nvib
        do j=1,3*Nat
            work(i,j) = 0.d0
            do k=1,3*Nat
                work(i,j) = work(i,j) + Aux(i,k)*(Hess(k,j)-AuxT(k,j))                
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

    !Now rotate F, F'=X^TFX. Store the rotated matrix in Aux3 (temporary array)
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
            Aux3(i,j) = 0.d0
            do k=1,Nvib
                Aux3(i,j) = Aux3(i,j) + Aux(i,k)*X(k,j)
            enddo
        enddo
    enddo

    !We can now diagonalize F
    call diagonalize_full(Aux3(1:Nvib,1:Nvib),Nvib,L(1:Nvib,1:Nvib),Freq(1:Nvib),"lapack")
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

end subroutine gf_method_V


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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine zmat2cart(molec,bond_s,angle_s,dihed_s,S,verbose)

    ! Transform z-matrix to cartesian coordinates
    ! Run with atomic coordinates (for S), but output in AMS!

    use structure_types
    use constants
    use MatrixMod

    integer,parameter :: NDIM=400
    logical,intent(in) :: verbose

    real(8) :: det

    type(str_resmol),intent(inout) :: molec
    integer,dimension(1:NDIM,1:4),intent(in) :: bond_s,angle_s,dihed_s
    real(8),dimension(NDIM),intent(in) :: S
   
    real(8) :: R,angle,dihed, xcom,ycom,zcom, MASS,xaux,yaux,zaux

    real(8),dimension(3,3) :: PIner,Aux, Aux3
    real(8),dimension(3) :: Vec

    real(8),dimension(NDIM,NDIM) :: work
    integer,dimension(NDIM) :: ipiv
    integer :: info

    integer :: Nat

    integer :: i, i_at, i_b,i_a,i_d

    Nat = molec%natoms

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
    if (verbose) then
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


!     print*, ""
!     print*, "COORDINATES IN XYZ"
!     print*, "==================="
!     print*, "Written to intermed.xyz"
    open(66,file="intermed.xyz")
    write(66,*) molec%natoms
    write(66,*) "final"
    do i=1,molec%natoms
        write(66,*) molec%atom(i)%name, molec%atom(i)%x,molec%atom(i)%y,molec%atom(i)%z
        !Convert back to a.u.
        molec%atom(i)%x = molec%atom(i)%x/BOHRtoAMS
        molec%atom(i)%y = molec%atom(i)%y/BOHRtoAMS
        molec%atom(i)%z = molec%atom(i)%z/BOHRtoAMS
    enddo
    close(66)
!     print*, "==================="
!     print*, ""
    return


end subroutine zmat2cart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine zmat2cart_ori(molec,bond_s,angle_s,dihed_s,S,T,info2,verbose)

    ! Transform z-matrix to cartesian coordinates
    ! Run with atomic coordinates (for S), but output in AMS!
    ! This version reuse the Rotation to std frame if requested

    use structure_types
    use constants
    use MatrixMod

    integer,parameter :: NDIM=400
    logical,intent(in) :: verbose

    real(8) :: det

    type(str_resmol),intent(inout) :: molec
    integer,dimension(1:NDIM,1:4),intent(in) :: bond_s,angle_s,dihed_s
    real(8),dimension(NDIM),intent(in) :: S
    real(8),dimension(1:3,1:3),intent(inout) :: T
    integer,intent(in) :: info2
   
    real(8) :: R,angle,dihed, xcom,ycom,zcom, MASS,xaux,yaux,zaux

    real(8),dimension(3,3) :: PIner,Aux, Aux3
    real(8),dimension(3) :: Vec

    real(8),dimension(NDIM,NDIM) :: work
    integer,dimension(NDIM) :: ipiv
    integer :: info

    integer :: Nat

    integer :: i, i_at, i_b,i_a,i_d

    Nat = molec%natoms

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
print*, "connected to i1", R
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
    if (info2 == 0) then
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
        if (verbose) then
        write(6,*) ""
        write(6,*) "I"
        do i=1,3
            write(6,'(100(F9.5,X))') PIner(i,1:3)
        enddo
        endif
        T(1:3,1:3) = Aux(1:3,1:3)
    else
        Aux(1:3,1:3) = T(1:3,1:3)
    endif

    write(6,*) ""
    write(6,*) "Xrot"
    do i=1,3
        write(6,'(100(F9.5,X))') Aux(i,1:3)
    enddo

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


!     print*, ""
!     print*, "COORDINATES IN XYZ"
!     print*, "==================="
!     print*, "Written to intermed.xyz"
    open(66,file="intermed.xyz")
    write(66,*) molec%natoms
    write(66,*) "final"
    do i=1,molec%natoms
        write(66,*) molec%atom(i)%name, molec%atom(i)%x,molec%atom(i)%y,molec%atom(i)%z
        !Convert back to a.u.
        molec%atom(i)%x = molec%atom(i)%x/BOHRtoAMS
        molec%atom(i)%y = molec%atom(i)%y/BOHRtoAMS
        molec%atom(i)%z = molec%atom(i)%z/BOHRtoAMS
    enddo
    close(66)
!     print*, "==================="
!     print*, ""
    return


end subroutine zmat2cart_ori

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine addcart(atom_1,atom_2,atom_3,&
                    bond,angle,dihed,atom_New)

! Angles in radians

    use structure_types
    use constants

    implicit none

    integer,parameter :: NDIM=400
    real(8),parameter :: ZERO = 1.d-10

    type(str_atom),intent(in) :: atom_1,atom_2,atom_3
    type(str_atom),intent(inout) :: atom_New 
    real(8),intent(in) :: bond,angle,dihed

    type(str_atom) :: aux_atom_1,aux_atom_2,aux_atom_3
    real(8) :: v1,v2,v3,sinTheta,Theta,baux,baux2,dihed_aux

    integer :: i

!     print*, " Input loactions"
!     print'(A,3F8.2)', "   Atom1 ", atom_1%x, atom_1%y, atom_1%z
!     print'(A,3F8.2)', "   Atom2 ", atom_2%x, atom_2%y, atom_2%z
!     print'(A,3F8.2)', "   Atom3 ", atom_3%x, atom_3%y, atom_3%z

    !New frame (only for the atoms involved)
    !1. Traslate to have atom i_b at origin
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

!     print*, " Std loactions"
!     print'(A,3F8.2)', "   C ", aux_atom_1%x, aux_atom_1%y, aux_atom_1%z
!     print'(A,3F8.2)', "   C ", aux_atom_2%x, aux_atom_2%y, aux_atom_2%z
!     print'(A,3F8.2)', "   C ", aux_atom_3%x, aux_atom_3%y, aux_atom_3%z
!     print'(A,3F8.2)', "   C ", atom_New%x, atom_New%y, atom_New%z

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

!     print*, " Final loactions"
!     print'(A,3F8.2)', "   C ", atom_1%x, atom_1%y, atom_1%z
!     print'(A,3F8.2)', "   C ", atom_2%x, atom_2%y, atom_2%z
!     print'(A,3F8.2)', "   C ", atom_3%x, atom_3%y, atom_3%z
!     print'(A,3F8.2)', "   C ", atom_New%x, atom_New%y, atom_New%z

    return

    contains

    subroutine rotation_3D(vx,vy,vz,tx,ty,tz,Theta) 

        !Description:
        ! Subroutine to rotate the vector (vx,vy,vz) around the axis
        ! defined by (tx,ty,tz) an angle Theta (rad).

        real(8), intent(inout) :: vx,vy,vz
        real(8), intent(in)    :: tx,ty,tz, Theta

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

       return

    end subroutine rotation_3D

end subroutine addcart

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine check_ori2(molec,molec2)

    ! Check if two molecules share the same orientation and, if not, reorinent molec to fit molec2
    ! This version "ensures" proper relative orientation (but is very time consuming!!)

    use structure_types

    implicit none

    type(str_resmol),intent(inout) :: molec
    type(str_resmol),intent(in)    :: molec2

    type(str_resmol)    :: molec_aux
    real(8),dimension(3,3) :: T
    real(8),dimension(3,3) :: T0
    real(8),dimension(3) :: Vec2

    integer :: iT, iTT, i, Nat
    real(8) :: xaux, yaux, zaux, det, rsum, rcheck, &
               xrot, yrot, zrot

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
        elseif (iTT ==6 ) then
        T(1,1) = Vec2(1)
        T(3,2) = Vec2(2)
        T(2,3) = Vec2(3)
        endif
        
    
        
!     print*, "T",iT, iTT
!     do i=1,3
!         print'(100(F8.3,2X))', T(i,1:3)
!     enddo

    det =       T(1,1)*T(2,2)*T(3,3)
    det = det + T(2,1)*T(3,2)*T(1,3)
    det = det + T(1,2)*T(2,3)*T(3,1)
    det = det - T(3,1)*T(2,2)*T(1,3)
    det = det - T(2,1)*T(1,2)*T(3,3)
    det = det - T(3,2)*T(2,3)*T(1,1)

!     print*, det

    if (det<0.d0) cycle

    !Rotate and check sum of distances
    molec_aux=molec
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
!     print*, "Sum of distance ", rsum

    if (rsum < rcheck) then
        T0=T
        rcheck=rsum
    endif

    enddo !
    enddo ! Possible T values

    print*, "T0"
    do i=1,3
        print'(100(F8.3,2X))', T0(i,1:3)
    enddo
    do i=1,Nat
        xaux = molec%atom(i)%x
        yaux = molec%atom(i)%y
        zaux = molec%atom(i)%z
        molec%atom(i)%x = T0(1,1)*xaux + T0(1,2)*yaux + T0(1,3)*zaux 
        molec%atom(i)%y = T0(2,1)*xaux + T0(2,2)*yaux + T0(2,3)*zaux
        molec%atom(i)%z = T0(3,1)*xaux + T0(3,2)*yaux + T0(3,3)*zaux 
    enddo
    print*, ""

    return

end subroutine check_ori2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine check_ori2b(molec,molec2,dist)

    ! Check if two molecules share the same orientation and, if not, reorinent molec to fit molec2
    ! This version "ensures" proper relative orientation (but is very time consuming!!)
    ! New: version b includes a distance threshold to avoid checking all the orientations

    use structure_types

    implicit none

    type(str_resmol),intent(inout) :: molec
    type(str_resmol),intent(in)    :: molec2
    real(8),intent(inout)          :: dist

    type(str_resmol)    :: molec_aux
    real(8),dimension(3,3) :: T
    real(8),dimension(3,3) :: T0
    real(8),dimension(3) :: Vec2

    integer :: iT, iTT, i, Nat
    real(8) :: xaux, yaux, zaux, det, rsum, rcheck, &
               xrot, yrot, zrot

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
        elseif (iTT ==6 ) then
        T(1,1) = Vec2(1)
        T(3,2) = Vec2(2)
        T(2,3) = Vec2(3)
        endif
        
    
        
!     print*, "T",iT, iTT
!     do i=1,3
!         print'(100(F8.3,2X))', T(i,1:3)
!     enddo

    det =       T(1,1)*T(2,2)*T(3,3)
    det = det + T(2,1)*T(3,2)*T(1,3)
    det = det + T(1,2)*T(2,3)*T(3,1)
    det = det - T(3,1)*T(2,2)*T(1,3)
    det = det - T(2,1)*T(1,2)*T(3,3)
    det = det - T(3,2)*T(2,3)*T(1,1)

!     print*, det

    if (det<0.d0) cycle

    !Rotate and check sum of distances
    molec_aux=molec
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
!     print*, "Sum of distance ", rsum

    if (rsum < rcheck) then
        T0=T
        rcheck=rsum
    endif

    if (rsum < dist) exit

    enddo !

    !If exited due to criterion met, exit here
    if (rsum < dist) exit

    enddo ! Possible T values

    !Update the distance criterion
    print*, "Distance criterion was: ", dist
    print*, "Updated to: ", rcheck
    dist=rcheck
 
    print*, "T0"
    do i=1,3
        print'(100(F8.3,2X))', T0(i,1:3)
    enddo
    do i=1,Nat
        xaux = molec%atom(i)%x
        yaux = molec%atom(i)%y
        zaux = molec%atom(i)%z
        molec%atom(i)%x = T0(1,1)*xaux + T0(1,2)*yaux + T0(1,3)*zaux 
        molec%atom(i)%y = T0(2,1)*xaux + T0(2,2)*yaux + T0(2,3)*zaux
        molec%atom(i)%z = T0(3,1)*xaux + T0(3,2)*yaux + T0(3,3)*zaux 
    enddo
    print*, ""

    return

end subroutine check_ori2b

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

subroutine check_ori3(molec,molec2,info)

    ! Check if two molecules share the same orientation and, if not, reorinent molec to fit molec2
    ! This version seeks effciency but may fail Needs the two structures to be only slightly different

    use structure_types

    implicit none

    type(str_resmol),intent(inout) :: molec
    type(str_resmol),intent(in)    :: molec2
    integer,intent(out) :: info

    type(str_resmol)    :: molec_aux
    real(8),dimension(3,3) :: T0, R
    real(8),dimension(3) :: Vec2

    integer :: i,j, Nat
    real(8) :: xaux, yaux, zaux, x2aux, y2aux, z2aux,&
               ux,uy,uz, umod, baux, baux2, Theta, costheta

    Nat = molec%natoms
    info=0

    !ASUMING THAT BOTH ARE A COM
    !Rotate atom1, so that it overlaps in both structures

    !The rotation axis is ...
    xaux  = (molec%atom(1)%x + &
             molec%atom(2)%x + &
             molec%atom(3)%x)/3.d0
    yaux  = (molec%atom(1)%y + &
             molec%atom(2)%y + &
             molec%atom(3)%y)/3.d0
    zaux  = (molec%atom(1)%z + &
             molec%atom(2)%z + &
             molec%atom(3)%z)/3.d0
    x2aux  = (molec2%atom(1)%x + &
             molec2%atom(2)%x + &
             molec2%atom(3)%x)/3.d0
    y2aux  = (molec2%atom(1)%y + &
             molec2%atom(2)%y + &
             molec2%atom(3)%y)/3.d0
    z2aux  = (molec2%atom(1)%z + &
             molec2%atom(2)%z + &
             molec2%atom(3)%z)/3.d0
    ux = yaux*z2aux - zaux*y2aux
    uy = zaux*x2aux - xaux*z2aux
    uz = xaux*y2aux - yaux*x2aux
    ! Vector u must be unitary
    umod = sqrt(ux**2 + uy**2 + uz**2)
    ux = ux/umod
    uy = uy/umod
    uz = uz/umod
    !And the angle is the one they form
    baux  = dsqrt(xaux**2+yaux**2+zaux**2)
    baux2 = dsqrt(x2aux**2+y2aux**2+z2aux**2)
    costheta = (xaux*x2aux + yaux*y2aux + zaux*z2aux)/baux/baux2
    Theta = dacos(costheta)

    !Set rotation matrix
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
    print*, "R"
    do i=1,3
        print'(100(F8.3,2X))', R(i,1:3)
    enddo

    !Then rotate to 

    !T rotation should be simular to R, but with 0/1 elements
    do i=1,3
        do j=1,3
!               if (dabs(R(i,j)) > 0.8d0) then
!                   T0(i,j) = dsign(1.d0,R(i,j))
!               else
!                   T0(i,j) = 0.d0
!               endif
            if (dabs(R(i,j)) > 0.04 .and. dabs(R(i,j)) < 0.96 .or. isnan(R(i,j))) then
                info=1
                print*, "Fast method failed"
                print*, ""
                return
            endif
            T0(i,j) = DFLOAT(NINT(R(i,j)))
        enddo
    enddo

    print*, "T0 (fast method)"
    do i=1,3
        print'(100(F8.3,2X))', T0(i,1:3)
    enddo
    do i=1,Nat
        xaux = molec%atom(i)%x
        yaux = molec%atom(i)%y
        zaux = molec%atom(i)%z
        molec%atom(i)%x = T0(1,1)*xaux + T0(1,2)*yaux + T0(1,3)*zaux 
        molec%atom(i)%y = T0(2,1)*xaux + T0(2,2)*yaux + T0(2,3)*zaux
        molec%atom(i)%z = T0(3,1)*xaux + T0(3,2)*yaux + T0(3,3)*zaux 
    enddo
    print*, ""

    return

end subroutine check_ori3

subroutine check_Ci(molec,molecP,molec2)

    !Compare which one, molec or molecP, is more similar to molec2
    ! The more similar is output on molec

    use structure_types

    type(str_resmol),intent(inout) :: molec
    type(str_resmol),intent(in)    :: molec2, molecP

    real(8) :: rmsd1,rmsd2,d, &
               tx,ty,tz

    
    !
    rmsd1=0.d0
    rmsd2=0.d0   
    do i=1,molec%natoms
        d=((molec%atom(i)%x-molec2%atom(i)%x)**2 + &
           (molec%atom(i)%y-molec2%atom(i)%y)**2 + &
           (molec%atom(i)%z-molec2%atom(i)%z)**2)
        rmsd1=rmsd1 + dsqrt(d)
        d=((molecP%atom(i)%x-molec2%atom(i)%x)**2 + &
           (molecP%atom(i)%y-molec2%atom(i)%y)**2 + &
           (molecP%atom(i)%z-molec2%atom(i)%z)**2)
        rmsd2=rmsd2 + dsqrt(d)
    enddo
!       print*, "RMSD", rmsd1, rmsd2
    if (rmsd2 < rmsd1) then
        print*, "INVERSION OCCURRED"
        molec=molecP  
    endif 

    return

end subroutine check_Ci

subroutine NumBDer(molec,S_sym,Bder)

!     use parameters
    use structure_types

    integer,parameter :: NDIM = 400
    real(8),parameter :: delta = 1.889726133d-3 !for numerical ders, in bohr(=10^-3 \AA, as Num freq in G09)

    type(str_resmol),intent(in) :: molec
    integer,dimension(NDIM),intent(in) :: S_sym
    real(8),dimension(:,:,:),intent(out) :: Bder     ! (Nvib x 3Nat x 3Nat) the last index is the second der


    type(str_resmol) :: molecB
    real(8),dimension(NDIM,NDIM) :: Bplus, Bmin
    !Unused stuff that go to the Wilson SR
    !Input-Output
    real(8),dimension(NDIM,NDIM) :: Asel !Used to manage nonredundant internal variables (disabled option)
    !Output
    real(8),dimension(NDIM,NDIM) :: G          !G matrix
    real(8),dimension(NDIM) :: S                  !Vector of internal coordinates
    character(len=100),dimension(NDIM) :: ModeDef !Definition of normal modes

    !input
    logical :: verbose=.false.
!     real(8),dimension(NDIM,NDIM) :: 
    
    !deactivate this option
    Asel(1,1) = 99.d0


    Nat = molec%natoms
    Nvib = 3*molec%natoms-6

    print*, ""
    print*, "COMPUTING NUMERICAL DERIVATIVES FOR B..."
    print*, ""

    do i=1,Nat

        !Displace X
        molecB = molec
        ii = 3*i-2
! print*, i, "X0", molecB%atom(i)%x
        molecB%atom(i)%x = molec%atom(i)%x + delta
! print*, i, "X0+d", molecB%atom(i)%x
        !Call B matrix at this geometry   
        call internal_Wilson(molecB,S,S_sym,ModeDef,Bplus,G,Asel,verbose)
        molecB%atom(i)%x = molec%atom(i)%x - delta
! print*, i, "X0-d", molecB%atom(i)%x
        !Call B matrix at this geometry   
        call internal_Wilson(molecB,S,S_sym,ModeDef,Bmin,G,Asel,verbose)

        do j=1,Nvib
        do k=1,3*Nat
           Bder(j,k,ii) = ( Bplus(j,k) - Bmin(j,k) ) / (2.d0*delta)
print*, j,k,ii, Bder(j,k,ii), Bplus(j,k), Bmin(j,k)
        enddo
        enddo

        !Displace Y
        molecB = molec
        ii = 3*i-1
        molecB%atom(i)%y = molec%atom(i)%y + delta
        !Call B matrix at this geometry   
        call internal_Wilson(molecB,S,S_sym,ModeDef,Bplus,G,Asel,verbose)
        molecB%atom(i)%y = molec%atom(i)%y - delta
        !Call B matrix at this geometry   
        call internal_Wilson(molecB,S,S_sym,ModeDef,Bmin,G,Asel,verbose)

        do j=1,Nvib
        do k=1,3*Nat
           Bder(j,k,ii) = ( Bplus(j,k) - Bmin(j,k) ) / (2.d0*delta)
print*, j,k,ii, Bder(j,k,ii), Bplus(j,k), Bmin(j,k)
        enddo
        enddo


        !Displace Z
        molecB = molec
        ii = 3*i
        molecB%atom(i)%z = molec%atom(i)%z + delta
        !Call B matrix at this geometry   
        call internal_Wilson(molecB,S,S_sym,ModeDef,Bplus,G,Asel,verbose)
        molecB%atom(i)%z = molec%atom(i)%z - delta
        !Call B matrix at this geometry   
        call internal_Wilson(molecB,S,S_sym,ModeDef,Bmin,G,Asel,verbose)

        do j=1,Nvib
        do k=1,3*Nat
           Bder(j,k,ii) = ( Bplus(j,k) - Bmin(j,k) ) / (2.d0*delta)
print*, j,k,ii, Bder(j,k,ii), Bplus(j,k), Bmin(j,k)
        enddo
        enddo

    enddo
stop

    return

end subroutine NumBDer


end module internal_module
