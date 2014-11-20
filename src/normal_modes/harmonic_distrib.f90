program har_dist


    !==============================================================
    ! This code uses of MOLECULAR_TOOLS (version 1.0/March 2012)
    !==============================================================
    !
    ! Description:
    ! -----------
    ! Program to analyse vibrations in term of internal coordinates.
    !
    ! Compilation instructions (for mymake script):
    !
    ! Change log:
    !
    ! TODO:
    ! ------
    !
    ! History
    ! v0: Adapted from nm_internal (v4.0.1.1)
    !
    !============================================================================    

!*****************
!   MODULE LOAD
!*****************
!============================================
!   Generic (structure_types independent)
!============================================
!     use alerts
!     use line_preprocess
    use constants
!   Matrix manipulation (i.e. rotation matrices)
!     use MatrixMod
!============================================
!   Structure types module
!============================================
!     use structure_types
!============================================
!   Structure dependent modules
!============================================
!     use gro_manage
!     use pdb_manage
!     use gaussian_manage
!     use gaussian_fchk_manage
!     use xyz_manage
!     use molcas_unsym_manage
!   Structural parameters
!     use molecular_structure
!     use ff_build
!   Bond/angle/dihed meassurement
!     use atomic_geom
!   Symmetry support
!     use symmetry_mod
!   For internal thingies
!     use internal_module
!     use zmat_manage

    implicit none

    real(8) :: KbT, T, g1, varB, varQ, freq
    integer :: Nvib, i

    !===================
    !CPU time 
    real(8) :: ti, tf
    !===================

! (End of variables declaration) 
!==================================================================================

    ! 0. GET COMMAND LINE ARGUMENTS
    read(*,*) T
    read(*,*) Nvib
    
    KbT=boltz*T / HARTtoJ
    do i=1,Nvib
        read(*,*) g1
        freq=g1
        g1 = g1 / autown
        varQ=dsqrt( 1.d0/(2.d0*dtanh(g1/2.d0/KbT)) )
        varB=dsqrt(KbT/g1)
        print*, i, varB, varQ
    enddo



    stop       

end program har_dist

! Try from SI to AU: 
!     do i=1,Nvib
!         read(*,*) g1
!         !Factor to adim (SI)
!         Factor = dsqrt(dabs(g1)*1.d2*clight*2.d0*PI/plankbar)
!         !====
!         g1 = g1 / autown
! !         varQ=1.d0/(dtanh(0.5d0*g1/KbT)*2.d0*g1)
!         varQ=1.d0/(dtanh(0.5d0*g1/KbT)*2.d0*dsqrt(KbT))
!         varB=dsqrt(KbT)/g1
!         !To SI
!         varQ=varQ*BOHRtoM*dsqrt(AUtoKG)
!         varB=varB*BOHRtoM*dsqrt(AUtoKG)
!         varQ=varQ*Factor
!         varB=varB*Factor
!         print*, i, varQ, varB
!     enddo

