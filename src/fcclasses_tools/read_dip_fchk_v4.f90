program read_dip_fchk

    !Compilation
    ! $FC ../modules/alerts.f90 ../modules/constants_mod.f90 ../modules/symmetry_mod_notypes.f90 ../modules/line_preprocess.f90 ../modules/gaussian_fchk_manage_v2.1.f90 read_dip_fchk_v2.f90 -o read_dip_fchk_v2.exe -cpp -DDOUBLE 

    !NOTES

    ! History
    ! Version 1: Initial Release
    ! Version 2: Reads geometry and guess symmetry in case symmetry was used to save computing
    !            numerical freqs. (so, some data would be missing)
    ! V4: addpated to V4 modules (no action is required, as it uses "no_types" versions)

    use gaussian_fchk_manage
    use alerts
    use line_preprocess
    use constants
    use symmetry_mod_notypes

    !Interesting info..
    integer :: Nat, Nes
    character(len=250) :: jobname
    character(len=20) :: jobtype,method,basis
    real(8),dimension(1:3) :: Dip
    real(8),dimension(:),allocatable :: DipD, geom
    integer,dimension(:),allocatable :: isym

    !Auxiliars
    integer :: N
    character :: dtype, cnull
    real(8),dimension(:),allocatable :: A
    integer,dimension(:),allocatable :: I
    integer :: error
    !Selection stuff
    integer :: Nsel 
    integer,dimension(1:1000) :: nm_select
    !Counters
    integer :: j,k, jj, j_sym, jj_sym, &
               at_num, at_sym, delta

    !I/O
    character(len=100) :: fchkfile, outfile
    integer :: ios
    integer :: I_FCHK=10,&
               O_DIP =20

    ! Get filenames from commandline or read standard input/output
    if ( iargc() == 0 ) then
        I_FCHK = 5
    else
        call getarg(1, fchkfile)
        open(I_FCHK,file=fchkfile,status="old",iostat=ios)
        if ( ios /= 0 ) call alert_msg("fatal","could not open the file "//trim(adjustl(fchkfile)))
    endif

    !INFORMATION IN THE FCHK FILE
    ! Natoms
    call read_fchk(I_FCHK,"Number of atoms",dtype,N,A,I,error)
    if (error == 0) then
        Nat = I(1)
        allocate(isym(1:Nat))
        deallocate(I)
    endif

    ! Cartesian geometry (a.u.)
    call read_fchk(I_FCHK,"Current cartesian coordinates",dtype,N,A,I,error)
    if (error == 0) then
        allocate(geom(1:N))
        geom(1:N) = A(1:N)
        deallocate(A)
    endif

    ! Number of excited states computed
    call read_fchk(I_FCHK,"ETran scalars",dtype,N,A,I,error)
    if (error == 0) then
        Nes = I(1)
        deallocate(I)
    else
        print*, "ERROR: 'ETran scalars' section not found"
        stop
    endif

    ! Read ETran state values
    call read_fchk(I_FCHK,"ETran state values",dtype,N,A,I,error)
    if (error /= 0) then
        print*, "ERROR: 'ETran state values' section not found"
        stop
    endif
    close(I_FCHK)

    !Set to zero very low values
!     do j=1,N 
!         if (abs(A(j)) .lt. 1d-10) A(j)=0.d0
!     enddo

    !===================================
    ! TRANSITION ELECTRIC DIPOLE MOMENT
    !===================================
    ! Get Dip (always present for TD/CIS calculations)
    ! NOTE: we get the first root: should be changed to allow other roots
    Dip(1:3) = A(2:4)

    ! If numerical freqs were computed, the size of the section is bigger:
        call symm_atoms_nt(geom,Nat,isym)
    if (N == Nes*16+48+3*Nat*16) then
        print*, "Derivatives of TrDip available"
        call symm_atoms_nt(geom,Nat,isym)
        allocate(DipD(1:3*3*Nat))
        !Note the order is not the same as needed for FCclasses
        ! FCHK: dmux/dx1, dmuy/dx1, dmuz/dx1,  dmux/dy1, dmuy/dy1, dmuz/dy1,  ...
        ! FCcl: dmux/dx1, dmux/dy1, dmux/dz1,  dmuy/dx1, dmuy/dy1, dmuy/dz1,  ...
        do j=1,3*Nat
            k = 16*Nes+48 + 16*(j-1)
            jj = j*3-2
            DipD(jj:jj+2) = A(k+2:k+4)

!             if (DipD(jj  ) == 0.d0 .and.&
!                 DipD(jj+1) == 0.d0 .and.&
!                 DipD(jj+2) == 0.d0) then
!                 !use symmetric values
!                 at_num = j/3
!                 if (at_num*3 /= j) at_num=at_num+1 
!                 delta=at_num*3 - j
!                 j_sym = 3*isym(at_num)-delta
! print'(X,4(A,I3))', "Taking symmetric for atom ", at_num, " from atom ", isym(at_num), "  |  Ders for geom: ", j, " <-->", j_sym
!                 jj_sym = j_sym*3-2
!                 DipD(jj:jj+2) = DipD(jj_sym:jj_sym+2)
!             endif
            if (DipD(jj  ) == 0.d0 .and.&
                DipD(jj+1) == 0.d0 .and.&
                DipD(jj+2) == 0.d0) then
                !Symmetry is used
                call alert_msg("warning","Looks like the computation was done with symmetry. If so, dipole ders. are not reliable!")
            endif
        enddo
    else
        print*, "Derivatives of TrDip NOT available"
    endif
    !WRITE ELDIP FILE
    ! outfile
    if ( I_FCHK /= 5 ) then
        call split_line(fchkfile,".",outfile,cnull)
        outfile = "eldip_"//trim(adjustl(outfile))//"_fchk"
    else
        outfile = "eldip_fchk"
    endif
    open(O_DIP,file=outfile,status="replace")
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)
    if (N == Nes*16+48+3*Nat*16) then
        do j=1,3*Nat
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))') DipD(jj:jj+2)
        enddo
    endif
    close(O_DIP)


    !===================================
    ! TRANSITION MAGNETIC DIPOLE MOMENT
    !===================================
    ! Get Dip (always present for TD/CIS calculations)
    ! NOTE: we get the first root: should be changed to allow other roots
    !Note that MagDip should be divided by -2 (see FCclasses manual)
    Dip(1:3) = A(8:10)/(-2.d0)

    ! If numerical freqs were computed, the size of the section is bigger
    if (N == Nes*16+48+3*Nat*16) then
        do j=1,3*Nat
            k = 16*Nes+48 + 16*(j-1)
            jj = j*3-2
            !Note that MagDip should be divided by -2 (see FCclasses manual)
            DipD(jj:jj+2) = A(k+8:k+10)/(-2.d0)
        enddo
    endif
    !WRITE MAGDIP FILE
    ! outfile
    if ( I_FCHK /= 5 ) then
        call split_line(fchkfile,".",outfile,cnull)
        outfile = "magdip_"//trim(adjustl(outfile))//"_fchk"
    else
        outfile = "magdip_fchk"
    endif
    open(O_DIP,file=outfile,status="replace")
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)
    write(O_DIP,'(3(X,E18.9))') Dip(1:3)
    if (N == Nes*16+48+3*Nat*16) then
        do j=1,3*Nat
            jj=j*3-2
            write(O_DIP,'(3(X,E18.9))') DipD(jj:jj+2)
        enddo
        !Its time to deallocate
        deallocate(DipD)
    endif
    close(O_DIP)

    deallocate(A)

    stop

    contains


end program read_dip_fchk

