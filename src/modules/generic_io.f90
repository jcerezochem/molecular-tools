module generic_io

    !==============================================================
    ! This code is part of FCC_TOOLS 
    !==============================================================
    ! Description
    !  This MODULE contains subroutines to manage output files 
    !  from different QM codes. Relies in specific modules for
    !  each package.
    !    
    !==============================================================

    !Common declarations:
    !===================
    use constants
    use line_preprocess
    use gaussian_manage
    use gamess_manage
    use psi4_manage
    use molcas_manage
    use molpro_manage
    use gmx_manage
    use pdb_manage
    use fcc_manage
    use xyz_manage
    implicit none

    contains

    subroutine assign_masses(Nat,AtName,Mass)

        integer,intent(in) :: Nat
        character(len=*),dimension(:),intent(in) ::  AtName
        real(8),dimension(:),intent(out) ::  Mass

        !Local
        integer :: i
        
        do i=1,Nat
            Mass(i) = atmass_from_atname(AtName(i)) 
        enddo

        return

    end subroutine assign_masses

    subroutine generic_natoms_reader(unt,filetype,Nat,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic natoms reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (out)  int /scalar   Number of atoms
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)           :: unt
        character(len=*),intent(in)  :: filetype
        integer,intent(out)          :: Nat
        integer,intent(out),optional :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        integer                          :: error_local

        error_local = 0
        select case (adjustl(filetype))
            case("log")
             call read_gausslog_natoms(unt,Nat,error_local)
            case("fchk")
             call read_fchk(unt,'Number of atoms',data_type,N,A,IA,error_local)
             Nat = IA(1)
             deallocate(IA)
            case("gms")
             call read_gamess_natoms(unt,Nat,error_local)
            case("psi4")
             call read_psi4_natoms(unt,Nat,error_local)
            case("molcas")
             call read_molcasUnSym_natoms(unt,Nat,error_local)
            case("molpro")
             call read_molpro_natoms(unt,Nat,error_local)
            case("g96")
             call read_g96_natoms(unt,Nat)
            case("xyz")
             call read_xyz_natoms(unt,Nat)
            case default
             call alert_msg("fatal","Unsupported filetype:"//trim(adjustl(filetype)))
!              call supported_filetype_list('freq')
             error_local = 99
         end select

         if (present(error_flag)) error_flag=error_local

         return

    end subroutine generic_natoms_reader


    subroutine generic_structure_reader(unt,filetype,Nat,X,Y,Z,Mass,AtName,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic geometry reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (io )  int /scalar   Number of atoms
        ! X,Y,Z   (out)  real/vectors  Coordinates
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(inout)           :: Nat
        real(8),dimension(:),intent(out):: X,Y,Z
        ! Made inout to handle incomplete files (fcc)
        character(len=*),dimension(:),intent(inout) :: AtName
        real(8),dimension(:),intent(inout):: Mass
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        integer,dimension(size(AtName))  :: AtNum
        !Other local
        integer                          :: i,j
        integer                          :: error_local

        error_local = 0
        select case (adjustl(filetype))
            case("log")
             call read_gauslog_geom(unt,Nat,AtName,X,Y,Z,error_local)
             call assign_masses(Nat,AtName,Mass)
            case("log-inpori")
             call read_gauslog_stdori(unt,Nat,AtNum,X,Y,Z,"Input orientation",error_local)
             do i=1,Nat
                  AtName(i) = atname_from_atnum(AtNum(i))
                  Mass(i)   = atmass_from_atnum(AtNum(i))
             enddo
            case("log-stdori")
             call read_gauslog_stdori(unt,Nat,AtNum,X,Y,Z,"Standard orientation",error_local)
             do i=1,Nat
                  AtName(i) = atname_from_atnum(AtNum(i))
                  Mass(i)   = atmass_from_atnum(AtNum(i))
             enddo
            case("fchk")
             call read_fchk(unt,'Current cartesian coordinates',data_type,N,A,IA,error_local)
             do i=1,N,3
                 j = (i-1)/3+1
                 X(j) = A(i)  *BOHRtoANGS
                 Y(j) = A(i+1)*BOHRtoANGS
                 Z(j) = A(i+2)*BOHRtoANGS
             enddo
             deallocate(A)
             Nat = N/3
             call read_fchk(unt,'Real atomic weights',data_type,N,A,IA,error_local)
             Mass(1:N) = A(1:N)
             deallocate(A)
             call read_fchk(unt,'Atomic numbers',data_type,N,A,IA,error_local)
             do i=1,N
                 AtName(i) = atname_from_atnum(IA(i))
                 AtName(i) = adjustl(AtName(i))
             enddo
             deallocate(IA)
            case("gms")
             call read_gamess_geom(unt,Nat,AtName,X,Y,Z,error_local)
             call assign_masses(Nat,AtName,Mass)
            case("psi4")
             call read_psi4_geom(unt,Nat,AtName,X,Y,Z,error_local)
             call assign_masses(Nat,AtName,Mass)
            case("molcas")
             call read_molcasUnSym_geom(unt,Nat,AtName,X,Y,Z,error_local)
             call assign_masses(Nat,AtName,Mass)
            case("molpro")
             call read_molpro_geom(unt,Nat,AtName,X,Y,Z,error_local)
             call assign_masses(Nat,AtName,Mass)
            case("g96")
             call read_g96_geom(unt,Nat,AtName,X,Y,Z)
             call assign_masses(Nat,AtName,Mass)
            case("xyz")
             call read_xyz_geom(unt,Nat,AtName,X,Y,Z)
             call assign_masses(Nat,AtName,Mass)
            case("pdb")
             call read_pdb_geom(unt,Nat,AtName,X,Y,Z)
             call assign_masses(Nat,AtName,Mass)
            case("fcc")
             call read_fccstate_geom(unt,Nat,X,Y,Z)
            case default
             call alert_msg("fatal","Unsupported filetype:"//trim(adjustl(filetype)))
!              call supported_filetype_list('freq')
             error_local = 99
         end select

         if (present(error_flag)) error_flag=error_local

         return

    end subroutine generic_structure_reader

    subroutine generic_structure_writer(unt,filetype,Nat,X,Y,Z,Mass,AtName,error_flag,title)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic geometry reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (io )  int /scalar   Number of atoms
        ! X,Y,Z   (out)  real/vectors  Coordinates
        ! error_flag (out) flag        0: Success
        !                              1: 
        ! title   (inp,opt) char/scalar title of the file
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(in)              :: Nat
        real(8),dimension(:),intent(in) :: X,Y,Z
        real(8),dimension(:),intent(in) :: Mass
        character(len=*),dimension(Nat),intent(in) :: AtName
        integer,intent(out),optional    :: error_flag
        character(len=*),intent(in),optional :: title

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        !Other local
        integer                          :: i,j
        integer                          :: error_local

        error_local = 0
        select case (adjustl(filetype))
            case("pdb")
             call write_pdb_geom(unt,Nat,AtName,X,Y,Z)
            case("fchk")
             if (present(title)) then
                 write(unt,'(A)') trim(adjustl(title))
             else
                 write(unt,'(A)') "FCHK generated with generic_io module"
             endif
             call write_fchk(unt,"Number of atoms","I",0,A,(/Nat/),error_local)
             ! Get atom number from atom names 
             allocate(IA(1:Nat))
             do i=1,Nat
                 IA(i) = atnum_from_atname(AtName(i))
             enddo
             call write_fchk(unt,"Atomic numbers","I",Nat,A,IA,error_local)
             call write_fchk(unt,"Nuclear charges","R",Nat,dfloat(IA),IA,error_local)
             deallocate(IA)
             call write_fchk(unt,"Real atomic weights","R",Nat,Mass,IA,error_local)
             ! Get coordinates as a vector
             allocate(A(1:3*Nat))
             j=0
             do i=1,Nat
                 j=j+1
                 A(j)=X(i)/BOHRtoANGS
                 j=j+1
                 A(j)=Y(i)/BOHRtoANGS
                 j=j+1
                 A(j)=Z(i)/BOHRtoANGS
             enddo   
             call write_fchk(unt,"Current cartesian coordinates","R",3*Nat,A,IA,error_local)
             deallocate(A)
            case default
             call alert_msg("fatal","Unsupported filetype:"//trim(adjustl(filetype)))
!              call supported_filetype_list('freq')
             error_local = 99
         end select

         if (present(error_flag)) error_flag=error_local

         return

    end subroutine generic_structure_writer


    subroutine generic_gradient_reader(unt,filetype,Nat,Grad,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic gradient reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (int)  int /scalar   Number of atoms
        ! Grad    (out)  real/vector   Gradient (AU)
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(in)              :: Nat
        real(8),dimension(:),intent(out):: Grad
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        !For gausslog read
        character(12*3*Nat)              :: section
        !Other auxiliar
        integer                          :: i
        integer                          :: error_local

        error_local = 0
        select case (adjustl(filetype))
            case("log")
             call summary_parser(unt,7,section,error_local)
             if (error_local /= 0) then
                 call alert_msg("warning","Gradient could not be read")
                 if (present(error_flag)) error_flag=error_local
                 return
             endif
             read(section,*) Grad(1:3*Nat)
            case("fchk")
             call read_fchk(unt,'Cartesian Gradient',data_type,N,A,IA,error_local)
             if (error_local /= 0) then
                 if (present(error_flag)) error_flag=error_local
                 return
             endif
             do i=1,N
                 Grad(i) = A(i)
             enddo
             deallocate(A)
            case default
             call alert_msg("fatal","Unsupported filetype:"//trim(adjustl(filetype)))
!              call supported_filetype_list('grad')
             error_local = 99
         end select

         if (present(error_flag)) error_flag=error_local

         return

    end subroutine generic_gradient_reader

    subroutine generic_Hessian_reader(unt,filetype,Nat,Hlt,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic Hessian reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (int)  int /scalar   Number of atoms
        ! Hlt     (out)  real/vector   Hessian (lower triangular form) (AU)
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(in)              :: Nat
        real(8),dimension(:),intent(out):: Hlt
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        !For gausslog read
        character(12*9*Nat*Nat)          :: section
        !Other auxiliar
        integer                          :: i
        integer                          :: error_local

        error_local = 0
        select case (adjustl(filetype))
            case("log")
             call summary_parser(unt,6,section,error_local)
             read(section,*) Hlt(1:3*Nat*(3*Nat+1)/2)
            case("fchk")
             call read_fchk(unt,'Cartesian Force Constants',data_type,N,A,IA,error_local)
             if (error_local /= 0) then
                 if (present(error_flag)) error_flag=error_local
                 return
             endif
             do i=1,N
                 Hlt(i) = A(i)
             enddo
             deallocate(A)
            case("gms")
             call read_gamess_hess(unt,Nat,Hlt,error_local)
            case("psi4")
             call read_psi4_hess(unt,Nat,Hlt,error_local)
            case("molcas")
             call read_molcasUnSym_hess(unt,Nat,Hlt,error_local)
            case("molpro")
             call read_molpro_hess(unt,Nat,Hlt,error_local)
            case("gmx")
             call read_gmx_hess(unt,Nat,Hlt,error_local)
            case default
             call alert_msg("fatal","Unsupported filetype:"//trim(adjustl(filetype)))
!              call supported_filetype_list('freq')
             error_local = 99
         end select

         if (present(error_flag)) error_flag=error_local

         return

    end subroutine generic_Hessian_reader


    subroutine generic_dip_reader(unt,filetype,Si,Sf,derivatives,dip_type,dx,Dip,DipD,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic dipole reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (io )  int /scalar   Number of atoms
        ! X,Y,Z   (out)  real/vectors  Coordinates
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(inout)           :: Si, Sf
        logical,intent(inout)           :: derivatives
        character(len=*),intent(in)     :: dip_type
        real(8),intent(inout)           :: dx !Only for numerical diffs
        real(8),dimension(:),intent(out):: Dip 
        real(8),dimension(:),intent(out):: DipD
        integer,intent(out),optional    :: error_flag
        !Local
        integer :: S
        integer :: error_local
        character(len=3) :: dummy_char

        error_local = 0
        select case (adjustl(filetype))
            case("log")
             !Get target state
             call read_gausslog_targestate(unt,S,error_local)
             if (present(error_flag)) error_flag=error_local
             if (Si == -1) Si = 0
             if (Sf == -1) Sf = S
             if (Si /= 0) then
                 write(dummy_char,'(I0)') Si
                 call alert_msg("fatal","TD-DFT calcs in G09 only provide trdip from/to GS, but requested S="//dummy_char)
                 error_local=-1
             else
                 !Need to rewind to read the first dip
                 rewind(unt)
                 !Now read dip
                 call read_gausslog_dip(unt,Si,Sf,dip_type,Dip,error_local)
                 if (derivatives) then
                     print*, " Computing derivatives.."
                     call read_gausslog_dipders(unt,Si,Sf,dip_type,dx,DipD,error_local)
                     if (error_local /= 0) derivatives=.false.
                     !Done this, we can safely reset error_flag to 0
                     error_local = 0
                 endif
             endif
            case("fchk")
             call read_gaussfchk_dip(unt,Si,Sf,derivatives,dip_type,Dip,DipD,error_local)
            case("gms")
             call alert_msg("fatal","Filetype not supported")
            case("psi4")
             call read_psi4_dip(unt,Si,Sf,dip_type,Dip,error_local)
             if (derivatives) then
                 print*, " Computing derivatives.."
                 call read_psi4_dipders(unt,Si,Sf,dip_type,dx,DipD,error_local)
                 if (error_local /= 0) derivatives=.false.
                 !Done this, we can safely reset error_flag to 0
                 error_local = 0
             endif
            case("molcas")
             call alert_msg("fatal","Filetype not supported")
            case("molpro")
             call alert_msg("fatal","Filetype not supported")
            case default
             call alert_msg("fatal","Unsupported filetype:"//trim(adjustl(filetype)))
!              call supported_filetype_list('freq')
             error_local = 99
         end select
         if (present(error_flag)) error_flag=error_local

         return

    end subroutine generic_dip_reader


    subroutine generic_nm_reader(unt,filetype,Nat,Nvib,Freq,L,error_flag)

        !==============================================================
        ! This code is part of FCC_TOOLS
        !==============================================================
        !Description
        ! Generic Hessian reader, using the modules for each QM program
        !
        !Arguments
        ! unt     (inp)  int /scalar   Unit of the file
        ! filetype(inp)  char/scalar   Filetype  
        ! Nat     (int)  int /scalar   Number of atoms
        ! Freq     (out)  real/vector  Frequencies (cm-1)
        ! L        (out)  real/matrix  Normal modes (Cartesian normalized, adim)
        ! error_flag (out) flag        0: Success
        !                              1: 
        !
        !==============================================================

        integer,intent(in)              :: unt
        character(len=*),intent(in)     :: filetype
        integer,intent(in)              :: Nat
        ! gauss files get Nvib here, fcc need it externally
        integer,intent(inout)           :: Nvib
        real(8),dimension(:),intent(out):: Freq
        real(8),dimension(:,:),intent(out):: L
        integer,intent(out),optional    :: error_flag

        !Local
        !Variables for read_fchk
        real(8),dimension(:),allocatable :: A
        integer,dimension(:),allocatable :: IA
        character(len=1)                 :: data_type
        integer                          :: N
        !Other auxiliar
        integer                          :: i, j, k
        integer                          :: error_local

        rewind(unt)

        error_local = 0
        select case (adjustl(filetype))
            case("log")
             call read_glog_nm(unt,Nvib,Nat,Freq,L,error_local)
            case("fchk")
             call read_fchk(unt,'Number of Normal Modes',data_type,N,A,IA,error_local)
             if (error_local /= 0) then
                 call alert_msg("fatal","The fchk does not contain normal modes")
                 if (present(error_flag)) error_flag=error_local
                 return
             endif
             Nvib=IA(1)
             deallocate(IA)
             call read_fchk(unt,'Vib-Modes',data_type,N,A,IA,error_local)
             if (error_local /= 0) then
                 if (present(error_flag)) error_flag=error_local
                 return
             endif
             ! Reconstruct Lcart 
             k=0
             do i=1,Nvib
                 do j=1,3*Nat
                     k=k+1
                     L(j,i) = A(k)
                 enddo
             enddo
             deallocate(A)
             call read_fchk(unt,'Vib-E2',data_type,N,A,IA,error_local)
             Freq(1:Nvib) = A(1:Nvib)
             deallocate(A)
            case("fcc")
             call read_fccstate_nm(unt,Nvib,Nat,Freq,L)
            case default
             call alert_msg("fatal","Unsupported filetype:"//trim(adjustl(filetype)))
!              call supported_filetype_list('freq')
             error_local = 99
         end select

         if (present(error_flag)) error_flag=error_local

         return

    end subroutine generic_nm_reader


end module generic_io

