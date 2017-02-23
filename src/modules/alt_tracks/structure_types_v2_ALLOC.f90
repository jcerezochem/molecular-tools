module structure_types

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    !Description
    ! This module contains SIZES of the program and the required
    ! derived types to manage molecular systems
    !==============================================================
    !---------------------------------------------------------
    ! Shoul work for GRO and PDB formats even if
    ! character lengths are different (the larger
    ! is used). If needed, use ADJUSTL(string).
    !
    ! NOTES
    ! -----
    ! Limits:
    ! MaxAtoms/system = 100 000 (maximum since counter is I5)
    ! MaxAtoms/residue = 200 (is safe)
    !
    ! HISTORY
    ! -------
    ! (16/11/11) Added charge to str_atom
    ! (21/11/11) Added connections and nbonds to str_atoms (beta)
    ! (19/01/12) Added normal modes and electronic properties stuff
    ! **VERSION 2**
    ! (27/01/12) New version release (v2): major chage
    !            Atom coordinates are stored in an array R(1:3) to 
    !            enhance algebra operations. Scalar x,y,z are
    !            kept for backward compability (to be deleted in the 
    !            next release) 
    ! (29/02/12) Allocatable arrays
    !            prop and nm are integrated in str_system
    !            bonded params are moved from atoms to geom (except environment)
    ! (03/05/12) Reverted change: atom coordinates are better stored
    !            as x,y,z. Tip: X(real array) = atom(:)%x
    ! (11/06/12) Splitted v2 in allocatable (_ALLOC) and not allocatable
    !            versions (this is allocatable)
    ! (13/06/12) Added fftype (character(len=6)) to atom. This stores the atom
    !            type in the force field, while attype stores the PDB attype
    ! (14/11/12) Added PG (point group) to system
    ! (09/07/13) Initialize string variables (long term goal)
    !---------------------------------------------------------

    !SIZES
    integer,parameter:: MAX_ATOMS = 100000, &
                        MAX_ATM_RES = 200,  &
                        MAX_CONNEXIONS = 6, &
                        MAX_DERIVED_CNX = 5000


    !DERIVED TYPES

    !============================================
    ! STR_ATOM
    !============================================
    type str_atom
        character::attype*6
        character::fftype*6
        character::name*5
        character::resname*5
        character::chain*1,alter_loc*1, ins_code*1
        integer::resseq, AtNum
#ifdef DOUBLE
        !-------------------------------------------------------------
        double precision:: x, y, z
        double precision,dimension(1:3) :: R !deprecated feature. not a good idea
        !-------------------------------------------------------------
        double precision:: mass, q
#else
        !-------------------------------------------------------------
        real:: x, y, z 
        real,dimension(1:3) :: R !deprecated feature. not a good idea
        !-------------------------------------------------------------
        real:: mass, q
#endif
        integer,dimension(:),allocatable :: connect
        integer :: nbonds
        character(len=6),dimension(:),allocatable :: env 
    end type str_atom

    type str_geom
        integer,dimension(:,:),allocatable :: bond, pair, angle, dihed
        integer :: npairs, nangles, ndihed, nbonds
    end type str_geom

    !============================================
    ! STR_SYSTEM
    !============================================
    type str_system
        !Atoms, fragments and system attributes
        integer::natoms
        integer::nres
        type(str_atom),dimension(:),allocatable::atom
        character(len=100) :: title=""
        character::PG*5=""
#ifdef DOUBLE
        double precision::boxX,boxY,boxZ
        !Electronic properties (GS and TD)
        double precision,dimension(1:3) :: dip, trans_dip
        double precision,dimension(:,:),allocatable :: ddip, trans_ddip
        double precision :: Energy, td_En
        double precision,dimension(1:100) :: energy_scan  !It's allocated (not very big)
        !Normal mode analysis
        double precision,dimension(:),allocatable :: freq
        double precision,dimension(:,:),allocatable :: L
        double precision,dimension(:),allocatable :: H
#else
        real::boxX,boxY,boxZ
        !Electronic properties (GS and TD)
        real,dimension(1:3) :: dip, trans_dip
        real,dimension(:,:),allocatable :: ddip, trans_ddip
        real :: Energy, td_En
        real,dimension(1:100) :: energy_scan  !It's allocated (not very big)
        !Normal mode analysis
        real,dimension(:),allocatable :: freq
        real,dimension(:,:),allocatable :: L
        real,dimension(:),allocatable :: H
#endif
        integer :: Nnm
    end type str_system

    !============================================
    ! STR_RESIDUE
    !============================================
    type str_residue
        character::name*5=""
        integer::natoms
        type(str_atom),dimension(:),allocatable::atom
        type(str_geom) :: geom
        integer::frst_atom,lst_atom
    end type str_residue



end module structure_types
