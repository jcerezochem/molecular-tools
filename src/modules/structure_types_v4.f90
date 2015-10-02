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
    !            versions (this is NOT allocatable)
    ! (13/06/12) Added fftype (character(len=6)) to atom. This stores the atom
    !            type in the force field, while attype stores the PDB attype
    ! (14/11/12) Added PG (point group) to system
    ! (10/01/13) Version3: str_residue and str_system to be merged. Residues will be  
    !             managed as str_system (TODO: recursively, is possible??). This is 
    !             transition version that make str_system fully featured (incl. geom) 
    !   NOTE: this version is under development!  -- see V2 to add here the same improvements (inititalisers)
    ! (11/02/14) Version4: New simplified types for:
    !                * atom -> resmol -x-> XsystX (system type is deleted)
    !            and additional types for:
    !                * bonded -> resmol // Not standalone type since it is useful 
    !                               to have residue and bonded to read topology
    !                * molprops  -x-> // Now is a standalone type
    ! (14/02/14) TESTING: new str_type: str_internal to store internal coordinates
    ! (27/02/14) Added str_job and included molec%job and molec%props(allocatable)
    ! (20/06/14) Version3.2 (in internal project): add impropers to geom
    !            Version4: as 4 but with added str_atom_light
    !---------------------------------------------------------

    !SIZES
    integer,parameter:: MAX_ATOMS = 100000, &
                        MAX_ATM_RES = 36000,  &
                        MAX_CONNEXIONS = 6, &
                        MAX_DERIVED_CNX = 5000

    !DERIVED TYPES
!   Auxiliar new types (need to be defined first)
!--------------------------
    type str_bonded
!--------------------------
        integer,dimension(1:MAX_DERIVED_CNX,2) :: bond
        integer,dimension(1:MAX_DERIVED_CNX,2) :: pair
        integer,dimension(1:MAX_DERIVED_CNX,3) :: angle
        integer,dimension(1:MAX_DERIVED_CNX,4) :: dihed
        integer,dimension(1:MAX_DERIVED_CNX,4) :: improp
        integer :: npairs, nangles, ndihed, nbonds, nimprop
    end type str_bonded

!--------------------------
    type str_molprops
!--------------------------
        character::PG*5
        integer :: Nnm
#ifdef DOUBLE
        !Electronic properties (GS and TD)
        double precision,dimension(1:3) :: dip, trans_dip
        double precision,dimension(1:3*MAX_ATOMS,1:3) :: ddip, trans_ddip
        double precision :: Energy, td_En
        double precision,dimension(1:100) :: energy_scan  !It's allocated (not very big)
        !Normal mode analysis
        double precision,dimension(1:1000) :: freq,grad
        double precision,dimension(1:1000,1:1000) :: L
        double precision,dimension(1:100000) :: H
#else
        !Electronic properties (GS and TD)
        real,dimension(1:3) :: dip, trans_dip
        real,dimension(1:3*MAX_ATOMS,1:3) :: ddip, trans_ddip
        real :: Energy, td_En
        real,dimension(1:100) :: energy_scan  !It's allocated (not very big)
        !Normal mode analysis
        real,dimension(1:1000) :: freq,grad
        real,dimension(1:1000,1:1000) :: L
        real,dimension(1:100000) :: H
#endif
    end type str_molprops

!--------------------------
    type str_job
!--------------------------
        character(len=100) :: title=""
        character(len=20)  :: type="XX"
        character(len=20)  :: method=""
        character(len=20)  :: basis=""
    end type str_job

!===========================
    type str_atom
!===========================
        character(len=6)::attype="ATOM  "
        character::fftype*6
        character::name*5
        character::resname*5="UNK"
        character(len=1)::chain=" ",alter_loc=" ", ins_code=" "
        character(len=2) :: element=""
        integer  ::resseq=1, AtNum
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
        integer,dimension(1:MAX_CONNEXIONS) :: connect
        integer :: nbonds
        character(len=6),dimension(1:MAX_CONNEXIONS+1) :: env 
    end type str_atom

!===========================
    type str_atom_light
!===========================
        character(len=6)::attype="ATOM  "
        character::fftype*6
        character::name*5
        character::resname*5="UNK"
        character(len=1)::chain=" ",alter_loc=" ", ins_code=" "
        character(len=2) :: element=""
        integer  ::resseq=1, AtNum
#ifdef DOUBLE
        !-------------------------------------------------------------
        double precision:: x, y, z
        !-------------------------------------------------------------
        double precision:: mass, q
#else
        !-------------------------------------------------------------
        real:: x, y, z 
        !-------------------------------------------------------------
        real:: mass, q
#endif
    end type str_atom_light

!===========================
    type str_resmol
!===========================
        ! This type can represent:
        ! * A residue (within a multiple residue system)
        ! * A molecule that can contain residues
        character::PG*5="XX" !Default for unknown
        character::name*5="UNK"
        integer::natoms
        type(str_atom),dimension(1:MAX_ATM_RES)::atom
        integer::frst_atom,lst_atom
        !If connected to other residues
!         integer,dimension(1:MAX_CONNEXIONS) :: connect
        integer::nbonds
        !If it contains several residues
        integer:: nres
        !Other info that can be read from structure files
        character(len=100) :: title=""
#ifdef DOUBLE
        double precision::boxX,boxY,boxZ
#else
        real::boxX,boxY,boxZ
#endif
!         !Include job info (allocatable, so it need to be "activated" before use)
!         type(str_molprops),allocatable :: props
        !Include job info (should be allocatable?) 
        type(str_job) :: job
        !Include bonded info (useful for force field management) (should be allocatable?)
        type(str_bonded) :: geom
        !Including com and cog (c. of geometry)
#ifdef DOUBLE
        double precision::cogX,cogY,cogZ, &
                          comX,comY,comZ
#else
        real::cogX,cogY,cogZ, & 
              comX,comY,comZ
#endif
!===================================
        !Make it recursive
!         type(str_resmol),dimension(:),allocatable :: residue
!.........................................................................
!       ERROR: gfortran complaints with:
!               type(str_resmol),dimension(:),allocatable :: residue
!                                                                   1
!               Error: El componente en (1) debe tener el atributo POINTER
!.........................................................................
!       This feature is not vital, though. The same behaviour is obtained 
!       using an array of resmol types

    end type str_resmol

!===========================
    type str_resmol_light
!===========================
        ! This type can represent:
        ! * A residue (within a multiple residue system)
        ! * A molecule that can contain residues
        character::PG*5="XX" !Default for unknown
        character::name*5="UNK"
        integer::natoms
        type(str_atom_light),dimension(1:MAX_ATM_RES)::atom
        integer::frst_atom,lst_atom
        !If connected to other residues
!         integer,dimension(1:MAX_CONNEXIONS) :: connect
        integer::nbonds
        !If it contains several residues
        integer:: nres
        !Other info that can be read from structure files
        character(len=100) :: title=""
#ifdef DOUBLE
        double precision::boxX,boxY,boxZ
#else
        real::boxX,boxY,boxZ
#endif
!         !Include job info (allocatable, so it need to be "activated" before use)
!         type(str_molprops),allocatable :: props
        !Include job info (should be allocatable?) 
        type(str_job) :: job
        !Include bonded info (useful for force field management) (should be allocatable?)
        type(str_bonded) :: geom
        !Including com and cog (c. of geometry)
#ifdef DOUBLE
        double precision::cogX,cogY,cogZ, &
                          comX,comY,comZ
#else
        real::cogX,cogY,cogZ, & 
              comX,comY,comZ
#endif
!===================================
        !Make it recursive
!         type(str_resmol),dimension(:),allocatable :: residue
!.........................................................................
!       ERROR: gfortran complaints with:
!               type(str_resmol),dimension(:),allocatable :: residue
!                                                                   1
!               Error: El componente en (1) debe tener el atributo POINTER
!.........................................................................
!       This feature is not vital, though. The same behaviour is obtained 
!       using an array of resmol types

    end type str_resmol_light


!   NEW INTERNAL STRUCTURE TYPE
!===========================
    type str_internal
!===========================
        character(len=1) :: coord="" ! identify type of internal: "b"(bond), "a"(angle), "d"(dihedral), "i"(improper), "x"(linear comb)
        integer,dimension(1:4) :: atoms !atom index involved
#ifdef DOUBLE
        double precision::value
#else
        real:: value
#endif
        character(len=3) :: units=""
    end type str_internal



end module structure_types
