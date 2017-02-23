module structure_types

    !==============================================================
    ! This code is part of MOLECULAR_TOOLS (version 0.4/February 2014)
    !==============================================================
    !Description
    ! Derived types handling molecular information. Only structural 
    ! information is included (NOT energy, gradients...), otherwise
    ! the types become too comples and not very usable.
    !
    ! NOTES
    ! -----
    ! Limits:
    ! MaxAtoms/system = 100 000 (maximum since counter is I5)
    ! MaxAtoms/residue = 200 (is safe)
    ! Shoul work for GRO and PDB formats even if
    ! character lengths are different (the larger
    ! is used). If needed, use ADJUSTL(string).
    !
    ! TODO
    !  - Remove job(str_job) and geom(str_bonded) from str_resmol.
    !    * Use job as standalone. Feed to required applications (e.g. 
    !      write_gcom as optional)
    !    * geom: replace by ic (internals) and use it as standalone
    !
    !  - Dynamic memory allocation
    !---------------------------------------------------------

    !SIZES
    integer,parameter:: MAX_ATOMS = 10000, &
                        MAX_ATM_RES = 6000,  &
                        MAX_CONNEXIONS = 6, &
                        MAX_DERIVED_CNX = 1000

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
        integer :: npairs=0, nangles=0, ndihed=0, nbonds=0, nimprop=0
    end type str_bonded

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
        integer  ::resseq=1, AtNum=0
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
    type str_resmol
!===========================
        ! This type can represent:
        ! * A residue (within a multiple residue system)
        ! * A molecule that can contain residues
        character(len=4)::units
        character::PG*5="XX" !Default for unknown
        character::name*5="UNK"
        integer::natoms
        type(str_atom),dimension(:),allocatable::atom
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


    contains

    subroutine allocate_atoms(molecule,natoms)
        
        type(str_resmol),intent(inout) :: molecule
        integer,intent(in),optional    :: natoms

        ! local
        integer :: n

        if (present(natoms)) then
            n = natoms
        else
            n = MAX_ATM_RES
        endif

        allocate(molecule%atom(1:n))

        return
        
    end subroutine allocate_atoms

end module structure_types
