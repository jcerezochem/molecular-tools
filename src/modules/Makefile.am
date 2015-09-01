AUTOMAKE_OPTIONS = subdir-objects

# Build a convenience library
# (but this migh also work as a real lib...)
# only to have all objects AND modules
# in this folder (otherwise, the modules
# would be by default on the folder of the
# first binary to call them.
# In order for the compilers to find the 
# module files, use I/path/to/here
noinst_LIBRARIES = libmoleculartools.a

# We compile only the modules that are
# used by the binaries out there
libmoleculartools_a_SOURCES =        \
    alerts.f90                        \
    line_preprocess.f90               \
    MatrixMod.f90                     \
    constants_mod.f90                 \
\
    symmetry_mod_notypes_v2.f90       \
    gaussian_manage_vNoTypes.f90      \
\
    structure_types_v4.f90            \
\
    atomic_geom_v2.f90                \
    fcc_manage.f90                    \
    ff_build_module_v4.f90            \
    gamess_manage.f90                 \
    gaussian_fchk_manage_v4.f90       \
    gaussian_manage_v4.f90            \
    gaussian_manage_lowlevel.f90      \
    gro_manage_v4.f90                 \
    internal_SR9_v4.f90               \
    molcas_unsym_manage_v4.f90        \
    molden_manage.f90                 \
    molecular_structure.f90           \
    molpro_manage_v4.f90              \
    pdb_manage_v4.f90                 \
    psi4_manage.f90                   \
    symmetry_mod_v4.f90               \
    xyz_manage_v4.f90                 \
    zmat_manage_v4.f90                

#    sizes.f90                         \
#    structure_types_v5.f90            

 