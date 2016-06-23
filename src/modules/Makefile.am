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
libmoleculartools_a_SOURCES =         \
    alerts.f90                        \
    line_preprocess.f90               \
    matrix.f90                        \
    matrix_print.f90                  \
    constants.f90                     \
    verbosity.f90                     \
\
    io.f90                            \
    version.f90                       \
    structure_types.f90               \
\
    metrics.f90                       \
    atomic_geom.f90                   \
    molecular_structure.f90           \
\
    symmetry.f90                      \
\
    manage_fcc.f90                    \
    manage_gamess.f90                 \
    manage_gaussian.f90               \
    manage_gmx.f90                    \
    manage_molcas.f90                 \
    manage_molden.f90                 \
    manage_molpro.f90                 \
    manage_turbomol.f90               \
    manage_pdb.f90                    \
    manage_psi4.f90                   \
\
    manage_pdb.f90                    \
    manage_xyz.f90                    \
    manage_xyz_molec.f90              \
    manage_gro_molec.f90              \
    zmat_manage.f90                   \
\
    generic_io.f90                    \
    generic_io_molec.f90              \
\
    ff_build.f90                      \
    internal.f90                      \
    thermochemistry.f90               \
    vibrational_analysis.f90          \
    vertical_model.f90

# Generate version file only if needed. To that end the
# dependency of version.f90 is a non-file-associated 
# instruction that only changes the file if requiered
# This avoids recompiling everything every time
version.f90: get_git_version
get_git_version:
	@bash ./get_git_version.sh $(FCFLAGS) $(LIBS) || echo "Reusing version file"


CLEANFILES =*.mod 

 