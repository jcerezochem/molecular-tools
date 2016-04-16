AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = get_isotopologues

# ============
# get_isotopologues
# ============
get_isotopologues_SOURCES = ../modules/version.f90                 \
                      ../modules/alerts.f90                  \
                      ../modules/line_preprocess.f90         \
                      ../modules/constants.f90               \
                      ../modules/verbosity.f90               \
                      ../modules/matrix.f90                  \
                      ../modules/matrix_print.f90            \
                      ../modules/structure_types.f90         \
                      ../modules/molecular_structure.f90     \
                      ../modules/metrics.f90                 \
                      ../modules/atomic_geom.f90             \
                      ../modules/manage_fcc.f90              \
                      ../modules/manage_gamess.f90           \
                      ../modules/manage_gaussian.f90         \
                      ../modules/manage_gmx.f90              \
                      ../modules/manage_molcas.f90           \
                      ../modules/manage_molden.f90           \
                      ../modules/manage_molpro.f90           \
                      ../modules/manage_pdb.f90              \
                      ../modules/manage_psi4.f90             \
                      ../modules/manage_xyz.f90              \
                      ../modules/manage_xyz_molec.f90        \
                      ../modules/manage_gro_molec.f90        \
                      ../modules/zmat_manage.f90             \
                      ../modules/generic_io.f90              \
                      ../modules/generic_io_molec.f90        \
                      ../modules/ff_build.f90                \
                      ../modules/symmetry.f90                \
                      get_isotopologues.f90
reform_fchk_LDFLAGS = @LIBS@
