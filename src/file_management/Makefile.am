AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = reform_fchk rmsd_fit fchk2gro read_irc
#get_gscan fchk2gro get_symatoms reform_fchk numderX

# ============
# reform_fchk
# ============
reform_fchk_SOURCES = ../modules/version.f90                 \
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
                      ../modules/internal.f90                \
                      ../modules/vibrational_analysis.f90    \
                      reform_fchk.f90
reform_fchk_LDFLAGS = @LIBS@

rmsd_fit_SOURCES = ../modules/version.f90                 \
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
                   ../modules/internal.f90                \
                   rmsd_fit.f90
rmsd_fit_LDFLAGS = @LIBS@

fchk2gro_SOURCES = ../modules/version.f90                 \
                   ../modules/alerts.f90                  \
                   ../modules/line_preprocess.f90         \
                   ../modules/constants.f90               \
                   ../modules/verbosity.f90               \
                   ../modules/matrix.f90                  \
                   ../modules/matrix_print.f90            \
                   ../modules/structure_types.f90         \
                   ../modules/molecular_structure.f90     \
                   ../modules/symmetry.f90                \
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
                   ../modules/internal.f90                \
                   fchk2gro.f90
fchk2gro_LDFLAGS = @LIBS@

read_irc_SOURCES = ../modules/version.f90                 \
                   ../modules/alerts.f90                  \
                   ../modules/line_preprocess.f90         \
                   ../modules/constants.f90               \
                   ../modules/verbosity.f90               \
                   ../modules/matrix.f90                  \
                   ../modules/matrix_print.f90            \
                   ../modules/structure_types.f90         \
                   ../modules/molecular_structure.f90     \
                   ../modules/symmetry.f90                \
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
                   ../modules/internal.f90                \
                   read_irc.f90
read_irc_LDFLAGS = @LIBS@

# 
# numderX_SOURCES = ../modules/alerts.f90\
#                       ../modules/MatrixMod.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       ../modules/structure_types_v4.f90\
#                       ../modules/gro_manage_v4.f90\
#                       ../modules/pdb_manage_v4.f90\
#                       ../modules/gaussian_manage_v4.f90\
#                       ../modules/gaussian_fchk_manage_v4.f90\
#                       ../modules/xyz_manage_v4.f90\
#                       ../modules/molden_manage.f90\
#                       ../modules/molcas_unsym_manage_v4.f90\
#                       ../modules/psi4_manage.f90\
#                       ../modules/gamess_manage.f90\
#                       ../modules/molecular_structure.f90\
#                       ../modules/ff_build_module_v4.f90\
#                       numderX.f90
# numderX_LDFLAGS = @LIBS@
# 
# get_gscan_SOURCES = ../modules/alerts.f90\
#                       ../modules/MatrixMod.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       ../modules/structure_types_v4.f90\
#                       ../modules/gaussian_manage_v4.f90\
#                       ../modules/gaussian_fchk_manage_v4.f90\
#                       ../modules/xyz_manage_v4.f90\
#                       ../modules/molden_manage.f90\
#                       ../modules/molecular_structure.f90\
#                       ../modules/ff_build_module_v4.f90\
#                       ../modules/atomic_geom_v2.f90\
#                       ../modules/symmetry_mod_v4.f90\
#                       get_gscan.f90
# get_gscan_LDFLAGS = @LIBS@
# 
# get_symatoms_SOURCES = ../modules/alerts.f90\
#                       ../modules/MatrixMod.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       ../modules/structure_types_v4.f90\
#                       ../modules/gro_manage_v4.f90\
#                       ../modules/pdb_manage_v4.f90\
#                       ../modules/gaussian_manage_v4.f90\
#                       ../modules/gaussian_fchk_manage_v4.f90\
#                       ../modules/xyz_manage_v4.f90\
#                       ../modules/molden_manage.f90\
#                       ../modules/molcas_unsym_manage_v4.f90\
#                       ../modules/molecular_structure.f90\
#                       ../modules/ff_build_module_v4.f90\
#                       ../modules/symmetry_mod_v4.f90\
#                       get_symatoms.f90
# get_symatoms_LDFLAGS = @LIBS@
# 

CLEANFILES =*.mod 

