AUTOMAKE_OPTIONS = subdir-objects

LDFLAGS = -L$(top_builddir)/libs/AnalyticalBders @LIBS@ 

bin_PROGRAMS = nm_internal_new numder_fc numder_L prj_nm_internal_new nm_cartesian_new gen_Lder_numx\
               check_derivatives numder_Lint

# freq2fc   vibrations_cart

# ============
# nm_internal
# ============
nm_internal_new_SOURCES =    ../modules/version.f90             \
                             ../modules/alerts.f90              \
                             ../modules/line_preprocess.f90     \
                             ../modules/constants.f90           \
                             ../modules/verbosity.f90           \
                             ../modules/matrix.f90              \
                             ../modules/matrix_print.f90        \
                             ../modules/structure_types.f90     \
                             ../modules/molecular_structure.f90 \
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90         \
                             ../modules/symmetry.f90            \
                             ../modules/manage_fcc.f90          \
                             ../modules/manage_gamess.f90       \
                             ../modules/manage_gaussian.f90     \
                             ../modules/manage_gmx.f90          \
                             ../modules/manage_molcas.f90       \
                             ../modules/manage_molden.f90       \
                             ../modules/manage_molpro.f90       \
                             ../modules/manage_pdb.f90          \
                             ../modules/manage_psi4.f90         \
                             ../modules/manage_xyz.f90          \
                             ../modules/manage_xyz_molec.f90    \
                             ../modules/manage_gro.f90          \
                             ../modules/zmat_manage.f90         \
                             ../modules/generic_io.f90          \
                             ../modules/generic_io_molec.f90    \
                             ../modules/ff_build.f90            \
                             ../modules/internal.f90            \
                             ../modules/vibrational_analysis.f90\
                             normal_modes_internal.f90

# ============
# nm_cartesian
# ============
nm_cartesian_new_SOURCES =   ../modules/version.f90             \
                             ../modules/alerts.f90              \
                             ../modules/line_preprocess.f90     \
                             ../modules/constants.f90           \
                             ../modules/verbosity.f90           \
                             ../modules/matrix.f90              \
                             ../modules/matrix_print.f90        \
                             ../modules/structure_types.f90     \
                             ../modules/molecular_structure.f90 \
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90         \
                             ../modules/manage_fcc.f90          \
                             ../modules/manage_gamess.f90       \
                             ../modules/manage_gaussian.f90     \
                             ../modules/manage_gmx.f90          \
                             ../modules/manage_molcas.f90       \
                             ../modules/manage_molden.f90       \
                             ../modules/manage_molpro.f90       \
                             ../modules/manage_pdb.f90          \
                             ../modules/manage_psi4.f90         \
                             ../modules/manage_xyz.f90          \
                             ../modules/manage_xyz_molec.f90    \
                             ../modules/manage_gro.f90          \
                             ../modules/zmat_manage.f90         \
                             ../modules/generic_io.f90          \
                             ../modules/generic_io_molec.f90    \
                             ../modules/symmetry.f90            \
                             ../modules/ff_build.f90            \
                             ../modules/internal.f90            \
                             ../modules/vibrational_analysis.f90\
                             normal_modes_Cartesian.f90
# ============
# gen_Lder_numx
# ============
gen_Lder_numx_SOURCES =../modules/alerts.f90                 \
                             ../modules/line_preprocess.f90     \
                             ../modules/constants.f90           \
                             ../modules/verbosity.f90           \
                             ../modules/matrix.f90              \
                             ../modules/matrix_print.f90        \
                             ../modules/structure_types.f90     \
                             ../modules/molecular_structure.f90 \
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90         \
                             ../modules/manage_fcc.f90          \
                             ../modules/manage_gamess.f90       \
                             ../modules/manage_gaussian.f90     \
                             ../modules/manage_gmx.f90          \
                             ../modules/manage_molcas.f90       \
                             ../modules/manage_molden.f90       \
                             ../modules/manage_molpro.f90       \
                             ../modules/manage_pdb.f90          \
                             ../modules/manage_psi4.f90         \
                             ../modules/manage_xyz.f90          \
                             ../modules/manage_xyz_molec.f90    \
                             ../modules/manage_gro.f90          \
                             ../modules/generic_io.f90          \
                             ../modules/generic_io_molec.f90    \
                             gen_Lder_numx.f90

# ============
# check_derivatives
# ============
check_derivatives_SOURCES = ../modules/alerts.f90                 \
                             ../modules/line_preprocess.f90     \
                             ../modules/constants.f90           \
                             ../modules/verbosity.f90           \
                             ../modules/matrix.f90              \
                             ../modules/matrix_print.f90        \
                             ../modules/structure_types.f90     \
                             ../modules/molecular_structure.f90 \
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90         \
                             ../modules/symmetry.f90            \
                             ../modules/manage_fcc.f90          \
                             ../modules/manage_gamess.f90       \
                             ../modules/manage_gaussian.f90     \
                             ../modules/manage_gmx.f90          \
                             ../modules/manage_molcas.f90       \
                             ../modules/manage_molden.f90       \
                             ../modules/manage_molpro.f90       \
                             ../modules/manage_pdb.f90          \
                             ../modules/manage_psi4.f90         \
                             ../modules/manage_xyz.f90          \
                             ../modules/manage_xyz_molec.f90    \
                             ../modules/manage_gro.f90          \
                             ../modules/zmat_manage.f90         \
                             ../modules/generic_io.f90          \
                             ../modules/generic_io_molec.f90    \
                             ../modules/ff_build.f90            \
                             ../modules/internal.f90            \
                             ../modules/vibrational_analysis.f90\
                             check_derivatives.f90

# 
# # ========
# # freq2fc
# # ========
# freq2fc_SOURCES =../modules/alerts.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       freq2fc.f90
# freq2fc_LDFLAGS = @LIBS@

# =========================
# numder_fc
# =========================
# numder_fc.f90
numder_fc_SOURCES =../modules/alerts.f90\
                   ../modules/line_preprocess.f90\
                   ../modules/constants.f90\
                   ../modules/manage_gaussian.f90     \
                   numder_fc.f90


# =========================
# numder_L
# =========================
# numder_fc.f90
numder_L_SOURCES = ../modules/alerts.f90\
                   ../modules/line_preprocess.f90\
                   ../modules/constants.f90\
                   ../modules/verbosity.f90           \
                   ../modules/matrix.f90              \
                   ../modules/matrix_print.f90        \
                   ../modules/structure_types.f90     \
                   ../modules/molecular_structure.f90 \
                   ../modules/manage_fcc.f90          \
                   ../modules/manage_gamess.f90       \
                   ../modules/manage_gaussian.f90     \
                   ../modules/manage_gmx.f90          \
                   ../modules/manage_molcas.f90       \
                   ../modules/manage_molden.f90       \
                   ../modules/manage_molpro.f90       \
                   ../modules/manage_pdb.f90          \
                   ../modules/manage_psi4.f90         \
                   ../modules/manage_xyz.f90          \
                   ../modules/manage_xyz_molec.f90    \
                   ../modules/manage_gro.f90          \
                   ../modules/generic_io.f90          \
                   ../modules/generic_io_molec.f90    \
                   ../modules/atomic_geom.f90         \
                   ../modules/symmetry.f90            \
                   ../modules/vibrational_analysis.f90\
                   numder_L.f90

# ============
# numder_Lint
# ============
numder_Lint_SOURCES = ../modules/alerts.f90                 \
                             ../modules/line_preprocess.f90     \
                             ../modules/constants.f90           \
                             ../modules/verbosity.f90           \
                             ../modules/matrix.f90              \
                             ../modules/matrix_print.f90        \
                             ../modules/structure_types.f90     \
                             ../modules/molecular_structure.f90 \
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90         \
                             ../modules/symmetry.f90            \
                             ../modules/manage_fcc.f90          \
                             ../modules/manage_gamess.f90       \
                             ../modules/manage_gaussian.f90     \
                             ../modules/manage_gmx.f90          \
                             ../modules/manage_molcas.f90       \
                             ../modules/manage_molden.f90       \
                             ../modules/manage_molpro.f90       \
                             ../modules/manage_pdb.f90          \
                             ../modules/manage_psi4.f90         \
                             ../modules/manage_xyz.f90          \
                             ../modules/manage_xyz_molec.f90    \
                             ../modules/manage_gro.f90          \
                             ../modules/zmat_manage.f90         \
                             ../modules/generic_io.f90          \
                             ../modules/generic_io_molec.f90    \
                             ../modules/ff_build.f90            \
                             ../modules/internal.f90            \
                             ../modules/vibrational_analysis.f90\
                              numder_Lint.f90


# ============
# prj_nm_internal
# ============
prj_nm_internal_new_SOURCES =../modules/alerts.f90              \
                             ../modules/line_preprocess.f90     \
                             ../modules/constants.f90           \
                             ../modules/verbosity.f90           \
                             ../modules/matrix.f90              \
                             ../modules/matrix_print.f90        \
                             ../modules/structure_types.f90     \
                             ../modules/molecular_structure.f90 \
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90         \
                             ../modules/manage_fcc.f90          \
                             ../modules/manage_gamess.f90       \
                             ../modules/manage_gaussian.f90     \
                             ../modules/manage_gmx.f90          \
                             ../modules/manage_molcas.f90       \
                             ../modules/manage_molden.f90       \
                             ../modules/manage_molpro.f90       \
                             ../modules/manage_pdb.f90          \
                             ../modules/manage_psi4.f90         \
                             ../modules/manage_xyz.f90          \
                             ../modules/manage_xyz_molec.f90    \
                             ../modules/manage_gro.f90          \
                             ../modules/zmat_manage.f90         \
                             ../modules/generic_io.f90          \
                             ../modules/generic_io_molec.f90    \
                             ../modules/ff_build.f90            \
                             ../modules/symmetry.f90            \
                             ../modules/internal.f90            \
                             ../modules/vibrational_analysis.f90\
                             project_normal_modes_internal.f90


CLEANFILES =*.mod 

