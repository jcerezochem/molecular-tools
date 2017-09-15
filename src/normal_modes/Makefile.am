AUTOMAKE_OPTIONS = subdir-objects

LDFLAGS = -L$(top_builddir)/libs/AnalyticalBders @LIBS@ 

bin_PROGRAMS = nm_internal numder_fc prj_nm_internal nm_cartesian \
               prj_grad_on_nm_int nm_internal_prjCrt nm_internal_prjInt \
               cartesian2internal_forces
# prj_grad_on_nm_cart 
# Deprecated:
# numder_L numder_Lint gen_Lder_numx check_derivatives

# freq2fc   vibrations_cart

# ============
# cartesian2internal_forces
# ============
cartesian2internal_forces_SOURCES =   \
                        ../modules/verbosity.f90\
                        ../modules/alerts.f90\
                        ../modules/line_preprocess.f90\
                        ../modules/constants.f90\
                        ../modules/io.f90                  \
                        ../modules/version.f90             \
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
                        ../modules/manage_cfour.f90        \
                        ../modules/manage_gmx.f90          \
                        ../modules/manage_molcas.f90       \
                        ../modules/manage_molden.f90       \
                        ../modules/manage_molpro.f90       \
                        ../modules/manage_turbomol.f90     \
                        ../modules/manage_pdb.f90          \
                        ../modules/manage_psi4.f90         \
                        ../modules/manage_xyz.f90          \
                        ../modules/manage_xyz_molec.f90    \
                        ../modules/manage_gro_molec.f90          \
                        ../modules/zmat_manage.f90         \
                        ../modules/generic_io.f90          \
                        ../modules/generic_io_molec.f90    \
                        ../modules/internal.f90            \
                        ../modules/thermochemistry.f90     \
                        ../modules/vibrational_analysis.f90\
                        ../modules/vertical_model.f90\
                        cartesian2internal_forces.f90

# ============
# nm_internal
# ============
nm_internal_SOURCES =   \
                        ../modules/verbosity.f90\
                        ../modules/alerts.f90\
                        ../modules/line_preprocess.f90\
                        ../modules/constants.f90\
                        ../modules/io.f90                  \
                        ../modules/version.f90             \
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
                        ../modules/manage_cfour.f90        \
                        ../modules/manage_gmx.f90          \
                        ../modules/manage_molcas.f90       \
                        ../modules/manage_molden.f90       \
                        ../modules/manage_molpro.f90       \
                        ../modules/manage_turbomol.f90     \
                        ../modules/manage_pdb.f90          \
                        ../modules/manage_psi4.f90         \
                        ../modules/manage_xyz.f90          \
                        ../modules/manage_xyz_molec.f90    \
                        ../modules/manage_gro_molec.f90          \
                        ../modules/zmat_manage.f90         \
                        ../modules/generic_io.f90          \
                        ../modules/generic_io_molec.f90    \
                        ../modules/internal.f90            \
                        ../modules/thermochemistry.f90     \
                        ../modules/vibrational_analysis.f90\
                        ../modules/vertical_model.f90\
                        normal_modes_internal.f90
# ============
# nm_internal
# ============
nm_internal_prjInt_SOURCES =   \
                        ../modules/verbosity.f90\
                        ../modules/alerts.f90\
                        ../modules/line_preprocess.f90\
                        ../modules/constants.f90\
                        ../modules/io.f90                  \
                        ../modules/version.f90             \
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
                        ../modules/manage_cfour.f90        \
                        ../modules/manage_gmx.f90          \
                        ../modules/manage_molcas.f90       \
                        ../modules/manage_molden.f90       \
                        ../modules/manage_molpro.f90       \
                        ../modules/manage_turbomol.f90     \
                        ../modules/manage_pdb.f90          \
                        ../modules/manage_psi4.f90         \
                        ../modules/manage_xyz.f90          \
                        ../modules/manage_xyz_molec.f90    \
                        ../modules/manage_gro_molec.f90          \
                        ../modules/zmat_manage.f90         \
                        ../modules/generic_io.f90          \
                        ../modules/generic_io_molec.f90    \
                        ../modules/internal.f90            \
                        ../modules/thermochemistry.f90     \
                        ../modules/vibrational_analysis.f90\
                        ../modules/vertical_model.f90\
                        normal_modes_internal_prjInt.f90
# ============
# nm_internal
# ============
nm_internal_prjCrt_SOURCES =   \
                        ../modules/verbosity.f90\
                        ../modules/alerts.f90\
                        ../modules/line_preprocess.f90\
                        ../modules/constants.f90\
                        ../modules/io.f90                  \
                        ../modules/version.f90             \
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
                        ../modules/manage_cfour.f90        \
                        ../modules/manage_gmx.f90          \
                        ../modules/manage_molcas.f90       \
                        ../modules/manage_molden.f90       \
                        ../modules/manage_molpro.f90       \
                        ../modules/manage_turbomol.f90     \
                        ../modules/manage_pdb.f90          \
                        ../modules/manage_psi4.f90         \
                        ../modules/manage_xyz.f90          \
                        ../modules/manage_xyz_molec.f90    \
                        ../modules/manage_gro_molec.f90          \
                        ../modules/zmat_manage.f90         \
                        ../modules/generic_io.f90          \
                        ../modules/generic_io_molec.f90    \
                        ../modules/internal.f90            \
                        ../modules/thermochemistry.f90     \
                        ../modules/vibrational_analysis.f90\
                        ../modules/vertical_model.f90\
                        normal_modes_internal_prjCrt.f90

# ============
# nm_cartesian
# ============
nm_cartesian_SOURCES =  \
                        ../modules/verbosity.f90\
                        ../modules/alerts.f90\
                        ../modules/line_preprocess.f90\
                        ../modules/constants.f90\
                        ../modules/io.f90                  \
                        ../modules/version.f90             \
                        ../modules/matrix.f90              \
                        ../modules/matrix_print.f90        \
                        ../modules/structure_types.f90     \
                        ../modules/molecular_structure.f90 \
                        ../modules/metrics.f90             \
                        ../modules/atomic_geom.f90         \
                        ../modules/manage_fcc.f90          \
                        ../modules/manage_gamess.f90       \
                        ../modules/manage_gaussian.f90     \
                        ../modules/manage_cfour.f90        \
                        ../modules/manage_gmx.f90          \
                        ../modules/manage_molcas.f90       \
                        ../modules/manage_molden.f90       \
                        ../modules/manage_molpro.f90       \
                        ../modules/manage_turbomol.f90     \
                        ../modules/manage_pdb.f90          \
                        ../modules/manage_psi4.f90         \
                        ../modules/manage_xyz.f90          \
                        ../modules/manage_xyz_molec.f90    \
                        ../modules/manage_gro_molec.f90          \
                        ../modules/zmat_manage.f90         \
                        ../modules/generic_io.f90          \
                        ../modules/generic_io_molec.f90    \
                        ../modules/symmetry.f90            \
                        ../modules/internal.f90            \
                        ../modules/thermochemistry.f90     \
                        ../modules/vibrational_analysis.f90\
                        ../modules/vertical_model.f90\
                        normal_modes_Cartesian.f90

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

# ============
# prj_nm_internal
# ============
prj_nm_internal_SOURCES =../modules/alerts.f90              \
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
                             ../modules/manage_cfour.f90        \
                             ../modules/manage_gmx.f90          \
                             ../modules/manage_molcas.f90       \
                             ../modules/manage_molden.f90       \
                             ../modules/manage_molpro.f90       \
                             ../modules/manage_turbomol.f90     \
                             ../modules/manage_pdb.f90          \
                             ../modules/manage_psi4.f90         \
                             ../modules/manage_xyz.f90          \
                             ../modules/manage_xyz_molec.f90    \
                             ../modules/manage_gro_molec.f90          \
                             ../modules/zmat_manage.f90         \
                             ../modules/generic_io.f90          \
                             ../modules/generic_io_molec.f90    \
                             ../modules/symmetry.f90            \
                             ../modules/internal.f90            \
                             ../modules/vibrational_analysis.f90\
                             project_normal_modes_internal.f90


prj_grad_on_nm_int_SOURCES = \
                             ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/io.f90                  \
                             ../modules/version.f90             \
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_cfour.f90                  \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_turbomol.f90     \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/manage_xyz_molec.f90              \
                             ../modules/manage_gro_molec.f90              \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             ../modules/vibrational_analysis.f90\
                             ../modules/vertical_model.f90\
                             project_gradient_on_nm_internal.f90

# prj_grad_on_nm_cart_SOURCES = ../modules/version.f90\
#                              ../modules/alerts.f90\
#                              ../modules/line_preprocess.f90\
#                              ../modules/constants.f90\
#                              ../modules/verbosity.f90\
#                              ../modules/io.f90\
#                              ../modules/matrix.f90\
#                              ../modules/matrix_print.f90\
#                              ../modules/structure_types.f90\
#                              ../modules/metrics.f90             \
#                              ../modules/atomic_geom.f90\
#                              ../modules/molecular_structure.f90\
#                              ../modules/manage_fcc.f90                    \
#                              ../modules/manage_gamess.f90                 \
#                              ../modules/manage_gaussian.f90               \
#                              ../modules/manage_gmx.f90                    \
#                              ../modules/manage_molcas.f90                 \
#                              ../modules/manage_molden.f90                 \
#                              ../modules/manage_molpro.f90                 \
#                              ../modules/manage_pdb.f90                    \
#                              ../modules/manage_psi4.f90                   \
#                              ../modules/manage_xyz.f90                    \
#                              ../modules/manage_xyz_molec.f90              \
#                              ../modules/manage_gro_molec.f90              \
#                              ../modules/zmat_manage.f90                   \
#                              ../modules/generic_io.f90\
#                              ../modules/generic_io_molec.f90\
#                              ../modules/symmetry.f90\
#                              ../modules/internal.f90\
#                              ../modules/vibrational_analysis.f90\
#                              ../modules/vertical_model.f90\
#                              project_gradient_on_nm_cartesian.f90


CLEANFILES =*.mod 

