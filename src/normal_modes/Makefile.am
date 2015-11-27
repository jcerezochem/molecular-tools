AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = nm_internal_new

#nm_cartesian freq2fc prj_nm_internal internal_numder_analysis vibrations_cart

# ============
# nm_internal
# ============
nm_internal_new_SOURCES = ../modules/alerts.f90                 \
                             ../modules/line_preprocess.f90     \
                             ../modules/constants.f90           \
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
                             ../modules/manage_gro.f90          \
                             ../modules/zmat_manage.f90         \
                             ../modules/generic_io.f90          \
                             ../modules/generic_io_molec.f90    \
                             ../modules/ff_build.f90            \
                             ../modules/atomic_geom.f90         \
                             ../modules/symmetry.f90            \
                             ../modules/internal.f90            \
                             ../modules/vibrational_analysis.f90\
                             normal_modes_internal.f90
nm_internal_new_LDFLAGS = @LIBS@


# # ============
# # nm_cartesian
# # ============
# nm_cartesian_SOURCES =../modules/alerts.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       ../modules/MatrixMod.f90\
#                       ../modules/matrix_print.f90\
#                       ../modules/vibrational_analysis.f90\
#                       ../modules/structure_types_v4.f90\
#                       ../modules/gro_manage_v4.f90\
#                       ../modules/pdb_manage_v4.f90\
#                       ../modules/gaussian_manage_v4.f90\
#                       ../modules/gaussian_manage_vNoTypes.f90\
#                       ../modules/gaussian_fchk_manage_v4.f90\
#                       ../modules/fcc_manage.f90\
#                       ../modules/xyz_manage_v4.f90\
#                       ../modules/molecular_structure.f90\
#                       ../modules/ff_build_module_v4.f90\
#                       ../modules/atomic_geom_v2.f90\
#                       normal_modes_Cartesian.f90
# nm_cartesian_LDFLAGS = @LIBS@
# 
# # ============
# # vibratios_cart
# # ============
# vibrations_cart_SOURCES =../modules/alerts.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       ../modules/MatrixMod.f90\
#                       ../modules/matrix_print.f90\
#                       ../modules/vibrational_analysis.f90\
#                       ../modules/structure_types_v4.f90\
#                       ../modules/gro_manage_v4.f90\
#                       ../modules/pdb_manage_v4.f90\
#                       ../modules/gaussian_manage_v4.f90\
#                       ../modules/gaussian_manage_vNoTypes.f90\
#                       ../modules/gaussian_manage_lowlevel.f90\
#                       ../modules/gaussian_fchk_manage_v4.f90\
#                       ../modules/fcc_manage.f90\
#                       ../modules/xyz_manage_v4.f90\
#                       ../modules/molecular_structure.f90\
#                       ../modules/ff_build_module_v4.f90\
#                       ../modules/atomic_geom_v2.f90\
#                       vibrations_cart.f90
# vibrations_cart_LDFLAGS = @LIBS@
# 
# # ========
# # freq2fc
# # ========
# freq2fc_SOURCES =../modules/alerts.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       freq2fc.f90
# freq2fc_LDFLAGS = @LIBS@
# 
# # =========================
# # internal_numder_analysis
# # =========================
# # internal_numder_analysis_v4.0.1.f90
# internal_numder_analysis_SOURCES =../modules/alerts.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       internal_numder_analysis_v4.0.1.f90
# internal_numder_analysis_LDFLAGS = @LIBS@
# 
# 
# # ============
# # prj_nm_internal
# # ============
# prj_nm_internal_SOURCES = ../modules/alerts.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       ../modules/MatrixMod.f90\
#                       ../modules/matrix_print.f90\
#                       ../modules/structure_types_v4.f90\
#                       ../modules/gro_manage_v4.f90\
#                       ../modules/pdb_manage_v4.f90\
#                       ../modules/gaussian_manage_v4.f90\
#                       ../modules/gaussian_fchk_manage_v4.f90\
#                       ../modules/xyz_manage_v4.f90\
#                       ../modules/molcas_unsym_manage_v4.f90\
#                       ../modules/molecular_structure.f90\
#                       ../modules/ff_build_module_v4.f90\
#                       ../modules/atomic_geom_v2.f90\
#                       ../modules/symmetry_mod_v4.f90\
#                       ../modules/internal_SR9_v4_NvibAdd.f90\
#                       ../modules/zmat_manage_v4.f90\
#                       ../modules/psi4_manage.f90\
#                       project_normal_modes_internal.f90
# prj_nm_internal_LDFLAGS = @LIBS@


CLEANFILES =*.mod 

