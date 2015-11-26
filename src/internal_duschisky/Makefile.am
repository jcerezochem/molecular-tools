AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = vertical2adiabatic_new 
# internal_duschinsky_new
# zmat2modredundant

vertical2adiabatic_new_SOURCES = ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/ff_build.f90\
                             ../modules/atomic_geom.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             ../modules/vibrational_analysis.f90\
                             vertical2adiabatic.f90
vertical2adiabatic_new_LDFLAGS = @LIBS@

# internal_duschinsky_new_SOURCES = ../modules/alerts.f90\
#                              ../modules/line_preprocess.f90\
#                              ../modules/constants.f90\
#                              ../modules/verbosity.f90\
#                              ../modules/matrix.f90\
#                              ../modules/matrix_print.f90\
#                              ../modules/structure_types.f90\
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
#                              ../modules/zmat_manage.f90                   \
#                              ../modules/generic_io.f90\
#                              ../modules/generic_io_molec.f90\
#                              ../modules/ff_build.f90\
#                              ../modules/atomic_geom.f90\
#                              ../modules/symmetry.f90\
#                              ../modules/internal.f90\
#                              internal_duschinsky.f90
# internal_duschinsky_new_LDFLAGS = @LIBS@

# 
# 
# zmat2modredundant_SOURCES = ../modules/alerts.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       ../modules/matrix.f90\
#                       ../modules/MatrixMod.f90\
#                       ../modules/matrix_print.f90\
#                       ../modules/structure_types_v4.f90\
#                       ../modules/gro_v4.f90\
#                       ../modules/pdb_v4.f90\
#                       ../modules/gaussian_v4.f90\
#                       ../modules/gaussian_fchk_v4.f90\
#                       ../modules/xyz_v4.f90\
#                       ../modules/molcas_unsym_v4.f90\
#                       ../modules/molecular_structure.f90\
#                       ../modules/ff_build_module_v4.f90\
#                       ../modules/atomic_geom_v2.f90\
#                       ../modules/symmetry_mod_v4.f90\
#                       ../modules/internal_SR9_v4_NvibAdd.f90\
#                       ../modules/zmat_v4.f90\
#                       zmat2modredundant.f90
# zmat2modredundant_LDFLAGS = @LIBS@

CLEANFILES =*.mod 
