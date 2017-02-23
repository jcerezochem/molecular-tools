AUTOMAKE_OPTIONS = subdir-objects

LDFLAGS = -L$(top_builddir)/libs/AnalyticalBders @LIBS@ 

bin_PROGRAMS = vertical2adiabatic cartesian_duschinsky internal_duschinsky generate_internals generate_connetivity state2hessian \
               rotate_disp substract_disps transpose_dusch NewtonRaphson_step get_DeltaX

# Simple tools:
rotate_disp_SOURCES = rotate_disp.f90
substract_disps_SOURCES = substract_disps.f90
transpose_dusch_SOURCES = transpose_dusch.f90


state2hessian_SOURCES = \
                        ../modules/alerts.f90\
                        ../modules/line_preprocess.f90\
                        ../modules/constants.f90\
                        ../modules/verbosity.f90\
                         ../modules/io.f90                      \
                         ../modules/version.f90                 \
                        ../modules/matrix.f90\
                        ../modules/matrix_print.f90\
                        ../modules/structure_types.f90\
                        ../modules/molecular_structure.f90\
                        ../modules/metrics.f90             \
                        ../modules/atomic_geom.f90\
                        ../modules/manage_fcc.f90                    \
                        ../modules/manage_gamess.f90                 \
                        ../modules/manage_gaussian.f90               \
                        ../modules/manage_gmx.f90                    \
                        ../modules/manage_molcas.f90                 \
                        ../modules/manage_molden.f90                 \
                        ../modules/manage_molpro.f90                 \
                        ../modules/manage_turbomol.f90               \
                        ../modules/manage_pdb.f90                    \
                        ../modules/manage_psi4.f90                   \
                        ../modules/manage_xyz.f90                    \
                        ../modules/manage_xyz_molec.f90              \
                        ../modules/manage_gro_molec.f90              \
                        ../modules/zmat_manage.f90                   \
                        ../modules/generic_io.f90\
                        ../modules/generic_io_molec.f90\
                        ../modules/ff_build.f90\
                        ../modules/symmetry.f90\
                        ../modules/internal.f90\
                        ../modules/vibrational_analysis.f90\
                        ../modules/vertical_model.f90\
                        state2hessian.f90

vertical2adiabatic_SOURCES = \
                             ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/io.f90                      \
                             ../modules/version.f90                 \
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_turbomol.f90               \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/manage_xyz_molec.f90              \
                             ../modules/manage_gro_molec.f90              \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/ff_build.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             ../modules/vibrational_analysis.f90\
                             ../modules/vertical_model.f90\
                             vertical2adiabatic.f90

NewtonRaphson_step_SOURCES = \
                             ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/io.f90                      \
                             ../modules/version.f90                 \
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_turbomol.f90               \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/manage_xyz_molec.f90              \
                             ../modules/manage_gro_molec.f90              \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/ff_build.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             ../modules/vibrational_analysis.f90\
                             ../modules/vertical_model.f90\
                             NewtonRaphson_step.f90

cartesian_duschinsky_SOURCES = \
                             ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/io.f90                      \
                             ../modules/version.f90                 \
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_turbomol.f90               \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/manage_xyz_molec.f90              \
                             ../modules/manage_gro_molec.f90              \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/ff_build.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             ../modules/vibrational_analysis.f90\
                             ../modules/vertical_model.f90\
                             cartesian_duschinsky.f90

internal_duschinsky_SOURCES = \
                             ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/io.f90                      \
                             ../modules/version.f90                 \
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_turbomol.f90               \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/manage_xyz_molec.f90              \
                             ../modules/manage_gro_molec.f90              \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/ff_build.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             ../modules/vibrational_analysis.f90\
                             ../modules/vertical_model.f90\
                             internal_duschinsky_prjCrt.f90

get_DeltaX_SOURCES = \
                             ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/io.f90                      \
                             ../modules/version.f90                 \
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_turbomol.f90               \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/manage_xyz_molec.f90              \
                             ../modules/manage_gro_molec.f90              \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/ff_build.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             ../modules/vibrational_analysis.f90\
                             ../modules/vertical_model.f90\
                             get_DeltaX.f90


generate_internals_SOURCES = ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/io.f90\
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_turbomol.f90               \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/manage_xyz_molec.f90              \
                             ../modules/manage_gro_molec.f90              \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/ff_build.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             generate_internals.f90 

generate_connetivity_SOURCES = ../modules/alerts.f90\
                             ../modules/line_preprocess.f90\
                             ../modules/constants.f90\
                             ../modules/verbosity.f90\
                             ../modules/io.f90\
                             ../modules/matrix.f90\
                             ../modules/matrix_print.f90\
                             ../modules/structure_types.f90\
                             ../modules/metrics.f90             \
                             ../modules/atomic_geom.f90\
                             ../modules/molecular_structure.f90\
                             ../modules/manage_fcc.f90                    \
                             ../modules/manage_gamess.f90                 \
                             ../modules/manage_gaussian.f90               \
                             ../modules/manage_gmx.f90                    \
                             ../modules/manage_molcas.f90                 \
                             ../modules/manage_molden.f90                 \
                             ../modules/manage_molpro.f90                 \
                             ../modules/manage_turbomol.f90               \
                             ../modules/manage_pdb.f90                    \
                             ../modules/manage_psi4.f90                   \
                             ../modules/manage_xyz.f90                    \
                             ../modules/manage_xyz_molec.f90              \
                             ../modules/manage_gro_molec.f90              \
                             ../modules/zmat_manage.f90                   \
                             ../modules/generic_io.f90\
                             ../modules/generic_io_molec.f90\
                             ../modules/ff_build.f90\
                             ../modules/symmetry.f90\
                             ../modules/internal.f90\
                             generate_connetivity.f90 

CLEANFILES =*.mod 
