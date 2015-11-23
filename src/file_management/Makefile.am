AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = get_gscan fchk2gro get_symatoms reform_fchk rmsd_fit numderX

fchk2gro_SOURCES = ../modules/alerts.f90\
                      ../modules/MatrixMod.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gro_manage_v4.f90\
                      ../modules/pdb_manage_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/xyz_manage_v4.f90\
                      ../modules/molden_manage.f90\
                      ../modules/molcas_unsym_manage_v4.f90\
                      ../modules/psi4_manage.f90\
                      ../modules/gamess_manage.f90\
                      ../modules/molecular_structure.f90\
                      ../modules/ff_build_module_v4.f90\
                      fchk2gro_v4.1.f90
fchk2gro_LDFLAGS = @LIBS@

numderX_SOURCES = ../modules/alerts.f90\
                      ../modules/MatrixMod.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gro_manage_v4.f90\
                      ../modules/pdb_manage_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/xyz_manage_v4.f90\
                      ../modules/molden_manage.f90\
                      ../modules/molcas_unsym_manage_v4.f90\
                      ../modules/psi4_manage.f90\
                      ../modules/gamess_manage.f90\
                      ../modules/molecular_structure.f90\
                      ../modules/ff_build_module_v4.f90\
                      numderX.f90
numderX_LDFLAGS = @LIBS@

rmsd_fit_SOURCES = ../modules/alerts.f90\
                      ../modules/MatrixMod.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gro_manage_v4.f90\
                      ../modules/pdb_manage_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/xyz_manage_v4.f90\
                      ../modules/molden_manage.f90\
                      ../modules/molcas_unsym_manage_v4.f90\
                      ../modules/psi4_manage.f90\
                      ../modules/gamess_manage.f90\
                      ../modules/molecular_structure.f90\
                      ../modules/ff_build_module_v4.f90\
                      rmsd_fit.f90
rmsd_fit_LDFLAGS = @LIBS@

get_gscan_SOURCES = ../modules/alerts.f90\
                      ../modules/MatrixMod.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/xyz_manage_v4.f90\
                      ../modules/molden_manage.f90\
                      ../modules/molecular_structure.f90\
                      ../modules/ff_build_module_v4.f90\
                      ../modules/atomic_geom_v2.f90\
                      ../modules/symmetry_mod_v4.f90\
                      get_gscan.f90
get_gscan_LDFLAGS = @LIBS@

get_symatoms_SOURCES = ../modules/alerts.f90\
                      ../modules/MatrixMod.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gro_manage_v4.f90\
                      ../modules/pdb_manage_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/xyz_manage_v4.f90\
                      ../modules/molden_manage.f90\
                      ../modules/molcas_unsym_manage_v4.f90\
                      ../modules/molecular_structure.f90\
                      ../modules/ff_build_module_v4.f90\
                      ../modules/symmetry_mod_v4.f90\
                      get_symatoms.f90
get_symatoms_LDFLAGS = @LIBS@

# ============
# reform_fchk
# ============
reform_fchk_SOURCES = ../modules/alerts.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/MatrixMod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gro_manage_v4.f90\
                      ../modules/pdb_manage_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/xyz_manage_v4.f90\
                      ../modules/molcas_unsym_manage_v4.f90\
                      ../modules/molecular_structure.f90\
                      reform_fchk.f90
reform_fchk_LDFLAGS = @LIBS@

# gen_oniom_SOURCES = ../modules/alerts.f90\
#                       ../modules/line_preprocess.f90\
#                       ../modules/constants_mod.f90\
#                       ../modules/structure_types_v4.f90\
#                       ../modules/gro_manage_v4.f90\
#                       ../modules/pdb_manage_v4.f90\
#                       ../modules/gaussian_manage_v4.f90\
#                       ../modules/gaussian_fchk_manage_v4.f90\
#                       ../modules/xyz_manage_v4.f90\
#                       ../modules/molecular_structure.f90\
#                       gen_oniom_v4.f90
# gen_oniom_LDFLAGS = @LIBS@

CLEANFILES =*.mod 

