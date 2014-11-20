bin_PROGRAMS = nm_internal nm_cartesian freq2fc prj_nm_internal

# ============
# nm_internal
# ============
nm_internal_SOURCES = ../modules/alerts.f90\
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
                      ../modules/ff_build_module_v4.f90\
                      ../modules/atomic_geom_v2.f90\
                      ../modules/symmetry_mod_v4.f90\
                      ../modules/internal_SR8_v4.0.2.f90\
                      ../modules/zmat_manage_v4.f90\
                      normal_modes_internal_v4.0.1.2.f90
nm_internal_LDFLAGS = @LIBS@
nm_internal_FCFLAGS = -cpp -DDOUBLE -DUSE_LAPACK

# ============
# nm_cartesian
# ============
nm_cartesian_SOURCES =../modules/alerts.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/MatrixMod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gro_manage_v4.f90\
                      ../modules/pdb_manage_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/gaussian_manage_vNoTypes.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/xyz_manage_v4.f90\
                      ../modules/molecular_structure.f90\
                      ../modules/ff_build_module_v4.f90\
                      ../modules/atomic_geom_v2.f90\
                      normal_modes_Cartesian1b_v4.0.1.f90
nm_cartesian_LDFLAGS = @LIBS@
nm_cartesian_FCFLAGS = -cpp -DDOUBLE

# ========
# freq2fc
# ========
freq2fc_SOURCES =../modules/alerts.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      freq2fc.f90
freq2fc_LDFLAGS = @LIBS@
freq2fc_FCFLAGS = -cpp -DDOUBLE

# =========================
# internal_numder_analysis
# =========================
# internal_numder_analysis_v4.0.1.f90

# ============
# prj_nm_internal
# ============
prj_nm_internal_SOURCES = ../modules/alerts.f90\
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
                      ../modules/ff_build_module_v4.f90\
                      ../modules/atomic_geom_v2.f90\
                      ../modules/symmetry_mod_v4.f90\
                      ../modules/internal_SR8_v4.0.2.f90\
                      ../modules/zmat_manage_v4.f90\
                      project_normal_modes_internal.f90
prj_nm_internal_LDFLAGS = @LIBS@
prj_nm_internal_FCFLAGS = -cpp -DDOUBLE -DUSE_LAPACK


CLEANFILES =*.mod 

