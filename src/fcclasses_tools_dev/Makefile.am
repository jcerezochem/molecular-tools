AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = read_state_hess

read_state_hess_SOURCES = ../modules/alerts.f90\
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
                       ../modules/gamess_manage.f90\
                       ../modules/molecular_structure.f90\
                       ../modules/ff_build_module_v4.f90\
                       ../modules/atomic_geom_v2.f90\
                       ../modules/symmetry_mod_v4.f90\
                       ../modules/internal_SR8_v4.0.1.f90\
                       ../modules/zmat_manage_v4.f90\
                       read_state_hess.f90
 read_state_hess_LDFLAGS = @LIBS@
 read_state_hess_FCFLAGS = -cpp -DDOUBLE -DUSE_LAPACK


CLEANFILES =*.mod 



