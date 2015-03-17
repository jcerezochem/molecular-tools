AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = check_diag 

check_diag_SOURCES =  ../modules/alerts.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/MatrixMod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/molecular_structure.f90\
                      ../modules/gro_manage_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/molcas_unsym_manage_v4.f90\
                      ../modules/psi4_manage.f90\
                      ../modules/gamess_manage.f90\
                      check_diag_v4.f90
check_diag_LDFLAGS = @LIBS@
check_diag_FCFLAGS = -cpp -DDOUBLE -DUSE_LAPACK


CLEANFILES =*.mod 

