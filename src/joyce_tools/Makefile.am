AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = jpp  modify_fchk

# ============
# jpp
# ============
jpp_SOURCES = ../modules/alerts.f90\
              ../modules/line_preprocess.f90\
              ../modules/constants_mod.f90\
             jpp.f90
jpp_LDFLAGS = @LIBS@


# ============
# modify_fchk
# ============
modify_fchk_SOURCES = ../modules/alerts.f90\
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
                      modify_fchk.f90
modify_fchk_LDFLAGS = @LIBS@


CLEANFILES =*.mod 



