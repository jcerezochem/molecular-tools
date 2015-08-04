AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = contribMO get_NTO_coef

contribMO_SOURCES = ../modules/alerts.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/atomic_geom_v2.f90\
                      contribMO_v4.f90
contribMO_LDFLAGS = @LIBS@

get_NTO_coef_SOURCES = ../modules/alerts.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      ../modules/structure_types_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      get_NTO_coef.f90
get_NTO_coef_LDFLAGS = @LIBS@

CLEANFILES =*.mod 



