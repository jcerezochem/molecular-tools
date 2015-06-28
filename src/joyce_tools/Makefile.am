AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = jpp 

jpp_SOURCES = ../modules/alerts.f90\
              ../modules/line_preprocess.f90\
              ../modules/constants_mod.f90\
             jpp.f90
jpp_LDFLAGS = @LIBS@
jpp_FCFLAGS = -cpp -DDOUBLE

CLEANFILES =*.mod 



