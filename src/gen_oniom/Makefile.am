AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = gen_oniom

gen_oniom_SOURCES = ../modules/alerts.f90\
                      ../modules/line_preprocess.f90\
                      ../modules/constants_mod.f90\
                      sizes.f90\
                      ../modules/structure_types_v5.f90\
                      ../modules/gro_manage_geno.f90\
                      ../modules/pdb_manage_v4.f90\
                      ../modules/gaussian_manage_v4.f90\
                      ../modules/gaussian_fchk_manage_v4.f90\
                      ../modules/xyz_manage_v4.f90\
                      ../modules/molecular_structure.f90\
                      gen_oniom_v5.f90
gen_oniom_LDFLAGS = @LIBS@
gen_oniom_FCFLAGS = -cpp -DDOUBLE

CLEANFILES =*.mod 

