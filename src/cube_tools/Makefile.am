AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS = charges_potential density_3D
# solvent_potential

#solvent_potential_SOURCES = ../modules/alerts.f90\
#                      ../modules/line_preprocess.f90\
#                      ../modules/constants_mod.f90\
#                      sizes.f90\
#                      ../modules/structure_types_v5.f90\
#                      ../modules/gro_manage_v4.f90\
#                      ../modules/pdb_manage_v4.f90\
#                      ../modules/gaussian_manage_v4.f90\
#                      ../modules/gaussian_fchk_manage_v4.f90\
#                      ../modules/xyz_manage_v4.f90\
#                      ../modules/molecular_structure.f90\
#                      solvent_potential.f90
#solvent_potential_LDFLAGS = @LIBS@
#solvent_potential_FCFLAGS = -cpp -DDOUBLE -O3


charges_potential_SOURCES = charges_potential.f90
charges_potential_FCFLAGS = -O3


density_3D_SOURCES = density_3D.f90
density_3D_FCFLAGS = -O3

