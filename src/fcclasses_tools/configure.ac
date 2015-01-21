#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([fcclasses_tools],[2.2.1],[jcb1@um.es])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR(./)
AC_LANG([Fortran])
# Javier specific things
AC_PREFIX_DEFAULT(/home/cerezo/Programas/MolecularTools)

# Checks for programs.
# Prevent autoconf from adding flags automatically  (e.g. -g)
if test -z $FFLAGS; then
    FFLAGS=''
fi
if test -z $FCFLAGS; then
    FCFLAGS=''
fi
# AC_PROG_F77([ifort,gfortran])
 # ifort compiler works, but gfortran does not!
AC_PROG_FC([])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

