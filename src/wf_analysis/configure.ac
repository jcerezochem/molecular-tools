#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([fcclasses_tools],[2.2.1],[jcb1@um.es])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR(contribMO_v4.f90)
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

# Checks for libraries: LAPACK
AC_ARG_WITH([mkl],
    AS_HELP_STRING([--with-mkl], [Use Intel MKL library for lapak subroutines]),
    [with_mkl=yes],
    [with_mkl=no])

AS_IF([test "x$with_mkl" = "xyes"],
      [AC_CHECK_LIB([mkl_intel],[zsytrf],
                    [use_lapack=no;LIBS="-lmkl_intel"],
                    [AC_MSG_WARN([MKL library requested but not found]);use_lapack=yes])],
      # MKL test needs some changes in the configure file:
      # sed -i "s/-lmkl_intel/-lmkl_intel_lp64 -lmkl_sequential -lmkl_core/g" configure
      [use_lapack=yes])

AS_IF([test "x$use_lapack" = "xyes"],
      [AC_CHECK_LIB([lapack],[zsytrf],
                    [have_lapack_funct=yes;LIBS=-llapack],
                    [have_lapack_funct=no])],
      [have_lapack_funct=yes])

AS_IF([test "x$have_lapack_funct" != "xyes"],
      [AC_MSG_ERROR([No lapack library found])],
      [have_lapack_funct=yes])

# Checks for libraries: FFT
AC_ARG_WITH([fftw],
    AS_HELP_STRING([--with-fftw], [Use the Fastes Fourier Transform in the West (FFTW)]),
    [with_fftw=yes],
    [with_fftw=no;use_fftw=no])

AS_IF([test "x$with_fftw" = "xyes"],
      [AC_CHECK_LIB([fftw3],[dfftw_plan_dft_1d],
                    [use_fftw=yes;LIBS="$LIBS -lfftw3"],
                    [AC_MSG_WARN([FFTW functionalities cannot be used]);use_fftw=no])],
      [use_fftw=no])

AS_IF([test "x$use_fftw" = "xyes"],
      [FFT_SR_SRC=fftw_calculator_SR.f90],
      [FFT_SR_SRC=fft_calculator_SR.f90 ])
# Make the substitution in Makefile
AC_SUBST([FFT_SR_SRC])
# AC_SUBST([FFT_SRC])
# Set a conditional for FFTW
AM_CONDITIONAL([FFTW],[test x$use_fftw = xyes])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

