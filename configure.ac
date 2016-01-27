#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([MolecularTools],[2.2.1],[jcb1@um.es])
AM_INIT_AUTOMAKE([foreign])
AM_EXTRA_RECURSIVE_TARGETS([clean-libs])
AC_CONFIG_SRCDIR(src/)
AC_LANG([Fortran])

# # Using Libtool to pack all module objects
# AC_CONFIG_MACRO_DIR([m4])
# LT_INIT
# Or using ranlig
AC_PROG_RANLIB([ranlib])

# Javier specific things
AC_PREFIX_DEFAULT(/home/cerezo/Programas/MolecularTools)

# Checks for programs.
# Prevent autoconf from adding flags automatically  (e.g. -g)
if test -z $FFLAGS; then
    FFLAGS='-I../modules/'
fi
if test -z $FCFLAGS; then
    FCFLAGS='-cpp -DDOUBLE -DUSE_LAPACK -I../modules/'
fi

LDFLAGS="$LDFLAGS -L$srcdir/libs/AnalyticalBders"

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
                    [use_lapack=no;LIBS="$LIBS -lmkl_intel"],
                    [AC_MSG_WARN([MKL library requested but not found]);use_lapack=yes])],
      # MKL test needs some changes in the configure file:
      # sed -i "s/-lmkl_intel/-lmkl_intel_lp64 -lmkl_sequential -lmkl_core/g" configure
      [use_lapack=yes])

AS_IF([test "x$use_lapack" = "xyes"],
      [AC_CHECK_LIB([lapack],[zsytrf],
                    [have_lapack_funct=yes;LIBS="$LIBS -llapack -lblas"],
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

# Check libbders
AC_ARG_WITH([libbders],
    AS_HELP_STRING([--with-libbders], [B analytical derivatives library (libbders)]),
    [with_libbders=yes],
    [with_libbders=no; use_libbders=no])

AS_IF([test "x$with_libbders" = "xyes"],
      [AC_CHECK_LIB([bders],[derBstre],
                    [use_libders=yes;
                     LIBS="$LIBS -lbders -L${abs_srcdir}/libs/AnalyticalBders";
                     FCFLAGS="$FCFLAGS -DWITH_LIBBDERS"],
                    [AC_MSG_WARN([libders disabled]);use_libders=no])],
      [use_libders=yes])
# Pass info to automake
AM_CONDITIONAL([WITH_LIBDERS], [test "$use_libders" = "yes"])


AS_IF([test "x$use_fftw" = "xyes"],
      [FFT_SR_SRC=fftw_calculator_SR.f90],
      [FFT_SR_SRC=fft_calculator_SR.f90 ])
# Make the substitution in Makefile
AC_SUBST([FFT_SR_SRC])
# AC_SUBST([FFT_SRC])
# Set a conditional for FFTW
AM_CONDITIONAL([FFTW],[test x$use_fftw = xyes])


# LIBS="$LIBS -L/home/cerezo/dev/MY_PROJECTS/GIT/molecular-tools/AnalyticalBders -lbders"
# LIBS="$LIBS -L$top_srcdir/libs/AnalyticalBders -lbders"

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 libs/AnalyticalBders/Makefile
                 src/modules/Makefile
                 src/normal_modes/Makefile          
                 src/calc_duschinsky/Makefile
                 src/file_management/Makefile 
                 src/bonded_info/Makefile    
])        
AC_OUTPUT    
               

#                 src/vibrations/Makefile           
#                 src/wf_analysis/Makefile
#                 src/gen_oniom/Makefile
#                 src/cube_tools/Makefile                                    
#                 src/build_top/Makefile    
#                 src/joyce_tools/Makefile    
  
