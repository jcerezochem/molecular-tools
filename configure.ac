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

###################################
#
# COMPILER CHECKS
#
###################################
AC_LANG([Fortran])

#--------------
# Flags f90
#--------------
if test -z $FFLAGS; then
    FFLAGS='-I../modules/'
fi
if test -z $FCFLAGS; then
    if test "$FC" == "lf95"; then
        FCFLAGS='-DDOUBLE -DUSE_LAPACK -I../modules/ -Cpp --chk --maxfatals 1'
# -g -fcheck=all -fbounds-check' 
    else
        FCFLAGS='-DDOUBLE -DUSE_LAPACK -I../modules/ -cpp'
    fi
fi

#--------------
# Set compiler
#--------------
AC_PROG_FC([])


#########################
#
# LIBRARY CHECKS 
#
#########################
#----------------------
# Analytical Bders library
#----------------------
LDFLAGS="$LDFLAGS -L$srcdir/libs/AnalyticalBders"
# Check libbders
AC_ARG_WITH([libbders],
    AS_HELP_STRING([--with-libbders], [B analytical derivatives library (libbders)]),
    [],
    [with_libbders=no])

# We need to define the relative path to the libs. To that end we use
# \${top_srcdir}, that will be expanded at make step
# top_srcdir: is the top source, where the main Makefile file was called
# srcdr: is the folder where the current Makefile was called
AS_IF([test "x$with_libbders" = "xyes"],
      [AC_CHECK_LIB([bders],[derBstre],
                    [LIBS="$LIBS -L\${top_srcdir}/libs/AnalyticalBders -lbders";
                     FCFLAGS="$FCFLAGS -DWITH_LIBBDERS"],
                    [AC_MSG_WARN([libders disabled])])],
      [AC_MSG_NOTICE([libbders (analytical derivatives) not enabled])])
# Pass info to automake
AM_CONDITIONAL([WITH_LIBDERS], [test "$with_libbders" = "yes"])

#--------------------------------
# MKL support (LAPACK and BLAS)
#--------------------------------
AC_ARG_WITH([mkl],
    AS_HELP_STRING([--with-mkl], [Use Intel MKL library for lapak subroutines]),
    ,
    [with_mkl=no])

AS_IF([test "x$with_mkl" = "xyes"],
      # Using macro:
      # AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
      # where we use the last (optional) argument [other-libraries] to enter the additional libs
      # needed to use mkl lapack. (Previously this was done post-processing the configure: not nice)
      [AC_CHECK_LIB([mkl_intel_lp64],[zsytrf],
                    [LIBS="$LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"], # default: add libraries to LIBS
                    [# If it did not work, try addind pthread library (mkl compile with threads?)
                     AC_MSG_NOTICE([Trying MKL adding -lpthread]);
                     AC_CHECK_LIB([mkl_intel_lp64],[ssytrf],
                                  [LIBS="$LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"], # default: add libraries to LIBS
                                  [AC_MSG_WARN([MKL library requested but not found (LAPACK)]);use_lapack=yes],
                                  [-lmkl_sequential -lmkl_core -lpthread])
                    ],
                    [-lmkl_sequential -lmkl_core])],
      [use_lapack=yes])

AS_IF([test "x$with_mkl" = "xyes"],
      # Using macro:
      # AC_CHECK_LIB (library, function, [action-if-found], [action-if-not-found], [other-libraries])
      # where we use the last (optional) argument [other-libraries] to enter the additional libs
      # needed to use mkl lapack. (Previously this was done post-processing the configure: not nice)
      [AC_CHECK_LIB([mkl_intel_lp64],[dgemm],
                    [LIBS="$LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"], # default: add libraries to LIBS
                    [# If it did not work, try addind pthread library (mkl compile with threads?)
                     AC_MSG_NOTICE([Trying MKL adding -lpthread]);
                     AC_CHECK_LIB([mkl_intel_lp64],[dgemm],
                                  [LIBS="$LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"], # default: add libraries to LIBS
                                  [AC_MSG_WARN([MKL library requested but not found (BLAS)]);use_blas=yes],
                                  [-lmkl_sequential -lmkl_core -lpthread])
                    ],
                    [-lmkl_sequential -lmkl_core])],
      [use_blas=yes])
      
      
#-------------
# LAPACK
#-------------
AS_IF([test "x$use_lapack" = "xyes"],
      [AC_CHECK_LIB([lapack],[zsytrf],
                    [have_lapack_funct=yes;LIBS="$LIBS -llapack"],
                    [have_lapack_funct=no])],
      [have_lapack_funct=yes])

AS_IF([test "x$have_lapack_funct" != "xyes"],
      [AC_MSG_ERROR([No lapack library found])],
      [have_lapack_funct=no])
      
      
#-------------
# BLAS
#-------------
AS_IF([test "x$use_blas" = "xyes"],
      [AC_CHECK_LIB([blas],[dgemm],
                    [have_blas_funct=yes;LIBS="$LIBS -lblas"],
                    [have_blas_funct=no])],
      [have_blas_funct=yes])

AS_IF([test "x$have_blas_funct" != "xyes"],
      [AC_MSG_ERROR([No blas library found])],
      [have_blas_funct=no])
          
      
#-------------
# FFTW3
#-------------
AC_ARG_WITH([fftw],
    AS_HELP_STRING([--with-fftw], [Use the Fastes Fourier Transform in the West (FFTW)]),
    ,
    [with_fftw=no;use_fftw=no])
# FFTW is mandatory since FFT is broken for the corrected DFT -> FT manipulation (to be resolved..)
# with_fftw=yes

AS_IF([test "x$with_fftw" = "xyes"],
      [AC_CHECK_LIB([fftw3],[dfftw_plan_dft_1d],
                    [use_fftw=yes;LIBS="$LIBS -lfftw3"],
                    [AC_MSG_ERROR([No fftw library found])])],
      [use_fftw=no])
#The WARNING instead of the error when FFT is fixed will be:
#                    [AC_MSG_WARN([FFTW functionalities cannot be used]);use_fftw=no])],

# Set conditional to be used by automake 
AM_CONDITIONAL([FFTW], [test x$use_fftw = xyes])


################################
#
# PROCESS MAKEFILES IN FOLDERS
#
################################
AC_CONFIG_FILES([Makefile
                 libs/AnalyticalBders/Makefile
                 src/modules/Makefile
                 src/normal_modes/Makefile          
                 src/calc_duschinsky/Makefile
                 src/file_management/Makefile 
                 src/bonded_info/Makefile    
                 src/gen_oniom/Makefile    
                 src/wf_analysis/Makefile
                 src/build_top/Makefile    
                 src/joyce_tools/Makefile    
])        
AC_OUTPUT    
               

#                 src/vibrations/Makefile           
#                 src/cube_tools/Makefile                                    
#                 src/joyce_tools/Makefile    
  
