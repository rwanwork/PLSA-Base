##  Prelude
AC_PREREQ([2.59])
AC_INIT([PLSA-BASE], [1.0.1], [r.wan@aist.go.jp])

AC_CONFIG_AUX_DIR(config-aux)
AM_INIT_AUTOMAKE

##  Set flags for C compiler and substitute them
##  Allowing the C99 standard so that we can declare variables in for loops
CFLAGS="-Wall -O3 -g -std=c99"
AC_SUBST(CFLAGS)

##  Checks for programs
AC_LANG(C)
AC_PROG_CC

##  Checks for header files
##  Checks for typedefs, structures, and compiler characteristics
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

##  Checks for library functions
AC_CHECK_HEADERS([float.h limits.h])
AC_FUNC_MALLOC
AC_FUNC_REALLOC

##  Output files
AC_CONFIG_HEADER([config.h])
AC_CONFIG_FILES([
  Makefile
])
AC_OUTPUT

