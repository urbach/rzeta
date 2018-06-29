
AC_INIT([src/qzeta_function.c])
AC_PROG_CXX

AC_CHECK_LIB(gsl, gsl_sf_airy_Ai_e,
             [LIBS="${LIBS} -lgsl -lgslcblas"], 
             [AC_MSG_ERROR([Cannot find Gnu Scientific Library.])],
             [-lgslcblas])
AC_CHECK_HEADER([gsl/gsl_version.h], ,
             [AC_MSG_ERROR([Cannot find GSL headers.])])
AC_MSG_CHECKING([if GSL version >= 1.8])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <string.h>
#include <gsl/gsl_version.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit ( major == 1 && minor < 8 );
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 1.8])
else
  AC_MSG_RESULT([yes])
fi

# if the user has not specified a library path for the quadmath library
# we try to auto-detect it
# currenty, if we fail to find it, we abort compilation
# in the future one could disable the quad precision
# implementation in this case
if test -z "${QUADMATH_LIB_DIR}"; then
  AC_SEARCH_LIBS([expq],
                 [quadmath],
                 [],
                 [AC_MSG_ERROR([unable to find the quadmath libary!])]
                 )
fi

gsl_custom_lib=""
gsl_custom_include=""
if test -n "${GSL_DIR}"; then
  gsl_custom_lib="-L${GSL_DIR}/lib -lgsl -lgslcblas"
  gsl_custom_include="-I${GSL_DIR}/include"
fi

quadmath_custom_lib=""
quadmath_custom_include=""
if test -n "${QUADMATH_LIB_DIR}"; then
  quadmath_custom_lib="-L${QUADMATH_LIB_DIR}/ -lquadmath"
fi
if test -n "${QUADMATH_INCLUDE_DIR}"; then
  quadmath_custom_include="-I${QUADMATH_INCLUDE_DIR}/"
fi

# construct src/Makevars with the manually or automatically detected library paths and flags
echo "PKG_LIBS=-lm ${LIBS} ${quadmath_custom_lib} ${gsl_custom_lib}" > src/Makevars
echo "PKG_CFLAGS=-march=native -mtune=native -O3 -Wall ${quadmath_custom_include} ${gsl_custom_include}" >> src/Makevars
echo "PKG_LDFLAGS=-lm ${LIBS} ${gsl_custom_lib} ${quadmath_custom_lib}" >> src/Makevars

AC_OUTPUT()

