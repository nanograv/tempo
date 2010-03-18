AC_INIT([tempo], [12.000])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CC

AC_PROG_F77
if test "$F77" = gfortran; then
    FFLAGS+=" -ffixed-form -ffixed-line-length-none -std=legacy"
fi

AC_C_BIGENDIAN([BIGENDIAN=".true."], [BIGENDIAN=".false."])
AC_SUBST(BIGENDIAN)

AC_CHECK_HEADERS([malloc.h])

AC_CONFIG_FILES([
	Makefile
	src/Makefile
	tempo.cfg
	src/config.h
	clock/Makefile
	ephem/Makefile
	test/Makefile
	tzpar/Makefile
	util/print_resid/Makefile
])
AC_OUTPUT

# Misc messages
echo 
echo "Tempo is ready to be compiled.  Please run 'make'."
echo 
echo "'make install' will copy 'tempo' to $prefix/bin."
echo
