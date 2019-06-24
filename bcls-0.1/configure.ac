#-------------------------------------------------------------------
# Process this file with autoconf to produce a configure script
#-------------------------------------------------------------------
AC_PREREQ(2.59)
AC_COPYRIGHT([Copyright (C) 2006 Michael P. Friedlander <mpf@cs.ubc.ca>])
AC_INIT([BCLS], [0.1], [mpf@cs.ubc.ca], [bcls])
AC_REVISION($Revision$)
AC_CONFIG_SRCDIR([src/bcls.c])
AC_CONFIG_AUX_DIR([config])
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE

# Shared library versioning
BCLS_LIBRARY_VERSION="1:0:0"
#                     | | |
#              +------+ | +---+
#              |        |     |
#           current:revision:age
#              |        |     |
#              |        |     +- increment if interfaces have been added;
#              |        |        set to zero if interfaces have been 
#              |        |        removed or changed
#              |        +- increment if source code has changed;
#              |           set to zero if current is incremented
#              +- increment if interfaces have been added, removed
#                 or changed
AC_SUBST(BCLS_LIBRARY_VERSION)

# Checks for programs
AC_LANG([C])
AC_PROG_MAKE_SET
save_CFLAGS="$CFLAGS"
AC_PROG_CC
CFLAGS="$save_CFLAGS"
AC_PROG_INSTALL

# Diable unnecessary tests for C++, Fortran, and Java
# before checking for libtool.
define([AC_LIBTOOL_LANG_CXX_CONFIG], [:])dnl
define([AC_LIBTOOL_LANG_F77_CONFIG], [:])dnl
define([AC_LIBTOOL_LANG_GCJ_CONFIG], [:])dnl
AC_PROG_LIBTOOL

# Check for debug switch
AC_ARG_ENABLE(debugging, 
              [AC_HELP_STRING([--enable-debugging],
                              [enable debugging [default=no]])])
AM_CONDITIONAL([ENABLE_DEBUGGING], [test x"$enable_debugging" = xyes])

# Check for htmldocs switch
AC_ARG_ENABLE(html-doc, 
              [AC_HELP_STRING([--enable-html-doc],
                              [enable HTML documentation [default=yes]])])
test "x$enable_html_doc" != xno && enable_html_doc=yes

# Debugging flags.
CFLAGS_DEBUG="-g $CFLAGS"
if test "x$GCC" = xyes; then
   CFLAGS_DEBUG="$CFLAGS_DEBUG -Wall -pedantic -std=c99"
fi
AC_SUBST(CFLAGS_DEBUG)

# If the user hasn't set CFLAGS, then we'll touch it!
if test "x$CFLAGS" = x; then
   CFLAGS="-O"
   # Only for GCC compiler
   if test "x$GCC" = xyes; then
      CFLAGS="$CFLAGS -std=c99"
   fi
fi

AC_CHECK_LIB(m, sqrt)

# Checks for header files
AC_HEADER_STDC
AC_HEADER_STDBOOL
AC_CHECK_HEADERS([stddef.h string.h sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([getrusage])

# Checks for MATLAB.
ACX_MATLAB

# Only link examples/bcsol with a foreign CBLAS if specifically requested.
AC_ARG_WITH(cblas,
  [AC_HELP_STRING([--with-cblas=<lib>], [use CBLAS library <lib>])])
case $with_cblas in
  no) acx_cblas_ok=no ;;
  -* | */* | *.a | *.so | *.so.* | *.o) CBLAS_LIBS="$with_cblas" ;;
  *) CBLAS_LIBS="-l$with_cblas" ;;
esac

AM_CONDITIONAL([ENABLE_LOCAL_CBLAS], [test x"$acx_cblas_ok" != xyes])

# Host-specific switches:
case $host in
     *-*-msdos* | *-*-go32* | *-*-mingw32* | *-*-cygwin* | *-*-windows*)
     AC_DEFINE(WINDOWS_HOST, 1, [Defined if host is Windows])
esac

# Checks for documentation build tools
AC_CHECK_PROG([doxygen_ok], [doxygen], [yes], [no])
AC_CHECK_PROG([latex_ok], [latex], [yes], [no])
AC_CHECK_PROG([dvips_ok], [dvips], [yes], [no])
AC_CHECK_PROG([gs_ok], [gs], [yes], [no])
if test "x$enable_html_doc" != xno; then
  test "x$doxygen_ok" != xyes && enable_html_doc=no
  test "x$latex_ok"   != xyes && enable_html_doc=no
  test "x$dvips_ok"   != xyes && enable_html_doc=no
  test "x$gs_ok"      != xyes && enable_html_doc=no
  AC_PATH_PROG([PERL], [perl], [])
fi
AM_CONDITIONAL([HTML_DOC], [test "x$enable_html_doc" = xyes])

# Define configured files
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 examples/Makefile
                 matlab/Makefile
                 doc/Makefile
                 doc/Doxyfile
                 doc/local/installation.doc])

AC_OUTPUT

echo "
------------------------------------------------------------------------------
$PACKAGE-$VERSION library configuration:
------------------------------------------------------------------------------

Directories:
  - prefix ......... : ${prefix}
  - exec_prefix .... : ${exec_prefix}
  - includedir ..... : ${includedir}
  - libdir ......... : ${libdir}

Switches:
  - shared ......... : ${enable_shared}
  - static ......... : ${enable_static}

Matlab: ............. : ${have_matlab}
  - Root ............ : ${MATLAB_DIR}
  - MEX extension ... : ${MEX_EXT}

Documentation tools:
  - doxygen ........ : ${doxygen_ok}
  - latex .......... : ${latex_ok}
  - dvips .......... : ${dvips_ok}
  - ghostscript .... : ${gs_ok}

External libs:
  - CBLAS .......... : ${with_cblas}

Compiler/linker flags/libs/defs:
  - CC ............. : ${CC}
  - CFLAGS ......... : ${CFLAGS}
  - CFLAGS_DEBUG ... : ${CFLAGS_DEBUG}
  - LDFLAGS ........ : ${LDFLAGS}
  - LIBS ........... : ${LIBS}

Host information:
  - Vendor.......... : ${build_vendor}
  - OS.............. : ${build_os}
  - CPU............. : ${build_cpu}
------------------------------------------------------------------------------
$PACKAGE-$VERSION successfully configured.
To build an optimized library,                 type  'make'
To build a library that can be debugged,       type  'make debug' 
To install the library,                        type  'make install'
To compile a Matlab interface,                 type  'make mex'
To build html docs (if Doxygen is installed),  type  'make html-doc'
To blow your nose,                             type  'make honk'
------------------------------------------------------------------------------
"
