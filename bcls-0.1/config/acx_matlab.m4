dnl Look for Matlab.
AC_DEFUN([ACX_MATLAB], [
AC_PREREQ(2.50)

AC_ARG_WITH(matlab,
  [AC_HELP_STRING([--with-matlab=DIR], [the DIR where Matlab is installed])],
  MATLAB_DIR=${withval},
  MATLAB_DIR=)

if test -n "${MATLAB_DIR}"
then
   AC_MSG_CHECKING(for a Matlab installation)
   if test ! -e "${MATLAB_DIR}/bin/mex"
   then
      AC_MSG_RESULT(mex not found in $MATLAB_DIR/bin)
      have_matlab=no
   else

dnl   Determine architecture and mex extension.
dnl   Most of this logic was lifted from Matlab's own
dnl   $MATLABROOT/bin/util/arch.sh.
      MEX_EXT="unknown"
      case $build_os in
          solaris)	
              case $build_cpu in
                  sparc)
                      if [ "$SOLARIS_64" -eq "1" ]; then
                          MEX_EXT="sol64"
                      else
                          MEX_EXT="sol2"
                      fi
                      ;;
                  i386)
                      MEX_EXT="sola64"
                      ;;
              esac
              ;;
          linux)
              case $build_cpu in
                  i*86)
                      MEX_EXT="glnx86"
                      ;;
                  ia64)
                      MEX_EXT="glnxi64"
                      ;;
                  x86_64)
                      MEX_EXT="glnxa64"
                      ;;
              esac
              ;;
          Darwin)					
              case $build_cpu in
                  powerpc)
                      MEX_EXT="mac"
                      ;;
                  i386)
                      MEX_EXT="maci"
                      ;;
              esac
              ;;
          *)
              :
              ;;
      esac
dnl   END

      AC_MSG_RESULT(yes)
      MEX=$MATLAB_DIR/bin/mex

      AC_SUBST(MATLAB_DIR)
      AC_SUBST(MEX_EXT)
      AC_SUBST(MEX)
      have_matlab=yes
   fi
else
   have_matlab=no
fi
AM_CONDITIONAL(HAVE_MATLAB, test "x$have_matlab" = "xyes")
])dnl ACX_MATLAB
