dnl     Adds to the include-path
dnl
dnl     Autoconf 1.11 should have provided a way to add include path options to
dnl     the cpp-tests.
dnl
AC_DEFUN([CF_INCLUDE_PATH],
[
for cf_path in $1
do
        cf_result="no"
        AC_MSG_CHECKING(for directory $cf_path)
        if test -d $cf_path
        then
                INCLUDES="$INCLUDES -I$cf_path"
                ac_cpp="${ac_cpp} -I$cf_path"
                CFLAGS="$CFLAGS -I$cf_path"
                cf_result="yes"
                case $cf_path in
                /usr/include|/usr/include/*)
                        ;;
                *)
changequote(,)dnl
                        cf_temp=`echo $cf_path | sed -e s'%/[^/]*$%%'`
changequote([,])dnl
                        case $cf_temp in
                        */include)
                                INCLUDES="$INCLUDES -I$cf_temp"
                                ac_cpp="${ac_cpp} -I$cf_temp"
                                CFLAGS="$CFLAGS -I$cf_temp"
                                ;;
                        esac
                esac
        fi
        AC_MSG_RESULT($cf_result)
done
])dnl

AC_DEFUN(AC_PROG_CC2,
[AC_BEFORE([$0], [AC_PROG_CPP2])dnl
AC_CHECK_PROG(CC, pgcc, pgcc, , , /usr/ucb/cc)
AC_CHECK_PROG(CC, cc, cc, , , /usr/ucb/cc)
AC_CHECK_PROG(CC, c89, c89, , , /usr/ucb/cc)
AC_CHECK_PROG(CC, gcc, gcc)

test -z "$CC" && AC_MSG_ERROR([no acceptable cc found in \$PATH])

AC_PROG_CC_WORKS
AC_PROG_CC_GNU

if test $ac_cv_prog_gcc = yes; then
  GCC=yes
dnl Check whether -g works, even if CFLAGS is set, in case the package
dnl plays around with CFLAGS (such as to build both debugging and
dnl normal versions of a library), tasteless as that idea is.
  ac_test_CFLAGS="${CFLAGS+set}"
  ac_save_CFLAGS="$CFLAGS"
  CFLAGS=
  AC_PROG_CC_G
  if test "$ac_test_CFLAGS" = set; then
    CFLAGS="$ac_save_CFLAGS"
  elif test $ac_cv_prog_cc_g = yes; then
    CFLAGS="-g"
  else
    CFLAGS=""
  fi
else
  GCC=
fi
])
dnl
dnl  Replacement for AC_PROG_CPP which tries rather harder to find a good
dnl  cpp.  
dnl
AC_DEFUN(AC_PROG_CPP2,
[AC_PATH_PROG(cpp_tmp,cpp, ,/usr/ccs/lib:/usr/lang:/opt/ctl/bin:/lib)
 AC_MSG_CHECKING(how to run the C preprocessor)
 AC_PROVIDE(AC_PROG_CPP)
# On Suns, sometimes $CPP names a directory.
if test -n "$CPP" && test -d "$CPP"; then
  CPP=
fi
if test -z "$CPP"; then
AC_CACHE_VAL(ac_cv_prog_CPP,
[  # This must be in double quotes, not single quotes, because CPP may get
  # substituted into the Makefile and "${CC-cc}" will confuse make.
  CPP="${CC-cc} -E"
  # On the NeXT, cc -E runs the code through the compiler's parser,
  # not just through cpp.
dnl Use a header file that comes with gcc, so configuring glibc
dnl with a fresh cross-compiler works.
  AC_TRY_CPP([#include <assert.h>
Syntax Error], ,
  CPP="${CC-cc} -E -traditional-cpp"
  AC_TRY_CPP([#include <assert.h>
Syntax Error], ,
  if test ${ac_cv_prog_gcc-no} = yes; then CPP="gcc -E"; else CPP="cc -E"; fi
  AC_TRY_CPP([#include <assert.h>
Syntax Error], ,CPP=$cpp_tmp)))
  ac_cv_prog_CPP="$CPP"])dnl
  CPP="$ac_cv_prog_CPP"
else
  ac_cv_prog_CPP="$CPP"
fi
AC_MSG_RESULT($CPP)
AC_SUBST(CPP)dnl
])
