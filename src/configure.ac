dnl Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
AC_INIT([Moldy],[3.6])
AC_CONFIG_SRCDIR([accel.c])

if test "x$CC" != "x" ; then
   for c in cc gcc c89 ; do
     ($c dummy.c -c -o $dummy.o) >/dev/null 2>&1
     if test $? = 0 ; then
         CC_FOR_BUILD="$c"; break
     fi
   done
   export CC_FOR_BUILD
fi

AC_CANONICAL_TARGET([])

AM_INIT_AUTOMAKE()
AM_CONFIG_HEADER(config.h)

dnl Checks for programs.

dnl Select parallel compilation
dnl
dnl Set cross-compiling flag because we can't guarantee to be able to run
dnl a parallel prog from configure.
AC_ARG_ENABLE(parallel, --enable-parallel: Compile for parallelism using one of mpi/bsp/shmem/tcgmsg (default is mpi), 
if test $enableval = yes ; then enableval="mpi" ; fi
cross_compiling="yes"
parlib=$enableval, parlib=none
)

if test $parlib = "mpi"; then
   compilers="mpicc mpcc hcc icc pgcc cc c89 gcc egcs"
elif test $parlib = "bsp"; then
   compilers="bspcc pgcc icc cc c89 gcc egcs"
else
   compilers="icc pgcc cc c89 gcc egcs"
fi

AC_PROG_CC($compilers)

if test "${CC-cc}" != gcc; then
dnl 
dnl Check for C compiler. DEFCC is preferred compiler, ALTCC is a list
dnl  of alternatives *which take the same flags*
dnl
   case "${target}" in
      alpha*-dec-osf*)    DEFCC="cc -migrate" ALTCC="mpicc bspcc hcc cc c89";;
      hppa?.*-hp-*)       DEFCC=c89 ALTCC="mpicc bspcc hcc cc";;
      rs6000-ibm-aix*)    DEFCC=cc ALTCC="mpcc mpicc bspcc hcc xlc c89" ;;
      sparc-sun-sunos*)   DEFCC=cc ALTCC="mpicc bspcc hcc acc";;
      sparc-sun-solaris*) DEFCC=cc ALTCC="mpicc bspcc hcc c89";;
      mips-sgi-irix*)     DEFCC=cc ALTCC="mpicc bspcc hcc c89";;
dnl   
dnl Unfortunately the Bourne-shell "env=value command" construct breaks the 
dnl configure script.
dnl However we can rely on the "env" command in the case of the crays.
dnl   
      alpha-cray-unicosmk*)   DEFCC="env TARGET=cray-t3d cc" ALTCC="mpicc bspcc cc";;
      alphaev*-cray-unicosmk*) DEFCC="env TARGET=cray-t3e cc" ALTCC="mpicc bspcc cc";;
dnl
dnl Alternative compilers for Linux
dnl
      i?86-pc-linux*)     DEFCC=pgcc ALTCC="icc mpicc bspcc hcc";;
      x86_64-*-linux*)    DEFCC=pgcc ALTCC="icc mpicc bspcc hcc";;
      ia64-*-linux*)      DEFCC=pgcc ALTCC="icc mpicc bspcc hcc";;
                     *)   ALTCC="cc mpicc bspcc hcc";;
   esac
dnl
dnl Here we test whether the specified (or default) C compiler is one of those
dnl we know a good set of optimization flags for.
dnl
   CC="${CC=$DEFCC}"
   cc_matches=0
   for trial_cc in "$DEFCC" $ALTCC ; do
        if test "$trial_cc" = "$CC" ; then cc_matches=1 ; fi
   done
   if test $cc_matches -eq 1 ; then
    case "${target}" in
     alpha*-dec-osf*)   : ${OPT=-tune\ host\ -O4}  \
                          ${OPT2=-tune\ host\ -ansi_alias\ -O5} \
                          ${OPTAUX=$OPT2};;
     #sparc-sun-sunos*) : ${OPT=-fast\ -xO4}  \
     #                    ${OPT2=-fast\ -xO4\ -fsimple}  \
     #                    ${LDFLAGS=-fast} ;;
     sparc-sun-sunos*)  : ${OPT=-fast\ -O4} ${OPT2=$OPT} ${OPTAUX=$OPT2}\
                          ${LDFLAGS=-fast} ;;
     sparc-sun-solaris*): ${OPT=-fast\ -xO4}  \
                          ${OPT2=${OPT}\ -fsimple\ -xrestrict=%all} \
                          ${OPTAUX=$OPT2} ${LDFLAGS=-fast} ;;
     rs6000-ibm-aix*)   : ${OPT=-qansialias\ -O3}  \
                          ${OPT2=$OPT} ${OPTAUX=$OPT2};;
     mips-sgi-irix5*)   : ${OPT=-O2\ -mips2} ${OPT2=$OPT} ${OPTAUX=$OPT2}\
                          ${LIBS=-lmalloc\ -lfastm} ;;
     mips-sgi-irix6*)   : ${EXTRA_CFLAGS=-woff\ 1174,84,1552,1680} ${LDFLAGS=-woff\ 84} \
                          ${OPT=-O3} \
                          ${OPT2=${OPT}\ -OPT:alias=restrict:fast_sqrt=on:roundoff=3:IEEE_arithmetic=3}} \
                          ${OPTAUX=$OPT2};;
     hppa1.*-hp-*)      : ${OPT=+O3\ +Olibcalls}  \
                          ${OPT2=${OPT}\ +Onoparmsoverlap} ${OPTAUX=$OPT2} ;;
     hppa2.*-hp-*)      : ${OPT=+O3\ +Olibcalls}  \
                          ${OPT2=${OPT}\ +Onoparmsoverlap\ +Odataprefetch\ +Oentrysched\ +Oregionsched} ${OPTAUX=$OPT2} ;;
     c*-convex-bsd*)    : ${OPT=-O2\ -na\ -nv\ -fi}  \
                          ${OPT2=${OPT}\ -alias\ standard\ -alias\ ptr_args\ -alias\ array_args} \
                          ${OPTAUX=$OPT2} ${LIBS=-fi\ -lveclib}  ;;
     *mp-cray-unicos*|*90-cray-unicos*|sv1-cray-unicos*) 
                        : ${OPT=-h\ nostdc,vector3,scalar3,ivdep}  \
                          ${OPT2=$OPT} ${OPTAUX=$OPT2} \
                          ${LDFLAGS=-Wl,-D,HEAP=50000+50000\ -Wl,-D,STACK=10000+10000} \
                          ${LIBS=-lsci}  ;;
     alpha-cray-unicosmk*)   : ${OPT=-O2} ${OPT2=-O3\ -hunroll,aggress,restrict=a} \
                          ${OPTAUX=$OPT2} ${LIBS=-lsci} ;;
     alphaev*-cray-unicosmk*) : ${OPT=-O2} ${OPTAUX=-O3\ -hunroll,aggress,restrict=a} \
                          ${OPT2=$OPTAUX,inlinefrom=auxil.c} \
                          ${LIBS=-lsci} ;;
     i?86-pc-linux*|ia64*linux*|x86_64*linux*)
	                  if test $CC = "pgcc" ; then
                           :  ${OPT=-O2} ${OPT2=$OPT} ${OPTAUX=$OPT2 -Munroll}
                          fi
                          if test $CC = "mpicc" ; then
                           :  ${OPT=-O} ${OPT2=$OPT} ${OPTAUX=$OPT2}
                          fi
                          if test $CC = "icc" ; then
                           :  ${OPT=-O3} ${OPT2=$OPT} ${OPTAUX=$OPT2}
                          fi ;;

                    *)  : ${OPT=-O} ${OPT2=$OPT} ${OPTAUX=$OPT2}
    esac
   fi
fi

dnl 
dnl Intel C compiler lies that it's gcc.  Morons
dnl
if test $CC = "icc" ; then   ac_cv_c_compiler_gnu="no"; fi
dnl

dnl Extra compiler options.  We know whether we are using gcc at this point.
if test "x$ac_cv_c_compiler_gnu" = xyes; then
dnl -O2 is added to CFLAGS by AC_PROG_CC so we don't ad it again here.
   case "${target}" in
        i586*|i686*)  OPT=-mcpu=pentium\ -malign-double; OPT2=$OPT;;
                  *)  OPT= ; OPT2=$OPT
   esac
   OPT=${OPT}\ -ffast-math
   OPT2=${OPT2}\ -ffast-math
   OPTAUX=${OPT2}\ -funroll-loops
else
dnl For Non-gcc compilers we turn OFF -g
   CFLAGS=${EXTRA_CFLAGS}
fi
AC_PROG_CPP

AC_PROG_INSTALL
AC_PROG_LN_S
dnl AC_PROG_MAKE_SET

dnl Checks for specific systems
AC_AIX

dnl Checks for libraries.
dnl Replace `main' with a function in -lm:
AC_SEARCH_LIBS(sin, m)
dnl AC_CHECK_LIB(m, main)

dnl Checks for header files.
dnl ANSI headers.  
AC_HEADER_STDC

AC_CHECK_HEADERS(sys/param.h sys/resource.h sys/time.h sys/times.h sys/types.h rpc/xdr.h)
AC_HEADER_TIME

dnl  
dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_CACHE_CHECK(whether $CC supports pointers to void,md_cv_prog_cc_voidptr,
AC_TRY_COMPILE(,void * func(); void *ip; ip=func(ip);, md_cv_prog_cc_voidptr=yes,md_cv_prog_cc_voidptr=no))
if test $md_cv_prog_cc_voidptr = yes; then
   AC_DEFINE(HAVE_VOID_PTRS,,[Define if compiler supports pointers to void])
fi
AC_CACHE_CHECK(whether $CC supports stdarg.h with function prototypes,md_cv_prog_cc_stdarg,
AC_TRY_COMPILE([#include <stdarg.h>], int argtst(int x, ...);, md_cv_prog_cc_stdarg=yes,md_cv_prog_cc_stdarg=no))
if test $md_cv_prog_cc_stdarg = yes; then
   AC_DEFINE(HAVE_STDARG_H,,[Define if compiler supports ANSI stdarg.h functionality])
fi

AC_CHECK_SIZEOF(float, 4)
AC_CHECK_SIZEOF(int, 4)

AC_TYPE_SIZE_T

dnl Checks for library functions.
AC_TYPE_SIGNAL
AC_FUNC_STRFTIME
AC_FUNC_VPRINTF

AC_SEARCH_LIBS(xdr_int, nsl sun rpclib tirpc,
AC_DEFINE(HAVE_XDR_INT,,[Define if sun XDR routines are available]))

AC_CHECK_FUNCS(gettimeofday getrusage strdup strerror strstr strchr memcpy memset remove strtol popen raise times)

AC_REPLACE_FUNCS(getopt)

if test $ac_cv_func_times = "yes" -a $ac_cv_header_sys_times_h = "yes" ; then
AC_CACHE_CHECK(whether times returns the wall-clock time,md_cv_func_times_value,
AC_TRY_RUN([#include <sys/times.h>
#include<time.h>
int main() {struct tms tm; if( times(&tm) <= 0 ) exit(-1); else exit(0);}],
md_cv_func_times_value=yes,md_cv_func_times_value=no,md_cv_func_times_value=yes))

if test $md_cv_func_times_value != yes; then
   AC_DEFINE(TIMES_RETURNS_STATUS,,[Define if the times system call returns a status flag instead of the RT clock])
fi
fi

dnl  Parallel compilation
if test $parlib = "mpi" ; then
dnl Don't check for headers and libraries if compiling with mpicc
dnl Many "mpicc's" are broken and don't support the -E flag needed
dnl to run the preprocessor in the test for headers.
   if test "$CC" = "mpicc" ; then
       ac_cv_header_mpi_h="yes"
       AC_CHECK_FUNC(MPI_Init)
   else
      CF_INCLUDE_PATH(/usr/local/include)
      AC_CHECK_HEADERS(mpi.h)
      AC_SEARCH_LIBS(MPI_Init, $parlib)
      ac_cv_func_MPI_Init="no"
      test "x$ac_cv_search_MPI_Init" = xno || ac_cv_func_MPI_Init="yes"
   fi   
   if test $ac_cv_func_MPI_Init = yes && test $ac_cv_header_mpi_h = yes; then 
      AC_DEFINE(SPMD,,[Define to compile for parallelism]) 
      AC_DEFINE(MPI,,[Define to use MPI parallel interface])
   else
      AC_MSG_WARN(Unable to compile and link with MPI libraries)
    fi
fi
if test $parlib = "bsp" ; then
dnl Don't check for headers  if compiling with mpicc
dnl bspcc v1.3 is broken and doesn't support the -E flag needed
dnl to run the preprocessor in the test for headers.
   if test "$CC" = "bspcc" ; then
       ac_cv_header_bsp_h="yes"
       ac_cv_header_bsp_level1_h="yes"
       AC_CHECK_FUNC(bsp_begin)
   else
      CF_INCLUDE_PATH(/usr/local/include)
      AC_CHECK_HEADERS(bsp.h bsp_level1.h)
      AC_SEARCH_LIBS(bsp_begin, $parlib)
      ac_cv_func_bsp_begin="no"
      test "x$ac_cv_search_bsp_begin" = xno || ac_cv_func_bsp_begin="yes"
   fi
   if  test $ac_cv_func_bsp_begin = "yes" && test $ac_cv_header_bsp_h = "yes" && test $ac_cv_header_bsp_level1_h = "yes"; then
      AC_DEFINE(SPMD) 
      AC_DEFINE(BSP,,[Define to use BSP parallel interface])
   else
      AC_MSG_WARN(Unable to compile and link with BSP libraries)
   fi
fi
if test $parlib = "oldbsp" ; then
   AC_CHECK_LIB(bsp, bspstart)
   if test $ac_cv_lib_bsp_bspstart = "yes"; then
      AC_DEFINE(SPMD) 
      AC_DEFINE(BSP0,,[Define to use old BSP parallel interface])
   else
      AC_MSG_WARN(Unable to compile and link with BSP libraries)
   fi
fi
if test $parlib = "tcgmsg" ; then
   CF_INCLUDE_PATH(/usr/local/include)
   AC_CHECK_HEADERS(sndrcv.h, tcgmsg_sndrcv_h=yes, tcgmsg_sndrcv_h=no)
   AC_SEARCH_LIBS(socket, socket)
   AC_CHECK_LIB($parlib, PBEGIN_, tcgmsg_PBEGIN=yes, tcgmsg_PBEGIN=no)
   if test $tcgmsg_PBEGIN = yes && test $tcgmsg_sndrcv_h=yes ; then 
      AC_DEFINE(SPMD) 
      AC_DEFINE(TCGMSG,,[Define to use TCGMSG parallel interface])
   else
      AC_MSG_WARN(Unable to compile and link with TCGMSG libraries)
   fi
fi
if test $parlib = "shmem" ; then
   AC_CHECK_HEADERS(mpp/shmem.h)
   AC_SEARCH_LIBS(shmem_barrier_all, sma)
   if test "x$ac_cv_search_shmem_barrier_all" != "xno" && test $ac_cv_header_mpp_shmem_h = yes; then 
      AC_DEFINE(SPMD) 
      AC_DEFINE(SHMEM,,[Define to use Cray SHMEM parallel interface])
      AC_CHECK_FUNCS(shmem_long_get globalexit sbreak)
   else
      AC_MSG_WARN(Unable to compile and link with SHMEM libraries)
   fi
fi
dnl Defines
dnl Substitutions
AC_SUBST(CC)
AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(OPT)
AC_SUBST(OPT2)
AC_SUBST(OPTAUX)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)


AC_CONFIG_FILES([Makefile])
AC_CONFIG_COMMANDS([default],[[echo timestamp > stamp-h]],[[]])
AC_OUTPUT
