/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1988, 1992, 1993 Keith Refson
 
This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 
In other words, you are welcome to use, share and improve this program.
You are forbidden to forbid anyone else to use, share and improve
what you give them.   Help stamp out software-hoarding!  */
/*
 * $Header: /home/eeyore_data/keith/md/moldy/RCS/defs.h,v 2.8 1994/07/07 17:03:39 keith Exp $
 *
 * $Log: defs.h,v $
 * Revision 2.8  1994/07/07  17:03:39  keith
 * Fixed up missing xdr_vector to be compiled in only if NEED_XDR_VECTOR defined.
 *
 * Revision 2.8  1994/07/07  17:03:39  keith
 * Fixed up missing xdr_vector to be compiled in only if NEED_XDR_VECTOR defined.
 *
 * Revision 2.7  1994/06/08  13:10:43  keith
 * New macro "balloc(n,size)" for cases when type isn't explicit.
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Got rid of all global (external) data items except for
 * "control" struct and constant data objects.  The latter
 * (pot_dim, potspec, prog_unit) are declared with CONST
 * qualifier macro which evaluates to "const" or nil
 * depending on ANSI/K&R environment.
 * Also moved as many "write" instantiations of "control"
 * members as possible to "startup", "main" leaving just
 * "dump".
 *
 * Declared as "static"  all functions which should be.
 *
 * Added CONST qualifier to (re-)declarations of ANSI library
 * emulation routines to give reliable compilation even
 * without ANSI_LIBS macro. (#define's away for K&R
 * compilers)
 *
 *
 * Now recognises _UNICOS macro and defines "unix".
 * Specific MSDOS file names added on __MSDOS__ macro.
 * Evaluates CONST macro to const or nil depending on
 * ANSI/KR environment.
 * Added size_mt typedef (ulong) for interfacing with lib fns.
 *
 * Revision 2.5  94/01/25  16:49:41  keith
 * Incorporated all portability experience to multiple platforms since 2.2.
 * Including ports to VAX/VMS and Open VMS on Alpha AXP and Solaris.
 * 
 * Revision 2.4  94/01/18  13:13:42  keith
 * Workaround for bugs and defined symbol _HPUX_SOURCE needed to compile xdr.
 * 
 * Revision 2.3.1.1  93/12/21  19:02:17  keith
 * Mods to allow HP's ANSI compiler to work.  I think it's broken.
 * Not tested on other architectures yet so don't incorporate
 * into main line.
 * NB Added _POSIX_SOURCE and _XOPEN_SOURCE symbols for *all* architectures.
 * This might or might no avoid problems.
 * 
 * Revision 2.3  93/10/28  17:41:01  keith
 * Added __unix to the list of recognised macros
 * 
 * Revision 2.3  93/10/28  10:28:29  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.2  93/10/14  18:17:56  keith
 * Fixed prortability problems to IBM RS6000
 * 
 * Revision 2.1  93/08/18  20:54:03  keith
 * Tidied up clashes over ABS, MIN, MAX macros.
 * 
 * Revision 2.1  93/07/19  13:27:15  keith
 * Added XDR capability for backup and dump files.
 * 
 * Revision 2.0  93/03/15  14:49:28  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.24  93/03/12  12:26:12  keith
 * Reorganized defines to recognise all ANSI (__type__) forms.
 * Fixed up Cray by defining old symbol.
 * 
 * Revision 1.23  93/03/09  15:59:24  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.22  92/06/15  16:52:00  keith
 * Put parens round arg on xfree macro to make it safe with expr.
 * 
 * Revision 1.21  92/06/12  12:55:56  keith
 * Mods to make it work on VMS again.  Ugh.
 * 
 * Revision 1.20  92/06/11  21:40:41  keith
 * Added LOCKEX macro for system-dependent lock extension.
 * 
 * Revision 1.19  92/06/10  15:53:33  keith
 * Added new potential type "generic" for Neal.
 * 
 * Revision 1.18  92/02/26  14:29:04  keith
 * Updated vectorization directive substitution for Convex C vsn 4.3
 * 
 * Revision 1.17  91/08/19  16:49:34  keith
 * Moved #if so that errno.h is included for system V.
 * 
 * Revision 1.16  91/08/17  13:58:59  keith
 * Added "__unix__" symbol for ANSI unix compilers.
 * 
 * Revision 1.15  91/08/15  18:12:22  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_mt and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.14  91/03/12  15:43:31  keith
 * Tidied up typedefs size_mt and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.13  90/09/28  13:29:45  keith
 * Inserted braces around VECTORIZE directives and changed include files
 * for STARDtardent 3000 series (via cond. comp symbol "ardent").
 * 
 * Revision 1.12  90/09/05  10:30:57  keith
 * Support for cray scc added - directives and ANSI_LIBS macro set.
 * 
 * Revision 1.11  90/08/22  10:58:42  keith
 * Changed ANSI libraries conditional compilation to rely on own
 * symbol ANSI_LIBS rather than __STDC__.
 * Corrected test for cray scc compiler in vectorization directives.
 * 
 * Revision 1.10  90/05/02  15:43:52  keith
 * Got rid of typedefs time_t and size_mt. 
 * 
 * Revision 1.9  90/04/12  16:23:44  keith
 * Used <errno.h> to check for Berkeley unix and define symbol BSD
 * 
 * Revision 1.8  90/04/06  11:09:22  keith
 * Aquired definition of NPOTP from structs.
 * 
 * Revision 1.7  90/03/27  17:36:12  keith
 * Moved O/S dependent conditionals to here, esp VPRINTF.
 * Reorganised configuration conditionals into one block.
 * 
 * Revision 1.6  90/03/26  18:04:33  keith
 * Tidied up system dependant includes.
 * Added system-dependant backup and temp file names (for input.c).
 * 
 * Revision 1.5  90/03/09  17:34:45  keith
 * Added preprocessor directives to define USG (ie system V) for unicos.
 * 
 * Revision 1.4  89/09/04  18:41:49  keith
 * Added conversion constants for charges.
 * 
 * Revision 1.3  89/06/14  14:16:35  keith
 * Added vectorisation for stellar and recognised sysV & SysV macros
 * 
 * Revision 1.2  89/05/22  14:05:51  keith
 * Added rescale-separately option, changed 'contr_t' format.
 * 
 * Revision 1.1  89/05/02  10:51:58  keith
 * Initial revision
 * 
 * 
 */

#ifndef	DEFS_ALREADY
#define DEFS_ALREADY

/*
 * Version ID strings
 */
#define          REVISION         "$Revision: 2.8 $"
#define		 REVISION_DATE    "$Date: 1994/07/07 17:03:39 $"
#define		 REVISION_STATE   "$State: Exp $"
/******************************************************************************
 *  Configurational information.  Edit this to tailor to your machine	      *
 ******************************************************************************/
/*
 * See if we can detect IBM RS6000. No preprocessor macro by default, so
 * user MUST define "-DRS6000" on command line".
 * _ALL_SOURCE is necessary to make XDR stuff work.
 */
#ifdef RS6000
#   define __unix__
#   define _ALL_SOURCE
#   define ANSI_LIBS
#endif
/*
 * To allow XDR stuff to work on HP.  Surely there's a more general
 * way of doing this? I think that HPs header files are broken.
 */
#ifdef __hpux
#define _HPUX_SOURCE
#endif

#if (defined(__unix__) || defined(__unix) || defined(_unix_) || defined(_unix) || defined(_UNICOS)) && !defined(unix)
#   define unix
#endif
#ifdef __vms
#   ifndef vms
#      define vms
#   endif
#   ifndef VMS
#      define VMS
#   endif
#endif
#ifdef ardent
#   define ARDENT
#   define NEED_XDR_VECTOR
#endif
/*
 * Define operating-system dependant default filenames
 */
#ifdef __MSDOS__
#   define BACKUP_FILE	"MDBCK"
#   define TEMP_FILE	"MDTEMPX"
#   define LOCKEX		"$LK"
#endif
#ifdef vms
#   define BACKUP_FILE	"MDBACKUP.DAT"
#   define TEMP_FILE	"MDTEMPXXXX.DAT"
#   define LOCKEX		"$LCK"
#endif
#ifdef CMS
#   define BACKUP_FILE	"MDBACKUP MOLDY A1"
#   define TEMP_FILE	"MDTEMP XXXXXXXX A1"
#endif
#ifdef unix
#   define LOCKEX		".lck"
#endif
#ifndef LOCKEX
#   define LOCKEX		"LK"
#endif
/*
 * Set ANSI_LIBS only if you have the standard ANSI headers and libraries
 */
#ifdef _CRAY
#ifndef CRAY
#   define CRAY
#endif
#endif
#if defined(CRAY) && defined(__STDC__)	/* scc compiler comes with libraries*/
#   define ANSI_LIBS
#endif
/*
 * New convex compiler is ANSI, but doesn't define __STDC__ or  convexvc.
 * It has a silly macro __stdc__ which we will use instead.  convexvc
 * may not be necessary as <fastmath.h> is no longer required.  But there
 * is still veclib.
 */
#if defined(__convexc__)
#   if defined(__stdc__) /* Anything but "-pcc" mode */
#      define ANSI
#      define ANSI_LIBS
#   endif
#endif
#if defined(vms)
#   define ANSI
#   define ANSI_LIBS
#endif
/*
 * Set HAVE_VPRINTF if this function is in target machine's library.
 */
#ifdef ANSI_LIBS			/* ANSI has it			*/
#define HAVE_VPRINTF
#endif
#if defined(sun) || defined(stellar) || defined(titan) /* So do these     */
#define HAVE_VPRINTF
#endif
#if defined(cray) && (defined(unix) || defined(__unix__)) /* ie UNICOS   	*/
#define HAVE_VPRINTF
#endif
/*
 *  Otherwise does machine have _doprnt?
 */
#if defined(convex) || defined(sequent)
#define HAVE_DOPRNT
#endif
/*
 *  Get rid of "const" keyword for non-ANSI compilers.
 */
#if defined(__STDC__) ||  defined(ANSI)
#   define CONST const
#else
#   define CONST /* */
#endif
/* 
 * Vectorisation directive translation.  N.B. Most preprocessors munge 
 * directives so the #define must substiture the preprocessor OUTPUT.
 */
#ifdef CRAY
#   ifdef __STDC__
#      define VECTORIZE
#      define NOVECTOR
#   else
#      define VECTORIZE ## ivdep
#      define NOVECTOR  ## novector
#   endif
#else
#ifdef __convexc__
#   define VECTORIZE ;/* $dir no_recurrence */;
#   define NOVECTOR  ;/* $dir scalar */;
#else
#ifdef stellar
#   define VECTORIZE __dir NO_RECURRENCE :
#   define NOVECTOR  __dir SCALAR :
#else
#ifdef ardent
#   define VECTORIZE # pragma ivdep
#   define NOVECTOR  # pragma novector
#else
#   define VECTORIZE /* Canny  vectorise on this machine!*/
#   define NOVECTOR  /* */
#endif
#endif
#endif
#endif
/*
 *  Set symbol USG to identify system V variant of unix, BSD for Berkeley.
 */
#include <errno.h>
/*
 * Berkeley error numbers appear to be very regular, so the following is
 * pretty likely to get it right.  If it doesn't, do the define by hand.
 * BSD is used to signal use of getrusage() rather than times() (unix only)
 * and absence of mem*() and strchr() functions from library.
 */
#if defined(unix) && !defined(USG)
#if EWOULDBLOCK==35 && EINPROGRESS==36 && EALREADY==37
# define BSD
#else
# define USG
#endif
#endif
/******************************************************************************
 *  End of machine/OS configuration.					      *
 ******************************************************************************/
#define NPOTP 7                 /* Must be number of doubles in pot_mt       */

#ifndef SEEK_END
#define	SEEK_SET	0
#define	SEEK_CUR	1
#define SEEK_END	2
#endif

#define	L_name		128			/* Max Length of file names  */
#define	NPE		2			/* real & Ewald PE's	     */
#undef  MAX
#define MAX(x,y)	((x) > (y) ? (x) : (y))
#define MAX3(x,y,z)	MAX(x, MAX(y,z))
#undef  MIN
#define MIN(x,y)	((x) < (y) ? (x) : (y))
#define	MIN3(x,y,z)	MIN(((x) < (y) ? (x) : (y)),z)
#define SUMSQ(x)	(x[0]*x[0] + x[1]*x[1] +x[2]*x[2])
#define SUMSQ2(x)	(x[1]*x[1] + x[2]*x[2] +x[3]*x[3])
#define	SQR(x)		((x) * (x))
#define CUBE(x)		((x) * (x) * (x))
#define DISTANCE(x,y)	sqrt(SQR(x[0]-y[0])+SQR(x[1]-y[1])+SQR(x[2]-y[2]))
#define	SIGN(x, y)	((y) > 0 ? fabs(x) : -fabs(x))
#define ALPHAMIN	1e-7
/* Fundamental constants '_' denotes MKS units.  SHOULD NEVER BE ALTERED      */
#define	PI		3.14159265358979323846
#define DTOR		(PI/180.0)
#define	_kB		1.380662e-23
#define _kcal		4184.0
#define	_EPS0		8.85418782e-12
#define	_ELCHG		1.6021892e-19
#define	_ROOT_4_PI_EPS	1.0548222865e-5
#define AVOGAD		6.022045e23
/* Program units relative to MKS units ONLY M, L, T, Q UNIT SHOULD BE CHANGED */
#define MUNIT		1.6605655e-27			/* atomic mass units  */
#define	LUNIT		1.0e-10				/* Angstrom	      */
#define	TUNIT		1.0e-12				/* Picosecond	      */
#define	EUNIT		(MUNIT * (LUNIT/TUNIT) * (LUNIT/TUNIT))
#define	QUNIT		/*sqrt(EUNIT*LUNIT)*/ (4.075e-17*_ROOT_4_PI_EPS)
#define MUNIT_N		"amu"
#define LUNIT_N		"A"
#define RLUNIT_N	"A(-1)"
#define	TUNIT_N		"ps"
#define	IUNIT_N		"amuA**2"
#define DIPUNIT_N	"Prog units"
/* Constants in program units						     */
#define	kB		(_kB / EUNIT)
#define EPS0		(0.25 / PI)
/* Output conversion units */
#define	CONV_E		(0.001*AVOGAD*EUNIT)			/* kJ/mol     */
#define	CONV_T		(kB / CONV_E)				/* kB / E unit*/
#define	CONV_P		(MUNIT/(LUNIT*TUNIT*TUNIT)/1.0e6) 	/* MPa	      */
#define	CONV_V		CONV_E
#define	CONV_F		(MUNIT*LUNIT/(TUNIT*TUNIT)*sqrt(AVOGAD))/* N/mol      */
#define	CONV_N		(CONV_F*LUNIT)				/* Nm/mol     */
#define	CONV_D		(QUNIT*LUNIT*4.8e10/_ELCHG)
#define CONV_Q		(QUNIT/_ELCHG)
#define	CONV_E_N	"kJ/mol"
#define	CONV_T_N	"K"
#define	CONV_P_N	"Mpa"
#define	CONV_F_N	"N**2/mol"
#define	CONV_N_N	"(Nm)**2/mol"
#define	CONV_D_N	"D"
#define CONV_Q_N	"Qe"


#define	false			0
#define	true			1
typedef	unsigned long int time_mt;/* Larger than any possible time_t */
typedef unsigned long size_mt;    /* Wide type for passing sizeof	      */

#if (defined(ANSI) || defined(__STDC__)) && !defined(VMS)
typedef void	gptr;
#else
typedef char	gptr;
#endif

typedef	double			real;
typedef	int			boolean;
typedef enum {inv, noinv} 	invrot;
typedef enum {tke_n, rke_n, pe_n, e_n, tt_n, rt_n, t_n, 
              h0_n, h1_n, h2_n, stress0_n, stress1_n,  stress2_n, press_n, 
              vir_n, msqf_n, msqt_n, dip_n, end} av_n;
           
typedef real    vec_mt[3];
typedef	vec_mt	*vec_mp;
typedef	real	quat_mt[4];
typedef quat_mt	*quat_mp;
typedef real    mat_mt[3][3];
typedef vec_mt	*mat_mp;

#define balloc(n,size)  talloc((int)(n),(size_mt)(size), __LINE__, __FILE__)
#define aalloc(n, type) (type *)balloc((n), sizeof(type))
				       
#define ialloc(n) aalloc(n, int)
#define dalloc(n) aalloc(n, real)
#define ralloc(n) aalloc(n, vec_mt)
#define palloc(n) aalloc(n, vec_mt *)
#define qalloc(n) aalloc(n, quat_mt)

#define	xfree( ptr )	tfree( (gptr *) (ptr) )
#define lsizeof         (size_mt)sizeof

#ifdef ANSI_LIBS
#   define memcp(s1,s2,n) (void)memcpy( (gptr*)(s1), (gptr*)(s2), (size_mt)(n))
#   define memst(s,c,n)   (void)memset( (gptr*)(s), (c), (size_mt)(n))
#else
#   define memcp(s1,s2,n) (void)memcpy( (gptr*)(s1), (gptr*)(s2), (int)(n))
#   define memst(s,c,n)   (void)memset( (gptr*)(s), (c), (int)(n))
#endif
#endif
