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
 * $Header: /home/eeyore/keith/md/moldy/RCS/defs.h,v 2.3 93/10/28 10:28:29 keith Stab $
 *
 * $Log:	defs.h,v $
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
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.14  91/03/12  15:43:31  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
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
 * Got rid of typedefs time_t and size_t. 
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
#define          REVISION         "$Revision: 2.3 $"
#define		 REVISION_DATE    "$Date: 93/10/28 10:28:29 $"
#define		 REVISION_STATE   "$State: Stab $"
/******************************************************************************
 *  Configurational information.  Edit this to tailor to your machine	      *
 ******************************************************************************/
/*
 * See if we can detect IBM RS6000. No preprocessor by default, so
 * user MUST user "-DRS6000" on command line".
 * _ALL_SOURCE is necessary to make XDR stuff work.
 */
#ifdef RS6000
#   define __unix__
#   define _ALL_SOURCE
#   define ANSI_LIBS
#endif
#if (defined(__unix__) || defined(__unix)) && !defined(unix)
#   define unix
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
/*
 * Define operating-system dependant default filenames
 */
#ifdef VMS
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
#ifdef __cray__
#   define CRAY
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
#define	MIN3(x,y,z)	MIN(x, MIN(y,z))
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

typedef	double			real;
typedef	int			boolean;
#define	false			0
#define	true			1
typedef	unsigned long int time_mt;/* Larger than any possible time_t */
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

#if (defined(ANSI) || defined(__STDC__)) && !defined(VMS)
typedef void	gptr;
#else
typedef char	gptr;
#endif

#define aalloc(n, type) (type *)talloc((int)(n),sizeof(type),__LINE__, __FILE__)
#define ialloc(n) aalloc(n, int)
#define dalloc(n) aalloc(n, real)
#define ralloc(n) aalloc(n, vec_mt)
#define palloc(n) aalloc(n, vec_mt *)
#define qalloc(n) aalloc(n, quat_mt)

#define	xfree( ptr )	tfree( (gptr *) (ptr) )

#ifdef ANSI_LIBS
#   define memcp(s1,s2,n) (void)memcpy( (gptr*)(s1), (gptr*)(s2), (size_t)(n))
#else
#   define memcp(s1,s2,n) (void)memcpy( (gptr*)(s1), (gptr*)(s2), (int)(n))
#endif
#endif
