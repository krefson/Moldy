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

#ifndef	DEFS_ALREADY
#define DEFS_ALREADY

/*
 * Version ID strings
 */
#include "version.h"

#ifdef HAVE_CONFIG_H
#   include "config.h"
#else
#   include "defconf.h"
#endif
/*
 * Post-config.h processing
 */
#ifdef _AIX
#   define __unix__
#endif

/*
 * To allow XDR stuff to work on HP.  Surely there's a more general
 * way of doing this? I think that HPs header files are broken.
 */
#if defined(__hpux) && ! defined(__convex_spp)
#   define _HPUX_SOURCE
#   define __unix__
#endif

#ifdef __convexc__              /* For"-std" and "-str" ANSI modes. Sigh. */
#   ifndef _POSIX_SOURCE
#      define _POSIX_SOURCE
#      define _CONVEX_SOURCE
#   endif
#endif

#if defined(HAVE_XDR_INT) && defined(HAVE_RPC_XDR_H)
#   define USE_XDR
#endif

/* 
 * Vectorisation directive translation.  N.B. Most preprocessors munge 
 * directives so the #define must substitute the preprocessor OUTPUT.
 */
#if defined(_CRAY1)	/* This seems to be defined for all CRI PVP systems */
#   define VECTOR    /* To choose vector vsn of site_neighbour_list */
#endif
#if defined(__convexc__) && ! defined(__convex_spp)
#   define VECTOR    /* To choose vector vsn of site_neighbour_list */
#endif
#ifdef stellar
#   define VECTORIZE __dir NO_RECURRENCE :
#   define NOVECTOR  __dir SCALAR :
#   define VECTOR    /* To choose vector vsn of site_neighbour_list */
#endif
#ifdef ardent
#   define VECTORIZE # pragma ivdep
#   define NOVECTOR  # pragma novector
#   define VECTOR    /* To choose vector vsn of site_neighbour_list */
#endif

#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__CYGWIN__)
#   define UNRESTRICTED_FILE_NAMES
#   define BACKUP_FILE	"MDBACKUP"
#   define LOCKEX		".lck"
#   define ALLOC_ALIGN
#endif

#include <errno.h>

#ifndef BACKUP_FILE
#   define BACKUP_FILE	"MDBCK"
#endif
#ifdef TEMP_FILE
#   define TEMP_FILE	"MDTEMPX"
#endif
#ifndef LOCKEX
#   define LOCKEX	"LK"
#endif

/*
 * Hardware configuration.  Specify cache-tuning params NCACHE and NLINE
 * NCACHE is minimum size of "sites" arrays and should be a sub-multiple
 *   of the cache length in WORDS (sizeof(real)).  It MUST be a power
 *   of two since the alignment algorithm relies on this.
 * NLINE is padding (in words) between arrays to avoid cache conflicts.
 * Theoretically NCACHE should equal the size of the cache and NLINE
 * the size of the cache line.  But too large an NCACHE wastes memory
 * (arrays are padded to a multiple of this) and is actually slower
 * for small systems presumably because it causes TLB misses. 512
 * works well for HP/PA and IBM RISC/POWER.
 * NLINE ought really to be a multiple of the cache line size but it
 *   doesn't seem to matter as long as it's not zero.
 * 512/4 works well on HP/PA/IBM RISC/POWER/SPARC/ALPHA/MIPS so use as
 * default.  Hardware-specific values could be set here conditionally.
 */
#ifndef NLINE
#   ifdef _CRAYT3E
#      define NLINE 8
#   endif
#   if __mips == 4 
#      define NLINE 16
#   endif
#endif
#ifndef NCACHE
#   define NCACHE 256
#endif
#ifndef NLINE
#   define NLINE 4
#endif
/* 
 * Vectorisation directive translation.  N.B. Most preprocessors munge 
 * directives so the #define must substitute the preprocessor OUTPUT.
 */
#ifndef VECTORIZE
#   define VECTORIZE
#   define NOVECTOR
#endif
/******************************************************************************
 *  End of machine/OS configuration.					      *
 ******************************************************************************/
#ifndef NPOTP
#define NPOTP 8                 /* Must be number of doubles in pot_mt       */
#endif

#ifndef SEEK_END
#define	SEEK_SET	0
#define	SEEK_CUR	1
#define SEEK_END	2
#endif

#define MOLPBC 1
#define SITEPBC 0
/*
 * N.B.  If these are changed the restart and dump formats will be incompatible.
 */
#define         SFORM   "%127[^#]" /* Format for scanf to read strings safely */
#define	L_name		128			/* Max Length of file names  */
#define L_spec           32			/* Max length of species name*/
#define L_site            8			/* Max length of site name   */
#define L_vsn            16			/* Max length of RCS version */
#define DLEN             28			/* Length of date/time string*/
#define	NPE		2			/* real & Ewald PE's	     */
#ifdef MAX
#undef  MAX
#endif
#define MAX(x,y)	((x) > (y) ? (x) : (y))
#define MAX3(x,y,z)	MAX(x, MAX(y,z))
#ifdef MIN
#undef  MIN
#endif
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
#define PISQ		(PI*PI)
#define ROOTPI		1.77245385090551602729
#define PI32		(PI*ROOTPI)
#define DTOR		(PI/180.0)
#ifdef OLDCONSTS
#define	_kB		1.380662e-23
#define _kcal		4184.0
#define	_EPS0		8.85418782e-12
#define	_ELCHG		1.6021892e-19
#define ECSTR           "1.6021892e-19"
#define	_ROOT_4_PI_EPS	1.0548222865e-5
#define AVOGAD		6.022045e23
#define AMU		1.6605655e-27
#define AMUSTR          "1.6605655e-27"
#define RTAMU           4.0745e-14		/*sqrt(amu) */
#else
/* Values from CODATA 1986 */
#define	_kB		1.380658e-23
#define _kcal		4184.0
#define	_EPS0		8.854187817e-12
#define	_ELCHG		1.60217733e-19
#define ECSTR           "1.60217733e-19"
#define	_ROOT_4_PI_EPS	1.05482230112e-05
#define AVOGAD		6.0221367e23
#define AMU		1.6605402e-27
#define AMUSTR          "1.6605402e-27"
#define RTAMU           4.07497263794495e-14	/*sqrt(amu) */
#endif
#define RGAS		(AVOGAD*_kB)
/* Program units relative to MKS units ONLY M, L, T, Q UNIT SHOULD BE CHANGED */
#define MUNIT		AMU				/* atomic mass units  */
#define	LUNIT		1.0e-10				/* Angstrom	      */
#define	TUNIT		1.0e-12				/* Picosecond	      */
#define	EUNIT		(MUNIT * (LUNIT/TUNIT) * (LUNIT/TUNIT))
#define	QUNIT		/*sqrt(EUNIT*LUNIT)*/ (RTAMU*1.e-3*_ROOT_4_PI_EPS)
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
#define CONV_TM		0.01
#define	CONV_E_N	"kJ/mol"
#define	CONV_T_N	"K"
#define	CONV_P_N	"MPa"
#define	CONV_F_N	"N**2/mol"
#define	CONV_N_N	"(Nm)**2/mol"
#define	CONV_D_N	"D"
#define CONV_Q_N	"Qe"
#define CONV_TM_N	"kJ/mol*ps**2"


#define	false			0
#define	true			1
typedef	unsigned long int time_mt;/* Larger than any possible time_t */
typedef unsigned long size_mt;    /* Wide type for passing sizeof	      */

#ifdef HAVE_VOID_PTRS
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

#ifdef STDC_HEADERS
#   define memcp(s1,s2,n) (void)memcpy( (gptr*)(s1), (gptr*)(s2), (size_mt)(n))
#   define memst(s,c,n)   (void)memset( (gptr*)(s), (c), (size_mt)(n))
#else
#   define memcp(s1,s2,n) (void)memcpy( (gptr*)(s1), (gptr*)(s2), (int)(n))
#   define memst(s,c,n)   (void)memset( (gptr*)(s), (c), (int)(n))
#endif
#endif
/*========================== Dump ==========================================*/
#define DUMP_SIZE(level, n, n_r)  \
   			 (( (level & 1)    * (3*n + 4*n_r + 9 + 2)+ \
			    (level>>1 & 1) * (3*n + 3*n_r + 9 + 1)+ \
                            (level>>3 & 1) * (3*n + 3*n_r + 9) ))
#define INERTIA_MIN  1.0e-14       
