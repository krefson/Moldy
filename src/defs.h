/*
 * $Header: /usr/data/keith/moldy/RCS/defs.h,v 1.11 90/08/10 17:47:29 keith Exp $
 *
 * $Log:	defs.h,v $
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
#define          REVISION         "$Revision: 1.11 $"
#define		 REVISION_DATE    "$Date: 90/08/10 17:47:29 $"
#define		 REVISION_STATE   "$State: Exp $"
/******************************************************************************
 *  Configurational information.  Edit this to tailor to your machine	      *
 ******************************************************************************/
/*
 *  Set symbol USG to identify system V variant of unix, BSD for Berkeley.
 */
#if defined(unix) && !defined(NOBSD)
#include <errno.h>
/*
 * Berkeley error numbers appear to be very regular, so the following is
 * pretty likely to get it right.  If it doesn't, do the define by hand.
 * BSD is used to signal use of getrusage() rather than times() (unix only)
 * and absence of mem*() and strchr() functions from library.
 */
#if EWOULDBLOCK==35 && EINPROGRESS==36 && EALREADY==37
# define BSD
#endif
#endif
/*
 * Define operating-system dependant default filenames
 */
#ifdef VMS
#define BACKUP_FILE	"MDBACKUP.DAT"
#define TEMP_FILE	"MDTEMPXXXX.DAT"
#endif
#ifdef CMS
#define BACKUP_FILE	"MDBACKUP MOLDY A1"
#define TEMP_FILE	"MDTEMP XXXXXXXX A1"
#endif
/*
 * Set ANSI_LIBS only if you have the standard ANSI headers and libraries
 */
#undef ANSI_LIBS
/*
 * Set HAVE_VPRINTF if this function is in target machine's library.
 */
#ifdef ANSI_LIBS			/* ANSI has it			*/
#define HAVE_VPRINTF
#endif
#if defined(sun) || defined(vms) || defined(stellar) /* So do these     */
#define HAVE_VPRINTF
#endif
#if defined(cray) && defined(unix)	/* ie UNICOS			*/
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
#if defined(CRAY) && !defined(__STDC__)
#define VECTORIZE ## ivdep
#define NOVECTOR  ## novector
#else
#ifdef convexvc
#define VECTORIZE __dir no_recurrence :
#define NOVECTOR  __dir scalar :
#else
#ifdef stellar
#define VECTORIZE __dir NO_RECURRENCE :
#define NOVECTOR  __dir SCALAR :
#else
#define VECTORIZE /* Canny  vectorise on this machine!*/
#define NOVECTOR  /* */
#endif
#endif
#endif
/******************************************************************************
 *  End of machine/OS configuration.					      *
 ******************************************************************************/
#define NPOTP 5                 /* Must be number of doubles in pot_t        */

#ifndef SEEK_END
#define	SEEK_SET	0
#define	SEEK_CUR	1
#define SEEK_END	2
#endif

#define	L_name		128			/* Max Length of file names  */
#define	NPE		2			/* real & Ewald PE's	     */
#undef	NULL
#define NULL		0
#define ABS(x)		((x) > 0 ? (x) : -(x))
#define MAX(x,y)	((x) > (y) ? (x) : (y))
#define MAX3(x,y,z)	MAX(x, MAX(y,z))
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
typedef	unsigned long int mytime_t;	/* Larger than any possible time_t */
typedef enum {inv, noinv} 	invrot;
typedef enum {tke_n, rke_n, pe_n, e_n, tt_n, rt_n, t_n, 
              h0_n, h1_n, h2_n, stress0_n, stress1_n,  stress2_n, press_n, 
              vir_n, msqf_n, msqt_n, dip_n, end} av_n;
           
typedef real    vec_t[3];
typedef	vec_t	*vec_p;
typedef	real	quat_t[4];
typedef quat_t	*quat_p;
typedef real    mat_t[3][3];
typedef vec_t	*mat_p;

char	*talloc();
#define aalloc(n, type) (type *)talloc((long)n,sizeof(type),__LINE__, __FILE__)
#define ialloc(n) aalloc(n, int)
#define dalloc(n) aalloc(n, real)
#define ralloc(n) aalloc(n, vec_t)
#define palloc(n) aalloc(n, vec_t *)
#define qalloc(n) aalloc(n, quat_t)

#endif
