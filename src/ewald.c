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
/******************************************************************************
 * Ewald	The reciprocal-space part of the standard Ewald sum technique *
 ******************************************************************************
 *      Revision Log
 *       $Log: ewald.c,v $
 *       Revision 2.23.4.1  2002/03/13 10:27:52  kr
 *       Trial version incorporating reciprocal-space summation for r^-2 and r^-6
 *       interactions.  This version implements a new potential "genpot46" to activate.
 *
 *       Revision 2.23  2001/03/02 11:43:30  keith
 *       Corrected fix for cache re-use bug.  The fix meant cache was never re-used!
 *
 *       Revision 2.22  2001/02/15 15:55:09  keith
 *       Fixed bug where qsincos  could incorrectly attempt to use
 *       cache coshxky etc on first iteration when unassigned.
 *
 *       Revision 2.21  2000/12/06 17:45:29  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.20  2000/04/27 17:57:07  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.19  1998/12/03 15:45:17  keith
 *       Hand unrolled stress loops to avoid attempts by compilers to optimise
 *       loops with 1 or 2 iterations.
 *
 *       Revision 2.18  1998/11/26 17:08:02  keith
 *       Performance improvements.
 *        a) Cache values of sin and  cos(hx+ky) between calls of qsincos().
 *        b) Eliminate use of chg[] array in inner loops by multiplying it
 *           into sin(lz) and cos(lz). Also avoids need for separae, cache
 *           -alighed copy.
 *
 *       Revision 2.17  1998/11/25 14:44:32  keith
 *       Modified loops in "qsincos" adding explicit scalar temporaries.
 *       This gives a speedup of 20-40% on a T3E when used with the "-hsplit"
 *       compiler option as it prevents thrashing the stream buffers.
 *
 *       Revision 2.16  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 *
 *       Revision 2.15  1998/01/27 15:46:45  keith
 *       Got rid of __MSDOS__ and "unix" macros.
 *       Replaced with ALLOC_SEPARATELY and ALLOC_ALIGN.
 *
 *       Revision 2.14  1997/10/06 09:01:27  keith
 *       Changed text of self-energy message to read kJ/mol.
 *
 *       Revision 2.13  1996/11/05 16:45:49  keith
 *       - Reorganized code by extracting generation of sin and cosine tables
 *         into new function trig_recur()
 *       - Optimized main loops and trig_recur to avoid cache conflicts.  Small
 *         gain on the T3D, HUGE gain (x4) on the IBM RS6000.
 *       - New function allocate_arrays() sets out arrays in memory.
 *       _ Preprocessor macros NCACHE and NLINE are used to tune this. NCACHE
 *         is a (possibly sub-multiple) of the cache size, and array dimensions
 *         are rounded up to this. NLINE is larger than the cache-line size and
 *         this space is inserted between arrays as padding.  Both are
 *         specified in WORDS (sizeof(real)).  Also used in accel.c for
 *         declaration of site force arrays.
 *
 *       Revision 2.12  1996/03/19 12:27:48  keith
 *       Parallelized trig function generation loops.  This is conditional
 *       on macro MPPMANY and is only a gain on the T3D. Requires
 *       par_collect_all() in parallel.c.
 *
 *       Revision 2.11  1996/02/07 18:26:58  keith
 *       Restructured for convergence of std and RIL versions
 *        - eliminated ewald-inner.
 *
 *       Revision 2.10  1994/12/30 11:46:08  keith
 *       Fixed bug which caused core dump for very small k-cutoff (hmax=0)
 *
 * Revision 2.9  1994/07/07  16:58:14  keith
 * Versions for SPMD parallel machines.
 * These are based on the *_parallel.c versions for shared-memory but
 * with the parallel loops and compiler directives deleted.
 *
 * Ewald.c has some more optimization of the serial bits.
 *
 * Revision 2.7  1994/06/08  13:13:59  keith
 * New version of array allocator which breaks up requests for DOS.
 * Now must use specific "afree()" paired with arralloc().
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Got rid of all global (external) data items except for
 * "control" struct and constant data objects.  The latter
 * (pot_dim, potspec, prog_unit) are declared with const
 * qualifier macro which evaluates to "const" or nil
 * depending on ANSI/K&R environment.
 * Also moved as many "write" instantiations of "control"
 * members as possible to "startup", "main" leaving just
 * "dump".
 *
 * Declared as "static"  all functions which should be.
 *
 * Added const qualifier to (re-)declarations of ANSI library
 * emulation routines to give reliable compilation even
 * without ANSI_LIBS macro. (#define's away for K&R
 * compilers)
 *
 * Revision 2.5  1994/01/26  16:34:36  keith
 * Fixed non-ansi #endif.
 *
 * Revision 2.3  93/10/28  10:28:59  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.1  93/09/02  12:34:44  keith
 * Optimized qsincos() -- should give up to 25% speed improvement on
 * compilers without assert no aliasing options.
 * 
 * Revision 2.0  93/03/15  14:49:49  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.22  93/03/12  12:23:05  keith
 * Reorganized defines to recognise all ANSI (__type__) forms.
 * Moved spxpy() from aux.c to force.c and force_parallel.c
 * 
 * 
 * Revision 1.21  93/03/09  15:59:58  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.20  93/03/05  15:01:53  keith
 * Added CRAY parallelising directives
 * 
 * Revision 1.19  92/08/13  17:56:58  keith
 * Modified nprocessors to limit execution to 1 proc
 * unless env var THREADS explicitly set.
 * 
 * Revision 1.18  92/06/26  17:03:02  keith
 * Got rid of assumption that memory returned by talloc() or
 * arralloc() is zeroed.  This enhances ANSI compatibility.
 * Removed memory zeroing from alloc.c() in consequence.
 * 
 * Revision 1.17  92/02/26  14:33:48  keith
 * Got rid of pstrip pragmas for convex -- they just broke things
 * 
 * Revision 1.16  91/11/27  15:15:56  keith
 * Corrected calculation of sheet energy term for charged framework.
 * Split force loop so as to omit frame-frame force (and stress) terms.
 * Added spaces so that Stellix compiler doesn't choke on "=*"
 * Replaced main loop variable with pointer to work round Stellix 2.3 bug
 * 
 * Revision 1.16  91/11/26  12:48:43  keith
 * Put #ifdefs around machine-specific pragmas for portability.  Now
 * supports Stardent 1000,2000 & 3000 (titan) series and convex.
 * 
 * Revision 1.13  91/08/24  16:55:18  keith
 * Added pragmas for convex C240 parallelization
 * 
 * Revision 1.12  91/08/15  18:13:17  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_mt and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.11  91/05/29  17:02:00  keith
 * Modified for minor speed improvement on Stardent titan
 * 
 * Revision 1.10  91/03/12  16:30:08  keith
 * Stardent Titan (ST3000) version.
 * 
 * Revision 1.9  90/09/28  13:29:19  keith
 * Inserted braces around VECTORIZE directives and changed include files
 * for STARDtardent 3000 series (via cond. comp symbol "ardent").
 * 
 * Revision 1.8  90/08/29  11:00:49  keith
 * Speeded up loop at 231 to improve parallel efficiency.
 * 
 * Revision 1.7  90/08/01  19:11:40  keith
 * Modified to exclude framework-framework interactions.
 * N.B. Excluded from pe and stress but NOT forces (as they sum to 0).
 * 
 * Revision 1.6  90/05/16  18:40:57  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.5  90/05/16  14:20:47  keith
 * *** empty log message ***
 * 
 * Revision 1.4  90/05/02  15:37:31  keith
 * Removed references to size_mt and time_t typedefs, no longer in "defs.h"
 * 
 * Revision 1.3  90/04/26  15:29:48  keith
 * Changed declaration of arralloc back to char*
 * 
 * Revision 1.2  90/04/25  10:37:06  keith
 * Fixed bug which led to ihkl[] accessing outside bounds of chx[] etc arrays
 * when in const-pressure mode and MD cell expanded between steps.
 * 
 * Revision 1.1  90/01/31  13:18:59  keith
 * Initial revision
 * 
 * Revision 1.6  89/12/15  12:56:26  keith
 * Added conditional ionclusion of <fastmath.h> for stellar
 * 
 * Revision 1.5  89/11/01  17:29:10  keith
 * Sin and cos loop vectorised - mat_vec_mul extracted from loop.
 * 'Uniform charge sheet' term added in case of electrically charged system.
 * 
 * Revision 1.5  89/10/26  11:29:26  keith
 * Sin and cos loop vectorised - mat_vec_mul extracted from loop.
 * 'Uniform charge sheet' term added in case of electrically charged layer.
 * 
 * Revision 1.5  89/10/26  11:27:31  keith
 * Sin and cos loop vectorised - mat_vec_mul extracted from loop.
 * 
 * Revision 1.4  89/10/02  11:39:14  keith
 * Fixed error in *star macros which assumed rlv's were cols of h(-1) r.t. rows.
 * 
 * Revision 1.3  89/06/09  13:38:17  keith
 * Older code for computation of q cos/sin k.r restored conditionally by
 * use of macro OLDEWALD.  This is for more primitive vectorising compilers.
 * 
 * Revision 1.2  89/06/08  10:42:51  keith
 * Modified to circumvent compiler bug in VMS/VAXC 2.4-026
 * 
 * Revision 1.1  89/04/20  16:00:39  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/kr/CVS/moldy/src/ewald.c,v 2.23.4.1 2002/03/13 10:27:52 kr Exp $";
#endif
/*========================== Program include files ===========================*/
#include 	"defs.h"
/*========================== Library include files ===========================*/
#ifdef stellar
#   include 	<fastmath.h>
#else
#ifdef titan
#   include 	<vmath.h>
#else
#   include 	<math.h>
#endif
#endif
#include	"stddef.h"
#include 	"stdlib.h"
#include        "string.h"
/*========================== Program include files ===========================*/
#include 	"structs.h"
#include 	"messages.h"
/*========================== External function declarations ==================*/
gptr            *talloc(int n, size_mt size, int line, char *file);
                                       /* Interface to memory allocator       */
void            tfree(gptr *p);	       /* Free allocated memory	      	      */
void            afree(gptr *pp);       /* Free allocated array	      	      */
double	err_fn(double x);              /* Error function		      */
double	det(real (*a)[3]);	       /* Determinant of 3x3 matrix	      */
void	invert(real (*a)[3], real (*b)[3]);/* Inverts a 3x3 matrix	      */
void	mat_vec_mul(real (*m)[3], vec_mp in_vec, vec_mp out_vec, int number);
                                        /* Multiplies a 3x3 matrix by 3xN vect*/
void	mat_sca_mul(register real s, real (*a)[3], real (*b)[3]);
					/* Multiplies a 3x3 matrix by scalar  */
void	transpose(real (*a)[3], real (*b)[3]); /* Transposes a 3x3 matrix     */
void    zero_real(real *r, int n);      /* Initialiser                        */
void    zero_double(double *r, int n); 	/* Initialiser                        */
double	sum(register int n, register double *x, register int ix);
					/* Sum of elements of 'real' vector   */
gptr	*arralloc(size_mt,int,...); 	/* Array allocator		      */
void	note(char *,...);		/* Write a message to the output file */
void	message(int *,...);		/* Write a warning or error message   */
/*========================== External data references ========================*/
extern	contr_mt	control;       	/* Main simulation control record     */
extern int		ithread, nthreads;
/*========================== Macros ==========================================*/
#define astar hinvp[0]
#define bstar hinvp[1]
#define cstar hinvp[2]
#define moda(hmat) (hmat[0][0])
#define modb(hmat) sqrt(SQR(hmat[0][1]) + SQR(hmat[1][1]))
#define modc(hmat) sqrt(SQR(hmat[0][2]) + SQR(hmat[1][2]) + SQR(hmat[2][2]))

#define FOURTH(x) SQR(SQR(x))
#define SIXTH(x) CUBE(SQR(x))
#define EIGTH(x) SQR(FOURTH(x))
#define POTZERO_TOL 1.0e-16
/*========================== Cache Parameters=================================*/
/* The default values are for the Cray T3D but are probably good enough 
 *  for most other systems too. */
#ifndef NCACHE
#   define NCACHE (256*sizeof(double)/sizeof(real))
#endif
#ifndef NLINE
#   define NLINE  (4*sizeof(double)/sizeof(real))
#endif
/*============================================================================*/
   struct s_hkl {double kx, ky, kz; int h,k,l;};
/*****************************************************************************
 * qsincos().  Evaluate q sin(k.r) and q cos(k.r).  This is in a separate    *
 * function because some compilers (notably Stellar's) generate MUCH better  *
 * vector code this way. 						     *
 * It is designed to be called with an arbitrary combination h,k,l, so can   *
 * not rely on the order of the h,k,l loops.  Instead it caches the partial  *
 * result for sin and cos(hx+ky) and reuses it if appropriate.  In practice  *
 * hit/miss rates vary from 6:1 to 12:1 on a uniprocessor.  To take          *
 * advantage in  parallel the k-vector loop should be blocked  not sliced.   *
 *****************************************************************************/
/*   static int hits=0, misses=0;*/
static
void      qsincos(real *coshx, real *sinhx, real *cosky, 
		  real *sinky, real *coslz, real *sinlz, 
		  real *qcoskr, real *qsinkr, 
		  real *coshxky, real *sinhxky, 
		  int h, int k, int l, int hlast, int klast, int nsites)
{
   int is;
   real qckr, chxky;
   
   if( h != hlast || k != klast ) 
   {
      /*      misses++;*/
      if( k >= 0 )
	 for(is = 0; is < nsites; is++)
	 {	
	    chxky       = coshx[is]*cosky[is] - sinhx[is]*sinky[is];
	    sinhxky[is] = sinhx[is]*cosky[is] + coshx[is]*sinky[is];
	    coshxky[is] = chxky;
	 }
      else
	 for(is = 0; is < nsites; is++)
	 {
	    chxky       = coshx[is]*cosky[is] + sinhx[is]*sinky[is];
	    sinhxky[is] = sinhx[is]*cosky[is] - coshx[is]*sinky[is];
	    coshxky[is] = chxky;
	 }
   }
   /*hits++;*/
   if( l >= 0 )
   {
      for(is = 0; is < nsites; is++)
      {	
	 qckr =       (coshxky[is]*coslz[is] - sinhxky[is]*sinlz[is]);
	 qsinkr[is] = (sinhxky[is]*coslz[is] + coshxky[is]*sinlz[is]);
	 qcoskr[is] = qckr;
      }
   }
   else
   {
      for(is = 0; is < nsites; is++)
      {
	 qckr =       (coshxky[is]*coslz[is] + sinhxky[is]*sinlz[is]);
	 qsinkr[is] = (sinhxky[is]*coslz[is] - coshxky[is]*sinlz[is]);
	 qcoskr[is] = qckr;
      }
   }
}

/******************************************************************************
 * Calculate cos and sin of astar*x, bstar*y & cstar*z for each charged site  *
 * Use addition formulae to get sin(h*astar*x)=sin(Kx*x) etc for each site    *
 ******************************************************************************/
static
void trig_recur(real **chx, real **shx, 
		real *sin1x, real *cos1x, 
		real *site0, real *site1, real *site2, 
		real *kstar, 
		int hmax, int ns0, int ns1)
{
   int h, is;
   real *coshx, *sinhx, *cm1, *sm1, coss, kr;
   real ksx = kstar[0], ksy = kstar[1], ksz = kstar[2];

   sinhx = shx[0];    coshx = chx[0]; 
   VECTORIZE
      for(is = ns0; is < ns1; is++)
      {
	 coshx[is] = 1.0;
	 sinhx[is] = 0.0;
      }      

   if( hmax >= 1 )
   {
      coshx = chx[1]; sinhx = shx[1];
      VECTORIZE
	 for(is = ns0; is < ns1; is++)
	 {
	    kr = ksx*site0[is]+ksy*site1[is]+ksz*site2[is];
	    coshx[is] = cos(kr); sinhx[is] = sin(kr);
	 }

      memcp(cos1x+ns0, coshx+ns0, (ns1-ns0)*sizeof(real));
      memcp(sin1x+ns0, sinhx+ns0, (ns1-ns0)*sizeof(real));

      for(h = 2; h <= hmax; h++)
      {
         coshx = chx[h]; sinhx = shx[h];
         cm1 = chx[h-1]; sm1 = shx[h-1];
	 VECTORIZE
	    for(is = ns0; is < ns1; is++)
	    {
	       coss      = cm1[is]*cos1x[is] - sm1[is]*sin1x[is];
	       sinhx[is] = sm1[is]*cos1x[is] + cm1[is]*sin1x[is];
	       coshx[is] = coss;
	    }
      }
   }
}

/******************************************************************************
 * allocate_arrays(). Memory allocation for Ewald's arrays.  We malloc one    *
 * large block and divide it out so as to precisely control the relative      *
 * offsets of the arrays.  This is to avoid cache conflicts.                  *
 * Revert to bog-standard method for MSDOS				      *
 ******************************************************************************/
static real*  allocate_arrays(int nsarray, int hmax, int kmax, int lmax, 
			      real ***chx, real ***cky, real ***clz, 
			      real ***shx, real ***sky, real ***slz, 
			      real **cos1x, real **cos1y, real **cos1z, 
			      real **sin1x, real **sin1y, real **sin1z, 
			      real **qcoskr, real **qsinkr, 
			      real **coshxky, real **sinhxky)
{
   int h, k,l; 
   real *csp, *base;

#ifndef ALLOC_SEPARATELY
   /*
    * Attempt to cache-align these arrays for optimum performance.
    * This requires NON-ANSI pointer conversion and arithmetic.
    * It should work with any UNIX address space, so make it conditional.
    */
   base = dalloc((2*(hmax+kmax+lmax)+16)*nsarray+11*NLINE+NCACHE);
#if defined(ALLOC_ALIGN) && !defined(__BOUNDS_CHECKING_ON)
   csp = (real*)0 + (((base - (real*)0) - 1 | NCACHE-1) + 1);
#else
   csp = base;
#endif

   *chx = (real**)arralloc(sizeof(real*),1,0,hmax);
   *cky = (real**)arralloc(sizeof(real*),1,0,kmax);
   *clz = (real**)arralloc(sizeof(real*),1,0,lmax);
   *shx = (real**)arralloc(sizeof(real*),1,0,hmax);
   *sky = (real**)arralloc(sizeof(real*),1,0,kmax);
   *slz = (real**)arralloc(sizeof(real*),1,0,lmax);
 
   *qcoskr = csp; csp += nsarray + NLINE;
   *qsinkr = csp; csp += nsarray + NLINE;
   *coshxky = csp; csp += nsarray + NLINE;
   *sinhxky = csp; csp += nsarray + NLINE;
   for(h = 0; h <= hmax; h++, csp += nsarray)
      (*chx)[h] = csp;
   csp += NLINE;
   for(h = 0; h <= hmax; h++, csp += nsarray)
      (*shx)[h] = csp;
   csp += NLINE;
   for(k = 0; k <= kmax; k++, csp += nsarray)
      (*cky)[k] = csp;
   csp += NLINE;
   for(k = 0; k <= kmax; k++, csp += nsarray)
      (*sky)[k] = csp;
   csp += NLINE;
   for(l = 0; l <= lmax; l++, csp += nsarray)
      (*clz)[l] = csp;
   csp += NLINE;
   for(l = 0; l <= lmax; l++, csp += nsarray)
      (*slz)[l] = csp;
   csp += NLINE;
   *cos1x = csp;  csp += nsarray;
   *cos1y = csp;  csp += nsarray;
   *cos1z = csp;  csp += nsarray + NLINE;
   *sin1x = csp;  csp += nsarray;
   *sin1y = csp;  csp += nsarray;
   *sin1z = csp;  /*csp += nsarray;*/
#else
   *chx = (real**)arralloc((size_mt)sizeof(real),2, 0, hmax, 0, nsarray-1);
   *cky = (real**)arralloc((size_mt)sizeof(real),2, 0, kmax, 0, nsarray-1);
   *clz = (real**)arralloc((size_mt)sizeof(real),2, 0, lmax, 0, nsarray-1);
   *shx = (real**)arralloc((size_mt)sizeof(real),2, 0, hmax, 0, nsarray-1);
   *sky = (real**)arralloc((size_mt)sizeof(real),2, 0, kmax, 0, nsarray-1);
   *slz = (real**)arralloc((size_mt)sizeof(real),2, 0, lmax, 0, nsarray-1);
   *cos1x = (*chx)[1]; *cos1y = (*cky)[1]; *cos1z = (*clz)[1];
   *sin1x = (*shx)[1]; *sin1y = (*sky)[1]; *sin1z = (*slz)[1];

   csp = base = dalloc(2*nsarray);
   *qcoskr = csp; csp += nsarray;
   *qsinkr = csp; /*csp += nsarray;*/
#endif

   return base;
}
/******************************************************************************
 *  Ewald  Calculate reciprocal-space part of coulombic forces		      *
 ******************************************************************************/
void	ewald(real **site,              /* Site co-ordinate arrays       (in) */
	      real **site_force,        /* Site force arrays            (out) */ 
	      system_mp system,         /* System record                 (in) */
	      spec_mt *species,         /* Array of species records      (in) */
	      real *chg, 	        /* Array of site charges         (in) */
	      double *pe, 	        /* Potential energy             (out) */
	      real (*stress)[3])        /* Stress virial                (out) */
{
   mat_mt	hinvp;			/* Matrix of reciprocal lattice vects*/
   int		h, k, l;		/* Recip. lattice vector indices     */
   int		i, is, ssite;		/* Counters.			     */
   spec_mp	spec;			/* species[ispec]		     */
   int		nsites = system->nsites;
   double	pe_k,			/* Pot'l energy for current K vector */
                coeff, coeff2;		/* 2/(e0V) * A(K) & similar	     */
   double	r_4_alpha = -1.0/(4.0 * control.alpha * control.alpha);
   double	sqcoskr,sqsinkr,	/* Sum q(i) sin/cos(K.r(i))          */
                sqcoskrn, sqsinkrn,
                sqcoskrf, sqsinkrf;
   double	ksq,			/* Squared magnitude of K vector     */
                kcsq = SQR(control.k_cutoff);
   double	kx,ky,kz,kzt;
   vec_mt	kv;			/* (Kx,Ky,Kz)  			     */
   real		force_comp, kv0, kv1, kv2, sfx, sfy;
   struct	s_hkl *hkl, *phkl;
   int		nhkl = 0, nhklp;
   /*
    * Maximum values of h, k, l  s.t. |k| < k_cutoff
    */
   int		hmax = floor(control.k_cutoff/(2*PI)*moda(system->h)),
                kmax = floor(control.k_cutoff/(2*PI)*modb(system->h)),
                lmax = floor(control.k_cutoff/(2*PI)*modc(system->h));
   /*
    * Saved previous h and k to signal validity of coshxky, sinhxky caches.
    */
   int		hlast=-1000000,klast=-1000000;
   /*
    * lower and upper limits for parallel loops.   
    */
#if defined(SPMD) && defined(MPPMANY)
   int  	nsnode = (nsites+nthreads-1)/nthreads;
   int  	nsarr0 = nsnode * nthreads;
   int  	ns0 = MIN(nsites, nsnode * ithread), 
        	ns1 = MIN(nsites, nsnode * (ithread+1));
#else	    
   int		ns0 = 0, ns1 = nsites;
   int		nsarr0 = nsites;
#endif
   /*
    * Round up size of arrays to cache sub-multiple size.
    */
#if 1
   int		nsarray = (nsarr0 - 1 | NCACHE - 1) + 1;
#else
   int		nsarray = nsarr0;
#endif
   /*
    * Arrays for cos & sin (h x(i)), (k y(i)) and (l z(i)) eg chx[h][isite]
    * and pointers to a particular h,k or l eg coshx[is] = chh[2][is]
    */
   real		**chx, **cky, **clz, **shx, **sky, **slz;
   real		*cshkl;
   real		*coshx, *cosky, *coslz, *sinhx, *sinky, *sinlz;
   real		*cos1x, *cos1y, *cos1z, *sin1x, *sin1y, *sin1z;
   real		q;
   real		*site_fx = site_force[0],
                *site_fy = site_force[1],
                *site_fz = site_force[2];
   real		*qcoskr,		/* q(i) cos(K.R(i))	      */
                *qsinkr;		/* q(i) sin(K.R(i))	      */
   real		*sinhxky, *coshxky;
   double	vol = det(system->h);	/* Volume of MD cell		      */
   static	double	self_energy,	/* Constant self energy term	      */
                sheet_energy;	        /* Correction for non-neutral system. */
   static	boolean init = true;	/* Flag for the first call of function*/
   static       int     nsitesxf;       /* Number of non-framework sites.     */
   /*
    * Trig array set-up.  The precise storage layout is to avoid cache conflicts.
    */
   cshkl = allocate_arrays(nsarray, hmax, kmax, lmax,
			   &chx, &cky, &clz, &shx, &sky, &slz,
			   &cos1x, &cos1y, &cos1z, &sin1x, &sin1y, &sin1z, 
			   &qcoskr, &qsinkr,&sinhxky,&coshxky);
   /*
    * First call only - evaluate self energy term and store for subsequent calls
    * Self energy includes terms for non-framework sites only.
    */
   if(init)
   {
      double	sqsq = 0, sq = 0, sqxf, intra, r;
      int	js, frame_flag;

      self_energy = sheet_energy = 0;
      ssite = 0;
      spec = species; 
      while( spec < species+system->nspecies && ! spec->framework)
      {
	 intra = 0.0;
	 for(is = 0; is < spec->nsites; is++)
	    for(js = is+1; js < spec->nsites; js++)
	    {
	       r = DISTANCE(spec->p_f_sites[is], spec->p_f_sites[js]);
	       intra += chg[ssite+is] * chg[ssite+js]
		  *err_fn(control.alpha * r) / r;
	    }
	 self_energy += spec->nmols * intra;
	 ssite += spec->nsites*spec->nmols;
	 spec++;
      }
      nsitesxf = ssite;
      frame_flag = (spec != species+system->nspecies);

      for(is = 0; is < nsitesxf; is++)
      {
	 sq += chg[is];
	 sqsq += SQR(chg[is]);
      }
      self_energy += control.alpha / sqrt(PI) * sqsq;
      /*
       * Sqxf is total non-framework charge.  Calculate grand total in sq.
       */
      sqxf = sq;
      for(; is < nsites; is++)
	 sq += chg[is];
      /*
       *  Charged-system/uniform sheet correction terms (really in direct
       *  space term but included here for convenience).
       *  1) For charged framework only.
       */
      if( frame_flag )
      {
	 sheet_energy = PI*SQR(sq-sqxf) / (2.0*SQR(control.alpha));
	 message(NULLI, NULLP, INFO, FRACHG, 
		 (sq-sqxf)*CONV_Q, sheet_energy/vol*CONV_E);
      }
      /* 
       *  2) Case of entire system non-neutral.
       */
      if( fabs(sq)*CONV_Q > 1.0e-5)
      {
	 sheet_energy -= intra = PI*SQR(sq) / (2.0*SQR(control.alpha));
	 message(NULLI, NULLP, WARNING, SYSCHG, sq*CONV_Q, intra/vol*CONV_E);
      }

      note("Ewald self-energy = %f kJ/mol",self_energy*CONV_E);
   }

   if( ithread == 0 )
   {
      *pe -= self_energy;		/* Subtract self energy term	      */
      *pe += sheet_energy/vol;		/* Uniform charge correction	      */
      for(i=0; i<3; i++)
	 stress[i][i] += sheet_energy/vol;
   }
   
   invert(system->h, hinvp);		/* Inverse of h is matrix of r.l.v.'s */
   mat_sca_mul(2*PI, hinvp, hinvp);
   /*
    * Build array hkl[] of k vectors within cutoff. 
    * N. B. This code presupposes upper-triangular H matrix.
    */
   hkl = aalloc(4*(hmax+1)*(kmax+1)*(lmax+1), struct s_hkl);
   for(h = 0; h <= hmax; h++)
      for(k = (h==0 ? 0 : -kmax); k <= kmax; k++)
      {
	 kx = h*astar[0] + k*bstar[0];
	 ky = h*astar[1] + k*bstar[1];
	 kzt = h*astar[2] + k*bstar[2];
	 ksq = SQR(kx) + SQR(ky);
	 for(l = (h==0 && k==0 ? 1 : -lmax); l <= lmax; l++)
	 {
	    kz = kzt + l*cstar[2];
	    if( SQR(kz)+ksq < kcsq )
	    {
	       hkl[nhkl].h = h; hkl[nhkl].k = k; hkl[nhkl].l = l;
	       hkl[nhkl].kx = kx;
	       hkl[nhkl].ky = ky;
	       hkl[nhkl].kz = kz;
	       nhkl++;
	    }
	 }
      }
   if( init )
      note("%d K-vectors included in reciprocal-space sum",nhkl);
   init = false;
   /*
    * Calculate cos and sin of astar*x, bstar*y & cstar*z for each charged site
    * Use addition formulae to get sin(h*astar*x)=sin(Kx*x) etc for each site
    */
   trig_recur(chx,shx,sin1x,cos1x,site[0],site[1],site[2],astar,hmax,ns0,ns1);
   trig_recur(cky,sky,sin1y,cos1y,site[0],site[1],site[2],bstar,kmax,ns0,ns1);
   trig_recur(clz,slz,sin1z,cos1z,site[0],site[1],site[2],cstar,lmax,ns0,ns1);
   /*
    * Pre-multiply cos(lz) and sin(lz) by q to avoid doing it in the inner loop
    */
   for(l=0; l <= lmax; l++)
      for(is=ns0; is < ns1; is++)
      {
	 q = chg[is];
	 clz[l][is] *= q;
	 slz[l][is] *= q;
      }
#if defined(SPMD) && defined(MPPMANY)
   par_collect_all(chx[0]+ns0, chx[0], nsnode, nsarray, hmax+1);
   par_collect_all(shx[0]+ns0, shx[0], nsnode, nsarray, hmax+1);
   par_collect_all(cky[0]+ns0, cky[0], nsnode, nsarray, kmax+1);
   par_collect_all(sky[0]+ns0, sky[0], nsnode, nsarray, kmax+1);
   par_collect_all(clz[0]+ns0, clz[0], nsnode, nsarray, lmax+1);
   par_collect_all(slz[0]+ns0, slz[0], nsnode, nsarray, lmax+1);
#endif
   /*
    * Start of main loops over K vector indices h, k, l between -*max, *max etc.
    * To avoid calculating K and -K, only half of the K-space box is covered. 
    * Points on the axes are included once and only once. (0,0,0) is omitted.
    */
   nhklp = (nhkl+nthreads-1)/nthreads;
   for(phkl = hkl+ithread*nhklp; phkl < hkl+MIN((ithread+1)*nhklp,nhkl); phkl ++)
   {
      /*
       * Calculate actual K vector and its squared magnitude.
       */
      h  = phkl->h;	    k     = phkl->k;  l     = phkl->l;
      kv0 = kv[0] = phkl->kx; 
      kv1 = kv[1] = phkl->ky; 
      kv2 = kv[2] = phkl->kz;

      ksq = SUMSQ(kv);

      /*
       * Calculate pre-factors A(K) etc
       */
      coeff  = 2.0 / (EPS0 * vol) * exp(ksq * r_4_alpha) / ksq;
      coeff2 = 2.0 * (1.0 - ksq * r_4_alpha) / ksq;
      
      /*
       * Set pointers to array of cos (h*astar*x) (for efficiency & vectorisation)
       */
      coshx = chx[h];  cosky = cky[abs(k)]; coslz = clz[abs(l)];
      sinhx = shx[h];  sinky = sky[abs(k)]; sinlz = slz[abs(l)];

      /*
       * Calculate q(i)*cos/sin(K.R(i)) by addition formulae. Note handling of
       * negative k and l by using sin(-x) = -sin(x), cos(-x) = cos(x). For
       * efficiency & vectorisation there is a loop for each case.
       */
      qsincos(coshx,sinhx,cosky,sinky,coslz,sinlz,
	      qcoskr,qsinkr,sinhxky, coshxky, h, k, l, hlast, klast, nsites);
      hlast = h; klast = k;
      sqcoskrn = sum(nsitesxf, qcoskr, 1);
      sqsinkrn = sum(nsitesxf, qsinkr, 1);
      sqcoskrf = sum(nsites-nsitesxf, qcoskr+nsitesxf, 1);
      sqsinkrf = sum(nsites-nsitesxf, qsinkr+nsitesxf, 1);
      sqcoskr = sqcoskrn + sqcoskrf;
      sqsinkr = sqsinkrn + sqsinkrf;
      
      /*
       * Evaluate potential energy contribution for this K and add to total.
       * Exclude frame-frame interaction terms.
       */
      pe_k = 0.5 * coeff * (sqcoskrn*(sqcoskrn+sqcoskrf+sqcoskrf) +
			    sqsinkrn*(sqsinkrn+sqsinkrf+sqsinkrf));
      *pe += pe_k;

      sqsinkr *= coeff; sqcoskr *= coeff;
      sqsinkrn *= coeff; sqcoskrn *= coeff;
      /*
       * Calculate long-range coulombic contribution to stress tensor
       */
      stress[0][0] += pe_k - pe_k * coeff2 * kv[0] * kv[0];
      stress[0][1] -=        pe_k * coeff2 * kv[0] * kv[1];
      stress[0][2] -=        pe_k * coeff2 * kv[0] * kv[2];
      stress[1][1] += pe_k - pe_k * coeff2 * kv[1] * kv[1];
      stress[1][2] -=        pe_k * coeff2 * kv[1] * kv[2];
      stress[2][2] += pe_k - pe_k * coeff2 * kv[2] * kv[2];
      /*
       * Evaluation of site forces.   Non-framework sites interact with all others
       */
VECTORIZE
      for(is = 0; is < nsitesxf; is++)
      {
	 force_comp = qsinkr[is]*sqcoskr - qcoskr[is]*sqsinkr;
	 sfx =  site_fx[is] + kv0 * force_comp;
	 sfy =  site_fy[is] + kv1 * force_comp;
	 site_fz[is] +=       kv2 * force_comp;
	 site_fx[is] = sfx;
	 site_fy[is] = sfy;
      }
      /*
       *  Framework sites -- only interact with non-framework sites
       */   
VECTORIZE
      for(is = nsitesxf; is < nsites; is++)
      {
	 force_comp = qsinkr[is]*sqcoskrn - qcoskr[is]*sqsinkrn;
	 sfx =  site_fx[is] + kv0 * force_comp;
	 sfy =  site_fy[is] + kv1 * force_comp;
	 site_fz[is] +=       kv2 * force_comp;
	 site_fx[is] = sfx;
	 site_fy[is] = sfy;
      }
   /*
    * End of loop over K vectors.
    */
   }
   afree((gptr*)chx); afree((gptr*)cky); afree((gptr*)clz); 
   afree((gptr*)shx); afree((gptr*)sky); afree((gptr*)slz);
   xfree(cshkl);
   xfree(hkl);
   /*   message(NULLI, NULLP, INFO, "Trig Caching: %d hits and %d misses", hits, misses);*/
}
/******************************************************************************
 *  site_id_sort(). Construct permutation arrays for site list sorted by id.  *
 ******************************************************************************/
site_id_sort(system_mt *system, spec_mt *species, int* site_idx, int *site_permute)
{
   int id, imol, ims, is, js;
   spec_mt *spec;

   is = 0;
   site_idx[0] = 0;
   for ( id=1; id < system->max_id; id++)
   {
      js = 0;
      for(spec = species; spec < &species[system->nspecies]; spec++)
	 for(imol=0; imol<spec->nmols; imol++)
	    for(ims=0; ims < spec->nsites; ims++)
	    {
	       if(spec->site_id[ims] == id)
		  site_permute[is++] = js;
	       js++;
	    }
      site_idx[id] = is;
   }
#ifdef DEBUG
   {
      int suma=0;
      for(is = 0; is < system->nsites; is++)
	 suma += site_permute[is];

      fprintf(stderr,"Nsites = %d, n(n-1)/2 = %d, sum(site_permute) = %d\n",system->nsites,system->nsites*(system->nsites-1)/2,suma);
   }
#endif
}
/******************************************************************************
 *  ewald46  Calculate reciprocal-space part of 1/r^4 and 1/r^6 forces	      *
 ******************************************************************************/
void	ewald46(real **site,            /* Site co-ordinate arrays       (in) */
		real **site_force,      /* Site force arrays            (out) */ 
		system_mp system,       /* System record                 (in) */
		spec_mt *species,       /* Array of species records      (in) */
		real  **pot4,           /* Potential parameters for r^-4 (in) */
		real  **pot6,           /* Potential parameters for r^-6 (in) */
		double *pe, 	        /* Potential energy             (out) */
		real (*stress)[3])      /* Stress virial                (out) */
{
   mat_mt	hinvp;			/* Matrix of reciprocal lattice vects*/
   int		h, k, l;		/* Recip. lattice vector indices     */
   int		i, is, js, idi, idj, nid;	/* Counters.			     */
   spec_mp	spec;			/* species[ispec]		     */
   int		nsites = system->nsites;
   double	struct_factor_4_sq, struct_factor_6_sq ;
   double	pe_k,			/* Pot'l energy for current K vector */
                coeff, coeff2;		/* 2/(e0V) * A(K) & similar	     */
   double	alpha = control.alpha46;
   double	*scoskr = arralloc(sizeof(double), 1, 0, system->max_id-1),
                *ssinkr = arralloc(sizeof(double), 1, 0, system->max_id-1); /* Sum of sin/cos(K.r(i))          */
   double	ksq,			/* Squared magnitude of K vector     */
                kmod,
      		kover2a,
      		erfc_k2a,
      		exp_k24a2,
      		ak4, ak6,
      		dak4dk, dak6dk,
                kcsq = SQR(control.k_cutoff46);
   double	kx,ky,kz,kzt;
   vec_mt	kv;			/* (Kx,Ky,Kz)  			     */
   real		force_comp, kv0, kv1, kv2, sfx, sfy;
   struct	s_hkl *hkl, *phkl;
   int		nhkl = 0, nhklp;
   /*
    * Maximum values of h, k, l  s.t. |k| < k_cutoff
    */
   int		hmax = floor(control.k_cutoff46/(2*PI)*moda(system->h)),
                kmax = floor(control.k_cutoff46/(2*PI)*modb(system->h)),
                lmax = floor(control.k_cutoff46/(2*PI)*modc(system->h));
   /*
    * Saved previous h and k to signal validity of coshxky, sinhxky caches.
    */
   int		hlast=-1000000,klast=-1000000;
   /*
    * lower and upper limits for parallel loops.   
    */
#if defined(SPMD) && defined(MPPMANY)
   int  	nsnode = (nsites+nthreads-1)/nthreads;
   int  	nsarr0 = nsnode * nthreads;
   int  	ns0 = MIN(nsites, nsnode * ithread), 
        	ns1 = MIN(nsites, nsnode * (ithread+1));
#else	    
   int		ns0 = 0, ns1 = nsites;
   int		nsarr0 = nsites;
#endif
   /*
    * Round up size of arrays to cache sub-multiple size.
    */
#if 1
   int		nsarray = (nsarr0 - 1 | NCACHE - 1) + 1;
#else
   int		nsarray = nsarr0;
#endif
   /*
    * Arrays for cos & sin (h x(i)), (k y(i)) and (l z(i)) eg chx[h][isite]
    * and pointers to a particular h,k or l eg coshx[is] = chh[2][is]
    */
   real		**chx, **cky, **clz, **shx, **sky, **slz;
   real		*cshkl;
   real		*coshx, *cosky, *coslz, *sinhx, *sinky, *sinlz;
   real		*cos1x, *cos1y, *cos1z, *sin1x, *sin1y, *sin1z;
   real		q;
   real		*site_x = dalloc(nsarray),
                *site_y = dalloc(nsarray),
                *site_z = dalloc(nsarray),
   		*site_fx = dalloc(nsarray),
                *site_fy = dalloc(nsarray),
                *site_fz = dalloc(nsarray);
   static int	*site_idx;
   static int   *site_permute;
   boolean	dopot4, dopot6;
   int	        rel;
   real		*coskr,		/* q(i) cos(K.R(i))	      */
                *sinkr;		/* q(i) sin(K.R(i))	      */
   real		*sinhxky, *coshxky;
   double	vol = det(system->h);	/* Volume of MD cell		      */
   static	double	self_energy,	/* Constant self energy term	      */
                sheet_energy;	        /* Correction for non-neutral system. */
   static	boolean init = true;	/* Flag for the first call of function*/

   zero_real(site_fx, nsites);   zero_real(site_fy, nsites);   zero_real(site_fz, nsites);
   /*
    * Trig array set-up.  The precise storage layout is to avoid cache conflicts.
    */
   cshkl = allocate_arrays(nsarray, hmax, kmax, lmax,
			   &chx, &cky, &clz, &shx, &sky, &slz,
			   &cos1x, &cos1y, &cos1z, &sin1x, &sin1y, &sin1z, 
			   &coskr, &sinkr,&sinhxky,&coshxky);
   /*
    * First call only - evaluate self energy term and store for subsequent calls
    */
   if(init)
   {
      double	sq4 = 0, sq6 = 0, sqq4 = 0, sqq6 = 0, intra, r, alphar, expar;
      int	js, idi, idj, ni, nj;

      site_idx = ialloc(system->max_id);
      site_permute = ialloc(nsites);

      site_id_sort(system, species, site_idx, site_permute);

      self_energy = sheet_energy = 0;
      spec = species; 
      while( spec < species+system->nspecies && ! spec->framework)
      {
	 intra = 0.0;
	 for(is = 0; is < spec->nsites; is++)
	 {
	    idi = spec->site_id[is];
	    for(js = is+1; js < spec->nsites; js++)
	    {
	       idj = spec->site_id[js];
	       r = DISTANCE(spec->p_f_sites[is], spec->p_f_sites[js]);
	       alphar = alpha*r;
	       expar=exp(-SQR(alphar));
	       intra += pot4[idi][idj]
		  *(1.0-(1.0+SQR(alphar))*expar)/FOURTH(r);
	       intra += pot6[idi][idj]
		  *(1.0-(1.0+SQR(alphar)+0.5*FOURTH(alphar))*expar)/SIXTH(r);
	    }
	 }
	 self_energy += spec->nmols * intra;

	 spec++;
      }

      for(idi=1; idi < system->max_id; idi++)
      {
	 ni = site_idx[idi] - site_idx[idi-1];
         sq4 += ni*pot4[idi][idi];
         sq6 += ni*pot6[idi][idi];
         for(idj=1; idj < system->max_id; idj++)
         {
	    nj = site_idx[idj] - site_idx[idj-1];
	    sqq4 += ni*nj*pot4[idi][idj];
            sqq6 += ni*nj*pot6[idi][idj];
         }
      }

      sheet_energy += PI32*(alpha*sqq4 + CUBE(alpha)/6.0*sqq6);
      self_energy += 0.25*FOURTH(alpha)*sq4 + SIXTH(alpha)/12.0*sq6;
      note("Ewald-4-6 G=0 term = %f kJ/mol",sheet_energy/vol*CONV_E);
      note("Ewald-4-6 self-energy = %f kJ/mol",self_energy*CONV_E);
   }

   if( ithread == 0 )
   {
      *pe -= self_energy;		/* Subtract self energy term	      */
      *pe += sheet_energy/vol;		/* Uniform charge correction	      */
      for(i=0; i<3; i++)
	 stress[i][i] += sheet_energy/vol;
   }
   
   invert(system->h, hinvp);		/* Inverse of h is matrix of r.l.v.'s */
   mat_sca_mul(2*PI, hinvp, hinvp);
   /*
    * Build array hkl[] of k vectors within cutoff. 
    * N. B. This code presupposes upper-triangular H matrix.
    */
   hkl = aalloc(4*(hmax+1)*(kmax+1)*(lmax+1), struct s_hkl);
   for(h = 0; h <= hmax; h++)
      for(k = (h==0 ? 0 : -kmax); k <= kmax; k++)
      {
	 kx = h*astar[0] + k*bstar[0];
	 ky = h*astar[1] + k*bstar[1];
	 kzt = h*astar[2] + k*bstar[2];
	 ksq = SQR(kx) + SQR(ky);
	 for(l = (h==0 && k==0 ? 1 : -lmax); l <= lmax; l++)
	 {
	    kz = kzt + l*cstar[2];
	    if( SQR(kz)+ksq < kcsq )
	    {
	       hkl[nhkl].h = h; hkl[nhkl].k = k; hkl[nhkl].l = l;
	       hkl[nhkl].kx = kx;
	       hkl[nhkl].ky = ky;
	       hkl[nhkl].kz = kz;
	       nhkl++;
	    }
	 }
      }
   if( init )
      note("%d K-vectors included in reciprocal-space sum",nhkl);
   init = false;

   for(is=0; is < nsites; is++)
   {
      rel = site_permute[is];
      site_x[is] = site[0][rel];
      site_y[is] = site[1][rel];
      site_z[is] = site[2][rel];
   }
   /*
    * Calculate cos and sin of astar*x, bstar*y & cstar*z for each charged site
    * Use addition formulae to get sin(h*astar*x)=sin(Kx*x) etc for each site
    */
   trig_recur(chx,shx,sin1x,cos1x,site_x,site_y,site_z,astar,hmax,ns0,ns1);
   trig_recur(cky,sky,sin1y,cos1y,site_x,site_y,site_z,bstar,kmax,ns0,ns1);
   trig_recur(clz,slz,sin1z,cos1z,site_x,site_y,site_z,cstar,lmax,ns0,ns1);
#if defined(SPMD) && defined(MPPMANY)
   par_collect_all(chx[0]+ns0, chx[0], nsnode, nsarray, hmax+1);
   par_collect_all(shx[0]+ns0, shx[0], nsnode, nsarray, hmax+1);
   par_collect_all(cky[0]+ns0, cky[0], nsnode, nsarray, kmax+1);
   par_collect_all(sky[0]+ns0, sky[0], nsnode, nsarray, kmax+1);
   par_collect_all(clz[0]+ns0, clz[0], nsnode, nsarray, lmax+1);
   par_collect_all(slz[0]+ns0, slz[0], nsnode, nsarray, lmax+1);
#endif
   /*
    * Start of main loops over K vector indices h, k, l between -*max, *max etc.
    * To avoid calculating K and -K, only half of the K-space box is covered. 
    * Points on the axes are included once and only once. (0,0,0) is omitted.
    */
   nhklp = (nhkl+nthreads-1)/nthreads;
   for(phkl = hkl+ithread*nhklp; phkl < hkl+MIN((ithread+1)*nhklp,nhkl); phkl ++)
   {
      /*
       * Calculate actual K vector and its squared magnitude.
       */
      h  = phkl->h;	    k     = phkl->k;  l     = phkl->l;
      kv0 = kv[0] = phkl->kx; 
      kv1 = kv[1] = phkl->ky; 
      kv2 = kv[2] = phkl->kz;

      ksq = SUMSQ(kv);
      kmod = sqrt(ksq);
      kover2a = kmod/(2*alpha);
      erfc_k2a = erfc(kover2a);
      exp_k24a2 = exp( - ksq/(4.0*SQR(alpha)));

      /*
       * Calculate pre-factors A(K) etc
       */
      ak4 = (-0.5*PISQ*kmod*erfc_k2a + alpha*PI32*exp_k24a2) * 2.0/vol;
      dak4dk = (-0.5*PISQ*erfc_k2a) * 2.0/vol;
      ak6 = (PISQ/24.0*kmod*ksq*erfc_k2a + PI32/6.0*alpha*(SQR(alpha)-0.5*ksq)*exp_k24a2) * 2.0/vol;
      dak6dk = ((1.0/24.0)*PI32*(3.0*ROOTPI*erfc_k2a*ksq - 6.0*alpha*kmod*exp_k24a2)) * 2.0/vol;
      
      /*
       * Set pointers to array of cos (h*astar*x) (for efficiency & vectorisation)
       */
      coshx = chx[h];  cosky = cky[abs(k)]; coslz = clz[abs(l)];
      sinhx = shx[h];  sinky = sky[abs(k)]; sinlz = slz[abs(l)];

      /*
       * Calculate q(i)*cos/sin(K.R(i)) by addition formulae. Note handling of
       * negative k and l by using sin(-x) = -sin(x), cos(-x) = cos(x). For
       * efficiency & vectorisation there is a loop for each case.
       */
      qsincos(coshx,sinhx,cosky,sinky,coslz,sinlz,
	      coskr,sinkr,sinhxky, coshxky, h, k, l, hlast, klast, nsites);
      hlast = h; klast = k;

      for(idi=1; idi < system->max_id; idi++)
      {
	 nid = site_idx[idi] - site_idx[idi-1];
	 js = site_idx[idi-1];
	 scoskr[idi] = sum(nid, &coskr[js], 1);
	 ssinkr[idi] = sum(nid, &sinkr[js], 1);
      }
     
      /*
       * Evaluate potential energy contribution for this K and add to total.
       */
      struct_factor_4_sq = 0.0; struct_factor_6_sq = 0.0; 
      for(idi=1; idi < system->max_id; idi++)
      {
	 for(idj=1; idj < system->max_id; idj++)
	 {
	    dopot4 = fabs(pot4[idi][idj]) > POTZERO_TOL;
	    dopot6 = fabs(pot6[idi][idj]) > POTZERO_TOL;
	    coeff = 0.0;
	    if(dopot4) 
	    {
	       struct_factor_4_sq += pot4[idi][idj]*(scoskr[idi]*scoskr[idj]+ssinkr[idi]*ssinkr[idj]);
	       coeff += pot4[idi][idj]*ak4;
	    }

	    if(dopot6) 
	    {
	       struct_factor_6_sq += pot6[idi][idj]*(scoskr[idi]*scoskr[idj]+ssinkr[idi]*ssinkr[idj]);
	       coeff += pot6[idi][idj]*ak6;
	    }
	    if( dopot4 || dopot6 )
	    {
	       /*
		* Evaluation of site forces. 
		*/
	       for(js = site_idx[idj-1]; js < site_idx[idj]; js++)
	       {
		  force_comp = coeff*(sinkr[js]*scoskr[idi]-coskr[js]*ssinkr[idi]);
		  sfx =  site_fx[js] + kv0 * force_comp;
		  sfy =  site_fy[js] + kv1 * force_comp;
		  site_fz[js] +=       kv2 * force_comp;
		  site_fx[js] = sfx;
		  site_fy[js] = sfy;
	       }
	    }
	 }
      }

      pe_k = ak4*struct_factor_4_sq+ak6*struct_factor_6_sq;
      *pe += pe_k;

      /*
       * Calculate contribution to stress tensor
       */
      coeff = (dak4dk*struct_factor_4_sq+dak6dk*struct_factor_6_sq)/kmod;
      stress[0][0] += pe_k + coeff * kv[0] * kv[0];
      stress[0][1] +=        coeff * kv[0] * kv[1];
      stress[0][2] +=        coeff * kv[0] * kv[2];
      stress[1][1] += pe_k + coeff * kv[1] * kv[1];
      stress[1][2] +=        coeff * kv[1] * kv[2];
      stress[2][2] += pe_k + coeff * kv[2] * kv[2];
   /*
    * End of loop over K vectors.
    */
   }
   /*
    * Scatter/add forces to main force arrays
    */
   for(is=0; is < nsites; is++)
   {
      rel = site_permute[is];
      site_force[0][rel] += site_fx[is];
      site_force[1][rel] += site_fy[is];
      site_force[2][rel] += site_fz[is];
   }

   afree((gptr*)chx); afree((gptr*)cky); afree((gptr*)clz); 
   afree((gptr*)shx); afree((gptr*)sky); afree((gptr*)slz);
   afree((gptr*)scoskr);    afree((gptr*)ssinkr);
   xfree(site_x); xfree(site_y); xfree(site_z);
   xfree(site_fx); xfree(site_fy); xfree(site_fz);
   xfree(cshkl);
   xfree(hkl);
   /*   message(NULLI, NULLP, INFO, "Trig Caching: %d hits and %d misses", hits, misses);*/
}
