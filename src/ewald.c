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
 *       Revision 2.8.1.5  1996/03/25 17:36:15  keith
 *       Fixed memory leak introduced in previous revision.
 *       Now prints # k-vectors just like RKL version
 *
 *       Revision 2.8.1.4  1996/02/07 20:51:44  keith
 *       Restructured for convergence of std and RIL versions
 *        - Implemented pre-computed k-vector set
 *        - Re-based trig arrays to start at zero for ANSI conformance.
 *
 *       Revision 2.8.1.3  1996/01/25 21:01:41  keith
 *       Fixed bug in allocation of sites to processors which caused
 *       crash on large # procs.
 *
 *       Revision 2.8.1.2  1995/12/06 15:07:50  keith
 *       Fixed bug which caused core dump for small k-cutoff (hmax=0)
 *
 *       Revision 2.8.1.2  1994/12/30 11:53:51  keith
 *       Finxed bug which caused core dump for small k-cutoff (hmax=0)
 *
 * Revision 2.8.1.1  1994/07/19  10:40:57  keith
 * Implementation of W. Smith's RIL parallelization strategy
 * (Comp Phys Commun, 67, (1992) 293-406
 * This involves distributing memory better but does communication
 * in the inner k-vector loop.  On most machines this seems to
 * defeat the parallelization altogether - it runs more slowly
 * in parallel on the Titan than in serial mode even for 32772 sites.
 *
 * Revision 2.8  1994/06/22  09:37:02  keith
 * Performance optimization of "trig rules" loops.
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
 * (pot_dim, potspec, prog_unit) are declared with CONST
 * qualifier macro which evaluates to "const" or nil
 * depending on ANSI/K&R environment.
 * Also moved as many "write" instantiations of "control"
 * members as possible to "startup", "main" leaving just
 * "dump".
 *
 * Declared as "static"  all functions which should be.
 *
 * Revision 2.5  1994/01/18  13:32:27  keith
 * Null update for XDR portability release
 *
 * Revision 2.3  93/10/28  10:27:48  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.1  93/09/02  12:31:55  keith
 * Optimized qsincos() -- should give up to 25% speed improvement on
 * compilers without assert no aliasing options.
 * 
 * Revision 2.1  93/05/17  10:42:22  keith
 * Optimized qsincos() -- should give up to 25% speed improvement on
 * compilers without assert no aliasing options.
 * 
 * Revision 2.0  93/03/15  14:49:02  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.23  93/03/12  12:22:38  keith
 * Reorganized defines to recognise all ANSI (__type__) forms.
 * Moved spxpy() from aux.c to force.c and force_parallel.c
 * 
 * 
 * Revision 1.23  93/03/12  12:21:50  keith
 * *** empty log message ***
 * 
 * Revision 1.22  93/03/09  15:58:28  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.21  92/06/26  17:02:58  keith
 * Got rid of assumption that memory returned by talloc() or
 * arralloc() is zeroed.  This enhances ANSI compatibility.
 * Removed memory zeroing from alloc.c() in consequence.
 * 
 * Revision 1.20  91/11/26  10:26:34  keith
 * Corrected calculation of sheet energy term for charged framework.
 * Split force loop so as to omit frame-frame force (and stress) terms.
 * 
 * Revision 1.19  91/08/15  18:11:52  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_mt and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.18  91/05/29  16:33:01  keith
 * Modified code for speed improvement in TITAN
 * 
 * Revision 1.17  91/03/12  15:42:31  keith
 * Tidied up typedefs size_mt and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.16  91/02/07  16:52:18  keith
 * Rewrote trig identity loops for better vectorization on Titan.
 * Finally deleted ancient commented-out code (#if OLDEWALD and VCALLS).
 * 
 * Revision 1.15  90/09/28  13:29:15  keith
 * Inserted braces around VECTORIZE directives and changed include files
 * for STARDtardent 3000 series (via cond. comp symbol "ardent").
 * 
 * Revision 1.14  90/08/29  11:01:19  keith
 * Modified to keep consistency with ewald_parallel.c r1.8
 * 
 * Revision 1.13  90/08/02  15:50:17  keith
 * Modified to exclude framework-framework interactions.
 * N.B. Excluded from pe and stress but NOT forces (as they sum to 0).
 * 
 * Revision 1.12  90/05/16  18:40:04  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.11  90/05/02  15:33:27  keith
 * Make declaration of saxpy() conditional along with use.
 * 
 * Revision 1.10  90/01/15  12:24:05  keith
 * Corrected declaration of arralloc from void* to char* to keep lint happy.
 * 
 * Revision 1.9  90/01/01  20:07:20  keith
 * Parcelled up generation of qcoskr etc into separate function and
 * created temp's site_fx etc to point at site_force[0] etc.
 * - Generates substabtially better code on Stellar.
 * 
 * Revision 1.8  89/12/22  19:31:53  keith
 * New version of arralloc() orders memory so that pointers come FIRST.
 * This means you can simply free() the pointer returned (if l.b. = 0).
 * 
 * Revision 1.7  89/12/21  16:29:47  keith
 * Reversed indices in 'site' and 'site_force' to allow stride of 1 in ewald.
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
static char *RCSid = "$Header: /home/eeyore_data/keith/md/moldy/RCS/ewald.c,v 2.8.1.5 1996/03/25 17:36:15 keith Exp $";
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
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
void            afree();	       /* Free allocated array	      	      */
double	err_fn();			/* Error function		      */
double	det();				/* Determinant of 3x3 matrix	      */
void	invert();			/* Inverts a 3x3 matrix		      */
void	mat_vec_mul();			/* Multiplies a 3x3 matrix by 3xN vect*/
void	mat_sca_mul();			/* Multiplies a 3x3 matrix by scalar  */
void	transpose();			/* Transposes a 3x3 matrix	      */
void    zero_real();            	/* Initialiser                        */
void    zero_double();          	/* Initialiser                        */
double	sum();				/* Sum of elements of 'real' vector   */
#if defined(ANSI) || defined(__STDC__)
gptr	*arralloc(size_mt,int,...); 	/* Array allocator		      */
void	note(char *,...);		/* Write a message to the output file */
void	message(int *,...);		/* Write a warning or error message   */
#else
gptr	*arralloc();	        	/* Array allocator		      */
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
#endif
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
 *****************************************************************************/
static
void      qsincos(coshx,sinhx,cosky,sinky,coslz,sinlz,chg,
		  qcoskr,qsinkr,k,l,nsites)
real coshx[], sinhx[], cosky[], sinky[], coslz[], sinlz[],
     chg[], qcoskr[], qsinkr[];
int  k,l,nsites;
{
   int is;
   real qckr;
   
   if( k >= 0 )
      if( l >= 0 )
      {
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qckr = chg[is]*(
		  (coshx[is]*cosky[is] - sinhx[is]*sinky[is])*coslz[is] 
                - (sinhx[is]*cosky[is] + coshx[is]*sinky[is])*sinlz[is]);
	    qsinkr[is] = chg[is]*(
                  (sinhx[is]*cosky[is] + coshx[is]*sinky[is])*coslz[is] 
		+ (coshx[is]*cosky[is] - sinhx[is]*sinky[is])*sinlz[is]);
	    qcoskr[is] = qckr;
	 }
      }
      else
      {
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qckr = chg[is]*(
		  (coshx[is]*cosky[is] - sinhx[is]*sinky[is])*coslz[is] 
                + (sinhx[is]*cosky[is] + coshx[is]*sinky[is])*sinlz[is]);
	    qsinkr[is] = chg[is]*(
                  (sinhx[is]*cosky[is] + coshx[is]*sinky[is])*coslz[is] 
		- (coshx[is]*cosky[is] - sinhx[is]*sinky[is])*sinlz[is]);
	    qcoskr[is] = qckr;
	 }
      }
   else
      if( l >= 0 )
      {
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qckr = chg[is]*(
		  (coshx[is]*cosky[is] + sinhx[is]*sinky[is])*coslz[is] 
                - (sinhx[is]*cosky[is] - coshx[is]*sinky[is])*sinlz[is]);
	    qsinkr[is] = chg[is]*(
                  (sinhx[is]*cosky[is] - coshx[is]*sinky[is])*coslz[is] 
		+ (coshx[is]*cosky[is] + sinhx[is]*sinky[is])*sinlz[is]);
	    qcoskr[is] = qckr;
	 }
      }
      else
      {
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qckr = chg[is]*(
		  (coshx[is]*cosky[is] + sinhx[is]*sinky[is])*coslz[is] 
                + (sinhx[is]*cosky[is] - coshx[is]*sinky[is])*sinlz[is]);
	    qsinkr[is] = chg[is]*(
                  (sinhx[is]*cosky[is] - coshx[is]*sinky[is])*coslz[is] 
		- (coshx[is]*cosky[is] + sinhx[is]*sinky[is])*sinlz[is]);
	    qcoskr[is] = qckr;
	 }
     }
}

/******************************************************************************
 * Calculate cos and sin of astar*x, bstar*y & cstar*z for each charged site  *
 * Use addition formulae to get sin(h*astar*x)=sin(Kx*x) etc for each site    *
 ******************************************************************************/
static
void trig_recur(chx,shx,sin1x,cos1x,site0,site1,site2,kstar,hmax,ns0,ns1)
real **chx, **shx, sin1x[],cos1x[], *site0, *site1, *site2, *kstar;
int ns0, ns1, hmax;
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
static real*  allocate_arrays(nsarray, hmax, kmax, lmax,
			      chx, cky, clz, shx, sky, slz, chg,
			      cos1x, cos1y, cos1z, sin1x, sin1y, sin1z, 
			      qcoskr, qsinkr)
int	nsarray, hmax, kmax, lmax;
real	***chx, ***cky, ***clz, ***shx, ***sky, ***slz;
real	**cos1x, **cos1y, **cos1z, **sin1x, **sin1y, **sin1z;
real	**chg, **qcoskr, **qsinkr;
{
   int h, k,l; 
   real *csp, *base;

#ifndef __MSDOS__
   /*
    * Attempt to cache-align these arrays for optimum performance.
    * This requires NON-ANSI pointer conversion and arithmetic.
    * It should work with any UNIX address space, so make it conditional.
    */
   base = dalloc((2*(hmax+kmax+lmax)+15)*nsarray+10*NLINE+NCACHE);
#if defined(unix) && !defined(__BOUNDS_CHECKING_ON)
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
   *chg = csp;    csp += nsarray + NLINE;
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
   *sin1z = csp;  csp += nsarray;
#else
   *chx = (real**)arralloc((size_mt)sizeof(real),2, 0, hmax, 0, nsarray-1);
   *cky = (real**)arralloc((size_mt)sizeof(real),2, 0, kmax, 0, nsarray-1);
   *clz = (real**)arralloc((size_mt)sizeof(real),2, 0, lmax, 0, nsarray-1);
   *shx = (real**)arralloc((size_mt)sizeof(real),2, 0, hmax, 0, nsarray-1);
   *sky = (real**)arralloc((size_mt)sizeof(real),2, 0, kmax, 0, nsarray-1);
   *slz = (real**)arralloc((size_mt)sizeof(real),2, 0, lmax, 0, nsarray-1);
   *cos1x = (*chx)[1]; *cos1y = (*cky)[1]; *cos1z = (*clz)[1];
   *sin1x = (*shx)[1]; *sin1y = (*sky)[1]; *sin1z = (*slz)[1];

   csp = base = dalloc(3*nsarray);
   *chg = csp;    csp += nsarray;
   *qcoskr = csp; csp += nsarray;
   *qsinkr = csp; csp += nsarray;
#endif

   return base;
}
/******************************************************************************
 *  Ewald  Calculate reciprocal-space part of coulombic forces		      *
 ******************************************************************************/
void	ewald(site,site_force,system,species,chgx,pe,stress)
real		**site,			/* Site co-ordinate arrays	 (in) */
		**site_force;		/* Site force arrays		(out) */
system_mp	system;			/* System record		 (in) */
spec_mt	species[];			/* Array of species records	 (in) */
real		chgx[];			/* Array of site charges	 (in) */
double		*pe;			/* Potential energy		(out) */
mat_mt		stress;			/* Stress virial		(out) */
{
   mat_mt	hinvp;			/* Matrix of reciprocal lattice vects*/
   int		h, k, l;		/* Recip. lattice vector indices     */
   int		i, j, is, ssite;	/* Counters.			     */
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
   int		nhkl = 0;
/*
 * Maximum values of h, k, l  s.t. |k| < k_cutoff
 */
   int		hmax = floor(control.k_cutoff/(2*PI)*moda(system->h)),
		kmax = floor(control.k_cutoff/(2*PI)*modb(system->h)),
		lmax = floor(control.k_cutoff/(2*PI)*modc(system->h));
   real		sqexpkr[4];
   mat_mt	stress_ew;
/*
 * lower and upper limits for parallel loops.  This doles out the sites
 * in parcels of "nsnode0" sites on "nns" threads and "nsnode0+1" sites
 * on the rest. "ns0" marks the beginning of the sites for this thread
 * and "nsarr0" the number of sites allocated to the node.  
 * N. B. The parcelling algorithm in W. Smith's paper FAILS if
 *       nthreads**2 > nsites.
 */
   int		nsnode0 = nsites/nthreads;
   int		nns = nthreads*(nsnode0 + 1 ) - nsites;
   int		nsarr0 = nsnode0 + ((ithread < nns)?0:1);
   int		ns0 = MIN(ithread, nns)*nsnode0 + MAX(ithread-nns,0)*(nsnode0+1);
   int		ns1 = ns0 + nsarr0, ns0f, ns1f;
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
   real		*chg;
   real		*site_fx = site_force[0]+ns0,
   		*site_fy = site_force[1]+ns0,
   		*site_fz = site_force[2]+ns0;
   real		*qcoskr,		/* q(i) cos(K.R(i))	      */
		*qsinkr;		/* q(i) sin(K.R(i))	      */
   double	vol = det(system->h);	/* Volume of MD cell		      */
   static	double	self_energy,	/* Constant self energy term	      */
   			sheet_energy;	/* Correction for non-neutral system. */
   static	boolean init = true;	/* Flag for the first call of function*/
   static	int	nsitesxf;	/* Number of non-framework sites.     */
/*
 * Trig array set-up.  The precise storage layout is to avoid cache conflicts.
 */
   cshkl = allocate_arrays(nsarray, hmax, kmax, lmax,
			   &chx, &cky, &clz, &shx, &sky, &slz, &chg,
			   &cos1x, &cos1y, &cos1z, &sin1x, &sin1y, &sin1z, 
			   &qcoskr, &qsinkr);
   memcp(chg, chgx+ns0, nsarr0*lsizeof(real));
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
	       intra += chgx[ssite+is] * chgx[ssite+js]
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
	 sq += chgx[is];
	 sqsq += SQR(chgx[is]);
      }
      self_energy += control.alpha / sqrt(PI) * sqsq;
      /*
       * Sqxf is total non-framework charge.  Calculate grand total in sq.
       */
      sqxf = sq;
      for(; is < nsites; is++)
	 sq += chgx[is];
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

      note("Ewald self-energy = %f Kj/mol",self_energy*CONV_E);
   }

   *pe -= self_energy;			/* Subtract self energy term	      */
   *pe += sheet_energy/vol;		/* Uniform charge correction	      */
   zero_real(stress_ew, 9);
   for(i=0; i<3; i++)
      stress_ew[i][i] += sheet_energy/vol;

   invert(system->h, hinvp);		/* Inverse of h is matrix of r.l.v.'s */
   mat_sca_mul(2*PI, hinvp, hinvp);
   /*
    * Build array hkl[] of k vectors within cutoff
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
   trig_recur(chx,shx,sin1x,cos1x,
	      site[0]+ns0,site[1]+ns0,site[2]+ns0,astar,hmax,0,nsarr0);
   trig_recur(cky,sky,sin1y,cos1y,
	      site[0]+ns0,site[1]+ns0,site[2]+ns0,bstar,kmax,0,nsarr0);
   trig_recur(clz,slz,sin1z,cos1z,
	      site[0]+ns0,site[1]+ns0,site[2]+ns0,cstar,lmax,0,nsarr0);
/*
 * Start of main loops over K vector indices h, k, l between -*max, *max etc.
 * To avoid calculating K and -K, only half of the K-space box is covered. 
 * Points on the axes are included once and only once. (0,0,0) is omitted.
 */
   for(phkl = hkl; phkl < hkl+nhkl; phkl++)
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
      qsincos(coshx,sinhx,cosky,sinky,coslz,sinlz,chg, 
	      qcoskr,qsinkr,k,l,nsarr0);
      ns1f = MIN(ns1, nsitesxf); ns0f = MAX(ns0, nsitesxf);
      sqexpkr[0] = sum(ns1f-ns0, qcoskr, 1);
      sqexpkr[1] = sum(ns1f-ns0, qsinkr, 1);
      sqexpkr[2] = sum(ns1-ns0f, qcoskr+(ns0f-ns0), 1);
      sqexpkr[3] = sum(ns1-ns0f, qsinkr+(ns0f-ns0), 1);
#ifdef SPMD
      par_rsum(sqexpkr, 4);
#endif
      sqcoskrn = sqexpkr[0]; 
      sqsinkrn = sqexpkr[1];
      sqcoskrf = sqexpkr[2];
      sqsinkrf = sqexpkr[3];
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
NOVECTOR
      for(i = 0; i < 3; i++)
      {
	 stress_ew[i][i] += pe_k;
NOVECTOR
	 for(j = i; j < 3; j++)
	    stress_ew[i][j] -= pe_k * coeff2 * kv[i] * kv[j];
      }
/*
 * Evaluation of site forces.   Non-framework sites interact with all others
 */
VECTORIZE
      for(is = 0; is < ns1f-ns0; is++)
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
      for(is = ns0f-ns0; is < ns1-ns0; is++)
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
   *pe /= nthreads;
   for(i=0; i<3; i++)
      for(j=0; j<3; j++)
	 stress[i][j] += stress_ew[i][j] / nthreads;

   afree((gptr*)chx); afree((gptr*)cky); afree((gptr*)clz); 
   afree((gptr*)shx); afree((gptr*)sky); afree((gptr*)slz);
   xfree(cshkl);
   xfree(hkl);
}
