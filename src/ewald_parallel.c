/******************************************************************************
 * Ewald	The reciprocal-space part of the standard Ewald sum technique *
 ******************************************************************************
 *      Revision Log
 *       $Log:	ewald_parallel.c,v $
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
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
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
 * Removed references to size_t and time_t typedefs, no longer in "defs.h"
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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/ewald_parallel.c,v 1.21 93/03/09 15:59:58 keith Exp $";
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
/*========================== Program include files ===========================*/
#include 	"structs.h"
#include 	"messages.h"
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
double	err_fn();			/* Error function		      */
double	det();				/* Determinant of 3x3 matrix	      */
void	invert();			/* Inverts a 3x3 matrix		      */
void	mat_vec_mul();			/* Multiplies a 3x3 matrix by 3xN vect*/
void	mat_sca_mul();			/* Multiplies a 3x3 matrix by scalar  */
void	transpose();			/* Transposes a 3x3 matrix	      */
double	sum();				/* Sum of elements of 'real' vector   */
gptr	*arralloc();			/* Allocates a dope vector array      */
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
void	ewald_inner();			/* Inner loop forward reference       */
int	nprocessors();			/* Return no. of procs to execute on. */
/*========================== External data references ========================*/
extern	contr_mt	control;		/* Main simulation control record     */
/*========================== Macros ==========================================*/
#define astar hinvp[0]
#define bstar hinvp[1]
#define cstar hinvp[2]
#define moda(hmat) (hmat[0][0])
#define modb(hmat) sqrt(SQR(hmat[0][1]) + SQR(hmat[1][1]))
#define modc(hmat) sqrt(SQR(hmat[0][2]) + SQR(hmat[1][2]) + SQR(hmat[2][2]))
/*============================================================================*/
   struct _hkl {double kx, ky, kz; int h,k,l;};
/******************************************************************************
 *  Ewald  Calculate reciprocal-space part of coulombic forces		      *
 ******************************************************************************/
void ewald_inner();
#ifdef titan
#ifdef PARALLEL
#pragma opt_level 3
#pragma pproc ewald_inner
#endif
#endif
void	ewald(site,site_force,system,species,chg,pe,stress)
real		**site,			/* Site co-ordinate arrays	 (in) */
		**site_force;		/* Site force arrays		(out) */
system_mp	system;			/* System record		 (in) */
spec_mt	species[];			/* Array of species records	 (in) */
real		chg[];			/* Array of site charges	 (in) */
double		*pe;			/* Potential energy		(out) */
mat_mt		stress;			/* Stress virial		(out) */
{
   mat_mt	hinvp;			/* Matrix of reciprocal lattice vects*/
   register	int	h, k, l;	/* Recip. lattice vector indices     */
		int	i, j, is, ssite;/* Counters.			     */
   		spec_mp	spec;		/* species[ispec]		     */
   register	int	nsites = system->nsites;
   		double	kx,ky,kz;
   		vec_mt	kv;		/* (Kx,Ky,Kz)  			     */
   	 	struct _hkl *hkl;
   		int	nhkl = 0;
/*
 * Maximum values of h, k, l  s.t. |k| < k_cutoff
 */
		int	hmax = floor(control.k_cutoff/(2*PI)*moda(system->h)),
			kmax = floor(control.k_cutoff/(2*PI)*modb(system->h)),
			lmax = floor(control.k_cutoff/(2*PI)*modc(system->h));
/*
 * Arrays for cos & sin (h x(i)), (k y(i)) and (l z(i)) eg chx[h][isite]
 * and pointers to a particular h,k or l eg coshx[is] = chh[2][is]
 */
   real		**chx = (real**)arralloc(sizeof(real),2, 0, hmax, 0, nsites-1),
		**cky = (real**)arralloc(sizeof(real),2, 0, kmax, 0, nsites-1),
		**clz = (real**)arralloc(sizeof(real),2, 0, lmax, 0, nsites-1),
		**shx = (real**)arralloc(sizeof(real),2, 0, hmax, 0, nsites-1),
		**sky = (real**)arralloc(sizeof(real),2, 0, kmax, 0, nsites-1),
		**slz = (real**)arralloc(sizeof(real),2, 0, lmax, 0, nsites-1);
   real		*coshx, *cosky, *coslz, *sinhx, *sinky, *sinlz;
   real               *c1, *s1, *cm1, *sm1;
   real		*sf0, *sf1, *sf2, *ssf0, *ssf1, *ssf2;
   real		*site0, *site1, *site2;
   int		nthreads = nprocessors(),
   		ithread;
   double	*pe_n = aalloc(nthreads, double);
   mat_mt	*stress_n = aalloc(nthreads, mat_mt);
   real		***s_f_n;
   double	r_4_alpha = -1.0/(4.0 * control.alpha * control.alpha);
   double	vol = det(system->h);	/* Volume of MD cell		      */
   static	double	self_energy,	/* Constant self energy term	      */
   			sheet_energy;	/* Correction for non-neutral system. */
   static	boolean init = true;	/* Flag for the first call of function*/
   static	int	nsitesxf;	/* Number of non-framework sites.     */

   invert(system->h, hinvp);		/* Inverse of h is matrix of r.l.v.'s */
   mat_sca_mul(2*PI, hinvp, hinvp);
   s_f_n = aalloc(nthreads, real**);
   s_f_n[0] = site_force;
   for(ithread = 1; ithread < nthreads; ithread++)
   {
      s_f_n[ithread] = (real**)arralloc(sizeof(real), 2, 0, 2, 0, nsites-1);
      zero_real(s_f_n[ithread][0],3*nsites);
   }
   zero_real(stress_n,9*nthreads);
   zero_double(pe_n, nthreads);

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

#ifdef titan
#pragma no_parallel
#endif
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
#ifdef titan
#pragma no_parallel
#endif
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

      note("Ewald self-energy = %f Kj/mol",self_energy*CONV_E);
      init = false;
   }
   
   /*
    * Build array hkl[] of k vectors within cutoff
    */
   hkl = aalloc(4*(hmax+1)*(kmax+1)*(lmax+1), struct _hkl);
   for(h = 0; h <= hmax; h++)
      for(k = (h==0 ? 0 : -kmax); k <= kmax; k++)
      {
	 kv[0] = h*astar[0] + k*bstar[0];
	 kv[1] = h*astar[1] + k*bstar[1];
	 kz = h*astar[2] + k*bstar[2];
	 for(l = (h==0 && k==0 ? 1 : -lmax); l <= lmax; l++)
	 {
	  /*  kv[0] = kx + l*cstar[0]; 
	    kv[1] = ky + l*cstar[1]; */
	    kv[2] = kz + l*cstar[2];
	    if( SUMSQ(kv) < SQR(control.k_cutoff) )
	    {
	       hkl[nhkl].h = h; hkl[nhkl].k = k; hkl[nhkl].l = l;
	       hkl[nhkl].kx = kv[0];
	       hkl[nhkl].ky = kv[1];
	       hkl[nhkl].kz = kv[2];
	       nhkl++;
	    }
	 }
      }

   *pe -= self_energy;			/* Subtract self energy term	      */
   *pe += sheet_energy/vol;		/* Uniform charge correction	      */
   for(i=0; i<3; i++)
      stress[i][i] += sheet_energy/vol;
/*
 * Calculate cos and sin of astar*x, bstar*y & cstar*z for each charged site
 */
   coshx = chx[0]; cosky = cky[0]; coslz = clz[0];
   sinhx = shx[0]; sinky = sky[0]; sinlz = slz[0];
#ifdef titan
#pragma no_parallel
#endif
VECTORIZE
   for(is = 0; is < nsites; is++)
   {
      coshx[is] = cosky[is] = coslz[is] = 1.0;
      sinhx[is] = sinky[is] = sinlz[is] = 0.0;
   }      

   coshx = chx[1]; cosky = cky[1]; coslz = clz[1];
   sinhx = shx[1]; sinky = sky[1]; sinlz = slz[1];
   site0 = site[0]; site1 = site[1]; site2 = site[2];
#ifdef titan
#pragma no_parallel
#endif
VECTORIZE
   for(is = 0; is < nsites; is++)
   {
      kx = astar[0]*site0[is]+astar[1]*site1[is]+astar[2]*site2[is];
      ky = bstar[0]*site0[is]+bstar[1]*site1[is]+bstar[2]*site2[is];
      kz = cstar[0]*site0[is]+cstar[1]*site1[is]+cstar[2]*site2[is];
      coshx[is] = cos(kx);
      sinhx[is] = sin(kx);
      cosky[is] = cos(ky);
      sinky[is] = sin(ky);
      coslz[is] = cos(kz);
      sinlz[is] = sin(kz);
   }
/*
 * Use addition formulae to get sin(h*astar*x)=sin(Kx*x) etc for each site
 */
   for(h = 2; h <= hmax; h++)
   {
      coshx = chx[h];
      sinhx = shx[h];
      cm1 = chx[h-1]; sm1 = shx[h-1];
      c1  = chx[1];   s1  = shx[1];
#ifdef titan
#pragma no_parallel
#endif
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 coshx[is] = cm1[is]*c1[is] - sm1[is]*s1[is];
	 sinhx[is] = sm1[is]*c1[is] + cm1[is]*s1[is];
      }
   }
   for(k = 2; k <= kmax; k++)
   {
      cosky = cky[k];
      sinky = sky[k];
      cm1 = cky[k-1]; sm1 = sky[k-1];
      c1  = cky[1];   s1  = sky[1];
#ifdef titan
#pragma no_parallel
#endif
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 cosky[is] = cm1[is]*c1[is] - sm1[is]*s1[is];
	 sinky[is] = sm1[is]*c1[is] + cm1[is]*s1[is];
      }
   }
   for(l = 2; l <= lmax; l++)
   {
      coslz = clz[l];
      sinlz = slz[l];
      cm1 = clz[l-1]; sm1 = slz[l-1];
      c1  = clz[1];   s1  = slz[1];
#ifdef titan
#pragma no_parallel
#endif
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 coslz[is] = cm1[is]*c1[is] - sm1[is]*s1[is];
	 sinlz[is] = sm1[is]*c1[is] + cm1[is]*s1[is];
      }
   }
/*
 * Start of main loops over K vector indices h, k, l between -*max, *max etc.
 * To avoid calculating K and -K, only half of the K-space box is covered. 
 * Points on the axes are included once and only once. (0,0,0) is omitted.
 */
#ifdef PARALLEL
#ifdef stellar
/*$dir parallel*/
#endif /* stellar */
#ifdef titan
#pragma ipdep
#pragma pproc ewald_inner
#endif titan
#ifdef __convexc__
#pragma _CNX force_parallel
#endif /* --convexc__ */
#ifdef CRAY
#pragma _CRI taskloop private(ithread) shared(nthreads, nhkl, hkl, nsites, \
		  nsitesxf, chx, cky, clz, shx, sky, slz, chg, vol, r_4_alpha, \
		  stress_n, pe_n, s_f_n)
#endif /* CRAY */
#endif /*PARALLEL */
   for(ithread = 0; ithread < nthreads; ithread++)
      ewald_inner(ithread, nthreads, nhkl, hkl, nsites, nsitesxf, 
		  chx, cky, clz, shx, sky, slz, chg, &vol, &r_4_alpha,
		  stress_n[ithread], pe_n+ithread, s_f_n[ithread]);
/*
 *  Sum Pot, energies, forces and stress from each parallel invocation
 */
   for(ithread = 0; ithread < nthreads; ithread++)
   {
      *pe += pe_n[ithread];
#ifdef titan
#pragma asis
#endif
      for(i = 0; i < 3; i++)
	 for(j = 0; j < 3; j++)
	    stress[i][j] += stress_n[ithread][i][j];
   }
   sf0 = site_force[0]; sf1 = site_force[1]; sf2 = site_force[2];
   for(ithread = 1; ithread < nthreads; ithread++)
   {
      ssf0 = s_f_n[ithread][0];
      ssf1 = s_f_n[ithread][1];
      ssf2 = s_f_n[ithread][2];
#ifdef titan
#pragma ipdep
#endif
#ifdef __convexc__
#pragma _CNX vstrip (64)
#pragma _CNX force_vector
#pragma _CNX force_parallel_ext
#endif
#ifdef CRAY
#pragma _CRI ivdep
#endif
      for(is = 0; is < nsites; is++)
      {
	 sf0[is] += ssf0[is];
	 sf1[is] += ssf1[is];
	 sf2[is] += ssf2[is];
      }
   }
   
   xfree(chx); xfree(cky); xfree(clz); 
   xfree(shx); xfree(sky); xfree(slz);
   xfree(pe_n);   xfree(stress_n);
   xfree(hkl);
   for( ithread = 1; ithread < nthreads; ithread++)
      xfree(s_f_n[ithread]);
   xfree(s_f_n);
}
#ifdef titan
#ifdef PARALLEL
#pragma opt_level 2
#endif
#endif
/*****************************************************************************
 * qsincos().  Evaluate q sin(k.r) and q cos(k.r).  This is in a separate    *
 * function because some compilers (notably Stellar's) generate MUCH better  *
 * vector code this way. 						     *
 *****************************************************************************/
void      qsincos(coshx,sinhx,cosky,sinky,coslz,sinlz,chg,
		  qcoskr,qsinkr,k,l,nsites)
real coshx[], sinhx[], cosky[], sinky[], coslz[], sinlz[],
     chg[], qcoskr[], qsinkr[];
int  k,l,nsites;
{
   int is;
   
   if( k >= 0 )
      if( l >= 0 )
      {
#ifdef __convexc__
#pragma _CNX no_parallel
#endif
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qcoskr[is] = chg[is]*(
		  (coshx[is]*cosky[is] - sinhx[is]*sinky[is])*coslz[is] 
                - (sinhx[is]*cosky[is] + coshx[is]*sinky[is])*sinlz[is]);
	    qsinkr[is] = chg[is]*(
                  (sinhx[is]*cosky[is] + coshx[is]*sinky[is])*coslz[is] 
		+ (coshx[is]*cosky[is] - sinhx[is]*sinky[is])*sinlz[is]);
	 }
      }
      else
      {
#ifdef __convexc__
#pragma _CNX no_parallel
#endif
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qcoskr[is] = chg[is]*(
		  (coshx[is]*cosky[is] - sinhx[is]*sinky[is])*coslz[is] 
                + (sinhx[is]*cosky[is] + coshx[is]*sinky[is])*sinlz[is]);
	    qsinkr[is] = chg[is]*(
                  (sinhx[is]*cosky[is] + coshx[is]*sinky[is])*coslz[is] 
		- (coshx[is]*cosky[is] - sinhx[is]*sinky[is])*sinlz[is]);
	 }
      }
   else
      if( l >= 0 )
      {
#ifdef __convexc__
#pragma _CNX no_parallel
#endif
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qcoskr[is] = chg[is]*(
		  (coshx[is]*cosky[is] + sinhx[is]*sinky[is])*coslz[is] 
                - (sinhx[is]*cosky[is] - coshx[is]*sinky[is])*sinlz[is]);
	    qsinkr[is] = chg[is]*(
                  (sinhx[is]*cosky[is] - coshx[is]*sinky[is])*coslz[is] 
		+ (coshx[is]*cosky[is] + sinhx[is]*sinky[is])*sinlz[is]);
	 }
      }
      else
      {
#ifdef __convexc__
#pragma _CNX no_parallel
#endif
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qcoskr[is] = chg[is]*(
		  (coshx[is]*cosky[is] + sinhx[is]*sinky[is])*coslz[is] 
                + (sinhx[is]*cosky[is] - coshx[is]*sinky[is])*sinlz[is]);
	    qsinkr[is] = chg[is]*(
                  (sinhx[is]*cosky[is] - coshx[is]*sinky[is])*coslz[is] 
		- (coshx[is]*cosky[is] + sinhx[is]*sinky[is])*sinlz[is]);
	 }
     }
}
/*****************************************************************************
 *  Ewald_inner().  Part of Ewald sum to run in parallel on multi-stream or  *
 *  multi-processor computers.  It splits up the loop over k-vectors by using*
 *  a stride of the number of threads in use.  The loop starts from a value  *
 *  unique to the thread.    Summable quantities, pe, stress, site_force are *
 *  accumulated where specified by the arguments so it is the callers	     *
 *  responsibility to provide a separte area and accumulate the grand totals *
 *****************************************************************************/
void
ewald_inner(ithread, nthreads, nhkl, hkl, nsites, nsitesxf, 
            chx, cky, clz, shx, sky, slz,
	    chg, volp, r_4_alphap, stress, pe, site_force)
int ithread, nthreads, nhkl;
struct _hkl hkl[];
int nsites;
int nsitesxf;			/* N sites excluding framework sites.	      */
real **chx, **cky, **clz, **shx, **sky, **slz;
double *volp, *r_4_alphap;
mat_mt	stress;
double *pe;
real	chg[];
real	**site_force;
{
   vec_mt	kv;
   double ksq, coeff, coeff2, pe_k;
   double sqcoskr, sqsinkr, sqcoskrn, sqsinkrn, sqcoskrf, sqsinkrf;
   real *coshx, *cosky, *coslz, *sinhx, *sinky, *sinlz;
   int is, i, j, h, k, l;
   struct _hkl *phkl;
   double vol = *volp, r_4_alpha = *r_4_alphap;
   real		force_comp;
   real *qcoskr = dalloc(nsites), *qsinkr = dalloc(nsites);
   real		*site_fx = site_force[0],
   		*site_fy = site_force[1],
   		*site_fz = site_force[2];

   for(phkl = hkl+ithread; phkl < hkl+nhkl; phkl += nthreads)
   {
      h  = phkl->h;	    k     = phkl->k;  l     = phkl->l;
      kv[0] = phkl->kx; kv[1] = phkl->ky; kv[2] = phkl->kz;
/*
 * Calculate pre-factors A(K) etc
 */
      ksq = SUMSQ(kv);
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
	      qcoskr,qsinkr,k,l,nsites);
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
NOVECTOR
      for(i = 0; i < 3; i++)
      {
	 stress[i][i] += pe_k;
NOVECTOR
	 for(j = i; j < 3; j++)
	    stress[i][j] -= pe_k * coeff2 * kv[i] * kv[j];
      }

/*
 * Evaluation of site forces.   Non-framework sites interact with all others
 */
VECTORIZE
      for(is = 0; is < nsitesxf; is++)
      {
	 force_comp = qsinkr[is]*sqcoskr - qcoskr[is]*sqsinkr;
	 site_fx[is] += kv[0] * force_comp;
	 site_fy[is] += kv[1] * force_comp;
	 site_fz[is] += kv[2] * force_comp;
      }
#if 1
/*
 *  Framework sites -- only interact with non-framework sites
 */
VECTORIZE
      for(is = nsitesxf; is < nsites; is++)
      {
	 force_comp = qsinkr[is]*sqcoskrn - qcoskr[is]*sqsinkrn;
	 site_fx[is] += kv[0] * force_comp;
	 site_fy[is] += kv[1] * force_comp;
	 site_fz[is] += kv[2] * force_comp;
      }
#endif
/*
 * End of loop over K vectors.
 */
   }
xfree(qcoskr); xfree(qsinkr);
}

char *getenv();
#ifdef titan
int nprocessors()
{
   char *env;
   static int n=0;
   int    nphys;

   if( n <= 0 )
   {
      nphys = MT_NUMBER_OF_PROCS();
      if( ( env = getenv("THREADS") ) == NULL )
	 n = 1;
      else
      {
	 n = atoi(env);
	 if ( n <= 0 || n > nphys)
	    n = nphys;
      }
   }
   return n;
}
#else
#ifdef CRAY
int nprocessors()
{
   char *env;
   static int n = 0;
   int nphys = 8;

   if( n <= 0 )
   {
      if( ( env = getenv("NCPUS") ) == NULL )
	 n = 1;
      else
      {
	 n = atoi(env);
	 if ( n <= 0 || n > nphys)
	    n = nphys;
      }
   }
   return n;
}
#else			/* GS1000/2000 but should compile on any unix */
int nprocessors()
{
   char *env;
   static int n = 0;
   int nphys = 4;

   if( n <= 0 )
   {
      if( ( env = getenv("THREADS") ) == NULL )
	 n = 1;
      else
      {
	 n = atoi(env);
	 if ( n <= 0 || n > nphys)
	    n = nphys;
      }
   }
   return n;
}
#endif
#endif
