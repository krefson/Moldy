/******************************************************************************
 * Ewald	The reciprocal-space part of the standard Ewald sum technique *
 ******************************************************************************
 *      Revision Log
 *       $Log:	ewald_parallel.c,v $
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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/ewald_parallel.c,v 1.6 90/05/16 18:40:57 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#if  defined(convexvc) || defined(stellar)
#include <fastmath.h>
#else
#include <math.h>
#endif
#include "stdlib.h"
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
void	tfree();
double	err_fn();			/* Error function		      */
double	det();				/* Determinant of 3x3 matrix	      */
void	invert();			/* Inverts a 3x3 matrix		      */
void	mat_vec_mul();			/* Multiplies a 3x3 matrix by 3xN vect*/
void	transpose();			/* Transposes a 3x3 matrix	      */
void	mat_sca_mul();			/* Multiplies 3x3 matrix by scalar    */
double	sum();				/* Sum of elements of 'real' vector   */
#ifdef VCALLS
void	saxpy();			/* A*x+y, x, y are long vectors	      */
#endif
char	*arralloc();			/* Allocates a dope vector array      */
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
void	ewald_inner();			/* Inner loop forward reference       */
int	nprocessors();			/* Return no. of procs to execute on. */
/*========================== External data references ========================*/
extern	contr_t	control;		/* Main simulation control record     */
/*========================== Macros ==========================================*/
#define astar hinv[0]
#define bstar hinv[1]
#define cstar hinv[2]
#define moda(hmat) (hmat[0][0])
#define modb(hmat) sqrt(SQR(hmat[0][1]) + SQR(hmat[1][1]))
#define modc(hmat) sqrt(SQR(hmat[0][2]) + SQR(hmat[1][2]) + SQR(hmat[2][2]))
/*============================================================================*/
   struct _hkl {double kx, ky, kz; int h,k,l;};
/******************************************************************************
 *  Ewald  Calculate reciprocal-space part of coulombic forces		      *
 ******************************************************************************/
void	ewald(site,site_force,system,species,chg,pe,stress)
real		**site,			/* Site co-ordinate arrays	 (in) */
		**site_force;		/* Site force arrays		(out) */
system_p	system;			/* System record		 (in) */
spec_t	species[];			/* Array of species records	 (in) */
real		chg[];			/* Array of site charges	 (in) */
double		*pe;			/* Potential energy		(out) */
mat_t		stress;			/* Stress virial		(out) */
{
   mat_t	hinv;			/* Matrix of reciprocal lattice vects*/
   register	int	h, k, l;	/* Recip. lattice vector indices     */
		int	i, j, is, ssite;/* Counters.			     */
   		spec_p	spec;		/* species[ispec]		     */
   register	int	nsites = system->nsites;
   		double	kx,ky,kz;
   		vec_t	kv;		/* (Kx,Ky,Kz)  			     */
   	 	struct _hkl *hkl;
   		int	nhkl = 0;
   mat_t	htrinv;
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
   int		nthreads = nprocessors(),
   		ithread;
   double	*pe_n = aalloc(nthreads, double);
   mat_t	*stress_n = aalloc(nthreads, mat_t);
   real		***s_f_n;
   double	r_4_alpha = -1.0/(4.0 * control.alpha * control.alpha);
   double	vol = det(system->h);	/* Volume of MD cell		      */
   static	double	self_energy,	/* Constant self energy term	      */
   			sheet_energy;	/* Correction for non-neutral system. */
   static	boolean init = true;	/* Flag for the first call of function*/
   static	int	nsitesxf;	/* Number of non-framework sites.     */

   invert(system->h, hinv);		/* Inverse of h is matrix of r.l.v.'s */
   if (nthreads > 1)
       s_f_n = (real***)arralloc(sizeof(real), 3,
				 0, nthreads-2, 0, 2, 0, nsites-1);
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
	 sheet_energy += PI*(sq-sqxf)*sqxf / (2.0*SQR(control.alpha));
	 message(NULLI, NULLP, INFO, FRACHG, 
		 sqxf*CONV_Q, sheet_energy/vol*CONV_E);
      }
      /* 
       *  2) Case of entire system non-neutral.
       */
      if( fabs(sq)*CONV_Q > 1.0e-5)
      {
	 sheet_energy += intra = PI*SQR(sq) / (2.0*SQR(control.alpha));
	 message(NULLI, NULLP, WARNING, SYSCHG, sq*CONV_Q, intra/vol*CONV_E);
      }

      note("Ewald self-energy = %f Kj/mol",self_energy*CONV_E);
      init = false;
   }
   
   /*
    * Build array hkl[] of k vectors within cutoff
    */
   transpose(hinv, htrinv);
   mat_sca_mul(2*PI, htrinv, htrinv);
   hkl = aalloc(4*(hmax+1)*(kmax+1)*(lmax+1), struct _hkl);
   for(h = 0; h <= hmax; h++)
      for(k = (h==0 ? 0 : -kmax); k <= kmax; k++)
	 for(l = (h==0 && k==0 ? 1 : -lmax); l <= lmax; l++)
	 {
	    kv[0] = h; kv[1] = k; kv[2] = l;
	    mat_vec_mul(htrinv, (vec_t*)kv,(vec_t*)kv,1);
	    if( SUMSQ(kv) < SQR(control.k_cutoff) )
	    {
	       hkl[nhkl].h = h; hkl[nhkl].k = k; hkl[nhkl].l = l;
	       hkl[nhkl].kx = kv[0];
	       hkl[nhkl].ky = kv[1];
	       hkl[nhkl].kz = kv[2];
	       nhkl++;
	    }
	 }

   *pe -= self_energy+sheet_energy/vol;	/* Subtract self energy term	      */
      

/*
 * Calculate cos and sin of astar*x, bstar*y & cstar*z for each charged site
 */
VECTORIZE
   for(is = 0; is < nsites; is++)
      chx[0][is] = cky[0][is] = clz[0][is] = 1.0;

VECTORIZE
   for(is = 0; is < nsites; is++)
   {
      kx = hinv[0][0]*site[0][is]+hinv[0][1]*site[1][is]+hinv[0][2]*site[2][is];
      ky = hinv[1][0]*site[0][is]+hinv[1][1]*site[1][is]+hinv[1][2]*site[2][is];
      kz = hinv[2][0]*site[0][is]+hinv[2][1]*site[1][is]+hinv[2][2]*site[2][is];
      chx[1][is] = cos(2.0 * PI * kx);
      shx[1][is] = sin(2.0 * PI * kx);
      cky[1][is] = cos(2.0 * PI * ky);
      sky[1][is] = sin(2.0 * PI * ky);
      clz[1][is] = cos(2.0 * PI * kz);
      slz[1][is] = sin(2.0 * PI * kz);
   }
/*
 * Use addition formulae to get sin(h*astar*x)=sin(Kx*x) etc for each site
 */
   for(h = 2; h <= hmax; h++)
   {
      coshx = chx[h];
      sinhx = shx[h];
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 coshx[is] = chx[h-1][is]*chx[1][is] - shx[h-1][is]*shx[1][is];
	 sinhx[is] = shx[h-1][is]*chx[1][is] + chx[h-1][is]*shx[1][is];
      }
   }
   for(k = 2; k <= kmax; k++)
   {
      cosky = cky[k];
      sinky = sky[k];
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 cosky[is] = cky[k-1][is]*cky[1][is] - sky[k-1][is]*sky[1][is];
	 sinky[is] = sky[k-1][is]*cky[1][is] + cky[k-1][is]*sky[1][is];
      }
   }
   for(l = 2; l <= lmax; l++)
   {
      coslz = clz[l];
      sinlz = slz[l];
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 coslz[is] = clz[l-1][is]*clz[1][is] - slz[l-1][is]*slz[1][is];
	 sinlz[is] = slz[l-1][is]*clz[1][is] + clz[l-1][is]*slz[1][is];
      }
   }
/*
 * Start of main loops over K vector indices h, k, l between -*max, *max etc.
 * To avoid calculating K and -K, only half of the K-space box is covered. 
 * Points on the axes are included once and only once. (0,0,0) is omitted.
 */
/*$dir parallel*/
   for(ithread = 0; ithread < nthreads; ithread++)
      ewald_inner(ithread, nthreads, nhkl, hkl, nsites, nsitesxf, 
		  chx, cky, clz, shx, sky, slz, chg, vol, r_4_alpha,
		  stress_n[ithread], pe_n+ithread,
		  ithread ? s_f_n[ithread-1] : site_force);
/*
 *  Sum Pot, energies, forces and stress from each parallel invocation
 */
   for(ithread = 0; ithread < nthreads; ithread++)
   {
      *pe += pe_n[ithread];
      for(i = 0; i < 3; i++)
	 for(j = 0; j < 3; j++)
	    stress[i][j] += stress_n[ithread][i][j];
   }
   for(ithread = 0; ithread < nthreads-1; ithread++)
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 site_force[0][is] += s_f_n[ithread][0][is];
	 site_force[1][is] += s_f_n[ithread][1][is];
	 site_force[2][is] += s_f_n[ithread][2][is];
      }

   
   tfree((char*)chx); tfree((char*)cky); tfree((char*)clz); 
   tfree((char*)shx); tfree((char*)sky); tfree((char*)slz);
   tfree((char*)pe_n);   tfree((char*)stress_n);
   tfree((char*)hkl);
   if( nthreads > 1)
      tfree((char*)s_f_n);
}
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
   real chxky, shxky;
   
   if( k >= 0 )
      if( l >= 0 )
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    chxky = coshx[is]*cosky[is] - sinhx[is]*sinky[is];
	    shxky = sinhx[is]*cosky[is] + coshx[is]*sinky[is];
	    qcoskr[is] = chg[is]*(chxky*coslz[is] - shxky*sinlz[is]);
	    qsinkr[is] = chg[is]*(shxky*coslz[is] + chxky*sinlz[is]);
	 }
      else
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    chxky = coshx[is]*cosky[is] - sinhx[is]*sinky[is];
	    shxky = sinhx[is]*cosky[is] + coshx[is]*sinky[is];
	    qcoskr[is] = chg[is]*(chxky*coslz[is] +shxky*sinlz[is]);
	    qsinkr[is] = chg[is]*(shxky*coslz[is] - chxky*sinlz[is]);
	 }
   else
      if( l >= 0 )
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    chxky = coshx[is]*cosky[is] + sinhx[is]*sinky[is];
	    shxky = sinhx[is]*cosky[is] - coshx[is]*sinky[is];
	    qcoskr[is] = chg[is]*(chxky*coslz[is] - shxky*sinlz[is]);
	    qsinkr[is] = chg[is]*(shxky*coslz[is] + chxky*sinlz[is]);
	 }
      else
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    chxky = coshx[is]*cosky[is] + sinhx[is]*sinky[is];
	    shxky = sinhx[is]*cosky[is] - coshx[is]*sinky[is];
	    qcoskr[is] = chg[is]*(chxky*coslz[is] + shxky*sinlz[is]);
	    qsinkr[is] = chg[is]*(shxky*coslz[is] - chxky*sinlz[is]);
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
	    chg, vol, r_4_alpha, stress, pe, site_force)
int ithread, nthreads, nhkl;
struct _hkl hkl[];
int nsites;
int nsitesxf;			/* N sites excluding framework sites.	      */
real **chx, **cky, **clz, **shx, **sky, **slz;
double vol, r_4_alpha;
mat_t	stress;
double *pe;
real	chg[];
real	**site_force;
{
   vec_t	kv;
   double ksq, coeff, coeff2, pe_k;
   double sqcoskr, sqsinkr, sqcoskrn, sqsinkrn, sqcoskrf, sqsinkrf;
   real *coshx, *cosky, *coslz, *sinhx, *sinky, *sinlz;
   int is, i, j, h, k, l, ihkl;
#ifdef OLDEWALD
   real		*chxky	= dalloc(nsites),
		*shxky	= dalloc(nsites);
#else
#endif
#ifdef VCALLS
   real		*temp	= dalloc(nsites);
#else
   real		force_comp;
#endif
   real *qcoskr = dalloc(nsites), *qsinkr = dalloc(nsites);
   real		*site_fx = site_force[0],
   		*site_fy = site_force[1],
   		*site_fz = site_force[2];

   for(ihkl = ithread; ihkl < nhkl; ihkl += nthreads)
   {
      h  = hkl[ihkl].h;	    k     = hkl[ihkl].k;  l     = hkl[ihkl].l;
      kv[0] = hkl[ihkl].kx; kv[1] = hkl[ihkl].ky; kv[2] = hkl[ihkl].kz;
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
#ifdef OLDEWALD
      if(k >= 0)
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    chxky[is] = coshx[is]*cosky[is] - sinhx[is]*sinky[is];
	    shxky[is] = sinhx[is]*cosky[is] + coshx[is]*sinky[is];
	 }
      else
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    chxky[is] = coshx[is]*cosky[is] + sinhx[is]*sinky[is];
	    shxky[is] = sinhx[is]*cosky[is] - coshx[is]*sinky[is];
	 }

      if(l >= 0)
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qcoskr[is] = chg[is]*(chxky[is]*coslz[is] - shxky[is]*sinlz[is]);
	    qsinkr[is] = chg[is]*(shxky[is]*coslz[is] + chxky[is]*sinlz[is]);
	 }
      else
VECTORIZE
	 for(is = 0; is < nsites; is++)
	 {
	    qcoskr[is] = chg[is]*(chxky[is]*coslz[is] + shxky[is]*sinlz[is]);
	    qsinkr[is] = chg[is]*(shxky[is]*coslz[is] - chxky[is]*sinlz[is]);
	 }
#else			/* Use of scalar temps is faster if compiler up to it*/
      qsincos(coshx,sinhx,cosky,sinky,coslz,sinlz,chg,
	      qcoskr,qsinkr,k,l,nsites);
#endif
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
/*
 * Calculate long-range coulombic contribution to stress tensor
 */
      for(i = 0; i < 3; i++)
      {
	 stress[i][i] += pe_k;
NOVECTOR
	 for(j = i; j < 3; j++)
	    stress[i][j] -= pe_k * coeff2 * kv[i] * kv[j];
      }

/*
 * Old and less efficient calculation of site forces. Retained for machines
 * with poor vectorising compilers.
 */
#ifdef VCALLS
VECTORIZE
      for(is = 0; is < nsites; is++)
	 temp[is] = qsinkr[is]*sqcoskr - qcoskr[is]*sqsinkr;
      saxpy(nsites, kv[0], temp, 1, site_force[0], 1);
      saxpy(nsites, kv[1], temp, 1, site_force[1], 1);
      saxpy(nsites, kv[2], temp, 1, site_force[2], 1);
#else
/*
 * Evaluation of site forces. Vectorises under CRAY CC 4.0 & Convex VC 2.0
 */
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 force_comp = qsinkr[is]*sqcoskr - qcoskr[is]*sqsinkr;
	 site_fx[is] += kv[0] * force_comp;
	 site_fy[is] += kv[1] * force_comp;
	 site_fz[is] += kv[2] * force_comp;
      }
#endif
/*
 * End of loop over K vectors.
 */
    }
tfree((char*)qcoskr); tfree((char*)qsinkr);
#ifdef OLDEWALD
   tfree((char*)chxky);	 tfree((char*)shxky); 
#endif
#ifdef VCALLS
   tfree((char*)temp);
#endif
}

char *getenv();
int nprocessors()
{
   char *env;
   int n;
   if( ( env = getenv("THREADS") ) == NULL || (n = atoi(env)) <= 0 )
      return 4;
   else
      return n;
}
