/******************************************************************************
 * Ewald	The reciprocal-space part of the standard Ewald sum technique *
 ******************************************************************************
 *      Revision Log
 *       $Log:	ewald.c,v $
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
static char *RCSid = "$Header: /mnt/keith/moldy/RCS/ewald.c,v 1.14 90/08/22 10:32:44 keith Exp $";
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
void	mat_sca_mul();			/* Multiplies a 3x3 matrix by scalar  */
double	sum();				/* Sum of elements of 'real' vector   */
#ifdef VCALLS
void	saxpy();			/* A*x+y, x, y are long vectors	      */
#endif
char	*arralloc();			/* Allocates a dope vector array      */
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
/*========================== External data references ========================*/
extern	contr_t	control;		/* Main simulation control record     */
/*========================== Macros ==========================================*/
#define astar hinvp[0]
#define bstar hinvp[1]
#define cstar hinvp[2]
#define moda(hmat) (hmat[0][0])
#define modb(hmat) sqrt(SQR(hmat[0][1]) + SQR(hmat[1][1]))
#define modc(hmat) sqrt(SQR(hmat[0][2]) + SQR(hmat[1][2]) + SQR(hmat[2][2]))
/*============================================================================*/
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
   mat_t	hinvp;			/* Matrix of reciprocal lattice vects*/
   register	int	h, k, l;	/* Recip. lattice vector indices     */
		int	i, j, is, ssite;/* Counters.			     */
   		spec_p	spec;		/* species[ispec]		     */
   register	int	nsites = system->nsites;
   register	real	pe_k,		/* Pot'l energy for current K vector */
		        coeff, coeff2;	/* 2/(e0V) * A(K) & similar	     */
   register	real	sqcoskr,sqsinkr,/* Sum q(i) sin/cos(K.r(i))          */
   			sqcoskrn, sqsinkrn,
			sqcoskrf, sqsinkrf;
		double	ksq;		/* Squared magnitude of K vector     */
   		double	kx,ky,kz;
   		vec_t	kv;		/* (Kx,Ky,Kz)  			     */
#ifdef OLDEWALD
   real		*chxky	= dalloc(nsites),
		*shxky	= dalloc(nsites);
#endif
   real		*site_fx = site_force[0],
   		*site_fy = site_force[1],
   		*site_fz = site_force[2];
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
   real		*qcoskr = dalloc(nsites),	/* q(i) cos(K.R(i))	      */
		*qsinkr = dalloc(nsites);	/* q(i) sin(K.R(i))	      */
#ifdef VCALLS
   real		*temp	= dalloc(nsites);
#else
   real		force_comp;
#endif
   double	r_4_alpha = -1.0/(4.0 * control.alpha * control.alpha);
   double	vol = det(system->h);	/* Volume of MD cell		      */
   static	double	self_energy,	/* Constant self energy term	      */
   			sheet_energy;	/* Correction for non-neutral system. */
   static	boolean init = true;	/* Flag for the first call of function*/
   static	int	nsitesxf;	/* Number of non-framework sites.     */

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

   *pe -= self_energy+sheet_energy/vol;	/* Subtract self energy term	      */

   invert(system->h, hinvp);		/* Inverse of h is matrix of r.l.v.'s */
   mat_sca_mul(2*PI, hinvp, hinvp);

/*
 * Calculate cos and sin of astar*x, bstar*y & cstar*z for each charged site
 */
   for(is = 0; is < nsites; is++)
      chx[0][is] = cky[0][is] = clz[0][is] = 1.0;

VECTORIZE
   for(is = 0; is < nsites; is++)
   {
      kx = astar[0]*site[0][is]+astar[1]*site[1][is]+astar[2]*site[2][is];
      ky = bstar[0]*site[0][is]+bstar[1]*site[1][is]+bstar[2]*site[2][is];
      kz = cstar[0]*site[0][is]+cstar[1]*site[1][is]+cstar[2]*site[2][is];
      chx[1][is] = cos(kx);
      shx[1][is] = sin(kx);
      cky[1][is] = cos(ky);
      sky[1][is] = sin(ky);
      clz[1][is] = cos(kz);
      slz[1][is] = sin(kz);
   }
/*
 * Use addition formulae to get sin(h*astar*x)=sin(Kx*x) etc for each site
 */
   for(h = 2; h <= hmax; h++)
   {
      coshx = chx[h];
      sinhx = shx[h];
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
   for(h = 0; h <= hmax; h++)
   for(k = (h==0 ? 0 : -kmax); k <= kmax; k++)
   for(l = (h==0 && k==0 ? 1 : -lmax); l <= lmax; l++)
   {
/*
 * Calculate actual K vector and its squared magnitude.
 */
      kv[0] = h*astar[0] + k*bstar[0] + l*cstar[0]; 
      kv[1] = h*astar[1] + k*bstar[1] + l*cstar[1]; 
      kv[2] = h*astar[2] + k*bstar[2] + l*cstar[2];
      
      ksq = SUMSQ(kv);
      
/*
 * Test whether K is within the specified cut-off and skip rest of loop if not
 */
      if(ksq >= SQR(control.k_cutoff)) continue;
      
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
   }

/*
 * End of loop over K vectors.
 */
   tfree((char*)chx); tfree((char*)cky); tfree((char*)clz); 
   tfree((char*)shx); tfree((char*)sky); tfree((char*)slz);
   tfree((char*)qcoskr); tfree((char*)qsinkr);
#ifdef OLDEWALD
   tfree((char*)chxky);	 tfree((char*)shxky); 
#endif
#ifdef VCALLS
   tfree((char*)temp);
#endif
}
