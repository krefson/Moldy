/******************************************************************************
 * Ewald	The reciprocal-space part of the standard Ewald sum technique *
 ******************************************************************************
 *      Revision Log
 *       $Log:	ewald.c,v $
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
static char *RCSid = "$Header: /home/tigger/keith/md/RCS/ewald.c,v 1.5 89/10/26 11:29:26 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#ifdef	convexvc
#include <fastmath.h>
#else
#include <math.h>
#endif
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== Library declarations ============================*/
int	abs();
void	cfree();
/*========================== External function declarations ==================*/
double	err_fn();			/* Error function		      */
double	det();				/* Determinant of 3x3 matrix	      */
void	invert();			/* Inverts a 3x3 matrix		      */
void	mat_vec_mul();			/* Multiplies a 3x3 matrix by 3xN vect*/
double	sum();				/* Sum of elements of 'real' vector   */
void	saxpy();			/* A*x+y, x, y are long vectors	      */
void	*arralloc();			/* Allocates a dope vector array      */
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
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

/******************************************************************************
 *  Ewald  Calculate reciprocal-space part of coulombic forces		      *
 ******************************************************************************/
void	ewald(site,site_force,system,species,chg,pe,stress)
vec_t		site[],			/* Site co-ordinate arrays	 (in) */
		site_force[];		/* Site force arrays		(out) */
system_p	system;			/* System record		 (in) */
spec_t	species[];			/* Array of species records	 (in) */
real		chg[];			/* Array of site charges	 (in) */
double		*pe;			/* Potential energy		(out) */
mat_t		stress;			/* Stress virial		(out) */
{
   mat_t	hinv;			/* Matrix of reciprocal lattice vects*/
   register	int	h, k, l;	/* Recip. lattice vector indices     */
		int	i, j, is;	/* Counters.			     */
		int	ispec, ssite;	/* Species counter		     */
   		spec_p	spec;		/* species[ispec]		     */
   register	int	nsites = system->nsites;
   register	real	pe_k,		/* Pot'l energy for current K vector */
		        coeff, coeff2;	/* 2/(e0V) * A(K) & similar	     */
   register	real	sqcoskr,sqsinkr;/* Sum q(i) sin/cos(K.r(i))          */
		double	ksq;		/* Squared magnitude of K vector     */
		double	pe_local = 0.0;	/* Local accumulator for pot. energy */
   		vec_t	kv;		/* (Kx,Ky,Kz)  			     */
#ifdef OLDEWALD
   real		*chxky	= dalloc(nsites),
		*shxky	= dalloc(nsites);
#else
   		real	chxky, shxky;	/* Temporaries for sin(hx+ky) etc     */
#endif
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
   vec_t	*kr = ralloc(nsites);		/* K.Ri			      */
   real		*qcoskr,			/* q(i) cos(K.R(i))	      */
		*qsinkr;			/* q(i) sin(K.R(i))	      */
#ifdef VCALLS
   real		*temp	= dalloc(nsites);
#else
   real		force_comp;
#endif
   double	r_4_alpha = -1.0/(4.0 * control.alpha * control.alpha);
   double	vol = det(system->h);	/* Volume of MD cell		      */
   static	double	self_energy;	/* Constant self energy term	      */
   static	boolean init = true;	/* Flag for the first call of function*/

/*
 * First call only - evaluate self energy term and store for subsequent calls
 */
   if(init)
   {
      double	sqsq = 0, sq = 0, intra, r;
      int	js;
      self_energy = 0;
      for(is = 0; is < nsites; is++)
      {
	 sq += chg[is];
	 sqsq += SQR(chg[is]);
      }
      if( fabs(sq)*CONV_Q > 1.0e-5)
      {
	 self_energy += PI*SQR(sq) / (2.0*vol*SQR(control.alpha));
	 message(NULLI, NULLP, WARNING, SYSCHG, sq*CONV_Q, self_energy*CONV_E);
      }

      self_energy += control.alpha / sqrt(PI) * sqsq;
      ssite = 0;
      for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
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
      }
      note("Ewald self-energy = %f Kj/mol",self_energy*CONV_E);
      init = false;
   }
   *pe -= self_energy;			/* Subtract self energy term	      */
      
   invert(system->h, hinv);		/* Inverse of h is matrix of r.l.v.'s */

/*
 * Calculate cos and sin of astar*x, bstar*y & cstar*z for each charged site
 */
   for(is = 0; is < nsites; is++)
      chx[0][is] = cky[0][is] = clz[0][is] = 1.0;

   mat_vec_mul(hinv, site, kr, nsites);
VECTORIZE
   for(is = 0; is < nsites; is++)
   {
      chx[1][is] = cos(2.0 * PI * kr[is][0]);
      shx[1][is] = sin(2.0 * PI * kr[is][0]);
      cky[1][is] = cos(2.0 * PI * kr[is][1]);
      sky[1][is] = sin(2.0 * PI * kr[is][1]);
      clz[1][is] = cos(2.0 * PI * kr[is][2]);
      slz[1][is] = sin(2.0 * PI * kr[is][2]);
   }
/*
 *  Finished with kr[].  Re-assign space to qcoskr/qsinkr
 */
   qcoskr = kr[0];
   qsinkr = qcoskr + nsites;
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
      kv[0] = 2.0*PI*(h*astar[0] + k*bstar[0] + l*cstar[0]); 
      kv[1] = 2.0*PI*(h*astar[1] + k*bstar[1] + l*cstar[1]); 
      kv[2] = 2.0*PI*(h*astar[2] + k*bstar[2] + l*cstar[2]);
      
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
#endif
      sqcoskr = sum(nsites, qcoskr, 1);	 sqsinkr = sum(nsites, qsinkr, 1);
      
/*
 * Evaluate potential energy contribution for this K and add to total.
 */
      pe_k = 0.5 * coeff * (SQR(sqcoskr) + SQR(sqsinkr));
      pe_local += pe_k;

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
      saxpy(nsites, kv[0]*coeff, temp, 1, site_force[0], 3);
      saxpy(nsites, kv[1]*coeff, temp, 1, site_force[0]+1, 3);
      saxpy(nsites, kv[2]*coeff, temp, 1, site_force[0]+2, 3);
#else
/*
 * Evaluation of site forces. Vectorises under CRAY CC 4.0 & Convex VC 2.0
 */
VECTORIZE
      for(is = 0; is < nsites; is++)
      {
	 force_comp = coeff * (qsinkr[is]*sqcoskr - qcoskr[is]*sqsinkr);
	 site_force[is][0] += kv[0] * force_comp;
	 site_force[is][1] += kv[1] * force_comp;
	 site_force[is][2] += kv[2] * force_comp;
      }
#endif
   }

/*
 * End of loop over K vectors.  Add the Ewald potential energy to total.
 */
   *pe += pe_local;
   
   cfree((char*)chx[0]); cfree((char*)cky[0]); cfree((char*)clz[0]); 
   cfree((char*)shx[0]); cfree((char*)sky[0]); cfree((char*)slz[0]);
   cfree((char*)kr);
#ifdef OLDEWALD
   cfree((char*)chxky);	 cfree((char*)shxky); 
#endif
#ifdef VCALLS
   cfree((char*)temp);
#endif
}
