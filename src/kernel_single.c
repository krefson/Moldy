#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/kernel_single.c,v 1.3 90/09/28 13:29:41 keith Exp $";
#endif

/*
 * $Log:	kernel_single.c,v $
 * Revision 1.3  90/09/28  13:29:41  keith
 * Inserted braces around VECTORIZE directives and changed include files
 * for STARDtardent 3000 series (via cond. comp symbol "ardent").
 * 
 * Revision 1.2  89/07/04  18:43:11  keith
 * Fixed error in kernel and force which led to sites being allocated the
 * wrong potential parameters.  Needed extra parameter to kernel.
 * 
 * Revision 1.1  89/04/20  16:01:53  keith
 * Initial revision
 * 
 * Revision 1.1  89/04/11  15:05:08  keith
 * Initial revision
 * 
 */

#if  defined(convexvc) || defined(stellar)
#   include <fastmath.h>
#else
#ifdef ardent
#   include <vmath.h>
#else
#   include <math.h>
#endif
#endif
#include "defs.h"
#include "messages.h"
#define a0 -1.26551223
#define a1  1.00002368
#define a2  0.37409196
#define a3  0.09678418
#define a4 -0.18628806
#define a5  0.27886807
#define a6 -1.13520398
#define a7  1.48851587
#define a8 -0.82215223
#define a9  0.17087277

#define LJPOT 0
#define E6POT 1
#define MCYPOT 2

#define	POLY(t)	(a0+(t)*(a1+(t)*(a2+(t)*(a3+(t)*(a4+(t)*(a5+(t)*(a6+\
                    (t)*(a7+(t)*(a8+(t)*a9)))))))))

void message();

#define SINGLE
#if defined(SINGLE) && defined(sun)	/* Single precision functions for SUN	      */
static float tmp;
static union {FLOATFUNCTIONTYPE i; float f;} res;
#define sqrt(x) (res.i = r_sqrt_(&x), res.f)
#define exp(x)  (res.i = r_exp_((tmp=(x),&tmp)), res.f)
#endif
#if defined(SINGLE) && defined(convex)
#define sqrt(x) ssqrt(x)
#define exp(x)  sexp(x)
#endif

/******************************************************************************
 *  dist_pot   return attractive part of potential integrated from cutoff to  *
 *  infinity for distant pressure calculation.				      *
 ******************************************************************************/
double	dist_pot(potpar, cutoff, ptype)
real	potpar[];			/* Array of potential parameters      */
double	cutoff;				/* Cutoff distance		      */
int	ptype;				/* Potential type selector	      */
{
   switch(ptype)
   {
    default:
      message(NULLI, NULLP, FATAL, UNKPTY, ptype);
    case LJPOT:
      return(potpar[0]*CUBE(SQR(potpar[1])/cutoff) / 3.0);
    case E6POT:
      return(potpar[0] / ( 3.0*CUBE(cutoff)));
    case MCYPOT:
      if( potpar[2] != 0.0 )
         return( potpar[2] * (SQR(cutoff)/potpar[3] + 2*cutoff/SQR(potpar[3])
	   		   + 2.0 / CUBE(potpar[3])) * exp(-potpar[3]*cutoff));
      else
         return( 0.0 );
   }
}
/******************************************************************************
 *  kernel   Innermost loop of force calculation.  Takes a vector of squared  *
 *  atomic distances (r_sqr), charges (chg) , pot'l parameters (pot[which]),  *
 *  and returns a vector of forces (forceij) and the potential energy (pe)    *
 *  Norm is 2*pi/sqrt(alpha), ptype is potential type selector.		      *
 ******************************************************************************/
void	kernel(j0, nnab, forceij, pe, r_sqr, nab_chg, d_chg, d_norm,
	       d_alpha, ptype, pot)
int	ptype, j0, nnab;
double	*pe;
double	d_norm, d_alpha, d_chg;
real	forceij[], r_sqr[], nab_chg[], **pot;
{
   real
   		ppe = 0.0,
   		expa2r2, t,		/* Temporary for component of force   */
   		erfc_term,exp_f1, exp_f2,/* Temporary for b*exp(-cr) etc      */
   		r, r_r,	  		/* Magnitude of site-site vector      */
   		r_sqr_r, r_6_r, r_12_r;	/* Powers of above                    */
   real		norm = d_norm, alpha = d_alpha, chg = d_chg;
   register int	jsite;
   real one = 1.0, two = 2.0, half = 0.5, six = 6.0;

   switch(ptype)
   {
    default:
      message(NULLI, NULLP, FATAL, UNKPTY, ptype);
    case LJPOT:
      if(alpha > 0.0)
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = one / r;
	    r_sqr_r = r_r * r_r;
	    r_6_r = SQR(pot[1][jsite])* r_sqr_r;
	    r_6_r   = r_6_r * r_6_r * r_6_r;
	    r_12_r  = SQR(r_6_r);
	    expa2r2 = nab_chg[jsite]* chg * exp(-SQR(alpha*r));
	    t = one/(one+half*(alpha*r));
	    erfc_term = t * expa2r2 * exp(POLY(t)) * r_r;
	    ppe += erfc_term + pot[0][jsite]*(r_12_r - r_6_r);
	    forceij[jsite] = r_sqr_r*(six*pot[0][jsite]*(two*r_12_r - r_6_r)
				      + erfc_term + norm * expa2r2);
	 }
      else
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = one / r;
	    r_sqr_r = r_r * r_r;
	    r_6_r = SQR(pot[1][jsite])* r_sqr_r;
	    r_6_r   = r_6_r * r_6_r * r_6_r;
	    r_12_r  = SQR(r_6_r);
	    ppe += pot[0][jsite]*(r_12_r - r_6_r);
	    forceij[jsite] = r_sqr_r*six*pot[0][jsite]*(two*r_12_r - r_6_r);
	 }
      break;
    case E6POT:
      if(alpha >= 0.0)
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = one / r;
	    r_sqr_r = r_r * r_r;
	    r_6_r   = pot[0][jsite] * r_sqr_r * r_sqr_r * r_sqr_r;
	    expa2r2 = nab_chg[jsite]* chg * exp(-SQR(alpha*r));
	    t = one/(one+half*(alpha*r));
	    erfc_term = t * expa2r2 * exp(POLY(t)) * r_r;
	    exp_f1 = pot[1][jsite] * exp(-pot[2][jsite] * r);
	    ppe += erfc_term - r_6_r + exp_f1;
	    forceij[jsite] = r_sqr_r*(-six* r_6_r+ erfc_term + norm * expa2r2)
	    + pot[2][jsite]*exp_f1 * r_r;
	 }
      else
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = one / r;
	    r_sqr_r = r_r * r_r;
	    r_6_r   = pot[0][jsite] * r_sqr_r * r_sqr_r * r_sqr_r;
	    exp_f1 = pot[1][jsite] * exp(-pot[2][jsite] * r);
	    ppe +=  - r_6_r + exp_f1;
	    forceij[jsite] = -r_sqr_r * six * r_6_r
	                     + pot[2][jsite]*exp_f1 * r_r;
	 }
      
      break;      
    case MCYPOT:
      if(alpha >= 0.0)
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = one / r;
	    r_sqr_r = r_r * r_r;
	    expa2r2 = nab_chg[jsite]* chg * exp(-SQR(alpha*r));
	    t = one/(one+half*(alpha*r));
	    erfc_term = t * expa2r2 * exp(POLY(t)) * r_r;
	    exp_f1 =  pot[0][jsite] * exp(-pot[1][jsite]*r);
	    exp_f2 = -pot[2][jsite] * exp(-pot[3][jsite]*r);
	    ppe += erfc_term + exp_f1 + exp_f2;
	    forceij[jsite] = (pot[1][jsite]*exp_f1 + pot[3][jsite]*exp_f2) * r_r
	    + (erfc_term + norm * expa2r2) * r_sqr_r;
	 }
       else
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = one / r;
	    r_sqr_r = r_r * r_r;
	    exp_f1 =  pot[0][jsite] * exp(-pot[1][jsite]*r);
	    exp_f2 = -pot[2][jsite] * exp(-pot[3][jsite]*r);
	    ppe += exp_f1 + exp_f2;
	    forceij[jsite] = (pot[1][jsite]*exp_f1 + pot[3][jsite]*exp_f2) *r_r;
	 }
      break;      
   }
   *pe += ppe;
}
