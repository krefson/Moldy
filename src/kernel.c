/******************************************************************************
 * kernel	Functions to calculate the forces, potential and distant pot'l*
 *		correction for various potentials.  ALL POTENTIAL-DEPENDANT   *
 *		CODE is in this module.					      *
 * dist_pot()		Return distant-potential correction		      *
 * kernel()		Calculate pe and forces				      *
 * types[]		Array of names of potential function types.	      *
 * npotp[]		Array containing number of parameters for each type.  *
 * npott		size of above two arrays.			      *
 *		Since CRAY CC 4.0 won't vectorise library function calls, a   *
 *		FORTRAN equivalent to kernel() is provided in kernel.f. To    *
 *		use it, compile force.c with macro FKERNEL defined and link   *
 *		with it AND THIS MODULE.				      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	kernel.c,v $
 * Revision 1.6  89/12/15  12:56:53  keith
 * Added conditional ionclusion of <fastmath.h> for stellar
 * 
 * Revision 1.5  89/11/20  13:29:26  keith
 * Replaced separate arrays "types" and "npotp" with array of structs "potspec"
 * 
 * Revision 1.4  89/10/12  16:42:37  keith
 * Eleminated a few inefficiencies in non-coulombic LJ and MCY loops.
 * 
 * Revision 1.3  89/07/04  18:42:07  keith
 * Fixed error in kernel and force which led to sites being allocated the
 * wrong potential parameters.  Needed extra parameter to kernel.
 * 
 * Revision 1.2  89/06/20  18:22:31  keith
 * Moved pot. par defs arrays 'types', 'npotp' and npott from input.c to kernel.c
 * 
 * Revision 1.1  89/04/20  16:00:45  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore1/keith/md/moldy/RCS/kernel.c,v 1.6 89/12/15 12:56:53 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#if  defined(convexvc) || defined(stellar)
#include <fastmath.h>
#else
#include <math.h>
#endif
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
void message();
/*========================== Potential type specification ====================*/
pots_t	potspec[]  = {{"lennard-jones",2},
		      {"buckingham",3},
                      {"mcy",4}};
int	npott=(sizeof potspec / sizeof(pots_t));
/*========================== Macros ==========================================*/

#define E1	 0.254829592		/* Polynomial Constants used in	      */
#define E2	-0.284496736		/* Evaluation of the complementary    */
#define E3	 1.421413741		/* Error function.  		      */
#define E4	-1.453152027		/* Approximation used is that of      */
#define E5	 1.061405429		/* Abramowitz & Stegun p299.	      */

#define PP	 0.3275911

#define POLY5(t)   ((t)*(E1 + (t)*(E2 + (t)*(E3 + (t)*(E4 + (t)*E5)))))

#define LJPOT 0
#define E6POT 1
#define MCYPOT 2
/*============================================================================*/
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
void	kernel(j0, nnab, forceij, pe, r_sqr, nab_chg, chg, norm,
	       alpha, ptype, pot)
int	ptype, j0, nnab;
double	*pe;
double	norm, alpha, chg;
real	forceij[], r_sqr[], nab_chg[];
register real	*pot[];
{
   register real t, r, r_r, expa2r2, erfc_term, r_6_r;
   real		ppe = 0.0,
   		exp_f1, exp_f2,		/* Temporary for b*exp(-cr) etc      */
   		r_sqr_r, r_12_r;	/* Powers of above                    */
   register int	jsite;
   real *p0 = pot[0], *p1 = pot[1], *p2 = pot[2], *p3 = pot[3];

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
	    r_r	 = 1.0 / r;
	    r_sqr_r = SQR(r_r);
	    r_6_r = SQR(p1[jsite])* r_sqr_r;
	    r_6_r   = CUBE(r_6_r);
	    r_12_r  = SQR(r_6_r);
	    expa2r2 = nab_chg[jsite]* chg * exp(-SQR(alpha*r));
	    t = 1.0/(1.0+PP*(alpha*r));
	    erfc_term = POLY5(t) * expa2r2 * r_r;
	    ppe += erfc_term + p0[jsite]*(r_12_r - r_6_r);
	    forceij[jsite] = r_sqr_r*(6.0*p0[jsite]*(2*r_12_r - r_6_r)
				      + erfc_term + norm * expa2r2);
	 }
      else
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r_sqr_r = 1.0 / r_sqr[jsite];
	    r_6_r = SQR(p1[jsite])* r_sqr_r;
	    r_6_r   = CUBE(r_6_r);
	    r_12_r  = SQR(r_6_r);
	    ppe += p0[jsite]*(r_12_r - r_6_r);
	    forceij[jsite] = r_sqr_r*6.0*p0[jsite]*(2*r_12_r - r_6_r);
	 }
      break;
    case E6POT:
      if(alpha >= 0.0)
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = 1.0 / r;
	    r_sqr_r = SQR(r_r);
	    r_6_r   = p0[jsite] * CUBE(r_sqr_r);
	    expa2r2 = nab_chg[jsite]* chg * exp(-SQR(alpha*r));
	    t = 1.0/(1.0+PP*(alpha*r));
	    erfc_term = POLY5(t) * expa2r2 * r_r;
	    exp_f1 = p1[jsite] * exp(-p2[jsite] * r);
	    ppe += erfc_term - r_6_r + exp_f1;
	    forceij[jsite] = r_sqr_r*(-6.0* r_6_r+ erfc_term + norm * expa2r2)
	    + p2[jsite]*exp_f1 * r_r;
	 }
      else
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = 1.0 / r;
	    r_sqr_r = SQR(r_r);
	    r_6_r   = p0[jsite] * r_sqr_r * r_sqr_r * r_sqr_r;
	    exp_f1 = p1[jsite] * exp(-p2[jsite] * r);
	    ppe +=  - r_6_r + exp_f1;
	    forceij[jsite] = -r_sqr_r * 6.0 * r_6_r
	                     + p2[jsite]*exp_f1 * r_r;
	 }
      
      break;      
    case MCYPOT:
      if(alpha >= 0.0)
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = 1.0 / r;
	    r_sqr_r = SQR(r_r);
	    expa2r2 = nab_chg[jsite]* chg * exp(-SQR(alpha*r));
	    t = 1.0/(1.0+PP*(alpha*r));
	    erfc_term = POLY5(t) * expa2r2 * r_r;
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    exp_f2 = -p2[jsite] * exp(-p3[jsite]*r);
	    ppe += erfc_term + exp_f1 + exp_f2;
	    forceij[jsite] = (p1[jsite]*exp_f1 + p3[jsite]*exp_f2) * r_r
	    + (erfc_term + norm * expa2r2) * r_sqr_r;
	 }
       else
VECTORIZE
         for(jsite=j0; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = 1.0 / r;
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    exp_f2 = -p2[jsite] * exp(-p3[jsite]*r);
	    ppe += exp_f1 + exp_f2;
	    forceij[jsite] = (p1[jsite]*exp_f1 + p3[jsite]*exp_f2) *r_r;
	 }
      break;      
   }
   *pe += ppe;
}
