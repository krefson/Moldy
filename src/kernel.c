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
 * Revision 1.15  92/09/22  14:45:20  keith
 * A few efficiency improvements for the Titan
 * 
 * Revision 1.14  92/06/10  15:53:16  keith
 * Added new potential type "generic" for Neal.
 * 
 * Revision 1.13  91/08/16  15:25:35  keith
 * Checked error returns from fread, fwrite, fseek and fclose more
 * rigourously.   Called strerror() to report errors.
 * 
 * Revision 1.12  91/08/15  18:12:03  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.11  91/02/19  14:51:26  keith
 * Minor changes to get rid of misleading compiler warnings.
 * 
 * Revision 1.10  90/09/28  13:29:39  keith
 * Inserted braces around VECTORIZE directives and changed include files
 * for STARDtardent 3000 series (via cond. comp symbol "ardent").
 * 
 * Revision 1.9  90/05/21  15:29:00  keith
 * Moved definition of struct pot_dim[][] from convert.c to kernel.c.
 * 
 * Revision 1.8  90/05/02  13:04:52  keith
 * Tidied up and restructured code to improve readability.
 * Reduced number of scalar temporaries in vector loops - this is
 * necessary to vectorize under CRAY CC 5.0.
 * 
 * Revision 1.7  89/12/22  19:33:44  keith
 * Added p0-p3 - temporaries for pot[0..3] pointers within loops.
 * Some machines (eg stellar) generate better code this way.
 * 
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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/kernel.c,v 1.15 92/09/22 14:45:20 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#if  defined(convexvc) || defined(stellar)
#   include <fastmath.h>
#else
#ifdef ardent
#   include <vmath.h>
#else
#   include <math.h>
#endif
#endif
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
void message();
/*========================== Potential type specification ====================*/
pots_mt	potspec[]  = {{"lennard-jones",2},
		      {"buckingham",3},
                      {"mcy",4},
		      {"generic",6}};
int	npott=(sizeof potspec / sizeof(pots_mt));
/*
 *  Array of dimensions of pot'l parameters.  Powers of {m,l,t} per parameter.
 */
dim_mt   pot_dim[][NPOTP]= {{{1,2,-2},{0,1,0}},
                           {{1,8,-2},{1,2,-2},{0,-1,0}},
		           {{1,2,-2},{0,-1,0},{1,2,-2},{0,-1,0}},
			   {{1,2,-2},{0,-1,0},{1,14,-2},
			    {1,6,-2},{1,8,-2},{1,10,-2}}};

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
#define GENPOT 3
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
      if( potpar[3] != 0.0 )
         return( potpar[2] * (SQR(cutoff)/potpar[3] + 2*cutoff/SQR(potpar[3])
	   		   + 2.0 / CUBE(potpar[3])) * exp(-potpar[3]*cutoff));
      else
         return( 0.0 );
    case GENPOT:
      return 0.0;
   }
}
/******************************************************************************
 *  kernel   Innermost loop of force calculation.  Takes a vector of squared  *
 *  atomic distances (r_sqr), charges (chg) , pot'l parameters (pot[which]),  *
 *  and returns a vector of forces (forceij) and the potential energy (pe)    *
 ******************************************************************************/
void	kernel(jmin, nnab, forceij, pe, r_sqr, nab_chg, chg, norm,
	       alpha, ptype, pot)
int	jmin, nnab;		/* Lower and upper limits for vectors.   (in) */
int	ptype;			/* Index of potential type in potspec[]. (in) */
double	*pe;			/* Potential energy accumulator.     (in/out) */
double	alpha, norm;		/* Ewald parameter and 2*alpha/sqrt(pi). (in) */
double	chg;			/* Electric charge of reference site.    (in) */
real	forceij[];		/* Vector of dU(r)/dr for each in r_sqr.(out) */
real	r_sqr[];		/* Vector of site-site distances (**2).  (in) */
real	nab_chg[];		/* Vector of charges of neighbour sites. (in) */
real	*pot[];			/* Vectors of potential parameters.	 (in) */
{
   register real t, ar;			/* Argument of erfc() polynomial.     */
   register real r;			/* Site-site distance.		      */
   register real r_r, r_6_r, r_sqr_r, r_12_r,	/* Reciprocal powers of r.    */
                 r_4_r, r_8_r;
   register real erfc_term;		/* Intermediates in erfc calculation. */
   	    real ppe = 0.0;		/* Local accumulator of pot. energy.  */
   	    real exp_f1, exp_f2;	/* Temporary for b*exp(-cr) etc       */
   register int	jsite;			/* Loop counter for vectors.	      */
   real *p0 = pot[0], *p1 = pot[1],     /* Local bases for arrays of pot'l    */
        *p2 = pot[2], *p3 = pot[3],     /* parameters.			      */
        *p4 = pot[4], *p5 = pot[5];

   if(alpha > 0.0)
      switch(ptype)
      {
       default:
	 message(NULLI, NULLP, FATAL, UNKPTY, ptype);
       case LJPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    r       = sqrt(r_sqr[jsite]);
	    ar	    = alpha*r;
	    r_r	 = 1.0 / r;
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    t = 1.0/(1.0+PP*ar);
	    t = POLY5(t) * erfc_term * r_r;
	    erfc_term = t + norm * erfc_term;
	    r_sqr_r = SQR(r_r);
	    /*
	     * Non-coulombic ie potential-specific part
	     */
	    r_6_r = SQR(p1[jsite])* r_sqr_r;
	    r_6_r   = CUBE(r_6_r);
	    r_12_r  = SQR(r_6_r);
	    ppe += t + p0[jsite]*(r_12_r - r_6_r);
	    forceij[jsite] = r_sqr_r*(6.0*p0[jsite]*(2*r_12_r - r_6_r)
				      + erfc_term);
	 }
	 break;
       case E6POT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    r       = sqrt(r_sqr[jsite]);
	    ar	    = alpha*r;
	    r_r	 = 1.0 / r;
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    t = 1.0/(1.0+PP*ar);
	    t = POLY5(t) * erfc_term * r_r;
	    erfc_term = t + norm * erfc_term;
	    r_sqr_r = SQR(r_r);
	    /*
	     * Non-coulombic ie potential-specific part
	     */
	    exp_f1 = p1[jsite] * exp(-p2[jsite] * r);
	    r_6_r   = p0[jsite] * CUBE(r_sqr_r);
	    ppe += t - r_6_r + exp_f1;
	    forceij[jsite] = r_sqr_r*(-6.0* r_6_r+ erfc_term)
	                             + p2[jsite]*exp_f1 * r_r;
	 }
	 break;
       case MCYPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    r       = sqrt(r_sqr[jsite]);
	    ar	    = alpha*r;
	    r_r	 = 1.0 / r;
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    t = 1.0/(1.0+PP*ar);
	    t = POLY5(t) * erfc_term * r_r;
	    erfc_term = t + norm * erfc_term;
	    r_sqr_r = SQR(r_r);
            /*
	     * Non-coulombic ie potential-specific part
	     */
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    exp_f2 = -p2[jsite] * exp(-p3[jsite]*r);
	    ppe +=t + exp_f1 + exp_f2;
	    forceij[jsite] = (p1[jsite]*exp_f1 + p3[jsite]*exp_f2) * r_r
	    + erfc_term * r_sqr_r;
	 }
	 break;      
       case GENPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    r       = sqrt(r_sqr[jsite]);
	    ar	    = alpha*r;
	    r_r	 = 1.0 / r;
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    t = 1.0/(1.0+PP*ar);
	    t = POLY5(t) * erfc_term * r_r;
	    erfc_term = t + norm * erfc_term;
	    r_sqr_r = SQR(r_r);
            /*
	     * Non-coulombic ie potential-specific part
	     */
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    r_4_r = SQR(r_sqr_r);
	    r_6_r = r_sqr_r * r_4_r;
	    r_8_r = p5[jsite] * SQR(r_4_r);
	    r_12_r = p2[jsite] * SQR(r_6_r);
	    r_4_r *= p3[jsite];
	    r_6_r *= p4[jsite];

	    ppe += t + exp_f1 + r_12_r -r_4_r - r_6_r - r_8_r;
	    forceij[jsite] = r_sqr_r*( 12.0*r_12_r - 4.0*r_4_r - 6.0*r_6_r 
				      - 8.0*r_8_r + erfc_term)
	                   + p1[jsite]*exp_f1 * r_r;
	 }
	 break;      
      }
   else
      switch(ptype)
      {
       default:
	 message(NULLI, NULLP, FATAL, UNKPTY, ptype);
       case LJPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
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
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
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
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    exp_f2 = -p2[jsite] * exp(-p3[jsite]*r);
	    ppe += exp_f1 + exp_f2;
	    r_r	 = 1.0 / r;
	    forceij[jsite] = (p1[jsite]*exp_f1 + p3[jsite]*exp_f2) *r_r;
	 }
	 break; 
       case GENPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = 1.0 / r;
	    r_sqr_r = SQR(r_r);
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    r_4_r = SQR(r_sqr_r);
	    r_6_r = r_sqr_r * r_4_r;
	    r_8_r = p5[jsite] * SQR(r_4_r);
	    r_12_r = p2[jsite] * SQR(r_6_r);
	    r_4_r *= p3[jsite];
	    r_6_r *= p4[jsite];

	    ppe += exp_f1 + r_12_r -r_4_r - r_6_r - r_8_r;
	    forceij[jsite] = r_sqr_r*( 12.0*r_12_r - 4.0*r_4_r - 6.0*r_6_r 
				      - 8.0*r_8_r)
	                   + p1[jsite]*exp_f1 * r_r;
	 }
	 break;      
      }     
   *pe += ppe;
}
