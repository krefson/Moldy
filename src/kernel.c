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
 *       $Log: kernel.c,v $
 *       Revision 2.14.2.4.2.2  2002/03/15 15:17:02  kr
 *       Fixed some bugs in forces
 *
 *       Revision 2.14.2.4.2.1  2002/03/13 10:27:52  kr
 *       Trial version incorporating reciprocal-space summation for r^-2 and r^-6
 *       interactions.  This version implements a new potential "genpot46" to activate.
 *
 *       Revision 2.14.2.4  2002/02/18 16:05:12  kr
 *       Fixed dist-pot term to include all terms for Generic and Buckingham potentials.
 *       This makes a difference for pathalogical cases like Floris ion-water potentials.
 *       Renamed "hiw7+win" to "hiwfl+win" for compatibility with Rafael's version.
 *
 *       Revision 2.14.2.4  2002/02/18 15:29:36  kr
 *       Fixed dist-pot term to include all terms for Generic and Buckingham potentials.
 *       This makes a difference for pathalogical cases like Floris ion-water potentials.
 *
 *       Revision 2.14.2.3  2001/05/28 17:48:36  keith
 *       Added "hiw7+win" potential for flexible Cr-6H2O cluster potential
 *
 *       Revision 2.14.2.2  2001/02/22 11:51:47  keith
 *       Added "generic+win" potential.
 *
 *       Revision 2.14.2.1  2000/12/08 17:11:22  keith
 *       Added Window potential back onto HIW
 *
 *       Revision 2.14  2000/12/08 15:22:07  keith
 *       Reorganized order of potentials to maintain existing numberinf of HIW
 *       etc -- for restart file compatibility.
 *
 *       Revision 2.13  2000/12/08 12:22:33  keith
 *       Incorporated Morse and HIW potentials
 *
 *       Revision 2.12  2000/12/07 10:11:08  keith
 *       Tidied up prototypes and incorporated mods of 2.10.4.3.
 *
 *       Revision 2.10.4.3  2000/12/06 17:45:30  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.10.4.2  2000/08/31 11:26:34  keith
 *       Fixed bugs in non-coulombic code which cause FPE errors on systems which trap
 *       Revision 2.11  2000/04/27 17:57:08  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.10  1999/09/20 10:27:30  keith
 *       Updated comments to assist with adding a new potential.
 *
 *       Revision 2.9  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 * Revision 1997/10/24 Wunderlich/Fisher 
 * New Morse Type potential with 7 para,eters (B-H) included
 * This potential is compatibel with the Mxdorto progarm written by Kawabata Tokyo
 * Remark: Colomb Term A/r is already includued in Moldy. 
 * Phi(r_ij) =         B * exp ((C-R)*D)   - E/r^6 +
 *		    +  F * exp (-2*G(r-H)) - 2F*exp (-G(r-H))
 *
 *       Revision 2.8  1996/10/04 17:27:24  keith
 *       Rescheduled line order to overlap divide/computation on DEC Alpha/T3D.
 *
 *       Revision 2.7  1994/06/08 13:22:31  keith
 *       Null update for version compatibility
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
 * Revision 2.5  94/01/18  13:32:38  keith
 * Null update for XDR portability release
 * 
 * Revision 2.4  94/01/13  12:46:53  keith
 * Aedded distant porential correction for GENPOT potential (NTS 13/1/94.
 * 
 * 
 * Revision 2.3  93/10/28  10:27:57  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:49:11  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.17  93/03/12  12:24:46  keith
 * Got rid of unneccesary convex special case.
 * 
 * Revision 1.16  93/03/09  15:58:41  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
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
static char *RCSid = "$Header: /home/kr/CVS/moldy/src/kernel.c,v 2.14.2.4.2.2 2002/03/15 15:17:02 kr Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include <math.h>
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
void	message(int *,...);		/* Write a warning or error message   */
/*========================== Potential type specification ====================*/
#define LJPOT 0		/* Lennard Jones.   U=p0*((p1/r)^6 - (p1/r)^12)       */
#define E6POT 1		/* 6-exp potential  U=-p0/r^6 + p1*exp(-p2*r)         */
#define MCYPOT 2	/* MCY water pot.  J.Chem.Phys 64,1351(1976)          */
			/*		    U= p0*exp(-p1*r) - p2*exp(-p3*r)  */
#define GENPOT 3	/* "generic" potential for multipurpose use.          */
			/* U= p0*exp(-p1*r) + p2/r^12 - p3/r^4 -p4/r^6 -p5/r^8*/
#define HIWPOT 4	/* HIW pot. R.R.Pappalardo, J.Phys.Chem 97,4500(1993) */
			/*		    U= p0/r^4 + p1/r^6 + p2/r^12      */
#define HIWWIN  5	/* HIW flexible + windows potential                   */
#define MORPOT 6	/* Busing-Ida-Gilbert plus Morse-potential of Mdxorto */
			/* Material Design using Personal computer, Ed        */
                        /* Kazuyuki Hirao,  (1994) ISBN 4-7853-6803-9         */
			/*     U=    B * exp ((C-R)*D)   - E/r^6 +            */
			/*           F * exp (-2*G(r-H)) - 2F*exp (-G(r-H))   */
#define GENWIN 7	/* "generic" potential for multipurpose use. +        */
                        /* harmonic window potential                          */
#define HIW7WIN 8       /* Version of HIW potential including r**-7 term      */
#define GENPOT46 9	/* "generic" potential with Ewald summation           */
			/* U= p0*exp(-p1*r) + p2/r^12 - p3/r^4 -p4/r^6 -p5/r^8*/

const pots_mt	potspec[]  = {{"lennard-jones",2},  /* Name, index & # parms  */
		              {"buckingham",3},
                              {"mcy",4},
		              {"generic",6},
			      {"hiw",3},
		              {"hiw+win",5},
		              {"morse",7},
		              {"generic+win",8},
			      {"hiwfl+win",6},
		              {"generic46",6},
		              {0,0}};	            /* MUST be null-terminated*/
/*
 *  Array of dimensions of pot'l parameters.  Triplets contain powers
 *  of {m,l,t} for each parameter and are used to convert from input
 *  to program units.   E.g LJ, e is an energy: kgm**2s-2 => {1,2,-2},
 *  sigma ls a length => {0,1,0}.
 */
const dim_mt   pot_dim[][NPOTP]= {
   /* Lennard-Jones */	{{1,2,-2},{0,1,0}},
   /* Buckingham    */  {{1,8,-2},{1,2,-2},{0,-1,0}},
   /* MCY           */  {{1,2,-2},{0,-1,0},{1,2,-2},{0,-1,0}},
   /* Generic       */  {{1,2,-2},{0,-1,0},{1,14,-2},{1,6,-2},{1,8,-2},{1,10,-2}}, 
   /* HIW           */  {{1,6,-2},{1,8,-2},{1,14,-2}},
   /* HIW+Win       */  {{1,6,-2},{1,8,-2},{1,14,-2},{1,0,-2},{0,1,0}},
   /* Morse         */  {{1,2,-2},{0, 1,0},{0,-1, 0},
			          {1,8,-2},{1,2,-2},{0,-1,0},{0, 1, 0}},
   /* Generic+Win   */  {{1,2,-2},{0,-1,0},{1,14,-2},{1,6,-2},
                         {1,8,-2},{1,10,-2},{1,0,-2},{0,1,0}}, 
   /* HIW7+Win      */  {{1,6,-2},{1,8,-2},{1,9,-2},{1,14,-2},{1,0,-2},{0,1,0}},
   /* Generic-4-6   */  {{1,2,-2},{0,-1,0},{1,14,-2},{1,6,-2},{1,8,-2},{1,10,-2}}, 
                                  };

/*========================== Macros ==========================================*/

#define E1	 0.254829592		/* Polynomial Constants used in	      */
#define E2	-0.284496736		/* Evaluation of the complementary    */
#define E3	 1.421413741		/* Error function.  		      */
#define E4	-1.453152027		/* Approximation used is that of      */
#define E5	 1.061405429		/* Abramowitz & Stegun p299.	      */

#define PP	 0.3275911

#define POLY5(t)   ((t)*(E1 + (t)*(E2 + (t)*(E3 + (t)*(E4 + (t)*E5)))))
#define BPAR_TOL 1.0e-7
/*============================================================================*/
/******************************************************************************
 *  mkpot46(). Build 2d potential perameter arrays for Ewald -4-6 sum.        *
 *  The code depends on the potential type, hence inclusion here.             *
 ******************************************************************************/
void mkpot46(real **pot4, real **pot6, int max_id, int ptype, pot_mt *potpar)
{
   int id, jd;
   switch(ptype)
   {
    default:
      message(NULLI, NULLP, FATAL, UNKPTY, ptype);
      /*FALLTHRU*/
   case GENPOT46:
      for( id = 1; id < max_id; id++)
	 for( jd = 1; jd < max_id; jd++)
	 {
	    pot4[id][jd] = -potpar[id*max_id+jd].p[3];
	    pot6[id][jd] = -potpar[id*max_id+jd].p[4];
	 }
   }
   
}
/******************************************************************************
 *  dist_pot   return attractive part of potential integrated outside cutoff. *
 *  dist_pot = - int_{r_c}^{infty} r^2 U(r) dr                                *
 ******************************************************************************/
double	dist_pot(real *potpar,          /* Array of potential parameters      */
		 double cutoff,         /* Cutoff distance                    */ 
		 int ptype)             /* Potential type selector            */
{
   switch(ptype)
   {
    default:
      message(NULLI, NULLP, FATAL, UNKPTY, ptype);
      /*FALLTHRU*/
    case LJPOT:
      return(potpar[0]*CUBE(SQR(potpar[1])/cutoff) / 3.0);
    case E6POT:
       if( potpar[2] > BPAR_TOL ) 
	  return( potpar[0] / ( 3.0*CUBE(cutoff))
		  - potpar[1] * exp(-potpar[2]*cutoff)
		  * (SQR(cutoff)/potpar[2] + 2*cutoff/SQR(potpar[2]) + 2.0 / CUBE(potpar[2])));
       else
	  return( potpar[0] / ( 3.0*CUBE(cutoff))
		  - potpar[1] * exp(-potpar[2]*cutoff)
		  * (SQR(cutoff)/potpar[2] + 2*cutoff/SQR(potpar[2]) + 2.0 / CUBE(potpar[2])));
    case MCYPOT:
      if( potpar[3]  > BPAR_TOL )
         return( potpar[2] * (SQR(cutoff)/potpar[3] + 2*cutoff/SQR(potpar[3])
	   		   + 2.0 / CUBE(potpar[3])) * exp(-potpar[3]*cutoff));
      else
         return( 0.0 );
    case GENWIN:
       /*FALLTHRU*/
    case GENPOT:
       if( potpar[1] > BPAR_TOL ) 
	  return ( - potpar[0] * exp(-potpar[1]*cutoff) *
		   (SQR(cutoff)/potpar[1] + 2*cutoff/SQR(potpar[1]) + 2.0 / CUBE(potpar[1])) 
		   -potpar[2] / ( 9.0*CUBE(CUBE(cutoff))) + potpar[3] / cutoff 
		   + potpar[4] / ( 3.0*CUBE(cutoff)) + potpar[5] / ( 5.0*SQR(cutoff)*CUBE(cutoff)));
       else
	  return ( -potpar[2] / ( 9.0*CUBE(CUBE(cutoff))) + potpar[3] / cutoff 
		   + potpar[4] / ( 3.0*CUBE(cutoff)) + potpar[5] / ( 5.0*SQR(cutoff)*CUBE(cutoff)));
    case GENPOT46:
       if( potpar[1] > BPAR_TOL ) 
	  return ( - potpar[0] * exp(-potpar[1]*cutoff) *
		   (SQR(cutoff)/potpar[1] + 2*cutoff/SQR(potpar[1]) + 2.0 / CUBE(potpar[1])) 
		   -potpar[2] / ( 9.0*CUBE(CUBE(cutoff))) 
		    + potpar[5] / ( 5.0*SQR(cutoff)*CUBE(cutoff)));
       else
	  return ( -potpar[2] / ( 9.0*CUBE(CUBE(cutoff))) 
		    + potpar[5] / ( 5.0*SQR(cutoff)*CUBE(cutoff)));
	  
    case MORPOT:
      if( potpar[5] != 0.0 )
         return( potpar[3] / ( 3.0*CUBE(cutoff))
            +2.0*potpar[4] * (SQR(cutoff)/potpar[5] + 2*cutoff/SQR(potpar[5])
	    +2.0 / CUBE(potpar[5])) * exp(-potpar[5]*(cutoff-potpar[6])));
      else
         return( potpar[3] / ( 3.0*CUBE(cutoff)) );
    case HIWWIN:
       /*FALLTHRU*/
    case HIWPOT:
         return( - potpar[0] /cutoff - potpar[1] /CUBE(cutoff)/3.0 -
                   potpar[2] /CUBE(CUBE(cutoff))/9.0);
    case HIW7WIN:
         return( - potpar[0] /cutoff - potpar[1] /CUBE(cutoff)/3.0
                 - potpar[2] /(4.0*cutoff*CUBE(cutoff))
                 - potpar[3] /CUBE(CUBE(cutoff))/9.0);
   }
}
/******************************************************************************
 *  kernel   Innermost loop of force calculation.  Takes a vector of squared  *
 *  atomic distances (r_sqr), charges (chg) , pot'l parameters (pot[which]),  *
 *  and returns a vector of forces (forceij) and the potential energy (pe)    *
 ******************************************************************************/
void	kernel(int jmin,     
	       int nnab,        /* Lower and upper limits for vectors.   (in) */
	       real *forceij,   /* Vector of dU(r)/dr for each in r_sqr.(out) */
	       double *pe,      /* Potential energy accumulator.     (in/out) */
	       real *r_sqr,     /* Vector of site-site distances (**2).  (in) */
	       real *nab_chg,   /* Vector of charges of neighbour sites. (in) */
	       double chg,      /* Electric charge of reference site.    (in) */
	       double alpha,    /* Ewald parameter and 2*alpha/sqrt(pi). (in) */
	       double alpha46,  /* Ewald parameter and 2*alpha/sqrt(pi). (in) */
	       int ptype,       /* Index of potential type in potspec[]. (in) */ 
	       real **pot)      /* Vectors of potential parameters.      (in) */
{
   register real t, ar;			/* Argument of erfc() polynomial.     */
   register real r;			/* Site-site distance.		      */
   register real r_r, r_6_r, r_sqr_r, r_12_r,	/* Reciprocal powers of r.    */
                 r_4_r, r_7_r, r_8_r, rsq;
   register real erfc_term, arfac, arfacsq;/* Intermediates in erfc calculation. */
   	    real ppe = 0.0;		/* Local accumulator of pot. energy.  */
   	    real exp_f1, exp_f2, exp_f3; /* Temporary for b*exp(-cr) etc      */
   register real rmr0, fwin;            /* Temporaries for window potential   */
   register int	jsite;			/* Loop counter for vectors.	      */
	    double norm = 2.0*alpha/ROOTPI;   /* Coulombic prefactor*/
	    double alpha46sq = SQR(alpha46);
   real *p0 = pot[0], *p1 = pot[1],     /* Local bases for arrays of pot'l    */
        *p2 = pot[2], *p3 = pot[3],     /* parameters.			      */
        *p4 = pot[4], *p5 = pot[5],
        *p6 = pot[6], *p7 = pot[7];

   if(alpha > 0.0)
      switch(ptype)
      {
       default:
	 message(NULLI, NULLP, FATAL, UNKPTY, ptype);
	 /*FALLTHRU*/
       case LJPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    r       = sqrt(r_sqr[jsite]);
	    ar	    = alpha*r;
	    t = 1.0/(1.0+PP*ar);
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    r_r	 = 1.0 / r;
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
	    t = 1.0/(1.0+PP*ar);
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    r_r	 = 1.0 / r;
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
	    t = 1.0/(1.0+PP*ar);
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    r_r	 = 1.0 / r;
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
	    t = 1.0/(1.0+PP*ar);
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    r_r	 = 1.0 / r;
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
       case GENPOT46:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    rsq = r_sqr[jsite];
	    r       = sqrt(rsq);
	    ar	    = alpha*r;
	    t = 1.0/(1.0+PP*ar);
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    r_r	 = 1.0 / r;
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
	    /*
	     * Ewald-like summation of r^-4 and r^-6
	     */
	    arfac = alpha46sq*rsq;
	    arfacsq = SQR(arfac);
	    exp_f2 = exp(-arfac);

	    ppe += t + exp_f1 + r_12_r - ((1.0+arfac)*r_4_r + (1.0+arfac+0.5*arfacsq)*r_6_r)*exp_f2 - r_8_r;
	    forceij[jsite] = r_sqr_r*( 12.0*r_12_r 
			   - ((4.0+4.0*arfac+2.0*arfacsq)*r_4_r+(6.0+6.0*arfac+3.0*arfacsq+arfacsq*arfac)*r_6_r)*exp_f2
				      - 8.0*r_8_r + erfc_term)
	                   + p1[jsite]*exp_f1 * r_r;
	 }
	 break;      
       case GENWIN:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    r       = sqrt(r_sqr[jsite]);
	    ar	    = alpha*r;
	    t = 1.0/(1.0+PP*ar);
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    r_r	 = 1.0 / r;
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
	    rmr0 = (r - p7[jsite]);
	    fwin = p6[jsite]*rmr0;

	    ppe += t + exp_f1 + r_12_r -r_4_r - r_6_r - r_8_r + 0.5*fwin*rmr0;
	    forceij[jsite] = r_sqr_r*( 12.0*r_12_r - 4.0*r_4_r - 6.0*r_6_r 
				      - 8.0*r_8_r + erfc_term - r*fwin)
	                   + p1[jsite]*exp_f1 * r_r;
	 }
	 break;      
       case MORPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    r         = sqrt(r_sqr[jsite]);
	    ar	      = alpha*r;
	    t         = 1.0/(1.0+PP*ar);
	    erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
	    r_r	      = 1.0 / r;
	    t         = POLY5(t) * erfc_term * r_r;
	    erfc_term = t + norm * erfc_term;
	    r_sqr_r   = SQR(r_r);
            /*
	     * Non-coulombic ie potential-specific part
	     */
	    exp_f1 = p0[jsite] * exp(( p1[jsite] - r)*p2[jsite]);
	    r_6_r  = p3[jsite] * CUBE(r_sqr_r);
	    exp_f2 = p4[jsite] * exp(-2.0 * p5[jsite] * ( r - p6[jsite]));
	    exp_f3 =-p4[jsite] * 2.0 * exp(-p5[jsite] * ( r - p6[jsite]));
	    ppe += t + exp_f1 - r_6_r + exp_f2  + exp_f3;
	    forceij[jsite] = r_sqr_r*( -6.0*r_6_r + erfc_term) 
	                   + r_r    *(p2[jsite] *exp_f1    
	                   +(2.0    * p5[jsite])*exp_f2 
	                   +          p5[jsite] *exp_f3 );
	 }
	 break;      
       case HIWPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
         {
            /*
             * Calculate r and coulombic part
             */
            r       = sqrt(r_sqr[jsite]);
            ar      = alpha*r;
            t = 1.0/(1.0+PP*ar);
            erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
            r_r  = 1.0 / r;
            t = POLY5(t) * erfc_term * r_r;
            erfc_term = t + norm * erfc_term;
            r_sqr_r = SQR(r_r);
            /*
             * Non-coulombic ie potential-specific part
             */
            r_4_r = SQR(r_sqr_r);
            r_6_r = r_sqr_r * r_4_r;
            r_12_r = p2[jsite] * SQR(r_6_r);
            r_6_r *= p1[jsite];
            r_4_r *= p0[jsite];

            ppe += t + r_4_r + r_6_r + r_12_r;

            forceij[jsite] =   r_sqr_r * ( 4.0 * r_4_r + 6.0 * r_6_r  
                               + 12.0 * r_12_r + erfc_term );  
         }
         break; 
       case HIWWIN:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
         {
            /*
             * Calculate r and coulombic part
             */
            r       = sqrt(r_sqr[jsite]);
            ar      = alpha*r;
            t = 1.0/(1.0+PP*ar);
            erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
            r_r  = 1.0 / r;
            t = POLY5(t) * erfc_term * r_r;
            erfc_term = t + norm * erfc_term;
            r_sqr_r = SQR(r_r);
            /*
             * Non-coulombic ie potential-specific part
             */
            r_4_r = SQR(r_sqr_r);
            r_6_r = r_sqr_r * r_4_r;
            r_12_r = p2[jsite] * SQR(r_6_r);
            r_6_r *= p1[jsite];
            r_4_r *= p0[jsite];
	    rmr0 = (r - p4[jsite]);
	    fwin = p3[jsite]*rmr0;

            ppe += t + r_4_r + r_6_r + r_12_r + 0.5*fwin*rmr0;

            forceij[jsite] =   r_sqr_r * ( 4.0 * r_4_r + 6.0 * r_6_r  
                               + 12.0 * r_12_r + erfc_term - r*fwin);  
         }
         break; 
        case HIW7WIN:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
         {
            /*
             * Calculate r and coulombic part
             */
            r       = sqrt(r_sqr[jsite]);
            ar      = alpha*r;
            t = 1.0/(1.0+PP*ar);
            erfc_term = nab_chg[jsite]* chg * exp(-SQR(ar));
            r_r  = 1.0 / r;
            t = POLY5(t) * erfc_term * r_r;
            erfc_term = t + norm * erfc_term;
            r_sqr_r = SQR(r_r);
            /*
             * Non-coulombic ie potential-specific part
             */
            r_4_r = SQR(r_sqr_r);
            r_6_r = r_sqr_r * r_4_r;
	    r_7_r = r_6_r*r_r;
            r_12_r = p3[jsite] * SQR(r_6_r);
            r_7_r *= p2[jsite];
            r_6_r *= p1[jsite];
            r_4_r *= p0[jsite];
	    rmr0 = (r - p5[jsite]);
	    fwin = p4[jsite]*rmr0;

            ppe += t + r_4_r + r_6_r + r_7_r + r_12_r + 0.5*fwin*rmr0;

            forceij[jsite] =   r_sqr_r * ( 4.0 * r_4_r + 6.0 * r_6_r + 7.0 * r_7_r 
                               + 12.0 * r_12_r + erfc_term - r*fwin);  
         }
         break; 
      }
   else
      switch(ptype)
      {
       default:
	 message(NULLI, NULLP, FATAL, UNKPTY, ptype);
	 /*FALLTHRU*/
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
	    exp_f1 = p1[jsite] * exp(-p2[jsite] * r);
	    r_sqr_r = SQR(r_r);
	    r_6_r   = p0[jsite] * r_sqr_r * r_sqr_r * r_sqr_r;
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
	    r_r	 = 1.0 / r;
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    exp_f2 = -p2[jsite] * exp(-p3[jsite]*r);
	    ppe += exp_f1 + exp_f2;
	    forceij[jsite] = (p1[jsite]*exp_f1 + p3[jsite]*exp_f2) *r_r;
	 }
	 break; 
       case GENPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = 1.0 / r;
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    r_sqr_r = SQR(r_r);
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
       case GENPOT46:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    rsq = r_sqr[jsite];
	    r       = sqrt(rsq);
	    r_r	 = 1.0 / r;
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
	    /*
	     * Ewald-like summation of r^-4 and r^-6
	     */
	    arfac = alpha46sq*rsq;
	    arfacsq = SQR(arfac);
	    exp_f2 = exp(-arfac);

	    ppe += exp_f1 + r_12_r - ((1.0+arfac)*r_4_r + (1.0+arfac+0.5*arfacsq)*r_6_r)*exp_f2 - r_8_r;
	    forceij[jsite] = r_sqr_r*( 12.0*r_12_r 
			   - ((4.0+4.0*arfac+2.0*arfacsq)*r_4_r+(6.0+6.0*arfac+3.0*arfacsq+arfacsq*arfac)*r_6_r)*exp_f2
				      - 8.0*r_8_r)
	                   + p1[jsite]*exp_f1 * r_r;
	 }
	 break;      
       case GENWIN:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    r       = sqrt(r_sqr[jsite]);
	    r_r	 = 1.0 / r;
	    exp_f1 =  p0[jsite] * exp(-p1[jsite]*r);
	    r_sqr_r = SQR(r_r);
	    r_4_r = SQR(r_sqr_r);
	    r_6_r = r_sqr_r * r_4_r;
	    r_8_r = p5[jsite] * SQR(r_4_r);
	    r_12_r = p2[jsite] * SQR(r_6_r);
	    r_4_r *= p3[jsite];
	    r_6_r *= p4[jsite];
	    rmr0 = (r - p7[jsite]);
	    fwin = p6[jsite]*rmr0;

	    ppe += exp_f1 + r_12_r -r_4_r - r_6_r - r_8_r + 0.5*fwin*rmr0;
	    forceij[jsite] = r_sqr_r*( 12.0*r_12_r - 4.0*r_4_r - 6.0*r_6_r 
				      - 8.0*r_8_r - r*fwin)
	                   + p1[jsite]*exp_f1 * r_r;

	 }
	 break;      
       case MORPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
	 {
	    /*
	     * Calculate r and coulombic part
	     */
	    r       = sqrt(r_sqr[jsite]);
	    r_r	    = 1.0 / r;
	    r_sqr_r = SQR(r_r);
	    exp_f1  = p0[jsite] * exp(( p1[jsite] - r)*p2[jsite]);
	    r_6_r   = p3[jsite] * CUBE(r_sqr_r);
	    exp_f2  = p4[jsite] * exp(-2.0 * p5[jsite] * ( r - p6[jsite]));
	    exp_f3  =-p4[jsite] * 2.0 * exp(-p5[jsite] * ( r - p6[jsite]));
	    ppe    +=     exp_f1 - r_6_r + exp_f2  + exp_f3;
	    forceij[jsite] = r_sqr_r*( -6.0*r_6_r            ) 
	                   + r_r    *(p2[jsite] *exp_f1    
	                   +(2.0    * p5[jsite])*exp_f2 
	                   +          p5[jsite] *exp_f3 );
	 }
	 break;      
       case HIWPOT:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
         {
            r       = sqrt(r_sqr[jsite]);
            r_r  = 1.0 / r;
            r_sqr_r = SQR(r_r);
            r_4_r = SQR(r_sqr_r);
            r_6_r = r_sqr_r * r_4_r;
            r_12_r = SQR(r_6_r) * p2[jsite];
            r_6_r *= p1[jsite]; 
            r_4_r *= p0[jsite];

            ppe += r_4_r + r_6_r + r_12_r;    

            forceij[jsite] = r_sqr_r * ( 4.0 * r_4_r + 6.0 * r_6_r  
                               + 12.0 * r_12_r);  
	 }  
         break;      
       case HIWWIN:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
         {
            r       = sqrt(r_sqr[jsite]);
            r_r  = 1.0 / r;
            r_sqr_r = SQR(r_r);
            r_4_r = SQR(r_sqr_r);
            r_6_r = r_sqr_r * r_4_r;
            r_12_r = SQR(r_6_r) * p2[jsite];
            r_6_r *= p1[jsite]; 
            r_4_r *= p0[jsite];
	    rmr0 = (r - p4[jsite]);
	    fwin = p3[jsite]*rmr0;

            ppe += r_4_r + r_6_r + r_12_r + 0.5*fwin*rmr0;

            forceij[jsite] =   r_sqr_r * ( 4.0 * r_4_r + 6.0 * r_6_r  
                               + 12.0 * r_12_r  - r*fwin); 
         }
         break; 
       case HIW7WIN:
VECTORIZE
         for(jsite=jmin; jsite < nnab; jsite++)
         {
            r       = sqrt(r_sqr[jsite]);
            r_r  = 1.0 / r;
            r_sqr_r = SQR(r_r);
            r_4_r = SQR(r_sqr_r);
            r_6_r = r_sqr_r * r_4_r;
	    r_7_r = r_6_r*r_r;
            r_12_r = SQR(r_6_r) * p3[jsite];
            r_7_r *= p2[jsite];
            r_6_r *= p1[jsite]; 
            r_4_r *= p0[jsite];
	    rmr0 = (r - p5[jsite]);
	    fwin = p4[jsite]*rmr0;

            ppe += r_4_r + r_6_r + r_7_r + r_12_r + 0.5*fwin*rmr0;

            forceij[jsite] =   r_sqr_r * ( 4.0 * r_4_r + 6.0 * r_6_r + 7.0 * r_7_r  
                               + 12.0 * r_12_r  - r*fwin); 
         }
         break; 
      }     
   *pe += ppe;
}
