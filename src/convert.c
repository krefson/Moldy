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
 * Convert	Functions for conversion of units of input parameters.	      *
 ******************************************************************************
 *      Revision Log
 *       $Log: convert.c,v $
 *       Revision 2.13  2000/12/06 17:45:28  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.12  2000/11/09 16:25:48  keith
 *       Fixed const argument consistency warning
 *
 *       Revision 2.11  2000/04/27 17:57:06  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.10  1999/11/01 17:22:10  keith
 *       Got rid of useless (and lint-unfriendly) _mp declarations.
 *       Made "const" qualifier conditionally excluded if using lint.
 *
 *       Revision 2.9  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 *
 *       Revision 2.8  1995/12/07 17:47:59  keith
 *       Reworked V. Murashov's thermostat code.
 *       Convert mass params from kJ/mol ps^2 to prog units. Defaults=1.
 *
 *       Revision 2.8  1995/12/06 17:49:03  keith
 *       Updated conversion table for thermostat parameters.
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
 * Revision 2.5  94/01/18  13:32:14  keith
 * Null update for XDR portability release
 * 
 * Revision 2.3  93/10/28  10:27:46  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:49:00  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.5  93/03/09  15:58:24  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.4  92/10/28  14:10:01  keith
 * Changed "site_[tp]" typedefs to avoid name clash on HP.
 * 
 * Revision 1.3  91/08/15  18:11:49  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.2  90/05/21  15:29:28  keith
 * Moved definition of struct pot_dim[][] from convert.c to kernel.c.
 * 
 * Revision 1.1  89/04/20  16:00:36  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/kr/CVS/moldy/src/convert.c,v 2.13 2000/12/06 17:45:28 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<math.h>
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External data references ========================*/
extern	      contr_mt	control;            /* Main simulation control parms. */
extern  const dim_mt	pot_dim[][NPOTP];    /* Pot'l dimension specification */
/*========================== External function declarations ==================*/
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
/*========================== Structs local to module =========================*/
/* This struct array contains all details for conversions in 'control' struct */
typedef struct
{
  double	*d;			/* Pointer to item to be converted    */
  dim_mt	dim;			/* Dimensions of quantity (eg MLT(-2))*/
  unit_mt	unit;			/* Units (in MKS) for conversion from */
}  conv_mt;
/*========================== Global variables ================================*/
#define TMUNIT (MUNIT/CONV_TM)
static const conv_mt 
	conv[] = {  {&control.step,	{0,0,1,0},	{1,1,TUNIT,1}},
                    {&control.pressure,	{1,-1,-2,0},	{1e6,1,1,1}},
		    {&control.pmass,	{1,0,0,0},	{MUNIT,1,1,1}},
		    {&control.ttmass,	{1,2,0,0},	{TMUNIT,LUNIT,1,1}},
		    {&control.rtmass,	{1,2,0,0},	{TMUNIT,LUNIT,1,1}},
		    {&control.cutoff,	{0,1,0,0},	{1,LUNIT,1,1}},
		    {&control.density,	{1,-3,0,0},	{.001,.01,1,1}},
		    {&control.alpha,	{0,-1,0,0},	{1,LUNIT,1,1}},
		    {&control.limit,	{0,-1,0,0},	{1,LUNIT,1,1}}};
#define	NCONV	(sizeof conv / sizeof(conv_mt))
/*============================================================================*/
/******************************************************************************
 *  unit_scale  Takes the dimensions of a quantity (dim), the units it is in  *
 *  (unit_from) and returns the scale by which it must be multiplied to       *
 *  convert it to new units (unit_to).  Uses logarithms to avoid overflow.    *
 ******************************************************************************/
#define	MAX_SCALE	80			/* 1e35 - safe for any machine*/
static double	unit_scale(const dim_mt *dim,   /* Dimensions    	      */
			   const unit_mt *unit_from, const unit_mt *unit_to)
{
   double	lnscale = 	dim->m*(log(unit_from->m) - log(unit_to->m))
			      + dim->l*(log(unit_from->l) - log(unit_to->l))
			      + dim->t*(log(unit_from->t) - log(unit_to->t))
			      + dim->q*(log(unit_from->q) - log(unit_to->q));
   if(lnscale < -MAX_SCALE || lnscale > MAX_SCALE)
      message(NULLI,NULLP,FATAL,BADUNI,lnscale);
   return(exp(lnscale));
}
/******************************************************************************
 *  Convert_potentials   Scale potential parameters, site masses and charges  *
 *  from input units to program units.					      *
 ******************************************************************************/
void	
conv_potentials(const unit_mt *unit_from,
		const unit_mt *unit_to, /* Values of units for conversion     */ 
		pot_mt *potpar, 	/* Array of potpar records[max_id**2] */ 
		int npotpar, 		/* Number of 'active' parameters      */ 
		int ptype, 		/* Potential type                     */ 
		site_mt *site_info,	/* Site specification array[max_id]   */ 
		int max_id)		/* How many site id's                 */
{
   int		idi, idj, ip;		/* Counters for id's and potpar	      */
   static dim_mt	mass_dim   = {1,0,0,0},	/* Dimensions of mass	      */
                charge_dim = {0,0,0,1};	/* Dimensions of charge		      */
   double	mscale = unit_scale(&mass_dim,   unit_from, unit_to),
                qscale = unit_scale(&charge_dim, unit_from, unit_to);
					/* Scale factors for mass and charge  */
   double	potscale[NPOTP];	/* Scale factors for pot'l parameters */

   for(ip = 0; ip < npotpar;ip++)	/* Work out scale factors for pot'l   */
      potscale[ip] = unit_scale(&pot_dim[ptype][ip],unit_from, unit_to);

   for(idi = 0; idi < max_id; idi++)
   {
      site_info[idi].mass   *= mscale;	/* Scale mass of each site	      */
      site_info[idi].charge *= qscale;	/* Scale charge of each site	      */
      for(idj = 0; idj < max_id; idj++)
          for(ip = 0; ip < npotpar; ip++)	/* Scale potential parameters */
             potpar[idi + max_id*idj].p[ip] *= potscale[ip];
   }
}
/******************************************************************************
 *   convert_control  Convert various quantities in the control record between*
 *   input and program units.  Which to convert are indicated by the array    *
 *   'conv' which contains a pointer to the data and a dimension struct.      *
 ******************************************************************************/
void	conv_control(const unit_mt *unit,/* Units conversion is to/from       */ 
		     boolean direction)	 /* True=from input, false=to input   */
{
   int		ic;			/* Counter			      */
   for(ic = 0; ic < NCONV; ic++)
      if(direction)
         *conv[ic].d *= unit_scale(&conv[ic].dim, &conv[ic].unit, unit);
      else
         *conv[ic].d *= unit_scale(&conv[ic].dim, unit, &conv[ic].unit);
}
