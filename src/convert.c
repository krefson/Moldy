/******************************************************************************
 * Convert	Functions for conversion of units of input parameters.	      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	convert.c,v $
 * Revision 1.2  90/05/21  15:29:28  keith
 * Moved definition of struct pot_dim[][] from convert.c to kernel.c.
 * 
 * Revision 1.1  89/04/20  16:00:36  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/convert.c,v 1.3 91/08/14 14:23:25 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<math.h>
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External data references ========================*/
extern	contr_t	control;
/*========================== External function declarations ==================*/
void		message();
/*========================== Structs local to module =========================*/
/* This struct array contains all details for conversions in 'control' struct */
typedef struct
{
  double	*d;			/* Pointer to item to be converted    */
  dim_t		dim;			/* Dimensions of quantity (eg MLT(-2))*/
  unit_t	unit;			/* Units (in MKS) for conversion from */
}  conv_t,  *conv_p;
/*========================== Global variables ================================*/
extern dim_t    pot_dim[][NPOTP];
#define	NPOTT	(sizeof npotp / sizeof(int))
		/* Specification array for conversion between units 	      */

static	conv_t	conv[] = {  {&control.step,	{0,0,1,0},	{1,1,TUNIT,1}},
                            {&control.pressure,	{1,-1,-2,0},	{1e6,1,1,1}},
			    {&control.pmass,	{1,0,0,0},	{MUNIT,1,1,1}},
			    {&control.cutoff,	{0,1,0,0},	{1,LUNIT,1,1}},
			    {&control.density,	{1,-3,0,0},	{.001,.01,1,1}},
			    {&control.alpha,	{0,-1,0,0},	{1,LUNIT,1,1}},
			    {&control.limit,	{0,-1,0,0},	{1,LUNIT,1,1}}};
#define	NCONV	(sizeof conv / sizeof(conv_t))
/*============================================================================*/
/******************************************************************************
 *  unit_scale  Takes the dimensions of a quantity (dim), the units it is in  *
 *  (unit_from) and returns the scale by which it must be multiplied to       *
 *  convert it to new units (unit_to).  Uses logarithms to avoid overflow.    *
 ******************************************************************************/
#define	MAX_SCALE	80			/* 1e35 - safe for any machine*/
static double	unit_scale(dim, unit_from, unit_to)
dim_p	dim;				/* Dimensions			      */
unit_p	unit_from, unit_to;		/* Units to convert from/to	      */
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
void	conv_potentials(unit_from, unit_to, potpar, npotpar, ptype,
			   site_info, max_id)
unit_p		unit_from, unit_to;	/* Values of units for conversion     */
pot_t		potpar[];		/* Array of potpar records[max_id**2] */
int		npotpar;		/* Number of 'active' parameters      */
int		ptype;			/* Potential type		      */
site_t		site_info[];		/* Site specification array[max_id]   */
int		max_id;			/* How many site id's		      */
{
   int		idi, idj, ip;		/* Counters for id's and potpar	      */
   static dim_t	mass_dim   = {1,0,0,0},	/* Dimensions of mass		      */
                charge_dim = {0,0,0,1};	/* Dimensions of charge		      */
   double	mscale = unit_scale(&mass_dim,   unit_from, unit_to),
                qscale = unit_scale(&charge_dim, unit_from, unit_to);
					/* Scale factors for mass and charge  */
   double	potscale[NPOTP];	/* Scale factors for pot'l parameters */

   for(ip = 0; ip < npotpar;ip++)	/* Work out scale factors for pot'l   */
      potscale[ip] = unit_scale(pot_dim[ptype]+ip,unit_from, unit_to);

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
void	conv_control(unit, direction)
unit_p		unit;			/* Units conversion is to/from	      */
boolean		direction;		/* True=from input, false=to input    */
{
   int		ic;			/* Counter			      */
   for(ic = 0; ic < NCONV; ic++)
      if(direction)
         *conv[ic].d *= unit_scale(&conv[ic].dim, &conv[ic].unit, unit);
      else
         *conv[ic].d *= unit_scale(&conv[ic].dim, unit, &conv[ic].unit);
}
