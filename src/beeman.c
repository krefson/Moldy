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
 *   Beeman - Routines to implement the modified Beeman algorithm for         *
 *            stepping the co-ordinates in a MD simulation including rigid    *
 *            molecules                                                       *
 *   Original algorithm - D. Beeman J. Comp. Phys 20, 131-139 (1976)          *
 *   Modified with velocity predictor to handle velocity-dependant forces     *
 *                      - K. Refson Physica 131B 256-266 (1985)               *
 *                                                                            *
 *   External routines:          step_1(system), step_2(system), beeman_2()   *
 *   External references:        none                                         *
 *   External data:              none                                         *
 *   External data references:   control                                      *
 *                                                                            *
 *   Note:  The routines beeman_1, beeman_2 and predict treat their arguments *
 *          as one dimensional arrays for simplicity.  Thus the 'number'      *
 *          argument is 3 (for c_of_m) or 4 (for quaternions) times the number*
 *          of particles.  This may cause 'consistency' messages from lint.   *
 ******************************************************************************
 *      Revision Log
 *       $Log: beeman.c,v $
 * Revision 2.7  1994/06/08  13:09:29  keith
 * Protected against possible bus error for systems with no
 * rotational freedom by making all references to "quat" etc conditional
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Got rid of all global (external) data items except for
 * "control" struct and constant data objects.  The latter
 * (pot_dim, potspec, prog_unit) are declared with CONST
 * qualifier macro which evaluates to "const" or nil
 * depending on ANSI/K&R environment.
 * Also moved as many "write" instantiations of "control"
 * members as possible to "startup", "main" leaving just
 * "dump".
 *
 * Declared as "static"  all functions which should be.
 *
 * Revision 2.5  1994/01/18  13:32:13  keith
 * Null update for XDR portability release
 *
 * Revision 2.3  93/10/28  10:27:45  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:48:58  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.5  93/03/15  14:41:37  keith
 * Added GPL copyleft notice to permit release and distribution.
 * N.B.  Previous versions were copyright (C) by author and 
 * only licensed by explicit permission.
 * 
 * Revision 1.4  93/03/09  15:58:21  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.3  91/08/15  18:11:47  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.2  89/10/24  17:18:37  keith
 * Modified pbc algorithm to use floor() library function.
 * Now works with non-orthorhombic cell.
 * 
 * Revision 1.1  89/04/20  16:00:35  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore_data/keith/md/moldy/RCS/beeman.c,v 2.7 1994/06/08 13:09:29 keith stab $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include <math.h>
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
#if defined(ANSI) || defined(__STDC__)
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
#else
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
#endif
/*========================== External data references ========================*/
extern	contr_mt	control;            /* Main simulation control parms. */
/*========================== Macros ==========================================*/
#define	TOLERANCE	1.0e-4
/*============================================================================*/
/******************************************************************************
 *   Normalise the new quaternions                                            *
 ******************************************************************************/
static void normalise(quat,n)
int		n;			/* Number of quaternions  (in)        */
quat_mp		quat;			/* Quaternions       (update)         */
{
   register int		i, j;
   register double	norm;

   for(i=0; i<n; i++)
   {
      norm = 0.0;
      for(j=0; j<4; j++)  norm += quat[i][j]*quat[i][j];
      norm = sqrt(norm);
      if( fabs(norm - 1.0) > TOLERANCE)
         message(NULLI, NULLP, FATAL, QNORM2, i,
		 quat[i][0],quat[i][1],quat[i][2],quat[i][3]);
      for(j=0; j<4; j++)  quat[i][j] /= norm;
   }
}
/******************************************************************************
 *   Constrain the quaternion derivatives				      *
 ******************************************************************************/
void constrain(quat, qdot ,n)
int		n;			/* Number of quaternions  (in)        */
quat_mp		quat;			/* Quaternions            (in)        */
quat_mp		qdot;			/* Quaternion derivatives (update)    */
{
   register int		i, j;
   register double	delta;

   for(i=0; i<n; i++)
   {
      delta = 0;
      for(j = 0; j < 4; j++)
         delta += quat[i][j]*qdot[i][j];
      if( fabs(delta)*control.step > TOLERANCE)
         message(NULLI, NULLP, FATAL, QCONST, i, delta);
      for(j=0; j<4; j++)  qdot[i][j] -= delta * quat[i][j];
   }
}
/******************************************************************************
 *   Apply periodic boundary conditions to put particles back in MD box       *
 ******************************************************************************/
static
void escape(c_of_m, nmols)
int		nmols;		/* First dimension of c-of-m                  */
vec_mp		c_of_m;		/* Centre of mass co-ordinates (updat)        */

{
   int	imol;			/* Molecule counter			      */

   for(imol = 0; imol < nmols; imol++)
   {
      c_of_m[imol][0] -= floor(c_of_m[imol][0] + 0.5);
      c_of_m[imol][1] -= floor(c_of_m[imol][1] + 0.5);
      c_of_m[imol][2] -= floor(c_of_m[imol][2] + 0.5);
   }
}
/******************************************************************************
 *   Step co-ordinates as first  (P(r)) stage of Beeman algorithm             *
 ******************************************************************************/
static void beeman_1(x, xdot, xddot, xddoto, n) 
int		n;			/* Length of vectors       (in)       */
real		x[],			/* Co-ordinate array   (update)       */
		xdot[],			/* First derivatives       (in)       */
		xddot[],		/* Accelerations           (in)       */
		xddoto[];		/* Old accelerations       (in)       */
{
     register int	i;
     register real	step = control.step,
			step_sq = step * step / 6.0;
            
      for(i=0; i<n; i++)
         x[i] += step*xdot[i] + step_sq*(4.0*xddot[i] - xddoto[i]);
}

/******************************************************************************
 *   Predict velocities at next timestep using 3-rd order predictor           *
 *   (3-rd order not yet implemented - for now use second order)              *
 ******************************************************************************/
/*ARGSUSED*/
static void predict(v, vp, a, ao, avo, n)
int		n;			/* Length of vectors       (in)       */
real		v[],			/* Velocities              (in)       */
		vp[],			/* Predicted velocities   (out)       */
		a[],			/* Accelerations           (in)       */
		ao[],			/* Old accelerations       (in)       */
		avo[];			/* Very old accelerations  (in)       */
{
   register int		i;
   register real	step = control.step/2.0;
   /*			step6 = step/3.0 */
            
   for(i=0; i<n; i++)
      vp[i] = v[i] + step * (3.0*a[i] - ao[i]);
/*     vp[i] = v[i] + step6 *(10.0*a[i] - 5.0*ao[i] + avo[i]); */
}

/******************************************************************************
 *   Step velocities as second (c(v)) stage of Beeman algorithm               *
 ******************************************************************************/
void beeman_2(v_in, v_out, a, ao, avo, n)
int		n;			/* Length of vectors       (in)       */
real		v_in[],			/* Velocities              (in)       */
                v_out[],		/* Velocities		  (out)	      */
		a[],			/* Accelerations           (in)       */
		ao[],			/* Old accelerations       (in)       */
		avo[];			/* Very old accelerations  (in)       */
{
   register int		i;
   register real	step = control.step/6.0;
            
   for(i=0; i<n; i++)
      v_out[i] = v_in[i] + step * (2.0*a[i] + 5.0*ao[i] - avo[i]);
}

/******************************************************************************
 *   Apply first stage of modified beeman algorithm to the whole system       *
 *   Step the centre of mass co-ordinates, quaternions and Parinello & Rahman *
 *   h matrix (if const pressure simulation), and predict the corresponding   *
 *   velocities.                                                              *
 ******************************************************************************/
void step_1(sys)
system_mp	sys;			/* pointer to whole-system record     */
{
   beeman_1(sys->c_of_m[0],sys->vel[0],sys->acc[0],sys->acco[0], 3*sys->nmols);
   escape(sys->c_of_m, sys->nmols);
   if( sys->nmols_r > 0 )
   {
      beeman_1(sys->quat[0],sys->qdot[0],sys->qddot[0],sys->qddoto[0],
	       4*sys->nmols_r);
      normalise(sys->quat, sys->nmols_r);
   
   }
   if(control.const_pressure)
      beeman_1(sys->h[0], sys->hdot[0], sys->hddot[0], sys->hddoto[0], 9);

   predict(sys->vel[0], sys->velp[0], sys->acc[0], 
           sys->acco[0], sys->accvo[0], 3*sys->nmols);
   if( sys->nmols_r > 0 )
   {
      predict(sys->qdot[0], sys->qdotp[0], sys->qddot[0], 
	      sys->qddoto[0], sys->qddotvo[0], 4*sys->nmols_r);
      constrain(sys->quat, sys->qdotp, sys->nmols_r);
   
      if (control.const_temp)
         predict(sys->ra, sys->rap, sys->radot, sys->radoto,
                 sys->radotvo, sys->nspecies);
   }
   if(control.const_pressure)
      predict(sys->h[0], sys->hdotp[0], sys->hdot[0],
              sys->hddoto[0], sys->hddotvo[0], 9);

   if (control.const_temp)
      predict(sys->ta, sys->tap, sys->tadot, sys->tadoto,
              sys->tadotvo, sys->nspecies);
#ifdef DEBUG3
   printf("Step1\n");
#endif 
}
   
/******************************************************************************
 *   Apply second stage of modified beeman algorithm to the whole system      *
 *   Step the centre of mass velocities, and the derivatives of the           *
 *   quaternions and Parinello & Rahman h matrix (if const pressure           *
 *   simulation), and predict the corresponding velocities.                   *
 ******************************************************************************/
void step_2(sys)
system_mp	sys;			/* pointer to whole-system record     */
{
   beeman_2(sys->vel[0], sys->vel[0], sys->acc[0], sys->acco[0],sys->accvo[0],
	    3*sys->nmols);
   if( sys->nmols_r > 0 )
   {
      beeman_2(sys->qdot[0], sys->qdot[0], sys->qddot[0], sys->qddoto[0],
	    sys->qddotvo[0], 4*sys->nmols_r);
      constrain(sys->quat, sys->qdot, sys->nmols_r);

      if (control.const_temp)
         beeman_2(sys->ra, sys->ra, sys->radot, sys->radoto,
                  sys->radotvo, sys->nspecies);

   }
   if(control.const_pressure)
      beeman_2(sys->hdot[0], sys->hdot[0], sys->hddot[0], sys->hddoto[0],
	       sys->hddotvo[0], 9);

   if (control.const_temp)
      beeman_2(sys->ta, sys->ta, sys->tadot, sys->tadoto,
               sys->tadotvo, sys->nspecies);
#ifdef DEBUG3
   printf("Step2\n");
#endif 
}
