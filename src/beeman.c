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
 *       $Log:	beeman.c,v $
 * Revision 1.1  89/04/20  16:00:35  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/tigger/keith/md/RCS/beeman.c,v 1.1 89/04/20 16:00:35 keith Stab $";
#endif
/*========================== Library include files ===========================*/
#include <stdio.h>
#include <math.h>
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
void	message();			/* Error  and exit handler	      */
/*========================== External data references ========================*/
extern	contr_t	control;
/*========================== Macros ==========================================*/
#define	TOLERANCE	1.0e-4
/*============================================================================*/
/******************************************************************************
 *   Normalise the new quaternions                                            *
 ******************************************************************************/
static void normalise(quat,n)
int		n;			/* Number of quaternions  (in)        */
quat_p		quat;			/* Quaternions       (update)         */
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
static void constrain(quat, qdot ,n)
int		n;			/* Number of quaternions  (in)        */
quat_p		quat;			/* Quaternions            (in)        */
quat_p		qdot;			/* Quaternion derivatives (update)    */
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
void escape(c_of_m, nmols)
int		nmols;		/* First dimension of c-of-m                  */
vec_p		c_of_m;		/* Centre of mass co-ordinates (updat)        */

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
system_p	sys;			/* pointer to whole-system record     */
{
   beeman_1(sys->c_of_m[0],sys->vel[0],sys->acc[0],sys->acco[0], 3*sys->nmols);
   escape(sys->c_of_m, sys->nmols);
   beeman_1(sys->quat[0],sys->qdot[0],sys->qddot[0],sys->qddoto[0],
                                                               4*sys->nmols_r);
   normalise(sys->quat, sys->nmols_r);
   if(control.const_pressure)
      beeman_1(sys->h[0], sys->hdot[0], sys->hddot[0], sys->hddoto[0], 9);
   
   predict(sys->vel[0], sys->velp[0], sys->acc[0], 
           sys->acco[0], sys->accvo[0], 3*sys->nmols);
   predict(sys->qdot[0], sys->qdotp[0], sys->qddot[0], 
           sys->qddoto[0], sys->qddotvo[0], 4*sys->nmols_r);
   constrain(sys->quat, sys->qdotp, sys->nmols_r);
   if(control.const_pressure)
      predict(sys->h[0], sys->hdotp[0], sys->hdot[0],
              sys->hddoto[0], sys->hddotvo[0], 9);
}
   
/******************************************************************************
 *   Apply second stage of modified beeman algorithm to the whole system      *
 *   Step the centre of mass velocities, and the derivatives of the           *
 *   quaternions and Parinello & Rahman h matrix (if const pressure           *
 *   simulation), and predict the corresponding velocities.                   *
 ******************************************************************************/
void step_2(sys)
system_p	sys;			/* pointer to whole-system record     */
{
   beeman_2(sys->vel[0], sys->vel[0], sys->acc[0], sys->acco[0],sys->accvo[0],
	    3*sys->nmols);
   beeman_2(sys->qdot[0], sys->qdot[0], sys->qddot[0], sys->qddoto[0],
	    sys->qddotvo[0], 4*sys->nmols_r);
   constrain(sys->quat, sys->qdot, sys->nmols_r);
   if(control.const_pressure)
      beeman_2(sys->hdot[0], sys->hdot[0], sys->hddot[0], sys->hddoto[0],
	       sys->hddotvo[0], 9);
}
