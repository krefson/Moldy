/******************************************************************************
 * quaterns	Functions for manipulating quaternions			      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	quaterns.c,v $
 * Revision 1.1  89/04/20  16:00:52  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/quaterns.c,v 1.1 89/04/20 16:00:52 keith Stab $";
#endif
/*========================== Library include files ===========================*/
#include <math.h>
/*========================== Program include files ===========================*/
#include "defs.h"
/*========================== External function declarations ==================*/
double precision();
/*============================================================================*/
/******************************************************************************
 * Quaternion multiplier.  Multiplies arrays of quaternions p by q to give r  *
 * Can be called with  r the same as p or q.                                  *
 ******************************************************************************/
void q_mul(p, q, r, n)
quat_mp	p,			/* First Quaternion array [n][4]        (in)  */
	q,			/* Second quaternion array [n][4]       (in)  */
	r;			/* Resultant quaternions [n][4]        (out)  */
int	n;			/* Number of quaternions in the arrays  (in)  */
{
   register	int	i;
   register	real	p0, p1, p2, p3;
   register	real	q0, q1, q2, q3;
   
   for(i = 0; i < n; i++)
   {
      p0 =  p[i][0]; p1 =  p[i][1]; p2 =  p[i][2]; p3 =  p[i][3];
      q0 =  q[i][0]; q1 =  q[i][1]; q2 =  q[i][2]; q3 =  q[i][3];

      r[i][0] = p0*q0 - p1*q1 - p2*q2 - p3*q3;
      r[i][1] = p1*q0 + p0*q1 - p3*q2 + p2*q3;
      r[i][2] = p2*q0 + p3*q1 + p0*q2 - p1*q3;
      r[i][3] = p3*q0 - p2*q1 + p1*q2 + p0*q3;
   }
}
/******************************************************************************
 * Quaternion multiplier.  Multiplies quaternions p by q to give r  	      *
 * Can be called with  r the same as p or q.                                  *
 ******************************************************************************/
void q_mul_1(p, q, r)
quat_mt	p,			/* First Quaternion array [n][4]        (in)  */
	q,			/* Second quaternion array [n][4]       (in)  */
	r;			/* Resultant quaternions [n][4]        (out)  */
{
   register	real	p0, p1, p2, p3;
   register	real	q0, q1, q2, q3;
   
   p0 =  p[0]; p1 =  p[1]; p2 =  p[2]; p3 =  p[3];
   q0 =  q[0]; q1 =  q[1]; q2 =  q[2]; q3 =  q[3];

   r[0] = p0*q0 - p1*q1 - p2*q2 - p3*q3;
   r[1] = p1*q0 + p0*q1 - p3*q2 + p2*q3;
   r[2] = p2*q0 + p3*q1 + p0*q2 - p1*q3;
   r[3] = p3*q0 - p2*q1 + p1*q2 + p0*q3;
}
/******************************************************************************
 * Quaternion multiplier.  Multiplies arrays of quaternions p(-1) by q to     *
 * give r. Can be called with  r the same as p or q.                          *
 ******************************************************************************/
void q_conj_mul(p, q, r, n)
quat_mp	p,			/* First Quaternion array [n][4]        (in)  */
	q,			/* Second quaternion array [n][4]       (in)  */
	r;			/* Resultant quaternions [n][4]        (out)  */
int	n;			/* Number of quaternions in the arrays  (in)  */
{
   register	int	i;
   register	real	p0, p1, p2, p3;
   register	real	q0, q1, q2, q3;
   
   for(i = 0; i < n; i++)
   {
      p0 =  p[i][0]; p1 = -p[i][1]; p2 = -p[i][2]; p3 = -p[i][3];
      q0 =  q[i][0]; q1 =  q[i][1]; q2 =  q[i][2]; q3 =  q[i][3];

      r[i][0] = p0*q0 - p1*q1 - p2*q2 - p3*q3;
      r[i][1] = p1*q0 + p0*q1 - p3*q2 + p2*q3;
      r[i][2] = p2*q0 + p3*q1 + p0*q2 - p1*q3;
      r[i][3] = p3*q0 - p2*q1 + p1*q2 + p0*q3;
   }
}
/******************************************************************************
 *  q_to_mat  Make the rotation matrix corresponding to quaternion q          *
 ******************************************************************************/
void q_to_rot(quat, rot)
quat_mt	quat;				/* Input quaternion              (in) */
mat_mt	rot;				/* Rotation matrix		(out) */
{
   register	real	q0, q1, q2, q3;
   register	real	a01, a02, a03,
			     a12, a13,
			          a23;

   q0 = quat[0];  q1 = quat[1]; q2 = quat[2];  q3 = quat[3];

   a01 = 2.0*q0*q1;  a02 = 2.0*q0*q2;  a03 = 2.0*q0*q3;
                     a12 = 2.0*q1*q2;  a13 = 2.0*q1*q3;
                                       a23 = 2.0*q2*q3;
   rot[0][1] = a12 - a03;
   rot[0][2] = a13 + a02;
   rot[1][0] = a12 + a03;
   rot[1][2] = a23 - a01;
   rot[2][0] = a13 - a02;
   rot[2][1] = a23 + a01;

   q0 = q0*q0;  q1 = q1*q1;  q2 = q2*q2;  q3 = q3*q3;

   rot[0][0] = q0 + q1 - q2 - q3;
   rot[1][1] = q0 - q1 + q2 - q3;
   rot[2][2] = q0 - q1 - q2 + q3;
}
/******************************************************************************
 * rot_to_q. Inverse of above.  Will fall over badly If rot is not orthogonal.*
 ******************************************************************************/
void	rot_to_q(rot, quat)
mat_mt	rot;				/* Rotation matrix		 (in) */
quat_mt	quat;				/* Input quaternion             (out) */
{
   int i, j, k;
   real sign;
   double eps = 100*precision();	/* Safety margin for square roots     */
   
   quat[0] = 0.5 * sqrt(1.0 + eps + rot[0][0] + rot[1][1] + rot[2][2]);

   for( i = 0, j = 1, k = 2; i < 3; i++, j=(i+1)%3, k=(j+1)%3 )
   {
      if( quat[0] > 0 )
         sign = (rot[j][k] - rot[k][j]) >= 0.0 ? 1.0 : -1.0;
      else
         sign = (rot[i][j] + rot[j][i] >= 0) ^ (rot[i][k] + rot[k][i] >= 0)
	        ? -1.0 : 1.0;

      quat[i+1] = sign*0.5*sqrt(1.0+eps + rot[i][i] - rot[j][j] - rot[k][k]);
   }
}
