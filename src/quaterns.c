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
 * quaterns	Functions for manipulating quaternions			      *
 ******************************************************************************
 *      Revision Log
 *       $Log$
 *       Revision 2.12.4.1  2004/03/01 04:53:56  moldydv
 *       Syswrite now treats non-periodic data (from XYZ and some CSSR files) as single species with initial configuration to be set using skew start.
 *       Options -n and -l added for no of particles and species label, respectively, for such systems.
 *
 *       Revision 2.12  2002/09/19 09:26:30  kr
 *       Tidied up header declarations.
 *       Changed old includes of string,stdlib,stddef and time to <> form
 *
 *       Revision 2.11  2000/12/06 17:45:32  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.10  2000/04/27 17:57:10  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.9  2000/04/26 16:01:02  keith
 *       Dullweber, Leimkuhler and McLachlan rotational leapfrog version.
 *
 *       Revision 2.8  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 *
 *       Revision 2.7  1994/06/08 13:22:31  keith
 *       Null update for version compatibility
 *
 * Revision 2.5  1994/01/18  13:32:54  keith
 * Null update for XDR portability release
 *
 * Revision 2.3  93/10/28  10:28:05  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:49:18  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.2  93/03/09  15:59:06  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.1  89/04/20  16:00:52  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header$";
#endif
/*========================== Library include files ===========================*/
#include <math.h>
/*========================== Program include files ===========================*/
#include "defs.h"
/*========================== External function declarations ==================*/
double precision(void);
/*============================================================================*/
/******************************************************************************
 * Quaternion multiplier.  Multiplies arrays of quaternions p by q to give r  *
 * Can be called with  r the same as p or q.                                  *
 ******************************************************************************/
void q_mul(quat_mp p, quat_mp q, quat_mp r, int n)
       	  			/* First Quaternion array [n][4]        (in)  */
	  			/* Second quaternion array [n][4]       (in)  */
	  			/* Resultant quaternions [n][4]        (out)  */
   	  			/* Number of quaternions in the arrays  (in)  */
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
void q_mul_1(real *p,           /* First Quaternion array [n][4]        (in)  */ 
	     real *q,           /* Second quaternion array [n][4]       (in)  */ 
	     real *r)           /* Resultant quaternions [n][4]        (out)  */
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
void q_conj_mul(quat_mp p,      /* First Quaternion array [n][4]        (in)  */ 
		quat_mp q,      /* Second quaternion array [n][4]       (in)  */ 
		quat_mp r,      /* Resultant quaternions [n][4]        (out)  */ 
		int n)	        /* Number of quaternions in the arrays  (in)  */
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
 * Quaternion multiplier.  Multiplies arrays of quaternions p by q(-1) to     *
 * give r. Can be called with  r the same as p or q.                          *
 ******************************************************************************/
void q_mul_conj(quat_mp p,      /* First Quaternion array [n][4]        (in)  */ 
		quat_mp q,      /* Second quaternion array [n][4]       (in)  */ 
		quat_mp r,      /* Resultant quaternions [n][4]        (out)  */ 
		int n)	        /* Number of quaternions in the arrays  (in)  */
{
   register	int	i;
   register	real	p0, p1, p2, p3;
   register	real	q0, q1, q2, q3;
   
   for(i = 0; i < n; i++)
   {
      p0 =  p[i][0]; p1 =  p[i][1]; p2 =  p[i][2]; p3 =  p[i][3];
      q0 =  q[i][0]; q1 = -q[i][1]; q2 = -q[i][2]; q3 = -q[i][3];

      r[i][0] = p0*q0 - p1*q1 - p2*q2 - p3*q3;
      r[i][1] = p1*q0 + p0*q1 - p3*q2 + p2*q3;
      r[i][2] = p2*q0 + p3*q1 + p0*q2 - p1*q3;
      r[i][3] = p3*q0 - p2*q1 + p1*q2 + p0*q3;
   }
}
/******************************************************************************
 *  q_to_mat  Make the rotation matrix corresponding to quaternion q          *
 ******************************************************************************/
void q_to_rot(real *quat,               /* Input quaternion              (in) */ 
	      mat_mt rot)               /* Rotation matrix              (out) */
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
void	rot_to_q(mat_mt rot,            /* Rotation matrix               (in) */ 
		 real *quat)            /* Input quaternion             (out) */
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
