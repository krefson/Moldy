#define DEBUG_THERMOSTAT
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
 *   Leapfrog - Routines to implement the leapfrog algorithm for              *
 *              stepping the co-ordinates in a MD simulation including rigid  *
 *              molecules                                                     *
 *              The rotational part is an implementation of the symplectic    *
 *              integrator described in Dullweber et al (1977) J. Chem. Phys  *
 *              107, 5840-5851.  Goodness knows why they thought you couldn't *
 *              continue to use quaternions to represent the rotations..      *
 *   External routines:          step_1(system), step_2(system), beeman_2()   *
 *   External references:        none                                         *
 *   External data:              none                                         *
 *   External data references:   control                                      *
 *                                                                            *
 ******************************************************************************
 *      Revision Log
 *       $Log: leapfrog.c,v $
 *       Revision 2.2  2000/04/27 17:57:08  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.1  2000/04/26 16:03:54  keith
 *       Dullweber, Leimkuhler and McLachlan rotational leapfrog version.
 *
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore_data/keith/CVS/moldy/src/leapfrog.c,v 2.2 2000/04/27 17:57:08 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include <math.h>
#ifdef DEBUG_THERMOSTAT
#   define DEBUG
#endif
#ifdef DEBUG
#   include <stdio.h>
#endif
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
gptr	*talloc(int n, size_mt size, int line, char *file);
void	q_mul(quat_mp p, quat_mp q, quat_mp r, int n);
void	q_conj_mul(quat_mp p, quat_mp q, quat_mp r, int n);
void	q_mul_conj(quat_mp p, quat_mp q, quat_mp r, int n);
void	vscale(register int n, register double s, register real *x, register int ix);
void	tfree(gptr *p);

void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
/*========================== External data references ========================*/
extern	contr_mt	control;            /* Main simulation control parms. */
/*========================== Macros ==========================================*/
#define	TOLERANCE	1.0e-4
#define INERTIA_MIN	1.0e-14		/* Tolerance for zero mom of I	      */
#ifndef DBL_MAX
#   define DBL_MAX 1.0e37
#endif
#ifndef FLT_MIN
#   define FLT_MIN 1.0e-37
#endif
/*============================================================================*/
/******************************************************************************
 *   Normalise the new quaternions                                            *
 ******************************************************************************/
static void normalise(quat_mp quat, int n)
   		  			/* Number of quaternions  (in)        */
       		     			/* Quaternions       (update)         */
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
 *   Apply periodic boundary conditions to put particles back in MD box       *
 ******************************************************************************/
void escape(vec_mp c_of_m, int nmols)
   		      		/* First dimension of c-of-m                  */
      		       		/* Centre of mass co-ordinates (updat)        */

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
 * leapf_coml(). Perform the centre-of-mass update of the leapfrog integration *
 ******************************************************************************/
void leapf_com(real step, vec_mt (*c_of_m), vec_mt (*vel), int nmols)
{
   int imol;
   real *com = c_of_m[0], *vm=vel[0];

   for(imol = 0; imol < 3*nmols; imol++)
      com[imol] += step * vm[imol];
   escape(c_of_m, nmols);
}
/******************************************************************************
 * leapf_vel().  Perform the velocity update of the leapfrog integration      *
 ******************************************************************************/
void leapf_vel(real step, real (*hinv)[3], vec_mt (*vel), vec_mt (*force), real mass, int nmols)
{
   int	imol;
   real	srmass = step/mass;
   real hi00=hinv[0][0], hi01=hinv[0][1], hi02=hinv[0][2];
   real hi10=hinv[1][0], hi11=hinv[1][1], hi12=hinv[1][2];
   real hi20=hinv[2][0], hi21=hinv[2][1], hi22=hinv[2][2];

   for(imol = 0; imol < nmols; imol++)
   {
      vel[imol][0] += srmass*(hi00*force[imol][0] + hi01*force[imol][1] 
                                                  + hi02*force[imol][2]);
      vel[imol][1] += srmass*(hi10*force[imol][0] + hi11*force[imol][1] 
                                                  + hi12*force[imol][2]);
      vel[imol][2] += srmass*(hi20*force[imol][0] + hi21*force[imol][1] 
                                                  + hi22*force[imol][2]);
   }
}
/******************************************************************************
 * leapf_s().  Perform the thermostat co-ordinate update of the leapfrog.     *
 ******************************************************************************/
void leapf_s(double step, real *s, real smom, double Q)
{
   double r = 0.5*step*smom/Q;

   *s *= (1.0 + r)/(1.0 - r);
#ifdef DEBUG_THERMOSTAT1
   fprintf(stderr,"s = %f\n",*s);
#endif
}
/******************************************************************************
 * leapf_smom().  Perform the thermostat momentum update of the leapfrog.     *
 *   This version performs a whole-step update, combining Eq. 35e and 35b.    *
 ***************************************x***************************************/
void leapf_smom(double step, real s, real *smom, double kepold, double pe, double Q, double gkt, double H_0)
{
   double C, smomo=*smom, smomn;

   C = -smomo-0.5*step*(kepold - 2.0*(pe + gkt*(log(s)+1) - H_0) - 0.5*SQR(smomo)/Q);

   smomn = -2.0*C/(1.0 + sqrt(1.0-C*step/Q));
#ifdef DEBUG_THERMOSTAT1
   fprintf(stderr,"leapf_smom: C=%f\t(1-hC/Q)=%f\tSmom %f -> %f\tDelta=%20.16g\tstep=%f\n", 
	   C, 1.0-C*step/Q, smomo, smomn, 0.25*step/Q*SQR(smomn)+(smomn)+C,step);
#endif
   *smom = smomn;
}
/******************************************************************************
 * leapf_smom_a().  Perform the thermostat momentum update of the leapfrog.   *
 *    Eq. 35b in Bond, Leimkuhler et al (1999)				      *
 ******************************************************************************/
void leapf_smom_a(double step, real s, real *smom, double ke, double pe, double Q, double gkt, double H_0)
{
   double C, smomo=*smom, smomn;

   C = -smomo-0.5*step*(ke - pe - gkt*(log(s)+1) + H_0);

   smomn = -2.0*C/(1.0 + sqrt(1.0-C*step/Q));
#ifdef DEBUG_THERMOSTAT1
   fprintf(stderr,"leapf_smom: C=%f\t(1-hC/Q)=%f\tSmom %f -> %f\tDelta=%20.16g\tstep=%f\n", 
	   C, 1.0-C*step/Q, smomo, smomn, 0.25*step/Q*SQR(smomn)+(smomn)+C,step);
#endif
   *smom = smomn;
}
/******************************************************************************
 * leapf_smom_b().  Perform the thermostat momentum update of the leapfrog.   *
 *    Eq. 35e in Bond, Leimkuhler et al (1999)				      *
 ******************************************************************************/
void leapf_smom_b(double step, real s, real *smom, double ke, double pe, double Q, double gkt, double H_0)
{
   double  smomo=*smom, smomn;

   smomn = smomo + 0.5*step*(2.0*ke-gkt)
      -0.5*step*(ke+pe+SQR(smomo)/(2.0*Q)+gkt*log(s) - H_0);

   *smom = smomn;
}
/******************************************************************************
 * leapf_smom().  Perform the thermostat momentum update of the leapfrog.     *
 *                in the straighforward Nose Hamiltonian case.		      *
 ***************************************x***************************************/
void leapf_smom_simple(double step, real *smom, double ke, double gkt)
{
   *smom += step*(ke - gkt);
}
/******************************************************************************
 * make_rot()  Construct the quaternion representing a rotation about a given *
 *    axis (x,y,z) by angle h*omega_i.                                        *
 ******************************************************************************/
void make_rot(double step, int axis, quat_mt (*amom), real rinertia, quat_mt (*rot), int nmols)
{
   int imol;
   double angle, ca, sa;
   
   for(imol = 0; imol < nmols; imol++)
   {
      angle = 0.5*step*rinertia*amom[imol][axis+1];
      ca = cos(angle); sa = sin(angle);
      rot[imol][0] = ca;
      rot[imol][1] = rot[imol][2] = rot[imol][3] = 0.0;
      rot[imol][axis+1] = sa;
   }
}
/******************************************************************************
 * make_rot_amom().  Construct the quaternions representing rotation about the*
 *   axis parallel to the angular momentum by angle h*|L|/I                   *
 ******************************************************************************/
void make_rot_amom(double step, quat_mt (*amom), real rinertia, quat_mt (*rot), int nmols)
{
   int imol;
   double samom, ramom, angle, ca, sa;
   
   for(imol = 0; imol < nmols; imol++)
   {
      samom=sqrt(SUMSQ2(amom[imol]));
      ramom = 1.0/(samom + FLT_MIN);	/* Don't fail is amom == 0           */
      angle = 0.5*step*rinertia*samom;
      ca = cos(angle); sa = sin(angle);
      rot[imol][0] = ca;
      rot[imol][1] = sa * ramom * amom[imol][1];
      rot[imol][2] = sa * ramom * amom[imol][2];
      rot[imol][3] = sa * ramom * amom[imol][3];
   }
}
/******************************************************************************
 * rot_substep().  Perform a combined rotation of the angular momentum and    *
 *     quaternions, implementing Equation 6 of Dullweber et al.  N.B. We use  *
 *     the opposite conversion and so our rotations are inverse to theirs.    *
 ******************************************************************************/
void rot_substep(quat_mt (*rot), quat_mt (*amom), quat_mt (*quat), int nmols)
{
   q_conj_mul(rot, amom, amom, nmols);
   q_mul(amom, rot, amom, nmols);
   
   q_mul(quat, rot, quat, nmols);
}
/******************************************************************************
 * leapf_quat().  Step quaternions using leapfrog integrator.  Angular momenta*
 *    are also updated.  This routine uses the method outlined in Appendix C  *
 *    of Dullweber et al (1997) J. Chem. Phys 107, 5840-5851..                *
 *    The Hamiltonion is split into a symmetric and asymmetric part, and is   *
 *    exact for the case of a free rotor with two equal moments of inertia.   *
 *    Adapted for use with quaternions by K.R.				      *
 ******************************************************************************/
void leapf_quat_b(real step, quat_mt (*quat), quat_mt (*avel), real *inertia, int nmols)
{
   int imol,i;
   double idmin, idiff;
   quat_mt	*rot  = qalloc(nmols);
   real	inx = inertia[0], iny = inertia[1], inz = inertia[2];	
   real rinertia[3];
   static int firsttime = 0, saxis, orthaxis1,orthaxis2;

   for(i=0; i<3; i++)
   {
      if( inertia[i]/(inertia[(i+1)%3]+inertia[(i+2)%3]) < INERTIA_MIN )
	 rinertia[i] = 0.0;
      else
	 rinertia[i] = 1.0/inertia[i];
   }
   /*
    * Find near-symmetry axis
    */
   if( firsttime == 0) 
   {
      idmin=DBL_MAX;
      for( i=0; i<3; i++)
      {
	 idiff = fabs(rinertia[(i+1)%3]-rinertia[(i+2)%3]);
	 if (idiff < idmin )
	 {
	    idmin = idiff;
	    saxis = i;
	 }
	 /*note("1/I%d-1/I%d = %f",(i+1)%3,(i+2)%3,idiff);*/
      }
      /* note("Chose %d",saxis);*/
      firsttime++;
      orthaxis1 = (saxis+1)%3;
      orthaxis2 = (saxis+2)%3;
   }
   /*
    *  Compute angular momentum.
    */
   for(imol = 0; imol < nmols; imol++)
   {
      avel[imol][1] *= inx;
      avel[imol][2] *= iny;
      avel[imol][3] *= inz;
   }
   /*
    * First, apply half of delta-from-symmetry rotation
    */
   make_rot(0.5*step, orthaxis1, avel, 
	    rinertia[orthaxis1]-rinertia[orthaxis2], rot, nmols);
   rot_substep(rot, avel, quat, nmols);
   /*
    * Apply rotations to angular momentum and to quaternions.
    * 1. Rotation about symmetry axis due to precession.
    */
   make_rot(step, saxis, avel, rinertia[saxis]-rinertia[orthaxis2], 
	    rot, nmols);
   rot_substep(rot, avel, quat, nmols);
   /*
    * 2. Rotation about angular momentum vector due to spin.
    */
   make_rot_amom(step, avel, rinertia[orthaxis2],rot, nmols);
   q_mul(quat, rot, quat, nmols);
   /*
    * Finally, apply half of delta-from-symmetry rotation
    */
   make_rot(0.5*step, orthaxis1, avel, 
	    rinertia[orthaxis1]-rinertia[orthaxis2], rot, nmols);
   rot_substep(rot, avel, quat, nmols);

   /*
    * Compute updated angular velocity from angular momentum.
    */
   for(imol = 0; imol < nmols; imol++)
   {
      avel[imol][1] *= rinertia[0];
      avel[imol][2] *= rinertia[1];
      avel[imol][3] *= rinertia[2];
   }
   normalise(quat, nmols);
   xfree(rot);
}
/******************************************************************************
 * leapf_quat().  Step quaternions using leapfrog integrator.  Angular momenta*
 *    are also updated.  This routine uses the symplectic splitting method    *
 *    outlined in Dullweber et al (1997) J. Chem. Phys 107, 5840-5851.        *
 *    Adapted for use with quaternions by K.R.				      *
 ******************************************************************************/
void leapf_quat_a(real step, quat_mt (*quat), quat_mt (*avel), real *inertia, int nmols)
{
   int imol,i;
   quat_mt	*rot  = qalloc(nmols);
   real	inx = inertia[0], iny = inertia[1], inz = inertia[2];	
   real rinertia[3];

   for(i=0; i<3; i++)
   {
      if( inertia[i]/(inertia[(i+1)%3]+inertia[(i+2)%3]) < INERTIA_MIN )
	 rinertia[i] = 0.0;
      else
	 rinertia[i] = 1.0/inertia[i];
   }

   /*
    *  Compute angular momentum.
    */
   for(imol = 0; imol < nmols; imol++)
   {
      avel[imol][1] *= inx;
      avel[imol][2] *= iny;
      avel[imol][3] *= inz;
   }
   /*
    * Apply rotations to angular momentum and to quaternions.
    */
   make_rot(0.5*step, 0, avel, rinertia[0], rot, nmols);
   rot_substep(rot, avel, quat, nmols);
   make_rot(0.5*step, 1, avel, rinertia[1], rot, nmols);
   rot_substep(rot, avel, quat, nmols);
   make_rot(    step, 2, avel, rinertia[2], rot, nmols);
   rot_substep(rot, avel, quat, nmols);
   make_rot(0.5*step, 1, avel, rinertia[1], rot, nmols);
   rot_substep(rot, avel, quat, nmols);
   make_rot(0.5*step, 0, avel, rinertia[0], rot, nmols);
   rot_substep(rot, avel, quat, nmols);
   /*
    * Compute updated angular velocity from angular momentum.
    */
   for(imol = 0; imol < nmols; imol++)
   {
      avel[imol][1] *= rinertia[0];
      avel[imol][2] *= rinertia[1];
      avel[imol][3] *= rinertia[2];
   }
   normalise(quat, nmols);
   xfree(rot);
}
/******************************************************************************
 * leapf_quat().  Chooses which symplectic splitting method to use.           *
 ******************************************************************************/
void leapf_quat(real step, quat_mt (*quat), quat_mt (*avel), real *inertia, int nmols)
{
#define SYMM_BODY
#ifdef SYMM_BODY
      leapf_quat_b(step, quat, avel, inertia, nmols);
#else
      leapf_quat_a(step, quat, avel, inertia, nmols);
#endif
}
/******************************************************************************
 * leapf_avel().  Perform the angular velocity update step of the leapfrog    *
 *    integration algorithm.                                                  *
 ******************************************************************************/
void leapf_avel(real step, quat_mt (*avel), vec_mt (*torque), real *inertia, int nmols)
{
   int imol, i;
   real rinertia[3];
   
   for(i=0; i<3; i++)
   {
      if( inertia[i]/(inertia[(i+1)%3]+inertia[(i+2)%3]) < INERTIA_MIN )
	 rinertia[i] = 0.0;
      else
	 rinertia[i] = 1.0/inertia[i];
   }

   for(imol = 0; imol < nmols; imol++)
   {
      avel[imol][1] += step*rinertia[0]*torque[imol][0];
      avel[imol][2] += step*rinertia[1]*torque[imol][1];
      avel[imol][3] += step*rinertia[2]*torque[imol][2];
   }
   
}