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
 *       Revision 2.12  2001/05/24 16:34:41  keith
 *       Eliminated some redundant variables
 *
 *       Revision 2.11  2001/05/24 16:26:43  keith
 *       Updated program to store and use angular momenta, not velocities.
 *        - added conversion routines for old restart files and dump files.
 *       Got rid of legacy 2.0 and lower restart file reading code.
 *
 *       Revision 2.10  2001/05/22 14:52:45  keith
 *       Added control param "dont-use-symm-rot" to switch between rotational
 *       leapfrog versions at runtime.
 *
 *       Revision 2.9  2001/02/19 19:36:44  keith
 *       First working version of combined isothermic/isobaric ensemble.
 *       (Previous version was faulty).
 *       Also includes uniform-dilation, constant-pressure mode.
 *
 *       Revision 2.8  2001/02/15 17:25:21  keith
 *       Tidied up code and corrected one error in combined barostat
 *       and thermostat. (The combination is still not tested though).
 *
 *       Revision 2.7  2001/02/13 17:45:08  keith
 *       Added symplectic Parrinello-Rahman constant pressure mode.
 *
 *       Revision 2.6  2000/12/06 17:45:31  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.5  2000/11/09 16:54:12  keith
 *       Updated utility progs to be consistent with new dump format
 *
 *       Revision 2.4  2000/11/06 16:02:06  keith
 *       First working version with a Nose-Poincare thermostat for rigid molecules.
 *
 *       System header updated to include H_0.
 *       Dump performs correct scaling  of angular velocities, but dumpext still
 *          needs to be updated to read this.
 *       XDR functions corrected to work with new structs.
 *       Parallel broadcast of config also updated.
 *       Some unneccessary functions and code deleted.
 *
 *       Revision 2.3  2000/05/23 15:23:08  keith
 *       First attempt at a thermostatted version of the Leapfrog code
 *       using either a Nose or a Nose-Poincare thermostat
 *
 *       Revision 2.2  2000/04/27 17:57:08  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.1  2000/04/26 16:03:54  keith
 *       Dullweber, Leimkuhler and McLachlan rotational leapfrog version.
 *
 */
#ifndef lint
static char *RCSid = "$Header: /home/minphys2/keith/CVS/moldy/src/leapfrog.c,v 2.12 2001/05/24 16:34:41 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include <math.h>
#ifdef DEBUG_THERMOSTAT1
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
double  vdot(int n, real *x, int ix, real *y, int iy); /* Fast  dot product   */
void	vscale(int n, double s, real *x, int ix);
void    transpose(mat_mt a, mat_mt b)    ;
void    invert(mat_mt a, mat_mt b)    ;
void    mvaxpy(int n, mat_mt a, vec_mt (*x), vec_mt (*y));
void	mat_mul(mat_mt a, mat_mt b, mat_mt c); /* 3 x 3 matrix multiplier     */
void	mat_add(mat_mt a, mat_mt b, mat_mt c); /* Add 2 3x3 matrices          */
void	mat_sca_mul(register real s, mat_mt a, mat_mt b); 
					/* Multiply 3x3 matrix by scalar      */
double rot_ke(quat_mt (*amom), real s, real *inertia, int nmols); 
                                    /* Compute rotational kinetic energy     */
double  ke_cell(mat_mt hmom, real w);
double  det(mat_mt );			/* Returns matrix determinant	     */
double  trace(mat_mt);
double  trace_sqr(mat_mt);
void    mk_sigma(mat_mt h, mat_mt sigma);
void	tfree(gptr *p);

void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
/*========================== External data references ========================*/
extern	contr_mt	control;            /* Main simulation control parms. */
/*========================== Macros ==========================================*/
#define	TOLERANCE	1.0e-4
#ifndef INERTIA_MIN
#define INERTIA_MIN	1.0e-14		/* Tolerance for zero mom of I	      */
#endif
#ifndef DBL_MAX
#   define DBL_MAX 1.0e37
#endif
#ifndef FLT_MIN
#   define FLT_MIN 1.0e-37
#endif
#ifndef DBL_MIN
#   define DBL_MIN 1.0e-37
#endif
#define REAL_MIN   (sizeof(real) == sizeof(double)?DBL_MIN:FLT_MIN)
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
 * leapf_com(). Perform the centre-of-mass update of the leapfrog integration *
 ******************************************************************************/
void leapf_com(double step, vec_mt (*c_of_m), vec_mt (*mom), 
	       mat_mt h, real s, real mass, int nmols)
{
   mat_mt G, G_inv,  h_tr;
   transpose(h,h_tr);
   mat_mul(h_tr,h,G);                           /* We now have the G matrix   */
   invert(G, G_inv);                            /* G (-1) done                */
   mat_sca_mul(step/(mass*s), G_inv, G_inv);

   mvaxpy(nmols, G_inv, mom, c_of_m);

   escape(c_of_m, nmols);
}
/******************************************************************************
 * leapf_mom().  Perform the linear momentum  update of the integration       *
 ******************************************************************************/
void leapf_mom(double step, real (*h)[3], vec_mt (*mom), vec_mt (*force), 
	       int nmols)
{
   mat_mt h_tr;

   transpose(h, h_tr);
   mat_sca_mul(step, h_tr, h_tr);
   mvaxpy(nmols, h_tr, force, mom);
}
/******************************************************************************
 * leapf_s().  Perform the thermostat co-ordinate update of the leapfrog.     *
 ******************************************************************************/
double leapf_s(double step, real s, real smom, double Q)
{
   double r = 0.5*step*smom/Q;

   s *= (1.0 + r)/(1.0 - r);
#ifdef DEBUG_THERMOSTAT1
   fprintf(stderr,"s = %f\n",*s);
#endif
   return s;
}
/******************************************************************************
 * leapf_smom_a().  Perform the thermostat momentum update of the leapfrog.   *
 *    Eq. 35b in Bond, Leimkuhler et al (1999)				      *
 ******************************************************************************/
double leapf_smom_a(double step, real s, real smomo,  double Q, double gkt)
{
   double C, smomn;

   C = -smomo+0.5*step*gkt*(log(s)+1);

   smomn = -2.0*C/(1.0 + sqrt(1.0-C*step/Q));
#ifdef DEBUG_THERMOSTAT1
   fprintf(stderr,"leapf_smom: C=%f\t(1-hC/Q)=%f\tSmom %f -> %f\tDelta=%20.16g\tstep=%f\n", 
	   C, 1.0-C*step/Q, smomo, smomn, 0.25*step/Q*SQR(smomn)+(smomn)+C,step);
#endif
   return smomn;
}
/******************************************************************************
 * leapf_smom_b().  Perform the thermostat momentum update of the leapfrog.   *
 *    Eq. 35e in Bond, Leimkuhler et al (1999)				      *
 ******************************************************************************/
double leapf_smom_b(double step, real s, real smomo, double Q, double gkt)
{
   double  smomn;

   smomn = smomo  -0.5*step*(SQR(smomo)/(2.0*Q)+gkt*(log(s)+1));

  return smomn;
}
/******************************************************************************
 * make_rot()  Construct the quaternion representing a rotation about a given *
 *    axis (x,y,z) by angle h*amom_i.                                         *
 ******************************************************************************/
void make_rot(double step, int axis, quat_mt (*amom), real rinertia, 
	      quat_mt (*rot), int nmols)
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
void make_rot_amom(double step, quat_mt (*amom), real rinertia, 
		   quat_mt (*rot), int nmols)
{
   int imol;
   double samom, ramom, angle, ca, sa;
   
   for(imol = 0; imol < nmols; imol++)
   {
      samom=sqrt(SUMSQ2(amom[imol]));
      ramom = 1.0/(samom + (8*DBL_MIN));	/* Don't fail is amom == 0           */
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
void leapf_quat_b(double step, quat_mt (*quat), quat_mt (*amom), 
		  real *inertia, real *smom, real ts, int nmols)
{
   int i;
   double idmin, idiff;
   double       stepdts = step/ts;

   quat_mt	*rot  = qalloc(nmols);
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
    * First, apply half of delta-from-symmetry rotation
    */
   make_rot(0.5*stepdts, orthaxis1, amom, 
	    rinertia[orthaxis1]-rinertia[orthaxis2], rot, nmols);
   rot_substep(rot, amom, quat, nmols);
   /*
    * Update thermostat momentum if necessary.
    */
   if( control.const_temp )
     *smom += 0.5*step*rot_ke(amom, ts, inertia, nmols);
   /*
    * Apply rotations to angular momentum and to quaternions.
    * 1. Rotation about symmetry axis due to precession.
    */
   make_rot(stepdts, saxis, amom, rinertia[saxis]-rinertia[orthaxis2], 
	    rot, nmols);
   rot_substep(rot, amom, quat, nmols);
   /*
    * 2. Rotation about angular momentum vector due to spin.
    */
   make_rot_amom(stepdts, amom, rinertia[orthaxis2],rot, nmols);
   q_mul(quat, rot, quat, nmols);
   /*
    * Update thermostat momentum if necessary.
    */
   if( control.const_temp )
     *smom += 0.5*step*rot_ke(amom, ts, inertia, nmols);
   /*
    * Finally, apply half of delta-from-symmetry rotation
    */
   make_rot(0.5*step, orthaxis1, amom, 
	    rinertia[orthaxis1]-rinertia[orthaxis2], rot, nmols);
   rot_substep(rot, amom, quat, nmols);

   normalise(quat, nmols);
   xfree(rot);
}
/******************************************************************************
 * leapf_quat().  Step quaternions using leapfrog integrator.  Angular momenta*
 *    are also updated.  This routine uses the symplectic splitting method    *
 *    outlined in Dullweber et al (1997) J. Chem. Phys 107, 5840-5851.        *
 *    Adapted for use with quaternions by K.R.				      *
 ******************************************************************************/
void leapf_quat_a(double step, quat_mt (*quat), quat_mt (*amom), 
		  real *inertia, real *smom, real ts, int nmols)
{
   int i;
   quat_mt	*rot  = qalloc(nmols);
   double       stepdts = step/ts, stepdts2 = step/SQR(ts);

   real rinertia[3];

   for(i=0; i<3; i++)
   {
      if( inertia[i]/(inertia[(i+1)%3]+inertia[(i+2)%3]) < INERTIA_MIN )
	 rinertia[i] = 0.0;
      else
	 rinertia[i] = 1.0/inertia[i];
   }

   /*
    * Apply rotations to angular momentum and to quaternions.
    */
   make_rot(0.5*stepdts, 0, amom, rinertia[0], rot, nmols);
   rot_substep(rot, amom, quat, nmols);
   if( control.const_temp )
   {
     *smom += 0.25*stepdts2*rinertia[0]*vdot(nmols, amom[0]+1, 4, amom[0]+1, 4);
     *smom += 0.25*stepdts2*rinertia[1]*vdot(nmols, amom[0]+2, 4, amom[0]+2, 4);
   }
   make_rot(0.5*stepdts, 1, amom, rinertia[1], rot, nmols);
   rot_substep(rot, amom, quat, nmols);
   if( control.const_temp )
   {
     *smom += 0.5*stepdts2*rinertia[2]*vdot(nmols, amom[0]+3, 4, amom[0]+3, 4);
   }
   make_rot(    stepdts, 2, amom, rinertia[2], rot, nmols);
   rot_substep(rot, amom, quat, nmols);
   make_rot(0.5*stepdts, 1, amom, rinertia[1], rot, nmols);
   rot_substep(rot, amom, quat, nmols);
   if( control.const_temp )
   {
     *smom += 0.25*stepdts2*rinertia[1]*vdot(nmols, amom[0]+2, 4, amom[0]+2, 4);
     *smom += 0.25*stepdts2*rinertia[0]*vdot(nmols, amom[0]+1, 4, amom[0]+1, 4);
   }
   make_rot(0.5*stepdts, 0, amom, rinertia[0], rot, nmols);
   rot_substep(rot, amom, quat, nmols);

   normalise(quat, nmols);
   xfree(rot);
}
/******************************************************************************
 * leapf_quat().  Chooses which symplectic splitting method to use.           *
 ******************************************************************************/
void leapf_quat(double step, quat_mt (*quat), quat_mt (*amom), 
		real *inertia, real *smom, real ts,  int nmols)
{
  if( control.nosymmetric_rot )
      leapf_quat_a(step, quat, amom, inertia, smom, ts, nmols);
  else
      leapf_quat_b(step, quat, amom, inertia, smom, ts, nmols);
}
/******************************************************************************
 * leapf_amom().  Perform the angular velocity update step of the leapfrog    *
 *    integration algorithm.                                                  *
 ******************************************************************************/
void leapf_amom(double step, quat_mt (*amom), vec_mt (*torque), int nmols)
{
   int imol;

   for(imol = 0; imol < nmols; imol++)
   {
      amom[imol][1] += step*torque[imol][0];
      amom[imol][2] += step*torque[imol][1];
      amom[imol][3] += step*torque[imol][2];
   }
   
}
/******************************************************************************
 * leapf_hmom().  Perform the h-matrix momentum update step.	              *
 ******************************************************************************/
void leapf_hmom(double step, mat_mt hmom, mat_mt sigma, real s, real pressure, 
		int mask)
{
   int i,j;
   double r = 0.5 * step * s * pressure;

   for(i=0; i<3; i++)
      for (j=0; j<3; j++)
      {
	 if( mask & 1 )
	    hmom[i][j] = 0;
	 else
	    hmom[i][j] -= r *sigma[i][j];
	 mask  >>= 1;    	 
      }
}
/******************************************************************************
 * leapf_h().  Perform the h-matrix update step.			      *
 ******************************************************************************/
void leapf_h(double step, mat_mt h, mat_mt hmom, real s, real w)
{
   int i,j;
   double r = step * s / w;

   for(i=0; i<3; i++)
      for (j=0; j<3; j++)
	 h[i][j] += r * hmom[i][j];
}
/******************************************************************************
 * gleap_therm.  Update thermostat variable and momenta using generalized leap*
 ******************************************************************************/
void gleap_therm(double step, real mass, real gkt, real *s, real *smom)
{
   *smom = leapf_smom_a(step, *s, *smom, mass, gkt );
   *s    = leapf_s     (step, *s, *smom, mass);
   *smom = leapf_smom_b(step, *s, *smom, mass, gkt);
}
/******************************************************************************
 * gleap_cell. Update cell variable and momenta using generalized leapfrog    *
 ******************************************************************************/
void gleap_cell(double step, real pmass, real s, real pressure, int strain_mask,
		mat_mt h, mat_mt hmom, real *smom, int uniform)
{
  double vol;
  mat_mt sigma;

  vol = det(h);

  if( uniform ) 
     leapf_hmom(step, hmom, h, 3.0*vol*s/trace_sqr(h), pressure, strain_mask);
  else
  {
     mk_sigma(h, sigma);
     leapf_hmom(step, hmom, sigma, s, pressure, strain_mask);
  }

  if( control.const_temp )
     *smom -= 0.5*step*(ke_cell(hmom, pmass) + pressure*vol);
  
  leapf_h(step, h, hmom, s, pmass);
  
  vol = det(h);
  if( control.const_temp )
     *smom -= 0.5*step*(ke_cell(hmom, pmass) + pressure*vol);

  if( uniform ) 
     leapf_hmom(step, hmom, h, 3.0*vol*s/trace_sqr(h), pressure, strain_mask);
  else
  {
     mk_sigma(h, sigma);
     leapf_hmom(step, hmom, sigma, s, pressure, strain_mask);
  }
}
/******************************************************************************
 * update_hmom().  Perform "simple" update of unit cell momenta.  This is used*
 *                 to include all of the stress terms.			      *
 * N.B. It is assumed that stress_part is assocuated with an actual parameter *
 *      matric containing either the kinetic or virial part of the stress     *
 *      tensor MULTIPLIED BY THE VOLUME.				      *
 ******************************************************************************/
void update_hmom(double step, real s, mat_mt h,
		 mat_mt stress_part, mat_mt hmom, int uniform)
{
   double vol;
   mat_mt sigma;
   mat_mt tmp_mat;

   if( uniform ) {
      mat_sca_mul(step*s*trace(stress_part)/trace_sqr(h), h, tmp_mat);
   } else {
      vol = det(h);
      mk_sigma(h, sigma);
      
      mat_mul(stress_part, sigma, tmp_mat);
      mat_sca_mul(step*s/vol, tmp_mat, tmp_mat);
   }
   mat_add(hmom, tmp_mat, hmom);
}
