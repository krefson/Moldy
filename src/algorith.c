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
 * Algorith  This file contains functions to implement the simulation algor - *
 *           ithms and other functions which do not need access to the system *
 *           and species structs.  Contents:                     	      *
 * rotate()		Perform co-ordinate transformation from quaternions   *
 * Vec_dist()		Return 'distance' between 2 long vectors	      *
 * mol_force()		Calculate molecular centre of mass forces	      *
 * mol_torque()		Calculate molecular torques			      *
 * make_sites()		Generate atomic site co-ordinates from c-of-m etc     *
 * newton()		Calculate accelerations from forces		      *
 * euler()		Calculate quaternion accelerations from torques	      *
 * parinello()		Calculate P&R c-of-m acceleration term		      *
 * rahman()		Calculate unit cell matrix accelerations from stress  *
 * trans_ke()		Return translational kinetic energy.		      *
 * rot_ke()		Return rotational kinetic energy		      *
 * energy_dyad()	Calculate kinetic energy part of stress tensor	      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	algorith.c,v $
 * Revision 2.0  93/03/15  14:48:51  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.1.1.12  93/03/15  14:41:28  keith
 * Added GPL copyleft notice to permit release and distribution.
 * N.B.  Previous versions were copyright (C) by author and 
 * only licensed by explicit permission.
 * 
 * Revision 1.1.1.11  93/03/09  15:58:12  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.1.1.10  91/08/15  18:11:39  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.1.1.9  91/03/12  15:42:03  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.1.1.8  90/10/22  16:41:43  keith
 * Make vec_dist() robust in case of all-zero vectors.
 * 
 * Revision 1.1.1.7  90/09/28  13:28:09  keith
 * Inserted braces around VECTORIZE directives and changed include files
 * for STARDtardent 3000 series (via cond. comp symbol "ardent").
 * 
 * Revision 1.1.1.6  90/07/16  15:55:25  keith
 * Fixed bugs in constant-stress code
 * 
 * Revision 1.1.1.5  90/05/16  18:39:06  keith
 * *** empty log message ***
 * 
 * Revision 1.1.1.4  90/04/16  18:19:29  keith
 * Modified rahman() to arbitrarily constrain h matrix by new parameter "mask".
 * 
 * Revision 1.1.1.3  89/12/21  16:29:43  keith
 * Reversed indices in 'site' and 'site_force' to allow stride of 1 in ewald.
 * 
 * Revision 1.1.1.2  89/10/24  17:17:25  keith
 * Modified pbc algorithm to use floor() library function.
 * Now works with non-orthorhombic cell.
 * 
 * Revision 1.1.1.1  89/10/06  16:23:57  keith
 * Make_sites() modified to wrap sites of framework back into MD box.
 * 
 * Revision 1.1  89/04/20  16:00:19  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/algorith.c,v 2.0 93/03/15 14:48:51 keith Rel $";
#endif
/*========================== program include files ===========================*/
#include 	"defs.h"
#include 	"messages.h"
/*========================== Library include files ===========================*/
#include 	<math.h>
#include 	"string.h"
#include	"stddef.h"
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
void	mat_vec_mul();			/* 3 x 3 Matrix by Vector multiplier  */
void	mat_mul();	          	/* 3 x 3 matrix multiplier	      */
void	mat_add();			/* Add 2 3x3 matrices                 */
void	mat_sca_mul();			/* Multiply 3x3 matrix by scalar      */
void	transpose();			/* transpose a 3x3 matrix	      */
void	invert();			/* invert a 3x3 matrix		      */
double	det();				/* Determinant of 3x3 matrix	      */
void	q_to_rot();			/* Make rotation matrix from quat'n   */
void	q_mul();
void	q_conj_mul();
double	vdot();				/* Vector dot product		      */
double	sum();				/* Vector sum			      */
void	vscale();
#if defined(ANSI) || defined(__STDC__)
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
#else
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
#endif
/*========================== Macros ==========================================*/
#define MATMUL(i, m, r, o) (m[i][0]*r[0][o] + m[i][1]*r[1][o] + m[i][2]*r[2][o])
/*============================================================================*/
/******************************************************************************
 *  rotate        Perform the rotation described by the quaternions in the    *
 *  second parameter on each of the co-ordinates in the first, putting the    *
 *  result in the third.  (Pawley,  Mol. Phys. 43, 1321-1330 (1981))          *
 *  NB this is different to Evans' formulation.                               *
 *  Apply each rotation to nvec/nquat vectors.                                *
 ******************************************************************************/
void rotate(r_in, r_out, nvec, quat, nquat, inv_mat)
vec_mp		r_in,		/* Co-ordinates to be rotated [n][3] (in)     */
		r_out;		/* Resulting co-ordinates [n][3]    (out)     */
quat_mp		quat;		/* Quaternions for the rotation.     (in)     */
int		nvec,		/* Number of co-ordinates.           (in)     */
		nquat;		/* Number of quaternions             (in)     */
invrot		inv_mat;	/* Flag to do inverse rotations      (in)     */
{
   mat_mt	rot;
   int		iquat;

   if(nvec % nquat != 0)
      message(NULLI, NULLP, FATAL, ROTLEN, nvec, nquat);

   for(iquat = 0; iquat < nquat; iquat++)
   {
      q_to_rot(quat[iquat], rot);
      if(inv_mat == inv)    transpose(rot, rot);
      mat_vec_mul(rot, r_in, r_out, nvec / nquat);
      r_in += nvec/nquat;   r_out += nvec/nquat;
   }
}
/******************************************************************************
 *  mean_square  Calculate the mean square of list of vectors for each cmpnt  *
 ******************************************************************************/
void	mean_square(x, meansq, nmols)
vec_mt	x[];
vec_mt	meansq;
int	nmols;
{
   int		i;
   for(i = 0; i < 3; i++)
      meansq[i] = vdot(nmols, x[0]+i, 3, x[0]+i, 3) / nmols;
}
/******************************************************************************
 *   vec_dist      Return the normalised distance between 2 vectors           *
 ******************************************************************************/
double	vec_dist(v1, v2, n)
real	*v1, *v2;		/* Input vectors			      */
int	n;			/* Length ie v1[n], v2[n]		      */
{
   register	double s=0, s1=0,s2=0;	/* Accumulators for sums	      */
   while(n-- > 0)
   {
      s  += (*v1 - *v2) * (*v1 - *v2);
      s1 += *v1 * *v1;
      s2 += *v2 * *v2;
      v1++; v2++;
   }
   s1 = MAX(s1,s2);
   if( s1 == 0.0 )
      return s1;
   else
      return(sqrt(s / s1));
}
/******************************************************************************
 *  molecule_force     Calculate the centre of mass forces on a number of     *
 *  molecules given the site forces and the site co-ordinates.                *
 ******************************************************************************/
void mol_force(site_force, force, nsites, nmols)
real		**site_force;	/* Site forces [nsites*nmols][3]        (in)  */
vec_mp		force;		/* Centre of mass forces [nmols][3]    (out)  */
int		nsites,		/* Number of sites on one molecule      (in)  */
		nmols;		/* Number of molecules                  (in)  */
{
   int	i,  imol;

   for(imol = 0; imol < nmols; imol++)
      for(i = 0; i < 3; i++)
         force[imol][i] = sum(nsites, site_force[i]+imol*nsites, 1);
}
/******************************************************************************
 *  molecule_torque    Calculate the torque on a number of identical          *
 *  molecules given the space frame site forces and co-ordinates.             *
 ******************************************************************************/
void mol_torque(site_force, site, torque, quat, nsites, nmols)
real		**site_force;	/* Principal frame site forces          (in)  */
vec_mp		site,		/* Principal frame site co-ordinates    (in)  */
		torque;		/* Molecular torques [nmols][3]        (out)  */
quat_mp		quat;		/* Molecular quaternions [nmols][4]     (in)  */
int		nsites,		/* Number of sites on one molecule      (in)  */
		nmols;		/* Number of molecules                  (in)  */
{
   vec_mp	princ_force = ralloc(nsites);
   int	i, j, k, imol, isite;
   register     double torq;

   for(imol = 0; imol < nmols; imol++)
   {
      for(i = 0; i < 3; i++)
      {
VECTORIZE
	 for(isite = 0; isite < nsites; isite++)
	    princ_force[isite][i] = site_force[i][isite+imol*nsites];
      }
      rotate(princ_force, princ_force, nsites, quat+imol, 1, inv);
      for(i = 0, j = 1, k = 2; i < 3; i++, j=(j+1)%3, k=(k+1)%3)
      {
         torq = 0.0;
VECTORIZE
	 for(isite = 0; isite < nsites; isite++)
	    torq += site[isite][j]*princ_force[isite][k]
	           -site[isite][k]*princ_force[isite][j];
         torque[imol][i] = torq;
      }
   }
   xfree(princ_force);
}
/******************************************************************************
 *  make_sites     Calculate the atomic site co-ordinates for nmols identical *
 *  molecules from the principal-frame sites, the quaternions and the centre  *
 *  of mass co-ordinates.  Called once for each molecular species.            *
 ******************************************************************************/
void make_sites(h, c_of_m_s , quat, p_f_sites, framework, site, nmols, nsites)
mat_mt		h;		/* Unit cell matrix h		     (in)     */
vec_mp		c_of_m_s,	/* Centre of mass co-ords [nmols][3] (in)     */
		p_f_sites;	/* Principal-frame sites [nsites][3] (in)     */
real		**site;		/* Sites [nmols*nsites][3]          (out)     */
int		framework;	/* Flag to signal framework structure (in)    */
quat_mp		quat;		/* Quaternions [nmols][4]            (in)     */
int		nmols,		/* Number of molecules                        */
		nsites;		/* Number of sites on each molecule           */
{
   int		imol, isite, i;	/* Counters				      */
   vec_mt	c_of_m;		/* Unscaled centre of mass co-ordinates       */
   vec_mt	*ssite = ralloc(nsites);
   register double	t;
   mat_mt	hinv;
   double	lx   = h[0][0], lxy  = h[0][1],
		ly   = h[1][1], lxz  = h[0][2],
		lz   = h[2][2], lyz  = h[1][2];
   invert(h,hinv);

   for(imol = 0; imol < nmols; imol++)
   {
      mat_vec_mul(h,c_of_m_s+imol,(vec_mp)c_of_m, 1);/* Get real c-of-m co-ords*/
      if(quat)
      {
         rotate(p_f_sites,ssite,nsites,quat+imol,1,noinv);
	 for(i = 0; i < 3; i++)
	    for(isite = 0; isite < nsites; isite++)
	       site[i][imol*nsites+isite] = ssite[isite][i] + c_of_m[i];
      }
      else
      {
	 for(i = 0; i < 3; i++)
	    for(isite = 0; isite < nsites; isite++)
	       site[i][imol*nsites+isite] = p_f_sites[isite][i] + c_of_m[i];
      }
   }

   if( framework )			/* Apply pbc's to put sites into cell */
      for( isite = 0; isite < nmols*nsites; isite++ )
      {
          site[0][isite] -= lx  *      floor(MATMUL(0,hinv,site,isite) + 0.5);
          site[0][isite] -= lxy * (t = floor(MATMUL(1,hinv,site,isite) + 0.5));
          site[1][isite] -= ly  * t;
          site[0][isite] -= lxz * (t = floor(MATMUL(2,hinv,site,isite) + 0.5));
          site[1][isite] -= lyz * t;
          site[2][isite] -= lz  * t;
      }
   xfree(ssite);
}
/******************************************************************************
 *  newton   Apply newton's equation to calculate the acceleration of a       *
 *  number of molecules given the force.                                      *
 ******************************************************************************/
void newton(force, acc, mass, nmols)
vec_mp		force,		/* Centre of mass forces [nmols][3]      (in) */
		acc;		/* Accelerations [nmols][3]             (out) */
double		mass;		/* Molecular mass                        (in) */
int		nmols;		/* Number of molecules                   (in) */
{
   int	imol, i;
   double	rmass = 1.0/mass;
   for(i=0; i < 3; i++)
      for(imol = 0; imol < nmols; imol++)
         acc[imol][i] = force[imol][i] * rmass;
}
/******************************************************************************
 * euler  Use the Euler equations and the second-order quaternion method to   *
 * calculate the second derivatives of the molecular quaternions from the     *
 * torques etc.  Test for zero moments of inertia is to handle case of linear *
 * molecule.                                                                  *
 ******************************************************************************/
void euler(torque, quat, qdot, qddot, inertia, nmols)
vec_mp		torque;		/* Space frame torques [nmols][3]        (in) */
quat_mp		quat, qdot,     /* Quaternions for this species and d/dt (in) */
		qddot;		/* Quaternion second derivatives        (out) */
vec_mt		inertia;	/* Principal moments of inertia          (in) */
int		nmols;		/* Number of molecules                   (in) */
{
   /* The following two quantities, though vectors, are stored in the last 3  *
    * components of a quaternion array to allow easy application of the       *
    * quaternion multiplication in the equations of motion.                   */
   quat_mp	omega,		/* Principal frame angular velocities         */
   		ang_acc=qalloc(nmols);	/* Principal frame ang accelerations  */
   register int	imol;
   int		i, j, k;
   register double	Iir, Ijk;	/* Temporaries for moments of inertia */

   if(quat == NULL) return;  /* Molecule is point atom or ion - no action     */

   omega = qddot;	/* Use the qddot array as workspace for angular vels  */

   q_conj_mul(quat, qdot, omega, nmols);
   vscale(4 * nmols, 2.0, omega[0], 1);

   for(imol = 0; imol < nmols; imol++)
      ang_acc[imol][0] = -2.0 * vdot(4,qdot[imol],1,qdot[imol],1);

   for(i=1, j=2, k=3; i<4; i++, j=i%3+1, k=j%3+1)
      if(inertia[i-1] != 0.0)
      {
         Iir = 1.0/inertia[i-1];  Ijk = inertia[j-1] - inertia[k-1];
         for(imol = 0; imol < nmols; imol++)
            ang_acc[imol][i] = Iir * 
                      (torque[imol][i-1] + Ijk * omega[imol][j]*omega[imol][k]);
      }
      else
         for(imol = 0; imol < nmols; imol++)
            ang_acc[imol][i] = 0.0;

   omega = NULL;		/* Omega not needed any more - re use as qddot*/

   q_mul(quat, ang_acc, qddot, nmols);
   vscale(4 * nmols, 0.5, qddot[0], 1);

   xfree(ang_acc);
}
/******************************************************************************
 *  Parinello   Calculate the correction to the scaled centre of mass         *
 *  accelerations in the Parinello and Rahman zero-stress method.             *
 *  Parinello M. and Rahman A. J. Appl. Phys. 52(12), 7182-7190 (1981)        *
 ******************************************************************************/
void parinello(h, h_dot, vel, acc, acc_out, nmols)
mat_mt	h,			/* P and R's unit cell matrix            (in) */
	h_dot;			/* Derivative of h matrix                (in) */
vec_mp	vel,			/* Centre of mass scaled velocities      (in) */
	acc, acc_out;		/* C of M accelerations              (in/out) */
int	nmols;			/* Size of vel and acc/ number molecules (in) */
{
   mat_mt	h_tr,		/* Transpose of h			      */
   		h_tr_dot,	/* Transpose of h_dot			      */
   		h_tmp_1,	/* Store for intermediate terms		      */
   		h_tmp_2,	/* Store for intermediate terms		      */
   		G,		/* h_tr * h	(metric tensor)		      */
   		G_inv,		/* Inverse of G				      */
   		G_dot,		/* Derivative of G			      */
   		G_i_d;		/* G_inv * G_dot			      */
   vec_mp	acc_corr=ralloc(nmols);	/* Correction term to accelerations   */
   int		i, imol;	/* Counters				      */

   transpose(h,h_tr);
   mat_mul(h_tr,h,G);				/* We now have the G matrix   */
   invert(G, G_inv);				/* G (-1) done                */
   transpose(h_dot, h_tr_dot);
   mat_mul(h_tr_dot, h, h_tmp_1);
   mat_mul(h_tr, h_dot, h_tmp_2);
   mat_add(h_tmp_1, h_tmp_2, G_dot);		/* G dot now complete         */
   mat_mul(G_inv, G_dot, G_i_d);		/* G_inv * G_dot              */

   mat_vec_mul(G_i_d, vel, acc_corr, nmols);    /* Calculate correction term  */

   for(i = 0; i < 3; i++)                       /* Add correction term        */
      for(imol = 0; imol < nmols; imol++)       /* to accelerations           */
         acc_out[imol][i] = acc[imol][i] - acc_corr[imol][i];

   xfree(acc_corr);
}
/******************************************************************************
 *  Trans_ke  calculate and return the translational kinetic energy           *
 ******************************************************************************/
double	trans_ke(h, vel_s, mass, nmols)
mat_mt	h;			/* Unit cell matrix			 (in) */
vec_mt	vel_s[];		/* Scaled c of m velocities		 (in) */
double	mass;			/* Mass of a molecule of this species	 (in) */
int	nmols;			/* Number of molecules			 (in) */
{
   double	ke;
   vec_mp	vel = ralloc(nmols);	/* Unscaled (real) velocities         */
   
   mat_vec_mul(h, vel_s, vel, nmols);   /* Calculate unscaled velocities      */

   ke = vdot(3*nmols, vel[0], 1, vel[0], 1);

   xfree(vel);
   return(0.5 * mass * ke);
}
   
/******************************************************************************
 *  rot_ke  calculate and return the rotational kinetic energy                *
 ******************************************************************************/
double	rot_ke(quat, qdot, inertia, nmols)
quat_mt	quat[],			/* Molecular quaternions		 (in) */
	qdot[];			/* Quaternion derivatives		 (in) */
vec_mt	inertia;		/* Principal moments of inertia		 (in) */
int	nmols;			/* Number of molecules			 (in) */
{
   double	ke = 0.0;
   quat_mp	omega_p = qalloc(nmols);   /* Principal angular velocities    */
   int		i;
   
   q_conj_mul(quat, qdot, omega_p, nmols); /* Calculate angular velocities    */
   vscale(4 * nmols, 2.0, omega_p[0], 1);  /* omega = 2*q~*qdot               */
   for(i = 0; i < 3; i++)
      ke += inertia[i] * vdot(nmols, omega_p[0]+i+1, 4, omega_p[0]+i+1, 4);

   xfree(omega_p);
   return(0.5 * ke);
}
/******************************************************************************
 * energy_dyad.  Calculate the dyadic sum m V V (dyad over V) for zero stress *
 ******************************************************************************/
void energy_dyad(ke_dyad, h, vels, mass, nmols)
mat_mt	ke_dyad,			/* Dyad is accumulated here  (in/out) */
	h;				/* Unit cell matrix		(in)  */
vec_mp	vels;				/* Scaled velocities		(in)  */
double	mass;				/* Mass of particles		(in)  */
int	nmols;				/* Number of molecules		(in)  */
{
   int		i, j;				/* Counters		      */
   vec_mp	vel = ralloc(nmols);		/* Real velocities	      */

   mat_vec_mul(h, vels, vel, nmols);	/* Calculate unscaled velocities      */

   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
      {
         ke_dyad[i][j] += mass * vdot(nmols, vel[0]+i, 3, vel[0]+j, 3);
      }     

   xfree(vel);
}
/******************************************************************************
 * Rahman   Calculate the unit cell matrix accelerations                      *
 ******************************************************************************/
void rahman(stress_vir, h, hddot, ke_dyad, press, W, mask)
mat_mt	stress_vir,			/* Stress virial		      */
	h,				/* Unit cell matrix		      */
	hddot,				/* Unit cell accelerations            */
	ke_dyad;			/* Translational kinetic energy dyad  */
double	press,				/* Externally applied pressure	      */
	W;				/* Piston mass parameter	      */
int	mask;				/* Mask constrained el's of h matrix  */
{
   double	vol = det(h);		/* Unit cell volume		      */
   mat_mt	stress,			/* Stress tensor		      */
   		h_tr,			/* Transpose of h		      */
   		h_tr_inv,		/* Inverse of transpose of h	      */
   		sigma;			/* P & R sigma matrix 		      */
   int		i, j;			/* Counters			      */

   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
         stress[i][j] = (ke_dyad[i][j] + stress_vir[i][j]) / vol;

   for(i = 0; i < 3; i++)
      stress[i][i] -= press;	/* Subtract applied pressure from diagonal    */

   transpose(h, h_tr);          /* Calculate sigma = vol*h transpose inverse  */
   invert(h_tr, h_tr_inv);
   mat_sca_mul(vol, h_tr_inv, sigma);

   mat_mul(stress, sigma, hddot);	/* Calculate unit cell accelerations  */
   mat_sca_mul(1.0/W, hddot, hddot);

   /* 
    * Zero unwanted degrees of freedom. Refson PhD Thesis (1986)
    */   
   for(i = 0; i < 9; i++)
   {
      if( mask & 1 )
	 hddot[0][i] = 0.0;		/* Access as [9] rather than [3][3]   */
      mask >>= 1;
   }
}
