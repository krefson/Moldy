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
 * parinello()          Calculate P&R c-of-m acceleration term                * 
 * rahman()		Calculate unit cell matrix accelerations from stress  *
 * trans_ke()		Return translational kinetic energy.		      *
 * rot_ke()		Return rotational kinetic energy		      *
 * energy_dyad()	Calculate kinetic energy part of stress tensor	      *
 * gaussiant()          Return alpha in Gaussian thermostat(trans)            *
 * gaussianr1()         Return alpha1 in Gaussian thermostat(rot)             *
 * gaussianr2()         Return alpha2 in Gaussian thermostat(rot)             *
 * hoover_tr()          Calculate correction to the forces(accelerations)     *
 * hoover_rot()         due to thermostat                                     *
 ******************************************************************************
 *      Revision Log
 *       $Log: algorith.c,v $
 *       Revision 2.14.2.2  2000/10/20 13:59:31  keith
 *       Incorporated new neightbour list stuff into accel.c.
 *       Removed old "poteval" from accel.c.  Now use one in
 *       force.
 *       Other errors corrected and declarations tidied somewhat.
 *
 *       Revision 2.14.2.1  2000/10/20 11:48:51  keith
 *       Incorporated new neighbour list indexing algorithm from 2.16
 *
 *       Revision 2.14  2000/05/23 15:23:07  keith
 *       First attempt at a thermostatted version of the Leapfrog code
 *       using either a Nose or a Nose-Poincare thermostat
 *
 *       Revision 2.13  2000/04/27 17:57:05  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.12  2000/04/26 16:01:01  keith
 *       Dullweber, Leimkuhler and McLachlan rotational leapfrog version.
 *
 *       Revision 2.11  1999/09/14 13:31:31  keith
 *       Experimental version implementing true const P (strain-mask=512).
 *
 *       Revision 2.10  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 *
 *       Revision 2.9  1995/12/04 11:45:49  keith
 *       Nose-Hoover and Gaussian (Hoover constrained) thermostats added.
 *       Thanks to V. Murashov.
 *
 * Revision 2.8  1994/07/07  16:52:14  keith
 * Performance optimization to mol_force.
 *
 * Revision 2.7  1994/06/08  13:22:31  keith
 * Null update for version compatibility
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
 * "dump".Changed size_t to own typedef size_mt == ulong.
 *
 * Revision 2.5  1994/01/18  13:32:05  keith
 * Null update for XDR portability release
 *
 * Revision 2.3  93/10/28  10:27:35  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
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
static char *RCSid = "$Header: /home/minphys2/keith/CVS/moldy/src/algorith.c,v 2.14.2.2 2000/10/20 13:59:31 keith Exp $";
#endif
/*========================== program include files ===========================*/
#include 	"defs.h"
#include 	"messages.h"
/*========================== Library include files ===========================*/
#include 	<math.h>
#include 	"string.h"
#include	"stddef.h"
/*========================== External function declarations ==================*/
gptr            *talloc(int n, size_mt size, int line, char *file);	       /* Interface to memory allocator       */
void            tfree(gptr *p);	       /* Free allocated memory	      	      */
void	mat_vec_mul(mat_mt m, vec_mp in_vec, vec_mp out_vec, int number);			/* 3 x 3 Matrix by Vector multiplier  */
void	mat_mul(mat_mt a, mat_mt b, mat_mt c);	          	/* 3 x 3 matrix multiplier	      */
void	mat_add(mat_mt a, mat_mt b, mat_mt c);			/* Add 2 3x3 matrices                 */
void	mat_sca_mul(register real s, mat_mt a, mat_mt b);			/* Multiply 3x3 matrix by scalar      */
void	transpose(mat_mt a, mat_mt b);			/* transpose a 3x3 matrix	      */
void	invert(mat_mt a, mat_mt b);			/* invert a 3x3 matrix		      */
double	det(mat_mt a);				/* Determinant of 3x3 matrix	      */
void	q_to_rot(real *quat, mat_mt rot);			/* Make rotation matrix from quat'n   */
void	q_mul(quat_mp p, quat_mp q, quat_mp r, int n);
void	q_conj_mul(quat_mp p, quat_mp q, quat_mp r, int n);
double	vdot(int n, real *x, int ix, real *y, int iy);				/* Vector dot product		      */
double	sum(register int n, register double *x, register int ix);				/* Vector sum			      */
void	vscale(register int n, register double s, register real *x, register int ix);
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
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
static
void rotate(vec_mp r_in, vec_mp r_out, int nvec, quat_mp quat, int nquat, invrot inv_mat)
      		     		/* Co-ordinates to be rotated [n][3] (in)     */
		      		/* Resulting co-ordinates [n][3]    (out)     */
       		     		/* Quaternions for the rotation.     (in)     */
   		     		/* Number of co-ordinates.           (in)     */
		      		/* Number of quaternions             (in)     */
      		        	/* Flag to do inverse rotations      (in)     */
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
void	mean_square(vec_mt (*x), real *meansq, int nmols)
{
   int		i;
   for(i = 0; i < 3; i++)
      meansq[i] = vdot(nmols, x[0]+i, 3, x[0]+i, 3) / nmols;
}
/******************************************************************************
 *   vec_dist      Return the normalised distance between 2 vectors           *
 ******************************************************************************/
double	vec_dist(real *v1, real *v2, int n)
    	         		/* Input vectors			      */
   	  			/* Length ie v1[n], v2[n]		      */
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
void mol_force(real **site_force, vec_mp force, int nsites, int nmols)
    		             	/* Site forces [nsites*nmols][3]        (in)  */
      		      		/* Centre of mass forces [nmols][3]    (out)  */
   		       		/* Number of sites on one molecule      (in)  */
		      		/* Number of molecules                  (in)  */
{
   int	i,  imol, isite;
   double	f;

   for(imol = 0; imol < nmols; imol++)
      for(i = 0; i < 3; i++)
      {
	 f = 0;
	 for(isite=0; isite < nsites; isite++)
	    f += site_force[i][isite+imol*nsites];
	 force[imol][i] = f;
      }
}
/******************************************************************************
 *  molecule_torque    Calculate the torque on a number of identical          *
 *  molecules given the space frame site forces and co-ordinates.             *
 ******************************************************************************/
void mol_torque(real **site_force, vec_mp site, vec_mp torque, quat_mp quat, int nsites, int nmols)
    		             	/* Principal frame site forces          (in)  */
      		     		/* Principal frame site co-ordinates    (in)  */
		       		/* Molecular torques [nmols][3]        (out)  */
       		     		/* Molecular quaternions [nmols][4]     (in)  */
   		       		/* Number of sites on one molecule      (in)  */
		      		/* Number of molecules                  (in)  */
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
/******************************************************************************
 *  make_sites     Calculate the atomic site co-ordinates for nmols identical *
 *  molecules from the principal-frame sites, the quaternions and the centre  *
 *  of mass co-ordinates.  Called once for each molecular species.            *
 ******************************************************************************/
void make_sites(mat_mt h, vec_mp c_of_m_s, quat_mp quat, vec_mp p_f_sites, real **site, int nmols, int nsites, int molflag)
      		  		/* Unit cell matrix h		     (in)     */
      		         	/* Centre of mass co-ords [nmols][3] (in)     */
		          	/* Principal-frame sites [nsites][3] (in)     */
    		       		/* Sites [nmols*nsites][3]          (out)     */
       		     		/* Quaternions [nmols][4]            (in)     */
   		      		/* Number of molecules                        */
		       		/* Number of sites on each molecule           */
   		        	/* Whether to apply pbc to sites or cofm      */
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

   if( molflag!=MOLPBC ) /* Apply pbc's to put sites into cell */
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
void newton(vec_mp force, vec_mp acc, double mass, int nmols)
      		      		/* Centre of mass forces [nmols][3]      (in) */
		    		/* Accelerations [nmols][3]             (out) */
      		     		/* Molecular mass                        (in) */
   		      		/* Number of molecules                   (in) */
{
   int	imol, i;
   double	rmass = 1.0/mass;
   for(i=0; i < 3; i++)
      for(imol = 0; imol < nmols; imol++)
      {
         acc[imol][i] = force[imol][i] * rmass;
#ifdef DEBUG2
   printf("Newton accelerations for %4i mol %4i orient %8.4f \n",imol,
           i, acc[imol][i]);
#endif
      }
}
/******************************************************************************
 *  Parinello   Calculate the correction to the scaled centre of mass         *
 *  accelerations in the Parinello and Rahman zero-stress method.             *
 *  Parinello M. and Rahman A. J. Appl. Phys. 52(12), 7182-7190 (1981)        *
 ******************************************************************************/
void parinello(mat_mt h, mat_mt h_dot, vec_mp vel, vec_mp acc, vec_mp acc_out, int nmols)
      	  			/* P and R's unit cell matrix            (in) */
	      			/* Derivative of h matrix                (in) */
      	    			/* Centre of mass scaled velocities      (in) */
	             		/* C of M accelerations              (in/out) */
   	      			/* Size of vel and acc/ number molecules (in) */
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
double	trans_ke(mat_mt h, vec_mt (*vel_s), real s, double mass, int nmols)
      	  			/* Unit cell matrix			 (in) */
      	        		/* Scaled c of m velocities		 (in) */
      	     			/* Mass of a molecule of this species	 (in) */
   	      			/* Number of molecules			 (in) */
{
   double	ke;
   vec_mp	vel = ralloc(nmols);	/* Unscaled (real) velocities         */
   
   mat_vec_mul(h, vel_s, vel, nmols);   /* Calculate unscaled velocities      */

   ke = vdot(3*nmols, vel[0], 1, vel[0], 1);

   xfree(vel);
   return(0.5 * mass * ke / SQR(s));
}
   
/******************************************************************************
 *  rot_ke  calculate and return the rotational kinetic energy                *
 ******************************************************************************/
double	rot_ke(quat_mt (*omega_p), real *inertia, int nmols)
       	          		/* Principal angular velocities    */
      	        		/* Principal moments of inertia		 (in) */
   	      			/* Number of molecules			 (in) */
{
   double	ke = 0.0;
   int		i;
   
   for(i = 0; i < 3; i++)
      ke += inertia[i] * vdot(nmols, omega_p[0]+i+1, 4, omega_p[0]+i+1, 4);

   return(0.5 * ke);
}
/******************************************************************************
 * energy_dyad.  Calculate the dyadic sum m V V (dyad over V) for zero stress *
 ******************************************************************************/
void energy_dyad(mat_mt ke_dyad, mat_mt h, double s, vec_mp vels, double mass, int nmols)
      	        			/* Dyad is accumulated here  (in/out) */
	  				/* Unit cell matrix		(in)  */
      	     				/* Scaled velocities		(in)  */
      	     				/* Mass of particles		(in)  */
   	      				/* Number of molecules		(in)  */
{
   int		i, j;				/* Counters		      */
   vec_mp	vel = ralloc(nmols);		/* Real velocities	      */

   mat_vec_mul(h, vels, vel, nmols);	/* Calculate unscaled velocities      */

   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
      {
         ke_dyad[i][j] += mass * vdot(nmols, vel[0]+i, 3, vel[0]+j, 3)/SQR(s);
      }     

   xfree(vel);
}
double trace(mat_mt mat)
{
  return(mat[0][0] + mat[1][1] + mat[2][2]);
}
/******************************************************************************
 * Rahman   Calculate the unit cell matrix accelerations                      *
 ******************************************************************************/
void rahman(mat_mt stress_vir, mat_mt h, mat_mt hddot, mat_mt ke_dyad, double press, double W, int mask)
      	           			/* Stress virial		      */
	  				/* Unit cell matrix		      */
	      				/* Unit cell accelerations            */
	        			/* Translational kinetic energy dyad  */
      	      				/* Externally applied pressure	      */
	  				/* Piston mass parameter	      */
   	     				/* Mask constrained el's of h matrix  */
{
   double	vol = det(h);		/* Unit cell volume		      */
   double       cpscale;                /* Scale factor for uniform scale case*/
   mat_mt	stress,			/* Stress tensor		      */
   		h_tr,			/* Transpose of h		      */
   		h_tr_inv,		/* Inverse of transpose of h	      */
                htrh,                   /* Product  h'h                       */
   		sigma;			/* P & R sigma matrix 		      */
   int		i, j;			/* Counters			      */

   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
         stress[i][j] = (ke_dyad[i][j] + stress_vir[i][j]) / vol;

   for(i = 0; i < 3; i++)
      stress[i][i] -= press;	/* Subtract applied pressure from diagonal    */

   transpose(h, h_tr);          /* Calculate sigma = vol*h transpose inverse  */
   if ( mask >> 9 & 1 )                         /* Uniform dilation case */
   {
     mat_mul(h_tr, h, htrh);
     cpscale = 3.0*vol*trace(stress)/(W*trace(htrh));
     mat_sca_mul(cpscale,h, hddot);
   }
   else                                         /* Full parrinello-rahman */
   {
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
	mask  >>= 1;
     }
   }

}
/******************************************************************************
 * Hoover_tr() function corrects forces to realize thermostat                 *
 * main formula: mass * accel = force - alpha * mass * vel                    *
 * Function is added by VVMurashov on 22/10/95                                *
 ******************************************************************************/
void hoover_tr(double alpha, vec_mp accel_in, vec_mp accel_out, vec_mp vel, int nmols)
                                /* Thermostat  multiplier                     */
      	    			/* Centre of mass scaled velocities      (in) */
                                /* Array of forces                   (in/out) */
                  
   	      			/* Size of vel and force/number molecules(in) */
{
   int		i, j;   	/* Counters				      */

   for(i = 0; i < 3; i++)                       /* Add correction term        */
      for(j = 0; j < nmols; j++)                /* to accelerations           */
         accel_out[j][i] = accel_in[j][i] - alpha*vel[j][i];
}
/******************************************************************************
 * Hoover_rot() function corrects forces to realize thermostat                *
 * main formula: inertia * accel = torque - alpha * inertia * omega           *
 * Function is added by VVMurashov on 22/10/95                                *
 ******************************************************************************/
void hoover_rot(double alpha, real *inertia, vec_mp force_in, vec_mp force_out, quat_mp omega, int nmols)
                                /* Nose-Hoover multiplier                     */
                                /* Inertia vector                             */
       	       		        /* Angular velocities                    (in) */
                                /* Array of forces                   (in/out) */
                   
   	       			/* Size of vel and force/number molecules(in) */
{
   int		i, imols;   	/* Counters				      */
   for(i = 0; i < 3; i++)                       /* Add correction term        */
      for(imols = 0; imols < nmols; imols++)    /* to accelerations           */
         force_out[imols][i] = force_in[imols][i] - alpha*inertia[i]*
                               omega[imols][i+1];
}
/*****************************************************************************
 * Gaussiant() function calculates alpha1 or alpha2 for the Gaussian         * 
 * thermostat. Principle formula: alpha = alpha1 / alpha2                    *
 * alpha1 = (SUM force * vel), and alpha2 = (SUM mass * vel^2)               *
 * Function is added by VVMurashov on 3/11/95                                *
 *****************************************************************************/
double gaussiant(vec_mp vec1, vec_mp vec2, int nmols)
                                /* First array                           (in) */
                                /* Second array                          (in) */
   	        		/* Size of arrays [nmols][n]             (in) */
{
   double   alpha = 0.0;        /* Sum of products                      (out) */
   int		i;      	/* Counters				      */
   for(i = 0; i < 3; i++)
      alpha += vdot(nmols, vec1[0]+i, 3, vec2[0]+i, 3);
   return(alpha);
}
/*****************************************************************************
 * Gaussianr1() function is the same as above, but allows for dealing with   *
 * quaternion array. Function is added by VVMurashov on 3/11/95              *
 *****************************************************************************/
double gaussianr1(vec_mp vec1, quat_mp vec2, int nmols)
                                /* First array                           (in) */
                                /* Second array                          (in) */
   	        		/* Size of arrays [nmols][n]             (in) */
{
   double   alpha = 0.0;        /* Sum of products                      (out) */
   int		i;      	/* Counters				      */
   for(i = 0; i < 3; i++)
      alpha += vdot(nmols, vec1[0]+i, 3, vec2[0]+i+1, 4);
   return(alpha);
}
/*****************************************************************************
 * Gaussianr2() function calculates alpha2 for the Gaussian thermostat       *
 * Principle formula: alpha = alpha1/alpha2                                  *
 * where alpha1 = (SUM torque * omega), and alpha2 = (SUM inertia * omega^2) *
 * Function is added by VVMurashov on 3/11/95                                *
 *****************************************************************************/
double gaussianr2(quat_mp omega, real *inertia, int nmols)
                                /* Inertia vector                             */
       	       		        /* Angular velocities                    (in) */
   	       			/* Size of vel and force/number molecules(in) */
{
   double   alpha = 0.0;        /* Gaussian multiplier                  (out) */
   int		i;      	/* Counters				      */
   for(i = 0; i < 3; i++)
      alpha += inertia[i] * vdot(nmols, omega[0]+i+1, 4, omega[0]+i+1, 4);
   return(alpha);
}
