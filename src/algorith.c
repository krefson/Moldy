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
 *       $Log$
 */
#ifndef lint
static char *RCSid = "$Header$";
#endif
/*========================== Library include files ===========================*/
#include 	<math.h>
#include 	"string.h"
/*========================== program include files ===========================*/
#include 	"defs.h"
#include 	"messages.h"
/*========================== Library declarations ============================*/
void	cfree();
/*========================== External function declarations ==================*/
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
void	message();			/* Error message and exit handler     */
/*========================== Macros ==========================================*/
#define	veccpy(v1,v2,n) (void)memcpy((char*)(v1),(char*)(v2),(n)*sizeof(vec_t))
/*============================================================================*/
/******************************************************************************
 *  rotate        Perform the rotation described by the quaternions in the    *
 *  second parameter on each of the co-ordinates in the first, putting the    *
 *  result in the third.  (Pawley,  Mol. Phys. 43, 1321-1330 (1981))          *
 *  NB this is different to Evans' formulation.                               *
 *  Apply each rotation to nvec/nquat vectors.                                *
 ******************************************************************************/
void rotate(r_in, r_out, nvec, quat, nquat, inv_mat)
vec_p		r_in,		/* Co-ordinates to be rotated [n][3] (in)     */
		r_out;		/* Resulting co-ordinates [n][3]    (out)     */
quat_p		quat;		/* Quaternions for the rotation.     (in)     */
int		nvec,		/* Number of co-ordinates.           (in)     */
		nquat;		/* Number of quaternions             (in)     */
invrot		inv_mat;	/* Flag to do inverse rotations      (in)     */
{
   mat_t	rot;
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
vec_t	x[];
vec_t	meansq;
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
   register	double s1=0,s2=0;	/* Accumulators for sums	      */
   while(n-- > 0)
   {
      s1 += (*v1 - *v2) * (*v1 - *v2);
      s2 += *v1 * *v1;
      v1++; v2++;
   }
   return(sqrt(s1 / s2));
}
/******************************************************************************
 *  molecule_force     Calculate the centre of mass forces on a number of     *
 *  molecules given the site forces and the site co-ordinates.                *
 ******************************************************************************/
void mol_force(site_force, force, nsites, nmols)
vec_p		site_force,	/* Site forces [nsites*nmols][3]        (in)  */
		force;		/* Centre of mass forces [nmols][3]    (out)  */
int		nsites,		/* Number of sites on one molecule      (in)  */
		nmols;		/* Number of molecules                  (in)  */
{
   int	i,  imol;

   for(imol = 0; imol < nmols; imol++)
      for(i = 0; i < 3; i++)
         force[imol][i] = sum(nsites, site_force[imol*nsites]+i, 3);
}
/******************************************************************************
 *  molecule_torque    Calculate the torque on a number of identical          *
 *  molecules given the space frame site forces and co-ordinates.             *
 ******************************************************************************/
void mol_torque(site_force, site, torque, quat, nsites, nmols)
vec_p		site_force,	/* Principal frame site forces          (in)  */
		site,		/* Principal frame site co-ordinates    (in)  */
		torque;		/* Molecular torques [nmols][3]        (out)  */
quat_p		quat;		/* Molecular quaternions [nmols][4]     (in)  */
int		nsites,		/* Number of sites on one molecule      (in)  */
		nmols;		/* Number of molecules                  (in)  */
{
   vec_p	princ_force = ralloc(nsites);
   int	i, j, k, imol, isite;

   for(imol = 0; imol < nmols; imol++)
   {
      rotate(site_force+imol*nsites, princ_force, nsites, quat+imol, 1, inv);
      for(i = 0, j = 1, k = 2; i < 3; i++, j=(j+1)%3, k=(k+1)%3)
      {
         torque[imol][i] = 0.0;
	 for(isite = 0; isite < nsites; isite++)
	    torque[imol][i] += site[isite][j]*princ_force[isite][k]
	                      -site[isite][k]*princ_force[isite][j];
      }
   }
   cfree((char *)princ_force);
}
/******************************************************************************
 *  make_sites     Calculate the atomic site co-ordinates for nmols identical *
 *  molecules from the principal-frame sites, the quaternions and the centre  *
 *  of mass co-ordinates.  Called once for each molecular species.            *
 ******************************************************************************/
void make_sites(h, c_of_m_s , quat, p_f_sites, site, nmols, nsites)
mat_t		h;		/* Unit cell matrix h		     (in)     */
vec_p		c_of_m_s,	/* Centre of mass co-ords [nmols][3] (in)     */
		p_f_sites,	/* Principal-frame sites [nsites][3] (in)     */
		site;		/* Sites [nmols*nsites][3]          (out)     */
quat_p		quat;		/* Quaternions [nmols][4]            (in)     */
int		nmols,		/* Number of molecules                        */
		nsites;		/* Number of sites on each molecule           */
{
   int		imol, isite, i;	/* Counters				      */
   vec_t	c_of_m;		/* Unscaled centre of mass co-ordinates       */

   for(imol = 0; imol < nmols; imol++)
   {
      mat_vec_mul(h,c_of_m_s+imol,(vec_p)c_of_m, 1);/* Get real c-of-m co-ords*/
      if(quat)
         rotate(p_f_sites,site+imol*nsites,nsites,quat+imol,1,noinv);
      else
         veccpy(site+imol*nsites, p_f_sites, nsites);
      for(i = 0; i < 3; i++)
         for(isite = 0; isite < nsites; isite++)
            site[imol*nsites+isite][i] += c_of_m[i];
   }
}
/******************************************************************************
 *  newton   Apply newton's equation to calculate the acceleration of a       *
 *  number of molecules given the force.                                      *
 ******************************************************************************/
void newton(force, acc, mass, nmols)
vec_p		force,		/* Centre of mass forces [nmols][3]      (in) */
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
vec_p		torque;		/* Space frame torques [nmols][3]        (in) */
quat_p		quat, qdot,     /* Quaternions for this species and d/dt (in) */
		qddot;		/* Quaternion second derivatives        (out) */
vec_t		inertia;	/* Principal moments of inertia          (in) */
int		nmols;		/* Number of molecules                   (in) */
{
   /* The following two quantities, though vectors, are stored in the last 3  *
    * components of a quaternion array to allow easy application of the       *
    * quaternion multiplication in the equations of motion.                   */
   quat_p	omega,		/* Principal frame angular velocities         */
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

   cfree((char*)ang_acc);
}
/******************************************************************************
 *  Parinello   Calculate the correction to the scaled centre of mass         *
 *  accelerations in the Parinello and Rahman zero-stress method.             *
 *  Parinello M. and Rahman A. J. Appl. Phys. 52(12), 7182-7190 (1981)        *
 ******************************************************************************/
void parinello(h, h_dot, vel, acc, acc_out, nmols)
mat_t	h,			/* P and R's unit cell matrix            (in) */
	h_dot;			/* Derivative of h matrix                (in) */
vec_p	vel,			/* Centre of mass scaled velocities      (in) */
	acc, acc_out;		/* C of M accelerations              (in/out) */
int	nmols;			/* Size of vel and acc/ number molecules (in) */
{
   mat_t	h_tr,		/* Transpose of h			      */
   		h_tr_dot,	/* Transpose of h_dot			      */
   		h_tmp_1,	/* Store for intermediate terms		      */
   		h_tmp_2,	/* Store for intermediate terms		      */
   		G,		/* h_tr * h	(metric tensor)		      */
   		G_inv,		/* Inverse of G				      */
   		G_dot,		/* Derivative of G			      */
   		G_i_d;		/* G_inv * G_dot			      */
   vec_p	acc_corr=ralloc(nmols);	/* Correction term to accelerations   */
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
         acc_out[imol][i] = acc[imol][i] + acc_corr[imol][i];

   cfree((char*)acc_corr);
}
/******************************************************************************
 *  Trans_ke  calculate and return the translational kinetic energy           *
 ******************************************************************************/
double	trans_ke(h, vel_s, mass, nmols)
mat_t	h;			/* Unit cell matrix			 (in) */
vec_t	vel_s[];		/* Scaled c of m velocities		 (in) */
double	mass;			/* Mass of a molecule of this species	 (in) */
int	nmols;			/* Number of molecules			 (in) */
{
   double	ke;
   vec_p	vel = ralloc(nmols);	/* Unscaled (real) velocities         */
   
   mat_vec_mul(h, vel_s, vel, nmols);   /* Calculate unscaled velocities      */

   ke = vdot(3*nmols, vel[0], 1, vel[0], 1);

   cfree((char*)vel);
   return(0.5 * mass * ke);
}
   
/******************************************************************************
 *  rot_ke  calculate and return the rotational kinetic energy                *
 ******************************************************************************/
double	rot_ke(quat, qdot, inertia, nmols)
quat_t	quat[],			/* Molecular quaternions		 (in) */
	qdot[];			/* Quaternion derivatives		 (in) */
vec_t	inertia;		/* Principal moments of inertia		 (in) */
int	nmols;			/* Number of molecules			 (in) */
{
   double	ke = 0.0;
   quat_p	omega_p = qalloc(nmols);   /* Principal angular velocities    */
   int		i;
   
   q_conj_mul(quat, qdot, omega_p, nmols); /* Calculate angular velocities    */
   vscale(4 * nmols, 2.0, omega_p[0], 1);  /* omega = 2*q~*qdot               */
   for(i = 0; i < 3; i++)
      ke += inertia[i] * vdot(nmols, omega_p[0]+i+1, 4, omega_p[0]+i+1, 4);

   cfree((char*)omega_p);
   return(0.5 * ke);
}
/******************************************************************************
 * energy_dyad.  Calculate the dyadic sum m V V (dyad over V) for zero stress *
 ******************************************************************************/
void energy_dyad(ke_dyad, h, vels, mass, nmols)
mat_t	ke_dyad,			/* Dyad is accumulated here  (in/out) */
	h;				/* Unit cell matrix		(in)  */
vec_p	vels;				/* Scaled velocities		(in)  */
double	mass;				/* Mass of particles		(in)  */
int	nmols;				/* Number of molecules		(in)  */
{
   int		i, j;				/* Counters		      */
   vec_p	vel = ralloc(nmols);		/* Real velocities	      */

   mat_vec_mul(h, vels, vel, nmols);	/* Calculate unscaled velocities      */

   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
      {
         ke_dyad[i][j] += mass * vdot(nmols, vel[0]+i, 3, vel[0]+j, 3);
      }     

   cfree((char*)vel);
}
/******************************************************************************
 * Rahman   Calculate the unit cell matrix accelerations                      *
 ******************************************************************************/
void rahman(stress_vir, h, hddot, ke_dyad, press, W)
mat_t	stress_vir,			/* Stress virial		      */
	h,				/* Unit cell matrix		      */
	hddot,				/* Unit cell accelerations            */
	ke_dyad;			/* Translational kinetic energy dyad  */
double	press,				/* Externally applied pressure	      */
	W;				/* Piston mass parameter	      */
{
   double	vol = det(h);		/* Unit cell volume		      */
   mat_t	stress,			/* Stress tensor		      */
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

   hddot[1][0] = hddot[2][0] = hddot[2][1] = 0.0;
   	/* Zero unwanted degrees of freedom. Refson PhD Thesis (1986)         */
}
