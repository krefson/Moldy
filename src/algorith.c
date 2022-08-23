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
 * mol_force()		Calculate molecular centre of mass forces	      *
 * mol_torque()		Calculate molecular torques			      *
 * make_sites()		Generate atomic site co-ordinates from c-of-m etc     *
 * parinello()          Calculate P&R c-of-m acceleration term                * 
 * rahman()		Calculate unit cell matrix accelerations from stress  *
 * trans_ke()		Return translational kinetic energy.		      *
 * rot_ke()		Return rotational kinetic energy		      *
 * energy_dyad()	Calculate kinetic energy part of stress tensor	      *
 ******************************************************************************
 */
/*========================== program include files ===========================*/
#include 	"defs.h"
#include 	"messages.h"
/*========================== Library include files ===========================*/
#include 	<math.h>
/*========================== External function declarations ==================*/
gptr            *talloc(int n, size_mt size, int line, char *file);
				       /* Interface to memory allocator       */
void    tfree(gptr *p);		/* Free allocated memory	      	      */
void	mat_vec_mul(mat_mt m, vec_mp in_vec, vec_mp out_vec, int number);
					/* 3 x 3 Matrix by Vector multiplier  */
void	mat_mul(mat_mt a, mat_mt b, mat_mt c); /* 3 x 3 matrix multiplier     */
void	mat_add(mat_mt a, mat_mt b, mat_mt c); /* Add 2 3x3 matrices          */
void	mat_sca_mul(register real s, mat_mt a, mat_mt b); 
					/* Multiply 3x3 matrix by scalar      */
void	transpose(mat_mt a, mat_mt b); /* transpose a 3x3 matrix	      */
void	invert(mat_mt a, mat_mt b); /* invert a 3x3 matrix		      */
double	det(mat_mt a);			/* Determinant of 3x3 matrix	      */
void	q_to_rot(real *quat, mat_mt rot); /* Make rotation matrix from quat'n */
void	q_mul(quat_mp p, quat_mp q, quat_mp r, int n);
void	q_conj_mul(quat_mp p, quat_mp q, quat_mp r, int n);
double	vdot(int n, real *x, int ix, real *y, int iy); /* Vector dot product */
double	sum(register int n, register double *x, register int ix); /* Vector sum */
void	vscale(register int n, register double s, register real *x, register int ix);
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
/*========================== Macros ==========================================*/
#define MATMUL(i, m, r, o) (m[i][0]*r[0][o] + m[i][1]*r[1][o] + m[i][2]*r[2][o])
#ifndef INERTIA_MIN
#define INERTIA_MIN	1.0e-14		/* Tolerance for zero mom of I	      */
#endif
/*============================================================================*/
/******************************************************************************
 *  rotate        Perform the rotation described by the quaternions in the    *
 *  third parameter on each of the co-ordinates in the first, putting the     *
 *  result in the second.  (Pawley,  Mol. Phys. 43, 1321-1330 (1981))         *
 *  NB this is different to Evans' formulation.                               *
 *  Apply each rotation to nvec/nquat vectors.                                *
 ******************************************************************************/
static
void rotate(vec_mp r_in,        /* Co-ordinates to be rotated [n][3] (in)     */ 
	    vec_mp r_out,       /* Resulting co-ordinates [n][3]    (out)     */ 
	    int nvec, 	        /* Number of co-ordinates.           (in)     */
	    quat_mp quat,       /* Quaternions for the rotation.     (in)     */
	    int nquat, 	        /* Number of quaternions             (in)     */
	    invrot inv_mat)     /* Flag to do inverse rotations      (in)     */
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
 *  molecule_force     Calculate the centre of mass forces on a number of     *
 *  molecules given the site forces and the site co-ordinates.                *
 ******************************************************************************/
void mol_force(real **site_force, /* Site forces [nsites*nmols][3]       (in) */
	       vec_mp force, 	  /* Centre of mass forces [nmols][3]   (out) */
	       int nsites, 	  /* Number of sites on one molecule     (in) */
	       int nmols)	  /* Number of molecules                 (in) */
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
void mol_torque(real **site_force, /* Principal frame site forces        (in) */
		vec_mp site, 	   /* Principal frame site co-ordinates  (in) */
		vec_mp torque, 	   /* Molecular torques [nmols][3]      (out) */
		quat_mp quat, 	   /* Molecular quaternions [nmols][4]   (in) */
		int nsites, 	   /* Number of sites on one molecule    (in) */
		int nmols)	   /* Number of molecules                (in) */
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
void make_sites(mat_mt h,        /* Unit cell matrix h                   (in) */
		vec_mp c_of_m_s, /* Centre of mass co-ords [nmols][3]    (in) */
		quat_mp quat, 	 /* Quaternions [nmols][4]               (in) */
		vec_mp p_f_sites,/* Principal-frame sites [nsites][3]    (in) */
		real **site, 	 /* Sites [nmols*nsites][3]             (out) */
		int nmols, 	 /* Number of molecules                  (in) */
		int nsites, 	 /* Number of sites on each molecule     (in) */
		int molflag)	 /* Whether to apply pbc to sites or cofm(in) */
{
   int		imol, isite, i;	/* Counters				      */
   vec_mt	c_of_m;		/* Unscaled centre of mass co-ordinates       */
   vec_mt	*ssite = ralloc(nsites);
   register double	tx, ty, tz;
   mat_mt	hinv;
   double	lx   = h[0][0], lxy  = h[0][1], lxz  = h[0][2],
                lyx  = h[1][0], ly   = h[1][1], lyz  = h[1][2],
                lzx  = h[2][0], lzy  = h[2][1], lz   = h[2][2];
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
	 tx = floor(MATMUL(0,hinv,site,isite) + 0.5);
	 ty = floor(MATMUL(1,hinv,site,isite) + 0.5);
	 tz = floor(MATMUL(2,hinv,site,isite) + 0.5);
	 site[0][isite] -= lx  * tx + lxy * ty + lxz * tz;
	 site[1][isite] -= lyx * tx + ly  * ty + lyz * tz;
	 site[2][isite] -= lzx * tx + lzy * ty + lz  * tz;
      }
   xfree(ssite);
}
/******************************************************************************
 *  Trans_ke  calculate and return the translational kinetic energy           *
 ******************************************************************************/
double	trans_ke(mat_mt h,       /* Unit cell matrix                      (in) */
		 vec_mt (*mom),  /* Scaled c of m momenta                 (in) */ 
		 real s,         /* Thermostat, time scaling variable     (in) */
		 double mass,    /* Mass of a molecule of this species    (in) */ 
		 int nmols)      /* Number of molecules                   (in) */
{
   double	ke;
   vec_mp	real_mom = ralloc(nmols);	/* Unscaled (real) momenta     */
   mat_mt	hinv;
   
   invert(h, hinv);
   transpose(hinv,hinv);
   mat_vec_mul(hinv, mom, real_mom, nmols);	/* Calculate unscaled momenta  */

   ke = vdot(3*nmols, real_mom[0], 1, real_mom[0], 1);

   xfree(real_mom);
   return(ke / (2.0*mass*SQR(s)));
}
   
/******************************************************************************
 *  rot_ke  calculate and return the rotational kinetic energy                *
 ******************************************************************************/
double	rot_ke(quat_mt (*amom),   /* Principal angular momenta           (in) */ 
	       real s,            /* Thermostat, time scaling variable   (in) */
	       real *inertia,     /* Principal moments of inertia        (in) */
	       int nmols)	  /* Number of molecules                 (in) */
{
   double	ke = 0.0;
   int		i;
   
   for(i = 0; i < 3; i++)
     if( inertia[i] > INERTIA_MIN )
       ke += vdot(nmols, amom[0]+i+1, 4, amom[0]+i+1, 4)/inertia[i];
   
   return(0.5 * ke/ SQR(s));
}
/******************************************************************************
 * energy_dyad.  Calculate the dyadic sum m V V (dyad over V) for zero stress *
 ******************************************************************************/
void energy_dyad(mat_mt ke_dyad,        /* Dyad is accumulated here  (in/out) */ 
		 mat_mt h, 	        /* Unit cell matrix             (in)  */
		 double s, 	        /* Thermostat variable           (in) */
		 vec_mp mom, 	        /* Scaled momenta               (in)  */
		 double mass, 	        /* Mass of particles            (in)  */
		 int nmols)		/* Number of molecules          (in)  */
{
   int		i, j;				/* Counters		      */
   vec_mp	real_mom = ralloc(nmols);	/* Real momenta	      */
   mat_mt       hinv;
   
   invert(h, hinv);
   transpose(hinv,hinv);
   mat_vec_mul(hinv, mom, real_mom, nmols);	/* Calculate real momenta     */

   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
      {
         ke_dyad[i][j] += vdot(nmols, real_mom[0]+i, 3, real_mom[0]+j, 3)
	    /(mass * SQR(s));
      }     

   xfree(real_mom);
}

/* Calculate sigma = vol*h transpose inverse  */
void mk_sigma(mat_mt h, mat_mt sigma)
{
   int i,j,k,l,m,n;

   for(i=0, j=1, k=2; i<3; i++, j=(j+1)%3, k=(k+1)%3)
      for (l=0, m=1, n=2; l<3; l++, m=(m+1)%3, n=(n+1)%3)
	 sigma[i][l] = h[j][m]*h[k][n] - h[k][m]*h[j][n];
}

double ke_cell(mat_mt hmom, real w)
{
   double s = 0;
   int    i,j;
   
   for(i=0; i < 3; i++)
      for(j=0; j < 3; j++)
	 s += hmom[i][j]*hmom[i][j];

   return 0.5/w*s;
}
