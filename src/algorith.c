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
 *      Revision Log
 *       $Log: algorith.c,v $
 *       Revision 2.19  2001/05/24 16:26:43  keith
 *       Updated program to store and use angular momenta, not velocities.
 *        - added conversion routines for old restart files and dump files.
 *       Got rid of legacy 2.0 and lower restart file reading code.
 *
 *       Revision 2.18  2001/02/13 17:45:07  keith
 *       Added symplectic Parrinello-Rahman constant pressure mode.
 *
 *       Revision 2.17  2000/12/06 17:45:27  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.16  2000/11/06 16:02:05  keith
 *       First working version with a Nose-Poincare thermostat for rigid molecules.
 *
 *       System header updated to include H_0.
 *       Dump performs correct scaling  of angular velocities, but dumpext still
 *          needs to be updated to read this.
 *       XDR functions corrected to work with new structs.
 *       Parallel broadcast of config also updated.
 *       Some unneccessary functions and code deleted.
 *
 *       Revision 2.15  2000/10/20 15:15:46  keith
 *       Incorporated all mods and bugfixes from Beeman branch up to Rel. 2.16
 *
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
static char *RCSid = "$Header: /home/kr/CVS/moldy/src/algorith.c,v 2.19 2001/05/24 16:26:43 keith Exp $";
#endif
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
 *  second parameter on each of the co-ordinates in the first, putting the    *
 *  result in the third.  (Pawley,  Mol. Phys. 43, 1321-1330 (1981))          *
 *  NB this is different to Evans' formulation.                               *
 *  Apply each rotation to nvec/nquat vectors.                                *
 ******************************************************************************/
static
void rotate(vec_mp r_in,        /* Co-ordinates to be rotated [n][3] (in)     */ 
	    vec_mp r_out        /* Resulting co-ordinates [n][3]    (out)     */, 
	    int nvec, 	        /* Quaternions for the rotation.     (in)     */
	    quat_mp quat        /* Number of co-ordinates.           (in)     */, 
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
		quat_mp quat, 	 /* Principal-frame sites [nsites][3]    (in) */
		vec_mp p_f_sites,/* Sites [nmols*nsites][3]             (out) */ 
		real **site, 	 /* Quaternions [nmols][4]               (in) */
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
