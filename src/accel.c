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
 * Accel     This file contains the 'driver' function for a single timestep - *
 *           "do_step" and other functions which need access to the system    *
 *           and species structs, "rescale" and "distant_const". 	      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	accel.c,v $
 * Revision 2.4  93/12/16  18:14:11  keith
 * Fixed divide-by-zero bug when rescaling on atomic systems.
 * (Only showed up on FPE trapping architectures.)
 * 
 * Revision 2.3  93/10/28  10:27:32  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:47:36  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.8.1.26  93/03/15  14:39:59  keith
 * Added GPL copyleft notice to permit release and distribution.
 * N.B.  Previous versions were copyright (C) by author and 
 * only licensed by explicit permission.
 * 
 * Revision 1.8.1.25  93/03/09  14:15:09  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.8.1.24  92/10/28  14:09:42  keith
 * Changed "site_[tp]" typedefs to avoid name clash on HP.
 * 
 * Revision 1.8.1.23  92/10/01  18:08:09  keith
 * Added mol_radius().  Function used in force.c for cutoff calculation.
 * 
 * Revision 1.8.1.22  92/09/22  14:48:19  keith
 * Tidied up calls to improve "lint" rating.
 * 
 * Revision 1.8.1.21  92/06/26  17:01:18  keith
 * Got rid of assumption that memory returned by talloc() or
 * arralloc() is zeroed.  This enhances ANSI compatibility.
 * Removed memory zeroing from alloc.c() in consequence.
 * 
 * Revision 1.8.1.20  92/03/19  15:46:28  keith
 * Removed spurious varlaibe errptr.
 * 
 * Revision 1.8.1.19  92/03/11  12:56:16  keith
 * Changed "scale-separately" parameter to "scale options"
 * 
 * Revision 1.8.1.18  91/11/26  10:22:58  keith
 * Corrections to calculate framework pressure/stress correctly.
 * Corrected calculation of distant-stress term  (Only 1 1/v necessary).
 * 
 * Revision 1.8.1.17  91/10/17  14:22:38  keith
 * Fixed bug in set up pointers in "torque" array.
 * 
 * Revision 1.8.1.16  91/08/15  18:09:16  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.8.1.15  91/03/12  15:41:07  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.8.1.14  91/02/19  14:50:13  keith
 * Rewrote loop to work around titan compiler bug
 * 
 * Revision 1.8.1.13  90/12/19  12:04:56  keith
 * Added test to protect against infinite looping for velocity non-convergence.
 * 
 * Revision 1.8.1.12  90/10/25  18:05:55  keith
 * Modified rescale() to correctly handle separate scaling in case
 * of framework with no dynamics.  Also call thermalise if T=0.
 * 
 * Revision 1.8.1.11  90/10/23  20:12:09  keith
 * Added dummy function call to inhibit vectorization.
 * This allows use of 'ivdep' compiler options and also
 * works round certain bugs in cray's scc compiler.
 * 
 * Revision 1.8.1.10  90/10/16  14:47:06  keith
 * Workaround added to inhibit (incorrect) vectorization of loop 
 * at line 411 under cray scc 1.0
 * 
 * Revision 1.8.1.9  90/08/02  15:51:53  keith
 * Modified to exclude framework-framework interactions.
 * N.B. Excluded from pe and stress but NOT forces (as they sum to 0).
 * 
 * Revision 1.8.1.8  90/07/16  15:54:51  keith
 * Fixed bugs in constant-stress code
 * 
 * Revision 1.8.1.7  90/05/16  18:38:30  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.8.1.6  90/05/16  14:18:58  keith
 * *** empty log message ***
 * 
 * Revision 1.8.1.5  90/04/16  18:18:58  keith
 * Changed call of "rahman" by adding "strain-mask" parameter.
 * 
 * Revision 1.8.1.4  90/01/15  12:23:27  keith
 * Corrected declaration of arralloc from void* to char* to keep lint happy.
 * 
 * Revision 1.8.1.3  89/12/22  19:31:56  keith
 * New version of arralloc() orders memory so that pointers come FIRST.
 * This means you can simply free() the pointer returned (if l.b. = 0).
 * 
 * Revision 1.8.1.2  89/12/22  11:14:38  keith
 * Reversed indices in 'site' and 'site_force' to allow stride of 1 in ewald.
 * 
 * Revision 1.8.1.1  89/10/06  16:22:51  keith
 * Make_sites() modified to wrap sites of framework back into MD box.
 * 
 * Revision 1.8  89/09/04  18:51:43  keith
 * Made De Leeuw surface dipole term optional (off by default).
 * This term SHOULD NOT be included for ionic systems.
 * 
 * Revision 1.7  89/08/10  17:28:05  keith
 * Fixed if statement so that rdf's started on rather than after 'begin-rdf'
 * 
 * Revision 1.6  89/07/03  18:17:25  keith
 * Made code to add dipole energy and force conditional, like Ewald.
 * 
 * Revision 1.5  89/06/22  15:42:17  keith
 * Tidied up loops over species to use one pointer as countes
 * 
 * Revision 1.4  89/06/14  18:18:12  keith
 * Removed call to SCILIB function VCOPY and equivalents - use memcpy instead.
 * 
 * Revision 1.3  89/05/22  18:37:04  keith
 * Added option to scale velocities of each species separately
 * 
 * Revision 1.2  89/04/20  17:49:08  keith
 * Added code for surface dipole part of Ewald sum (After De Leeuw et al).
 * 
 * Revision 1.1  89/04/20  15:58:41  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/accel.c,v 2.4 93/12/16 18:14:11 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<math.h>
#include	"string.h"
#include	"stddef.h"
#ifdef DEBUG6
#include	<stdio.h>
#endif
/*========================== program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
void            step_1();	       /* Step co-ordinates by Beeman algrthm */
void            step_2();	       /* Step velocities at above            */
void            beeman_2();	       /* As above for individual components  */
void            make_sites();	       /* Construct site coordinate arrays    */
void            mol_force();	       /* Calculare molecular from site force */
void            mol_torque();	       /* Calculate torques from site forces  */
void            rotate();	       /* Perform rotations given quaternions */
void            newton();	       /* Calculate accelerations from forces */
void            euler();	       /* Get quat 2nd derivatives            */
void            parinello();	       /* Get correction to c of m accns      */
void            rahman();	       /* Get h 2nd derivatives	              */
void            energy_dyad();	       /* Calculate mvv for stress term       */
void            force_calc();	       /* Calculate direct-space forces       */
double          dist_pot();	       /* Returns integrated potential fn     */
void            ewald();	       /* Get Ewald sum forces		      */
void            dump();		       /* Maintain and write dump data files  */
void            zero_real();	       /* Clear area of memory		      */
void            zero_double();	       /* Clear area of memory		      */
void            invert();	       /* Matrix inverter		      */
double          det();		       /* Returns matrix determinant	      */
void            mat_vec_mul();	       /* 3 x 3 Matrix by Vector multiplier   */
void            mean_square();	       /* Caluculates mean square of args     */
void            rdf_calc();	       /* Accumulate and bin rdf	      */
double 		value();	       /* Return thermodynamic average	      */
double 		roll_av();	       /* Return thermodynamic average	      */
double          vdot();		       /* Fast vector dot product	      */
void            vscale();	       /* Vector by constant multiply	      */
double          vec_dist();	       /* normalised vector distance	      */
void		thermalise();	       /* Randomize velocities to given temp  */
void	inhibit_vectorization();       /* Self-explanatory dummy              */
#if defined(ANSI) || defined(__STDC__)
gptr		*arralloc(size_t,int,...); /* Array allocator		      */
void		note(char *, ...);	/* Write a message to the output file */
void		message(int *, ...);	/* Write a warning or error message   */
#else
gptr		*arralloc();	        /* Array allocator		      */
void		note();			/* Write a message to the output file */
void		message();		/* Write a warning or error message   */
#endif
/*========================== External data references ========================*/
extern contr_mt control;
/*========================== Macros ==========================================*/
#define ITER_MAX 10
#define	CONVRG	1.0e-7
/*  Can't rely on ANSI yet. */
#ifndef DBL_MIN
#   define DBL_MIN 1.0e-36
#endif
/*============================================================================*/
/******************************************************************************
 *   rescale    rescale velocities and quaternion derivatives to desired temp.*
 *   Exact behaviour is controlled by flag "control.scale_options".	      *
 *   This is a bit flag with the following meanings:			      *
 *	bit 0:	scale temperature for each species separately.		      *
 *	bit 1:  scale rotational and translational velocities separately      *
 *      bit 2:	use rolling averages rather than instantaneous "temperatures" *
 *	bit 3:  don't scale at all, but re-inititlize from MB distribution.   *
 ******************************************************************************/
void
rescale(system, species)
system_mp	system;
spec_mp		species;
{
   spec_mp	spec;
   int		ispec;
   double 	*temp_value = dalloc(2*system->nspecies);
   double	min_temp=MIN(value(t_n,0),roll_av(t_n,0));
   double	rtemp = 0.0, ttemp = 0.0, scale;
   int		tdof=0, rdof=0;

   /*
    *  First get trans. and rot. "temperatures", either instantaneous
    *  or rolling averages for each species individually.
    */
   for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
   {
      if( control.scale_options & 0x4 )
      {
	 temp_value[2*ispec  ] = roll_av(tt_n,ispec);
	 temp_value[2*ispec+1] = roll_av(rt_n,ispec);
      }
      else
      {
	 temp_value[2*ispec  ] = value(tt_n,ispec);
	 temp_value[2*ispec+1] = value(rt_n,ispec);
      }
      if( ! spec->framework )
      {
	 if( temp_value[2*ispec  ] < min_temp )
	    min_temp = temp_value[2*ispec  ];
	 if( spec->rdof > 0 && temp_value[2*ispec+1] < min_temp )
	    min_temp = temp_value[2*ispec+1];
      }
   }

   /*
    *  Re initialize from Maxwell-Boltzmann distribution if explicitly
    *  requested, or if any temperature is zero.
    */
   if( min_temp < 1.0e5*DBL_MIN ||
       (control.scale_options & 0x8) )
   {
      thermalise(system, species);
      return;
   }

   /*
    *  Get average of translational and rotational temps (per species)
    */
   if( ! (control.scale_options & 0x2) )
      for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
	 if( ! spec->framework )
	    temp_value[2*ispec  ] = temp_value[2*ispec+1] = 
	       (3*temp_value[2*ispec  ] + spec->rdof*temp_value[2*ispec+1]) /
		  (3+spec->rdof);
	 
   /*
    *  Perform average over species if scaling together.
    */
   if( ! (control.scale_options & 0x1) )
   {
      for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
	 if( ! spec->framework )
	 {
	    ttemp += 3*spec->nmols*temp_value[2*ispec  ]; 
	    tdof += 3*spec->nmols;
	    rtemp += spec->rdof*spec->nmols*temp_value[2*ispec+1]; 
	    rdof += spec->rdof*spec->nmols;
	 }
      ttemp /= tdof;
      if( rdof > 0 )
	 rtemp /= rdof;
      for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
	 if( ! spec->framework )
	 {
	    temp_value[2*ispec  ] = ttemp;
	    temp_value[2*ispec+1] = rtemp;
	 }
   }

   /*
    *  Actually do the scaling
    */
#ifdef DEBUG6
   fprintf(stderr,"Trans T\t\tRot T\n");
#endif
   for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
   {
#ifdef DEBUG6
      fprintf(stderr,"%8.2f\t%8.2f\n",
	      temp_value[2*ispec],temp_value[2*ispec+1]);
#endif
      if( ! spec->framework )
      {
	 scale = sqrt(control.temp / temp_value[2*ispec]);
	 vscale(3 * spec->nmols,   scale, spec->vel[0], 1);
	 if( spec->rdof > 0 )
	 { 
	    scale = sqrt(control.temp / temp_value[2*ispec+1]);
	    vscale(4 * spec->nmols, scale, spec->qdot[0], 1);
	    vscale(4 * spec->nmols, scale*scale, spec->qddot[0], 1);
	    vscale(4 * spec->nmols, scale*scale, spec->qddoto[0], 1);
	 }
	 
      }
   }
   tfree((gptr*)temp_value);
}
/******************************************************************************
 *  distant_const     Return the constant part of the distant-potential term  *
 *  c = - 2 pi sum i sum j Ni Nj Aij(cutoff), where i,j are site types, Ni,Nj *
 *  are their populations and Aij(r) is pot'l integrated from r to infinity.  *
 ******************************************************************************/
static double 
distant_const(system, species, potpar, cutoff)
system_mp       system;
spec_mt         species[];
pot_mt          potpar[];
double          cutoff;
{
   int             isite, id, jd;	/* Counters		      */
   spec_mp         spec;	       		/* pointer to current species */
   int            *site_count = ialloc(system->max_id);	/* Numbers of each site
							 * type */
   double          c = 0.0;	       		/* Accumulator for result     */
/*
 * Count the sites
 */
   (void)memset((char*)site_count,0,system->max_id*sizeof(int));
   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
NOVECTOR
      for (isite = 0; isite < spec->nsites; isite++)
      {
	 inhibit_vectorization();
	 site_count[spec->site_id[isite]] += spec->nmols;
      }
   }

   for (id = 1; id < system->max_id; id++)	/* Loops for sum over i,j     */
      for (jd = 1; jd < system->max_id; jd++)
	 c -= 2 * PI * site_count[id] * site_count[jd]
	    * dist_pot(potpar[id + system->max_id*jd].p, cutoff, system->ptype);

   xfree(site_count);
   return (c);
}
/******************************************************************************
 *  Shuffle    move down the 'acceleration' co-ordinates                      *
 *  current->old, old->very old, very old->oblivion to make room for the new  *
 *  ones at the next timestep.  Only the pointers are actually moved          *
 ******************************************************************************/
static vec_mp   v_tmp;
static quat_mp  q_tmp;
#define shuffle(a, ao, avo, tmp)	{tmp = avo; \
					avo = ao; \
					ao  = a; \
					a   = tmp; }

/******************************************************************************
 *   do_step       This routine controls the main part of the calculation for *
 *   each timestep.   It performs the following actions:  It                  *
 *   1)  Allocates space for the dynamic arrays for site co-ordinates,        *
 *       site forces, molecular forces and torques.                           *
 *   2)  Builds the array, site, of site co-ordinates from the principal-frame*
 *       co-ordinates, c of m positions and molecular quaternions.            *
 *   3)  Steps the co-ordinates using the Beeman algorithm		      *
 *   4)  Calls the main inter-site force calculating routine, force_calc.     *
 *   5)  Calls the Ewald sum routine for charged systems.                     *
 *   6)  Evaluates the molecular forces and torques from the site forces      *
 *   7)  'shuffles' all the accelerations down to make room for the new ones  *
 *        about to be evaluated.                                              *
 *   8)  Applies Newton's equations of motion to calculate accelerations,     *
 *       Euler's equations and 2nd order quaternion method to get the         *
 *       quaternion accelerations and calculates the h-matrix accelerations   *
 *       in the Parinello and Rahman constant-stress method.                  *
 *   9)  Applies the P & R correction to the c of m accelerations.            *
 *   10) Steps the velocities using the Beeman algorithm,		      *
 *       and iterates steps 8-10 until convergence.			      *
 *   11) Deallocates all the dynamic arrays and returns the space to the heap *
 ******************************************************************************/
void 
do_step(sys, species, site_info, potpar, meansq_f_t, pe, dip_mom, stress)
system_mp       sys;		       /* Pointer to system info	 (in) */
spec_mt         species[];	       /* Array of species info	 (in) */
site_mt         site_info[];	       /* Array of site info structures (in) */
pot_mt          potpar[];	       /* Array of potential parameters (in) */
vec_mt          meansq_f_t[][2];       /* Mean square force and torque (out) */
double          pe[];		       /* Potential energy real/Ewald  (out) */
vec_mt          dip_mom;	       /* Total system dipole moment   (out) */
mat_mt          stress;		       /* Virial part of stress	(out) */
{
/*
 * The following declarations are arrays of pointers to the forces
 * etc for each species.  That is force[i] is a pointer to the force on
 * molecule 0 of species i
 */
   vec_mp	*force = palloc(sys->nspecies),
		*torque = palloc(sys->nspecies);
   real		***site_sp = (real***)arralloc(sizeof(real*), 2,
					       0, sys->nspecies-1, 0, 2),
  		***force_sp = (real***)arralloc(sizeof(real*), 2,
					       0, sys->nspecies-1, 0, 2); 
/*
 * The following declarations are pointers to the force etc for the whole
 * system, and are set equal to (eg) force[0]
 */
   vec_mp          force_base = ralloc(sys->nmols),
		   torque_base = ralloc(sys->nmols_r);
   real		**site = (real**)arralloc(sizeof(real), 2,
					  0, 2, 0, sys->nsites-1),
   		**site_force = (real**)arralloc(sizeof(real), 2,
						0, 2, 0, sys->nsites-1);
/*
 * Other local variables
 */
   real           *chg = dalloc(sys->nsites), *chg_ptr;
   vec_mp		   c_of_m = ralloc(sys->nmols);
   spec_mp         spec, fspec /*Framework species */;
   int             nspecies = sys->nspecies;
   int             ispec, imol, imol_r, isite;
   int             i, j;
   static boolean  init = true;
   static double   dist;
   double          vol = det(sys->h);
   int             iter;
   mat_mt          ke_dyad, hinv;
   quat_mp         qd_tmp;	       /* Temporary for velocities   	      */
   vec_mp          acc_tmp, vel_tmp;   /* Temporaries for iteration	      */
   int		   nsitesxf, nmolsxf;  /* Count of non-framework sites, mols. */
/*
 * Initialisation of distant potential constant - executed first time only
 */
   if (init)
   {
      dist = distant_const(sys, species, potpar, control.cutoff);
      note("Distant potential correction = %f, Pressure correction = %f",
	   CONV_E * dist / vol, CONV_P * dist / (vol * vol));
      init = false;
   }
/*
 * The next chunk of code sets up the dynamic arrays for the centre of mass
 * forces, torques, site vectors and site forces.  It is a little complex but
 * results in great simplicity in all of the called routines. The arrays
 * themselves are indicated by pointers *_base (*=force etc), and have
 * dimensions [n][3], n= number of molecules, polyatomics and sites. They are
 * subdivided into the appropriate segments for each molecular species eg
 * site_base[nspecies][nsites][3], BUT nsites depends on the species so it is
 * NOT a 3 dimensional array.  Instead the array of pointer 'site' is used.
 * Site[ispec] is a pointer to the sub-array for species ispec of dimension
 * [nsites][3]. If any or all of the molecular species are monatomic, then no
 * space is allocated for the torque array, and all the pointers are NULL.
 */
/*
 *  Count non-framework sites and molecules.
 */
   fspec = species; nsitesxf = 0; nmolsxf = 0;
   while (fspec < species+nspecies && ! fspec->framework)
   {
      nsitesxf += fspec->nsites * fspec->nmols;
      nmolsxf  += fspec->nmols;
      fspec++;
   }   
/*
 * Set up arrays of pointers to sites, forces etc for each species
 */
   isite = 0; imol = 0; imol_r = 0;
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {
      force[ispec] = force_base+imol;
      if (spec->quat)
	 torque[ispec] = torque_base+imol_r;
      for( i = 0; i < 3; i++ )
      {
	 site_sp[ispec][i] = site[i]+isite;
	 force_sp[ispec][i] = site_force[i]+isite;
      }
      imol += spec->nmols;
      if (spec->quat)
	 imol_r += spec->nmols;
      isite += spec->nmols * spec->nsites;
   }
/*
 * Set up array of site charges
 */
   chg_ptr = chg;
   for (spec = species; spec < species+nspecies; spec++)
      for (imol = 0; imol < spec->nmols; imol++)
	 for (isite = 0; isite < spec->nsites; isite++)
	    *chg_ptr++ = site_info[spec->site_id[isite]].charge;
   mat_vec_mul(sys->h, sys->c_of_m, c_of_m, sys->nmols);
/*
 * Set some accumulators to zero
 */
   zero_real(stress[0], 9);	       /* Initialise stress tensor   */
   zero_real(meansq_f_t[0][0], 6 * sys->nspecies);
   zero_real(site_force[0], 3*sys->nsites);
   zero_double(pe, NPE);
/*
 * Initial co-ordinate step of Beeman algorithm.
 */
   step_1(sys);
/*
 * Calculate the site positions at this timestep - loop over species
 */
   for (spec = species; spec < &species[nspecies]; spec++)
   {
      make_sites(sys->h, spec->c_of_m, spec->quat, spec->p_f_sites,
		spec->framework,site_sp[spec-species],spec->nmols,spec->nsites);
#ifdef DEBUG1
   { int is;
     printf("%s co-ordinates\n",spec->name);
     for(is = 0; is < spec->nsites*spec->nmols; is++)
	printf("%24.15f %24.15f %24.15f\n", site_sp[spec-species][0][is],
	       		                   site_sp[spec-species][1][is],
	                                   site_sp[spec-species][2][is]);
  }
#endif
   }


/*
 * Real-space part of force evaluation - no loop over species for efficiency
 */
   force_calc(site, site_force, sys, species, chg, potpar, pe, stress);

/*
 * Reciprocal-space part of Ewald sum
 */
   if (control.alpha > ALPHAMIN)
   {
      ewald(site, site_force, sys, species, chg, pe + 1, stress);
/*
 * Dipole moment contribution to forces and potential (De Leeuw, Perram
 * and Smith Proc Roy Soc A373, 27-56 (1980)
 */
      for (i = 0; i < 3; i++)
	 dip_mom[i] = vdot(sys->nsites, site[i], 1, chg, 1);
      if( control.surface_dipole )
      {
	 for (i = 0; i < 3; i++)
	    for ( isite = 0; isite < sys->nsites; isite++ )
	       site_force[i][isite] -= 4.0*PI/(3.0*vol)*dip_mom[i] * chg[isite];
	 
	 pe[1] += 2.0*PI/(3.0*vol) * SUMSQ(dip_mom);
      }
   }

/*
 * Calculate the centre of mass forces and torques from the site forces
 */
   for (spec = species; spec < &species[nspecies]; spec++)
   {
      ispec = spec-species;
      mol_force(force_sp[ispec], force[ispec], spec->nsites, spec->nmols);
      if (spec->rdof > 0)
	 mol_torque(force_sp[ispec], spec->p_f_sites,
		    torque[ispec], spec->quat, spec->nsites, spec->nmols);
   }

/*
 * Add correction term to convert from site to molecular virial
 */
   for (i = 0; i < 3; i++)
   {
      for (j = i + 1; j < 3; j++)
	 stress[j][i] = stress[i][j];
      for (j = 0; j < 3; j++)
      {
        /*
	 * Non-framework sites sum f.s = sum f.r - sum F.R
	 */
	 stress[i][j] -=
		  vdot(nsitesxf, site_force[i], 1, site[j], 1)
	        - vdot(nmolsxf, force_base[0] + i, 3, c_of_m[0] + j, 3);
	/*
         *  Framework sites -- can't use "site" array as sites have been
         *  relocated by pbc's back into md cell.
	 */
	if( fspec < species+nspecies )
	   for(imol=0; imol < fspec->nmols; imol++)
	      stress[i][j] -= vdot(fspec->nsites,
				   site_force[i]+nsitesxf+imol*fspec->nsites,1,
				   fspec->p_f_sites[0]+j,3);
      }
		                            
   }

/*
 * Calculate distant stress and potential terms and add
 */
   pe[0] += dist / vol;
   for (i = 0; i < 3; i++)
      stress[i][i] += dist / vol;

/*
 * Shuffle the accelerations (acc -> acco, acco-> accvo, accvo -> acc) Don't
 * actually move the data - just the pointers
 */
   shuffle(sys->acc, sys->acco, sys->accvo, v_tmp);
   shuffle(sys->qddot, sys->qddoto, sys->qddotvo, q_tmp);
   if (control.const_pressure)
      shuffle(sys->hddot, sys->hddoto, sys->hddotvo, v_tmp);

   for (spec = species; spec < &species[nspecies]; spec++)
   {
      inhibit_vectorization();      /* Inhibits (incorrect) vectorization */
      shuffle(spec->acc, spec->acco, spec->accvo, v_tmp);
      if (spec->rdof > 0)
	 shuffle(spec->qddot, spec->qddoto, spec->qddotvo, q_tmp);
   }

/*
 * Now apply the Newton/Euler equations to find the accelerations and
 * quaternion second derivatives.
 */
   for (spec = species; spec < &species[nspecies]; spec++)
      newton(force[spec-species], spec->acc, spec->mass, spec->nmols);

/*
 * Correction to centre of mass accelerations for constant pressure algorithm
 * First get scaled accelerations by multiplying by the inverse h matrix and
 * add P&R term. Then calculate the 'accelerations' of the unit cell matrix
 */
   invert(sys->h, hinv);
   mat_vec_mul(hinv, sys->acc, sys->acc, sys->nmols);
   if (control.const_pressure)
   {
      zero_real(ke_dyad[0], 9);
      for (spec = species; spec < &species[nspecies]; spec++)
	 energy_dyad(ke_dyad, sys->h, spec->velp, spec->mass, spec->nmols);
      rahman(stress, sys->h, sys->hddot, ke_dyad,
	     control.pressure, control.pmass, control.strain_mask);
   }

/*
 * Iterate velocity dependant parts with beeman step 2 until convergence
 */
   if (sys->nmols_r > 0)
   {
      iter = 0;
      qd_tmp = qalloc(sys->nmols_r);
      do
      {
	 iter++;
	 if(iter > ITER_MAX)
	    message(NULLI, NULLP, FATAL, NCNVRG, iter, 
		    vec_dist(qd_tmp[0], sys->qdotp[0], 4 * sys->nmols_r));
	 for (spec = species; spec < &species[nspecies]; spec++)
	    if (spec->rdof > 0)
	       euler(torque[spec-species], spec->quat, spec->qdotp,
		     spec->qddot, spec->inertia, spec->nmols);
	 memcp(qd_tmp[0], sys->qdotp[0], 4 * sys->nmols_r * sizeof(real));
	 beeman_2(sys->qdot[0], sys->qdotp[0], sys->qddot[0], sys->qddoto[0],
		  sys->qddotvo[0], 4 * sys->nmols_r);
      } while (vec_dist(qd_tmp[0], sys->qdotp[0], 4 * sys->nmols_r) > CONVRG);
#ifdef DEBUG
      printf("Quaternion derivatives converged in %d iterations \n", iter);
#endif
      xfree(qd_tmp);
   }
   if (control.const_pressure)
   {
      acc_tmp = ralloc(sys->nmols);
      vel_tmp = ralloc(sys->nmols);
      beeman_2(sys->hdot[0], sys->hdotp[0], sys->hddot[0], sys->hddoto[0],
	       sys->hddotvo[0], 9);
      do
      {
	 parinello(sys->h, sys->hdotp, sys->velp, sys->acc, acc_tmp,
		   sys->nmols);
	 memcp(vel_tmp[0], sys->velp[0], 3 * sys->nmols * sizeof(real));
	 beeman_2(sys->vel[0], sys->velp[0], acc_tmp[0], sys->acco[0],
		  sys->accvo[0], 3 * sys->nmols);
      } while (vec_dist(vel_tmp[0], sys->velp[0], 3 * sys->nmols) > CONVRG);
      memcp(sys->acc, acc_tmp, 3*sys->nmols * sizeof(real));
      xfree(acc_tmp);
      xfree(vel_tmp);
   }
/*
 *  Apply constraint to any framework molecules.
 */
   for (spec = species; spec < &species[nspecies]; spec++)
      if( spec->framework )
      {
	 zero_real(spec->acc[0], 3*spec->nmols);
	 zero_real(spec->vel[0], 3*spec->nmols);
      }
/*
 * Final MD update step
 */
   step_2(sys);

/*
 * Calculate mean-square forces and torques
 */
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {
      mean_square(force[ispec], meansq_f_t[ispec][0], spec->nmols);
      if (spec->rdof > 0)
	 mean_square(torque[ispec], meansq_f_t[ispec][1], spec->nmols);
   }

/*
 * Accumulate radial distribution functions
 */
   if (control.rdf_interval > 0 && control.istep >= control.begin_rdf &&
	 control.istep % control.rdf_interval == 0)
      rdf_calc(site, sys, species);

/*
 * Perform periodic dump of dynamic data
 */
   if (control.dump_interval > 0 && control.dump_level != 0 &&
	 control.istep >= control.begin_dump &&
	 (control.istep - control.begin_dump) % control.dump_interval == 0)
      dump(sys, force_base, torque_base, stress, pe[0] + pe[1]);

/*
 * Deallocate the dynamic storage before exiting
 */
   xfree(force);
   xfree(torque);
   xfree(site);
   xfree(site_force);
   xfree(force_base);
   if (torque_base != NULL)
      xfree(torque_base);
   xfree(site_sp);
   xfree(force_sp);
   xfree(chg);
   xfree(c_of_m);
}
/*******************************************************************************
 *  Mol_radius.  Determine and return the greatest molecular radius of any     *
 *  species in the system.  That is, the largest c-of-mass - site distance.    *
 *******************************************************************************/
double mol_radius(species, nspecies)
spec_mt	species[];
int	nspecies;
{
   spec_mp spec;
   static double radius = -1.0;
   double r;
   int	isite;

   if( radius >= 0.0 )
      return radius;

   radius = 0.0;
   for(spec = species; spec < species+nspecies; spec++)
      if( !spec->framework )
      {
	 for( isite = 0; isite < spec->nsites; isite++ )
	 {
	    r = SUMSQ(spec->p_f_sites[isite]);
	    radius = MAX(radius, r);
	 }
      }
   return radius;
}

   
