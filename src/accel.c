/******************************************************************************
 * Accel     This file contains the 'driver' function for a single timestep - *
 *           "do_step" and other functions which need access to the system    *
 *           and species structs, "rescale" and "distant_const". 	      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	accel.c,v $
 * Revision 1.2  89/04/20  17:49:08  keith
 * Added code for surface dipole part of Ewald sum (After De Leeuw et al).
 * 
 * Revision 1.1  89/04/20  15:58:41  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: accel.c,v 1.3 89/05/22 14:04:42 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#include	<math.h>
/*========================== program include files ===========================*/
#include "structs.h"
/*========================== Library declarations ============================*/
void            cfree();	       /* Free allocated memory	      */
/*========================== External function declarations ==================*/
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
void            invert();	       /* Matrix inverter		      */
double          det();		       /* Returns matrix determinant	      */
void            mat_vec_mul();	       /* 3 x 3 Matrix by Vector multiplier   */
void            mean_square();	       /* Caluculates mean square of args     */
void            rdf_calc();	       /* Accumulate and bin rdf	      */
double 		value();	       /* Return thermodynamic average	      */
void            vcopy();	       /* Vector copy.			      */
double          vdot();		       /* Fast vector dot product	      */
void            vscale();	       /* Vector by constant multiply	      */
double          vec_dist();	       /* normalised vector distance	      */
void            note();		       /* Write message to output file	      */
/*========================== External data references ========================*/
extern contr_t  control;
/*========================== Macros ==========================================*/
#define	CONVRG	1.0e-7
/*============================================================================*/
/******************************************************************************
 *   Scale_spec	rescale for an individual species			      *
 ******************************************************************************/
void 
scale_spec(scale, spec)
double		scale;		       	/* Scaling factor		 (in) */
spec_p		spec;		      	/* System info struct 	     (in/out) */
{
   int		nmols   = spec->nmols,
   		nmols_r = (spec->rdof ? spec->nmols : 0);
   vscale(3 * nmols,   scale, spec->vel[0], 1);
   vscale(4 * nmols_r, scale, spec->qdot[0], 1);
   vscale(4 * nmols_r, scale*scale, spec->qddot[0], 1);
   vscale(4 * nmols_r, scale*scale, spec->qddoto[0], 1);
}
/******************************************************************************
 *   rescale    rescale velocities and quaternion derivatives to desired temp.*
 ******************************************************************************/
void
rescale(system, species)
system_p	system;
spec_p		species;
{
   spec_p	spec;
   double scale_factor;

   scale_factor = control.temp / value(t_n, 0);	/* `Global' scale factor      */

   for(spec = species; spec < &species[system->nspecies]; spec++)
   {
      if ( control.scale_separately )		/* Per-species scale factor   */
         scale_factor = (3+spec->rdof) * control.temp / 
                        (  3.0 *      value(tt_n, spec-species)
		         + spec->rdof*value(rt_n, spec-species));
      scale_spec(sqrt(scale_factor), spec);
   }
}
/******************************************************************************
 *  distant_const     Return the constant part of the distant-potential term  *
 *  c = - 2 pi sum i sum j Ni Nj Aij(cutoff), where i,j are site types, Ni,Nj *
 *  are their populations and Aij(r) is pot'l integrated from r to infinity.  *
 ******************************************************************************/
static double 
distant_const(system, species, potpar, cutoff)
system_p        system;
spec_t          species[];
pot_t           potpar[];
double          cutoff;
{
   int             ispec, isite, id, jd;	/* Counters		      */
   spec_p          spec;	       		/* pointer to current species */
   int            *site_count = ialloc(system->max_id);	/* Numbers of each site
							 * type */
   double          c = 0.0;	       		/* Accumulator for result     */
/*
 * Count the sites
 */
   for (ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
      for (isite = 0; isite < spec->nsites; isite++)
	 site_count[spec->site_id[isite]] += spec->nmols;

   for (id = 1; id < system->max_id; id++)	/* Loops for sum over i,j     */
      for (jd = 1; jd < system->max_id; jd++)
	 c -= 2 * PI * site_count[id] * site_count[jd]
	    * dist_pot(potpar[id + system->max_id*jd].p, cutoff, system->ptype);

   cfree((char *) site_count);
   return (c);
}
/******************************************************************************
 *  Shuffle    move down the 'acceleration' co-ordinates                      *
 *  current->old, old->very old, very old->oblivion to make room for the new  *
 *  ones at the next timestep.  Only the pointers are actually moved          *
 ******************************************************************************/
static vec_p    v_tmp;
static quat_p   q_tmp;
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
system_p        sys;		       /* Pointer to system info	 (in) */
spec_t          species[];	       /* Array of species info	 (in) */
site_t          site_info[];	       /* Array of site info structures (in) */
pot_t           potpar[];	       /* Array of potential parameters (in) */
vec_t           meansq_f_t[][2];       /* Mean square force and torque (out) */
double          pe[];		       /* Potential energy real/Ewald  (out) */
vec_t           dip_mom;	       /* Total system dipole moment   (out) */
mat_t           stress;		       /* Virial part of stress	(out) */
{
/*
 * The following declarations are arrays of pointers to the   * sites, forces
 * etc for each species.  That is force[i]      * is a pointer to the force on
 * molecule 0 of species i
 */
   vec_p          *force = palloc(sys->nspecies),
		  *torque = palloc(sys->nspecies),
		  *site = palloc(sys->nspecies),
		  *site_force = palloc(sys->nspecies);
/*
 * The following declarations are pointers to the force etc   * for the whole
 * system, and are set equal to (eg) force[0]
 */
   vec_p           force_base = ralloc(sys->nmols),
		   torque_base = ralloc(sys->nmols_r),
		   site_base = ralloc(sys->nsites),
		   s_f_base = ralloc(sys->nsites);
/*
 * Temporary pointers corresponding to above 
 */
   vec_p           force_ptr, torque_ptr, site_ptr, s_f_ptr;

   real           *chg = dalloc(sys->nsites), *chg_ptr;
   vec_p		   c_of_m = ralloc(sys->nmols);
   spec_p          spec;
   int             nspecies = sys->nspecies;
   int             ispec, imol, isite;
   int             i, j;
   static boolean  init = true;
   static double   dist;
   double          vol = det(sys->h);
   double         *ppe;
#ifdef DEBUG
   int             iter;
#endif
   mat_t           ke_dyad, hinv;
   quat_p          qd_tmp;	       /* Temporary for velocities   	      */
   vec_p           acc_tmp, vel_tmp;   /* Temporaries for iteration	      */
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
 * Set up arrays of pointers to sites, forces etc for each species
 */
   force_ptr = force_base;
   torque_ptr = torque_base;
   site_ptr = site_base;
   s_f_ptr = s_f_base;
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {
      force[ispec] = force_ptr;
      torque[ispec] = torque_ptr;
      site[ispec] = site_ptr;
      site_force[ispec] = s_f_ptr;
      force_ptr += spec->nmols;
      if (species[ispec].quat)
	 torque_ptr += spec->nmols;
      site_ptr += spec->nmols * spec->nsites;
      s_f_ptr += spec->nmols * spec->nsites;
   }
/*
 * Set up array of site charges
 */
   chg_ptr = chg;
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
      for (imol = 0; imol < spec->nmols; imol++)
	 for (isite = 0; isite < spec->nsites; isite++)
	    *chg_ptr++ = site_info[spec->site_id[isite]].charge;
   mat_vec_mul(sys->h, sys->c_of_m, c_of_m, sys->nmols);
/*
 * Set some accumulators to zero
 */
   zero_real(stress[0], 9);	       /* Initialise stress tensor   */
   zero_real(meansq_f_t[0][0], 6 * sys->nspecies);
   for (ppe = pe; ppe < pe + NPE; ppe++)
      *ppe = 0.0;
/*
 * Initial co-ordinate step of Beeman algorithm.
 */
   step_1(sys);
/*
 * Calculate the site positions at this timestep - loop over species
 */
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
      make_sites(sys->h, spec->c_of_m, spec->quat, spec->p_f_sites,
		 site[ispec], spec->nmols, spec->nsites);

/*
 * Real-space part of force evaluation - no loop over species for efficiency
 */
   force_calc(site_base, s_f_base, sys, species, chg, potpar, pe, stress);

/*
 * Reciprocal-space part of Ewald sum
 */
   if (control.alpha > ALPHAMIN)
      ewald(site_base, s_f_base, sys, species, chg, pe + 1, stress);
/*
 * Dipole moment contribution to forces and potential (De Leeuw, Perram
 * and Smith Proc Roy Soc A373, 27-56 (1980)
 */
   for (i = 0; i < 3; i++)
   {
      dip_mom[i] = vdot(sys->nsites, site_base[0] + i, 3, chg, 1);
      for ( isite = 0; isite < sys->nsites; isite++ )
         s_f_base[isite][i] -= 4.0*PI/(3.0*vol) * dip_mom[i] * chg[isite];
   }
   pe[1] += 2.0*PI/(3.0*vol) * SUMSQ(dip_mom);

/*
 * Calculate the centre of mass forces and torques from the site forces
 */
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {
      mol_force(site_force[ispec], force[ispec], spec->nsites, spec->nmols);
      if (spec->rdof > 0)
	 mol_torque(site_force[ispec], spec->p_f_sites,
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
	 stress[i][j] -=
		  vdot(sys->nsites, s_f_base[0] + i, 3, site_base[0] + j, 3)
	        - vdot(sys->nmols, force_base[0] + i, 3, c_of_m[0] + j, 3);
   }

/*
 * Calculate distant stress and potential terms and add
 */
   pe[0] += dist / vol;
   for (i = 0; i < 3; i++)
      stress[i][i] += dist / (vol * vol);

/*
 * Shuffle the accelerations (acc -> acco, acco-> accvo, accvo -> acc) Don't
 * actually move the data - just the pointers
 */
   shuffle(sys->acc, sys->acco, sys->accvo, v_tmp);
   shuffle(sys->qddot, sys->qddoto, sys->qddotvo, q_tmp);
   if (control.const_pressure)
      shuffle(sys->hddot, sys->hddoto, sys->hddotvo, v_tmp);
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {
      shuffle(spec->acc, spec->acco, spec->accvo, v_tmp);
      if (spec->rdof > 0)
	 shuffle(spec->qddot, spec->qddoto, spec->qddotvo, q_tmp);
   }

/*
 * Now apply the Newton/Euler equations to find the accelerations and
 * quaternion second derivatives.
 */
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
      newton(force[ispec], spec->acc, spec->mass, spec->nmols);

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
      for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
	 energy_dyad(ke_dyad, sys->h, spec->velp, spec->mass, spec->nmols);
      rahman(stress, sys->h, sys->hddot, ke_dyad,
	     control.pressure, control.pmass);
   }

/*
 * Iterate velocity dependant parts with beeman step 2 until convergence
 */
   if (sys->nmols_r > 0)
   {
#ifdef DEBUG
      iter = 0;
#endif
      qd_tmp = qalloc(sys->nmols_r);
      do
      {
#ifdef DEBUG
	 iter++;
#endif
	 for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
	    if (spec->rdof > 0)
	       euler(torque[ispec], spec->quat, spec->qdotp,
		     spec->qddot, spec->inertia, spec->nmols);
	 vcopy(4 * sys->nmols_r, sys->qdotp[0], 1, qd_tmp[0], 1);
	 beeman_2(sys->qdot[0], sys->qdotp[0], sys->qddot[0], sys->qddoto[0],
		  sys->qddotvo[0], 4 * sys->nmols_r);
      } while (vec_dist(qd_tmp[0], sys->qdotp[0], 4 * sys->nmols_r) > CONVRG);
#ifdef DEBUG
      printf("Quaternion derivatives converged in %d iterations \n", iter);
#endif
      cfree((char *) qd_tmp);
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
	 vcopy(3 * sys->nmols, sys->velp[0], 1, vel_tmp[0], 1);
	 beeman_2(sys->vel[0], sys->velp[0], acc_tmp[0], sys->acco[0],
		  sys->accvo[0], 3 * sys->nmols);
      } while (vec_dist(vel_tmp[0], sys->velp[0], 3 * sys->nmols) > CONVRG);
      cfree((char *) acc_tmp);
      cfree((char *) vel_tmp);
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
   if (control.rdf_interval > 0 && control.istep > control.begin_rdf &&
	 control.istep % control.rdf_interval == 0)
      rdf_calc(site_base, sys, species);

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
   cfree((char *) force);
   cfree((char *) torque);
   cfree((char *) site);
   cfree((char *) site_force);
   cfree((char *) force_base);
   if (torque_base != NULL)
      cfree((char *) torque_base);
   cfree((char *) site_base);
   cfree((char *) s_f_base);
   cfree((char *) chg);
   cfree((char *) c_of_m);
}
