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
 * Accel     This file contains the 'driver' function for a single timestep - *
 *           "do_step" and other functions which need access to the system    *
 *           and species structs, "rescale" and "distant_const".               *
 ******************************************************************************
 *      Revision Log
 *       $Log: accel.c,v $
 *       Revision 2.41  2002/09/19 09:26:27  kr
 *       Tidied up header declarations.
 *       Changed old includes of string,stdlib,stddef and time to <> form
 *
 *       Revision 2.40  2002/09/18 09:03:20  kr
 *       Merged H_0 NVT bugfix for ewald46 branch
 *
 *       Revision 2.39  2001/07/31 09:51:02  keith
 *       Merged H0halfstep branch.
 *       Calculates H_0 at half-step, which corresponds to on-step in VTV
 *       ordering.
 *       Now uses Nose's splitting integrator by default.
 *
 *       Revision 2.39  2001/07/31 09:40:57  keith
 *       Merged H0halfstep branch.
 *       Calculates H_0 at half-step, which corresponds to on-step in VTV
 *       ordering.
 *       Now uses Nose's splitting integrator by default.
 *
 *       Revision 2.38.2.2  2001/07/26 17:30:17  keith
 *       Version using Nose's splitting with exact solution
 *
 *       Revision 2.38.2.1  2001/07/25 14:26:59  keith
 *       Moved calculation of H_0 to half-step - still uses velocity verlet
 *       momenta.
 *
 *       Revision 2.38  2001/07/20 11:43:24  keith
 *       Added code to print out 1/2 step momenta-based kinetic energy
 *
 *       Revision 2.37  2001/07/13 16:52:53  keith
 *       Fixed calculation of H_0 to use co-ords and momenta at consistent timestep.
 *
 *       Revision 2.36  2001/07/13 14:54:17  keith
 *       Fixed bug in DEBUG_THERMOSTAT code which calculayed KE and PE
 *       at different half-steps.
 *
 *       Revision 2.35  2001/06/28 15:48:22  keith
 *       Fixed bug in implementation of Nose-poincare thermostat to do with
 *       update of smom in rotational update stage.
 *
 *       Revision 2.34  2001/05/24 16:26:43  keith
 *       Updated program to store and use angular momenta, not velocities.
 *        - added conversion routines for old restart files and dump files.
 *       Got rid of legacy 2.0 and lower restart file reading code.
 *
 *       Revision 2.33  2001/02/22 10:24:35  keith
 *       Reinstated capability of "molecular" cutoff, but still using Bekker stress.
 *
 *       Revision 2.32  2001/02/21 09:45:11  keith
 *       Modified rescale() function to reset thermostat variable to 1 and
 *       associated momentum to 0.  Every scaling now resets the thermostat
 *       to the default condition.
 *
 *       Revision 2.31  2001/02/19 19:36:44  keith
 *       First working version of combined isothermic/isobaric ensemble.
 *       (Previous version was faulty).
 *       Also includes uniform-dilation, constant-pressure mode.
 *
 *       Revision 2.30  2001/02/15 19:00:37  keith
 *       Tidied up code and corrected one error in combined barostat  and
 *       themostat.
 *
 *       Revision 2.29  2001/02/13 17:45:06  keith
 *       Added symplectic Parrinello-Rahman constant pressure mode.
 *
 *       Revision 2.28  2000/12/06 17:45:27  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.27  2000/11/15 17:51:57  keith
 *       Changed format of dump files.
 *       Added second struct with sufficient information
 *       about the simulation that most utility programs
 *       (namely those which do not need site co-ordinates)
 *       should not need to read sys-spec or restart files.
 *
 *       New options "-c -1" to dumpext prints header info.
 *       -- dumpanal removed.
 *
 *       Revision 2.26  2000/11/08 18:36:24  keith
 *       Moved calculation of H0 from beginning to midpoint.  This gives a much more
 *       accurate value -- any error will bias the temperature.
 *       Recalculates H0 after every velocity rescaling.  The thermostat should not
 *       fight with rescaling as much now.
 *
 *       Revision 2.25  2000/11/06 16:02:04  keith
 *       First working version with a Nose-Poincare thermostat for rigid molecules.
 *
 *       System header updated to include H_0.
 *       Dump performs correct scaling  of angular velocities, but dumpext still
 *          needs to be updated to read this.
 *       XDR functions corrected to work with new structs.
 *       Parallel broadcast of config also updated.
 *       Some unneccessary functions and code deleted.
 *
 *       Revision 2.22.2.3  2000/10/20 15:17:38  keith
 *       Incorporated fix og 2.15c into this branch
 *
 *       Revision 2.22.2.2  2000/10/20 13:59:31  keith
 *       Incorporated new neightbour list stuff into accel.c.
 *       Removed old "poteval" from accel.c.  Now use one in
 *       force.
 *       Other errors corrected and declarations tidied somewhat.
 *
 *       Revision 2.22.2.2  2000/10/20 11:48:50  keith
 *       Incorporated new neighbour list indexing algorithm from 2.16
 *
 *       Revision 2.22.2.1  2000/05/24 11:53:08  keith
 *       Split momentum update into 0.5*step at beginning of cycle\nand 0.5*step at end as written in the papers
 *
 *       Revision 2.21  2000/04/27 17:57:05  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.19.2.1  2000/04/24 17:05:41  keith
 *       Dullweber et al Leapfrog version
 *
 *       Revision 3.1  2000/04/14 15:00:47  keith
 *       Dullweber et al Leapfrog version
 *
 *       Revision 2.19  1999/10/08 15:49:58  keith
 *       Fully implemented new constant-pressure algorithm.
 *       Select by "const-pressure=2" in control.
 *
 *       Revision 2.18  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 *
 *       Revision 2.17  1996/11/05 16:47:19  keith
 *       Optimized site and site_force arrays to avoid cache conflicts.
 *
 *       Revision 2.16  1996/08/23 15:06:01  keith
 *       Fixed bug whereby rot elements of temp_value[] were uninitialized.
 *       This caused a crash on non-ieee machines.
 *
 *       Revision 2.15  1996/08/14 16:23:24  keith
 *       Fixed error in thermoststat implementation and integration.
 *
 *       Revision 2.14  1996/03/14 14:42:27  keith
 *       Altered "rescale()" to suntract net velocity in case of
 *       separate-species rescaling, since that doesn't conserve momentum.
 *
 *       Revision 2.13  1996/01/15 15:14:00  keith
 *       De "lint"-ed the code.
 *       Removed "old" RDF code call.
 *
 *       Revision 2.12  1995/12/07 17:43:53  keith
 *       Reworked V. Murashov's thermostat code.  Created new functions
 *       nhtherm() and gtherm() to calculate alphadot.
 *       Changed to Program units.
 *
 *       Revision 2.11  1995/12/05 11:24:57  keith
 *       Added new function "poteval" which calls "kernel" to evaluate
 *       potential at single point and modified "distant_const" to
 *       correctly evaluate long-range pressure correction.
 *
 *       Revision 2.10  1995/12/04 11:45:49  keith
 *       Nose-Hoover and Gaussian (Hoover constrained) thermostats added.
 *       Thanks to V. Murashov.
 *
 * Revision 2.9  1994/07/12  16:20:26  keith
 * Fixed bug whereby "dip_mom" left uninitialized for non-coulomb system.
 *
 * Revision 2.8  1994/07/07  16:57:01  keith
 * Updated for parallel execution on SPMD machines.
 * Interface to MP library routines hidden by par_*() calls.
 * Compile with -DSPMD to activate
 *
 * Revision 2.7  1994/06/08  13:07:39  keith
 * New version of array allocator which breaks up requests for DOS.
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
 * Revision 2.5  94/01/18  13:31:46  keith
 * Null update for XDR portability release
 * 
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
 * --Modified declaration of size_mt and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.8.1.15  91/03/12  15:41:07  keith
 * Tidied up typedefs size_mt and include file <sys/types.h>
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
static char *RCSid = "$Header: /home/moldy/CVS/moldy/src/accel.c,v 2.41 2002/09/19 09:26:27 kr Exp $";
#endif
/*========================== Library include files ===========================*/
#include        "defs.h"
/*========================== Library include files ===========================*/
#include        <math.h>
#include        <string.h>
#if defined(DEBUG10) || defined(DEBUG2) || defined(DEBUG_THERMOSTAT)
#include        <stdio.h>
#endif
/*========================== program include files ===========================*/
#include "structs.h"
/*========================== External function declarations ==================*/
gptr   *talloc(int n, size_mt size, int line, char *file); 
                              /* Interface to memory allocator       */
void   tfree(gptr *p);               /* Free allocated memory                            */
void   afree(gptr *p);               /* Free allocated array                            */
void   leapf_com(double step, vec_mt (*c_of_m), vec_mt (*mom), 
                 mat_mt h, real s, real mass, int nmols);
void   leapf_mom(double step, mat_mt, vec_mt (*mom), 
                 vec_mt (*force),  int nmols);
void   leapf_quat(double step, quat_mt (*quat), quat_mt (*amom), 
                  real *inertia, real *smom, real s, int nmols);
void   leapf_amom(double step, quat_mt (*amom), vec_mt (*torque), int nmols);
double leapf_s(double step, real s, real smom, double Q);
double leapf_smom_a(double step, real s, real smomo,  double Q, double gkt);
double leapf_smom_b(double step, real s, real smomo, double Q, double gkt);
void   leapf_hmom(double step, mat_mt hmom, mat_mt sigma, real s,
                  real pressure, int mask);
void   leapf_h(double step, mat_mt h, mat_mt hmom, real s, real W);
void   gleap_therm(double step, real mass, real gkt, real *s, real *smom);
void   leapf_nose_therm(double step, real mass, real *s, real *smom);
void   gleap_cell(double step, real pmass, real s, real pressure, int strain_mask, 
                  mat_mt h, mat_mt hmom, real *smom, boolean uniform);
void   update_hmom(double step, real s, mat_mt h,
                   mat_mt stress_part, mat_mt hmom, boolean uniform);

void   make_sites(mat_mt, vec_mp c_of_m_s, quat_mp quat, 
                  vec_mp p_f_sites, real **site, int nmols, int nsites, 
                  int molflag);                /* Construct site coordinate arrays  */
void   mol_force(real **site_force, vec_mp force, int nsites, int nmols);
                                      /* Calculate molecular from site force */
void   mol_torque(real **site_force, vec_mp site, vec_mp torque, quat_mp quat, 
                  int nsites, int nmols);    /* Calculate torques from site forces  */
                                        /* Get h 2nd derivatives             */
void   escape(vec_mt *, int);
void   energy_dyad(mat_mt ke_dyad, mat_mt, double s, vec_mp moms, double mass,  
                   int nmols);                        /* Calculate mvv for stress term     */
double ke_cell(mat_mt hmom, real w);
void   force_calc(real **site, real **site_force, system_mt *system, 
                  spec_mt *species, real *chg, pot_mt *potpar, double *pe, 
                  mat_mt stress);                /* Calculate direct-space forces     */
double poteval(real *potpar, double r, int ptype, double chgsq);
double dist_pot(real *potpar, double cutoff, int ptype); 
                                      /* Returns integrated potential fn     */
void   ewald(real **site, real **site_force, system_mp system, spec_mt *species, 
             real *chg, double *pe, mat_mt stress); 
void   dump(system_mp system, spec_mt *species, vec_mt (*force), vec_mt (*torque), 
            mat_mt stress, double pe, restrt_mt *restart_header, 
            int backup_restart);        /*  Write dump data files            */
void   zero_real(real *r, int n);      /* Clear area of memory                     */
void   zero_double(double *r, int n);  /* Clear area of memory                     */
void   invert(mat_mt , mat_mt );        /* Matrix inverter                     */
double det(mat_mt );                        /* Returns matrix determinant             */
void   mat_vec_mul(mat_mt , vec_mp, vec_mp, int);
void   mat_mul(mat_mt a, mat_mt b, mat_mt c); /* 3 x 3 matrix multiplier     */
void   mat_add(mat_mt a, mat_mt b, mat_mt c); /* Add 2 3x3 matrices          */
void   mat_sca_mul(real s, mat_mt a, mat_mt b); 
void   mk_sigma(mat_mt h, mat_mt sigma);
void   mean_square(vec_mt (*x), real *meansq, int nmols);
                                             /* Calculates mean square of args    */
void   rdf_calc(real **site, system_mp system, spec_mt *species);
                                           /* Accumulate and bin rdf             */
double value(av_n type, int comp); /* Return thermodynamic average             */
double roll_av(av_n type, int comp); /* Return thermodynamic average             */
double vdot(int n, real *x, int ix, real *y, int iy); /* Fast  dot product   */
double sum(int n, real *x, int ix); /* Fast sum */
void   vscale(int n, double s, real *x, int ix); /* Vector by const multiply */
void   thermalise(system_mp system, spec_mt *species);              
                                      /* Randomize velocities to given temp  */
double trans_ke(mat_mt, vec_mt (*mom_s), real s, double mass, int nmols);
                                     /* Compute translational kinetic energy */
double rot_ke(quat_mt (*amom), real s, real *inertia, int nmols); 
                                    /* Compute rotational kinetic energy     */
void   q_conj_mul(quat_mp p, quat_mp q, quat_mp r, int n);   
                                      /* Quat. conjugated x by quat. dot     */
void   inhibit_vectorization(void);       /* Self-explanatory dummy          */
void   kernel(int, int, real *, double *, real *, real *, double, double, 
              double, int, real **);   /* Potential function evaluation      */
#ifdef SPMD
void   par_rsum(real *buf, int n);
void   par_dsum(double *buf, int n);
#endif 
gptr   *arralloc(size_mt,int,...);        /* Array allocator                      */
void   note(char *, ...);                /* Write a message to the output file */
/*========================== External data references ========================*/
extern contr_mt control;                    /* Main simulation control parms. */
extern int         ithread, nthreads;
/*========================== Macros ==========================================*/
#define ITER_MAX 10
#define        CONVRG        1.0e-7
/*  Can't rely on ANSI yet. */
#ifndef DBL_MIN
#   define DBL_MIN 1.0e-36
#endif
/*========================== Cache Parameters=================================*/
/* The default values are for the Cray T3D but are probably good enough 
 *  for most other systems too. */
#ifndef NCACHE
#   define NCACHE (256*sizeof(double)/sizeof(real))
#endif
#ifndef NLINE
#   define NLINE  (4*sizeof(double)/sizeof(real))
#endif
/*============================================================================*/
/******************************************************************************
 *   rescale    rescale velocities and quaternion derivatives to desired temp.*
 *   Exact behaviour is controlled by flag "control.scale_options".              *
 *   This is a bit flag with the following meanings:                              *
 *        bit 0:        scale temperature for each species separately.                      *
 *        bit 1:  scale rotational and translational velocities separately      *
 *      bit 2:        use rolling averages rather than instantaneous "temperatures" *
 *        bit 3:  don't scale at all, but re-initialize from MB distribution.   *
 ******************************************************************************/
void
rescale(system_mp system, spec_mp species)
{
   spec_mp        spec;
   int                ispec, imol, i;
   double         *temp_value = (double*)aalloc(2*system->nspecies,double);
   double        min_temp=MIN(value(t_n,0),roll_av(t_n,0));
   double        rtemp = 0.0, ttemp = 0.0, scale;
   double        total_mass;
   vec_mt        momentum;
   int                tdof=0, rdof=0;

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
         scale = sqrt(control.temp / temp_value[2*ispec])/system->ts;
         vscale(3 * spec->nmols,   scale, spec->mom[0], 1);
         if( spec->rdof > 0 )
         { 
            scale = sqrt(control.temp / temp_value[2*ispec+1]);
            vscale(4 * spec->nmols, scale, spec->amom[0], 1);
         }
         
      }
   }
   /*
    * Reset thermostat variable.  N.B need to recalculate H_0 in do_step()
    */
   system->ts = 1.0;
   system->tsmom = 0.0;
   /* 
    * Subtract spurious net velocity introduced by scaling species
    * separately.  This will introduce an apparent error into the
    * instantaneous "temperature".  But the error was there anyway
    * since net linear velocity shouldn't be included in the
    * calculation.  Do nothing for joint rescaling since that does
    * conserve momentum or for a framework system.
    */
   if( control.scale_options & 0x1 )
   {
      total_mass = 0.0;
      zero_real(momentum, 3);
      for (spec = species; spec < species+system->nspecies && ! spec->framework;
           spec++)
      {
         total_mass += spec->mass*spec->nmols;
         for(i = 0; i < 3; i++)        
            momentum[i] += sum(spec->nmols, spec->mom[0]+i,3);
      }
      if(spec == species+system->nspecies)/* Normal loop exit => no framework */
         for (spec = species; spec < species+system->nspecies; spec++)
            for(i = 0; i < 3; i++)            
               for(imol = 0; imol < spec->nmols; imol++)
                  spec->mom[imol][i] -= momentum[i] * spec->mass / total_mass;
   }
   tfree((gptr*)temp_value);
}
/******************************************************************************
 *  distant_const     Return the constant part of the distant-potential term  *
 *  c = - 2 pi sum i sum j Ni Nj Aij(cutoff), where i,j are site types, Ni,Nj *
 *  are their populations and Aij(r) is pot'l integrated from r to infinity.  *
 *  If Iflag == 0, return potential correction, == 1, pressure correction.    *
 ******************************************************************************/
static double 
distant_const(system_mp system, spec_mt *species, pot_mt *potpar, 
              double cutoff, int iflag)
{
   int             isite, id, jd;                /* Counters                      */
   spec_mp         spec;                               /* pointer to current species */
   int            *site_count = ialloc(system->max_id);        /* Numbers of each site
                                                         * type */
   double          c = 0.0;                               /* Accumulator for result     */
/*
 * Count the sites
 */
   memst(site_count,0,system->max_id*sizeof(int));
   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
NOVECTOR
      for (isite = 0; isite < spec->nsites; isite++)
      {
         inhibit_vectorization();
         site_count[spec->site_id[isite]] += spec->nmols;
      }
   }

   for (id = 1; id < system->max_id; id++)        /* Loops for sum over i,j     */
      for (jd = 1; jd < system->max_id; jd++)
      {
         c -= 2 * PI * site_count[id] * site_count[jd]
            * dist_pot(potpar[id + system->max_id*jd].p, cutoff, system->ptype);
         if( iflag )
            c += 2.0/3.0 * PI * site_count[id] * site_count[jd]
            * CUBE(cutoff) * poteval(potpar[id + system->max_id*jd].p, cutoff, 
                                     system->ptype, 0.0);
      }

   xfree(site_count);
   return (c);
}
/******************************************************************************
 * tot_ke().           Evaluate total kinetic energy.                              *
 ******************************************************************************/
double tot_ke (system_mt *sys, spec_mt *species, boolean rot_flag)
{
   double   ke = 0.0;
   spec_mt  *spec;
   int            ispec;
   for(spec=species, ispec = 0; ispec < sys->nspecies; spec++, ispec++)
   {
      ke += trans_ke(sys->h, spec->mom, sys->ts, spec->mass, spec->nmols) ;
      if(rot_flag && spec->rdof > 0)         /* Only if polyatomic species */
         ke += rot_ke(spec->amom, sys->ts, spec->inertia, spec->nmols);
   }
   return ke;
}
/******************************************************************************
 * stress_kin ().  Evaluate the kinetic part of the stress.                      *
 ******************************************************************************/
void stress_kin (mat_mt ke_dyad, system_mt *sys, spec_mt *species)
{
   spec_mt *spec;
   zero_real(ke_dyad[0], 9);
      
   for (spec = species; spec < &species[sys->nspecies]; spec++)
      energy_dyad(ke_dyad, sys->h, sys->ts, spec->mom, spec->mass, spec->nmols); 
}
/******************************************************************************
 * leapf_all_coords(). Update centre-of-mass and quaternion co-ordinates      *
 ******************************************************************************/
void leapf_all_coords(double step, system_mt *sys, spec_mt *species)
{
   spec_mt *spec;

   for (spec = species; spec < &species[sys->nspecies]; spec++)
   {
      leapf_com(step, spec->c_of_m, spec->mom, sys->h, sys->ts, 
                spec->mass, spec->nmols);
      if( spec->rdof > 0 )
         leapf_quat(step, spec->quat, spec->amom, 
                    spec->inertia, &sys->tsmom, sys->ts, spec->nmols);
   }
}
/******************************************************************************
 * leapf_all_momenta(). Update linear and angular momenta.                      *
 ******************************************************************************/
void leapf_all_momenta(double step, system_mt *sys, spec_mt *species,
                       vec_mp *force, vec_mp *torque)
{
   spec_mt *spec;
   int     ispec;

   for (ispec = 0, spec = species; ispec < sys->nspecies; ispec++, spec++)
   {
      leapf_mom(step*sys->ts, sys->h, spec->mom,  force[ispec], spec->nmols);
      if( spec->rdof > 0 )
         leapf_amom(step*sys->ts, spec->amom, torque[ispec], spec->nmols);
   }
}
/******************************************************************************
 * eval_forces().  This is the master potential and force evaluation routine. *
 *                 It calls  force() and  ewald to evaluate the site forces   *
 *                 converts these to the molecular forces and torques which   *
 *                 are passed out as arguments.  All of the site_co-ordinates *
 *                 and force arrays are local to this routine.                *
 *                 Site->molecular virial, durface dipole and distant         *
 *                 pressure corrections are also applied here.                      *
 ******************************************************************************/
void 
eval_forces(system_mp sys,             /* Pointer to system info        (in) */
            spec_mt *species,          /* Array of species info         (in) */
            site_mt *site_info,        /* Array of site info structures (in) */ 
            pot_mt *potpar,            /* Array of potential parameters (in) */
            double *pe,                /* Potential energy real/Ewald  (out) */
            real *dip_mom,             /* Total system dipole moment   (out) */
            mat_mt stress,               /* Virial part of stress        (out) */
            vec_mp *force,             /* Array of molecular forces    (out) */
            vec_mp *torque)            /* Array of molecular torques   (out) */
{
/*
 * The following declarations are arrays of pointers to the forces
 * etc for each species.  That is force[i] is a pointer to the force on
 * molecule 0 of species i
 */
   real                ***site_sp = (real***)arralloc((size_mt)sizeof(real*), 2,
                                               0, sys->nspecies-1, 0, 2),
                  ***force_sp = (real***)arralloc((size_mt)sizeof(real*), 2,
                                               0, sys->nspecies-1, 0, 2); 
/*
 * The following declarations are pointers to the force etc for the whole
 * system, and are set equal to (eg) force[0]
 */
   int                nsarray = ((sys->nsites - 1) | (NCACHE - 1)) + 1+NLINE;
   real                **site = (real**)arralloc((size_mt)sizeof(real), 2,
                                          0, 2, 0, nsarray-1);
   real                **site_force = (real**)arralloc((size_mt)sizeof(real), 2,
                                                0, 2, 0, nsarray-1);
   /*
    * Other local variables
    */
   real            *chg = dalloc(sys->nsites), *chg_ptr;
   vec_mp          c_of_m = ralloc(sys->nmols);
   spec_mp         spec, fspec /*Framework species */;
   int             nspecies = sys->nspecies;
   int             ispec, imol, isite;
   int             i, j;
   double          vol = det(sys->h);
   int                   nsitesxf, nmolsxf;  /* Count of non-framework sites, mols. */
   static double   dist, distp;
   static boolean  init = true;

   /*
    * Initialisation of distant potential constant - executed first time only
    */
   if (init)
   {
      dist  = distant_const(sys, species, potpar, control.cutoff,0);
      distp = distant_const(sys, species, potpar, control.cutoff,1);

      if( ithread == 0 )
         note("Distant potential correction = %f, Pressure correction = %f",
              CONV_E * dist / vol, CONV_P * distp / (vol * vol));
      init = false;
   }
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
   isite = 0;
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {
      for( i = 0; i < 3; i++ )
      {
         site_sp[ispec][i] = site[i]+isite;
         force_sp[ispec][i] = site_force[i]+isite;
      }
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
/*
 * Set some accumulators to zero
 */
   zero_real(stress[0], 9);               /* Initialise stress tensor   */
   zero_real(dip_mom, 3);
   zero_real(site_force[0], nsarray);
   zero_real(site_force[1], nsarray);
   zero_real(site_force[2], nsarray);
   zero_double(pe, NPE);
   zero_real(force[0][0], sys->nmols);
   zero_real(torque[0][0], sys->nmols_r);
   
/*
 * Calculate the site positions at this timestep - loop over species
 */
   for (spec = species; spec < &species[nspecies]; spec++)
   {
      make_sites(sys->h, spec->c_of_m, spec->quat, spec->p_f_sites,
                 site_sp[spec-species],spec->nmols,spec->nsites,
                 control.molpbc?MOLPBC:SITEPBC);
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
   }
/*
 *  Sum Pot, energies, forces and stress from each parallel invocation
 */
#ifdef SPMD
   par_dsum(pe, NPE);
   par_rsum(stress[0], 9);
   par_rsum(site_force[0], 3*nsarray);
#endif
   for (spec = species; spec < &species[nspecies]; spec++)
   {
      make_sites(sys->h, spec->c_of_m, spec->quat, spec->p_f_sites,
                 site_sp[spec-species],spec->nmols,spec->nsites,
                 spec->framework?SITEPBC:MOLPBC);
   }
   if (control.alpha > ALPHAMIN)
   {
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
   mat_vec_mul(sys->h, sys->c_of_m, c_of_m, sys->nmols);
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
                - vdot(nmolsxf, force[0][0] + i, 3, c_of_m[0] + j, 3);
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
      stress[i][i] += distp / vol;

/*
 * Deallocate the dynamic storage before exiting
 */
   afree((gptr*)site);
   afree((gptr*)site_force);
   afree((gptr*)site_sp);
   afree((gptr*)force_sp);
   xfree(chg);
   xfree(c_of_m);
}
/******************************************************************************
 *   do_step       This routine controls the main part of the calculation for *
 *                 each timestep.  This is where the actual integration       *
 *                 algorithm is implemented calling the support functions     *
 *                 from leapfrog.c and eval_forces().                              *
 *                 The periodic data dump is called from here while the forces*
 *                 and torques are in scope.                                  *
 ******************************************************************************/
void
do_step(system_mt *sys,                 /* Pointer to system info        (in) */
        spec_mt *species,                 /* Array of species info         (in) */
        site_mt *site_info,                 /* Array of site info structures (in) */
        pot_mt *potpar,                 /* Array of potential parameters (in) */
        vec_mt (*meansq_f_t)[2],        /* Mean square force and torque (out) */
        double *pe,                         /* Potential energy real/Ewald  (out) */
        real *dip_mom,                         /* Total system dipole moment   (out) */
        mat_mt stress_vir,                 /* Virial part of stress        (out) */
        restrt_mt *restart_header,      /* What the name says.           (in) */
        int backup_restart,                /* Flag signalling backup restart (in)*/
        int init_H_0)                        /* Flag re-init of H_0 required. (in) */
{
 /*
 * The following declarations are arrays of pointers to the forces
 * etc for each species.  That is force[i] is a pointer to the force on
 * molecule 0 of species i
 */
   vec_mp   *force  = palloc(sys->nspecies),
            *torque = palloc(sys->nspecies);
/*
 * The following declarations are pointers to the force etc for the whole
 * system, and are set equal to (eg) force[0]
 */
   vec_mp   force_base = ralloc(sys->nmols),
            torque_base = sys->nmols_r?ralloc(sys->nmols_r):0L;
   mat_mt          ke_dyad;
   double           ke, tke;
   spec_mp         spec;
   int             ispec, imol, imol_r;
   int                   nspecies = sys->nspecies;
   boolean           uni = control.const_pressure> 0 && (control.const_pressure%2==0);
   static           double saved_pe;
   /*
    * Initialize force and torque arays.
    */
   imol = 0; imol_r = 0;
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {
      force[ispec] = force_base+imol;
      torque[ispec] = torque_base+imol_r;
      imol += spec->nmols;
      if (spec->quat)
         imol_r += spec->nmols;
   }

   /*
    * Half step for H1. 
    */
   if( control.const_temp )
   {
#ifdef GL_THERM
      gleap_therm(0.5*control.step, control.ttmass, sys->d_of_f*kB*control.temp, 
                  &sys->ts, &sys->tsmom);
#else
      leapf_nose_therm(0.5*control.step, control.ttmass, &sys->ts, &sys->tsmom);
      sys->tsmom -= 0.5*control.step*sys->d_of_f*kB*control.temp*(1.0+log(sys->ts));
#endif
   }
   /*
    * H2 half step
    */
   if( control.const_pressure )
      gleap_cell(0.5*control.step, control.pmass, sys->ts, control.pressure, 
                 control.strain_mask, sys->h, sys->hmom, &sys->tsmom, uni);
   /*
    * H3 half step
    */
   /*
    * Partial update of thermostat momentum.  Note trans. KE only as
    * rotational contribution added inside leapf_quat()
    */
   if( control.const_temp )
     sys->tsmom += 0.5*control.step*tot_ke(sys, species, 0);
   if( control.const_pressure )
   {
      stress_kin(ke_dyad, sys, species);
      update_hmom(0.5*control.step, sys->ts, sys->h, ke_dyad, sys->hmom, uni);
   }
   leapf_all_coords(0.5*control.step, sys, species);
   /*
    * H4 full step
    */
   eval_forces(sys, species, site_info, potpar,
               pe, dip_mom, stress_vir, force, torque);
   saved_pe = pe[0] + pe[1];

   leapf_all_momenta(0.5*control.step, sys, species, force, torque);
   /*
    * Evaluate initial Hamiltonian H_0 at start of run or after a velocity
    * scaling step.  In this version we evaluate it at the half-step.
    */ 
   if(control.istep == 1 || init_H_0 )
   {
      ke = tot_ke(sys, species,1);
      sys->H_0 = ke + saved_pe + SQR(sys->tsmom)/(2.0*control.ttmass) 
                            + sys->d_of_f*kB*control.temp * log(sys->ts)
                            + ke_cell(sys->hmom, control.pmass);
       if( control.const_pressure )
         sys->H_0 += control.pressure*det(sys->h);
   }
#ifdef DEBUG_THERMOSTAT
   if( control.istep%control.print_interval == 0 && ithread == 0 )
   {
      double H, HP, HS, HHP, HHS;
      HS = sys->d_of_f*kB*control.temp * log(sys->ts);
      HP = SQR(sys->tsmom)/(2.0*control.ttmass);
      HHP = ke_cell(sys->hmom, control.pmass);
      HHS = control.pressure*det(sys->h);
      ke = tot_ke(sys, species,1);
      H = ke + pe[0] + pe[1];
      if( control.const_temp )
         H += HP + HS;
      if( control.const_pressure )
         H += HHP + HHS;
      fprintf(stderr,
             "do_step:  s= %7.4g          ps= %9.5g      H= %12.7g  HK= %12.7g  HV= %12.7g  \n          HP= %12.7g    HS= %12.7g   HHP= %12.7g  HHS= %12.7g   (H-H_0)s= %12.7g\n",
             sys->ts,sys->tsmom, H,ke, pe[0] + pe[1], HP, HS,HHP, HHS, 
              (H-sys->H_0)*sys->ts); 
   }
#endif
   leapf_all_momenta(0.5*control.step, sys, species, force, torque);

   if( control.const_temp )
      sys->tsmom -= control.step*(pe[0]+pe[1] - sys->H_0);

   if( control.const_pressure )
      update_hmom(control.step, sys->ts, sys->h, stress_vir, sys->hmom, uni);
   /*
    * Second H3 half step.
    */
   if( control.const_temp )
     sys->tsmom += 0.5*control.step*tot_ke(sys, species, 0);
   if( control.const_pressure )
   {
      stress_kin(ke_dyad, sys, species);
      update_hmom(0.5*control.step, sys->ts, sys->h, ke_dyad, sys->hmom, uni);
   }

   leapf_all_coords(0.5*control.step, sys, species);
   /*
    * Second H2 half step.
    */
   if( control.const_pressure )
      gleap_cell(0.5*control.step, control.pmass, sys->ts, control.pressure, 
                 control.strain_mask, sys->h, sys->hmom, &sys->tsmom, uni);
   /*
    * Final half step for H1.
    */
   if( control.const_temp )
#ifdef GL_THERM
      gleap_therm(0.5*control.step, control.ttmass, sys->d_of_f*kB*control.temp, 
                  &sys->ts, &sys->tsmom);
#else
      sys->tsmom -= 0.5*control.step*sys->d_of_f*kB*control.temp*(1.0+log(sys->ts));
      leapf_nose_therm(0.5*control.step, control.ttmass, &sys->ts, &sys->tsmom);
#endif

#ifdef DEBUG_THERMOSTAT
   if( control.istep%control.print_interval == 0 && ithread == 0 )
   {
      ke = tot_ke(sys, species,1);
      fprintf(stderr, "do_step2:  HK2= %12.7g\n", ke);
   }
#endif
/*
 *  Apply constraint to any framework molecules.
 */
   for (spec = species; spec < &species[nspecies]; spec++)
      if( spec->framework )
      {
         zero_real(spec->mom[0], 3*spec->nmols);
      }
/*
 * Calculate mean-square forces and torques
 */
   zero_real(meansq_f_t[0][0], 6 * sys->nspecies);
   for (ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {
      mean_square(force[ispec], meansq_f_t[ispec][0], spec->nmols);
      if (spec->rdof > 0)
         mean_square(torque[ispec], meansq_f_t[ispec][1], spec->nmols);
   }

/*
 * Perform periodic dump of dynamic data
 */
   if( ithread == 0 )
   {
      if (control.dump_interval > 0 && control.dump_level != 0 &&
         control.istep >= control.begin_dump &&
          (control.istep - control.begin_dump) % control.dump_interval == 0)
       dump(sys, species, force_base, torque_base, stress_vir, pe[0] + pe[1], 
            restart_header, backup_restart);
      
   }
   xfree(force);
   xfree(torque);
   xfree(force_base);
   if (torque_base != 0)
      xfree(torque_base);
}
