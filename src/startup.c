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
 * startup	Functions to perform operations to set up and initialise the  *
 *		simulation and control the reading of control parameters,     *
 *		system specification and restart files.  Contents:	      *
 * start_up()		Master start-up control function		      *
 * allocate_dynamics()	Set up dynamic variable arrays by memory allocation   *
 * initialise_sysdef()	Calculate 'whole-system' quantities and more	      *
 * check_sysdef()	Check new system specification is consistent with old *
 * interp()		Scale derivatives when changing timestep	      *
 * skew_start()		Set up skew-cyclic initial configuration	      *
 * random_start()	Set up random initial configuration (not used)	      *
 *   			(N.B) lattice_start in "input.c"		      *
 * thermalise()		Set up Maxwell-Boltzmann velocity distribution	      *
 * select_spec()	Randomly choose a species			      *
 * random_quat()	Generate quaternion  for uniform random rotation      *
 * gauss_rand()		Return random sample from univariant gaussian         *
 ******************************************************************************
 *      Revision Log
 *      $Log: startup.c,v $
 *      Revision 2.19  2000/05/23 15:23:08  keith
 *      First attempt at a thermostatted version of the Leapfrog code
 *      using either a Nose or a Nose-Poincare thermostat
 *
 *      Revision 2.18  2000/04/27 17:57:11  keith
 *      Converted to use full ANSI function prototypes
 *
 *      Revision 2.17  2000/04/26 16:01:02  keith
 *      Dullweber, Leimkuhler and McLachlan rotational leapfrog version.
 *
 *      Revision 2.16  1999/12/20 15:19:26  keith
 *      Check for rdf-limit or nbinds changed on restart, and handle
 *      gracefully.
 *
 *      Revision 2.15  1999/10/11 09:51:07  keith
 *      Fully implemented new constant-pressure algorithm.
 *      Select by "const-pressure=2" in control.
 *
 *
 *      Revision 2.14  1999/10/08 10:49:39  keith
 *      Added checks to behave sensibly if rdf-interval changed during accumulation
 *      of RDF data upon restart and to discard excess data if begin-rdf changed.
 *      Added "validate_control()" to perform some simple checks on (usually sign of)
 *      input parameters.
 *
 *      Revision 2.13  1998/05/07 17:06:11  keith
 *      Reworked all conditional compliation macros to be
 *      feature-specific rather than OS specific.
 *      This is for use with GNU autoconf.
 *
 *      Revision 2.12  1996/08/15 14:35:39  keith
 *      Fixed restart structure correctly - broken in prev version.
 *      Thermostat parameters may not be properly read.
 *
 *      Revision 2.11  1996/01/15 15:23:30  keith
 *      Changed input to new units
 *      Reset defaults for cutoff, alpha to zero for "auto" select.
 *      (Murashov's code missed this)
 *      Added safety checks to auto-ewald routine.
 *
 * Revision 2.10  1995/12/04 11:45:49  keith
 * Nose-Hoover and Gaussian (Hoover constrained) thermostats added.
 * Thanks to V. Murashov.
 *
 * Revision 2.9  1995/01/03  13:56:39  keith
 * Added auto-generation of Ewald Sum parameters using Fincham's formulae.
 *
 * Revision 2.8  1994/07/07  16:57:01  keith
 * Updated for parallel execution on SPMD machines.
 * Interface to MP library routines hidden by par_*() calls.
 * Compile with -DSPMD to activate
 *
 * Revision 2.7  1994/06/08  13:16:34  keith
 * Changed all timestep-related parameters to type "long". This means
 * that 16-bit DOS compilers can do more than 32767 timesteps.
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
 * "dump".
 *
 * Declared as "static"  all functions which should be.
 *
 * Added const qualifier to (re-)declarations of ANSI library
 * emulation routines to give reliable compilation even
 * without ANSI_LIBS macro. (#define's away for K&R
 * compilers)
 *
 *
 * Moved declaration of "match" structure from input.c
 * Moved a few sanity tests & modifications of "control"
 * members to here from other modules.
 *
 * Revision 2.5  94/01/24  18:20:11  keith
 * Null checkin for release compatibility.
 * 
 * Revision 2.4  94/01/24  18:18:35  keith
 * Deleted commented-out code for NR jacobi() function. Eigens() is
 * now well-tested.  Eliminated compiler warnings.
 * 
 * Added external "backup_restart" to flag that run started from a
 * backup file.  For dump().
 * 
 * 
 * Revision 2.3  93/10/28  10:28:13  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:49:22  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.6.1.27  93/03/13  01:45:58  keith
 * Replaced NR "jacobi.c" with "eigens.c" for public release.
 * 
 * Revision 1.6.1.26  93/03/12  12:14:28  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.6.1.26  93/03/09  15:59:15  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.6.1.25  92/10/28  14:09:53  keith
 * Changed "site_[tp]" typedefs to avoid name clash on HP.
 * 
 * Revision 1.6.1.24  92/09/22  14:48:08  keith
 * Tidied up calls to improve "lint" rating.
 * 
 * Revision 1.6.1.23  92/06/26  17:03:28  keith
 * Got rid of assumption that memory returned by talloc() or
 * arralloc() is zeroed.  This enhances ANSI compatibility.
 * Removed memory zeroing from alloc.c() in consequence.
 * 
 * Revision 1.6.1.22  92/06/16  12:42:49  keith
 * Zeroed unit cell velocities/accns if constraint changed.
 * 
 * Revision 1.6.1.21  92/06/11  21:39:58  keith
 * Fixed bug which meant you couldn't change backup file name in restart if backup existed.
 * Added file locking against multiple runs using same dump or backup files.
 * 
 * Revision 1.6.1.20  92/04/21  17:50:08  keith
 * Fixed bug which compared char to NULL instead of 0.
 * 
 * Revision 1.6.1.19  92/03/19  15:45:54  keith
 * Added support for dynamic allocation of rolling average arrays,
 * conversion of existing restart files is done on fly.
 * 
 * Revision 1.6.1.18  92/03/11  12:56:11  keith
 * Changed "scale-separately" parameter to "scale options"
 * 
 * Revision 1.6.1.17  91/08/16  15:26:09  keith
 * Checked error returns from fread, fwrite, fseek and fclose more
 * rigourously.   Called strerror() to report errors.
 * 
 * Revision 1.6.1.16  91/08/15  18:12:16  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.6.1.15  91/03/12  15:43:19  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.6.1.14  91/02/04  18:20:08  keith
 * Call to inhibit vectorization() added on line 509 for titan.
 * Alloas use of vector_c option which otherwise caused crash.
 * 
 * Revision 1.6.1.13  90/10/25  19:00:26  keith
 * Made interpolate_derivatives() robust against interpolation from zero.
 * 
 * Revision 1.6.1.12  90/10/23  20:13:20  keith
 * Added dummy function call to inhibit vectorization.
 * This allows use of 'ivdep' compiler options and also
 * works round certain bugs in cray's scc compiler.
 * 
 * Revision 1.6.1.11  90/08/24  17:46:42  keith
 * Made 'reset-averages' parameter actually do something.
 * 
 * Revision 1.6.1.10  90/05/16  18:40:46  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.6.1.9  90/05/16  14:20:28  keith
 * *** empty log message ***
 * 
 * Revision 1.6.1.8  90/02/22  18:03:40  keith
 * Modified backup-restart code so as not to copy whole backup header into
 * restart header.  This erroneously wrote a timestamp to dumps after restrart.
 * 
 * Revision 1.6.1.7  89/11/21  16:32:38  keith
 * Removed member out_file from control and all uses. (Now command parameter).
 * Added new argument to specify output file.
 * Added default_control() to initialise control struct.
 * Fixed bug in start_up which tried to free unallocated var qpf.
 * 
 * Revision 1.6.1.6  89/11/01  17:06:24  keith
 * Made allocate_dynamics() externally visible to allow link with 'mdshak'.
 * 
 * Revision 1.6.1.5  89/10/05  18:18:02  keith
 * Modified initialise_sysdef() to stop rotation to princ frame of framework.
 * 
 * Revision 1.6.1.4  89/09/21  15:02:57  keith
 * Test for old timestep of zero when altering timestep, in which case do nothing.
 * 
 * Revision 1.6.1.3  89/09/19  14:50:15  keith
 * Fixed bug in check_sysdef which didn't increment species pointer.
 * 
 * Revision 1.6.1.2  89/09/05  10:28:10  keith
 * Fixed bug in code to disable Ewald sum for uncharged system (init.._sysdef()
 * Added calculation of species total charge.
 * 
 * Revision 1.6.1.1  89/08/31  10:14:40  keith
 * Mods to simulate framework structure
 * 
 * Revision 1.6  89/08/30  12:36:26  keith
 * Modified start_up to fix bug which only considered rotations
 * of one species.  In conjunction with change in lattice_start().
 * 
 * Revision 1.5  89/07/06  16:23:56  keith
 * Eliminated 'dump-offset' - renumbering starts at 1 for new dump run.
 * 
 * Revision 1.4  89/06/23  15:35:10  keith
 * print-control option deleted.
 * 
 * Revision 1.3  89/06/22  15:45:19  keith
 * Tidied up loops over species to use one pointer as counter.
 * 
 * Revision 1.2  89/06/14  17:57:48  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.1  89/04/27  15:16:41  keith
 * Initial revision
 * 
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/minphys2/keith/CVS/moldy/src/startup.c,v 2.19 2000/05/23 15:23:08 keith Exp $";
#endif
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<math.h>
#include 	"string.h"
#include	"stddef.h"
#include	"stdlib.h"
#include	<stdio.h>
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External function declarations ==================*/
gptr            *talloc(int n, size_mt size, int line, char *file);	       /* Interface to memory allocator       */
void            tfree(gptr *p);	       /* Free allocated memory	      	      */
void		read_control(FILE *file, const match_mt *match);
void		read_sysdef(FILE *file, system_mp system, spec_mp *spec_pp, site_mp *site_info, pot_mp *pot_ptr);
void		lattice_start(FILE *file, system_mp system, spec_mp species, quat_mt (*qpf));
void		conv_control(const unit_mt *unit, boolean direction);
void		conv_potentials(const unit_mt *unit_from, const unit_mt *unit_to, pot_mt *potpar, int npotpar, int ptype, site_mt *site_info, int max_id);
void		init_averages(int nspecies, char *vsn, long int roll_interval, long int old_roll_interval, int *av_convert);
void		convert_averages(long int roll_interval, long int old_roll_interval, int av_convert);
void		init_rdf(system_mp system);
void		re_re_header(FILE *restart, restrt_mt *header, contr_mt *contr);
void		re_re_sysdef(FILE *restart, char *vsn, system_mp system, spec_mp *spec_ptr, site_mp *site_info, pot_mp *pot_ptr);
void		read_restart(FILE *restart, char *vsn, system_mp system, int av_convert);
void		banner_page(system_mp system, spec_mt *species, restrt_mt *restart_header);
void		zero_real(real *r, int n);
void		eigens(real *A, real *RR, real *E, int N);
void		transpose(real (*a)[3], real (*b)[3]);
void		mat_vec_mul(real (*m)[3], vec_mp in_vec, vec_mp out_vec, int number);
double          det(real (*a)[3]);
void		q_mul(quat_mp p, quat_mp q, quat_mp r, int n);
void		rot_to_q(real (*rot)[3], real *quat);
char		*atime(void);
double		mdrand(void);
double		precision(void);
void		smdrand(long unsigned int seed);
void		inhibit_vectorization(void);	/* Self-explanatory dummy     */
void		note(char *, ...);	/* Write a message to the output file */
void		message(int *, ...);	/* Write a warning or error message   */
/*========================== External data references ========================*/
extern	contr_mt	control;       /* Main simulation control parms. */
extern int 		ithread, nthreads;
/*========================== GLOBAL variables ================================*/
const   unit_mt	prog_unit = {MUNIT, LUNIT, TUNIT, QUNIT};
	unit_mt	input_unit;		/* Unit specification (see Convert.c) */
static	char		backup_lockname[L_name];
static	char		dump_lockname[L_name];
#ifdef	DEBUG
static	char	afmt[] = "    %8s = %8X %8s = %8X %8s = %8X %8s = %8X\
 %8s = %8X %8s = %8X\n";
#endif
/*
 *  Default backup and temporary file names if not set in "defs.h"
 */
#ifndef BACKUP_FILE
#define BACKUP_FILE	"MDBACKUP"
#endif
#ifndef TEMP_FILE
#define TEMP_FILE	"MDTEMPX"
#endif
/*========================== Control file keyword template ===================*/
/*
 * format SFORM is defined as %NAMLENs in structs.h, to avoid overflow.
 */
const match_mt	match[] = {
{"title",            SFORM,  "Test Simulation",(gptr*) control.title},
{"nsteps",           "%ld",  "0",            (gptr*)&control.nsteps},
{"step",             "%lf",  "0.005",        (gptr*)&control.step},
{"text-mode-save",   "%d",   "0",            (gptr*)&control.print_sysdef},
{"new-sys-spec",     "%d",   "0",            (gptr*)&control.new_sysdef},
{"scale-options"   , "%d",   "0",            (gptr*)&control.scale_options},
{"therm-options"   , "%d",   "0",            (gptr*)&control.scale_options},
{"surface-dipole",   "%d",   "0",            (gptr*)&control.surface_dipole},
{"lattice-start",    "%d",   "0",            (gptr*)&control.lattice_start},
{"sys-spec-file",    SFORM,  "",             (gptr*)control.sysdef},
{"restart-file",     SFORM,  "",             (gptr*)control.restart_file},
{"save-file",        SFORM,  "",             (gptr*)control.save_file},
{"dump-file",        SFORM,  "",             (gptr*)control.dump_file},
{"backup-file",      SFORM,  BACKUP_FILE,    (gptr*)control.backup_file},
{"temp-file",        SFORM,  TEMP_FILE,      (gptr*)control.temp_file},
{"strict-cutoff",    "%d",   "0",            (gptr*)&control.strict_cutoff},
{"xdr",    	     "%d",   "1",            (gptr*)&control.xdr_write},
{"strain-mask",	     "%d",   "200",	     (gptr*)&control.strain_mask},
{"nbins",            "%d",   "100",          (gptr*)&control.nbins},
{"seed",             "%ld",  "1234567",      (gptr*)&control.seed},
{"page-width",       "%d",   "132",          (gptr*)&control.page_width},
{"page-length",      "%d",   "44",           (gptr*)&control.page_length},
{"scale-interval",   "%ld",  "10",           (gptr*)&control.scale_interval},
{"const-pressure",   "%d",   "0",            (gptr*)&control.const_pressure},
{"const-temp",       "%d",   "0",            (gptr*)&control.const_temp},
{"reset-averages",   "%d",   "0",            (gptr*)&control.reset_averages},
{"scale-end",        "%ld",  "1000000",      (gptr*)&control.scale_end},
{"begin-average",    "%ld",  "1001",         (gptr*)&control.begin_average},
{"average-interval", "%ld",  "5000",         (gptr*)&control.average_interval},
{"begin-dump",       "%ld",  "1",            (gptr*)&control.begin_dump},
{"dump-interval",    "%ld",  "20",           (gptr*)&control.dump_interval},
{"dump-level",       "%d",   "0",            (gptr*)&control.dump_level},
{"ndumps",           "%d",   "250",          (gptr*)&control.maxdumps},
{"backup-interval",  "%ld",  "500",          (gptr*)&control.backup_interval},
{"roll-interval",    "%ld",  "10",           (gptr*)&control.roll_interval},
{"print-interval",   "%ld",  "10",           (gptr*)&control.print_interval},
{"begin-rdf",        "%ld",  "1000000",      (gptr*)&control.begin_rdf},
{"rdf-interval",     "%ld",  "20",           (gptr*)&control.rdf_interval},
{"rdf-out",          "%ld",  "5000",         (gptr*)&control.rdf_out},
{"temperature",      "%lf",  "0.0",          (gptr*)&control.temp},
{"pressure",         "%lf",  "0.0",          (gptr*)&control.pressure},
{"w",                "%lf",  "100.0",        (gptr*)&control.pmass},
{"rtmass",           "%lf",  "100.0",        (gptr*)&control.rtmass},
{"ttmass",           "%lf",  "100.0",        (gptr*)&control.ttmass},
{"cutoff",           "%lf",  "0.0",          (gptr*)&control.cutoff},
{"subcell",          "%lf",  "0.0",          (gptr*)&control.subcell},
{"density",          "%lf",  "1.0",          (gptr*)&control.density},
{"alpha",            "%lf",  "0.0",          (gptr*)&control.alpha},
{"k-cutoff",         "%lf",  "0.0",          (gptr*)&control.k_cutoff},
{"rdf-limit",        "%lf",  "10.0",         (gptr*)&control.limit},
{"cpu-limit",        "%lf",  "1.0e20",       (gptr*)&control.cpu_limit},
{"mass-unit",        "%lf",  AMUSTR,         (gptr*)&input_unit.m},
{"length-unit",      "%lf",  "1.0e-10",      (gptr*)&input_unit.l},
{"time-unit",        "%lf",  "1.0e-13",      (gptr*)&input_unit.t},
{"charge-unit",      "%lf",  ECSTR,          (gptr*)&input_unit.q},
{0,0,0,0}	     }; 		/* Null termination essential.	*/

/*============================================================================*
 * Functions for external access to globals.		      *
 *============================================================================*/
void	rmlockfiles(void)
{
   if( backup_lockname[0] )
      (void)remove(backup_lockname);
   if( dump_lockname[0] )
      (void)remove(dump_lockname);
}
/*============================================================================*/
/******************************************************************************
 *  default_control.   Initialise 'control' with default values from 'match' *
 ******************************************************************************/
static
void	default_control(void)
{
   const match_mt	*match_p;
   char	tmp[64];

   for( match_p = match; match_p->key; match_p++)
      (void)sscanf(strncpy(tmp,match_p->defalt, sizeof tmp),
		   match_p->format, match_p->ptr);
}
/******************************************************************************
 * gauss_rand. Return a random variable from a gaussian distribution with uni *
 * variance.  After Press, Plannery, Teulkolsk & Vetterling p202.             *
 ******************************************************************************/
static
double	gauss_rand(void)
{
   static int		set = 1;
   static double	gset;
   double		v1, v2, r, fac;
   if(set)
   {
      do
      {
         v1 = 2.0*mdrand()-1.0;  v2 = 2.0*mdrand() - 1.0;
         r = v1*v1 + v2*v2;
      } while(r > 1.0);
      fac = sqrt(-2.0 * log(r) / r);
      gset = v1 * fac;
      set = 0;
      return(v2*fac);
   }
   else
   {
      set = 1;
      return(gset);
   }
}
/******************************************************************************
 * random_quat.   Generate quaternions to represent uniform random rotations. *
 * If q=(cos(psi/2), l sin(psi/2)) , l is a unit vector, uniformity requires  *
 * that q(0) be distributed linearly ie p(q(0)) = q(0). L is a random unit    *
 * vector, generated from spherical co-ordinates theta and phi with           *
 * distribution p(theta) = sin(theta) 					      *
 ******************************************************************************/
static
void	random_quat(quat_mp q, int n)
       	  				/* First quaternion		(out) */
   	  				/* Number to be generated.       (in) */
{
   double	phi, cos_theta, sin_theta, st2;
   while(n-- > 0)
   {
      phi = 2.0*PI*mdrand();		/* Phi is uniform on [0, 2pi)	      */
      cos_theta = 1.0 - 2.0*mdrand();	/* 0 <= theta < pi, p(theta)=sin()    */
      sin_theta = sqrt(1.0 - SQR(cos_theta));
      (*q)[0] = sqrt(mdrand());
      st2 = sqrt(1.0 - SQR((*q)[0]));
      (*q)[1] = st2*sin_theta*sin(phi);
      (*q)[2] = st2*sin_theta*cos(phi);
      (*q)[3] = st2*cos_theta;
      q++;
   }
}
/******************************************************************************
 *  select_spec  Choose a species at random, weighted by proportion of whole  *
 ******************************************************************************/
static
spec_mp	select_spec(system_mp system, spec_mt *species)
{
   int		sel = mdrand() * system->nmols;
   spec_mp	spec;

   for (spec = species; spec < species+system->nspecies; spec++)
   {
      if(sel < spec->nmols)
         return(spec);
      sel -= spec->nmols;
   }
   return((spec_mp)-1);
}
/******************************************************************************
 *  skew_start.  Make a starting configuration with all molecules arranged at *
 *  regular intervals on a line passing at an angle through the MD cell.  This*
 *  permits a reasonable spacing between molecules without restricting the    *
 *  number allowed as in a lattice start.                                     *
 ******************************************************************************/
static
void	skew_start(system_mp system, spec_mt *species)
{
   int		ispec, imol;		/* Counters for species, molecules etc*/
   spec_mp	spec;
   double	mass = 0.0;		/* Whole system mass		      */
   double	n_third = pow((double)system->nmols,1.0/3.0);
   int		nz = 1, ny = (int)(n_third+0.5), nx = (int)(SQR(n_third)+0.5);
   int		*nmols = ialloc(system->nspecies), nm;
   double	delta_x = (double)nx / system->nmols,
                delta_y = (double)ny / system->nmols,
                delta_z = (double)nz / system->nmols;

   for (spec = species; spec < species+system->nspecies; spec++)
   {
      inhibit_vectorization();      /* Circumvent cray scc bug  */
      mass += spec->mass * spec->nmols;
   }
      
   system->h[0][0] = system->h[1][1] = system->h[2][2] 
                   = pow(mass/control.density, 1.0/3.0);
 					/* L = cube root of mass/density      */
   memst(nmols, 0, system->nspecies*sizeof(int));
   for(imol = 0; imol < system->nmols; imol++)
   {
      do
      {
         spec = select_spec(system, species);		/* Choose species     */
	 ispec = spec-species;
      }
      while(nmols[ispec] >= spec->nmols);	/* Repeat if all set  */

      nm = nmols[ispec];
      spec->c_of_m[nm][0] = imol*delta_x;
      spec->c_of_m[nm][1] = imol*delta_y;
      spec->c_of_m[nm][2] = imol*delta_z;
      nmols[ispec]++;
   }

   random_quat(system->quat, system->nmols_r);
   xfree(nmols);
}
/******************************************************************************
 * random_start.  This function generates a completely random starting        *
 * configuration.  Molecules are placed at random locations and orientations  *
 * in md box, whose size is chosed to give the required density.	      *
 ******************************************************************************/
#ifdef NOT_FOR_NOW
void	random_start(system, species)
system_mp	system;
spec_mt	species[];
{
   int		imol, i;		/* Counters for species, molecules etc*/
   double	mass = 0.0;		/* Whole system mass		      */

   for (spec = species; spec < species+system->nspecies; spec++)
      mass += spec->mass * spec->nmols;
      
   system->h[0][0] = system->h[1][1] = system->h[2][2] 
                   = pow(mass/control.density, 1.0/3.0);
    					/* L = cube root of mass/density      */
   for(imol = 0; imol < system->nmols; imol++)
      for(i = 0; i < 3; i++)		/* Centre of mass co-ords -1 < x < 1  */
         system->c_of_m[imol][i] = mdrand() - 1.0;
   random_quat(system->quat, system->nmols_r);
}
#endif
/******************************************************************************
 *  thermalise  set velocities and quaternion derivatives to values sampled   *
 *  at random from the Boltzmann distribution for the required temperature.   *
 ******************************************************************************/
void	thermalise(system_mp system, spec_mt *species)
{
   int		imol, i;		/* Counters for species, molecules etc*/
   spec_mp	spec;			/* Pointer to species[ispec]	      */
   double	omega_sq;		/* |omega|squared / 4		      */
   double	root_ktm, root_kti[3];	/* Gaussian widths of MB distribution */
   double	total_mass = 0;
   vec_mt	momentum;	      	/* Whole system momentum	      */

   zero_real(momentum, 3);
   /*
    *  Set accelerations to zero.
    */
   zero_real(system->vel[0],   3*system->nmols);
   system->ts = 1.0;
   system->rs = 1.0;
   system->tsmom = 0.0;
   system->rsmom = 0.0;
   
   for (spec = species; spec < species+system->nspecies; spec++)
   {
      if( !spec->framework)
      {
	 root_ktm = sqrt(kB * control.temp / spec->mass) / system->h[0][0];
	 total_mass += spec->mass*spec->nmols;
	 for(imol = 0; imol < spec->nmols; imol++)
	    for(i = 0; i < 3; i++)	/* Centre of mass co-ords -1 < x < 1  */
	    {
	       spec->vel[imol][i]    = system->ts* root_ktm * gauss_rand();
	       momentum[i] += spec->mass*spec->vel[imol][i];
	    }
	 
	 if(spec->rdof > 0)
	 {
	    for(i = 0; i < 3; i++)
	       if(spec->inertia[i] != 0.0)
		  root_kti[i] = sqrt(kB * control.temp / spec->inertia[i]);
	       else
		  root_kti[i] = 0.0;
	    
	    for(imol = 0; imol < spec->nmols; imol++)
	    {
	       spec->avel[imol][0] = 0.0;
	       for(i = 0; i < 3; i++)	/* Centre of mass co-ords -1 < x < 1  */
		  spec->avel[imol][i+1] = root_kti[i] * gauss_rand();
	    }
	 }
      }
   }
   for (spec = species; spec < species+system->nspecies; spec++)
      if( !spec->framework)
      {
	 for(i = 0; i < 3; i++)
	    
	    for(imol = 0; imol < spec->nmols; imol++)
	       spec->vel[imol][i] -= momentum[i] / (system->ts * total_mass);
      }
}
/******************************************************************************
 *  initialise_sysdef    Computes various quanities and completes the set-up  *
 *  of the system once it has been read in.  Quantities computed are :-       *
 *  Total numbers of molecules, sites, degrees of freedom in whole system,    *
 *  mass and inertia tensor for each species.  The molecular site co-         *
 *  ordinates are re-expressed wrt the molecular centre of mass, the inertia  *
 *  tensor is diagonalised and the site co-ordinates rotated to the principal *
 *  frame.								      *
 ******************************************************************************/
#define LTR(i,j) (((i)*(i)+(i))/2+(j))
void	initialise_sysdef(system_mp system, spec_mt *species, site_mt *site_info, quat_mt (*qpf))
         	       
       		          
       		            
       		      			/* Quaternion rotation to princ.frame*/
{
   vec_mt	c_of_m;			/* Co-ordinates of centre of mass    */
   vec_mt	dipole;			/* Molecular dipole moment	     */
   real 	inertia[6];		/* Inertia tensor		     */
   mat_mt	v;			/* Transformation matrix to prin. fr.*/
   spec_mp	spec;			/* Used for looping over species     */
   double	mass;			/* Temporary for site mass	     */
   int		nz;			/* Count of zero moments of inertia  */
   int		i, j, isite, id; 	/* Various loop counters	     */
   boolean	flag;			/* Used to test for charges	     */
   double	imax;			/* Largest moment of inertia	     */
   double	eps = 10.0*precision(); /* Criterion for "zero" moment.	     */

   system->nsites  = 0;  system->nmols  = 0;
   system->nmols_r = 0;  system->d_of_f = 0;

   for (spec = species; spec < species+system->nspecies; spec++)
   {					/* Loop over molecular species       */
      system->nmols  += spec->nmols;
      system->nsites += spec->nmols * spec->nsites;
      if( !spec->framework )
	 system->d_of_f += spec->nmols * 3;

      zero_real(c_of_m,3);		/* Initialise C_of_M for this species*/
      zero_real(dipole,3);		/* And dipole moment		     */
      spec->mass = 0.0;			/* Mass				     */
      spec->charge = 0.0;		/* And total charge		     */
      for(isite=0; isite < spec->nsites; isite++) /* Calculate (sum m*r) and */
      {					/* molecular mass.		     */
         for(i=0; i<3; i++)
            c_of_m[i] += spec->p_f_sites[isite][i] 
                         * site_info[spec->site_id[isite]].mass;
         spec->mass += site_info[spec->site_id[isite]].mass;
	 spec->charge += site_info[spec->site_id[isite]].charge;
      }

      if(spec->mass < 1.0)		/* Lighter than 1 amu ?              */
         message(NULLI,NULLP,FATAL,ZMASS,spec->name,spec->mass);

      for(i=0; i < 3; i++)		/* Finish calculation of c. of mass. */
         c_of_m[i] /= spec->mass;

      for(isite=0; isite < spec->nsites; isite++)
         for(i=0; i < 3; i++)		/* Subtract c_of_m from co-ordinates */
            spec->p_f_sites[isite][i] -= c_of_m[i];

      if(spec->nsites > 1)		/* If molecule is polyatomic         */
      {
         zero_real(inertia,6);	/* Initialise inertia tensor	     */
         for(isite=0; isite < spec->nsites; isite++)
         {
            mass = site_info[spec->site_id[isite]].mass;
            for(i=0; i < 3; i++)	/* Calculate inertia tensor	     */
            {
               inertia[LTR(i,i)] += mass * SUMSQ(spec->p_f_sites[isite]);
               for(j=0; j <= i; j++)
         	  inertia[LTR(i,j)] -= mass * spec->p_f_sites[isite][i] 
         	                            * spec->p_f_sites[isite][j];
            }
         }
#ifdef	DEBUG
         printf(" *D* Molecule type %d, mass = %g, C of M = (%g,%g,%g)\n",
                spec-species, spec->mass, c_of_m[0], c_of_m[1], c_of_m[2]);
         print_mat(inertia, " *D* Inertia Tensor");
#endif
	 eigens(inertia,v[0],spec->inertia,3);
	 /*	 eigensort(v[0], spec->inertia, 3);*/
	 rot_to_q(v, qpf[spec-species]);	/* make equivalent quaternion*/
#ifdef	DEBUG
         print_mat(v," *D* Rotation Mat.");
#endif
	 imax = MAX3(spec->inertia[0],spec->inertia[1],spec->inertia[2]);
	 nz = 0;
	 for( i=0; i<3; i++)			 /* Count zero  moments.     */ 
	    if( spec->inertia[i] < eps*imax )
	    {
	       nz++;
	       spec->inertia[i] = 0.0;
	    }
	    
         spec->rdof = 3-nz;			/* Rotational deg. of freedom*/
	 if( spec->framework )			/* Frameworks can't rotate   */
	 {
	    spec->rdof = 0;
	    spec->quat = 0;
	 }
         if(spec->rdof > 0)			/* Count molecules with      */
	 {
            system->nmols_r += spec->nmols;     /* rotational freedom.       */
	    mat_vec_mul(v,spec->p_f_sites, spec->p_f_sites, spec->nsites);
	 }
         system->d_of_f += spec->rdof * spec->nmols;/* Count total d of f    */

	 for(i = 0; i < 3; i++)			/* Calculate molecular dipole*/
	    for(isite = 0; isite < spec->nsites; isite++)
	       dipole[i] += spec->p_f_sites[isite][i] * 
	       		    site_info[spec->site_id[isite]].charge;
	 spec->dipole = sqrt(SUMSQ(dipole));
      }
      else
      {
	 spec->dipole = 0;
	 spec->rdof = 0;
	 zero_real(spec->inertia,3);
	 spec->quat = 0;
      }
   }
#ifdef	DEBUG
   printf(" *D* Totals: nsites = %d, nmols = %d, nmols_r = %d, dof = %d\n",
          system->nsites, system->nmols, system->nmols_r, system->d_of_f);
#endif
   
   flag = false;			/* Test to see if any charges present */
   for(id = 1; id < system->max_id; id++)
      flag |= (site_info[id].charge != 0.0);
   if(!flag) control.alpha = -1.0;	/* Don't call Ewald sum if not	      */
}
/******************************************************************************
 *  allocate_dynamics	 Allocate memory for the dynamic MD variables         *
 ******************************************************************************/
void	allocate_dynamics(system_mp system, spec_mt *species)
{
   spec_mp	spec;			/* Alias for species[ispec]           */
   int		nmol_cum = 0,		/* Cumulative number of molecules     */
   		nmolr_cum = 0;		/* As above excluding point atoms     */

   system->c_of_m = ralloc(system->nmols);
   system->vel    = ralloc(system->nmols);
   system->velp   = ralloc(system->nmols);
#ifdef	DEBUG
   printf(" *D* System Dynamic variables (all %d x 3 reals)\n",system->nmols);
   printf(afmt,"c_of_m",system->c_of_m,"vel",system->vel,"velp",system->velp,
          "acc",system->acc,"acco",system->acco,"accvo",system->accvo);
#endif

   if(system->nmols_r > 0)
   {
      system->quat   = qalloc(system->nmols_r);
      system->avel   = qalloc(system->nmols_r);
      system->avelp  = qalloc(system->nmols_r);
   }
   else
      system->quat = system->avel = system->avelp = 0;


   system->h       = ralloc(3);
   system->hdot    = ralloc(3);
   system->hdotp   = ralloc(3);
   zero_real(system->h[0],       9);
   zero_real(system->hdot[0],    9);
   zero_real(system->hdotp[0],   9);
#ifdef	DEBUG
   printf(" *D* System Dynamic variables (all 9 reals)\n");
   printf(afmt,"h",system->h, "hdot",system->hdot, "hdotp",system->hdotp,
       "hddot",system->hddot,"hddoto",system->hddoto,"hddotvo",system->hddotvo);
#endif

   for (spec = species; spec < species+system->nspecies; spec++)
   {
      inhibit_vectorization();
      spec->c_of_m = system->c_of_m + nmol_cum;
      spec->vel    = system->vel    + nmol_cum;
      spec->velp   = system->velp   + nmol_cum;
#ifdef	DEBUG
      printf(" *D* Species %d Dynamic variables (all %d x 3 reals)\n",
	     spec-species, spec->nmols);
      printf(afmt,"c_of_m",spec->c_of_m,"vel",spec->vel,"velp",spec->velp,
          "acc",spec->acc,"acco",spec->acco,"accvo",spec->accvo);
#endif
      if(spec->rdof > 0)
      {
         spec->quat    = system->quat    + nmolr_cum;
         spec->avel    = system->avel    + nmolr_cum;
         spec->avelp   = system->avelp   + nmolr_cum;
         nmolr_cum += spec->nmols;
#ifdef	DEBUG
         printf(" *D* Species %d Dynamic variables (all %d x 4 reals)\n",
		spec-species, spec->nmols);
         printf(afmt,"quat",spec->quat, "avel",spec->avel, "avelp",spec->avelp,
             "qddot",spec->qddot,"qddoto",spec->qddoto,"qddotvo",spec->qddotvo);
#endif
      }
      else
	 spec->quat = spec->avel = spec->avelp = 0;
      nmol_cum += spec->nmols;
   }
   /*
    * These may not be initialised if reading an old restart file,
    * so zero them here for safety.
    */
   system->ts = 1.0;
   system->rs = 1.0;
   system->tsmom = 0.0;
   system->rsmom = 0.0;
}
/******************************************************************************
 *  Interpolate_derivatives & interp.    Interp is a quadratic interpolation  *
 *  routine for scaling derivatives for a new timestep.                       *
 *  Interpolate_derivatives calls it for all the dynamic variables.	      *
 ******************************************************************************/
static
void	interp(double ratio, real *x, real *xo, real *xvo, int n)
      	      				/* Between old and new timesteps      */
    	              			/* Pointers to 1st dynamic variable   */
   	  				/* Size of arrays x, xo, xvo          */
{
   double	c1 = 1.0 - 1.5*ratio + 0.5*SQR(ratio),
   		c2 = 2.0*ratio - SQR(ratio),
   		c3 = -0.5*ratio + 0.5*SQR(ratio),
   		c4 = 1.0 - 3.0*ratio + 2.0*SQR(ratio),
   		c5 = 4.0*ratio - 4.0*SQR(ratio),
   		c6 = -1.0*ratio + 2.0*SQR(ratio);
   double	tmp;
   
   while(n-- > 0)
   {
      tmp	= *xo;
      *xo	= c1 * *x + c2 * *xo + c3 * *xvo;
      *xvo	= c4 * *x + c5 * tmp + c6 * *xvo;
      x++; xo++; xvo++;
   }
}
static
void	interpolate_derivatives(system_mp sys, double step, double step1)
{
   double	ratio;
   if( step == 0.0 )
   {
      message(NULLI, NULLP, WARNING, ZEROTS, step1);
      return;
   }
   ratio = step1/step;
   message(NULLI, NULLP, INFO, NEWTS, step, step1);
#ifdef BEEMAN
   interp(ratio, sys->acc[0],   sys->acco[0],  sys->accvo[0],   3*sys->nmols);
   if( sys->nmols_r > 0 )
      interp(ratio, sys->qddot[0], sys->qddoto[0],sys->qddotvo[0], 
	     4*sys->nmols_r);
   interp(ratio, sys->hddot[0], sys->hddoto[0],sys->hddotvo[0], 9);
   interp(ratio, sys->tadot, sys->tadoto, sys->tadotvo, sys->nspecies);
   interp(ratio, sys->radot, sys->radoto, sys->radotvo, sys->nspecies);
#endif
}
/******************************************************************************
 *  check_sysdef.   		Read in system specification from the restart *
 *  file and check its consistency with the previously defined spec.  This is *
 *  for use when replacing the system spec in the middle of a run.  The system*
 *  spec from restart is discarded.  Only the minimum checks are applied for  *
 *  the dynamic variables stored in the restart file to make sense with the   *
 *  new system spec.  Checked are a) Number of species, b) Number of molecules*
 *  of each species and c) that the molecules have the same rotational degrees*
 *  of freedom.  The number of sites can change subject to c).                *
 ******************************************************************************/
static
void	check_sysdef(FILE *restart, char *vsn, system_mp system, spec_mt *species)
    		         		/* Restart file pointer		      */
                                        /* restart file version.              */
         	       			/* NEW 'system' struct		      */
       	          		        /* NEW 'species' struct array	      */
{
   system_mt	sys_tmp;		/* Local temporaries of system,       */
   spec_mp	spec_tmp, spec;		/* species, site_info and potpar      */
   site_mp	site_tmp;		/* used when overwriting restart      */
   pot_mp	pot_tmp;		/* sysdef.			      */

   re_re_sysdef(restart, vsn, &sys_tmp, &spec_tmp, &site_tmp, &pot_tmp);

   if(system->nspecies != sys_tmp.nspecies)
      message(NULLI, NULLP, FATAL, NSPCON, system->nspecies, sys_tmp.nspecies);
   for (spec = species; spec < species+system->nspecies; spec++, spec_tmp++)
   {
      if(spec->nmols != spec_tmp->nmols)
         message(NULLI, NULLP, FATAL, NMLCON, spec_tmp->name,
                 spec->nmols, spec_tmp->nmols);
      if(spec->rdof != spec_tmp->rdof)
         message(NULLI, NULLP, FATAL, NDFCON, spec_tmp->name,
                 spec->rdof, spec_tmp->rdof);
      xfree(spec_tmp->p_f_sites);
      xfree(spec_tmp->site_id);
   }
   xfree((spec_tmp-system->nspecies));
   xfree(site_tmp);
   xfree(pot_tmp);
}
/******************************************************************************
 * Init_cutoffs.  Set the initial values of the alpha parameter and cutoffs   *
 * for the Ewald Sum according to the formulae in Fincham, (1993) CCP5 38, p17*
 * alpha = (nsites*pi^3/V^2 tr/tf)^1/6 is the optmum value		      *
 * r_c = p^1/2 / alpha and k_c = 2 * alpha / p^1/2                            *
 * are the cutoffs which giva an accuracy of exp(-p). Here use 10^-5.         *
 * Any value explicity specified in the control file is left alone.           *
 ******************************************************************************/
#define LOGACC   11.5 /* p =11.5 <=> acc = 10^-5 */
#define TRoverTF 5.5  /* Ratio of indiv. interaction times - approx */
void      init_cutoffs(double *alpha, double *cutoff, double *k_cutoff, real (*h)[3], int nsites)
{
   double max_cutoff = MIN3(h[0][0], h[1][1], h[2][2]);
   double vol = det(h);

   if( *alpha == 0.0 )
      *alpha = pow(nsites*CUBE(PI)/SQR(vol)*TRoverTF, 1.0/6.0);
   if( *alpha > 0.0 )
   {
      if( *cutoff == 0.0 )
      {
	 *cutoff = sqrt(LOGACC)/ *alpha;
	 if( *cutoff > max_cutoff )
	 {
	    message( NULLI, NULLP, WARNING, MAXCUT, *cutoff, max_cutoff);
	    *cutoff = max_cutoff;
	 }
      }
      if( *k_cutoff == 0.0 )
	 *k_cutoff = 2.0* (*alpha) * sqrt(LOGACC);
   }
   else if ( *cutoff <= 0.0 )
      message(NULLI, NULLP, FATAL, NOCUT, *cutoff);
}
/******************************************************************************
 * validate_control.  Determine whether control parameters are sensible.      *
 ******************************************************************************/
void validate_control(void)
{   
   int nerrs = 0;
   if( control.nsteps < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.nsteps, "nsteps");
   if( control.nbins <= 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.nbins, "nbins");
   if( control.begin_average < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.begin_average, "begin-averages");
   if( control.begin_rdf < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.begin_rdf, "begin-rdf");
   if( control.begin_dump < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.begin_dump, "begin-dump");
   if( control.dump_interval < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.dump_interval, "dump-interval");
   if( control.rdf_interval < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.rdf_interval, "rdf-interval");
   if( control.rdf_out < 0 || (control.rdf_interval > 0 && control.rdf_out % control.rdf_interval != 0))
      message(&nerrs, NULLP, ERROR, INVVAL, control.rdf_out, "rdf-out");
   if( control.average_interval < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.average_interval, "average-interval");
   if( control.scale_interval < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.scale_interval, "scale-interval");
   if( control.print_interval <= 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.print_interval, "print-interval");
   if( control.roll_interval <= 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.roll_interval, "roll-interval");
   if( control.backup_interval < 0 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.backup_interval, "backup-interval");
   if( control.cpu_limit < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.cpu_limit, "cpu-limit");
   if( control.density < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.density, "density");
   if( control.pressure < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.pressure, "pressure");
   if( control.temp < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.temp, "temperature");
   if( control.k_cutoff < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.k_cutoff, "k-cutoff");
   if( control.cutoff < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.cutoff, "cutoff");
   if( control.limit <= 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.limit, "rdf-limit");
   if( control.subcell < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.subcell, "subcell");
   if( control.const_pressure < 0 || control.const_pressure > 2 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.const_pressure, "const-pressure");
   if( control.const_temp < 0 || control.const_temp > 2 )
      message(&nerrs, NULLP, ERROR, INVVAL, control.const_temp, "const-temp");
   if( control.rtmass < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.rtmass, "rtmass");
   if( control.ttmass < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.ttmass, "ttmass");
   if( control.pmass < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, control.pmass, "w");
   if( input_unit.m < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, input_unit.m, "mass-unit");
   if( input_unit.l < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, input_unit.l, "length-unit");
   if( input_unit.t < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, input_unit.t, "time-unit");
   if( input_unit.q < 0.0 )
      message(&nerrs, NULLP, ERROR, INVVLF, input_unit.q, "charge-unit");
   if( nerrs > 0 )
      message(&nerrs,NULLP,FATAL,ERRCON,nerrs,(nerrs>1)?'s':' ');
}
/******************************************************************************
 *  startup	This function sets up everything that is needed to start a    *
 *  run.  It controls the reading in of the control, system specification and *
 *  restart files, conversions to program units, calculation of molecular     *
 *  mass & moments of inertia and evaluation of 'whole system' quantities eg  *
 *  number of molecules.  						      *
 ******************************************************************************/
void	start_up(char *contr_name, char *out_name, system_mp system, spec_mp *species, site_mp *site_info, pot_mp *potpar, restrt_mt *restart_header, int *backup_restart)
    		            		/* Name of control file "" for stdin  */
   		          		/* Name of output file "" for stdout  */
         	       			/* Pointer to system struct	      */
       		         		/* Pointer to species array	      */
       		           		/* Pointer to site_info array	      */
      		        		/* Pointer to pot'l parameter array   */
         	                	/* Pointer to restart hdr info (out)  */
   		                	/* (ptr to) flag said purpose   (out) */
{
   FILE		*contr_file,		/* File pointer for control read      */
   		*sysdef,		/* File pointer for sysdef file read  */
   		*backup = NULL,		/* File pointer for backup file read  */
   		*restart = NULL,	/* File pointer for restart file read */
   		*lock;			/* File pointer for lockfile	      */
   double	old_step;		/* Timestep read from restart file    */
   long		old_dump_interval;	/* To check if altered on restart     */
   int		old_max_dumps;		/* To check if altered on restart     */
   long		old_roll_interval;	/* To check if altered on restart     */
   long		old_rdf_interval;	/* To check if altered on restart     */
   long		old_rdf_out;		/* To check if altered on restart     */
   long		old_begin_rdf;		/* To check if altered on restart     */
   int		old_const_pressure;     /* To check if altered on restart     */
   int		old_nbins;              /* To check if altered on restart     */
   double       old_limit;              /* To check if altered on restart     */
   boolean	flag;			/* Used to test 'fseek'		      */
   long		pos;			/* Where control info starts on input */
   restrt_mt	backup_header;		/* To read backup file header into    */
   contr_mt	backup_control;		/* Control struct from backup file    */
   quat_mt	*qpf=0;			/* Quat of rotation to princ. frame   */
   int		av_convert;		/* Flag for old-fmt averages in restrt*/
   int		i;
   *backup_restart = 0;
   (void)memst(restart_header,0,sizeof(*restart_header));

   if(contr_name[0] == '\0')		/* Null name - read control from      */
      contr_file = stdin;		/* standard input.		      */
   else					/* Open named file for reading control*/
   {
      contr_file = fopen(contr_name,"r");
      if(contr_file == NULL)
         message(NULLI, NULLP, FATAL, OCFAIL, contr_name, strerror(errno));
   }
   pos = ftell(contr_file);		/* Current file pos needed for cray   */
   default_control();			/* Set up default values of params    */
   read_control(contr_file, match);	/* Do keyword read of control params  */
   conv_control(&prog_unit,true);	/* Convert to program units           */
   validate_control();

   if(control.restart_file[0] != '\0')	/* Open restart file, get backup name */
   {
      if((restart = fopen(control.restart_file,"rb")) == NULL)
         message(NULLI, NULLP, FATAL, ORFAIL, control.restart_file, 
		 strerror(errno));
      re_re_header(restart, restart_header, &control);
      /*
       *  Now reread control file to override restart file defaults.
       *  Need to do this here in case backup file name changed.
       */
      old_step               = control.step; 		/* Needed for scaling */
      old_roll_interval	     = control.roll_interval;   /* In case it changed */
      old_dump_interval	     = control.dump_interval;	/* Check if these par */
      old_max_dumps	     = control.maxdumps;	/* -amaters altered.  */
      old_rdf_interval	     = control.rdf_interval;
      old_begin_rdf	     = control.begin_rdf;
      old_rdf_out	     = control.rdf_out;
      old_const_pressure     = control.const_pressure;
      old_nbins              = control.nbins;
      old_limit              = control.limit;
      conv_control(&prog_unit, false);
      control.scale_end     -= control.istep;		/* These parameters   */
      control.begin_average -= control.istep;		/* are respecified    */
      control.begin_rdf     -= control.istep;		/* RELATIVE to current*/
      control.begin_dump    -= control.istep;		/* time step.	      */
      control.nsteps	    -= control.istep;

      flag = fseek(contr_file,pos,0);			/* Rewind control file*/
      if( flag )
         message(NULLI,NULLP,FATAL,SEFAIL,contr_name[0]?contr_name:"stdin",
		 strerror(errno));
      read_control(contr_file, match);			/* Reread control file*/
      conv_control(&prog_unit,true);			/* Back to prog units */
   }
   /*
    *  At this point we have read the control file, read the control
    *  parameters from the restart file (if specified).  Now we
    *  know the name of the backup file (Whew!).  Now we have
    *  all the parameter values we can set up and test the lockfiles.
    */
   if( ithread == 0 && control.backup_file[0] )
      (void)strcat(strncpy(backup_lockname,control.backup_file,L_name-5),
		   LOCKEX);
   if( ithread == 0 && control.dump_file[0] )
   {
      (void)sprintf(dump_lockname, control.dump_file, 0);
      (void)strncat(dump_lockname, LOCKEX, L_name-1);
   }
   if( backup_lockname[0] )
   {
      if( fopen(backup_lockname,"r") )
	 message(NULLI, NULLP, -FATAL, LOCKED, "backup", backup_lockname);
      if( (lock=fopen(backup_lockname,"w")) != 0 )
	 (void)fclose(lock);
      else
	 message(NULLI, NULLP, WARNING, LOCKFL, backup_lockname);
   }
   if( dump_lockname[0] )
   {
      if( fopen(dump_lockname,"r") )
      {
	 if( backup_lockname[0] )
	    remove(backup_lockname);
	 message(NULLI, NULLP, -FATAL, LOCKED, "dump", dump_lockname);
      }
      if( (lock=fopen(dump_lockname,"w")) != 0 )
	 (void)fclose(lock);
      else
	 message(NULLI, NULLP, WARNING, LOCKFL, dump_lockname);
   }

   /*
    *  Check for the existance of a backup file and restart from it.
    */
   if( control.backup_file[0] != '\0' &&
      ( backup = fopen(control.backup_file,"rb")) != NULL )  /* Backup exists */
   {
      re_re_header(backup, &backup_header, &backup_control);
      if(restart && (backup_header.timestamp < restart_header->timestamp))
         backup = NULL;			/* Backup older than restart-don't use*/
   }

   if( backup )
   /* We are restarting from backup file. File is already open and header and *
    * control struct have been read into backup_header and backup_control resp*
    * 1) Set number of steps to number yet to do.			      *
    * 2) Copy backup header and control into the correct structs.	      *
    * 3) Read rest of backup file as for ordinary restart, omitting reread of *
    *    control and scaling since nothing can have changed.		      */
   {
      control = backup_control;
      (void)strcpy(restart_header->init_date, backup_header.init_date);
      (void)strcpy(restart_header->title,backup_header.title);
      re_re_sysdef(backup,backup_header.vsn,system,species,site_info,potpar);
      allocate_dynamics(system, *species);/* Memory for dynamic variables     */
      init_averages(system->nspecies, backup_header.vsn,
		    control.roll_interval, control.roll_interval,&av_convert);
      if(control.rdf_interval > 0)
         init_rdf(system);		/* Prepare to calculate rdf	      */
      read_restart(backup, backup_header.vsn, system, av_convert); 
                                        /* Saved dynamic vars etc             */
      convert_averages(control.roll_interval, control.roll_interval,av_convert);
      (void)fclose(backup);
      (*backup_restart)++;
      message(NULLI, NULLP, INFO, BACKUP, control.backup_file);
   }
   else if( !restart )			
   /* Initiate a new run from scratch, ie not restart			      *
    * 1) Read the system specification.  This is either in a file of its own  *
    *    or follows the control info.                                         *
    * 2) Convert potential parameters, site masses and charges to prog. units.*
    * 3) Calculate molecular masses, moments of inertia and nmols etc.        *
    * 4) Allocate memory for the MD dynamic variables.			      *
    * 5) Call the dynamic var initialisation routine.			      *
    * 6) Set up averages data structures.				      */
   {
      if( control.sysdef[0] == '\0' || strcmp(control.sysdef,contr_name) == 0 )
         sysdef = contr_file;		/* Sys def'n is tacked onto control   */
      else				/* Sys def'n is in separate file      */
      {					/* Open system specification file     */
         sysdef = fopen(control.sysdef,"r");
         if(sysdef == NULL)
            message(NULLI, NULLP, FATAL, ODFAIL,control.sysdef,strerror(errno));
      }
					/* Read system specification file     */
      read_sysdef(sysdef, system, species, site_info, potpar);

      conv_potentials(&input_unit, &prog_unit, *potpar, system->n_potpar,
                         system->ptype, *site_info, system->max_id);
      qpf = qalloc(system->nspecies);
      initialise_sysdef(system, *species, *site_info, qpf);
      allocate_dynamics(system, *species);	/* Allocate dynamic arrays    */

      smdrand(control.seed);			/* Seed random number generato*/
      if( control.lattice_start )		/* Choose startup method      */
	 lattice_start(sysdef, system, *species, qpf); /* Lattice - from file */
      else	
         skew_start(system, *species);		/* Start from skew lattice    */
      thermalise(system, *species);		/* Initialise velocities      */
      
      init_averages(system->nspecies, (char*)0,
		    control.roll_interval, control.roll_interval ,&av_convert);
      if( control.rdf_interval > 0 )
         init_rdf(system);			/* Prepare to calculate rdf   */
      if(control.limit <= 0.0)			/* Choose RDF limit           */
	 control.limit = 0.5*MIN3(system->h[0][0],system->h[1][1],
				  system->h[2][2]);

      init_cutoffs(&control.alpha, &control.cutoff, &control.k_cutoff, 
		   system->h, system->nsites);

      (void)strcpy(restart_header->init_date, atime());
      (void)strcpy(restart_header->title,control.title);

      if(sysdef != contr_file)
         (void)fclose(sysdef);			/* Finished with sys spec file*/
   }
   else
   /* Continue from a restart file.  The restart file contains a header, the  *
    * saved 'control' struct, system specification, dynamic variables and     *
    * averages.  The old 'control' values act as defaults for this run. To do *
    * this they are read into 'control', and the control file is rewound and  *
    * reread.  To keep the units consistent, the old control is converted to  *
    * the NEW input units before rereading, and then back to program units.   *
    * Also the scale_end and begin_average parameters are incremented by the  *
    * initial timestep only if they are explicitly specified.                 *
    * By the time we get here, this has all been done.			      *
    *									      *
    * The next thing to be read is the system spec.  Normally this is just    *
    * taken from the restart file, but a facility is provided to allow its    *
    * replacement with a setup from a new system spec 'source' file.  Its     *
    * consistency with the old sys-spec is only checked in so far as the size *
    * and shape of the dynamic variable arrays are identical. Beware.         *
    * Memory is allocated for the dynamic variables and averages which are    *
    * read from the input file.  Finally if the timestep has been altered the *
    * accelerations at previous timesteps are adjusted to the new step.       */
   {
      control.scale_end     += control.istep;
      control.begin_average += control.istep;
      control.begin_rdf     += control.istep;
      control.begin_dump    += control.istep;
      control.nsteps	    += control.istep;
      if( control.istep >= control.begin_dump && 
	  (control.maxdumps != old_max_dumps ||		/* Need to restart    */
	   control.dump_interval != old_dump_interval) )/* dump seq if changed*/
      {
	 message(NULLI, NULLP, WARNING, DMPAL2);
         control.begin_dump = control.istep + 1;	/* Set new beginning  */
      }

      if( !control.new_sysdef )		/* Usual case, get sysdef from restart*/
         re_re_sysdef(restart, restart_header->vsn, 
		      system, species, site_info, potpar);
      else				/* Get sysdef from new sys-spec file  */
      {
         if( control.sysdef[0]=='\0' || strcmp(control.sysdef,contr_name) == 0 )
            sysdef = contr_file;	/* Sys def'n is tacked onto control   */
         else				/* Sys def'n is in separate file      */
         {				/* Open system specification file     */
            sysdef = fopen(control.sysdef,"r");
            if( sysdef == NULL )
               message(NULLI,NULLP,FATAL,OSFAIL,control.sysdef,strerror(errno));
         }
	 /* Read in and set up new system spec (just as in case of new run)   */
         read_sysdef(sysdef, system, species, site_info, potpar);
         conv_potentials(&input_unit, &prog_unit, *potpar, system->n_potpar,
                            system->ptype, *site_info, system->max_id);
#ifdef	DEBUG
         printf(" *D* Read and converted new system specification\n");
#endif
	 qpf = qalloc(system->nspecies);
         initialise_sysdef(system, *species, *site_info, qpf);
	 /* Consistent with saved one?*/
         check_sysdef(restart, restart_header->vsn, system, *species);
         if(sysdef != contr_file)
            (void)fclose(sysdef);
         control.reset_averages = 1;	/* Averages invalid if sysdef changed */
	 control.new_sysdef = 0;        /* Don't store new-sys-spec in restart*/
      }
      allocate_dynamics(system, *species);/* Memory for dynamic variables     */
      init_averages(system->nspecies, restart_header->vsn,
		    control.roll_interval, old_roll_interval, &av_convert);
      if(control.rdf_interval > 0)
         init_rdf(system);		/* Prepare to calculate rdf	      */
      if( (control.rdf_interval != old_rdf_interval ||
	   control.nbins != old_nbins ||
	   control.limit != old_limit) &&
	   control.istep > control.begin_rdf )
      {  
	 message(NULLI, NULLP, WARNING, RDFALT);
	 control.begin_rdf = control.istep + 1;
      }
      if((control.istep-old_begin_rdf) % old_rdf_out > control.rdf_out) 
	 message(NULLI, NULLP, WARNING, RDFDIS);

      read_restart(restart, restart_header->vsn, system, av_convert);  
                                        /* Saved dynamic vars etc             */
      convert_averages(control.roll_interval, old_roll_interval, av_convert);
      control.reset_averages = 0;                /* This flag never propagated.*/

      if(control.step != old_step)
         interpolate_derivatives(system, old_step, control.step);
#ifdef BEEMAN
      for(i = 0; i < 9; i++)		/* Zap cell velocities if constrained */
	 if( (control.strain_mask >> i) & 1 )
	    system->hddot[0][i] = system->hddoto[0][i] = system->hdot[0][i] = 0;
      if(control.const_pressure == 2 && old_const_pressure == 1) 
      {                      /* Enforce consistency if const-p method changed */
	 for(i = 0; i < 9; i++)		           /* Zap cell velocities etc */
	    system->hddot[0][i] = system->hddoto[0][i] = system->hdot[0][i] = 0;
      }
#endif
      (void)fclose(restart);
      message(NULLI, NULLP, INFO, RESUCC, control.restart_file);
   }               
   (void)fclose(contr_file);

   if( out_name[0] != '\0' )	/* Open output file (or use stdout)   */
   {
      (void)fflush(stdout);		/* Purge buffer before opening file   */
      if( freopen(out_name, "a", stdout) == NULL )
         message(NULLI, NULLP, FATAL, OOFAIL, out_name, strerror(errno));
   }
   if( ithread == 0 )
      banner_page(system, *species, restart_header);
   if( qpf != NULL )
      xfree(qpf);
}
