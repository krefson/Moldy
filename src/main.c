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
 * MAIN		Driver program for MOLDY.  Calls startup to read input and    *
 *		control parameters and set up the simulation variables and    *
 *		executes the main cycle over timesteps.  It also contains     *
 *		the definition of the 'control' struct with the default	      *
 *		values of the simulation control parameters.		      *
 ******************************************************************************
 */
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== System include files ============================*/
#include	<signal.h>
#ifdef  SIGCPULIM			/* Alternative name to SIGXCPU.	      */
#define SIGXCPU SIGCPULIM
#endif		/* Unicos uses SIGCPULIM not SIGXCPU. */
#include	<stdio.h>
#include	<string.h>
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External function declarations ==================*/
void	start_up(char *contr_name, char *out_name, system_mp system, 
		 spec_mp *species, site_mp *site_info, pot_mp *potpar, 
		 restrt_mt *restart_header, int *backup_restart, boolean *);
void	do_step(system_mp sys, spec_mt *species, site_mt *site_info, 
		pot_mt *potpar, vec_mt (*meansq_f_t)[2], double *pe, 
		real *dip_mom, real (*stress)[3], restrt_mt *restart_header, 
		int backup_restart, boolean re_init_H0);
void	values(system_mp system, spec_mt *species, vec_mt (*meansq_f_t)[2], 
	       double *pe, real *dipole, real (*stress_vir)[3]);
double	value(av_n type, int comp);
void	averages(void);
void	output(void);
void	rescale(system_mp system, spec_mp species);
void	print_rdf(system_mt *system, spec_mt *species, site_mt *site_info);
void	print_config(char *save_name, system_mp system, spec_mp species, 
		     site_mp site_info, pot_mp potpar);
double	cpu(void);
void	write_restart(char *save_name, restrt_mt *header, system_mp system, 
		      spec_mp species, site_mp site_info, pot_mp potpar);
void	purge(char *file);
double  rt_clock(void);
gptr    *talloc(int n, size_mt size, int line, char *file);
				       /* Interface to memory allocator       */
void    tfree(gptr *p);		       /* Free allocated memory	      	      */
void	save_version_inc(char *save_file, int namesize);
#ifdef SPMD
void    par_begin(int *argc, char ***argv, int *ithread, int *nthreads);
void    par_sigintreset(void);
void    par_finish(void);
void    par_fsum(float *buf, int n);
void    par_imax(int *idat);
void    par_broadcast(gptr *buf, int n, size_mt size, int ifrom);
void    replicate(contr_mt *control, system_mt *system, spec_mt **spec_ptr, 
		  site_mt **site_info, pot_mt **pot_ptr, 
		  restrt_mt *restart_header);
#endif
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
void	rmlockfiles(void);		/* Delete all lock files.	      */
gptr    *rdf_ptr(int *size);            /* Return ptr to start of rdf data    */
/*========================== External data definition ========================*/
contr_mt control;                           /* Main simulation control parms. */
int ithread=0, nthreads=1;
/*============================================================================*/
/******************************************************************************
 *  Signal handler.  Just set flag and return.				      *
 ******************************************************************************/
static int	sig_flag = 0;
static
void	shutdown(int sig)
{
   sig_flag = sig;
}
static
void	siglock(int sig)
{
   rmlockfiles();
#ifdef SPMD
   if( sig == SIGINT )
      par_sigintreset();
   else   
#endif
   signal(sig, SIG_DFL);
   raise(sig);
}
/******************************************************************************
 *  Main program.							      *
 ******************************************************************************/
int main(int argc, char **argv)
{
   system_mt	system;
   spec_mt	*species;
   site_mt	*site_info;
   pot_mt	*potpar;
   restrt_mt	restart_header;
   int		backup_restart;
   boolean	init_H_0 = false;
   static mat_mt	stress_vir;
   static double	pe[NPE];
#ifdef SPMD
   double t0, t00;
#endif
   double	delta_cpu = 0.0, cpu_base = cpu();
   double	rt = rt_clock();
   vec_mt	(*meansq_f_t)[2];
   vec_mt	dip_mom;
   float        *rdf_base;
   int          rdf_size;

#ifdef SPMD
   par_begin(&argc, &argv, &ithread, &nthreads);
#endif
#ifdef PARALLEL
# ifdef ardent
   int nthreads = nprocessors();
   int stacksize = 65536;
# endif
#endif

#if defined(SPMD) && !defined(READALL)
   if( ithread == 0 )
#endif
      start_up((argc>1)?argv[1]:"", (argc>2)?argv[2]:"",
	       &system, &species, &site_info, &potpar, 
	       &restart_header, &backup_restart, &init_H_0);
#if defined(SPMD) && !defined(READALL)
   replicate(&control, &system, &species, &site_info, &potpar, 
	     &restart_header);
   par_broadcast( &init_H_0, 1, sizeof(boolean), 0);
#endif
   rdf_base = (float*)rdf_ptr(&rdf_size);

   meansq_f_t = (vec_mt (*)[2])ralloc(2*system.nspecies);
   
   /*
    *  Set signal handlers -- attempt clean shutdown
    */
   (void)signal(SIGTERM, shutdown);
#ifdef SIGXCPU
   (void)signal(SIGXCPU, shutdown);
#endif
   /*
    *  Set signal handlers -- remove lock files and exit
    *  Check that the signal is actually defined if it is non-ansi
    */
#ifdef SIGHUP
   (void)signal(SIGHUP, siglock);
#endif
   (void)signal(SIGINT, siglock);
#ifdef SIGQUIT
   (void)signal(SIGQUIT, siglock);
#endif
#ifdef SIGABRT
   (void)signal(SIGABRT, siglock);
#endif
#ifdef PARALLEL
# ifdef ardent
   MT_SET_THREAD_NUMBER(&nthreads);
   MT_INIT(&stacksize);
# endif
#endif
   /*
    *  Main MD timestep loop
    */
   while( control.istep < control.nsteps && sig_flag == 0)
   {
      control.istep++;

      if((control.istep-control.begin_rdf) % control.rdf_out == 0) 
	 memst((gptr*)rdf_base, 0.0, rdf_size*sizeof(float));

      do_step(&system, species, site_info, potpar,
	      meansq_f_t, pe, dip_mom, stress_vir, 
	      &restart_header, backup_restart, init_H_0);

      init_H_0 = false;
      
      values(&system, species, meansq_f_t, pe, dip_mom, stress_vir);
   
      if( ithread == 0 && control.istep % control.print_interval == 0)
	 output();

      if(control.scale_interval > 0)
      {
	 if(control.istep <= control.scale_end &&
	    control.istep % control.scale_interval == 0)
	 {
	    rescale(&system, species);
	    init_H_0 = true;
	 }
	 if( ithread == 0 && control.istep == control.scale_end)
	    note("Temperature scaling turned off after step %ld", control.istep);
      }

      if(control.average_interval > 0 && 
	 control.istep >= control.begin_average &&
	 ithread == 0)
      {
         if( control.istep == control.begin_average )
	 {
	    note("started accumulating thermodynamic averages on timestep %ld", 
		 control.istep);
	 }
         else if ( (control.istep-control.begin_average + 1) %
		    control.average_interval == 0)
            averages();
      }

      if(control.rdf_interval > 0 && control.istep >= control.begin_rdf && 
	 (control.istep-control.begin_rdf+1) % control.rdf_out == 0)
      {
#if defined(SPMD) && ! defined OLDRDF
	 par_fsum(rdf_base, rdf_size);
#endif
	 if( ithread == 0 )
	    print_rdf(&system, species, site_info);
      }

      if(control.backup_interval > 0 && control.backup_file[0] &&
	 control.istep % control.backup_interval == 0)
      {
#if defined(SPMD) && ! defined OLDRDF
	 par_fsum(rdf_base, rdf_size);
	 if( ithread != 0 )
	    memst((gptr*)rdf_base, 0, rdf_size*sizeof(int));
#endif
	 if( ithread == 0 )
	 {
	    write_restart(control.backup_file, &restart_header,
			  &system, species, site_info, potpar);
	    purge(control.backup_file);
	 }
      }
      if(delta_cpu == 0.0) delta_cpu = cpu() - cpu_base;/* Time for a timestep*/
      if(  cpu()-cpu_base+delta_cpu >= control.cpu_limit )
	 sig_flag++; 	/* Cheat */
#ifdef SPMD
      /*
       * Better make sure every process wants to stop at the same time!
       */
      par_imax(&sig_flag);
      /*
       * Vital consistency check that all threads are absolutely in sync.
       * Use Temperature as it's a function of all KEs & therefore trajs.
       */
      if( control.istep%10 == 0)
      {
	 t0 = t00 = value(t_n,0);
	 par_broadcast(&t0, 1, sizeof(double), 0);
	 if( t0 != t00)
	    message(NULLI, NULLP, FATAL, DESYNC, ithread, t0, t00);
      }
#endif
   }					/* End of main MD timestep loop	      */
   
#if defined(SPMD) && ! defined OLDRDF
   par_fsum(rdf_base, rdf_size);
#endif
   if( ithread == 0 )
   {
      if(control.istep < control.nsteps) /* Run ended prematurely	      */
      {
	 if(sig_flag == SIGTERM)
	    note("Run ended after step %ld - SIGTERM received", control.istep);
	 else
	    note("Run ended after step %ld - cpu limit exceeded", control.istep);
	 write_restart(control.backup_file, &restart_header, 
		       &system, species, site_info, potpar);
      }
      else if(control.save_file[0] != '\0')
      {
	 /*
	  * Save-file name management.  Save overwriting.
	  */
	 if( strcmp(control.restart_file,control.save_file) == 0)
	 {
	    save_version_inc(control.save_file, L_name);
	    message(NULLI, NULLP, WARNING, SAVINC, control.save_file);
	 }

	 if( control.print_sysdef )
	    print_config(control.save_file, &system, species, site_info, potpar);
	 else
	    write_restart(control.save_file, &restart_header,
			  &system, species, site_info, potpar);
	 (void)remove(control.backup_file);		/* Get rid of backup */
      }
      else
	 (void)remove(control.backup_file);		/* Get rid of backup */

      rmlockfiles();
#ifdef SPMD
      printf(" *I* Run used %.2fs of CPU time and %.2fs elapsed on %d processors\n", 
	     cpu()-cpu_base, rt_clock()-rt, nthreads);
#else
      printf(" *I* Run used %.2fs of CPU time and %.2fs elapsed\n", 
	     cpu()-cpu_base, rt_clock()-rt);
#endif
   }
#ifdef SPMD
   par_finish();
#endif
   return(0);
}
