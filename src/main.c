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
 *      Revision Log
 *       $Log: main.c,v $
 * Revision 2.5  94/01/18  13:32:42  keith
 * Null update for XDR portability release
 * 
 * Revision 2.3  93/10/28  10:27:59  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:49:13  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.20  93/03/09  15:58:58  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.19  92/10/28  14:09:30  keith
 * Changed "site_[tp]" typedefs to avoid name clash on HP.
 * 
 * Revision 1.18  92/06/12  12:55:58  keith
 * Mods to make it work on VMS again.  Ugh.
 * 
 * Revision 1.17  92/06/11  20:31:49  keith
 * Added file locking against multiple runs using same dump or backup files.
 * 
 * Revision 1.16  91/08/15  18:12:06  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.15  91/03/12  15:43:04  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.14  91/02/21  15:27:19  keith
 * Mods for parallel version for titan added
 * 
 * Revision 1.13  90/05/16  14:20:04  keith
 * *** empty log message ***
 * 
 * Revision 1.12  90/04/14  17:53:41  keith
 * Added signal handler to catch CPU exceeded and TERM signal.
 * 
 * Revision 1.11  89/12/15  12:57:00  keith
 * Now prints elapsed as well as cpu time.
 * 
 * Revision 1.10  89/11/20  18:04:15  keith
 * Moved initialisation of control and units to 'input.c'
 * Added 2nd command line arg to specify output file.
 * 
 * Revision 1.9  89/11/20  12:02:09  keith
 * Changed interface to print_rdf.  cf rdf.c 1.6
 * Modified write of restart and backup files - added 'purge' call.
 * 
 * Revision 1.8  89/09/04  18:48:31  keith
 * Chhanged initialisation of 'control' commensurate with structs 1.6.1.2
 * 
 * Revision 1.7  89/08/10  17:30:54  keith
 * Fixed if statement so that rdf's started on rather than after 'begin-rdf'
 * 
 * Revision 1.6  89/07/05  18:19:37  keith
 * Code to support the portable text mode save configuration added.
 * Calculation of when to output averages and rdf data fixed to print
 * when (istep-start+1) % interval == 0. 
 * 
 * Revision 1.5  89/06/01  21:24:38  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.4  89/05/24  11:08:19  keith
 * Fixed bug which called 'averages()' before and at begin-average.
 * Velocities are now rescaled up to and including scale_end.
 * Put better defaults for 'control' parameters.
 * 
 * Revision 1.3  89/05/22  14:05:34  keith
 * Added rescale-separately option, changed 'contr_t' format.
 * 
 * Revision 1.2  89/04/21  10:48:38  keith
 * Corrected bug which left step counter 1 too high at end of run
 * 
 * Revision 1.1  89/04/20  16:00:48  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/main.c,v 2.5.1.1 1994/02/03 18:36:12 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== System include files ============================*/
#include	<signal.h>
#ifdef  SIGCPULIM			/* Alternative name to SIGXCPU.	      */
#define SIGXCPU SIGCPULIM
#endif		/* Unicos uses SIGCPULIM not SIGXCPU. */
#include	<stdio.h>
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External function declarations ==================*/
void	start_up();
void	do_step();
void	values();
void	averages();
void	output();
void	rescale();
void	dump();
void	print_rdf();
void	print_config();
double	cpu();
void	write_restart();
void	purge();
double  rt_clock();
gptr    *talloc();		       /* Interface to memory allocator       */
void    tfree();		       /* Free allocated memory	      	      */
#if defined(ANSI) || defined(__STDC__)
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
#else
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
#endif
void	rmlockfiles();			/* Delete all lock files.	      */
/*========================== External data definition ========================*/
contr_mt control;                           /* Main simulation control parms. */
/*============================================================================*/
/******************************************************************************
 *  Signal handler.  Just set flag and return.				      *
 ******************************************************************************/
static int	sig_flag = 0;
static
void	shutdown(sig)
int sig;
{
   sig_flag = sig;
}
static
void	siglock(sig)
int sig;
{
   rmlockfiles();
   signal(sig, SIG_DFL);
   raise(sig);
}
/******************************************************************************
 *  Main program.							      *
 ******************************************************************************/
int main(argc, argv)
int	argc;
char	*argv[];
{
   system_mt	system;
   spec_mt	*species;
   site_mt	*site_info;
   pot_mt	*potpar;
   restrt_mt	restart_header;
   int		backup_restart;
   mat_mt	stress_vir;
   double	pe[NPE];
   double	delta_cpu = 0.0, cpu_base = cpu();
   double	rt = rt_clock();
   vec_mt	(*meansq_f_t)[2];
   vec_mt	dip_mom;

#ifdef PARALLEL
# ifdef ardent
   int nthreads = nprocessors();
   int stacksize = 65536;
# endif
#endif

   start_up((argc>1)?argv[1]:"", (argc>2)?argv[2]:"",
	    &system, &species, &site_info, &potpar, 
	    &restart_header, &backup_restart);
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
   while( control.istep < control.nsteps &&
	  cpu()-cpu_base+delta_cpu < control.cpu_limit &&
	  sig_flag == 0)
   {
      control.istep++;
      do_step(&system, species, site_info, potpar,
	      meansq_f_t, pe, dip_mom, stress_vir, 
	      &restart_header, backup_restart);
   
      values(&system, species, meansq_f_t, pe, dip_mom, stress_vir);
   
      if(control.istep % control.print_interval == 0)
         output();
      
      if(control.scale_interval > 0)
      {
         if(control.istep <= control.scale_end &&
            control.istep % control.scale_interval == 0)
            rescale(&system, species);
         if(control.istep == control.scale_end)
            note("Temperature scaling turned off after step %d", control.istep);
      }

      if(control.average_interval > 0 && control.istep >= control.begin_average)
      {
         if( control.istep == control.begin_average )
            note("started accumulating thermodynamic averages on timestep %d", 
		 control.istep);
         else if ( (control.istep-control.begin_average + 1) %
		    control.average_interval == 0)
            averages();
      }

      if(control.rdf_interval > 0 && control.istep >= control.begin_rdf && 
        (control.istep-control.begin_rdf+1) % control.rdf_out == 0)
         print_rdf(&system, species, site_info);

      if(control.backup_interval > 0 &&
	 control.istep % control.backup_interval == 0)
      {
	 write_restart(control.backup_file, &restart_header,
		       &system, species, site_info, potpar);
	 purge(control.backup_file);
      }

      if(delta_cpu == 0.0) delta_cpu = cpu() - cpu_base;/* Time for a timestep*/

   }					/* End of main MD timestep loop	      */
   
   if(control.istep < control.nsteps)	/* Run ended prematurely	      */
   {
      if(sig_flag == SIGTERM)
	 note("Run ended after step %d - SIGTERM received", control.istep);
      else
	 note("Run ended after step %d - cpu limit exceeded", control.istep);
      write_restart(control.backup_file, &restart_header, 
		    &system, species, site_info, potpar);
   }
   else if(control.save_file[0] != '\0')
   {
      if( control.print_sysdef )
	 print_config(control.save_file, &system, species, site_info, potpar);
      else
	 write_restart(control.save_file, &restart_header,
		       &system, species, site_info, potpar);
      (void)remove(control.backup_file);		/* Get rid of backup */
   }
   else
      (void)remove(control.backup_file);		/* Get rid of backup */
      
   note("Run used %.2fs of CPU time and %.2fs elapsed", cpu()-cpu_base,
	rt_clock()-rt);

   rmlockfiles();
   return(0);
}
