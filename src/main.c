/******************************************************************************
 * MAIN		Driver program for MOLDY.  Calls startup to read input and    *
 *		control parameters and set up the simulation variables and    *
 *		executes the main cycle over timesteps.  It also contains     *
 *		the definition of the 'control' struct with the default	      *
 *		values of the simulation control parameters.		      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	main.c,v $
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
static char *RCSid = "$Header: /home/tigger/keith/md/RCS/main.c,v 1.9 89/11/20 12:02:09 keith Exp $";
#endif
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
void	note();
void	dump();
void	print_rdf();
void	print_config();
void	message();
double	cpu();
void	write_restart();
void	purge();
/*========================== External data definition ========================*/
contr_t		control;
unit_t		input_unit;
/*============================================================================*/
main(argc, argv)
int	argc;
char	*argv[];
{
   system_t	system;
   spec_t	*species;
   site_t	*site_info;
   pot_t	*potpar;
   mat_t	stress_vir;
   double	pe[NPE];
   double	delta_cpu = 0.0, cpu_base = cpu();
   vec_t	(*meansq_f_t)[2];
   vec_t	dip_mom;

   start_up((argc>1)?argv[1]:"", (argc>2)?argv[2]:"",
	    &system, &species, &site_info, &potpar);
   meansq_f_t = (vec_t (*)[2])ralloc(2*system.nspecies);
   
   while( control.istep < control.nsteps &&
	  cpu()-cpu_base+delta_cpu < control.cpu_limit)
   {
      control.istep++;
      do_step(&system, species, site_info, potpar,
	      meansq_f_t, pe, dip_mom, stress_vir);
   
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
	 write_restart(control.backup_file,&system, species, site_info, potpar);
	 purge(control.backup_file);
      }

      if(delta_cpu == 0.0) delta_cpu = cpu() - cpu_base;/* Time for a timestep*/

   }					/* End of main MD timestep loop	      */
   
   if(control.istep < control.nsteps)	/* Run ended when CPU limit exceeded  */
   {
      note("Run ended after step %d - cpu limit exceeded", control.istep);
      write_restart(control.backup_file, &system, species, site_info, potpar);
   }
   else if(control.save_file[0] != '\0')
   {
      if( control.print_sysdef )
	 print_config(control.save_file, &system, species, site_info, potpar);
      else
	 write_restart(control.save_file, &system, species, site_info, potpar);
      (void)remove(control.backup_file);		/* Get rid of backup */
   }
   else
      (void)remove(control.backup_file);		/* Get rid of backup */
      
   note("Run used %.2fs of CPU time", cpu()-cpu_base);

   return(0);
}
