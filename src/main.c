/******************************************************************************
 * MAIN		Driver program for MOLDY.  Calls startup to read input and    *
 *		control parameters and set up the simulation variables and    *
 *		executes the main cycle over timesteps.  It also contains     *
 *		the definition of the 'control' struct with the default	      *
 *		values of the simulation control parameters.		      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	main.c,v $
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
static char *RCSid = "$Header: main.c,v 1.7 89/07/04 18:46:55 keith Exp $";
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
/*========================== External data definition ========================*/
/*
 * Control struct with default values for the parameters 
 */
contr_t	control = {
	"Test Simulation",	/* Job title				      */
	0,			/* Current timestep - used as loop counter    */
	0,			/* Number of timesteps to execute	      */
	0.005,			/* Value of timestep in program units	      */
	false,			/* Flag to print out system specification file*/
	false,			/* Read new sysdef instead of restart file one*/
	false,			/* Flag to turn on P&R CP method	      */
	false,			/* Flag to set average counters to zero       */
	false,			/* Flag to read initial state from sysdef file*/
	"",			/* Name of system specification file	      */
	"",			/* Name of file to read restart conf. from    */
	"",			/* Name of file to write restart conf. to     */
	"",			/* Name of file 'dump' writes to              */
	"MDBACKUP",		/* Name of file for periodic save of state    */
#ifdef CMS
	"MDTEMP XXXXXXXX A1",	/* Temporary file name for restart to write   */
#else
#ifdef vms
	"MDTEMPXXXX.DAT",
#else
	"MDTEMPX",
#endif
#endif
	"",			/* Name of main output file		      */
	0,			/* To preserve alignment of structure	      */
	100,			/* Number of bins for rdf calculation	      */
	1234567,		/* Seed for random number generator	      */
	132,			/* Line width for output file		      */
	44,			/* Length of page on output file	      */
	10,			/* Number of timesteps between scales	      */
	1000000,		/* Stop scaling after n timesteps	      */
	1001,			/* Number of 'equilibration' steps	      */
	0,			/* Whether to scale each species separately   */
	100,			/* Frequency of averages calculation	      */
	1,			/* Start of configuration dumps		      */
	0,			/* Dump filename offset (internal use only)   */
	20,			/* Frequency of configuration dumps	      */
	0,			/* Level of dump to perform		      */
	10,			/* How many dump records in a dump file	      */
	100,			/* Frequency to write save configuration      */
	10,			/* Number of timesteps for rolling avgs       */
	10,			/* Number of timesteps between printouts      */
        1001,			/* When to start accumulating rdf data        */
	20,			/* How frequently to perform binning          */
	1000,			/* How frequently to calculate & print rdf    */
	0,			/* Required temperature 		      */
	0.0,			/* Required pressure			      */
	50,			/* Parinello and Rahman W parameter	      */
	10.0,			/* Cut off radius			      */
	0.0,			/* Size of side of interaction cells	      */
        1.0,			/* Density 1g/cc			      */
	0.3,			/* Convergence parameter for Ewald sum	      */
        2.0,			/* K space cutoff for Ewald sum		      */
	10.0,			/* Limiting distance for RDF calculation      */
	1.0e20};		/* Default CPU limit - very large	      */
	

unit_t		input_unit = {MUNIT, LUNIT, TUNIT/10, _ELCHG};
				/* amu, A, ps/10 => energy unit = kJ/mol      */
/*============================================================================*/
main(argc, argv)
int	argc;
char	*argv[];
{
   system_t	system;
   spec_p	species;
   site_p	site_info;
   pot_p	potpar;
   mat_t	stress_vir;
   double	pe[NPE];
   double	delta_cpu = 0.0, cpu_base = cpu();
   vec_t	(*meansq_f_t)[2];
   vec_t	dip_mom;

   start_up((argc>1)?argv[1]:"", &system, &species, &site_info, &potpar);
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

      if(control.rdf_interval > 0 && control.istep > control.begin_rdf && 
        (control.istep-control.begin_rdf+1) % control.rdf_out == 0)
         print_rdf(&system, site_info);

      if(control.backup_interval > 0 &&
	 control.istep % control.backup_interval == 0)
	 write_restart(control.backup_file,&system, species, site_info, potpar);

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
