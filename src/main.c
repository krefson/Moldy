/******************************************************************************
 * MAIN		Driver program for MOLDY.  Calls startup to read input and    *
 *		control parameters and set up the simulation variables and    *
 *		executes the main cycle over timesteps.  It also contains     *
 *		the definition of the 'control' struct with the default	      *
 *		values of the simulation control parameters.		      *
 ******************************************************************************
 *      Revision Log
 *       $Log$
 */
#ifndef lint
static char *RCSid = "$Header$";
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
double	value();
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
	NULL,			/* Pointer to main output file		      */
	100,			/* Number of bins for rdf calculation	      */
	1234567,		/* Seed for random number generator	      */
	0,			/* Padding to align structure fields	      */
	132,			/* Line width for output file		      */
	44,			/* Length of page on output file	      */
	0,			/* Number of timesteps between scales	      */
	1000000,		/* Stop scaling after n timesteps	      */
	1000,			/* Number of 'equilibration' steps	      */
	100,			/* Frequency of averages calculation	      */
	0,			/* Start of configuration dumps		      */
	0,			/* Dump filename offset (internal use only)   */
	0,			/* Frequency of configuration dumps	      */
	0,			/* Level of dump to perform		      */
	10,			/* How many dump records in a dump file	      */
	100,			/* Frequency to write save configuration      */
	10,			/* Number of timesteps for rolling avgs       */
	10,			/* Number of timesteps between printouts      */
        0,			/* When to start accumulating rdf data        */
	0,			/* How frequently to perform binning          */
	100,			/* How frequently to calculate & print rdf    */
	0,			/* Required temperature 		      */
	0.0,			/* Required pressure			      */
	50,			/* Parinello and Rahman W parameter	      */
	10.0,			/* Cut off radius			      */
	0.0,			/* Size of side of interaction cells	      */
        1.0,			/* Density 1g/cc			      */
	0.5,			/* Convergence parameter for Ewald sum	      */
        3.5,			/* K space cutoff for Ewald sum		      */
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
   
   while( control.istep++ < control.nsteps &&
	  cpu()-cpu_base+delta_cpu < control.cpu_limit)
   {
      do_step(&system, species, site_info, potpar,
	      meansq_f_t, pe, dip_mom, stress_vir);
   
      values(&system, species, meansq_f_t, pe, dip_mom, stress_vir);
   
      if(control.istep % control.print_interval == 0)
         output();
      
      if(control.scale_interval > 0)
         if(control.istep < control.scale_end &&
            control.istep % control.scale_interval == 0)
            rescale(value(t_n,0), control.temp, &system);
         else if(control.istep == control.scale_end)
            note("temperature scaling turned off");

      if(control.average_interval > 0 && 
        (control.istep-control.begin_average) % control.average_interval == 0)
           averages();

      if(control.istep == control.begin_average)
         note("started accumulating thermodynamic averages on timestep %d", 
              control.istep);

      if(control.rdf_interval > 0 && control.istep > control.begin_rdf && 
        (control.istep-control.begin_rdf) % control.rdf_out == 0)
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
      write_restart(control.save_file, &system, species, site_info, potpar);
      (void)remove(control.backup_file);		/* Get rid of backup */
   }
   else
      (void)remove(control.backup_file);		/* Get rid of backup */
      
   note("Run used %fs of CPU time", cpu()-cpu_base);

   return(0);
}
