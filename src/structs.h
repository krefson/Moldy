/*
 * $Header: structs.h,v 1.3 89/05/15 15:48:13 keith Exp $
 *
 * $Log:	structs.h,v $
 * Revision 1.2  89/05/11  13:50:32  keith
 * Modified restrt_t to allow commensurate sun3/sun4 padding
 * 
 * Revision 1.1  89/04/27  14:44:58  keith
 * Initial revision
 * 
 * 
 */
#ifndef	STRUCT_ALREADY
#define	STRUCT_ALREADY

#include <stdio.h>
#include "defs.h"

#define		SFORM	"%127[^#]" /* Format for scanf to read strings safely */
typedef struct			/* Control parameters for simulation	      */
{
   char		title[L_name];	/* Job title				      */
   int		istep,		/* Current timestep - used as loop counter    */
		nsteps; 	/* Number of timesteps to execute	      */
   double 	step;		/* Value of timestep in program units	      */
   boolean	print_sysdef,	/* Flag to print out system specification file*/
		new_sysdef,	/* Read new sysdef instead of restart file one*/
		const_pressure, /* Flag to turn on P&R CP method	      */
		reset_averages,	/* Flag to set average counters to zero       */
   		lattice_start;	/* Flag to read starting state from sysdef    */
   char 	sysdef[L_name],		/* Name of system specification file  */
		restart_file[L_name],	/* Name of restart configuration file */
		save_file[L_name],	/* Name of file to write restart conf */
		dump_file[L_name],	/* Name of file 'dump' writes to      */
   		backup_file[L_name],	/* Name of backup save file	      */
   		temp_file[L_name],	/* Temporary file for writing restart */
   	 	out_file[L_name]; 	/* Name of main output file	      */
   FILE		*out;		/* Pointer to main output file		      */
   int		nbins;		/* Number of bins for rdf calculation         */
   unsigned long seed; 		/* Seed for random number generator	      */
   int		pad;		/* Padding for structure alignment	      */
   int		page_width,	/* Line width for output file		      */
		page_length,	/* Length of page on output file	      */
		scale_interval, /* Number of timesteps between scales	      */
		scale_end,	/* Stop scaling after n timesteps	      */
		begin_average,	/* Number of 'equilibration' steps	      */
		average_interval,/* Frequency of averages calculation	      */
   		begin_dump,	/* When to start storing dumps for analysis   */
   		dump_offset,	/* Used in for dump file names - internal only*/
		dump_interval,	/* Frequency of configuration dumps	      */
   		dump_level,	/* What to dump to file			      */
   		maxdumps,	/* How many dump records in a dump file	      */
   		backup_interval,/* Frequency to save state to backup file     */
		roll_interval,	/* Number of timesteps for rolling avgs       */
		print_interval,	/* Number of timesteps between printouts      */
		begin_rdf,	/* When to start RDF calculation              */
		rdf_interval,	/* Frequency to accumulate rdf data           */
		rdf_out;	/* Frequency to calculate and output rdf      */
   double 	temp,		/* Required temperature 		      */
		pressure,	/* Required pressure			      */
		pmass,		/* Parinello and Rahman W parameter	      */
		cutoff, 	/* Cut off radius			      */
   		subcell,	/* Size of side of interaction cells	      */
		density,	/* Initial density on set-up		      */
		alpha,		/* Convergence factor for Ewald sum	      */
   		k_cutoff,	/* Cutoff in k space for ewald sum	      */
		limit,		/* Limiting distance for rdf calculation      */
   		cpu_limit;	/* Maximum CPU allowed before run is stopped  */
} contr_t, *contr_p;

typedef struct			/* Whole system information		      */
{
   int		nsites, 	/* Total number of sites		      */
		nmols,		/* Total number of molecules/atoms	      */
		nmols_r,	/* Total number of polyatomics		      */
		nspecies,	/* Number of different molecule types	      */
		max_id,		/* Last dimension of potpar array	      */
		d_of_f;		/* Degrees of freedom of whole system	      */
   int		ptype,		/* 0 = LJ, 1= buckingham, 2= MCY	      */
		n_potpar;	/* # parameters for this potential	      */
		/* Dynamic variable arrays for whole system		      */
		/* Dimensions for C of M quantities are [nmols][3]	      */
		/* and for quaternions and derivatives, [nmols_r][4]	      */
   vec_p	c_of_m, 	/* Centre of mass positions		      */
		vel,		/* " " " velocities			      */
		velp,		/* Predicted C of M velocities		      */
		acc,		/* C of M accelerations 		      */
		acco,		/* " " at previous timestep		      */
		accvo;		/* " " two timesteps before		      */
   quat_p	quat,		/* Quaternions for this component	      */
		qdot,		/* Quaternion derivatives		      */
		qdotp,		/* Predicted quaternion derivatives	      */
		qddot,		/* Quaternion second derivatives	      */
		qddoto, 	/* Old quaternion second derivatives	      */
		qddotvo;	/* Second derivatives two timesteps before    */
   mat_p	h,		/* Unit cell for zero-stress simulation       */
		hdot,		/* Unit cell derivatives		      */
		hdotp,		/* Predicted unit cell derivatives	      */
		hddot,		/* Unit cell second derivatives 	      */
		hddoto, 	/* Old unit cell second derivatives	      */
		hddotvo;	/* Very old unit cell second derivatives      */
} system_t, *system_p;


typedef struct			/* Information for one species		      */
{
   real		inertia[3],	/* Principal moments of inertia 	      */
		mass,		/* Mass of whole molecule		      */
		dipole;		/* Dipole Moment			      */
   int		nsites, 	/* Number of sites on this species	      */
		nmols;		/* Number of molecules of this species	      */
   int		rdof;		/* Rotational degrees of freedom (2=linear)   */
   char 	name[32];	/* Name of this species 		      */
   int		*site_id;	/* site identifier array		      */
   vec_p	p_f_sites;	/* Site co-ordinates in principal frame       */
		/* Dynamic variable arrays for this species		      */
		/* These point to a subset of the whole-system arrays	      */
		/* Dimensions for C of M quantities are [nmols][3]	      */
		/* and for quaternions and derivatives, [nmols][4]	      */
		/* If species is monatomic, quaternion pointers are null      */
   vec_p	c_of_m, 	/* Centre of mass positions		      */
		vel,		/* " " " velocities			      */
		velp,		/* Predicted C of M velocities		      */
		acc,		/* C of M accelerations 		      */
		acco,		/* " " at previous timestep		      */
		accvo;		/* " " two timesteps before		      */
   quat_p	quat,		/* Quaternions for this species 	      */
		qdot,		/* Quaternion derivatives		      */
		qdotp,		/* Predicted quaternion derivatives	      */
		qddot,		/* Quaternion second derivatives	      */
		qddoto, 	/* Old quaternion second derivatives	      */
		qddotvo;	/* Second derivatives two timesteps before    */
   int		pad;		/* Padding to align on 8-byte boundary	      */
} spec_t, *spec_p;

typedef struct			/* site info template.			      */
{
   double	mass,
   		charge;
   char		name[8];
   int		flag;
   int		pad;
}	site_t,  *site_p;
   		
#define NPOTP 5                 /* Must be number of doubles in pot_t        */
typedef struct			/* Holds potential parameter information      */
{
   int		flag;
   int		pad;
   real		p[NPOTP];
} pot_t, *pot_p;

typedef struct			/* Units used for program input		      */
{
   double	m,		/* mass					      */
		l,		/* length				      */
		t,		/* time					      */
		q;		/* charge				      */

} unit_t, *unit_p;

typedef struct                  /* Record of dimensions of physical quantity  */
{
  int           m,              /* Number of powers of mass in unit           */
                l,
                t,
                q;
} dim_t, *dim_p;

#define DLEN	28		/* Length of date/time string		      */
typedef	struct			/* Restart file header format 		      */
{
   time_t	timestamp,	/* Date and time restart file was written     */
   		prev_timestamp;	/* Timestamp of preceding restart file	      */
   char		init_date[DLEN],/* Date run was initiated (propagated through)*/
		title[L_name],	/* Title when run was initiated		      */
		vsn[16];	/* Version SID of program that wrote restart  */
   int		seq;		/* Sequence NO.  eg 5th restart in run        */
}	restrt_t;

typedef struct			/* Dump file header format		      */
{
   char		title[L_name],	/* Run title at beginning of dump run	      */
   		vsn[16];	/* RCS Revision number			      */
   int		istep,		/* Timestep at beginning of this file	      */
   		dump_interval,	/* How many steps between dumps		      */
   		dump_level,	/* Parameter determining contents of record   */
   		maxdumps,	/* Maximum number of dump records in file     */
   		ndumps,		/* How many dump records in file	      */
   		dump_size;	/* Size of a dump record		      */
   time_t	timestamp,	/* Time file was written		      */
   		dump_init,	/* Time dump run was started (ie first file)  */
   		restart_timestamp;/* Time corresponding restart file written  */
}	dump_t;
#endif
