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
/*
 * $Header: /home/eeyore/keith/md/moldy/RCS/structs.h,v 2.1 93/07/19 13:28:21 keith Exp $
 *
 * $Log:	structs.h,v $
 * Revision 2.1  93/07/19  13:28:21  keith
 * Added XDR capability for backup and dump files.
 * 
 * Revision 2.0  93/03/15  14:49:26  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.6.1.15  93/03/12  12:14:34  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.6.1.15  93/03/09  15:59:22  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.6.1.14  92/10/28  14:09:40  keith
 * Changed "site_[tp]" typedefs to avoid name clash on HP.
 * 
 * Revision 1.6.1.13  92/09/22  14:55:12  keith
 * Added support for "strict cutoff" mode.
 * 
 * Revision 1.6.1.12  92/03/11  12:56:24  keith
 * Changed "scale-separately" parameter to "scale options"
 * 
 * Revision 1.6.1.11  91/08/23  11:34:49  keith
 * Added extra padding to struct spec_t to round size up to 8 byte
 * boundary.  This promotes portability of restart files.  Rounded
 * up rather than down for compatibility of existing restart files.
 * 
 * Revision 1.6.1.10  91/08/15  18:12:21  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.6.1.9  90/05/02  15:44:36  keith
 * Got rid of typedefs time_t and size_t. 
 * 
 * Revision 1.6.1.8  90/04/16  18:20:40  keith
 * Added new field "strain-mask" to control.
 * 
 * Revision 1.6.1.7  90/04/12  16:29:10  keith
 * removed unneccessary include of <stdio.h>
 * 
 * Revision 1.6.1.6  90/04/06  11:09:49  keith
 * Moved definition of NPOTP to defs.h
 * 
 * Revision 1.6.1.5  89/11/21  16:31:32  keith
 * Removed member out_file from control and all uses. (Now command parameter).
 * 
 * Revision 1.6.1.4  89/11/20  18:10:41  keith
 * Added "defalt" field to match_t.
 * 
 * Revision 1.6.1.3  89/11/20  13:30:14  keith
 * Replaced separate arrays "types" and "npotp" with array of structs "potspec"
 * 
 * Revision 1.6.1.2  89/09/04  18:40:19  keith
 * Added 'surface_dipole' to control_t (& removed pad), moved 'scale_separately'
 * Added field 'charge' to spec_t.
 * 
 * Revision 1.6.1.1  89/08/25  15:24:21  keith
 * Mods to add framework structures to simulation model
 * 
 * Revision 1.6  89/06/20  18:25:36  keith
 * Moved definition of match_t to structs.h 
 * 
 * Revision 1.5  89/06/01  21:25:33  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.4  89/05/22  14:05:48  keith
 * Added rescale-separately option, changed 'contr_t' format.
 * 
 * Revision 1.3  89/05/15  16:12:04  keith
 * Added new members, 'vsn' and 'dump_size' to dump_t.
 * * Must use with r1.3 or later of 'dump.c'.
 * 
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
   		scale_options,	/* Scale each species separately	      */
		surface_dipole,	/* Flag surface dipole term in Ewald sum      */
   		lattice_start;	/* Flag to read starting state from sysdef    */
   char 	sysdef[L_name],		/* Name of system specification file  */
		restart_file[L_name],	/* Name of restart configuration file */
		save_file[L_name],	/* Name of file to write restart conf */
		dump_file[L_name],	/* Name of file 'dump' writes to      */
   		backup_file[L_name],	/* Name of backup save file	      */
   		temp_file[L_name];	/* Temporary file for writing restart */
   int	 	spare[29];	/* Extra space for expansion 	              */
   boolean	xdr_write,	/* Write restart, dump files in portable way. */
                strict_cutoff;  /* Perform real-space cutoff rigorously       */
   int		strain_mask;	/* Mask of constrained elements of h matrix   */
   int		nbins;		/* Number of bins for rdf calculation         */
   unsigned long seed; 		/* Seed for random number generator	      */
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
} contr_mt, *contr_mp;

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
   vec_mp	c_of_m, 	/* Centre of mass positions		      */
		vel,		/* " " " velocities			      */
		velp,		/* Predicted C of M velocities		      */
		acc,		/* C of M accelerations 		      */
		acco,		/* " " at previous timestep		      */
		accvo;		/* " " two timesteps before		      */
   quat_mp	quat,		/* Quaternions for this component	      */
		qdot,		/* Quaternion derivatives		      */
		qdotp,		/* Predicted quaternion derivatives	      */
		qddot,		/* Quaternion second derivatives	      */
		qddoto, 	/* Old quaternion second derivatives	      */
		qddotvo;	/* Second derivatives two timesteps before    */
   mat_mp	h,		/* Unit cell for zero-stress simulation       */
		hdot,		/* Unit cell derivatives		      */
		hdotp,		/* Predicted unit cell derivatives	      */
		hddot,		/* Unit cell second derivatives 	      */
		hddoto, 	/* Old unit cell second derivatives	      */
		hddotvo;	/* Very old unit cell second derivatives      */
} system_mt, *system_mp;


typedef struct			/* Information for one species		      */
{
   real		inertia[3],	/* Principal moments of inertia 	      */
		mass,		/* Mass of whole molecule		      */
		dipole,		/* Dipole Moment			      */
   		charge;		/* Total charge				      */
   int		nsites, 	/* Number of sites on this species	      */
		nmols;		/* Number of molecules of this species	      */
   int		rdof,		/* Rotational degrees of freedom (2=linear)   */
   		framework;	/* Flag to signal this is a framework species */
   char 	name[32];	/* Name of this species 		      */
   int		*site_id;	/* site identifier array		      */
   vec_mp	p_f_sites;	/* Site co-ordinates in principal frame       */
		/* Dynamic variable arrays for this species		      */
		/* These point to a subset of the whole-system arrays	      */
		/* Dimensions for C of M quantities are [nmols][3]	      */
		/* and for quaternions and derivatives, [nmols][4]	      */
		/* If species is monatomic, quaternion pointers are null      */
   vec_mp	c_of_m, 	/* Centre of mass positions		      */
		vel,		/* " " " velocities			      */
		velp,		/* Predicted C of M velocities		      */
		acc,		/* C of M accelerations 		      */
		acco,		/* " " at previous timestep		      */
		accvo;		/* " " two timesteps before		      */
   quat_mp	quat,		/* Quaternions for this species 	      */
		qdot,		/* Quaternion derivatives		      */
		qdotp,		/* Predicted quaternion derivatives	      */
		qddot,		/* Quaternion second derivatives	      */
		qddoto, 	/* Old quaternion second derivatives	      */
		qddotvo;	/* Second derivatives two timesteps before    */
   int		pad[2];		/* Needed for compatibility of binary restart */
   				/* files due to historical cock-up.	      */ 
} spec_mt, *spec_mp;

typedef struct			/* site info template.			      */
{
   double	mass,
   		charge;
   char		name[8];
   int		flag;
   int		pad;
}	site_mt,  *site_mp;
   		
typedef struct			/* Holds potential parameter information      */
{
   int		flag;
   int		pad;
   real		p[NPOTP];
} pot_mt, *pot_mp;

typedef struct
{
   char	*name;
   int  npar;
} pots_mt;

typedef struct			/* Units used for program input		      */
{
   double	m,		/* mass					      */
		l,		/* length				      */
		t,		/* time					      */
		q;		/* charge				      */

} unit_mt, *unit_mp;

typedef struct                  /* Record of dimensions of physical quantity  */
{
  int           m,              /* Number of powers of mass in unit           */
                l,
                t,
                q;
} dim_mt, *dim_mp;

typedef	struct				/* Struct template for keyword	      */
{					/* in read_control.		      */
   char	*key,
	*format,
        *defalt;
   gptr	*ptr;
}	match_mt;

#define DLEN	28		/* Length of date/time string		      */
typedef	struct			/* Restart file header format 		      */
{
   time_mt	timestamp,	/* Date and time restart file was written     */
   		prev_timestamp;	/* Timestamp of preceding restart file	      */
   char		init_date[DLEN],/* Date run was initiated (propagated through)*/
		title[L_name],	/* Title when run was initiated		      */
		vsn[16];	/* Version SID of program that wrote restart  */
   int		seq;		/* Sequence NO.  eg 5th restart in run        */
}	restrt_mt;

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
   time_mt	timestamp,	/* Time file was written		      */
   		dump_init,	/* Time dump run was started (ie first file)  */
   		restart_timestamp;/* Time corresponding restart file written  */
}	dump_mt;

#define	MAX_ROLL_INTERVAL	100
typedef struct
{   double	value,
   		sum,
   		sum_sq,
		mean,
		sd,
   		roll[MAX_ROLL_INTERVAL],
   		roll_mean,
   		roll_sd;
} old_av_mt;

typedef	union
{
   old_av_mt		av;
   struct
   {	
      int	av, roll;
   }		cnt;
} old_av_u_mt;

typedef	struct
{	
   int		nav, 
   		nroll, 
		iroll, 
		pad;
   double align;
} av_head_mt;

#endif
