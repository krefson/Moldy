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

#ifndef STRUCT_ALREADY
#define STRUCT_ALREADY

#include "defs.h"

typedef struct                  /* Control parameters for simulation          */
{
   char         title[L_name];  /* Job title                                  */
   long         istep,          /* Current timestep - used as loop counter    */
                nsteps;         /* Number of timesteps to execute             */
   double       step;           /* Value of timestep in program units         */
   boolean      print_sysdef,   /* Flag to print out system specification file*/
                new_sysdef,     /* Read new sysdef instead of restart file one*/
                molpbc,		/* Switch to molecular from site cutoff scheme*/
                reset_averages; /* Flag to set average counters to zero       */
   int          scale_options;  /* Scale each species separately              */
   boolean      surface_dipole, /* Flag surface dipole term in Ewald sum      */
                lattice_start;  /* Flag to read starting state from sysdef    */
   char         sysdef[L_name],         /* Name of system specification file  */
                restart_file[L_name],   /* Name of restart configuration file */
                save_file[L_name],      /* Name of file to write restart conf */
                dump_file[L_name],      /* Name of file 'dump' writes to      */
                backup_file[L_name],    /* Name of backup save file           */
                temp_file[L_name];      /* Temporary file for writing restart */
   int          spare[20];      /* Extra space for expansion (should be EVEN) */
   boolean      nosymmetric_rot;/* Don't use symm. variant of rot'l leapfrog  */
   double	ewald_accuracy; /* Accuracy parameter for Ewald sum.          */
   double       ttmass,         /* Nose-Hoover trans temp mass parameter      */
                rtmass;         /* Nose-Hoover rotat temp mass parameter      */
   int		const_pressure; /* Flag to turn on P&R CP method              */
   int          const_temp;     /* Flag to turn on NP themostat method        */
   boolean      xdr_write,      /* Write restart, dump files in portable way. */
                strict_cutoff;  /* Perform real-space cutoff rigorously       */
   int          strain_mask;    /* Mask of constrained elements of h matrix   */
   int          nbins;          /* Number of bins for rdf calculation         */
   unsigned long seed;          /* Seed for random number generator           */
   int          page_width,     /* Line width for output file                 */
                page_length;    /* Length of page on output file              */
   long         scale_interval, /* Number of timesteps between scales         */
                scale_end,      /* Stop scaling after n timesteps             */
                begin_average,  /* Number of 'equilibration' steps            */
                average_interval,/* Frequency of averages calculation         */
                begin_dump,     /* When to start storing dumps for analysis   */
                dump_offset,    /* Used in for dump file names - internal only*/
                dump_interval;  /* Frequency of configuration dumps           */
   int          dump_level,     /* What to dump to file                       */
                maxdumps;       /* How many dump records in a dump file       */
   long         backup_interval,/* Frequency to save state to backup file     */
                roll_interval,  /* Number of timesteps for rolling avgs       */
                print_interval, /* Number of timesteps between printouts      */
                begin_rdf,      /* When to start RDF calculation              */
                rdf_interval,   /* Frequency to accumulate rdf data           */
                rdf_out;        /* Frequency to calculate and output rdf      */
   double       temp,           /* Required temperature                       */
                pressure,       /* Required pressure                          */
                pmass,          /* Parinello and Rahman W parameter           */
                cutoff,         /* Cut off radius                             */
                subcell,        /* Size of side of interaction cells          */
                density,        /* Initial density on set-up                  */
                alpha,          /* Convergence factor for Ewald sum           */
                k_cutoff,       /* Cutoff in k space for ewald sum            */
                limit,          /* Limiting distance for rdf calculation      */
                cpu_limit;      /* Maximum CPU allowed before run is stopped  */
} contr_mt, *contr_mp;

typedef struct                  /* Whole system information                   */
{
   int          nsites,         /* Total number of sites                      */
                nmols,          /* Total number of molecules/atoms            */
                nmols_r,        /* Total number of polyatomics                */
                nspecies,       /* Number of different molecule types         */
                max_id,         /* Last dimension of potpar array             */
                d_of_f;         /* Degrees of freedom of whole system         */
   int          ptype,          /* 0 = LJ, 1= buckingham, 2= MCY              */
                n_potpar;       /* # parameters for this potential            */
                /* Dynamic variable arrays for whole system                   */
                /* Dimensions for C of M quantities are [nmols][3]            */
                /* and for quaternions and derivatives, [nmols_r][4]          */
   vec_mp       c_of_m,         /* Centre of mass positions                   */
                mom,            /* " " " momenta                              */
                momp;           /* Predicted C of M momenta                   */
   quat_mp      quat,           /* Quaternions for this component             */
                amom,           /* Quaternion derivatives                     */
                amomp;          /* Predicted quaternion derivatives           */
   mat_mp       h,              /* Unit cell for zero-stress simulation       */
                hmom,           /* Unit cell derivatives                      */
                hmomp;          /* Predicted unit cell derivatives            */
                /*
                 *  Following variables ta..., ra.. have been introduced
                 *  by VVM and have dimensions of [nspecies]
                 */
   real        ts,             /* Nose-Poincare temperature co-ordinate      */
               tsmom;          /* Nose-Poincare temperature momentum         */
   real	       H_0;	       /* Original value of Nose Hamiltonian         */
   real        rs,             /* Reserved for future use.	 	     */
               rsmom;          /* Reserved for future use.	 	     */
} system_mt, *system_mp;


typedef struct                  /* Information for one species                */
{
   real         inertia[3],     /* Principal moments of inertia               */
                mass,           /* Mass of whole molecule                     */
                dipole,         /* Dipole Moment                              */
                charge;         /* Total charge                               */
   int          nsites,         /* Number of sites on this species            */
                nmols;          /* Number of molecules of this species        */
   int          rdof,           /* Rotational degrees of freedom (2=linear)   */
                framework;      /* Flag to signal this is a framework species */
   char         name[L_spec];   /* Name of this species                       */
   int          *site_id;       /* site identifier array                      */
   vec_mp       p_f_sites;      /* Site co-ordinates in principal frame       */
                /* Dynamic variable arrays for this species                   */
                /* These point to a subset of the whole-system arrays         */
                /* Dimensions for C of M quantities are [nmols][3]            */
                /* and for quaternions and derivatives, [nmols][4]            */
                /* If species is monatomic, quaternion pointers are null      */
   vec_mp       c_of_m,         /* Centre of mass positions                   */
                mom,            /* " " " momenta                              */
                momp;           /* Predicted C of M momenta                   */
   quat_mp      quat,           /* Quaternions for this species               */
                amom,           /* Quaternion derivatives                     */
                amomp;          /* Predicted quaternion derivatives           */
} spec_mt, *spec_mp;

typedef struct                  /* site info template.                        */
{
   double       mass,
                charge;
   char         name[L_site];
   int          flag;
   int          pad;
}       site_mt,  *site_mp;
                
typedef struct                  /* Holds potential parameter information      */
{
   int          flag;
   int          pad;
   real         p[NPOTP];
} pot_mt, *pot_mp;

typedef struct
{
   char *name;
   int  npar;
} pots_mt;

typedef struct                  /* Units used for program input               */
{
   double       m,              /* mass                                       */
                l,              /* length                                     */
                t,              /* time                                       */
                q;              /* charge                                     */

} unit_mt, *unit_mp;

typedef struct                  /* Record of dimensions of physical quantity  */
{
  int           m,              /* Number of powers of mass in unit           */
                l,
                t,
                q;
} dim_mt, *dim_mp;

typedef struct                          /* Struct template for keyword        */
{                                       /* in read_control.                   */
   char *key,
        *format,
        *defalt;
   gptr *ptr;
}       match_mt;

typedef struct                  /* Restart file header format                 */
{
   time_mt      timestamp,      /* Date and time restart file was written     */
                prev_timestamp; /* Timestamp of preceding restart file        */
   char         init_date[DLEN],/* Date run was initiated (propagated through)*/
                title[L_name],  /* Title when run was initiated               */
                vsn[L_vsn];     /* Version SID of program that wrote restart  */
   int          seq;            /* Sequence NO.  eg 5th restart in run        */
}       restrt_mt;

typedef struct                  /* Dump file header format                    */
{
   char         title[L_name],  /* Run title at beginning of dump run         */
                vsn[L_vsn];     /* RCS Revision number                        */
   long         istep,          /* Timestep at beginning of this file         */
                dump_interval;  /* How many steps between dumps               */
   int          dump_level,     /* Parameter determining contents of record   */
                maxdumps,       /* Maximum number of dump records in file     */
                ndumps,         /* How many dump records in file              */
                dump_size;      /* Size of a dump record                      */
   time_mt      timestamp,      /* Time file was written                      */
                dump_init,      /* Time dump run was started (ie first file)  */
                restart_timestamp;/* Time corresponding restart file written  */
   size_mt	sysinfo_size;
}       dump_mt;

typedef struct mol_mt		/* Per species information for dump header.   */
{
   float inertia[3];
   float mass;
   float dipole;
   float charge;
   int nmols;
   int rdof;
   int framework;
   char name[L_spec];
} mol_mt;

typedef struct dump_sysinfo_mt  /* Dump file system info format               */
{
   float      deltat;		/* Timestep*dump-interval		      */
   int	      nmols, nmols_r;   /* Number of molecules and rotations          */
   int	      nspecies;         /* Number of species in system.		      */
   struct     mol_mt mol[1];   /* Number of molecules per species         1   */
}       dump_sysinfo_mt;

#define MAX_ROLL_INTERVAL       100
typedef struct
{       
   int          nav, 
                nroll, 
                iroll, 
                pad;
   double align;
} av_head_mt;

#endif
