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
 * Input	Functions for reading and verifying the input files (except   *
 *		the restart file). Contents:				      *
 * Strlower();		Convert a string to lowercase (Internal use only)     *
 * Get_line();		Read next input line.         (Internal use only)     *
 * Read_sysdef()       	Read the system specification file     	       	      *
 * Lattice_start()	Read initial crystal structure and set it up	      *
 * Read_control()       Read control file				      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	input.c,v $
 * Revision 2.0  93/03/15  14:49:08  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.8.1.21  93/03/12  12:14:14  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.8.1.21  93/03/09  15:58:36  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.8.1.20  92/10/28  14:09:48  keith
 * Changed "site_[tp]" typedefs to avoid name clash on HP.
 * 
 * Revision 1.8.1.19  92/09/22  14:48:15  keith
 * Tidied up calls to improve "lint" rating.
 * 
 * Revision 1.8.1.18  92/06/26  17:03:10  keith
 * Got rid of assumption that memory returned by talloc() or
 * arralloc() is zeroed.  This enhances ANSI compatibility.
 * Removed memory zeroing from alloc.c() in consequence.
 * 
 * Revision 1.8.1.17  92/03/11  12:56:21  keith
 * Changed "scale-separately" parameter to "scale options"
 * 
 * Revision 1.8.1.16  91/08/19  16:46:39  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * 
 * Revision 1.8.1.15  91/08/16  15:25:30  keith
 * Checked error returns from fread, fwrite, fseek and fclose more
 * rigourously.   Called strerror() to report errors.
 * 
 * Revision 1.8.1.13  91/03/12  15:42:49  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.8.1.12  90/10/22  17:47:25  keith
 * Corrected conversion of unit cell angles and lengths to vectors
 * in lattice_start().
 * 
 * Revision 1.8.1.11  90/08/20  17:25:45  keith
 * Modified to order species so that frameworks are last. 
 * 
 * Revision 1.8.1.10  90/05/03  16:41:24  keith
 * Fixed sys-spec parsing to cope with GEC and other broken scanf's which
 * return too many items matched.
 * 
 * Revision 1.8.1.9  90/04/16  18:18:16  keith
 * Added "strain-mask" field to input parse table.
 * 
 * Revision 1.8.1.8  90/04/12  16:27:09  keith
 * Added include of <stdio.h> which was removed from "structs.h"
 * 
 * Revision 1.8.1.7  90/03/26  18:05:50  keith
 * Moved system-dependant backup and temp file names to "defs.h"
 * 
 * Revision 1.8.1.6  89/11/22  14:34:44  keith
 * Changed default values of begin-rdf and average-interval.
 * 
 * Revision 1.8.1.5  89/11/21  16:32:30  keith
 * Removed member out_file from control and all uses. (Now command parameter).
 * 
 * Revision 1.8.1.4  89/11/20  18:06:23  keith
 * Modified form of match[] to include default values.
 * 
 * Revision 1.8.1.3  89/11/20  13:30:00  keith
 * Replaced separate arrays "types" and "npotp" with array of structs "potspec"
 * 
 * Revision 1.8.1.2  89/09/04  18:56:08  keith
 * Added 'surface-dipole' keyword to control file.
 * 
 * Revision 1.8.1.1  89/08/30  12:55:44  keith
 * Mods to add framework structures to simulation model
 * 
 * Revision 1.8  89/08/30  12:51:41  keith
 * Fixed read_sysdef() to keep original input line when reading potentials
 * to make error message informative.
 * Modified lattice_start() to fix bug which only considered rotations
 * of one species.  In conjunction with change in startup.
 * 
 * Revision 1.7  89/07/07  10:49:56  keith
 * Fixed lattice_start() so as not to test quaternion normalisation for
 * monatomic sopecies.
 * 
 * Revision 1.6  89/06/26  13:55:34  keith
 * Tidied up loops over species to use one pointer as counter.
 * Incorrect code to print control params removed from read_control()
 * 
 * Revision 1.5  89/06/21  13:36:42  keith
 * Moved pot. par. defs arrays types[], npotp[] and npott to kernel.c
 * Moved print_sysdef() to output.c
 * Made match[] external and its dimension nmatch into an external int
 * 
 * Revision 1.4  89/06/16  16:56:08  keith
 * Corrected bug in lattice_start() which crashed for point atoms/ions
 * Added message for successful lattice start
 * 
 * Revision 1.3  89/06/01  21:24:24  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.2  89/05/22  14:05:38  keith
 * Added rescale-separately option, changed 'contr_t' format.
 * 
 * Revision 1.1  89/04/20  16:00:42  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/input.c,v 2.0 93/03/15 14:49:08 keith Rel $";
#endif
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<ctype.h>
#include	<math.h>
#include 	"string.h"
#include	"stddef.h"
#include	<stdio.h>
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
void		q_mul_1();
void		message();
void    	zero_real();            /* Initialiser                        */
/*========================== External data references ========================*/
extern	contr_mt	control;		/* Main simulation control record     */
extern	unit_mt	input_unit;		/* Unit specification (see Convert.c) */
extern	pots_mt	potspec[];		/* Potential type specification       */
extern	int	npott;			/* Dimensions of above arrays	      */
/*========================== Macros ==========================================*/
#define		LLEN		132
		/* Flags to indicate status of potpar and site_info records   */
#define		S_USED		0x01
#define		S_MASS		0x02
#define		S_CHARGE	0x04
#define		S_NAME		0x08
/*========================== Control file keyword template ===================*/
					/* format SFORM is defined as %NAMLENs*/
					/* in structs.h, to avoid overflow */
/*
 *  Default backup and temporary file names if not set in "defs.h"
 */
#ifndef BACKUP_FILE
#define BACKUP_FILE	"MDBACKUP"
#endif
#ifndef TEMP_FILE
#define TEMP_FILE	"MDTEMPX"
#endif

match_mt	match[] = {
{"title",            SFORM,  "Test Simulation",(gptr*) control.title},
{"nsteps",           "%d",   "0",            (gptr*)&control.nsteps},
{"step",             "%lf",  "0.005",        (gptr*)&control.step},
{"text-mode-save",   "%d",   "0",            (gptr*)&control.print_sysdef},
{"new-sys-spec",     "%d",   "0",            (gptr*)&control.new_sysdef},
{"scale-options"   , "%d",   "0",            (gptr*)&control.scale_options},
{"surface-dipole",   "%d",   "0",            (gptr*)&control.surface_dipole},
{"lattice-start",    "%d",   "0",            (gptr*)&control.lattice_start},
{"sys-spec-file",    SFORM,  "",             (gptr*)control.sysdef},
{"restart-file",     SFORM,  "",             (gptr*)control.restart_file},
{"save-file",        SFORM,  "",             (gptr*)control.save_file},
{"dump-file",        SFORM,  "",             (gptr*)control.dump_file},
{"backup-file",      SFORM,  BACKUP_FILE,    (gptr*)control.backup_file},
{"temp-file",        SFORM,  TEMP_FILE,      (gptr*)control.temp_file},
{"strict-cutoff",    "%d",   "0",            (gptr*)&control.strict_cutoff},
{"xdr",    	     "%d",   "0",            (gptr*)&control.xdr_write},
{"strain-mask",	     "%d",   "200",	     (gptr*)&control.strain_mask},
{"nbins",            "%d",   "100",          (gptr*)&control.nbins},
{"seed",             "%d",   "1234567",      (gptr*)&control.seed},
{"page-width",       "%d",   "132",          (gptr*)&control.page_width},
{"page-length",      "%d",   "44",           (gptr*)&control.page_length},
{"scale-interval",   "%d",   "10",           (gptr*)&control.scale_interval},
{"const-pressure",   "%d",   "0",            (gptr*)&control.const_pressure},
{"reset-averages",   "%d",   "0",            (gptr*)&control.reset_averages},
{"scale-end",        "%d",   "1000000",      (gptr*)&control.scale_end},
{"begin-average",    "%d",   "1001",         (gptr*)&control.begin_average},
{"average-interval", "%d",   "5000",         (gptr*)&control.average_interval},
{"begin-dump",       "%d",   "1",            (gptr*)&control.begin_dump},
{"dump-interval",    "%d",   "20",           (gptr*)&control.dump_interval},
{"dump-level",       "%d",   "0",            (gptr*)&control.dump_level},
{"ndumps",           "%d",   "250",          (gptr*)&control.maxdumps},
{"backup-interval",  "%d",   "500",          (gptr*)&control.backup_interval},
{"roll-interval",    "%d",   "10",           (gptr*)&control.roll_interval},
{"print-interval",   "%d",   "10",           (gptr*)&control.print_interval},
{"begin-rdf",        "%d",   "1000000",      (gptr*)&control.begin_rdf},
{"rdf-interval",     "%d",   "20",           (gptr*)&control.rdf_interval},
{"rdf-out",          "%d",   "5000",         (gptr*)&control.rdf_out},
{"temperature",      "%lf",  "0.0",          (gptr*)&control.temp},
{"pressure",         "%lf",  "0.0",          (gptr*)&control.pressure},
{"w",                "%lf",  "100.0",        (gptr*)&control.pmass},
{"cutoff",           "%lf",  "10.0",         (gptr*)&control.cutoff},
{"subcell",          "%lf",  "0.0",          (gptr*)&control.subcell},
{"density",          "%lf",  "1.0",          (gptr*)&control.density},
{"alpha",            "%lf",  "0.3",          (gptr*)&control.alpha},
{"k-cutoff",         "%lf",  "2.0",          (gptr*)&control.k_cutoff},
{"rdf-limit",        "%lf",  "10.0",         (gptr*)&control.limit},
{"cpu-limit",        "%lf",  "1.0e20",       (gptr*)&control.cpu_limit},
{"mass-unit",        "%lf",  "1.6605655e-27",(gptr*)&input_unit.m},
{"length-unit",      "%lf",  "1.0e-10",      (gptr*)&input_unit.l},
{"time-unit",        "%lf",  "1.0e-13",      (gptr*)&input_unit.t},
{"charge-unit",      "%lf",  "1.6021892e-19",(gptr*)&input_unit.q}
	     };
int	nmatch=(sizeof match / sizeof(match_mt));
/*=============================================================================
 |   Start of functions							      |
 =============================================================================*/
/******************************************************************************
 *  get_line  read an input line skipping blank and comment lines	      *
 ******************************************************************************/
char	*get_line(line, len, file)
char	*line;
int	len;
FILE	*file;
{
   char	*s, *t;
   do
   {
      s = fgets(line, len, file);		/* Read one line of input     */
      if(s == NULL) break;			/* exit if end of file        */
      t = s + strlen(s) - 1;
      while(t >= s && (*t == ' ' || *t == '\t' || *t == '\n'))
         *t-- = '\0';				/* Strip trailing white space */
   }
   while(*s == '\0' || *s == '#');		/* Repeat if blank or comment */
   if(s == NULL)
      *line = '\0';				/* Return null at eof         */
   return(line);
}
/******************************************************************************
 * strlower   convert a string to lowercase anr return a pointer to it        *
 ******************************************************************************/
char	*strlower(s)
char	*s;
{
   char	*t;
   for(t = s; *t != '\0'; t++)
      *t = isupper(*t) ? tolower(*t) : *t;
   return(s);
}
/******************************************************************************
 *  Sort array of species structs so frameworks are at end.		      *
******************************************************************************/
void sort_species(species, nspecies)
spec_mt	*species;
int	nspecies;
{
   spec_mt tmp, *lo=species, *hi=species+nspecies-1;

   while( lo < hi )
   {
      while( lo < hi && ! lo->framework)
	 lo++;
      while( lo < hi && hi->framework)
	 hi--;

      if( lo >= hi )
	 break;

      tmp = *hi;
      *hi = *lo;
      *lo = tmp;

      lo++;
      hi--;
   }
}

/******************************************************************************
 *  read_sysdef    Read the system specification file which must be open and  *
 *  pointed to by parameter 'file'.  Set up the structures system and species *
 *  and arrays site_info and potpar (allocating space) and read in and check  *
 *  the supplied values.  The reading is done in two passes.  Pass 1 simply   *
 *  counts the number of species, number of sites on each species and the     *
 *  largest site identifier index in order to allocate the dynamic arrays.    *
 *  Pass 2 does the actual reading and checking.                              *
 ******************************************************************************/
void	read_sysdef(file, system, spec_pp, site_info, pot_ptr)
FILE		*file;			/* File pointer to read info from     */
system_mp	system;			/* Pointer to system array (in main)  */
spec_mp		*spec_pp;		/* Pointer to be set to species array */
site_mp		*site_info;		/* To be pointed at site_info array   */
pot_mp		*pot_ptr;		/* To be pointed at potpar array      */
{
   int		nspecies = 0,		/* Number of distinct species         */
   		max_id = 0,		/* Largest site identifier index      */
   		id, idi, idj,		/* Temp. site identifier index        */
   		isite,			/* species and site counters	      */
		sflag,			/* Temporary flag		      */
		i,			/* Counter			      */
		n_potpar,		/* Number of parameters for this pot'l*/
   		n_items;		/* How many items scanf found in input*/
   struct list_mt {int n; struct list_mt *p;};/* Template for linked list nsites*/
   struct list_mt nsites_base,		/* Head of list (contains no datum)   */
   		*nsites,		/* List entry for current species     */
                *last = &nsites_base;	/* List entry for previous species    */
   int		nerrs = 0;		/* Accumulated error count	      */
   int		flag;			/* Used to test 'fseek' result        */
   long		start_pos = ftell(file);/* Rewind marker for second pass      */
   char		name[LLEN],		/* Species name temporary             */
   		keywd[LLEN],		/* Species attribute keywords	      */
   		line[LLEN],		/* Store for input line from file     */
   		pline[LLEN];		/* Used in pot'l paramater parsing    */
   double	mass, charge, p_tmp;	/* Local temporaries		      */
   double	p_f_sites[3];		/* Local temporary		      */
   pot_mp	pp1;			/* Used for acces to potpar ij and ji */
   spec_mp	species, spec;		/* Local pointer to species array     */
   site_mp	s_ptr;			/* Local pointer to site info array   */
   static pot_mt	pot = {S_USED};	/* Local storage for potentials       */

   message(&nerrs,NULLP,INFO,SYSRD);
   /* First pass - read system definition and count nspecies, nsites, max_id  */
   (void)get_line(line,LLEN,file);		/* Read first line.	      */
   while(sscanf(line, "%s", name) > 0 && strcmp(strlower(name), "end") != 0)
   {						/* Loop, parsing 'line' for   */
      nspecies++;				/* name of new species.       */
      nsites = aalloc(1, struct list_mt); 	/* Make new list element      */
      last->p = nsites;				/* Link it in		      */
      last = nsites;				/* Backwards pointer for link */
      nsites->p=NULL; nsites->n=0;
      while(sscanf(get_line(line,LLEN,file), "%d", &id) > 0)
      {						/* Loop, reading and parsing  */
         nsites->n++;				/* for integer ie new site id.*/
         max_id = MAX(max_id, id);		/* Count nsites, greatest id. */
      }						/* Leave 'line' if parse fails*/
   }
   if(nspecies == 0)				/* Empty file??		      */
      message(&nerrs,NULLP,FATAL,NOSPEC);

   /* Allocate arrays of species and site info records */
   max_id++;
   system->max_id = max_id;
   *spec_pp    = aalloc(nspecies, spec_mt );
   *site_info  = aalloc(max_id, site_mt );
   (void)memset((char*)(*site_info), 0, max_id*sizeof(site_mt));
   *pot_ptr    = aalloc(SQR(max_id), pot_mt );
   for( i = 0; i < SQR(max_id); i++)
   {
      (*pot_ptr)[i].flag = 0;
      zero_real((*pot_ptr)[i].p,NPOTP);
   }
   species = *spec_pp;   			/* Local pointer for neatness.*/
   system->nspecies = nspecies;

   flag = fseek(file, start_pos, 0);		/* Prepare to reread input.   */
   if(flag)
      message(NULLI, NULLP, FATAL, SEFAIL, "control file", strerror(errno));
   nsites = &nsites_base;
   /* Pass 2.  read system definition and set up species and site_info arrays */
   for (spec = species; spec < species+system->nspecies; spec++)
   {						/* Loop over all species.     */
      n_items = sscanf(get_line(line,LLEN,file),"%s %d %s",
		       name, &spec->nmols, keywd);
      name[sizeof spec->name-1] = '\0';		/* Truncate before copying    */
      (void)strcpy(spec->name, name);		/* to avoid overflow.         */
      nsites = nsites->p;			/* Find next element of list  */
      spec->nsites = nsites->n;			/* which contains nsites.     */
      switch(n_items)
      {				   /* Relies on fall-through: do not re-order */
       case 3:
	 if(! strcmp(strlower(keywd), "framework"))
	    spec->framework = true;
	 else if(*keywd != '\0')   /* Kludge for broken scanf's.              */
	    message(&nerrs, NULLP, ERROR, UNKEY, keywd);
	 break;
       case 2:
	 spec->framework = false;
	 break;					/* Normal exit from switch    */
       default:
         message(&nerrs,line, ERROR, NONUM, name);
      }
      if(spec->nmols <= 0)			/* Can't have <=0 molecules   */
         message(&nerrs,line,ERROR, NOMOLS, spec->nmols, name);
      if(spec->nsites <=0)			/* or ghost molecules!        */
         message(&nerrs,NULLP,ERROR,NOSITE,spec->nsites,name);
      spec->p_f_sites = ralloc(spec->nsites);	/* Allocate space and set     */
      spec->site_id   = ialloc(spec->nsites+1);	/* pointers for each species. */

      for(isite = 0; isite < spec->nsites; isite++)
      {						/* Loop over sites on molecule*/
        n_items =sscanf(get_line(line,LLEN,file), "%d %lf %lf %lf %lf %lf %s",
                        &id,			/* Get and parse line of input*/
                        p_f_sites,
                        p_f_sites + 1,
                        p_f_sites + 2,
                        &mass,  &charge,  name);
        if(id <= 0)				/* Test for valid site index. */
        {
           message(&nerrs,line, ERROR, INVSID, id);
           id = 0;
        }
        spec->site_id[isite] = id;		/* Put id into rightful place.*/
        name[sizeof s_ptr->name - 1] = '\0';	/* Truncate site name.        */
        s_ptr = *site_info + id;		/* Reference (*site_info)[id].*/
        s_ptr->flag |= S_USED;
        switch (n_items)			/* Handle input items in      */
        {					/* reverse order.             */
           case 7:				/* Site name supplied.	      */
              if(s_ptr->flag & S_NAME && strcmp(name, s_ptr->name) != 0)
                 message(&nerrs,line,ERROR,NCONF,id,s_ptr->name);
              else
                 (void)strcpy(s_ptr->name, name);
              s_ptr->flag |= S_NAME;
           case 6:				/* Site charge supplied.      */
              if(s_ptr->flag & S_CHARGE && charge != s_ptr->charge)
                 message(&nerrs,line,ERROR,CCONF,id, s_ptr->charge);
              else
                 s_ptr->charge = charge;
              s_ptr->flag |= S_CHARGE;
           case 5:				/* Site mass supplied.        */
              if(s_ptr->flag & S_MASS && mass != s_ptr->mass)
                 message(&nerrs,line,ERROR,MCONF,id, s_ptr->mass);
              else if(mass < 0.0)
                 message(&nerrs,NULLP,ERROR,INVMAS,id,mass);
              else
                 s_ptr->mass = mass;
              s_ptr->flag |= S_MASS;
           case 4:				/* All site co-ordinates      */
	      for( i = 0; i < 3; i++ )
	         spec->p_f_sites[isite][i] = p_f_sites[i];
              break;				/* Normal exit from switch    */
           case 3:                      	/* One or more site 	      */
           case 2:                      	/* co-ordinates are missing.  */
           case 1:
              message(&nerrs,line,ERROR,MISSCO,n_items-1,id);
              break;
           default:				/* Should never occur.	      */
              message(&nerrs,NULLP, FATAL, BROKEN, n_items);
        }
      }
   }
   /*
    *  Order species structs with frameworks last
    */
   sort_species(species, nspecies);
   /*
    * Check that all sites have been fully defined, and for gaps in ordering.
    */
   for(id = 1; id < max_id; id++)
   {
      sflag = (*site_info)[id].flag;
      if(sflag & S_USED)
      {
         if(! (sflag & S_MASS))
            message(&nerrs,NULLP,ERROR,NOMASS,id);
         if(! (sflag & S_CHARGE))
            message(&nerrs,NULLP,ERROR,NOCGRG,id);
         if(! (sflag & S_NAME))
            message(&nerrs,NULLP,ERROR,NONAME,id);
      }
      else
         message(&nerrs,NULLP, WARNING, NOTUSD,id);
   }
   (void)get_line(line,LLEN,file);		/* read "end" -left by pass 1 */

   /* Next line is keyword indicating type of potentials to be used	      */
   n_items = sscanf(get_line(line,LLEN,file), "%s", name);
   for(i = 0; i < npott; i++)			/* Is 'name' a known type?    */
      if(strcmp(strlower(name), potspec[i].name) == 0)
         break;
   if(i == npott)				/* Did the loop find 'name'?  */
      message(&nerrs,line,FATAL,UNKPOT,name);	/* no			      */
   system->ptype = i;				/* yes		              */
   n_potpar = system->n_potpar = potspec[i].npar;
   						/* Now read in parameters     */
   while(sscanf(get_line(line,LLEN,file),"%s",name) > 0
                    && strcmp(strlower(name), "end") != 0)
   {
      n_items = 0;
      if(sscanf(line,"%d %d %[^#]",&idi,&idj,pline) <= 2)/* Not enough values */
         message(&nerrs,line,ERROR,NOPAIR);
      else
      {						/* Parse potential parameters */
	 (void)strcat(pline, "$");		/* Add marker to end	      */
         while(n_items < NPOTP && sscanf(pline,"%lf %[^#]", &p_tmp, pline) > 1 )
	    pot.p[n_items++] = p_tmp;
      }
      if (n_items < n_potpar)
         message(&nerrs,line,ERROR,NOPOTP,n_potpar);
      else				/* Test id pair and if OK store values*/
      {
         if(idi < 1 || idi >= max_id)
            message(&nerrs,line,ERROR,IDOUTR,idi);
         if(idj < 1 || idj >= max_id)
            message(&nerrs,line,ERROR,IDOUTR,idj);
         if(!(  (*site_info)[idi].flag & S_USED
              &&(*site_info)[idj].flag & S_USED))
            message(&nerrs,line,WARNING,EXTPOT);
         pp1 = (*pot_ptr) + (idi + idj * system->max_id);
         if(pp1->flag & S_USED)		/* pot'l for this id pair already set?*/
            message(&nerrs,line,ERROR,DUPPOT);
         else				/* Put values into pp1  and  pp2      */
         {
            (*pot_ptr)[idi + idj * system->max_id] = pot;
            (*pot_ptr)[idj + idi * system->max_id] = pot;
         }
      }
   }
   /* Check whether potentials have been defined for all 'used' site id's     */
   for(idi = 0; idi < max_id; idi++)
      for(idj = idi; idj < max_id; idj++)
         if( (   (*site_info)[idi].flag & S_USED    /* True if sites idi, idj */
              && (*site_info)[idj].flag & S_USED)   /* both used but pot'l not*/
            && !((*pot_ptr)[idi + max_id * idj].flag & S_USED))
            message(&nerrs,NULLP,WARNING,NOPOT,idi,idj);

   if(nerrs > 0)			/* if any errors have been detected   */
      message(&nerrs,NULLP,FATAL,ERRS,nerrs);
   else
      message(&nerrs,NULLP,INFO,SUCCES);
}
/******************************************************************************
 * lattice_start   Initialse the simulation co-ordinates on a lattice. Read   *
 * from the end of the system specification file.  The format is one line     *
 * specifying the unit cell (6 x floating point + 3 x int # cells)	      *
 *    a  b  c  alpha  beta  gamma  nx ny nz                                   *
 * followed by one line for each molecule in the unit cell:                   *
 *    species  x  y  z  q0  q1  q2  q3                                        *
 * 'species'  is the name,  x, y, z are FRACTIONAL co-ords and 4 quaternions. *
 ******************************************************************************/
void	lattice_start(file, system, species, qpf)
FILE	*file; 				/* File to read info from	      */
system_mp system;			/* System info struct		      */
spec_mp	species;			/* Array of species info structs      */
quat_mt	qpf[];				/* Princ frame rotation quaternion    */
{
   typedef struct init_s {int species;  struct init_s *next;
                  double r[3], q[4];} init_mt; 	/* For linked list of coords  */
   init_mt	*cur, *init = NULL;		/* Current and header of list */
   double	a, b, c, calpha, cbeta, cgamma;	/* Unit cell lengths, angles  */
   int		ix, iy, iz, nx, ny, nz;		/* Number of unit cells in MDC*/
   spec_mp	spec;
   char		line[LLEN], name[LLEN];
   int		n_items, nerrs = 0, ispec, imol, i;
   int		*nmols = ialloc(system->nspecies);
   real		ca, cb, cg, sg;
   quat_mt	q;

   (void)memset((char*)nmols,0,system->nspecies*sizeof(int));
   n_items = sscanf(get_line(line,LLEN,file),"%lf%lf%lf%lf%lf%lf%d%d%d",
		    &a, &b, &c, &calpha, &cbeta, &cgamma, &nx, &ny, &nz);
   if(n_items < 9)
      message(&nerrs, line, ERROR, NOCELL);
   if( ! (a > 0 && b > 0 && c > 0 && nx > 0 && ny > 0 && nz > 0 &&
	  calpha > 0 && calpha < 180.0 && cbeta > 0 && cbeta < 180.0 &&
	  cgamma > 0 && cgamma < 180.0))
      message(&nerrs, line,  ERROR, INVCEL);

   ca = cos(calpha*DTOR); cb = cos(cbeta*DTOR); cg = cos(cgamma*DTOR);
   sg = sin(cgamma*DTOR);

   system->h[0][0] = nx*a;			/* Set up MD cell matrix      */
   system->h[0][1] = ny*b * cg;			/* from lengths and angles.   */
   system->h[1][1] = ny*b * sg;
   system->h[0][2] = nz*c * cb;
   system->h[1][2] = nz*c / sg * (ca - cb*cg);
   system->h[2][2] = nz*c / sg * sqrt(1 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg);

   while(sscanf(get_line(line,LLEN,file), "%s", name) > 0 &&
	 strcmp(strlower(name), "end") != 0)	/* Cycle over lines in file   */
   {
      cur =aalloc(1, init_mt );			/* Get struct for new molecule*/
      cur->next = init;  init = cur;		/* Link it into list	      */
      n_items = sscanf(line, "%s%lf%lf%lf%lf%lf%lf%lf",
		       name, &cur->r[0], &cur->r[1], &cur->r[2],
		       &cur->q[0], &cur->q[1], &cur->q[2], &cur->q[3]);
      if(n_items > 1)				/* Have name of molecule      */
      {
	 for (spec = species; spec < species+system->nspecies; spec++)
	    if(strcmp(strlower(name),strlower(spec->name)) == 0)
	       break;
	 if(spec >= species+system->nspecies)	/* Didn't find it	      */
	    message(&nerrs,NULLP,ERROR,UNKSPE,name);
	 else					/* Found it - check values    */
	 {
	    cur->species = spec-species;
	    if(n_items < 4)
	       message(&nerrs, line, ERROR, FEWCOO, name);
	    if(cur->r[0] < 0 || cur->r[0] >= 1 ||
	       cur->r[1] < 0 || cur->r[1] >= 1 ||
	       cur->r[2] < 0 || cur->r[2] >= 1)
	       message(&nerrs,NULLP,ERROR,FRACCO,cur->r[0],cur->r[1],cur->r[2]);

	    if(species[cur->species].rdof != 0)
	    {
	       if( n_items < 8 )
		  message(&nerrs, line, ERROR, FEWQUA, name);
	       if(fabs(1.0 - SQR(cur->q[0]) - SQR(cur->q[1])
		           - SQR(cur->q[2]) - SQR(cur->q[3])) > 1e-4)
		  message(&nerrs,NULLP,ERROR,QNORM,
			  cur->q[0],cur->q[1],cur->q[2],cur->q[3]);
	    }
	    nmols[cur->species] += nx*ny*nz;
	 }
      }
   }
   for (spec = species; spec < species+system->nspecies; spec++)
   {
      ispec = spec-species;
      if(nmols[ispec] != spec->nmols)
         message(&nerrs,NULLP,ERROR,NIMOLS,spec->name,
		 nmols[ispec],spec->nmols);
      nmols[ispec] = 0;
   }

   if(nerrs > 0)				/* Is file all correct?       */
      message(NULLI, NULLP, FATAL, INITER, nerrs);

   for(cur = init; cur != NULL; cur = cur->next)
   {
      spec = species + cur->species;
      for(ix = 0; ix < nx; ix++)
         for(iy = 0; iy < ny; iy++)
            for(iz = 0; iz < nz; iz++)
	    {
	       imol = nmols[cur->species];
	       spec->c_of_m[imol][0] = (cur->r[0]+ix)/nx - 0.5;
	       spec->c_of_m[imol][1] = (cur->r[1]+iy)/ny - 0.5;
	       spec->c_of_m[imol][2] = (cur->r[2]+iz)/nz - 0.5;
	       if(spec->rdof > 0 )
	       {
		  for( i = 0; i < 4; i++ )
		     q[i] = cur->q[i];		/* Convert type to 'real'     */
		  q_mul_1(q, qpf[cur->species], spec->quat[imol]);
	       }
	       nmols[cur->species]++;
	    }
   }
   message(NULLI, NULLP, INFO, LATTIC);
}
/*******************************************************************************
 * assign()  Convert string value by format and assign to pointer location.    *
 ******************************************************************************/
int assign(strval, fmt, ptr)
char	*strval, *fmt;
gptr	*ptr;
{
   int len = strlen(fmt);
   int code = fmt[MAX(0,len-1)];
   if( len > 3 && fmt[len-2] == 'l' ) code = toupper(code);

   switch(code)
   {
    case 's':
    case ']':
      return sscanf(strval, fmt, (char*)ptr);
    case 'd':
      return sscanf(strval, fmt, (int*)ptr);
    case 'D':
      return sscanf(strval, fmt, (long*)ptr);
    case 'f':
      return sscanf(strval, fmt, (float*)ptr);
    case 'F':
      return sscanf(strval, fmt, (double*)ptr);
    default:
      message(NULLI, NULLP, FATAL,
	      "Scanf code \"%s\" not catered for", fmt);
   }
   return -1;		/* This statement is never reached		*/
}
/******************************************************************************
 *  read_control.   Read the control parameters from the standard input.      *
 *  Input lines are of the form " key = value ", one per line.  The struct    *
 *  'match' is searched for a matching key, and if found converts the value   *
 *  according to the format string and stores it at the value of the pointer  *
 *  in 'match'.	"name=" with no value means assign a null string. 	      *
 ******************************************************************************/
void	read_control(file)
FILE	*file;
{
   char		line[LLEN],
   		name[LLEN],
   		value[LLEN];
   int		i, n_items;
   int		nerrs = 0;

   while( *get_line(line,LLEN,file) != '\0' )
   {
      n_items = sscanf(line, " %[^= ] = %127[^#]",name, value);
      if(!strcmp(strlower(name),"end"))
         break;
      if( n_items < 1 )
         message(&nerrs,line,ERROR,NOVAL,name);
      else
      {
         for( i = 0; i < nmatch; i++ )		/* Search table for key	      */
            if( !strcmp(strlower(name), match[i].key) )
               break;				/* Found it          	      */
	 if( i == nmatch )			/* Reached end without success*/
            message(&nerrs,line,ERROR,NOTFND,name);
         else					/* Found it, so convert value */
         {
            if( n_items == 1 && strcmp(match[i].format,SFORM) == 0 )
                *(char*)match[i].ptr = '\0';	/* name=<empty> - assign null */
	    else
	    {
               n_items = assign(value, match[i].format, match[i].ptr);
               if( n_items < 1 )
                  message(&nerrs,NULLP,ERROR,BADVAL,value,name);
	    }
         }
      }
   }
   if( nerrs > 0 )
      message(&nerrs,NULLP,FATAL,ERRCON,nerrs,(nerrs>1)?'s':' ');
   else
      message(&nerrs,NULLP,INFO,SUCCON);
}
