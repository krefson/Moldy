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
 * output	Contains various output and error handling functions, except  *
 *		for 'print_frame' and 'output' which, because of their	      *
 *		intimate connection with the averages/values database are     *
 *		located in "values.c".  Contents:			      *
 * new_line()		Write new line & manage page length		      *
 * new_page()		Start new output page				      *
 * put_line()		Write a line of symbols				      *
 * note()		Write a message to the output file		      *
 * message()		Write an error or warning message, possibly exiting   *
 * print_array()	\						      *
 * format_int()		 \   Internal (static) procedures for use by	      *
 * format_dbl()		 /		banner_page()			      *
 * format_vec()		/						      *
 * banner_page()	Write main startup banner and simulation parameters   *
 * print_sysdef()	Print system specification readable by read_sysdef()  *
 ******************************************************************************
 */
/*========================== Program include files ===========================*/
#include "defs.h"
/*========================== Library include files ===========================*/
#include 	<stdarg.h>
#include 	<math.h>
#include 	<stdlib.h>
#include 	<string.h>
#include        <stdio.h>
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
void	conv_potentials(const unit_mt *unit_from, const unit_mt *unit_to, 
			pot_mt *potpar, int npotpar, int ptype, 
			site_mt *site_info, int max_id);
void	conv_control(const unit_mt *unit, boolean direction);
					/* Unit conversion for 'control'      */
char	*atime(void);			/* Current date and time in ASCII     */
char	*cctime(time_mt *timeloc); /* Convert long time to ASCII.	      */
void	rmlockfiles(void);
void    par_abort(int);
/*========================== External data references ========================*/
extern	      contr_mt 	control;	    /* Main simulation control parms. */
extern	const match_mt	match[];	    /* Control file keyword table.    */
extern  const pots_mt	potspec[];	    /* Potential type specification   */
extern int ithread, nthreads;
/*========================== External data definitions  ======================*/
static  int	out_page = 1;		    /* Which page of output we are on */
static  int	out_line = 999999;	    /* Which line of output           */
/*========================== Macros ==========================================*/
#define		S_USED		0x01
/*========================== Special Control output cases ====================*/
static        int	one=1;
extern const  unit_mt	prog_unit;
extern	      unit_mt	input_unit;
static	const match_mt	special[] = {
        {"lattice-start",	"%d", "",	(gptr*)&one},
	{"restart-file",	"%s", "",	(gptr*)""},
	{"sys-spec-file",	"%s", "",	(gptr*)""},
	{"save-file",	        "%s", "",	(gptr*)""}
		      };
static	const int	nspecial = sizeof(special) / sizeof(match_mt);
/******************************************************************************
 * lines_left().  How many lines are left on page?			      *
 ******************************************************************************/
int lines_left(void)
{
   if( control.page_length > 0) 
      return MAX(0,control.page_length - out_line);
   else
      return 999999;
}
/******************************************************************************
 * new_line.   print a newline and update line counter                        *
 ******************************************************************************/
void	new_line(void)
{
   void	new_page(void);
   (void)putchar('\n');
   out_line++;
   if(out_line > control.page_length && control.page_length > 0)   new_page();
}
void	new_lins(int n)
{
   while(n-- > 0)
      new_line();
}
/******************************************************************************
 * new_page   Take a new page on the output and print a header                * 
 ******************************************************************************/
void	new_page(void)
{
   (void)putchar('\f');					/* Take new page      */
   out_line = 0;					/* Print page header  */
   (void)printf("\t%s\t%s\tPage %d", atime(), control.title, out_page++);
   new_line();
}
/******************************************************************************
 *  Banner line.							      *
 ******************************************************************************/
void	put_line(int c)
{
   int n = control.page_width;
   while(n-- > 0)
      (void)putchar(c);
   new_line();
}
/******************************************************************************
 *  message.   Deliver error message to possibly exiting.  It can be called   *
 *	       BEFORE output file is opened, in which case outt to stderr.    *
 ******************************************************************************/
/*VARARGS*/
void	message(int *nerrs, ...)

{
   va_list	ap;
   char		*buff;
   int		sev;
   char		*format;
   static char	*sev_txt[] = {" *I* "," *W* "," *E* "," *F* "};
   va_start(ap, nerrs);

   buff  = va_arg(ap, char *);
   sev   = va_arg(ap, int);
   format= va_arg(ap, char *);

   if( ithread == 0 || abs(sev) == FATAL)
   {
      (void)printf(sev_txt[abs(sev)]);
      (void)vprintf(format, ap);
      new_line();		      /* To maintain pagination	      */

      if(buff != 0)                /* null ptr means don't print buffer  */
      {
	 (void)printf("     buffer contents=\"%s\"",buff);
	 new_line();
      }
   }
   va_end(ap);
   if(sev >= ERROR && nerrs != 0)
      (*nerrs)++;
   if(sev == FATAL)
      rmlockfiles();
   if(abs(sev) == FATAL)
   {
     (void)fflush(stdout);
#ifdef SPMD
      par_abort(3);
#endif
      exit(3);
   }
}
/******************************************************************************
 *  note   write a message to the output file				      *
 ******************************************************************************/
/*VARARGS*/
void	note(char *text, ...)

{
   va_list	ap;
   va_start(ap, text);

   if( ithread > 0 )
      return;

   (void)printf(" *I* "); 
   (void)vprintf( text, ap);  new_line();
   va_end(ap);
}
/******************************************************************************
 *  Print_array    Print out an array of strings in a common format 	      *
 ******************************************************************************/
static void	print_array(char **text, size_mt n)
{
   int i;
   for(i=0; i<n; i++)
   {
      (void)printf("\t\t%s",text[i]);
      new_line();
   }
   new_line();
}
/******************************************************************************
 *   Format_int     Print the name and value of some parameter in same format *
 ******************************************************************************/
static void	format_int(char *text, int value)
{
   (void)printf("\t%-32s = %d",text,value);
   new_line();
}
/******************************************************************************
 *   Format_long     Print the name and value of some parameter in same format *
 ******************************************************************************/
static void	format_long(char *text, long int value)
{
   (void)printf("\t%-32s = %ld",text,value);
   new_line();
}
/******************************************************************************
 *   Format_dbl     Print the name and value of some parameter in same format *
 ******************************************************************************/
static void	format_dbl(char *text, double value, char *units)
{
   (void)printf("\t%-32s = %g %s",text,value,units);
   new_line();
}
/******************************************************************************
 *   Format_vec     Print the name and value of some parameter in same format *
 ******************************************************************************/
static void	format_vec(char *text, double value1, double value2, 
			   double value3, char *units)
{
   (void)printf("\t%-32s = %g %g %g %s",
		 text,value1,value2,value3,units);
   new_line();
}
/******************************************************************************
 *   Main banner, version string, name, address and copyright notice          *
 ******************************************************************************/
static char	*banner[] = {
		"#     # ####### #       ######  #     #",
		"##   ## #     # #       #     #  #   # ",
		"# # # # #     # #       #     #   # #  ",
		"#  #  # #     # #       #     #    #   ",
		"#     # #     # #       #     #    #   ",
		"#     # #     # #       #     #    #   ",
		"#     # ####### ####### ######     #   "};

static char	*Revision	= REVISION,
		*Revision_Date  = REVISION_DATE,
		*Revision_State = REVISION_STATE;

static char	*name_addr[] = {"Keith Refson",
				"Department of Earth Sciences",
				"Parks Road, Oxford OX1 3PR",
				"keith@earth.ox.ac.uk"};
static char	*copy_notice[] = {"Moldy Copyright (C) Keith Refson 1988, 1992, 1993",
				  "Moldy comes with ABSOLUTELY NO WARRANTY:",
				  "This is free software and you are welcome to",
				  "redistribute it under certain conditions.",
				  "For details see file COPYING included with source."};

	
/******************************************************************************
 *  banner_page   Write the banner and relevant system/run information        *
 ******************************************************************************/
void	banner_page(system_mp system, spec_mt *species, restrt_mt *restart_header)
{
   spec_mp	spec;
   mat_mp	h = system->h;
   real		chg;
   int rdof = 0;
   char		version[132], *vsn=version;

   new_page(); new_lins(2);
   print_array( banner, lsizeof banner / sizeof(char*));
   (void)sprintf(version, "Version %.*s (%.*s) %.*s",
		 	(int)strlen(Revision+7)-1,      Revision+7,
		 	(int)strlen(Revision_State+8)-1, Revision_State+8,
		 	(int)strlen(Revision_Date+7)-1,  Revision_Date+7);
   print_array( &vsn, (size_mt)1);
   print_array( name_addr, lsizeof name_addr / sizeof(char*));
   print_array( copy_notice, lsizeof copy_notice / sizeof(char*));
   if(control.restart_file[0] != '\0')
      if(control.new_sysdef)
         (void)printf( " New system specification read in from file %s",
		       control.sysdef);
      else
         (void)printf( " System specification read in from restart file %s",
		       control.restart_file);
   else
      (void)printf( " System specification read in from file %s",
		    control.sysdef);
   new_line();


   for(spec = species; spec < &species[system->nspecies]; spec++)
   {
      (void)printf(" %s", spec->name); new_line();
      if( spec->charge == 0.0 )
         format_int("Number of molecules",spec->nmols);
      else
         format_int("Number of ions",spec->nmols);
      format_int("Number of sites",spec->nsites);
      format_dbl("Mass",spec->mass,MUNIT_N);
      format_dbl("Electric charge", spec->charge*CONV_Q,CONV_Q_N);
      if(spec->nsites > 1 )
	 format_dbl("Dipole moment",spec->dipole*CONV_D,CONV_D_N);
      if(spec->rdof == 0)
      {
	 (void)printf(
	     "\t%s %s has no rotational degrees of freedom", spec->name, (spec->charge==0)?
                     "molecule":"ion");
	 new_line();
      }
      else
      {
	 rdof += spec->rdof;
	 if(spec->rdof == 2)
	 {
	    (void)printf("\t%s molecule is linear",spec->name);
	    new_line();
	 }
	 format_vec("Moments of inertia",
		    spec->inertia[0],spec->inertia[1],spec->inertia[2],IUNIT_N);
      }
   }
   new_line();
   (void)printf(" MD cell vectors"); new_line();
   format_vec("a",h[0][0],h[1][0],h[2][0],LUNIT_N);
   format_vec("b",h[0][1],h[1][1],h[2][1],LUNIT_N);
   format_vec("c",h[0][2],h[1][2],h[2][2],LUNIT_N);
   (void)printf(" Run parameters"); new_line();
   if(control.istep > 0)
      format_long("Initial step",control.istep); 
   format_long("Final step",control.nsteps);
   format_dbl("Size of step",control.step,TUNIT_N);
   format_dbl("CPU limit",control.cpu_limit,"s");
   if(control.scale_interval > 0 && control.istep < control.scale_end)
   {
      if( control.scale_options & 0x8 )
      {
	 (void)printf(" Velocities to be periodically RESET from MB distribution");
	 new_line();
      }
      else
      {
	 (void)printf(" Temperature will be scaled using %s kinetic energy",
		      control.scale_options & 0x4 ? 
		      "rolling average" : "instantaneous");
	 new_line();
	 if( control.scale_options & 0x3)
	 {
	    (void)printf(" (for ");
	    if( control.scale_options & 0x2 )
	    {
	       (void)printf("transl. and rotl.");
	       if( control.scale_options & 0x1 )
		  (void)printf(" and ");
	    }
	    if( control.scale_options & 0x1 )
	       (void)printf("each species");
	    
	    (void)printf(" individually)");
	    new_line();
	 }
      }
      format_long("No. steps between scalings",control.scale_interval);
      format_long("End scaling at step",control.scale_end);
   }
   if(control.const_temp)
   {
      (void)printf(" Nose-Poincare thermostat will be used");
      new_line();
      format_dbl("Temperature mass parameter ",
		 control.ttmass*CONV_TM, CONV_TM_N);
   }
   if((control.scale_interval > 0) || (control.const_temp != 0))
      format_dbl("Applied Temperature",control.temp,"K");

   if(control.const_pressure)
   {
      if(control.const_pressure == 1)
	 (void)printf(" Constant stress ensemble will be used");
      else if(control.const_pressure == 2)
	 (void)printf(" Constant pressure ensemble will be used");
      new_line();
      format_dbl("Applied pressure", CONV_P*control.pressure,CONV_P_N);
      format_dbl("Mass parameter W",control.pmass,MUNIT_N);
      if(control.const_pressure == 1)
	 format_int("h-matrix constraint mask",control.strain_mask);

      if(control.const_temp == 2)
         message(NULLI, NULLP, WARNING, GANDP);
   }
   (void)printf(" Parameters controlling accuracy of potential function evaluation");
   new_line();
   if( control.strict_cutoff )
      format_dbl("Interaction cut-off (strict)",control.cutoff,LUNIT_N);
   else
      format_dbl("Interaction cut-off (lazy)",control.cutoff,LUNIT_N);
   if(control.alpha != 0.0)
   {
      format_dbl("Ewald Sum accuracy", control.ewald_accuracy,"");
      format_dbl("Ewald Sum Alpha parameter",control.alpha,RLUNIT_N);
      format_dbl("Reciprocal space cut-off",control.k_cutoff,RLUNIT_N);
      if( control.surface_dipole )
      {
	 chg = 0.0;
	 for(spec = species; spec < &species[system->nspecies]; spec++)
	    chg += fabs(spec->charge);
	 
	 if( chg > 1.0e-6 )
	    message(NULLI, NULLP, WARNING,DPSCHG);
	 else
	    (void)printf("\tDe-Leeuw et al. surface dipole term included");
	 new_line();
      }
   }

   if ( rdof > 0 )
   {
     if( control.nosymmetric_rot )
       (void)printf(" Using General version of rotational leapfrog integrator");
     else
       (void)printf(" Using symmetry-adapted version of rotational leapfrog integrator");
     new_line();
   }

   if( control.rdf_interval > 0 && control.begin_rdf <= control.nsteps)
   {
      (void)printf(" Radial distribution functions will be calculated");
      new_line();
      format_dbl("Pair cutoff for RDF calculation",control.limit,LUNIT_N);
      format_long("Starting at timestep", control.begin_rdf);
      format_long("No. steps between binnings", control.rdf_interval);
      format_long("Calculate and print after", control.rdf_out);
   }

   if( control.dump_level > 0 && control.dump_interval > 0 )
   {
      (void)printf(" Configurational data will be dumped to file(s) %s",
		   control.dump_file);
      new_line();
      format_long("Starting at timestep", control.begin_dump);
      format_long("No. steps between dumps", control.dump_interval);
      format_int("Dump level", control.dump_level);
   }
   
   if(control.restart_file[0] == '\0')
   {
      (void)printf( " New run entitled \"%s\" started %s",
	      restart_header->title, restart_header->init_date);
      new_line();
   }
   else
   {
      (void)printf( " Run initialised from restart file %s written %s",
		    control.restart_file, cctime(&restart_header->timestamp));
      new_line();
      (void)printf( " This is restart No %d of run \"%s\" started %s",
	    restart_header->seq, restart_header->title, 
		                 restart_header->init_date);
      new_line();
   }
   (void)fflush(stdout);
}
/******************************************************************************
 *  print sysdef   Print out the definition of the system, in the format that *
 *  read_sysdef can interpret.                                                *
 ******************************************************************************/
static
void    print_sysdef(FILE *file,       
		     system_mp system,  /* Pointer to system array (in main)  */
		     spec_mt *species,  /* Pointer to species array           */
		     site_mp site_info, /* pointer to site_info array         */ 
		     pot_mt *potpar)    /* Potential parameter array          */
{
   spec_mp      spec;
   int  isite, idi, idj, idij, ip;
   int  n_potpar = potspec[system->ptype].npar;
   for(spec = species; spec < &species[system->nspecies]; spec++)
   {
      (void)fprintf(file, " %-16s  %d  %s\n", spec->name, spec->nmols,
		    spec->framework ? "framework" : "");
      for(isite=0; isite < spec->nsites; isite++)
         (void)fprintf(file, " %6d %12.12g %12.12g %12.12g %12.12g %12.12g %s\n",
                        spec->site_id[isite],
                        spec->p_f_sites[isite][0],
                        spec->p_f_sites[isite][1],
                        spec->p_f_sites[isite][2],
                        site_info[spec->site_id[isite]].mass,
                        site_info[spec->site_id[isite]].charge,
                        site_info[spec->site_id[isite]].name);
   }
   (void)fprintf(file, " end\n");
   (void)fprintf(file," %s potential parameters\n",potspec[system->ptype].name);
   for(idi = 1; idi < system->max_id; idi++)
      for(idj = idi; idj < system->max_id; idj++)
      {
         idij = idj + idi * system->max_id;
         if(potpar[idij].flag & S_USED)
         {
            (void)fprintf(file, " %6d %6d", idi, idj);
            for(ip = 0; ip < n_potpar; ip++)
               (void)fprintf(file, " %16.16g",potpar[idij].p[ip]);
            (void)fputc('\n',file);
         }
      }
   (void)fprintf(file, " end\n");
}
/******************************************************************************
 * Print_config()	Print out the configuration of the system in a text   *
 * format.  Control parameters system definition and 'lattice start' are      *
 * output allowing a portable restart.					      *
 ******************************************************************************/
void	print_config(char *save_name,   /* Name of save file to be written    */
		     system_mp system,  /* Pointer to system array (in main)  */ 
		     spec_mp species,   /* Pointer to be set to species array */
		     site_mp site_info, /* To be pointed at site_info array   */
		     pot_mp potpar)     /* To be pointed at potpar array      */
{
   FILE 	*out;
   const match_mt	*match_p, *cur, *special_p;
   spec_mp	spec;
   int		imol, code, i, j, k;
   double	cell_length[3], cell_angle[3];
   mat_mp	h = system->h;

   /*
    * Convert 'control' to input units for correct rereading.
    * Save current values and restore afterwards.
    */
   conv_control(&prog_unit, false);
   conv_potentials(&prog_unit, &input_unit, potpar, system->n_potpar,
                         system->ptype, site_info, system->max_id);

   if( (out = fopen(save_name, "w")) == 0 )
      message(NULLI, NULLP, FATAL, OSFAIL, save_name);

   for( match_p = match; match_p->key; match_p++)
   {
      for( special_p = special; special_p < &special[nspecial]; special_p++)
	 if( ! strcmp(match_p->key, special_p->key) )
	    break;

      if( special_p < &special[nspecial] )
	 cur = special_p;
      else
	 cur = match_p;

      code = cur->format[MAX(0, (long)strlen(cur->format)-1)];
      switch(code)
      {
       case 's':
       case ']':
	 (void)fprintf(out, "%s = %s\n", cur->key, (char*)cur->ptr);
	 break;
       case 'd':
	 (void)fprintf(out, "%s = %d\n", cur->key, *(int*)cur->ptr);
	 break;
       case 'f':
	 (void)fprintf(out, "%s = %.9g\n", cur->key, *(double*)cur->ptr);
	 break;
       default:
	 message(NULLI, NULLP, FATAL,
		 "Printf code \"%s\" not catered for", cur->format);
      }
   }
   (void)fprintf(out, "end\n");

   for( i = 0; i < 3; i++)
      cell_length[i] = sqrt(SQR(h[0][i]) + SQR(h[1][i]) + SQR(h[2][i]));
   for( i=0, j=1, k=2; i < 3; i++, j=(i+1)%3, k=(j+1)%3)
      cell_angle[i] = acos(
			 (h[0][j]*h[0][k] + h[1][j]*h[1][k] + h[2][j]*h[2][k])/
			 (cell_length[j]*cell_length[k])) / DTOR;
   
   print_sysdef(out, system, species, site_info, potpar);

   (void)fprintf(out, "%g %g %g %g %g %g 1 1 1\n",
		 cell_length[0], cell_length[1], cell_length[2],
		 cell_angle[0], cell_angle[1], cell_angle[2]);
   for( spec = species; spec < &species[system->nspecies]; spec++ )
      for( imol = 0; imol < spec->nmols; imol++ )
      {
	 (void)fprintf(out, "%s %g %g %g",
		       spec->name,                spec->c_of_m[imol][0]+0.5,
		       spec->c_of_m[imol][1]+0.5, spec->c_of_m[imol][2]+0.5);
	 if( spec->rdof > 0 )
	    (void)fprintf(out, " %g %g %g %g\n",
			  spec->quat[imol][0], spec->quat[imol][1],
			  spec->quat[imol][2], spec->quat[imol][3]);
	 else
	    (void)fputc('\n', out);
      }
   (void)fprintf(out, "end\n");
   
   if( ferror(out) || fclose(out) )
      message(NULLI,NULLP,FATAL,REWRT,strerror(errno));

   note("Configuration written to sys-spec plus lattice-start file \"%s\" ",
	save_name);

   conv_control(&prog_unit, true);
   conv_potentials(&input_unit, &prog_unit, potpar, system->n_potpar,
		   system->ptype, site_info, system->max_id);
}

	 
