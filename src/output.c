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
 *      Revision Log
 *       $Log:	output.c,v $
 * Revision 1.3  89/06/01  21:25:07  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.2  89/05/24  13:55:03  keith
 * Changed ifdef's to select on __STDC__ macro
 * Message() now prints to user specified output file after initial set up
 * 
 * Revision 1.1  89/04/27  16:52:19  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: output.c,v 1.3 89/06/01 21:25:07 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#if ANSI || __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
char	*atime();			/* Current date and time in ASCII     */
char	*cctime();			/* Convert long time to ASCII.	      */
/*========================== External data references ========================*/
extern  contr_t 	control;
extern	restrt_t	restart_header;
extern	char		*types[];
extern	int		npotp[];
/*========================== External data definitions  ======================*/
int	out_page = 1;			/* Which page of output we are on     */
int	out_line = 999999;	        /* Which line of output               */
/*========================== Macros ==========================================*/
#define		S_USED		0x01
/******************************************************************************
 * new_line.   print a newline and update line counter                        *
 ******************************************************************************/
void	new_line()
{
   void	new_page();
   (void)putchar('\n');
   out_line++;
   if(out_line > control.page_length)   new_page();
}
void	new_lins(n)
int	n;
{
   while(n-- > 0)
      new_line();
}
/******************************************************************************
 * new_page   Take a new page on the output and print a header                * 
 ******************************************************************************/
void	new_page()
{
   (void)putchar('\f');					/* Take new page      */
   out_line = 0;					/* Print page header  */
   (void)printf("\t%s\t%s\tPage %d", atime(), control.title, out_page++);
   new_line();
}
/******************************************************************************
 *  Banner line.							      *
 ******************************************************************************/
void	put_line(c)
char	c;
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
#if ANSI || __STDC__
#undef  va_alist
#define	va_alist int *nerrs, ...
#define va_dcl /* */
#endif
/*VARARGS*/
void	message(va_alist)
va_dcl
{
   va_list	ap;
   char		*buff;
   int		sev;
   char		*format;
   static char	*sev_txt[] = {" *I* "," *W* "," *E* "," *F* "};
#if ANSI || __STDC__
   va_start(ap, nerrs);
#else
   int		*nerrs;

   va_start(ap);
   nerrs = va_arg(ap, int *);
#endif

   buff  = va_arg(ap, char *);
   sev   = va_arg(ap, int);
   format= va_arg(ap, char *);

   (void)printf(sev_txt[sev]);
   (void)vprintf(format, ap);
   va_end(ap);
   new_line();			/* To maintain pagination	      */

   if(buff != NULL)                     /* null ptr means don't print buffer  */
   {
      (void)printf("     buffer contents=\"%s\"",buff);
      new_line();
   }
   if(sev >= ERROR && nerrs != NULL)
      (*nerrs)++;
   if(sev == FATAL)
      exit(3);
}
/******************************************************************************
 *  note   write a message to the output file				      *
 ******************************************************************************/
#if ANSI || __STDC__
#undef  va_alist
#define	va_alist char *text, ...
#define va_dcl /* */
#endif
/*VARARGS*/
void	note(va_alist)
va_dcl
{
   va_list	ap;
#if ANSI || __STDC__
   va_start(ap, text);
#else
   char		*text;

   va_start(ap);
   text = va_arg(ap, char *);
#endif

   (void)printf(" *I* "); 
   (void)vprintf( text, ap);  new_line();
   va_end(ap);
}
/******************************************************************************
 *  Print_array    Print out an array of strings in a common format 	      *
 ******************************************************************************/
static void	print_array(text, n)
char	*text[];
int	n;
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
static void	format_int(text,value)
char	*text;
int	value;
{
   (void)printf("\t%-32s = %d",text,value);
   new_line();
}
/******************************************************************************
 *   Format_dbl     Print the name and value of some parameter in same format *
 ******************************************************************************/
static void	format_dbl(text,value,units)
char	*text;
double	value;
char	*units;
{
   (void)printf("\t%-32s = %g %s",text,value,units);
   new_line();
}
/******************************************************************************
 *   Format_vec     Print the name and value of some parameter in same format *
 ******************************************************************************/
static void	format_vec(text,value1,value2,value3,units)
char	*text;
double	value1,value2,value3;
char	*units;
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
				"REFSON@UK.AC.OX.VAX  (JANET)",
				"REFSON@VAX.OX.AC.UK  (BITNET)"};
static char	*copy_notice[] = {"Copyright Keith Refson 1988",
				  ""};

	
/******************************************************************************
 *  banner_page   Write the banner and relevant system/run information        *
 ******************************************************************************/
void	banner_page(system, species)
system_p	system;
spec_t	species[];
{
   int		ispec;
   spec_p	spec;
   mat_p	h = system->h;
   char		version[132], *vsn=version;

   new_page(); new_lins(2);
   print_array( banner, sizeof banner / sizeof(char*));
   (void)sprintf(version, "Version %.*s (%.*s) %.*s",
		 	(int)strlen(Revision+11)-1,      Revision+11,
		 	(int)strlen(Revision_State+8)-1, Revision_State+8,
		 	(int)strlen(Revision_Date+7)-1,  Revision_Date+7);
   print_array( &vsn, 1);
   print_array( name_addr, sizeof name_addr / sizeof(char*));
   print_array( copy_notice, sizeof copy_notice / sizeof(char*));
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

   for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
   {
      (void)printf(" %s", spec->name); new_line();
      format_int("Number of molecules",spec->nmols);
      format_int("Number of sites",spec->nsites);
      format_dbl("Mass",spec->mass,MUNIT_N);
      if(spec->rdof == 0)
      {
	 (void)printf(
	     "\t%s molecule has no rotational degrees of freedom", spec->name);
	 new_line();
      }
      else
      {
	 if(spec->rdof == 2)
	 {
	    (void)printf("\t%s molecule is linear",spec->name);
	    new_line();
	 }
	 format_vec("Moments of inertia",
		    spec->inertia[0],spec->inertia[1],spec->inertia[2],IUNIT_N);
	 format_dbl("Dipole moment",spec->dipole*CONV_D,CONV_D_N);
      }
   }
   new_line();
   (void)printf(" MD cell vectors"); new_line();
   format_vec("a",h[0][0],h[1][0],h[2][0],LUNIT_N);
   format_vec("b",h[0][1],h[1][1],h[2][1],LUNIT_N);
   format_vec("c",h[0][2],h[1][2],h[2][2],LUNIT_N);
   (void)printf(" Run parameters"); new_line();
   if(control.istep > 0)
      format_int("Initial step",control.istep); 
   format_int("Final step",control.nsteps);
   format_dbl("Size of step",control.step,TUNIT_N);
   format_dbl("CPU limit",control.cpu_limit,"s");
   if(control.scale_interval > 0)
   {
      (void)printf(" Velocities will be scaled");
      if( control.scale_separately )
	 (void)printf(" (for each species individually)");
      new_line();
      format_dbl("Applied Temperature",control.temp,"K");
      format_int("No. steps between scalings",control.scale_interval);
      format_int("End scaling at step",control.scale_end);
   }
   if(control.const_pressure)
   {
      (void)printf(" Constant stress ensemble will be used");
      new_line();
      format_dbl("Applied pressure", CONV_P*control.pressure,CONV_P_N);
      format_dbl("Mass parameter W",control.pmass,MUNIT_N);
   }
   format_dbl("Interaction cut-off",control.cutoff,LUNIT_N);
   if(control.alpha != 0.0)
   {
      format_dbl("Alpha parameter for Ewald sum",control.alpha,RLUNIT_N);
      format_dbl("Reciprocal space cut-off",control.k_cutoff,RLUNIT_N);
   }

   if( control.rdf_interval > 0 )
   {
      (void)printf(" Radial distribution functions will be calculated");
      new_line();
      format_int("Starting at timestep", control.begin_rdf);
      format_int("No. steps between binnings", control.rdf_interval);
      format_int("Calculate and print after", control.rdf_out);
   }

   if( control.dump_level > 0 && control.dump_interval > 0 )
   {
      (void)printf(" Configurational data will be dumped to file(s) %s",
		   control.dump_file);
      new_line();
      format_int("Starting at timestep", control.begin_dump);
      format_int("No. steps between dumps", control.dump_interval);
      format_int("Dump level", control.dump_level);
   }
   
   if(control.restart_file[0] == '\0')
   {
      (void)printf( " New run entitled \"%s\" started %s",
	      restart_header.title, restart_header.init_date);
      new_line();
   }
   else
   {
      (void)printf( " Run initialised from restart file %s written %s",
		    control.restart_file, cctime(&restart_header.timestamp));
      new_line();
      (void)printf( " This is restart No %d of run \"%s\" started %s",
	    restart_header.seq, restart_header.title, restart_header.init_date);
      new_line();
   }
   (void)fflush(stdout);
}
/******************************************************************************
 *  print sysdef   Print out the definition of the system, in the format that *
 *  read_sysdef can interpret.                                                *
 ******************************************************************************/
void    print_sysdef(system, species, site_info, potpar)
system_p        system;                 /* Pointer to system array (in main)  */
spec_t          species[];              /* Pointer to species array           */
site_p          site_info;              /* pointer to site_info array         */
pot_t           potpar[];               /* Potential parameter array          */
{
   spec_p       spec;
   int  ispec, isite, idi, idj, idij, ip;
   int  n_potpar = npotp[system->ptype];
   for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
   {
      (void)printf(" %-16s  %d\n", spec->name, spec->nmols);
      for(isite=0; isite < spec->nsites; isite++)
         (void)printf(" %6d %9g %9g %9g %9g %9g %s\n",
                        spec->site_id[isite],
                        spec->p_f_sites[isite][0],
                        spec->p_f_sites[isite][1],
                        spec->p_f_sites[isite][2],
                        site_info[spec->site_id[isite]].mass,
                        site_info[spec->site_id[isite]].charge,
                        site_info[spec->site_id[isite]].name);
   }
   (void)printf(" end\n");
   (void)printf(" %s potential parameters\n",types[system->ptype]);
   for(idi = 1; idi < system->max_id; idi++)
      for(idj = idi; idj < system->max_id; idj++)
      {
         idij = idj + idi * system->max_id;
         if(potpar[idij].flag & S_USED)
         {
            (void)printf(" %6d %6d", idi, idj);
            for(ip = 0; ip < n_potpar; ip++)
               (void)printf("%9g",potpar[idij].p[ip]);
            (void)printf("\n");
         }
      }
   (void)printf(" end\n");
}
 
