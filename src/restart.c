/******************************************************************************
 * restart	functions for writing program configurations to a data file   *
 *		and rereading them.  Because of the way they are used, the    *
 *		writing function is unitary, but the reading part is split    *
 *		between three routines.		Contents:		      *
 * cread()		Read data record from file, check length	      *
 * cwrite()		Write length and data as record to file		      *
 * re_re_header()	Read the restart header and control structs	      *
 * re_re_sysdef()	Read system definition records from restart file      *
 * read_restart()	Read simulation dynamic variables, averages & rdf's   *
 ******************************************************************************
 *      Revision Log
 *       $Log:	restart.c,v $
 * Revision 1.16  92/06/02  10:38:27  keith
 * Added check of ferror() as well as return from fwrite().  Talk
 * about belt and braces.
 * 
 * Revision 1.15  92/03/24  12:41:07  keith
 * Moved reset-averages code to values.c and did it properly.
 * 
 * Revision 1.14  92/03/19  15:14:12  keith
 * Added support for dynamic allocation of rolling average arrays,
 * conversion of existing restart files is done on fly.
 * 
 * Revision 1.13  91/08/23  14:12:46  keith
 * Fixed declaration of av_ptr to agree with definition in values.c
 * 
 * Revision 1.12  91/08/16  15:59:48  keith
 * Checked error returns from fread, fwrite, fseek and fclose more
 * rigourously.   Called strerror() to report errors.
 * 
 * Revision 1.11  91/08/15  18:12:14  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.10  91/03/12  15:43:16  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.9  91/02/19  14:51:31  keith
 * Minor changes to get rid of misleading compiler warnings.
 * 
 * Revision 1.8  90/08/24  17:47:26  keith
 * Made 'reset-averages' parameter actually do something.
 * 
 * Revision 1.7  90/05/02  15:28:55  keith
 * Removed references to size_t and time_t typedefs, no longer in "defs.h"
 * 
 * Revision 1.6  89/11/01  17:40:36  keith
 * Read_restart modified to allow skip of read of averages data -
 * specified by asize=0.
 * 
 * Revision 1.4  89/06/22  15:45:14  keith
 * Tidied up loops over species to use one pointer as counter.
 * 
 * Revision 1.3  89/05/22  14:05:43  keith
 * Added rescale-separately option, changed 'contr_t' format.
 * 
 * Revision 1.2  89/05/22  13:53:20  keith
 * Fixed to put correct RCS Revision into restart file header
 * 
 * Revision 1.1  89/04/27  15:16:09  keith
 * Initial revision
 * 
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/restart.c,v 1.16 92/06/02 10:38:27 keith Exp $";
#endif
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<stdio.h>
#include	"stddef.h"
#include 	"string.h"
#include	"time.h"
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
int		replace();
void		message();
void		note();
gptr		*av_ptr();
char		*atime();
char		*cctime();
/*========================== External data references ========================*/
extern  contr_t control;
extern	int	***rdf;				/* Accumulated RDF bins       */
/*========================== External data definitions =======================*/
restrt_t		restart_header = {0L, 0L, "", "", "0.0$", 0};
/*============================================================================*/
static size_t stored_size = 0;
static int    size_flg = 0;
/******************************************************************************
 *   cnext.  Read the size of the next 'record' stored in the file.           *
 *   Leave file pointer unchanged.  ie do lookahead.			      *
 ******************************************************************************/
size_t cnext(file)
FILE	*file;
{
   (void)fread((gptr*)&stored_size, sizeof stored_size, 1, file); 
   if(ferror(file))
      message(NULLI,NULLP,FATAL,REREAD,ftell(file),strerror(errno));
   else if(feof(file))
      message(NULLI,NULLP,FATAL,REEOF);
   size_flg++;
   return stored_size;
}
/******************************************************************************
 *   cread.   read a word from the binary restart file.  This should be the   *
 *   size of the next 'record' stored in the file.  Check it against the      *
 *   expected value.  Finally call fread to read in this data and check for   *
 *   error and end of file						      *
 ******************************************************************************/
void	cread(file, ptr, size, nitems)
FILE	*file;
gptr	*ptr;
size_t	size;
int	nitems;
{
   if( size_flg )
      size_flg = 0;
   else
   {
      (void)fread((gptr*)&stored_size, sizeof stored_size, 1, file); 
      if(ferror(file))
	 message(NULLI,NULLP,FATAL,REREAD,ftell(file),strerror(errno));
      else if(feof(file))
	 message(NULLI,NULLP,FATAL,REEOF);
   }
   if(stored_size != size * nitems)
      message(NULLI,NULLP,FATAL,REFORM, ftell(file),
              stored_size, size * nitems);
   (void)fread(ptr, size, nitems, file);
   if(ferror(file))
      message(NULLI,NULLP,FATAL,REREAD,ftell(file),strerror(errno));
   else if(feof(file))
      message(NULLI,NULLP,FATAL,REEOF);
}
/******************************************************************************
 *   cskip.  Skip over  the next record in the file.			      *
 ******************************************************************************/
void	cskip(file)
FILE	*file;
{
   if( size_flg )
      size_flg = 0;
   else
   {
      (void)fread((gptr*)&stored_size, sizeof stored_size, 1, file); 
      if(ferror(file))
	 message(NULLI,NULLP,FATAL,REREAD,ftell(file),strerror(errno));
      else if(feof(file))
	 message(NULLI,NULLP,FATAL,REEOF);
   }
   (void)fseek(file, (long)stored_size, SEEK_CUR);
   if(ferror(file))
      message(NULLI,NULLP,FATAL,REREAD,ftell(file),strerror(errno));
   else if(feof(file))
      message(NULLI,NULLP,FATAL,REEOF);
}
/******************************************************************************
 *  cwrite.  opposite of cread.  write the length followed by the data	      *
 ******************************************************************************/
void	cwrite(file, ptr, size, nitems)
FILE	*file;
gptr	*ptr;
size_t	size;
int	nitems;
{
   size_t       length = size*nitems;
   if( fwrite((gptr*)&length, sizeof length, 1, file) == 0 )
      message(NULLI,NULLP,FATAL,REWRT,strerror(errno));
   if( fwrite(ptr, size, nitems, file) < nitems )
      message(NULLI,NULLP,FATAL,REWRT,strerror(errno));
}
/******************************************************************************
 *  re_re_header   Read from the (already opened) restart file, the	      *
 *  restart_header and control structs.                                       *
 ******************************************************************************/
void	re_re_header(restart, header, contr)
FILE	*restart;
restrt_t *header;
contr_t	*contr;
{
   cread(restart,  (gptr*)header, sizeof(restrt_t), 1);
   cread(restart,  (gptr*)contr, sizeof(contr_t), 1); 
}
/******************************************************************************
 *  conv_potsize    Convert potential parameters array if NPOTP has changed   *
 ******************************************************************************/
void conv_potsize(pot, restart_size, system)
pot_t	*pot;
size_t	restart_size;
system_p system;
{
   size_t	old_pot_size;
   int		old_npotp, i;
   char		*tmp_pot;

   old_npotp = ((long)restart_size / SQR(system->max_id) 
		- (long)sizeof(pot_t)) / (long)sizeof(real) + NPOTP;
   old_pot_size = sizeof(pot_t) + (old_npotp - NPOTP) * sizeof(real);
   if( old_npotp == NPOTP )
      return;
   else if ( system->n_potpar <= NPOTP )
   {
      note("Old potential parameter array size = %d, new = %d",old_npotp,NPOTP);
      tmp_pot = aalloc( restart_size, char);
      memcpy(tmp_pot, (gptr*)pot, restart_size);
      for( i=0; i < SQR(system->max_id); i++)
	 memcpy((gptr*)(pot+i), tmp_pot+old_pot_size*i, 
		MIN(old_pot_size,sizeof(pot_t)));
      xfree(tmp_pot);
   }
   else
      message(NULLI, NULLP, FATAL, CPOTFL, NPOTP, system->n_potpar);
}
/******************************************************************************
 *  re_re_sysdef    Read the system specification from the restart file       *
 *  which must be open and pointed to by parameter 'file'.  Set up the        *
 *  structures system and species and arrays site_info and potpar (allocating *
 *  space) and read in the supplied values.  				      *
 ******************************************************************************/
void	re_re_sysdef(file, system, spec_ptr, site_info, pot_ptr)
system_p	system;			/* Pointer to system array (in main)  */
spec_p		*spec_ptr;		/* Pointer to be set to species array */
site_p		*site_info;		/* To be pointed at site_info array   */
pot_p		*pot_ptr;		/* To be pointed at potpar array      */
FILE		*file;			/* File pointer to read info from     */
{
   spec_p	spec;
   size_t	potsize;

   cread(file,  (gptr*)system, sizeof(system_t), 1);/* Read in system structure*/

   /* Allocate space for species, site_info and potpar arrays and set pointers*/
   *spec_ptr  = aalloc(system->nspecies,                spec_t );
   *site_info = aalloc(system->max_id,                  site_t );
   *pot_ptr   = aalloc(system->max_id * system->max_id, pot_t );
   /*  read species array into allocated space				      */
   cread(file,  (gptr*)*spec_ptr, sizeof(spec_t), system->nspecies);
   
   for (spec = *spec_ptr; spec < &(*spec_ptr)[system->nspecies]; spec++)
   {
      spec->p_f_sites = ralloc(spec->nsites);	/* Allocate the species -     */
      spec->site_id   = ialloc(spec->nsites);	/* specific arrays	      */
      cread(file,  (gptr*)spec->p_f_sites, sizeof(vec_t), spec->nsites);
      cread(file,  (gptr*)spec->site_id, sizeof(int), spec->nsites);
   }
   cread(file,  (gptr*)*site_info, sizeof(site_t), system->max_id);
   potsize = cnext(file);
   *pot_ptr = aalloc(MAX(potsize/sizeof(pot_t)+1,SQR(system->max_id) ), pot_t);
   cread(file,  (gptr*)*pot_ptr, potsize, 1);
   conv_potsize(*pot_ptr, potsize, system);
}
/******************************************************************************
 *  read_restart.   read the dynamic simulation variables from restart file.  *
 ******************************************************************************/
void	read_restart(restart, system)
FILE		*restart;
system_p	system; 
{
   gptr		*ap;			/* Pointer to averages database       */
   size_t	asize;			/* Size of averages database	      */
   boolean	rdf_flag;		/* Indicates whether file contains rdf*/
   int		rdf_size = control.nbins*system->max_id*(system->max_id-1)/2;

   cread(restart,  (gptr*)system->c_of_m, sizeof(vec_t), system->nmols);
   cread(restart,  (gptr*)system->vel,    sizeof(vec_t), system->nmols);
   cread(restart,  (gptr*)system->velp,   sizeof(vec_t), system->nmols);
   cread(restart,  (gptr*)system->acc,    sizeof(vec_t), system->nmols);
   cread(restart,  (gptr*)system->acco,   sizeof(vec_t), system->nmols);
   cread(restart,  (gptr*)system->accvo,  sizeof(vec_t), system->nmols);
   if(system->nmols_r > 0)
   {
      cread(restart,  (gptr*)system->quat,    sizeof(quat_t), system->nmols_r);
      cread(restart,  (gptr*)system->qdot,    sizeof(quat_t), system->nmols_r);
      cread(restart,  (gptr*)system->qdotp,   sizeof(quat_t), system->nmols_r);
      cread(restart,  (gptr*)system->qddot,   sizeof(quat_t), system->nmols_r);
      cread(restart,  (gptr*)system->qddoto,  sizeof(quat_t), system->nmols_r);
      cread(restart,  (gptr*)system->qddotvo, sizeof(quat_t), system->nmols_r);
   }
   cread(restart,  (gptr*)system->h,       sizeof(vec_t), 3);
   cread(restart,  (gptr*)system->hdot,    sizeof(vec_t), 3);
   cread(restart,  (gptr*)system->hdotp,   sizeof(vec_t), 3);
   cread(restart,  (gptr*)system->hddot,   sizeof(vec_t), 3);
   cread(restart,  (gptr*)system->hddoto,  sizeof(vec_t), 3);
   cread(restart,  (gptr*)system->hddotvo, sizeof(vec_t), 3);

   ap = av_ptr(&asize);			/* get addr and size of database      */
   if( asize == 0 )			/* Don't read in any data	      */
      cskip(restart);
   else
      cread(restart, ap, asize, 1);

   cread(restart,  (gptr*)&rdf_flag, sizeof rdf_flag, 1); /* Stored RDF data?  */
   if(rdf_flag && control.rdf_interval>0)/* Only read if data there and needed*/
      cread(restart,  (gptr*)rdf[1][1], sizeof(int), rdf_size);
}
/******************************************************************************
 *  write_restart.  Write the restart file.  Included (in order) are          *
 *  the 'restart_header' and 'control' structs, the system-specification      *
 *  (system species, site_info and potpar data) and the dynamic variables.    *
 ******************************************************************************/
void	write_restart(save_name, system, species, site_info, potpar)
char		*save_name;		/* Name of save file to be written    */
system_p	system;			/* Pointer to system array (in main)  */
spec_p		species;		/* Pointer to be set to species array */
site_p		site_info;		/* To be pointed at site_info array   */
pot_p		potpar;			/* To be pointed at potpar array      */
{
   spec_p	spec;
   gptr		*ap;			/* Pointer to averages database       */
   size_t	asize;			/* Size of averages database	      */
   int		rdf_size = control.nbins*system->max_id*(system->max_id-1)/2;
   int		zero = 0, one = 1;
   restrt_t	save_header;
   FILE		*save;
   char		*vsn = "$Revision: 1.16 $"+11;

   save = fopen(control.temp_file, "wb");
   if(save == NULL)
   {
      message(NULLI, NULLP, ERROR, OSFAIL, control.temp_file);
      return;
   }

   control.reset_averages = 0;		/* This flag never propagated.	      */
   save_header = restart_header;
   (void)strncpy(save_header.vsn, vsn, sizeof save_header.vsn-1);
   save_header.vsn[strlen(save_header.vsn)-2] = '\0'; /* Strip trailing $     */
   save_header.prev_timestamp = restart_header.timestamp;
   save_header.timestamp = time((time_t*)0);		/* Update header      */
   save_header.seq++;
   
   cwrite(save,  (gptr*)&save_header, sizeof save_header, 1);
   cwrite(save,  (gptr*)&control, sizeof control, 1); 

   cwrite(save,  (gptr*)system, sizeof(system_t), 1);
   cwrite(save,  (gptr*)species, sizeof(spec_t), system->nspecies);
   
   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
      cwrite(save,  (gptr*)spec->p_f_sites, sizeof(vec_t), spec->nsites);
      cwrite(save,  (gptr*)spec->site_id, sizeof(int), spec->nsites);
   }
   cwrite(save,  (gptr*)site_info, sizeof(site_t), system->max_id);
   cwrite(save,  (gptr*)potpar, sizeof(pot_t), SQR(system->max_id));

   cwrite(save,  (gptr*)system->c_of_m, sizeof(vec_t), system->nmols);
   cwrite(save,  (gptr*)system->vel,    sizeof(vec_t), system->nmols);
   cwrite(save,  (gptr*)system->velp,   sizeof(vec_t), system->nmols);
   cwrite(save,  (gptr*)system->acc,    sizeof(vec_t), system->nmols);
   cwrite(save,  (gptr*)system->acco,   sizeof(vec_t), system->nmols);
   cwrite(save,  (gptr*)system->accvo,  sizeof(vec_t), system->nmols);
   if(system->nmols_r > 0)
   {
      cwrite(save,  (gptr*)system->quat,    sizeof(quat_t), system->nmols_r);
      cwrite(save,  (gptr*)system->qdot,    sizeof(quat_t), system->nmols_r);
      cwrite(save,  (gptr*)system->qdotp,   sizeof(quat_t), system->nmols_r);
      cwrite(save,  (gptr*)system->qddot,   sizeof(quat_t), system->nmols_r);
      cwrite(save,  (gptr*)system->qddoto,  sizeof(quat_t), system->nmols_r);
      cwrite(save,  (gptr*)system->qddotvo, sizeof(quat_t), system->nmols_r);
   }
   cwrite(save,  (gptr*)system->h,       sizeof(vec_t), 3);
   cwrite(save,  (gptr*)system->hdot,    sizeof(vec_t), 3);
   cwrite(save,  (gptr*)system->hdotp,   sizeof(vec_t), 3);
   cwrite(save,  (gptr*)system->hddot,   sizeof(vec_t), 3);
   cwrite(save,  (gptr*)system->hddoto,  sizeof(vec_t), 3);
   cwrite(save,  (gptr*)system->hddotvo, sizeof(vec_t), 3);

   ap = av_ptr(&asize);				/* get addr, size of database */
   cwrite(save, ap, asize, 1);
   
   if(control.rdf_interval > 0)			/* If we have rdf data	      */
   {
      cwrite(save,  (gptr*)&one, sizeof(int), 1);/* Flag rdf data in file      */
      cwrite(save,  (gptr*)rdf[1][1], sizeof(int), rdf_size);/* write data     */
   }
   else
      cwrite(save,  (gptr*)&zero, sizeof(int), 1);/* Otherwise flag no rdf data*/
      
   if( ferror(save) || fclose(save) )            /* Delayed write error?       */
      message(NULLI, NULLP, FATAL, REWRT, strerror(errno));

   if(replace(control.temp_file, save_name) == 0)/* Rename temp file to output*/
   {						/* Don't signal backup write  */
      if(strcmp(save_name, control.backup_file))
         note("Configuration written to save file \"%s\" at %s, Seq no %d",
	      save_name, cctime(&save_header.timestamp), save_header.seq);
   }
   else
      message(NULLI, NULLP, ERROR, RNFAIL, save_name, control.temp_file,
	      strerror(errno));
}
