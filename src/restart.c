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
static char *RCSid = "$Header: /home/tigger/keith/md/moldy/RCS/restart.c,v 1.6 89/11/01 17:40:36 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#include	<stdio.h>
#include 	"string.h"
#include	"time.h"
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External function declarations ==================*/
int		replace();
void		message();
void		note();
char		*av_ptr();
char		*atime();
char		*cctime();
/*========================== External data references ========================*/
extern  contr_t control;
extern	int	***rdf;				/* Accumulated RDF bins       */
/*========================== External data definitions =======================*/
restrt_t		restart_header = {0L, 0L, "", "", "0.0$", 0};
/*============================================================================*/
/******************************************************************************
 *   cread.   read a word from the binary restart file.  This should be the   *
 *   size of the next 'record' stored in the file.  Check it against the      *
 *   expected value.  Finally call fread to read in this data and check for   *
 *   error and end of file						      *
 ******************************************************************************/
void	cread(file, ptr, size, nitems)
FILE	*file;
char	*ptr;
int	size;
int	nitems;
{
   long	stored_size = 0;

   (void)fread((char*)&stored_size, sizeof stored_size, 1, file); 
   if(ferror(file))
      message(NULLI,NULLP,FATAL,REREAD,ferror(file),ftell(file));
   else if(feof(file))
      message(NULLI,NULLP,FATAL,REEOF);
   if(stored_size != (long)size * nitems)
      message(NULLI,NULLP,FATAL,REFORM, ftell(file),
              stored_size, (long)size * nitems);
   (void)fread(ptr, size, nitems, file);
   if(ferror(file))
      message(NULLI,NULLP,FATAL,REREAD,ferror(file),ftell(file));
   else if(feof(file))
      message(NULLI,NULLP,FATAL,REEOF);
}
/******************************************************************************
 *   cskip.  Skip over  the next record in the file.			      *
 ******************************************************************************/
void	cskip(file)
FILE	*file;
{
   long	stored_size = 0;

   (void)fread((char*)&stored_size, sizeof stored_size, 1, file); 
   if(ferror(file))
      message(NULLI,NULLP,FATAL,REREAD,ferror(file),ftell(file));
   else if(feof(file))
      message(NULLI,NULLP,FATAL,REEOF);
   (void)fseek(file, stored_size, 1);
   if(ferror(file))
      message(NULLI,NULLP,FATAL,REREAD,ferror(file),ftell(file));
   else if(feof(file))
      message(NULLI,NULLP,FATAL,REEOF);
}
/******************************************************************************
 *  cwrite.  opposite of cread.  write the length followed by the data	      *
 ******************************************************************************/
void	cwrite(file, ptr, size, nitems)
FILE	*file;
char	*ptr;
int	size;
int	nitems;
{
   long       length = (long)size*nitems;
   (void)fwrite((char*)&length, sizeof length, 1, file);
   (void)fwrite(ptr, size, nitems, file);
   if(ferror(file))
      message(NULLI,NULLP,FATAL,REWRT,ferror(file));
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
   cread(restart, (char*)header, sizeof(restrt_t), 1);
   cread(restart, (char*)contr, sizeof(contr_t), 1); 
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

   cread(file, (char*)system, sizeof(system_t), 1);/* Read in system structure*/

   /* Allocate space for species, site_info and potpar arrays and set pointers*/
   *spec_ptr  = aalloc(system->nspecies,                spec_t );
   *site_info = aalloc(system->max_id,                  site_t );
   *pot_ptr   = aalloc(system->max_id * system->max_id, pot_t );
   /*  read species array into allocated space				      */
   cread(file, (char*)*spec_ptr, sizeof(spec_t), system->nspecies);
   
   for (spec = *spec_ptr; spec < &(*spec_ptr)[system->nspecies]; spec++)
   {
      spec->p_f_sites = ralloc(spec->nsites);	/* Allocate the species -     */
      spec->site_id   = ialloc(spec->nsites);	/* specific arrays	      */
      cread(file, (char*)spec->p_f_sites, sizeof(vec_t), spec->nsites);
      cread(file, (char*)spec->site_id, sizeof(int), spec->nsites);
   }
   cread(file, (char*)*site_info, sizeof(site_t), system->max_id);
   cread(file, (char*)*pot_ptr, sizeof(pot_t), SQR(system->max_id));
}
/******************************************************************************
 *  read_restart.   read the dynamic simulation variables from restart file.  *
 ******************************************************************************/
void	read_restart(restart, system)
FILE		*restart;
system_p	system; 
{
   char		*ap;			/* Pointer to averages database       */
   int		asize;			/* Size of averages database	      */
   boolean	rdf_flag;		/* Indicates whether file contains rdf*/
   int		rdf_size = control.nbins*system->max_id*(system->max_id-1)/2;

   cread(restart, (char*)system->c_of_m, sizeof(vec_t), system->nmols);
   cread(restart, (char*)system->vel,    sizeof(vec_t), system->nmols);
   cread(restart, (char*)system->velp,   sizeof(vec_t), system->nmols);
   cread(restart, (char*)system->acc,    sizeof(vec_t), system->nmols);
   cread(restart, (char*)system->acco,   sizeof(vec_t), system->nmols);
   cread(restart, (char*)system->accvo,  sizeof(vec_t), system->nmols);
   if(system->nmols_r > 0)
   {
      cread(restart, (char*)system->quat,    sizeof(quat_t), system->nmols_r);
      cread(restart, (char*)system->qdot,    sizeof(quat_t), system->nmols_r);
      cread(restart, (char*)system->qdotp,   sizeof(quat_t), system->nmols_r);
      cread(restart, (char*)system->qddot,   sizeof(quat_t), system->nmols_r);
      cread(restart, (char*)system->qddoto,  sizeof(quat_t), system->nmols_r);
      cread(restart, (char*)system->qddotvo, sizeof(quat_t), system->nmols_r);
   }
   cread(restart, (char*)system->h,       sizeof(vec_t), 3);
   cread(restart, (char*)system->hdot,    sizeof(vec_t), 3);
   cread(restart, (char*)system->hdotp,   sizeof(vec_t), 3);
   cread(restart, (char*)system->hddot,   sizeof(vec_t), 3);
   cread(restart, (char*)system->hddoto,  sizeof(vec_t), 3);
   cread(restart, (char*)system->hddotvo, sizeof(vec_t), 3);

   ap = av_ptr(&asize);			/* get addr and size of database      */
   if( asize == 0 )			/* Don't read in any data	      */
      cskip(restart);
   else
      cread(restart, ap, asize, 1);

   cread(restart, (char*)&rdf_flag, sizeof rdf_flag, 1); /* Stored RDF data?  */
   if(rdf_flag && control.rdf_interval>0)/* Only read if data there and needed*/
      cread(restart, (char*)rdf[1][1], sizeof(int), rdf_size);
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
   char		*ap;			/* Pointer to averages database       */
   int		asize;			/* Size of averages database	      */
   int		rdf_size = control.nbins*system->max_id*(system->max_id-1)/2;
   int		zero = 0, one = 1;
   restrt_t	save_header;
   FILE		*save;

   save = fopen(control.temp_file, "wb");
   if(save == NULL)
   {
      message(NULLI, NULLP, ERROR, OSFAIL, control.temp_file);
      return;
   }

   (void)memcpy((char*)&save_header, (char*)&restart_header, sizeof(restrt_t));
   (void)strncpy(save_header.vsn, "$Revision: 1.6 $"+11,
		                  sizeof save_header.vsn-1);
   save_header.prev_timestamp = restart_header.timestamp;
   save_header.timestamp = time((time_t*)0);		/* Update header      */
   save_header.seq++;
   
   cwrite(save, (char*)&save_header, sizeof save_header, 1);
   cwrite(save, (char*)&control, sizeof control, 1); 

   cwrite(save, (char*)system, sizeof(system_t), 1);/* write system structure */
   cwrite(save, (char*)species, sizeof(spec_t), system->nspecies);
   
   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
      cwrite(save, (char*)spec->p_f_sites, sizeof(vec_t), spec->nsites);
      cwrite(save, (char*)spec->site_id, sizeof(int), spec->nsites);
   }
   cwrite(save, (char*)site_info, sizeof(site_t), system->max_id);
   cwrite(save, (char*)potpar, sizeof(pot_t), SQR(system->max_id));

   cwrite(save, (char*)system->c_of_m, sizeof(vec_t), system->nmols);
   cwrite(save, (char*)system->vel,    sizeof(vec_t), system->nmols);
   cwrite(save, (char*)system->velp,   sizeof(vec_t), system->nmols);
   cwrite(save, (char*)system->acc,    sizeof(vec_t), system->nmols);
   cwrite(save, (char*)system->acco,   sizeof(vec_t), system->nmols);
   cwrite(save, (char*)system->accvo,  sizeof(vec_t), system->nmols);
   if(system->nmols_r > 0)
   {
      cwrite(save, (char*)system->quat,    sizeof(quat_t), system->nmols_r);
      cwrite(save, (char*)system->qdot,    sizeof(quat_t), system->nmols_r);
      cwrite(save, (char*)system->qdotp,   sizeof(quat_t), system->nmols_r);
      cwrite(save, (char*)system->qddot,   sizeof(quat_t), system->nmols_r);
      cwrite(save, (char*)system->qddoto,  sizeof(quat_t), system->nmols_r);
      cwrite(save, (char*)system->qddotvo, sizeof(quat_t), system->nmols_r);
   }
   cwrite(save, (char*)system->h,       sizeof(vec_t), 3);
   cwrite(save, (char*)system->hdot,    sizeof(vec_t), 3);
   cwrite(save, (char*)system->hdotp,   sizeof(vec_t), 3);
   cwrite(save, (char*)system->hddot,   sizeof(vec_t), 3);
   cwrite(save, (char*)system->hddoto,  sizeof(vec_t), 3);
   cwrite(save, (char*)system->hddotvo, sizeof(vec_t), 3);

   ap = av_ptr(&asize);				/* get addr, size of database */
   cwrite(save, ap, asize, 1);
   
   if(control.rdf_interval > 0)			/* If we have rdf data	      */
   {
      cwrite(save, (char*)&one, sizeof(int), 1);/* Flag rdf data in file      */
      cwrite(save, (char*)rdf[1][1], sizeof(int), rdf_size);/* write data     */
   }
   else
      cwrite(save, (char*)&zero, sizeof(int), 1);/* Otherwise flag no rdf data*/
      
   (void)fclose(save);				/* Finished writing temp. file*/

   if(replace(control.temp_file, save_name) == 0)/* Rename temp file to output*/
   {						/* Don't signal backup write  */
      if(strcmp(save_name, control.backup_file))
         note("Configuration written to save file \"%s\" at %s, Seq no %d",
	      save_name, cctime(&save_header.timestamp), save_header.seq);
   }
   else
      message(NULLI, NULLP, ERROR, RNFAIL, save_name, control.temp_file);
}
