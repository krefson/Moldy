/******************************************************************************
 *  dump   Does a periodic dump of simulation data to a file.	In fact most  *
 *  of this function is concerned with the bookkeeping associated with 	      *
 *  maintaining a sequentially numbered set of contiguous dump files, 	      *
 *  ensuring contiguity and dealing with errors. The actual writing is done   *
 *  by 'write_dump_record'.
 * 									      *
 *    The user must supply a prototype filename (control.dump_file, read in   *
 *  by keyword "dump-name".  This should contain a numeric conversion code    *
 *  in 'printf' format (eg %3.3d ) which will be replaced by a sequential file*
 *  number.  (If it is absent the file number is added to the end).	      *
 *  Each dump file consists of a header record which is a struct (type dump_t)*
 *  followed by up to control.ndumps dump records.  A file may be incomplete- *
 *  the ndumps field of the header struct records the actual number.  It is   *
 *  not necessary to write one file per run - records are added to the end    *
 *  of an existing file and a new one is started only when it is full. The    *
 *  contents of the record are determined by the "dump-level" parameter and   *
 *  consist of binary data.  Each header struct contains a unique timestamp - *
 *  "dump_init" which is the same for all files in a sequence of dumps. (Note *
 *  that the previous dump file MUST be accessible in order to read and       *
 *  propagate this number - so when archiving dumps, always leave the very    *
 *  last one.)								      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	dump.c,v $
 * Revision 1.2  89/05/11  13:56:05  keith
 * Tidied up code using flags rather than GOTO's
 * Fixed bug in 'dump_convert' which wrote outside array bounds
 * Put PE in dump-level 1 rather than 8
 * 
 * Revision 1.1  89/04/27  15:13:56  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: dump.c,v 1.3 89/05/15 17:38:22 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#include	<stdio.h>
#include	<ctype.h>
#include 	"string.h"
/*========================== program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== Library declarations ============================*/
void	cfree();			/* Free allocated memory	      */
time_t	time();
/*========================== External function declarations ==================*/
static char	*mutate();
double		mdrand();
void		mat_vec_mul();
static void	dump_convert();
static void	real_to_float();
void		message();
void		note();
/*========================== External data references ========================*/
extern contr_t	control;
extern restrt_t restart_header;
/*========================== Macros ==========================================*/
#define DUMP_SIZE(level)	(((level & 1)+(level & 2) + \
			  (level & 4) ) * \
			 (3*system->nmols + 4*system->nmols_r + 9)+ \
			  (level & 8) * \
			 (3*system->nmols + 3*system->nmols_r + 9) +\
			  (level & 1))
/*============================================================================*/

void	dump(system, force, torque, stress, pe)
system_p	system;
vec_t		force[], torque[];
mat_t		stress;
double		pe;
{
   FILE		*dumpf;			/* File pointer to dump files	      */
   dump_t	dump_header;		/* Header record proforma	      */
   char		cur_file[L_name],	/* Names of current and previous      */
   		prev_file[L_name],	/* dump files.			      */
   		*fname;			/* Pointer to one of above filenames  */
   int		filenum=control.dump_offset+(control.istep-control.begin_dump) /
   		        (control.dump_interval * control.maxdumps),
		ndumps =(control.istep-control.begin_dump) /
   			 control.dump_interval % control.maxdumps; 
   int		dump_size = DUMP_SIZE(control.dump_level);
                                        /* Size in floats of each dump record */
   float	*dump_buf = aalloc(dump_size, float);	/* For converted data */
   long		file_pos,		/* Offset within file of current rec. */
   		file_len;		/* Length of file		      */
   		boolean errflg = false;
   static	boolean init = true;
   
   if( ! strchr(control.dump_file, '%') )
      	(void)strcat(control.dump_file, "%d");
   (void)sprintf(cur_file,  control.dump_file, filenum);
   (void)sprintf(prev_file, control.dump_file, filenum-1);

   if( control.istep > control.begin_dump )	/* Continue existing dump run */
   {						/* Get last dump's header     */
      if( ndumps == 0 )		fname = prev_file;
                        else    fname = cur_file;

      if( (dumpf = fopen(fname, "r+b")) == NULL)	/* Open dump file     */
         	errflg = true;

      if( !errflg &&
	   fread((char*)&dump_header, sizeof(dump_t), 1, dumpf) == 0 ) 
        	errflg = true;

      if( errflg )
	 message(NULLI, NULLP, WARNING, DUMPFI, fname);
      else if( control.dump_level != dump_header.dump_level )
      {
	 errflg = true;
         message(NULLI, NULLP, INFO, DMPALT);
      }

      if( !errflg )
      {
	 if( init &&
	    dump_header.restart_timestamp != restart_header.prev_timestamp )
	     message(NULLI, NULLP, FATAL, CONTIG, fname);

	 if( dump_header.ndumps < (ndumps ? ndumps : control.maxdumps ) )
            message(NULLI,NULLP,FATAL, SHTDMP,
		    fname, dump_header.ndumps, ndumps);
	 else if( dump_header.ndumps >  (ndumps ? ndumps : control.maxdumps) )
            message(NULLI,NULLP,WARNING,LNGDMP,
		    fname, dump_header.ndumps, ndumps);

	 dump_header.ndumps = ndumps;
      
	 (void)fseek(dumpf, 0L, SEEK_END);
	 file_len = ftell(dumpf);		/* Get length of file	      */
	 file_pos = sizeof(dump_t)
	            + ndumps*dump_size*sizeof(float);	/* Expected length    */
	 if( file_len < file_pos )
         	message(NULLI, NULLP, FATAL, CORUPT, fname, file_len, file_pos);
      }
      else
      {
	 (void)fclose(dumpf);
	 if( ndumps != 0 )
	 {
	    ndumps = 0;
	    control.begin_dump = control.istep;
	    control.dump_offset = ++filenum;
            (void)sprintf(cur_file, control.dump_file, filenum);
	 }
      }
   }
   if( control.istep == control.begin_dump )
   {
      (void)strcpy(dump_header.title, control.title);
      (void)strncpy(dump_header.vsn, "$Revision: 1.3 $"+11,
		                     sizeof dump_header.vsn-1);
      dump_header.dump_interval = control.dump_interval;
      dump_header.dump_level    = control.dump_level;
      dump_header.maxdumps	= control.maxdumps;
      dump_header.dump_size	= dump_size;
      dump_header.dump_init     = time((time_t *)0);

      while( (dumpf = fopen(cur_file, "r")) != 0 )
      {
	 (void)fclose(dumpf);
	 (void)mutate(control.dump_file);
	 message(NULLI, NULLP, WARNING, DMPEXS, cur_file, control.dump_file);
	 (void)sprintf(cur_file, control.dump_file, filenum);
      }
      note(DUMPST, cur_file, control.istep);
   }
      
   if( ndumps == 0)				/* Start of new dump file     */
   {
      if( dumpf ) (void)fclose(dumpf);		/* Finished with old file     */
      while( dumpf = fopen(cur_file, "r"))
      {
	 (void)fclose(dumpf);
	 (void)mutate(strcpy(prev_file, control.dump_file));
	 message(NULLI, NULLP, WARNING, DMPEXS, cur_file , prev_file);
	 (void)sprintf(cur_file, prev_file, filenum);
      }

      dump_header.ndumps = 0;	dump_header.istep = control.istep;
      dump_header.timestamp = time((time_t *)0);

      if( (dumpf = fopen(cur_file, "w+b")) == 0)
         	message(NULLI, NULLP, FATAL, DOFAIL, cur_file);
      if( fwrite((char*)&dump_header, sizeof(dump_t), 1, dumpf) == 0)
         	message(NULLI, NULLP, FATAL, DWFAIL, cur_file);
      file_pos = ftell(dumpf);
   }

   dump_convert(dump_buf, system, force, torque, stress, pe);
   dump_header.ndumps++;
   dump_header.restart_timestamp = restart_header.timestamp;

   (void)fseek(dumpf, file_pos, SEEK_SET);		/* Write data at end */
   if( fwrite((char*)dump_buf, dump_size, sizeof(float), dumpf) == 0 )
      	message(NULLI, NULLP, FATAL, DWFAIL, cur_file);

   (void)fseek(dumpf, 0L, SEEK_SET);			/* Write header      */
   if( fwrite((char*)&dump_header, sizeof(dump_t), 1, dumpf) == 0)
      	message(NULLI, NULLP, FATAL, DWFAIL, cur_file);

   (void)fclose(dumpf);
   cfree((char*)dump_buf);
   init = false;
}
/******************************************************************************
 *  mutate  Take a string defining a file name and randomly alter characters  *
 *  to generate another, hopefully valid file name.  Assumes that a printf    *
 *  conversion code is present and only alters contiguously alphabetic chars  *
 *  before the percent.							      *
 ******************************************************************************/
#define	N_MUTATE	3			/* Max number of mutated chars*/
static char	*mutate(name)
char	*name;
{
   char	*pc_pos = strchr(name,'%'),
   	*begin  = pc_pos;
   static char	alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZZ";

   while( begin > name && (pc_pos - begin) < N_MUTATE && isalnum(*(begin-1)) )
      begin--;

   while( begin < pc_pos )
      *begin++ = alpha[(int)(mdrand() * 26.0)];

   return(name);
}
/******************************************************************************
 *  dump_convert  format data and place in buffer for dumping to file         *
 *  Quantities dumped depend on bits set in 'control.dump_level' as follows:  *
 *	bit 0:	centre of mass co-ords, quaternions, unit cell matrix.	      *
 *	bit 1:	derivatives of above.					      *.
 *	bit 2:	second derivatives of above.				      *
 *	bit 3:	forces, torques, stress and potential energy.		      *
 *  All quantities are converted to floats (from real, whatever that is) and  *
 *  c_of_m and its derivatives are converted to unscaled co-ordinates.	      *
 ******************************************************************************/
static void	dump_convert(buf, system, force, torque, stress, pe)
float		*buf;
system_p	system;
vec_t		force[], torque[];
mat_t		stress;
double		pe;
{
   int		nmols   = system->nmols,
   		nmols_r = system->nmols_r;
   vec_t	*scale_buf = ralloc(nmols);
   real		ppe = pe;

   if( control.dump_level & 1)
   {
      mat_vec_mul(system->h, system->c_of_m, scale_buf, nmols);
      real_to_float(scale_buf[0],    buf, 3*nmols);	buf += 3*nmols;
      real_to_float(system->quat[0], buf, 4*nmols_r);	buf += 4*nmols_r;
      real_to_float(system->h[0],    buf, 9);		buf += 9;
      real_to_float(&ppe, buf, 1);		buf += 1;
   }
   if( control.dump_level & 2)
   {
      mat_vec_mul(system->h, system->vel, scale_buf, nmols);
      real_to_float(scale_buf[0],    buf, 3*nmols);	buf += 3*nmols;
      real_to_float(system->qdot[0], buf, 4*nmols_r);	buf += 4*nmols_r;
      real_to_float(system->hdot[0],    buf, 9);	buf += 9;
   }
   if( control.dump_level & 4)
   {
      mat_vec_mul(system->h, system->acc, scale_buf, nmols);
      real_to_float(scale_buf[0],     buf, 3*nmols);	buf += 3*nmols;
      real_to_float(system->qddot[0], buf, 4*nmols_r);	buf += 4*nmols_r;
      real_to_float(system->hddot[0], buf, 9);		buf += 9;
   }
   if( control.dump_level & 8)
   {
      real_to_float(force[0],  buf, 3*nmols);	buf += 3*nmols;
      real_to_float(torque[0], buf, 3*nmols_r);	buf += 3*nmols_r;
      real_to_float(stress[0], buf, 9);		buf += 9;
   }
   cfree((char*)scale_buf);
}
/******************************************************************************
 *  real_to_float  Copy data from one array to another, converting from       *
 *  typedef real to type float in the interests of space.  (may be null op'n) *
 ******************************************************************************/
static void	real_to_float(b, a, n)
float	a[];
real	b[];
int	n;
{
   int i;
   for( i=0; i<n; i++)
      a[i] = b[i];
}