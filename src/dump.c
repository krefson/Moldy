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
 *       $Log: dump.c,v $
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Got rid of all global (external) data items except for
 * "control" struct and constant data objects.  The latter
 * (pot_dim, potspec, prog_unit) are declared with CONST
 * qualifier macro which evaluates to "const" or nil
 * depending on ANSI/K&R environment.
 * Also moved as many "write" instantiations of "control"
 * members as possible to "startup", "main" leaving just
 * "dump".
 *
 * Declared as "static"  all functions which should be.
 *
 * Revision 2.5  1994/02/01  14:33:45  keith
 * Revised consistency checks to safeguard existing files:
 * Now checks timestep in dump header.
 * Reverted behaviour on dump-level change to start new run
 *  - acidentally botched in 2.4.
 * Changed behaviour on restart-file consistency check failure
 * to start new run too.
 * Improved messages.
 * 
 * Revision 2.4  1994/01/24  18:32:32  keith
 * Don't rename dump run if restarting from backup file and
 * dump file exists and matches.
 *
 * Revision 2.3  93/10/28  10:27:47  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.1  93/08/18  20:52:08  keith
 * Added support for dumps in XDR format.
 * 
 * Revision 2.0  93/03/15  14:49:01  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.19  93/03/09  15:58:26  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.18  92/06/02  10:37:52  keith
 * Added check of ferror() as well as return from fwrite().  Talk
 * about belt and braces.
 * 
 * Revision 1.17  92/04/21  17:49:16  keith
 * Corrected error message
 * 
 * Revision 1.16  92/03/02  15:28:32  keith
 * Fixed bug which caused dump/restart consistency check to erroneously
 * report an error on the second dump write.
 * Relaxed check to allow dump run to abort and restart from beginning.
 * 
 * Revision 1.15  91/08/16  15:23:57  keith
 * Checked error returns from fread, fwrite, fseek and fclose more
 * rigourously.   Called strerror() to report errors.
 * 
 * Revision 1.14  91/08/15  18:11:51  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.13  91/03/12  15:42:24  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.12  91/02/19  14:51:05  keith
 * Minor changes to get rid of misleading compiler warnings.
 * 
 * Revision 1.11  90/05/16  18:40:01  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.10  90/05/16  14:10:18  keith
 * Modified so that mutates always apply to rest of run.
 * 
 * Revision 1.9  90/05/02  15:28:52  keith
 * Removed references to size_t and time_t typedefs, no longer in "defs.h"
 * 
 * Revision 1.8  90/04/25  17:09:23  keith
 * Corrected test for mutation failure.
 * 
 * Revision 1.7  90/04/17  10:49:47  keith
 * Corrected test for dump and restart consistency.
 * 
 * Revision 1.6  90/04/09  14:49:33  keith
 * Now tests for failure of mutate() and gives up rather than lloping.
 * 
 * Revision 1.5  89/10/16  15:47:37  keith
 * Fixed DUMP_SIZE macro to be correct!
 * 
 * Revision 1.4  89/07/05  18:42:06  keith
 * Eliminated 'dump-offset' - renumbering starts at 1 for new dump run.
 * 
 * Revision 1.3  89/05/15  17:54:12  keith
 * Added members 'vsn' and 'dump_size' to header (cf structs.h r1.3)
 * Fixed bugs in call of write and ftell to write data correctly.
 * 
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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/dump.c,v 2.6 1994/02/17 16:38:16 keith Exp $";
#endif
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<ctype.h>
#include 	"string.h"
#include	"stddef.h"
#include	"time.h"
#include	<stdio.h>
/*========================== program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
#include	"xdr.h"
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
static char	*mutate();
double		mdrand();
void		mat_vec_mul();
static void	dump_convert();
static void	real_to_float();
#if defined(ANSI) || defined(__STDC__)
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
#else
void	note();				/* Write a message to the output file */
void	message();			/* Write a warning or error message   */
#endif
/*========================== External data references ========================*/
extern contr_mt	control;                    /* Main simulation control parms. */
#ifdef USE_XDR
static   XDR		xdrs;
#endif
/*========================== Macros ==========================================*/
#define DUMP_SIZE(level)  (( (level & 1) + (level>>1 & 1) + (level>>2 & 1) ) * \
			            (3*system->nmols + 4*system->nmols_r + 9)+ \
			     (level>>3 & 1) * \
			            (3*system->nmols + 3*system->nmols_r + 9) +\
			     (level & 1))
/*============================================================================*/
static
int read_dump_hdr(fname, dumpf, hdr_p, xdr_write)
char	*fname;
FILE	**dumpf;
dump_mt	*hdr_p;
boolean	*xdr_write;
{
   int      errflg = true;	/* Provisionally !!   */

   *xdr_write = false;
   if( (*dumpf = fopen(fname, "r+b")) == NULL)	/* Open dump file     */
      message(NULLI, NULLP, WARNING, DOERRR, fname, strerror(errno));
   else 
   {
#ifdef USE_XDR
      /*
       * Attempt to read dump header in XDR format
       */
      xdrstdio_create(&xdrs, *dumpf, XDR_DECODE);
      if( xdr_dump(&xdrs, hdr_p) )
      {
	 hdr_p->vsn[sizeof hdr_p->vsn - 1] = '\0';
	 if( strstr(hdr_p->vsn,"(XDR)") )
	 {
	    errflg = false;
	    *xdr_write = true;
	 }
      }
#endif
      /*
       * If we failed, try to read header as native struct image.
       */
      if( ! *xdr_write )
      {
	 if( fseek(*dumpf, 0L, 0) ) 
	    message(NULLI, NULLP, WARNING, SEFAIL, fname, strerror(errno));
	 else if( fread((gptr*)&*hdr_p, sizeof(dump_mt), 1, *dumpf) == 0 )
	    message(NULLI, NULLP, WARNING, DRERR, fname, strerror(errno));
	 else
	    errflg = false;
      }
   }
   return errflg;
}

/*============================================================================*/

void	dump(system, force, torque, stress, pe, restart_header, backup_restart)
restrt_mt	*restart_header;
int		backup_restart;
system_mp	system;
vec_mt		force[], torque[];
mat_mt		stress;
double		pe;
{
   FILE		*dumpf=NULL;		/* File pointer to dump files	      */
   dump_mt	dump_header,		/* Header record proforma	      */
   		hdr_tmp;
   char		cur_file[L_name],	/* Names of current and previous      */
   		prev_file[L_name],	/* dump files.			      */
   		*fname;			/* Pointer to one of above filenames  */
   int		filenum=(control.istep-control.begin_dump) /
   		        (control.dump_interval * control.maxdumps),
		ndumps =(control.istep-control.begin_dump) /
   			 control.dump_interval % control.maxdumps; 
   int		istep_hdr;		/* Timestep for dump header	      */
   unsigned     dump_size = DUMP_SIZE(control.dump_level);
                                        /* Size in floats of each dump record */
   float	*dump_buf=aalloc(dump_size,float);      /* For converted data */
   long		file_pos=0,		/* Offset within file of current rec. */
   		file_len;		/* Length of file		      */
   		boolean errflg = false;
#define		NMUTATES 10   		/* Max number of mutation attempts.   */
   int		nmutates = 0;   	/* Number of mutation attempts.	      */
   int		junk;
   boolean	xdr_write = false;	/* Is current dump in XDR format?     */
   static int	firsttime = 1;

   if( ! strchr(control.dump_file, '%') )
      	(void)strcat(control.dump_file, "%d");
   (void)sprintf(cur_file,  control.dump_file, filenum);
   (void)sprintf(prev_file, control.dump_file, filenum-1);

   if( control.istep > control.begin_dump )	/* Continue existing dump run */
   {						/* Get last dump's header     */
      if( ndumps == 0 )		fname = prev_file;
                        else    fname = cur_file;

      istep_hdr = control.istep - ndumps*control.dump_interval;
      if( ndumps == 0 )
	 istep_hdr -= control.maxdumps;
      errflg = true;
      /*
       * Attempt to read header and perform consistency checks which,
       * if failed, are "fatal" to dump run and initiate a *NEW* run
       * by setting errflg. 
       */
      if( read_dump_hdr(fname, &dumpf, &dump_header, &xdr_write) /* fails */ )
	 errflg = true;		                 /* message already printed */
      else if( control.dump_level != dump_header.dump_level )	/* Level    */
	 message(NULLI, NULLP, WARNING, DMPALT);
      else if( istep_hdr != dump_header.istep ) 		/* Timeslice */
	 message(NULLI, NULLP, WARNING, DUMPTS, fname, 
		 istep_hdr,dump_header.istep);
      else if( firsttime && 				  /* Matches Restart */
	    dump_header.timestamp < restart_header->timestamp &&
	    dump_header.restart_timestamp != restart_header->prev_timestamp &&
	    dump_header.restart_timestamp != restart_header->timestamp )
	 message(NULLI, NULLP, WARNING, CONTIG, fname);
      else
	 errflg = false;
      /*
       * At this point we think we have a matching dump header.  
       * Check to see if amount of data in file matches and whether
       * length of file is consistent with header.  Any errors now are
       * too serious to just start another dump sequence and so abort
       * the run. 
       */
      if( !errflg )
      {
	 if( dump_header.ndumps < (ndumps ? ndumps : control.maxdumps ) )
            message(NULLI,NULLP,FATAL, SHTDMP,
		    fname, dump_header.ndumps, ndumps);
	 else if( dump_header.ndumps >  (ndumps ? ndumps : control.maxdumps) )
            message(NULLI,NULLP,WARNING,LNGDMP,
		    fname, dump_header.ndumps, ndumps);

	 dump_header.ndumps = ndumps;
      
	 if( fseek(dumpf, 0L, SEEK_END) )
	    message(NULLI, NULLP, FATAL, SEFAIL, fname, strerror(errno));
	 file_len = ftell(dumpf);		/* Get length of file	      */
#ifdef USE_XDR
	 if( xdr_write )
	    file_pos = XDR_DUMP_SIZE
	       + ndumps*dump_size*XDR_FLOAT_SIZE;	/* Expected length    */
	 else
#endif
	    file_pos = sizeof(dump_mt)
	       + ndumps*dump_size*sizeof(float);	/* Expected length    */
	 if( file_len < file_pos )
         	message(NULLI, NULLP, FATAL, CORUPT, fname, file_pos, file_len);
      }
      else
      {
	 if( dumpf )
	 {
#ifdef USE_XDR
	    if( xdr_write )
	       xdr_destroy(&xdrs);
#endif
	    (void)fclose(dumpf);
	 }
	 if( ndumps != 0 )
	 {
	    ndumps = 0; filenum = 0;
	    control.begin_dump = control.istep;
            (void)sprintf(cur_file, control.dump_file, filenum);
	 }
	 message(NULLI,NULLP,WARNING,DRESET);
      }
   }
   if( errflg || control.istep == control.begin_dump )
   {
      (void)strcpy(dump_header.title, control.title);
      (void)strncpy(dump_header.vsn, "$Revision: 2.6 $"+11,
		                     sizeof dump_header.vsn-1);
#ifdef USE_XDR
      if( control.xdr_write )
      {
	 (void)strncat(dump_header.vsn, " (XDR)",
		                     sizeof dump_header.vsn-1);
	 xdr_write = true;
      }
#endif
      dump_header.dump_interval = control.dump_interval;
      dump_header.dump_level    = control.dump_level;
      dump_header.maxdumps	= control.maxdumps;
      dump_header.dump_size	= dump_size;
      dump_header.dump_init     = time((time_t *)0);

      note(DUMPST, cur_file, control.istep);
   }
      
   if( ndumps == 0)				/* Start of new dump file     */
   {
      if( dumpf ) 		/* Finished with old file     */
      {
#ifdef USE_XDR
	 if( xdr_write )
	    xdr_destroy(&xdrs);
#endif
	 (void)fclose(dumpf);
      }
      while( (dumpf = fopen(cur_file, "r")) != 0 )	/* File of that name */
      { 						/* already exists    */
	 (void)fclose(dumpf);
	 /*
	  *  Test whether existing file belongs to current dump
	  *  sequence in case of run restarted from backup.  In that
	  *  case, OVERWRITE it. 
          */
	 if( backup_restart )
	 {
	    errflg = read_dump_hdr(cur_file, &dumpf, &hdr_tmp, &junk);
	    if( errflg == false && (hdr_tmp.dump_init == dump_header.dump_init) )
	    {
	       (void)fclose(dumpf);
	       break;
	    }
	 }
	 
	 if( nmutates++ >= NMUTATES || mutate(control.dump_file) == NULL)
	    message(NULLI, NULLP, FATAL, MUFAIL, control.dump_file, nmutates);
	 message(NULLI, NULLP, WARNING, DMPEXS, cur_file, control.dump_file);
	 (void)sprintf(cur_file, control.dump_file, filenum);
      }

      dump_header.ndumps = 0;	dump_header.istep = control.istep;
      dump_header.timestamp = time((time_t *)0);
      dump_header.restart_timestamp = 0;

      if( (dumpf = fopen(cur_file, "w+b")) == 0)
	 message(NULLI, NULLP, FATAL, DOERRW, cur_file, strerror(errno));
#ifdef USE_XDR
      if( xdr_write )
      {
	 xdrstdio_create(&xdrs, dumpf, XDR_ENCODE);
	 if( ! xdr_dump(&xdrs, &dump_header) )
	    message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
	 if( (file_pos = xdr_getpos(&xdrs)) == -1)
	    message(NULLI, NULLP, FATAL, GPFAIL, cur_file, strerror(errno));
      }
      else
#endif
      {
	 if( fwrite((gptr*)&dump_header, sizeof(dump_mt), 1, dumpf) == 0)
	    message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
	 if( (file_pos = ftell(dumpf)) == -1)
	    message(NULLI, NULLP, FATAL, GPFAIL, cur_file, strerror(errno));
      }
   }
#ifdef USE_XDR
   else
   {
      if( xdr_write )
      {
	 xdr_destroy(&xdrs);
	 xdrstdio_create(&xdrs, dumpf, XDR_ENCODE);
      }
   }
#endif

   dump_convert(dump_buf, system, force, torque, stress, pe);
   dump_header.ndumps++;
   dump_header.restart_timestamp = restart_header->timestamp;

#ifdef USE_XDR
   if( xdr_write )
   {
      if( ! xdr_setpos(&xdrs, file_pos) )		/* Write data at end */
	 message(NULLI, NULLP, FATAL, SEFAIL, cur_file, strerror(errno));
      if( ! xdr_vector(&xdrs, (gptr*)dump_buf, dump_size, sizeof(float), 
		     xdr_float) )
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
   }
   else
#endif
   {
      if( fseek(dumpf, file_pos, SEEK_SET) )		/* Write data at end */
	 message(NULLI, NULLP, FATAL, SEFAIL, cur_file, strerror(errno));
      if( fwrite((gptr*)dump_buf, sizeof(float), dump_size, dumpf) < dump_size )
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
   }

#ifdef USE_XDR
   if( xdr_write )
   {
      (void)xdr_setpos(&xdrs, 0);
      if( ! xdr_dump(&xdrs, &dump_header) )
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
      xdr_destroy(&xdrs);
   }
   else
#endif
   {
      (void)fseek(dumpf, 0L, SEEK_SET);			/* Write header      */
      if( fwrite((gptr*)&dump_header, sizeof(dump_mt), 1, dumpf) == 0)
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
   }
   if( ferror(dumpf) || fclose(dumpf) )
      message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));

   firsttime = 0;
   xfree(dump_buf);
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

   if(begin >= pc_pos)				/* No characters to mutate    */
      return NULL;

   while( begin < pc_pos )
      *begin++ = alpha[(int)(mdrand() * 26.0)];

   return name;
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
system_mp	system;
vec_mt		force[], torque[];
mat_mt		stress;
double		pe;
{
   int		nmols   = system->nmols,
   		nmols_r = system->nmols_r;
   vec_mt	*scale_buf = ralloc(nmols);
   real		ppe = pe;

   if( control.dump_level & 1)
   {
      mat_vec_mul(system->h, system->c_of_m, scale_buf, nmols);
      real_to_float(scale_buf[0],    buf, 3*nmols);	buf += 3*nmols;
      if( system->nmols_r > 0 )
      {
	 real_to_float(system->quat[0], buf, 4*nmols_r);	
	 buf += 4*nmols_r;
      }
      real_to_float(system->h[0],    buf, 9);		buf += 9;
      real_to_float(&ppe, buf, 1);		buf += 1;
   }
   if( control.dump_level & 2)
   {
      mat_vec_mul(system->h, system->vel, scale_buf, nmols);
      real_to_float(scale_buf[0],    buf, 3*nmols);	buf += 3*nmols;
      if( system->nmols_r > 0 )
      {
	 real_to_float(system->qdot[0], buf, 4*nmols_r);
	 buf += 4*nmols_r;
      }
      real_to_float(system->hdot[0],    buf, 9);	buf += 9;
   }
   if( control.dump_level & 4)
   {
      mat_vec_mul(system->h, system->acc, scale_buf, nmols);
      real_to_float(scale_buf[0],     buf, 3*nmols);	buf += 3*nmols;
      if( system->nmols_r > 0 )
      {
	 real_to_float(system->qddot[0], buf, 4*nmols_r);
	 buf += 4*nmols_r;
      }
      real_to_float(system->hddot[0], buf, 9);		buf += 9;
   }
   if( control.dump_level & 8)
   {
      real_to_float(force[0],  buf, 3*nmols);	buf += 3*nmols;
      if( system->nmols_r > 0 )
      {
	 real_to_float(torque[0], buf, 3*nmols_r);
	 buf += 3*nmols_r;
      }
      real_to_float(stress[0], buf, 9);		buf += 9;
   }
   xfree(scale_buf);
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
