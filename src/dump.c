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
 */
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<ctype.h>
#include 	<string.h>
#include	<stddef.h>
#include	<time.h>
#include	<stdio.h>
/*========================== program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
#include	"xdr.h"
/*========================== External function declarations ==================*/
gptr            *talloc(int n, size_mt size, int line, char *file);
	       /* Interface to memory allocator       */
void            tfree(gptr *p);	       /* Free allocated memory	      	      */
static char	*mutate(char *name);
double		mdrand1(void);
void		mat_vec_mul(real (*m)[3], vec_mp , vec_mp , int );
void		mat_sca_mul(real s, mat_mt a, mat_mt b); 
void		invert(mat_mt , mat_mt );	/* Matrix inverter             */
void		transpose(mat_mt a, mat_mt b); /* transpose a 3x3 matrix       */
void            mvaxpy(int n, mat_mt a, vec_mt (*x), vec_mt (*y));
void            vscale( int,  double,  real *, int); /* Vector* const multiply */
static void	dump_convert(float *, system_mt (*), spec_mt (*), vec_mt (*), 
			     vec_mt (*), mat_mt, double );
static void	real_to_float(real *b, float *a, int n);
static void	amom_to_float(spec_mt *species, float *a, int n, double ts);
void		note(char *, ...);	/* Write a message to the output file */
void		message(int *, ...);	/* Write a warning or error message   */
/*========================== External data references ========================*/
extern contr_mt	control;                    /* Main simulation control parms. */
#ifdef USE_XDR
static   XDR		xdrs;
#endif
/*============================================================================*/
#ifndef INERTIA_MIN
#define INERTIA_MIN	1.0e-14		/* Tolerance for zero mom of I	      */
#endif
/*============================================================================*/
FILE  *open_dump(char *fname, char *mode)
{
   FILE *dumpf;

   dumpf = fopen(fname, mode);
   
#ifdef USE_XDR
   if( dumpf )
   {
      if( mode[0] == 'w' || (mode[0] && mode[1] == '+') ||  (mode[1] && mode[2] == '+'))
	 xdrstdio_create(&xdrs, dumpf, XDR_ENCODE);
      else
	 xdrstdio_create(&xdrs, dumpf, XDR_DECODE);
   }
#endif
    return dumpf;
}

int close_dump(FILE *dumpf)
{
#ifdef USE_XDR
   xdr_destroy(&xdrs);
#endif
   return fclose(dumpf);
}

/*ARGSUSED3*/
size_mt	dump_curpos(size_mt sysinfo_size, int dump_size, 
		    int ndumps, int nspecies,  boolean xdr_write)
{
#ifdef USE_XDR
   if( xdr_write )
      return XDR_DUMP_SIZE + XDR_SYSINFO_SIZE(nspecies) 
	                   + ndumps*dump_size*XDR_FLOAT_SIZE;
   else
#endif
      return sizeof(dump_mt) + sysinfo_size + ndumps*dump_size*sizeof(float);
}

/*ARGSUSED2*/
static
void dump_setpos(FILE *dumpf, size_mt file_pos, boolean xdr_write)
{
#ifdef USE_XDR
   if( xdr_write )
   {
      if( ! xdr_setpos(&xdrs, file_pos) )		/* Write data at end */
	 message(NULLI, NULLP, FATAL, SEFAIL, "", strerror(errno));
   }
   else
#endif
   {
      if( fseek(dumpf, file_pos, SEEK_SET) )		/* Write data at end */
	 message(NULLI, NULLP, FATAL, SEFAIL, "", strerror(errno));
   }
}


static
int read_dump_header(char *fname, FILE *dumpf, dump_mt *hdr_p, boolean *xdr_write,
		     int sysinfo_size, dump_sysinfo_mt *dump_sysinfo)
{
   int      errflg = true;	/* Provisionally !!   */
#ifdef USE_XDR
   char     vbuf[sizeof hdr_p->vsn + 1];
#endif
   int	    vmajor,vminor;

   *xdr_write = false;
#ifdef USE_XDR
   /*
       * Attempt to read dump header in XDR format
       */
   if( xdr_dump(&xdrs, hdr_p) )
   {
      strncpy(vbuf,hdr_p->vsn,sizeof hdr_p->vsn);
      vbuf[sizeof hdr_p->vsn] = '\0';
      if( strstr(vbuf,"(XDR)") )
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
      if( fseek(dumpf, 0L, 0) ) 
	 message(NULLI, NULLP, WARNING, SEFAIL, fname, strerror(errno));
      else if( fread((gptr*)&*hdr_p, sizeof(dump_mt), 1, dumpf) == 0 )
	 message(NULLI, NULLP, WARNING, DRERR, fname, strerror(errno));
      else
	 errflg = false;
   }
   if( ! errflg )
   {
      /*
       * Parse header version
       */
      errflg = true;
      if( sscanf(hdr_p->vsn, "%d.%d", &vmajor, &vminor) < 2 )
	 message(NULLI, NULLP, WARNING, INDVSN, hdr_p->vsn);
      if( vmajor < 2 || vminor <= 22)
	 message(NULLI, NULLP, WARNING, OLDVSN, hdr_p->vsn);
      else
	 errflg = false;
   }
   if( ! errflg && dump_sysinfo )
   {
#if 0
      if( hdr_p->sysinfo_size != sysinfo_size )
	 message(NULLI, NULLP, WARNING, CORUPT, fname, sysinfo_size, 
		 hdr_p->sysinfo_size);
#endif
      /*
       * Now check for sysinfo and read it
       */
#ifdef USE_XDR
      if( *xdr_write ) {
	 if( ! xdr_dump_sysinfo(&xdrs, dump_sysinfo, vmajor, vminor) )
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
	 errflg = false;
      } else
#endif
      {
	 if( fread((gptr*)dump_sysinfo, sysinfo_size, 1, dumpf) == 0)
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
	 errflg = false;
      }
   }
   return errflg;
}

static
void write_dump_header(FILE *dumpf, char *cur_file, dump_mt *dump_header, 
		      boolean xdr_write,
		      int sysinfo_size, dump_sysinfo_mt *dump_sysinfo)
{
   int vmajor, vminor;

   if( sscanf(dump_header->vsn, "%d.%d", &vmajor, &vminor) < 2 )
      message(NULLI, NULLP, WARNING, INDVSN, dump_header->vsn);

   dump_setpos(dumpf, 0L, xdr_write);
#ifdef USE_XDR
   if( xdr_write )
   {
      if( ! xdr_dump(&xdrs, dump_header) )
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
      if( dump_sysinfo ) 
      {
	 if( ! xdr_dump_sysinfo(&xdrs, dump_sysinfo, vmajor, vminor) )
	    message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
      }
   }
   else
#endif
   {
      if( fwrite((gptr*)dump_header, sizeof(dump_mt), 1, dumpf) == 0)
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
      if( dump_sysinfo ) 
      {
	 if( fwrite((gptr*)dump_sysinfo, sysinfo_size, 1, dumpf) == 0)
	    message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
      }
   }
}

/*ARGSUSED4*/
static
void write_dump_record(gptr *dump_buf, FILE *dumpf, size_mt dump_size, 
		       char *cur_file, boolean xdr_write)
{
#ifdef USE_XDR
   if( xdr_write )
   {
      if( ! xdr_vector(&xdrs, dump_buf, dump_size, sizeof(float), 
		     (xdrproc_t)xdr_float) )
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
   }
   else
#endif
   {
      if( fwrite(dump_buf, sizeof(float), dump_size, dumpf) < dump_size )
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
   }
}


/*============================================================================*/

void	dump(system_mp system, spec_mt *species, 
	     vec_mt (*force), vec_mt (*torque), 
	     mat_mt stress, double pe, 
	     restrt_mt *restart_header, int backup_restart)
{
   FILE		*dumpf;			/* File pointer to dump files	      */
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
   unsigned     dump_size = DUMP_SIZE(control.dump_level,
				      system->nmols, system->nmols_r);
                                        /* Size in floats of each dump record */
   float	*dump_buf=aalloc(dump_size,float);      /* For converted data */
   long		file_pos,		/* Offset within file of current rec. */
   		file_len;		/* Length of file		      */
   boolean      errflg = false;
   spec_mt	*spec;
   int		ispec;
   size_mt	sysinfo_size = sizeof(dump_sysinfo_mt) + 
                               sizeof(mol_mt) * (system->nspecies-1);
   dump_sysinfo_mt *dump_sysinfo =  (dump_sysinfo_mt*)balloc(1, sysinfo_size);
   		
#define		NMUTATES 10   		/* Max number of mutation attempts.   */
   int		nmutates = 0;   	/* Number of mutation attempts.	      */
   int		junk;
   boolean	xdr_write = false;	/* Is current dump in XDR format?     */
   static int	firsttime = 1;
#define REV_OFFSET 11
   char		*vsn = "$Revision: 2.27 $"+REV_OFFSET;
#define LEN_REVISION strlen(vsn)

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
	 istep_hdr -= control.maxdumps*control.dump_interval;
      errflg = true;
      /*
       * Attempt to read header and perform consistency checks which,
       * if failed, are "fatal" to dump run and initiate a *NEW* run
       * by setting errflg. 
       */
      if( (dumpf = open_dump(fname, "rb")) == 0 )
	 message(NULLI, NULLP, WARNING, DOERRR, fname, strerror(errno)); 
      else if( read_dump_header(fname, dumpf, &dump_header, 
				&xdr_write, sysinfo_size, dump_sysinfo)  )
	 message(NULLI, NULLP, WARNING, DHDERR, fname); 
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
	 file_pos = dump_curpos(sysinfo_size, dump_size, ndumps, 
				system->nspecies, xdr_write);
	 if( file_len < file_pos )
         	message(NULLI, NULLP, FATAL, CORUPT, fname, file_pos, file_len);
      }
      else
      {
	 if( ndumps != 0 )
	 {
	    ndumps = 0; filenum = 0;
	    control.begin_dump = control.istep;
            (void)sprintf(cur_file, control.dump_file, filenum);
	 }
	 message(NULLI,NULLP,WARNING,DRESET);
      }
      (void)close_dump(dumpf);
   }
   /*
    * We now know whether to start a new dump run, or are able to
    * successfully restart an existing run.
    * At this point all dump files should be closed.
    *
    * The next block is executed if we need to initialize a new dump sequence.
    */
   if( errflg || control.istep == control.begin_dump )
   {
      (void)strcpy(dump_header.title, control.title);
      (void)strncpy(dump_header.vsn, vsn, sizeof dump_header.vsn);
      dump_header.vsn[LEN_REVISION-2] = '\0';
#ifdef USE_XDR
#define LEN_XDR 5
      if( control.xdr_write )
      {
	 if(LEN_REVISION-2+LEN_XDR >= sizeof dump_header.vsn)
	    message(NULLI,NULLP,FATAL,
		    "Internal error: dump_mt header field VSN too small");
	 (void)strncat(dump_header.vsn, "(XDR)",sizeof dump_header.vsn-1);
	 dump_header.vsn[sizeof dump_header.vsn-1] = '\0';
	 xdr_write = true;
      }
#endif
      dump_header.dump_interval = control.dump_interval;
      dump_header.dump_level    = control.dump_level;
      dump_header.maxdumps	= control.maxdumps;
      dump_header.dump_size	= dump_size;
      dump_header.dump_init     = time((time_t *)0);
#ifdef USE_XDR
      if( xdr_write )
	 dump_header.sysinfo_size  =  XDR_SYSINFO_SIZE(system->nspecies);
      else
#endif      
	 dump_header.sysinfo_size  = sysinfo_size;
      dump_sysinfo->nmols        = system->nmols;
      dump_sysinfo->nmols_r      = system->nmols_r;
      dump_sysinfo->nspecies     = system->nspecies;
      dump_sysinfo->deltat       = control.step*control.dump_interval;
      for (ispec = 0, spec = species; spec < &species[system->nspecies]; 
	   ispec++, spec++)
      {
	 dump_sysinfo->mol[ispec].inertia[0] = spec->inertia[0];
	 dump_sysinfo->mol[ispec].inertia[1] = spec->inertia[1];
	 dump_sysinfo->mol[ispec].inertia[2] = spec->inertia[2];
	 dump_sysinfo->mol[ispec].mass = spec->mass;
	 dump_sysinfo->mol[ispec].dipole = spec->dipole;
	 dump_sysinfo->mol[ispec].charge = spec->charge;
	 dump_sysinfo->mol[ispec].nmols = spec->nmols;
	 dump_sysinfo->mol[ispec].rdof = spec->rdof;
	 dump_sysinfo->mol[ispec].framework = spec->framework;
	 strcpy(dump_sysinfo->mol[ispec].name, spec->name);
      }

      note(DUMPST, cur_file, control.istep);
   }
   /*
    * This block is executed if we are beginning a new dump file, either
    * continuing a sequence or starting a new one.
    *
    * No dump files should be open at this point.
    */
   if( ndumps == 0)
   {
      /*
       *  Test whether existing file belongs to current dump
       *  sequence in case of run restarted from backup.  In that
       *  case, OVERWRITE it. 
       */
      while( (dumpf = open_dump(cur_file, "r+b")) != 0 )/* File of that name */
      { 						/* already exists    */
	 if( backup_restart )
	 {   
	    errflg = read_dump_header(cur_file, dumpf, &hdr_tmp, &junk,
				      sysinfo_size, 0);
	    if( errflg == false && (hdr_tmp.dump_init == dump_header.dump_init) )
	       break;
	 }
	 (void)close_dump(dumpf);
	 
	 if( nmutates++ >= NMUTATES || mutate(control.dump_file) == NULL)
	    message(NULLI, NULLP, FATAL, MUFAIL, control.dump_file, nmutates,
                 (nmutates==1?")":"s)"));
	 message(NULLI, NULLP, WARNING, DMPEXS, cur_file, control.dump_file);
	 message(NULLI, NULLP, WARNING, DMPEXS, cur_file, control.dump_file);
	 (void)sprintf(cur_file, control.dump_file, filenum);
      }

      dump_header.ndumps = 0;	dump_header.istep = control.istep;
      dump_header.timestamp = time((time_t *)0);
      dump_header.restart_timestamp = 0;

      dumpf = open_dump(cur_file, "w+b");
      
      write_dump_header(dumpf,cur_file, &dump_header, xdr_write,
			sysinfo_size, dump_sysinfo);
   }
   else
      dumpf = open_dump(cur_file, "r+b");
   /*
    * Ready to convert and write data to dump file.
    *
    * Dump file should already be open for write.
    */
   dump_convert(dump_buf, system, species, force, torque, stress, pe);
   dump_header.ndumps++;
   dump_header.restart_timestamp = restart_header->timestamp;

   file_pos = dump_curpos(sysinfo_size, dump_size, ndumps, 
			  system->nspecies, xdr_write);
   dump_setpos(dumpf, file_pos, xdr_write);
   write_dump_record(dump_buf, dumpf, dump_size, cur_file, xdr_write);

   write_dump_header(dumpf, cur_file, &dump_header, xdr_write, sysinfo_size, 0);

   if( ferror(dumpf) || close_dump(dumpf) )
      message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));

   firsttime = 0;
   xfree(dump_buf);
   xfree(dump_sysinfo);
}
/******************************************************************************
 *  mutate  Take a string defining a file name and randomly alter characters  *
 *  to generate another, hopefully valid file name.  Assumes that a printf    *
 *  conversion code is present and only alters contiguously alphabetic chars  *
 *  before the percent.							      *
 ******************************************************************************/
#define	N_MUTATE	3			/* Max number of mutated chars*/
static char	*mutate(char *name)
{
   char	*pc_pos = strchr(name,'%'),
   	*begin  = pc_pos;
   static char	alpha[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZZ";

   while( begin > name && (pc_pos - begin) < N_MUTATE && isalnum(*(begin-1)) )
      begin--;

   if(begin >= pc_pos)				/* No characters to mutate    */
      return NULL;

   while( begin < pc_pos )
      *begin++ = alpha[(int)(mdrand1() * 26.0)];

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
static void	dump_convert(float *buf, system_mt *system, spec_mt *species,
			     vec_mt (*force), vec_mt (*torque), 
			     mat_mt stress, double pe)
{
   int		nmols   = system->nmols,
   		nmols_r = system->nmols_r;
   int		imol;
   vec_mt	*scale_buf = ralloc(nmols);
   spec_mt	*spec;
   mat_mt	hinvt, hmomrw;
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
      real_to_float(&system->ts, buf, 1);		buf += 1;
      real_to_float(&ppe, buf, 1);		buf += 1;
   }
   if( control.dump_level & 2)
   {
      /*
       * Evaluate dr/dt = 1/w*hmom*rho[i] + h^{-1}'p[i]/(m_i*s)
       * Second (main) term first.
       */
      invert(system->h, hinvt);
      transpose(hinvt, hinvt);
      mat_vec_mul(hinvt, system->mom, scale_buf, nmols);
      for (spec = species, imol=0; spec < &species[system->nspecies]; 
	   imol += spec->nmols,spec++)
	 vscale(3 * spec->nmols,   1.0/(spec->mass*system->ts), scale_buf[imol], 1);
      /*
       * 1/w*hmom*rho[i] term.
       */
      mat_sca_mul(1.0/control.pmass, system->hmom, hmomrw);
      mvaxpy(nmols, hmomrw, system->c_of_m, scale_buf);
      real_to_float(scale_buf[0],    buf, 3*nmols);	buf += 3*nmols;

      if( system->nmols_r > 0 )
      {
	 amom_to_float(species, buf, system->nspecies, 1.0/system->ts);
	 buf += 3*nmols_r;
      }
      real_to_float(system->hmom[0],    buf, 9);	buf += 9;
      real_to_float(&system->tsmom, buf, 1);		buf += 1;
   }
   if( control.dump_level & 8)
   {
      real_to_float(force[0],  buf, 3*nmols);	buf += 3*nmols;
      if( system->nmols_r > 0 )
      {
	 real_to_float(torque[0], buf, 3*nmols_r);
	 buf += 3*nmols_r;
      }
      real_to_float(stress[0], buf, 9);		/* buf += 9;*/
   }
   xfree(scale_buf);
}
/******************************************************************************
 *  real_to_float  Copy data from one array to another, converting from       *
 *  typedef real to type float in the interests of space.  (may be null op'n) *
 ******************************************************************************/
static void	real_to_float(real *b, float *a, int n)
{
   int i;
   for( i=0; i<n; i++)
      a[i] = b[i];
}
/******************************************************************************
 *  amom_to_float  Copy data from angular velocity array of quaternions to    *
 *  vector, float array.                                                      *
 ******************************************************************************/
static void	amom_to_float(spec_mt *species, float *a, 
			      int nspecies, double rts)
{
   int i,imol = 0,im;
   spec_mt *spec;
   real rinertia[3];
   
   for (spec = species; spec < &species[nspecies]; spec++)
   {
     if( spec-> rdof > 0 )
     {
	for(i=0; i<3; i++)
	{
	   if( spec->inertia[i]/(spec->inertia[(i+1)%3]+spec->inertia[(i+2)%3]) 
	       < INERTIA_MIN )
	      rinertia[i] = 0.0;
	   else
	      rinertia[i] = 1.0/spec->inertia[i];
	}

        for( im=0; im < spec->nmols; im++, imol++)
        {
	  for(i=0; i<3; i++)
	    a[3*imol+i] = spec->amom[im][i+1]*rts*rinertia[i];
	}
     }
   }
}
