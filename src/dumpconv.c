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
#ifndef lint
static char *RCSid = "$Header: /home/kr/CVS/moldy/src/dumpconv.c,v 2.16 2001/08/01 11:56:36 keith Exp $";
#endif

/*
 * $Log: dumpconv.c,v $
 * Revision 2.16  2001/08/01 11:56:36  keith
 * Incorporated all info from "species" struct into dump file headers.
 * - fixed utilities and a few bugs.
 *
 * Revision 2.15  2001/07/31 17:58:18  keith
 * Incorporated all info from "species" struct into dump file headers.
 *
 * Revision 2.14  2000/11/15 17:51:59  keith
 * Changed format of dump files.
 * Added second struct with sufficient information
 * about the simulation that most utility programs
 * (namely those which do not need site co-ordinates)
 * should not need to read sys-spec or restart files.
 *
 * New options "-c -1" to dumpext prints header info.
 * -- dumpanal removed.
 *
 * Revision 2.13  2000/11/10 12:16:27  keith
 * Tidied up some dubious cases to get rid of compiler warnings.
 * Updated configure scripts -- fix for non-pgcc linux case.
 * Got rid of redundant Makefile.w32 and Makefile.mak
 * make -f xmakefile Makefile.in now works under Linux
 *
 * Revision 2.12  2000/11/09 16:28:03  keith
 * Updated dump file format for new Leapfrog dynamics\nAdded #molecules etc to dump header format
 *
 * Revision 2.11  2000/04/27 17:57:06  keith
 * Converted to use full ANSI function prototypes
 *
 * Revision 2.10  1998/05/22 11:11:47  keith
 * Protected va_dcl redefinition.
 *
 * Revision 2.9  1998/05/07 17:06:11  keith
 * Reworked all conditional compliation macros to be
 * feature-specific rather than OS specific.
 * This is for use with GNU autoconf.
 *
 * Revision 2.8  1995/12/05 20:55:10  keith
 * Separated ANSI replacement routines from Auxil.c into Ansi.c
 * Removed all COS functionality.
 *
 * Revision 2.7  1994/06/08 13:12:34  keith
 * Changed "%g" scan format to "%f" - should be identical but
 * some systems, particularly VAX/VMS didn't grok "%g".
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Revision 2.5  94/01/21  12:51:04  keith
 * Corrected declaration of main()
 * 
 * Revision 2.4  94/01/21  12:35:03  keith
 * Incorporated all portability experience to multiple platforms since 2.2.
 * Rewrote varargs functions to use stdargs conditionally on __STDC__
 * 
 * Revision 2.5  1994/01/18  16:26:29  keith
 * Incorporated all portability experience to multiple platforms since 2.2.
 * Rewrote varargs functions to use stdargs conditionally on __STDC__
 *
 * Revision 2.5  94/01/18  13:13:09  keith
 * Incorporated all portability experience to multiple platforms since 2.2.
 * Rewrote varargs functions to use stdargs conditionally on __STDC__
 * 
 * Revision 2.4  93/12/21  18:49:40  keith
 * Portability improvements:
 * 1. Moved malloc etc declarations into header files
 * 2. Rewrote varargs functions to use stdargs conditionally on __STDC__
 * 
 * Revision 2.3  93/10/28  10:28:50  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.2  93/09/06  14:42:43  keith
 * Fixed portability problems/bugs in XDR code.
 * 
 * Revision 2.1  93/08/18  20:52:05  keith
 * Added support for dumps in XDR format.
 * 
 * Revision 2.0  93/03/15  14:49:44  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.5  93/03/09  15:59:51  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.4  91/08/15  18:13:10  keith
 * 
 * 
 * Revision 1.3  90/09/28  10:51:27  keith
 * Amended #idfefs for unicos
 * 
 * Revision 1.2  89/09/07  18:15:54  keith
 * checked in with -k by keith at 89.09.08.15.48.37.
 * 
 * Revision 1.2  89/09/07  18:15:54  keith
 * Rationalised command-line parameters.  Now takes option '-d', + two
 * file names for input and output.  A file name of '-' means use stdin/out
 * as does an absent name.   
 * This should work better on non-unix machines.
 * 
 * Revision 1.1  89/06/27  12:10:34  keith
 * Initial revision
 * 
 */

#include "defs.h"
#include "messages.h"
#include "stdlib.h"
#include "stddef.h"
#include "structs.h"
#include "string.h"
#include <stdarg.h>
#include <stdio.h>
#ifdef USE_XDR
#include "xdr.h"
XDR	 xdrs;
XDR	 xdrsw;
bool_t xdr_dump(XDR *xdrs, dump_mt *sp);
#endif

#undef OLDVSN
#define OLDVSN "Can not determine dump file version - \"%s\""

/******************************************************************************
 * strstr replacement for pre-ANSI machines which don't have it.              *
 ******************************************************************************/
#ifndef STDC_HEADERS
char *strstr(cs, ct)
char *cs, *ct;
{
   char *end = cs+strlen(cs)-strlen(ct);
   for(; cs <= end; cs++)
      if( !strcmp(cs,ct) )
	 return cs;
   return 0;      
}
#endif

/*VARARGS*/
void error(char *format, ...)

{
   va_list p;
   va_start(p, format);
   vfprintf(stderr,format,p);
   fputc('\n',stderr);
   va_end(p);
   exit(3);
}

/******************************************************************************
 *  message.   Deliver error message to possibly exiting.  It can be called   *
 *             BEFORE output file is opened, in which case output to stderr.  *
 ******************************************************************************/
char    *comm = "dumpconv";             
/*VARARGS*/
void    message(int *nerrs, ...)
{
   va_list      ap;
   char         *buff;
   int          sev;
   char         *format;
   static char  *sev_txt[] = {" *I* "," *W* "," *E* "," *F* "};
   va_start(ap, nerrs);
   buff  = va_arg(ap, char *);
   sev   = va_arg(ap, int);
   format= va_arg(ap, char *);
 
   (void)fprintf(stderr,"%s: ",comm);
   (void)vfprintf(stderr, format, ap);
   va_end(ap);
   fputc('\n',stderr);
 
   if(buff != NULL)                     /* null ptr means don't print buffer  */
   {
      (void)fprintf(stderr,"     buffer contents=\"%s\"",buff);
      fputc('\n',stderr);
   }
   if(sev >= ERROR && nerrs != NULL)
      (*nerrs)++;
   if(sev == FATAL)
      exit(3);
}

void read_text(float *buf, int buflen)
{
   int i;

   for( i = 0; i < buflen; i++ )
      scanf("%f", &buf[i]);
   if( ferror(stdin) )
      error("Read error on input file (Error code %d).",ferror(stdin));
}

void write_text(float *buf, int buflen)
{
   int i;

   for( i = 0; i < buflen; i++ )
      printf("%.8g%c",buf[i],(i+1) % 5 ? ' ' : '\n');
   if( i % 5 ) fputc('\n', stdout);
   if( ferror(stdout) )
      error("Write error on output file (Error code %d).",ferror(stdout));
}

static dump_sysinfo_mt *dump_sysinfo;

int	read_header(dump_mt *header, dump_sysinfo_mt **sysinfo)
{
   int num, ispec, nmols, nmols_r, nspecies;
   float deltat;
   char *c, next;
   size_mt sysinfo_size;
   mol_mt *mol_p;
   
   fgets(header->title, sizeof header->title, stdin);
   if((c = strchr(header->title, '\n')))
      *c = '\0';
   fgets(header->vsn, sizeof header->vsn, stdin);
   if((c = strchr(header->vsn, '\n')))
      *c = '\0';
   num  = 2;
   num += scanf("%ld %ld %d %d %d", &header->istep, &header->dump_interval,
	  &header->dump_level, &header->maxdumps, &header->ndumps);
   num += scanf("%ld %ld %ld %d %ld",&header->timestamp, &header->restart_timestamp,
		 &header->dump_init, &header->dump_size, &header->sysinfo_size);
   if( num < 12 ) 
      return -1;
   num = scanf("%f %d %d %d ", &deltat, &nmols, &nmols_r, &nspecies);
   sysinfo_size = sizeof(dump_sysinfo_mt) + sizeof(mol_mt) * (nspecies-1);
   if ( ( dump_sysinfo = *sysinfo = malloc(header->sysinfo_size)) == 0)
      error("Failed to malloc memory for sysinfo header");

   dump_sysinfo->deltat =   deltat;
   dump_sysinfo->nmols  =   nmols;
   dump_sysinfo->nmols_r=   nmols_r;
   dump_sysinfo->nspecies = nspecies;
   for(ispec = 0; ispec < dump_sysinfo->nspecies; ispec++)
   {
      mol_p = &dump_sysinfo->mol[ispec];
      fgets(dump_sysinfo->mol[ispec].name, L_spec, stdin);     
      if((c = strchr(dump_sysinfo->mol[ispec].name, '\n')))
	 *c = '\0';
      num ++;
      num += scanf("%d %d %d", &mol_p->nmols, &mol_p->framework, &mol_p->rdof);
      num += scanf("%f %f %f %f", &mol_p->mass, &mol_p->inertia[0],
		   &mol_p->inertia[1],&mol_p->inertia[2]);
      num += scanf("%f %f", &mol_p->charge, &mol_p->dipole);
      do 
      {
	 next = fgetc(stdin);
      } while( isspace(next) );
      ungetc(next, stdin);
   }
   if( num < 4+10*dump_sysinfo->nspecies) 
      return -1;
   return 0;
}

void	print_header(dump_mt *header, dump_sysinfo_mt *sysinfo)
{
   char *s;
   int ispec;
   if( (s = strstr(header->vsn,"(XDR") ) != 0 )
      *s = 0;
   printf("%s\n%s\n",header->title, header->vsn);
   printf("%ld %ld %d %d %d\n", header->istep, header->dump_interval,
	  header->dump_level, header->maxdumps, header->ndumps);
   printf("%ld %ld %ld %d %ld\n",header->timestamp, header->restart_timestamp,
	   header->dump_init, header->dump_size, header->sysinfo_size);
   printf("%.8g %d %d %d\n", sysinfo->deltat, sysinfo->nmols, 
	                     sysinfo->nmols_r, sysinfo->nspecies);
   for(ispec = 0; ispec < sysinfo->nspecies; ispec++)
   {
      printf("%s\n%d %d %d\n",sysinfo->mol[ispec].name,sysinfo->mol[ispec].nmols,
	                   sysinfo->mol[ispec].framework, sysinfo->mol[ispec].rdof);
      printf("%.8g %.8g %.8g %.8g\n", sysinfo->mol[ispec].mass,
	     sysinfo->mol[ispec].inertia[0],sysinfo->mol[ispec].inertia[1],
	     sysinfo->mol[ispec].inertia[2]);
      printf("%.8g %.8g\n",sysinfo->mol[ispec].charge,sysinfo->mol[ispec].dipole);
   }
}


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

FILE  *reopen_dump(FILE *dumpf, char *mode)
{

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

int rewind_dump(FILE *dumpf, int xdr)
{
#ifdef USE_XDR
   if( xdr )
      xdr_setpos(&xdrs, 0);
#endif
   return fseek(dumpf, 0L, SEEK_SET);
}

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
		     size_mt sysinfo_size, dump_sysinfo_mt *dump_sysinfo)
{
   int      errflg = true;	/* Provisionally !!   */
   char     vbuf[sizeof hdr_p->vsn + 1];
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
      if( vmajor < 2 || vminor <= 17)
	 message(NULLI, NULLP, WARNING, OLDVSN, hdr_p->vsn);
      else
	 errflg = false;
   }
   if( errflg ) return errflg;

   
   if( dump_sysinfo == 0)
      return errflg;
   else if ( sysinfo_size == sizeof(dump_sysinfo_mt) )
   {
      /*
       * Now check for sysinfo and read fixed part of it.  This is needed to
       * determine species count to read the whole thing.
       */
#ifdef USE_XDR
      if( *xdr_write ) {
	 if( ! xdr_dump_sysinfo_hdr(&xdrs, dump_sysinfo) )
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
	 errflg = false;
      } else
#endif
      {
	 if( fread((gptr*)dump_sysinfo,sizeof(dump_sysinfo_mt), 1, dumpf) == 0)
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
	 errflg = false;
      }
   }
   else
   {
      /*
       * Now check for sysinfo and read it all.  N.B.  Buffer must be 
       * allocated to full expected size by prior call to read_dump_header.
       */
#ifdef USE_XDR
      if( *xdr_write ) {
	 if( ! xdr_dump_sysinfo(&xdrs, dump_sysinfo, vmajor, vminor) )
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
      if (sizeof(dump_sysinfo_mt) 
	  + sizeof(mol_mt) * (dump_sysinfo->nspecies-1) > sysinfo_size)
      {
	 /*
	  * We have already overrun the end of the "dump_sysinfo" buffer.
	  * Perhaps we can exit gracefully before crashing?
          */
	 message(NULLI, NULLP, FATAL, RDHERR,  sizeof(dump_sysinfo_mt) 
		 + sizeof(mol_mt) * (dump_sysinfo->nspecies-1), sysinfo_size);
      }
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

static
void read_dump_record(gptr *dump_buf, FILE *dumpf, size_mt dump_size, 
		       char *cur_file, boolean xdr)
{
#ifdef USE_XDR
   if( xdr )
   {
      if( ! xdr_vector(&xdrs, dump_buf, dump_size, sizeof(float), 
		     (xdrproc_t)xdr_float) )
	 message(NULLI, NULLP, FATAL, DWERR, cur_file, strerror(errno));
   }
   else
#endif
   {
      if( fread(dump_buf, sizeof(float), dump_size, dumpf)== 0 )
	 message(NULLI, NULLP, FATAL, DRERR, cur_file, strerror(errno));
   }
}

int
main(int argc, char **argv)
{
   static char *ity[2] = {"rb","r"}, *oty[2] = {"w","wb"};
   dump_mt header;
   float *buf;
   int   idump;
   int   textin = 0, textout = 1, xdrin, xdrout = 1;
   dump_sysinfo_mt *dump_sysinfo;
   size_mt	sysinfo_size;

   while( argc > 0 && argv[1][0] == '-' )
   {
      switch( argv[1][1] ) {
       case 'b':	
	 textin  = 0;
	 textout = 0;
	 break;
       case 'd':	
	 textin++;
	 textout = 0;
	 break;
       case 'n':
	  xdrout=0;
	 break;
       default:
	 fprintf(stderr,"Usage: dumpconvert [-b] [-d] [-n] infile outfile\n");
	 exit(2);
      }
      argc--; argv++;
   }

   if( argc > 1 )
   {
      if( strcmp(argv[1],"-") && freopen(argv[1], ity[textin], stdin) == NULL )
	 error("Failed to open file \"%s\" for input", argv[1]);
   }
   if( argc > 2 )
   {
      if( strcmp(argv[2],"-") && freopen(argv[2], oty[textin], stdout) == NULL )
	 error("Failed to open file \"%s\" for writing", argv[2]);
   }
	    
      
   if( textin )
   {
      if( read_header(&header,&dump_sysinfo) )
	 error("Failed to read text header");
   }
   else
   {
      reopen_dump(stdin,"rb");
      sysinfo_size = sizeof(dump_sysinfo_mt);
      dump_sysinfo = (dump_sysinfo_mt*)malloc(sysinfo_size);
      if( read_dump_header("", stdin, &header, &xdrin, sysinfo_size, dump_sysinfo) )
	 error("Failed to read dump header");
      sysinfo_size = sizeof(dump_sysinfo_mt) + sizeof(mol_mt) * (dump_sysinfo->nspecies-1);
      (void)free(dump_sysinfo);
       dump_sysinfo = (dump_sysinfo_mt*)malloc(sysinfo_size);
      /*
       * Rewind and reread header, this time including all of sysinfo.
       */
      (void)rewind_dump(stdin, xdrin);
      if( read_dump_header("", stdin, &header, &xdrin, 
			   header.sysinfo_size, dump_sysinfo) 
	  || ferror(stdin) || feof(stdin) )
	 error("Failed to read dump header.");
   }
      

   if( textout )
      print_header(&header, dump_sysinfo);
   else
   {
      reopen_dump(stdout,"wb");
      if( xdrout )
	   strncat(header.vsn,"(XDR)",sizeof header.vsn);
      write_dump_header(stdout,"", &header, xdrout, header.sysinfo_size, dump_sysinfo);
   }
      

   if( ! (buf = (float *)calloc(header.dump_size, sizeof(float))))
      error("Failed to allocate memory (%d words requested)\n",
	    header.dump_size);
   
   for( idump = 0; idump < header.ndumps; idump++)
   {
      if( textin )
	 read_text(buf, header.dump_size);
      else
	 read_dump_record(buf, stdin, header.dump_size, "", xdrin);

      if( textout )
	 write_text(buf, header.dump_size);
      else
	 write_dump_record(buf, stdout, header.dump_size, "", xdrout);
   }
return 0;
}
