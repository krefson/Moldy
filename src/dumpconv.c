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
static char *RCSid = "$Header: /home/eeyore_data/keith/moldy/src/RCS/dumpconv.c,v 2.9 1998/05/07 17:06:11 keith Exp $";
#endif

/*
 * $Log: dumpconv.c,v $
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
#include "stdlib.h"
#include "stddef.h"
#include "structs.h"
#include "string.h"
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#include <stdio.h>
#ifdef USE_XDR
#include "xdr.h"
XDR	 xdrsr;
XDR	 xdrsw;
bool_t xdr_dump();
#endif

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

#ifdef HAVE_STDARG_H
#undef  va_alist
#define	va_alist char *format, ...
#ifdef va_dcl
#   undef va_dcl
#endif
#define va_dcl /* */
#endif
/*VARARGS*/
void error(va_alist)
va_dcl
{
   va_list p;
#ifdef HAVE_STDARG_H
   va_start(p, format);
#else
   char	*format;

   va_start(p);
   format = va_arg(p, char *);
#endif
   vfprintf(stderr,format,p);
   fputc('\n',stderr);
   va_end(p);
   exit(3);
}

void read_text(buf,buflen)
float  *buf;
int    buflen;
{
   int i;

   for( i = 0; i < buflen; i++ )
      scanf("%f", &buf[i]);
   if( ferror(stdin) )
      error("Read error on input file (Error code %d).",ferror(stdin));
}

void write_text(buf,buflen)
float  *buf;
int    buflen;
{
   int i;

   for( i = 0; i < buflen; i++ )
      printf("%.7g%c",buf[i],(i+1) % 5 ? ' ' : '\n');
   if( i % 5 ) fputc('\n', stdout);
   if( ferror(stdout) )
      error("Write error on output file (Error code %d).",ferror(stdout));
}

void read_binary(buf,buflen,xdr)
float  *buf;
int    buflen, xdr;
{
#ifdef USE_XDR
   if( xdr )
   {
      xdr_vector(&xdrsr, (char*)buf, buflen, XDR_FLOAT_SIZE, xdr_float);
   }
   else
#endif
   {
      fread(buf , sizeof(float), buflen, stdin);
   }
   if( ferror(stdin) )
      error("Read error on input file (Error code %d).",ferror(stdin));
}

void write_native(buf, buflen)
float  *buf;
int    buflen;
{
   fwrite(buf , sizeof(float), buflen, stdout);
   if( ferror(stdout) )
      error("Write error on output file (Error code %d).",ferror(stdout));
}

#ifdef USE_XDR
void write_xdr(buf, buflen)
float  *buf;
int    buflen;
{
   xdr_vector(&xdrsw, (char*)buf, buflen, XDR_FLOAT_SIZE, xdr_float);
   if( ferror(stdout) )
      error("Write error on output file (Error code %d).",ferror(stdout));
}
#endif

void read_bin_hdr(header, xdrw)
dump_mt *header;
int     *xdrw;
{
   int    xdr = 0, errflg = 0;

#ifdef USE_XDR
   /*
    * Attempt to read dump header in XDR format
    */
   xdrstdio_create(&xdrsr, stdin, XDR_DECODE);
   if( xdr_dump(&xdrsr, header) )
   {
      header->vsn[sizeof header->vsn - 1] = '\0';
      if( strstr(header->vsn,"(XDR)") )
      {
	 errflg = 0;
	 xdr = 1;
      }
   }
   else
      errflg = 1;
#endif
   /*
    * If we failed, try to read header as native struct image.
    */
   if( ! xdr )
   {
      if( fseek(stdin, 0L, 0) == 0 &&
	 fread((char*)header, sizeof(dump_mt), 1, stdin) == 1) 
	 errflg = 0;
   }
   if( errflg || ferror(stdin) || feof(stdin) )
   {
      error("Failed to read header record (Error code %d).",ferror(stdin));
      exit(2);
   }
   *xdrw = xdr;
}

int	read_header(header)
dump_mt	*header;
{
   int num;
   char *c;
   
   fgets(header->title, sizeof header->title, stdin);
   if((c = strchr(header->title, '\n')))
      *c = '\0';
   fgets(header->vsn, sizeof header->vsn, stdin);
   if((c = strchr(header->vsn, '\n')))
      *c = '\0';
   num  = 2;
   num += scanf("%d %d %d %d %d", &header->istep, &header->dump_interval,
	  &header->dump_level, &header->maxdumps, &header->ndumps);
   num += scanf("%ld %ld %ld %d",&header->timestamp,
		 &header->restart_timestamp,
		 &header->dump_init, &header->dump_size);
   if( num == 11 ) return 0;
   else            return -1;
}

void write_native_hdr(header)
dump_mt *header;
{
   char *s;
   if( (s = strstr(header->vsn,"(XDR)") ) != 0 )
      *s = 0;
   if( ! fwrite((char*)header, sizeof(dump_mt), 1, stdout) )
      error("Write error on output file (Error code %d).",ferror(stdout));
}

#ifdef USE_XDR
void write_xdr_hdr(header)
dump_mt *header;
{
   strncat(header->vsn,"(XDR)",sizeof header->vsn);
   xdrstdio_create(&xdrsw, stdout, XDR_ENCODE);
   if( ! xdr_dump(&xdrsw, header) )
      error("Write error on output file (Error code %d).",ferror(stdout));
}
#endif

void	print_header(header)
dump_mt	*header;
{
   char *s;
   if( (s = strstr(header->vsn,"(XDR)") ) != 0 )
      *s = 0;
   printf("%s\n%s\n",header->title, header->vsn);
   printf("%d %d %d %d %d\n", header->istep, header->dump_interval,
	  header->dump_level, header->maxdumps, header->ndumps);
   printf("%ld %ld %ld %d\n",header->timestamp, header->restart_timestamp,
	   header->dump_init, header->dump_size); 
}

int
main(argc, argv)
int	argc;
char	*argv[];
{
   static char *ity[2] = {"rb","r"}, *oty[2] = {"w","wb"};
   dump_mt header;
   float *buf;
   int   idump;
   int   textin = 0, xdrin, xdrout = 0;

   while( argc > 0 && argv[1][0] == '-' )
   {
      switch( argv[1][1] ) {
       case 'd':	
	 textin++;
	 break;
       case 'x':
	  xdrout++;
	 break;
       default:
	 fprintf(stderr,"Usage: dumpconvert [-d] [-x] infile outfile\n");
	 exit(2);
      }
      argc--; argv++;
   }

   if( argc > 0 )
   {
      if( strcmp(argv[1],"-") && freopen(argv[1], ity[textin], stdin) == NULL )
	 error("Failed to open file \"%s\" for input", argv[1]);
   }
   if( argc > 1 )
   {
      if( strcmp(argv[2],"-") && freopen(argv[2], oty[textin], stdout) == NULL )
	 error("Failed to open file \"%s\" for writing", argv[2]);
   }
	    
      
   if( textin )
   {
      read_header(&header);
#ifdef USE_XDR
      if( xdrout )
	 write_xdr_hdr(&header);
      else
#endif
	 write_native_hdr(&header);
   }
   else
   {
      read_bin_hdr(&header, &xdrin);
#ifdef USE_XDR
      if( xdrout )
	 write_xdr_hdr(&header);
      else
#endif
	 print_header(&header);
   }

   if( ! (buf = (float *)calloc(header.dump_size, sizeof(float))))
      error("Failed to allocate memory (%d words requested)\n",
	    header.dump_size);
   
   for( idump = 0; idump < header.ndumps; idump++)
   {
      if( textin )
      {
	 read_text(buf, header.dump_size);
#ifdef USE_XDR
	 if( xdrout )
	    write_xdr(buf, header.dump_size);
	 else
#endif
	    write_native(buf, header.dump_size);
      }
      else
      {
	 read_binary(buf, header.dump_size, xdrin);
#ifdef USE_XDR
	 if( xdrout )
	    write_xdr(buf, header.dump_size);
	 else
#endif
	    write_text(buf, header.dump_size);
      }
   }
return 0;
}
