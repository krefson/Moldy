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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/dumpconvert.c,v 1.5 93/03/09 15:59:51 keith Exp $";
#endif

/*
 * $Log:	dumpconvert.c,v $
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

#include "stddef.h"
#include "structs.h"
#include "string.h"
#include <varargs.h>
#include <stdio.h>
char	*calloc();
void	print_header();
int	read_header();
void	binary_to_text();
void	text_to_binary();

#if defined(CRAY) && ! defined(unix)
int     vfprintf (file, fmt, args)
FILE    *file;
char    *fmt;
va_list args;
{
    return(XPRINTF(&fmt, &args, &file));
}
#endif

/*VARARGS*/
void error(va_alist)
va_dcl
{
   char	*format;
   va_list p;
   va_start(p);
   format = va_arg(p, char *);
   vfprintf(stderr,format,p);
   fputc('\n',stderr);
   va_end(p);
   exit(3);
}

main(argc, argv)
int	argc;
char	*argv[];
{
   int	i = 1, detext = 0;
   static char *ity[2] = {"rb","r"}, *oty[2] = {"w","wb"};
   

   if ( argc > 1 && ! strcmp(argv[1],"-d") )
   {
      i++;
      detext++;
   }

   if( argc > i )
   {
      if( strcmp(argv[i],"-") && freopen(argv[i], ity[detext], stdin) == NULL )
	 error("Failed to open file \"%s\" for input", argv[i]);
   }
   if( argc > ++i )
   {
      if( strcmp(argv[i],"-") && freopen(argv[i], oty[detext], stdout) == NULL )
	 error("Failed to open file \"%s\" for writing", argv[i]);
   }
	    
      
   if( detext)
      text_to_binary();
   else
      binary_to_text();
}

void binary_to_text()
{
   dump_mt header;
   float  *buf;
   int    buflen, idump, i;

   if( fread((char*)&header, sizeof(dump_mt), 1, stdin) )
      print_header(&header);
   else
      error("Failed to read header record (Error code %d).",ferror(stdin));

   buflen = header.dump_size;
   if( ! (buf = (float *)calloc(buflen, sizeof(float))))
      error("Failed to allocate memory (%d words requested)\n",buflen);

   for( idump = 0; idump < header.ndumps; idump++ )
   {
      fread(buf , sizeof(float), buflen, stdin);
      if( ferror(stdin) )
	 error("Read error on input file (Error code %d).",ferror(stdin));

      for( i = 0; i < buflen; i++ )
         printf("%.7g%c",buf[i],(i+1) % 5 ? ' ' : '\n');
      if( i % 5 ) fputc('\n', stdout);
      if( ferror(stdout) )
	 error("Write error on output file (Error code %d).",ferror(stdout));

   }
}

void	print_header(header)
dump_mt	*header;
{
   printf("%s\n%s\n",header->title, header->vsn);
   printf("%d %d %d %d %d\n", header->istep, header->dump_interval,
	  header->dump_level, header->maxdumps, header->ndumps);
   printf("%ld %ld %ld %ld\n",header->timestamp, header->restart_timestamp,
	   header->dump_init, header->dump_size); 
}


void text_to_binary()
{
   dump_mt header;
   float  *buf;
   int    buflen, idump, i;

   if( read_header(&header) )
      error("Failed to read header record (Error code %d).",ferror(stdin));

   buflen = header.dump_size;
   if( ! (buf = (float *)calloc(buflen, sizeof(float))))
      error("Failed to allocate memory (%d words requested)\n",buflen);

   if( ! fwrite((char*)&header, sizeof(dump_mt), 1, stdout) )
      error("Write error on output file (Error code %d).",ferror(stdout));

   for( idump = 0; idump < header.ndumps; idump++ )
   {
      for( i = 0; i < buflen; i++ )
        scanf("%g", &buf[i]);
      if( ferror(stdin) )
	 error("Read error on input file (Error code %d).",ferror(stdin));

      fwrite(buf , sizeof(float), buflen, stdout);
      if( ferror(stdout) )
	 error("Write error on output file (Error code %d).",ferror(stdout));

   }
}

int	read_header(header)
dump_mt	*header;
{
   int num;
   char *c;
   
   fgets(header->title, sizeof header->title, stdin);
   if(c = strchr(header->title, '\n'))
      *c = '\0';
   fgets(header->vsn, sizeof header->vsn, stdin);
   if(c = strchr(header->vsn, '\n'))
      *c = '\0';
   num  = 2;
   num += scanf("%d %d %d %d %d", &header->istep, &header->dump_interval,
	  &header->dump_level, &header->maxdumps, &header->ndumps);
   num += scanf("%ld %ld %ld %ld",&header->timestamp,
		 &header->restart_timestamp,
		 &header->dump_init, &header->dump_size);
   if( num == 11 ) return 0;
   else            return -1;
}
