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
static char *RCSid = "$Header: /home/kr/CVS/moldy/src/manalyze.c,v 2.10 2000/12/06 10:47:32 keith Exp $";
#endif

/*
 * $Log: manalyze.c,v $
 * Revision 2.10  2000/12/06 10:47:32  keith
 * Fixed call of make_sites() in utlsup.c to be compatible with new version.
 * Tidied up declarations and added lint flags to reduce lint noise.
 *
 * Revision 2.9  2000/04/27 17:57:09  keith
 * Converted to use full ANSI function prototypes
 *
 * Revision 2.8  1996/11/05 09:53:50  keith
 * Fixed bug which reported last record twice.
 * Now prints offsets too for ease of debugging.
 *
 * Revision 2.7  1994/06/08 13:22:31  keith
 * Null update for version compatibility
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Revision 2.5  1994/01/26  11:57:07  keith
 * Tidied up lint/gcc warnings.
 * Fixed def'n of main to "int" coz it failed on broken (?) VMS compiler.
 * Rewrote varargs functions to use stdargs conditionally on __STDC__
 *
 * Revision 2.5  1994/01/21  12:46:01  keith
 * Lint/gcc -Wall tidying
 *
 * Revision 2.4  94/01/18  13:23:14  keith
 * Incorporated all portability experience to multiple platforms since 2.2.
 * Including ports to VAX/VMS and Open VMS on Alpha AXP and Solaris.
 * 
 * Revision 2.3  93/10/28  10:28:54  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:49:46  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.5  93/03/09  15:59:54  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.4  91/08/15  18:13:13  keith
 * 
 * 
 * Revision 1.3  91/03/07  18:10:59  keith
 * Changed message
 * 
 * Revision 1.2  90/04/25  14:47:28  keith
 * Declared malloc().
 * 
 * Revision 1.1  90/02/22  17:46:11  keith
 * Initial revision
 * 
 */

#include	<time.h>
#include	"defs.h"
#include	<stddef.h>
#include	"structs.h"
#include	<stdlib.h>
#include	<stdio.h>


char	*ctime(const time_t *);
int
main(int argc, char **argv)
{
   FILE		*f = stdin;
   restrt_mt	restart_header;
   unsigned		size, offset;
   char		*ptr;
   int		n, rec = 1;
   if(argc > 1)
   {
      f = fopen(argv[1],"r");
      if(f == NULL)
      {
         fprintf(stderr,"Couldn't open restart file \"%s\"\n",argv[1]);
         exit(1);
      }
   }
   (void)fread((gptr*)&size, sizeof size, 1, f);
   if(size != sizeof restart_header)
   {
   	fprintf(stderr,"This isn't a Moldy restart file\n");
   	exit(1);
   }
   fread(&restart_header, size, 1,f);
   offset = sizeof size;
   printf("Restart file was written at %s", ctime((time_t*)&restart_header.timestamp));
   printf("This is restart No %d of run \"%s\" started %s\n",
            restart_header.seq, restart_header.title, restart_header.init_date);
   printf("It was written by version %s of write_restart\n",restart_header.vsn);
   printf("\n\tHeader record\t%d\tbytes %9X %9X\n",size, offset, size);
   offset += sizeof restart_header + sizeof size;
   do
   {
      n = fread((gptr*)&size, sizeof size, 1, f);
      if( n < 1 )
	 break;
      printf("\tRecord %d \t%d\tbytes %9X %9X\n",rec++,size, offset, size);
      offset += size + sizeof size;
      ptr = malloc(size);
      fread(ptr,size,1,f);
      free(ptr);
   } while(!feof(f));
   return 0;
}
