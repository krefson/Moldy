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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/manalyze.c,v 2.6 1994/02/17 16:38:16 keith Exp $";
#endif

/*
 * $Log: manalyze.c,v $
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

#include	"defs.h"
#include	"stddef.h"
#include	"structs.h"
#include	"stdlib.h"
#include	<stdio.h>


char	*ctime();
int
main(argc, argv)
int	argc;
char	*argv[];
{
   FILE		*f = stdin;
   restrt_mt	restart_header;
   unsigned		size;
   char		*ptr;
   int		rec = 1;
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
   printf("Restart file was written at %s", ctime(&restart_header.timestamp));
   printf("This is restart No %d of run \"%s\" started %s\n",
            restart_header.seq, restart_header.title, restart_header.init_date);
   printf("It was written by version %s of write_restart\n",restart_header.vsn);
   printf("\n\tHeader record\t%d\tbytes\n",size);
   while(!feof(f))
   {
      (void)fread((gptr*)&size, sizeof size, 1, f);
      printf("\tRecord %d \t%d\tbytes\n",rec++,size);
      ptr = malloc(size);
      fread(ptr,size,1,f);
      free(ptr);
   }
   return 0;
}
