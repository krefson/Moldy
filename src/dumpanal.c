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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/dumpanalyze.c,v 1.3 93/03/09 15:59:49 keith Exp $";
#endif

/*
 * $Log:	dumpanalyze.c,v $
 * Revision 1.3  93/03/09  15:59:49  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.2  89/05/15  16:51:27  keith
 * Modified to work with structs.h 1.3 or later
 * 
 * Revision 1.1  89/04/11  15:04:51  keith
 * Initial revision
 * 
 */

#include <stdio.h>
#include "structs.h"
char	*ctime();
void	analyze();
void	print_header();
main(argc, argv)
int	argc;
char	*argv[];
{
   int	i;
   for(i = 1; i < argc; i++)
      analyze(argv[i]);
}

void analyze(file)
char *file;
{
   FILE *f;
   dump_mt header;

   if( f = fopen(file,"rb") )
   {
      printf("\n***** %s *****\n",file);
      if( fread((char*)&header, sizeof(dump_mt), 1, f) )
	 print_header(&header);
      else
         printf("Failed to read anything from %s",file);
      (void)fclose(f);
   }
}
void	print_header(header)
dump_mt	*header;
{
   printf("Title\t\t\t= \"%s\"\n",header->title);
   printf("RCS Revision\t\t= %.*s\n", strlen(header->vsn)-1, header->vsn);
   printf("Istep\t\t\t= %d\n",header->istep);
   printf("Dump_interval\t\t= %d\n", header->dump_interval);
   printf("Dump_level\t\t= %d\n", header->dump_level);
   printf("Max dumps\t\t= %d\n", header->maxdumps);
   printf("Dump Size\t\t= %d\n", header->dump_size);
   printf("Number of dumps\t\t= %d\n", header->ndumps);
   printf("Timestamp\t\t= %s", ctime(&header->timestamp));
   printf("Restart Timestamp\t= %s", ctime(&header->restart_timestamp));
   printf("Dump Start\t\t= %s", ctime(&header->dump_init));
}
