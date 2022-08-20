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

/*
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
