#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/moldyanalyze.c,v 1.4 91/08/13 15:31:40 keith Exp $";
#endif

/*
 * $Log:	moldyanalyze.c,v $
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

#include	<stdio.h>
#include	"stddef.h"
#include	"structs.h"

char	*malloc();

#undef cfree

char	*ctime();
main(argc, argv)
int	argc;
char	*argv[];
{
   FILE		*f = stdin;
   restrt_t	restart_header;
   int		size;
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
   size = getw(f);
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
      size = getw(f);
      printf("\tRecord %d \t%d\tbytes\n",rec++,size);
      ptr = malloc(size);
      fread(ptr,size,1,f);
      cfree(ptr);
   }
}
