#ifndef lint
static char *RCSid = "$Header: /home/tigger/keith/md/moldy/RCS/header.c,v 1.1 89/04/11 15:05:01 keith Stab $";
#endif

/*
 * $Log:	$
 */

#include	<stdio.h>
#include	"structs.h"

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
   	fprintf(stderr,"This isn't a MDKEITH restart file\n");
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