#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/dumpanalyze.c,v 1.2 89/05/15 16:51:27 keith Stab $";
#endif

/*
 * $Log:	dumpanalyze.c,v $
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
