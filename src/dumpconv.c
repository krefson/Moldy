#ifndef lint
static char *RCSid = "$Header: dumpanalyze.c,v 1.1 89/04/11 15:04:51 keith Exp $";
#endif

/*
 * $Log:	dumpanalyze.c,v $
 */

#include <stdio.h>
#include "structs.h"
#include "string.h"
#include <varargs.h>
char	*calloc();
void	print_header();
int	read_header();
void	binary_to_text();
void	text_to_binary();

#define EXTN ".T"

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

   if ( argc > 1 && ! strcmp(argv[1],"-d") )
   {
      i=2;
      detext++;
   }
   
   for(; i < argc; i++)
      if( detext)
         text_to_binary(strcat(strdup(argv[i]),EXTN), argv[i]);
      else
         binary_to_text(argv[i], strcat(strdup(argv[i]),EXTN));
}

void binary_to_text(in_file, out_file)
char *in_file, *out_file;
{
   FILE *in, *out;
   dump_t header;
   float  *buf;
   int    buflen, idump, i;

   if( ! (in = fopen(in_file,"rb") ) )
      error("Failed to open file \"%s\" for input", in_file);

   if( ! ( out = fopen(out_file,"w") ) )
      error("Failed to open file \"%s\" for writing", out_file);

   if( fread((char*)&header, sizeof(dump_t), 1, in) )
      print_header(&header, out);
   else
      error("Failed to read header record from %s",in_file);

   buflen = header.dump_size;
   if( ! (buf = (float *)calloc(buflen, sizeof(float))))
      error("Failed to allocate memory (%d words requested)\n",buflen);

   for( idump = 0; idump < header.ndumps; idump++ )
   {
      fread(buf , sizeof(float), buflen, in);
      if( ferror(in) )
	 error("Read error on file %s",in_file);

      for( i = 0; i < buflen; i++ )
         fprintf(out,"%.7g%c",buf[i],(i+1) % 5 ? ' ' : '\n');
      if( i % 5 ) fputc('\n', out);
      if( ferror(out) )
	 error("Write error on file %s",out_file);

   }
   fclose(in);
   fclose(out);
}

void	print_header(header, out)
dump_t	*header;
FILE	*out;
{
   fprintf(out,"%s\n%s\n",header->title, header->vsn);
   fprintf(out,"%d %d %d %d %d\n", header->istep, header->dump_interval,
	  header->dump_level, header->maxdumps, header->ndumps);
   fprintf(out,"%ld %ld %ld %ld\n",header->timestamp, header->restart_timestamp,
	   header->dump_init, header->dump_size); 
}


void text_to_binary(in_file, out_file)
char *in_file, *out_file;
{
   FILE *in, *out;
   dump_t header;
   float  *buf;
   int    buflen, idump, i;

   if( ! (in = fopen(in_file,"r") ) )
      error("Failed to open file \"%s\" for input", in_file);

   if( read_header(&header, in) )
      error("Failed to read header record from %s",in_file);

   buflen = header.dump_size;
   if( ! (buf = (float *)calloc(buflen, sizeof(float))))
      error("Failed to allocate memory (%d words requested)\n",buflen);

   if( ! ( out = fopen(out_file,"wb") ) )
      error("Failed to open file \"%s\" for writing", out_file);

   if( ! fwrite((char*)&header, sizeof(dump_t), 1, out) )
      error("Failed to write to file %s\n", out_file);

   for( idump = 0; idump < header.ndumps; idump++ )
   {
      for( i = 0; i < buflen; i++ )
        fscanf(in,"%g", &buf[i]);
      if( ferror(in) )
	 error("Read error on file %s",in_file);

      fwrite(buf , sizeof(float), buflen, out);
      if( ferror(out) )
	 error("Write error on file %s",out_file);

   }
   fclose(in);
   fclose(out);
}

int	read_header(header, in)
dump_t	*header;
FILE	*in;
{
   int num;
   fgets(header->title, sizeof header->title, in);
   *strchr(header->title, '\n') = '\0';
   fgets(header->vsn, sizeof header->vsn, in);
   *strchr(header->vsn, '\n') = '\0';
   num  = 2;
   num += fscanf(in,"%d %d %d %d %d", &header->istep, &header->dump_interval,
	  &header->dump_level, &header->maxdumps, &header->ndumps);
   num += fscanf(in,"%ld %ld %ld %ld",&header->timestamp,
		 &header->restart_timestamp,
		 &header->dump_init, &header->dump_size);
   if( num == 11 ) return 0;
   else            return -1;
}
