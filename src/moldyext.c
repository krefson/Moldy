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
#include "defs.h"
#include "string.h"
#include <stdarg.h>
#include "stdlib.h"
#include "stddef.h"
#include <stdio.h>

int	getopt(int, char *const *, const char *);

#define NSIGNAL 8
#define buf_inc 128

/*VARARGS*/
void error(char *format, ...)

{
   va_list p;
   va_start(p, format);
   vfprintf(stderr,format,p);
   fputc('\n',stderr);
   va_end(p);
   exit(3);
}
static char * mystrdup(char *s)
{
   char * t=malloc(strlen(s)+1);
   return t?strcpy(t,s):0;
}
/******************************************************************************
 *  Tokenise().  Parse the string of fields to be returned and return a mask  *
 *  in a char array.  Format is 1,3,6-9,3 . . . ie comma-separated with cont- *
 *  iguous range specified with hyphen.  Numbering starts at 1.		      *
 ******************************************************************************/
int	tokenise(char *fields, char *mask, int len)
{
   char	*s;
   int	lo, hi, i, n;

   for(i = 0; i < len; i++)
      mask[i] = 0;

   while( ( s = strtok(fields,",") ) != NULL )
   {
      n = sscanf(s, "%d-%d", &lo, &hi);
      if( n == 0 )
	 return 0;

      if( n == 1 )
	 hi = lo;

      if( lo < 1 || hi < lo || hi > len)
	 return 0;

      for( i = lo-1; i < hi; i++)
	 mask[i] = 1;
      fields = NULL;
   }
   return 1;
}
/******************************************************************************
 * read_record().    Read one 'record' of Moldy output into buffer.  A record *
 * starts on the line following 8 '=' characters and enbtd either at the first*
 * line containing 8 '-' chars or a linefeed. It takes a char* pointer to a   *
 * malloc'ed buffer area, and realloc's this if it needs more space.          *
 *****************************************************************************/
char	*read_record(void)
{
   static	char	*buf = NULL;
   static	int	buf_len = 132;
   int	c, nsymb, buf_cnt;

   if( buf == 0)
   {
      if( ( buf = malloc(buf_len) ) == 0)
	 error("Memory allocation fails: %d bytes requested", buf_len);
   }
      
   nsymb = 0;
   /*
    * Read input, discarding up to and including next 8 contiguous '='
    */
   while( nsymb < NSIGNAL && ( c = getchar()) != EOF)
   {
      if( c == '=' )
	 nsymb++;
      else
	 nsymb = 0;
   }
   /*
    *  Read up to and including next newline
    */
   while( (c = getchar()) != EOF && c != '\n' )
      ; /* Empty loop body */
   /*
    *  Now read input into buffer, extending it if necessary.
    *  Read up to (and including) next 8 '-' or formfeed
    */
   buf_cnt = 0;  nsymb = 0;
   while( nsymb < NSIGNAL && (c = getchar()) != EOF && c != '\f' )
   {
      if( buf_cnt >= buf_len )		/* Buffer too small */
      {
	 buf_len += buf_inc;		/* Make it bigger   */
	 if( ( buf = realloc(buf, buf_len) ) == 0)
	    error("Memory allocation fails: %d bytes requested", buf_len);
      }
      buf[buf_cnt++] = c;
      if( c == '-' )
	 nsymb++;
      else
	 nsymb = 0;
   }
   buf_cnt -= nsymb;			/* Don't return trailing '-'  */
   buf[buf_cnt] = '\0';			/* Terminate string	      */

   return( buf );
}
/******************************************************************************
 * main program								      *
 ******************************************************************************/
#define MAX_FIELDS 256
int
main(int argc, char **argv)
{
   char *buf, *fields;
   char mask[MAX_FIELDS];
   char sfield[64];
   int field, cnt, end, inc;

   extern char *optarg;
   extern int optind;
 
   if( getopt(argc, argv, "f:") == -1)
      fields = "1-256";
   else
      fields = optarg;
      
   if( tokenise(mystrdup(fields), mask, MAX_FIELDS) == 0 )
      error("Invalid field specification \"%s\": usage eg 1,3,5-9,4",fields);

#ifdef DEBUG
   {
      int i;
      for(i = 0; i < MAX_FIELDS; i++)
         if(mask[i])
	    putchar('1');
         else
	    putchar('0');
      putchar('\n');
   }
#endif
      
   while(optind < argc)
     {
       if(strcmp(argv[optind],"-") && freopen(argv[optind], "r", stdin) == NULL)
	 error("Failed to open file \"%s\" for reading\n", argv[optind]);  
       optind++;     
       while( ! feof( stdin) )
	 {
	   buf = read_record();
	   end = strlen(buf);
	   cnt = 0;
	   field = 0;
	   while( cnt < end )
	     {
	       if(sscanf(buf+cnt, "%63s%n", sfield, &inc) != EOF)
		 {
		   if(mask[field])
		     {
		       putchar('\t');
		       fputs(sfield,stdout);
		     }
		   cnt += inc;
		   field++;
		 }
	       else
		 cnt = end;			/* Flag exit from inner loop	*/
	     }
	   if( field > 0 )
	     putchar('\n');
	 }
     }
   return 0;
}


