#include <stdio.h>
#include "string.h"
#include <varargs.h>
#include "stddef.h"
char	*malloc(), *realloc();
char	*strtok();

#define NSIGNAL 8
#define buf_inc 128

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
/******************************************************************************
 *  Tokenise().  Parse the string of fields to be returned and return a mask  *
 *  in a char array.  Format is 1,3,6-9,3 . . . ie comma-separated with cont- *
 *  iguous range specified with hyphen.  Numbering starts at 1.		      *
 ******************************************************************************/
int	tokenise(fields, mask, len)
char	*fields, *mask;
int	len;
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
char	*read_record()
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
main(argc, argv)
int	argc;
char	*argv[];
{
   char *buf, *fields;
   char mask[MAX_FIELDS];
   char sfield[64];
   int field, cnt, end, inc;

   extern char *optarg;
   extern int optind, opterr;
 
   if( getopt(argc, argv, "f:") == -1)
      fields = "1-256";
   else
      fields = optarg;
      
   if( argc > optind )
   {
      if(strcmp(argv[optind],"-") && freopen(argv[optind], "r", stdin) == NULL)
	 error("Failed to open file \"%s\" for input", argv[optind]);
   }
   if( argc > ++optind )
   {
      if(strcmp(argv[optind],"-") && freopen(argv[optind], "w", stdout) == NULL)
	 error("Failed to open file \"%s\" for writing", argv[optind]);
   }

   if( tokenise(strdup(fields), mask, MAX_FIELDS) == 0 )
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
      putchar('\n');
   }
}


