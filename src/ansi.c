/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1988, 1992, 1993, 1995 Keith Refson
 
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
/******************************************************************************
 * Ansi		This file contains implementations of some ANSI C library     *
 *              functions for machines which don't have them.                 *
 ******************************************************************************/
/*========================== program include files ===========================*/
#include	"defs.h"
#ifndef _POSIX_SOURCE
#   define _POSIX_SOURCE
#endif
/*========================== Library include files ===========================*/
#include 	<stdarg.h>
#include 	<string.h>
#include	<stddef.h>
#include	<stdlib.h>
#include	<stdio.h>
#ifndef STDC_HEADERS
/******************************************************************************
 * strerror for pre-ANSI unix machines					      *
 ******************************************************************************/
#ifndef HAVE_STRERROR
char	*strerror(i)
int i;
{
   extern int sys_nerr;
   extern char *sys_errlist[];
   if( i >= 0 && i < sys_nerr )
      return sys_errlist[i];
   else
      return "invalid error code";
}   
#endif
/******************************************************************************
 * raise for pre-ANSI unix machines					      *
 ******************************************************************************/
#ifndef HAVE_RAISE
int 	raise(sig)
int sig;
{
   extern int getpid();
   extern int kill();

   return(kill(getpid(), sig));
}   
#endif
/******************************************************************************
 * strstr replacement for pre-ANSI machines which don't have it.              *
 ******************************************************************************/
#ifndef HAVE_STRSTR
char *strstr(cs, ct)
const char *cs, *ct;
{
   int  i, sl = strlen(cs)-strlen(ct);

   for(i = 0; i <= sl; i++)
      if( !strcmp(cs+i,ct) )
	 return (char*)cs+i;
   return 0;      
}
#endif
/******************************************************************************
 * mem{cpy,set} and strchr replacement for BSD machines which don't have them.*
 ******************************************************************************/
#ifndef HAVE_STRCHR
char *strchr(s, c)
const char	*s;
int	c;
{
   extern char	*index();
   return index(s,c);
}
#endif
#ifndef HAVE_MEMCPY
gptr *memcpy(s1, s2, n)
gptr *s1, *s2;
int n;
{
   int bcopy();
   (void)bcopy((char *)s2, (char *)s1, n);
   return(s1);
}
#endif
#ifndef HAVE_MEMSET
gptr *memset(s, c, n)
gptr *s;
int c, n;
{
   void bzero();
   char	*sp;
   if( c == 0 )
      (void)bzero(s, n);
   else
   {
      for( sp=s; sp < (char*)s+n; sp++)
	 *sp = c;
   }
   return(s);
}
#endif
/******************************************************************************
 * remove.  delete (unlink) a file.   (ANSI replacement)		      *
 ******************************************************************************/
#ifndef HAVE_REMOVE
int remove(file)
const char	*file;
{
   int unlink();
   return (unlink(file));
}
#endif
/******************************************************************************
 * raise for VMS.  Earlier VAX/VMS libs don't it have it.  Assume that any vsn*
 * late enough to have __DECC is OK. Don't believe ANSI_LIBS -- kludge.       *
 ******************************************************************************/
#if defined(vms) && ! defined(__DECC)
int 	raise(sig)
int sig;
{
   extern int getpid();
   extern int kill();

   return(kill(getpid(), sig));
}   
#endif
/******************************************************************************
 *  vprintf for those machines which don't have it			      *
 ******************************************************************************/
#ifndef HAVE_VPRINTF

#ifdef	HAVE_DOPRNT

int	vprintf (fmt, args)
char	*fmt;
va_list	args;
{
   return(_doprnt(fmt, args, stdout));
}


#else		/* No _doprnt or equivalent */

#include <ctype.h>
int	vprintf (format, ap)
const char	*format;
va_list	ap;
{
    int     pos, charsout, fpos, error, modflag;
    char    fmt[1024], temps[1024];
    char    ch;

    pos = charsout = error = 0;
    while (format[pos] != 0)
    {
        modflag = 0;
        if (format[pos] != '%')
        {
            putchar (format[pos++]);
        }
        else
        {
            fmt[0] = '%';
            pos++;
            fpos = 1;
            while (format[pos] != 0 && !islower (format[pos]) &&
                format[pos] != '%')
                fmt[fpos++] = format[pos++];
	    if( format[pos] == 'l' || format[pos] == 'L' || format[pos] == 'h')
	    {
	       modflag = format[pos];
	       fmt[fpos++] = format[pos++];
	    }
            if (format[pos] == 0)
            {
            /* error in format string */
                error++;
            }
            ch = fmt[fpos++] = format[pos++];
            fmt[fpos] = 0;
#ifdef DEBUG
            printf ("Format is \"%s\"\n", fmt);
#endif
        /* Now fmt contains the format for this part of the output and ch
           contains the conversion code. */
            temps[0] = 0;
            switch (ch)
            {
              case 'n':
                *va_arg (ap, int *) = charsout;
                break;
              case 's':
                sprintf (temps, fmt, va_arg (ap, char *));
                break;
              case 'p':
                sprintf (temps, fmt, va_arg (ap, void *));
                break;
              case 'c':
              case 'd':
              case 'i':
              case 'o':
              case 'u':
              case 'x':
              case 'X':
		if( modflag == 'l')
		   sprintf (temps, fmt, va_arg (ap, long));
		else
		   sprintf (temps, fmt, va_arg (ap, int));
                break;
              case 'e':
              case 'E':
              case 'f':
              case 'g':
              case 'G':
                sprintf (temps, fmt, va_arg (ap, double));
                break;
              default:
                (void)strcpy (temps, fmt);
                error++;
            }
#ifdef DEBUG
            printf ("temps is \"%s\"\n", temps);
#endif
            fputs (temps, stdout);
            charsout += strlen (temps);
#ifdef DEBUG
            printf ("still to interpret \"%s\"\n", (char *) &format[pos]);
#endif
        }
    }
    if (error)
        return (-1);
    return (charsout);
}
#endif		/* HAVE_DOPRNT       */
#endif		/* Vprintf needed    */

#endif
