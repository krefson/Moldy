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
/******************************************************************************
 * Aux		This file contains auxilliary service functions and (almost)  *
 *		all machine-dependant stuff.  The ANSI standard interface is  *
 *              used as much as possible for library functions; if a machine  *
 *		is missing something an ANSI interface routine is supplied.   *
 ******************************************************************************
 *      Revision Log
 *       $Log: auxil.c,v $
 *       Revision 2.8  1995/12/05 20:55:10  keith
 *       Separated ANSI replacement routines from Auxil.c into Ansi.c
 *       Removed all COS functionality.
 *
 *       Revision 2.7  1994/06/08 13:09:00  keith
 *       Added "zero_dbls()" function.
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Got rid of all global (external) data items except for
 * "control" struct and constant data objects.  The latter
 * (pot_dim, potspec, prog_unit) are declared with CONST
 * qualifier macro which evaluates to "const" or nil
 * depending on ANSI/K&R environment.
 * Also moved as many "write" instantiations of "control"
 * members as possible to "startup", "main" leaving just
 * "dump".
 *
 * Moved declaration of "match" structure from input.c
 * Moved a few sanity tests & modifications of "control"
 * members to here from other modules.
 *
 * Changed size_t to own typedef size_mt == ulong.
 *
 * Declared as "static"  all functions which should be.
 *
 * Added CONST qualifier to (re-)declarations of ANSI library
 * emulation routines to give reliable compilation even
 * without ANSI_LIBS macro. (#define's away for K&R
 * compilers)
 *
 * Renamed Aux.c to Auxil.c for DOS (ugh).
 * Revamped cpu time interface to use POSIX "times()" whenever
 * available - including for BSD systems - and got rid of
 * "getrusage()" version entirely.
 * Thus eliminating various marginally portable header
 * kludges.
 *
 * Added memset() emulation for *really* old BSD systems.
 *
 * Revision 2.5  94/01/29  11:09:53  keith
 * Fixed bug in "strstr()" replacement for non-ANSI libraries
 * 
 * Revision 2.4  94/01/18  13:23:04  keith
 * Incorporated all portability experience to multiple platforms since 2.2.
 * Including ports to VAX/VMS and Open VMS on Alpha AXP and Solaris.
 * 
 * Revision 2.3  93/10/28  10:27:38  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.2  93/09/06  14:43:28  keith
 * Fixed portability problems/bugs in XDR code.
 * 
 * Revision 2.0  93/03/15  14:48:56  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.7.1.32  93/03/15  14:41:33  keith
 * Added GPL copyleft notice to permit release and distribution.
 * N.B.  Previous versions were copyright (C) by author and 
 * only licensed by explicit permission.
 * 
 * Revision 1.7.1.31  93/03/12  12:19:22  keith
 * Reorganized defines to recognise all ANSI (__type__) forms.
 * Moved spxpy() from aux.c to force.c and force_parallel.c
 * 
 * Revision 1.7.1.30  93/03/09  15:58:18  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.7.1.29  93/03/05  15:02:58  keith
 * Now is POSIX - safe wrt declaration of times().
 * 
 * Revision 1.7.1.28  92/09/22  14:48:22  keith
 * Tidied up calls to improve "lint" rating.
 * 
 * Revision 1.7.1.27  92/06/26  17:02:46  keith
 * Got rid of assumption that memory returned by talloc() or
 * arralloc() is zeroed.  This enhances ANSI compatibility.
 * Removed memory zeroing from alloc.c() in consequence.
 * 
 * Revision 1.7.1.26  92/06/19  17:28:40  keith
 * Added faster version of "search_lt" for Titan.
 * 
 * Revision 1.7.1.25  92/06/12  12:55:41  keith
 * Mods to make it work on VMS again.  Ugh.
 * 
 * Revision 1.7.1.24  92/06/11  20:31:04  keith
 * Added file locking against multiple runs using same dump or backup files.
 * 
 * Revision 1.7.1.23  91/11/26  10:24:31  keith
 * Modified timing routines to work under POSIX.
 * Replaced saxpy() with sxpy() -- UNVECTORIZED as indices may be repeated.
 * 
 * Revision 1.7.1.22  91/08/17  16:24:15  keith
 * Added strerror() for pre-ANSI unix systems
 * Updated #defines with __unix__ symbol and protection of sys/time.h for convex.
 * 
 * Revision 1.7.1.21  91/08/15  18:11:44  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.7.1.20  91/05/29  10:58:00  keith
 * Modified wheneq to generate better code on Stardent TITAN
 * Added #ifdefe'd "const" qualifier in "remove" for ANSI compilers
 * without ANSI libraries.
 * 
 * Revision 1.7.1.19  91/03/12  15:42:12  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.7.1.18  91/02/07  16:48:36  keith
 * Modified BLAS functions to vectorise better undar Stardent Titan compiler.
 * 
 * Revision 1.7.1.17  90/10/23  20:10:47  keith
 * Added inhibit_vectorization() dummy function.
 * 
 * Revision 1.7.1.16  90/10/18  17:24:26  keith
 * Added extra arg to gather() for (conditional) bounds check.
 * C.f. force.c 1.8.1.17.
 * 
 * Revision 1.7.1.15  90/09/28  13:29:08  keith
 * Inserted braces around VECTORIZE directives and changed include files
 * for STARDtardent 3000 series (via cond. comp symbol "ardent").
 * 
 * Revision 1.7.1.14  90/08/30  16:00:21  keith
 * Replaced malloc() calls with calls to own allocator 'aalloc()'
 * Ifdefed out 'remove' and 'memset/memcpy' if ANSI_LIBS is set.
 * 
 * Revision 1.7.1.13  90/05/16  18:39:46  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.7.1.12  90/05/15  19:00:16  keith
 * Incorporate new ANSI CLOCKDS_PER_SEC instead of CLK_TCK.
 * 
 * Revision 1.7.1.11  90/05/02  18:15:12  keith
 * Tidied up include files, added "time.h".
 * 
 * Revision 1.7.1.10  90/04/16  13:04:08  keith
 * Improved sysV version of rt_clock() by using result of times().
 * Coded zero-real explicitly for ANSI compliance.
 * Added strchr() for BSD systems (calls index()).
 * 
 * Revision 1.7.1.9  90/03/27  17:34:46  keith
 * Moved selection of own `vprintf into defs.h - 
 * Define symbol HAVE_VPRINTF if not needed. (Also HAVE_DOPRNT).
 * 
 * Revision 1.7.1.8  90/03/09  17:32:35  keith
 * Modified preprocessor conditionals to handle unicos.
 * 
 * Revision 1.7.1.7  90/01/02  19:01:03  keith
 * Rewrote loops to go from zero to n again - substantially faster on Stellar
 * 
 * Revision 1.7.1.6  89/12/18  17:50:23  keith
 * Sum() tests for special case of stride=1 for efficiency.
 * Rt_clock() added (returns time in seconds).
 * 
 * Revision 1.7.1.5  89/11/21  15:51:10  keith
 * Added purge() for cray and fixed replace() for standard case.
 * 
 * Revision 1.7.1.4  89/10/25  10:05:24  keith
 * Added spaxpy() in support of change in force.c.
 * 
 * Revision 1.7.1.3  89/09/20  17:08:26  keith
 * Added wheneq() and gatheri() to veclib calls for convex.
 * Modified search_lt() for cray to return zero if called with n=0.
 * 
 * Revision 1.7.1.2  89/09/18  17:39:37  keith
 * Fixed error in vprintf which used puts() instead of fputs() giving extra \n.
 * Made 'eps' static in precision() to give correct answer on subsequent calls.
 * 
 * Revision 1.7.1.1  89/08/25  15:21:36  keith
 * Mods to add framework structures to simulation model
 * 
 * Revision 1.7  89/08/15  18:29:58  keith
 * Fixed bug which assumed clock tick for Sys V machines = 1/60s
 * 
 * Revision 1.6  89/08/10  17:29:42  keith
 * Fixed search_lt() to return correct result for non-unit stride.
 * 
 * Revision 1.5  89/06/14  14:14:44  keith
 * Fixed ifdef's for CRAY to handle case of unicos.
 * Fixed mistake in #define's for typedef clash in sysV CPU().
 * 
 * Revision 1.4  89/06/09  17:01:53  keith
 * Made zero_real() call memset/bzero
 * Fixed vprintf() for machines which use _doprnt()
 * Added alliant vector functions
 * Modified cpu() to avoid typedef clash on Stellar GS 1000
 * Rewrote 'scalar' versions of vector functions to allow vectorisation.
 * 
 * Revision 1.3  89/06/01  21:23:04  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.2  89/06/01  18:00:42  keith
 * Moved `vadd()' from aux.c to force.c for ease of vectorisation.
 * Now no need to compile aux.c with vectorisation.
 * 
 * Revision 1.1  89/04/27  15:01:14  keith
 * Initial revision
 * 
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore_data/keith/md/moldy/RCS/auxil.c,v 2.8 1995/12/05 20:55:10 keith Exp $";
#endif
/*========================== program include files ===========================*/
#include	"defs.h"
#ifndef _POSIX_SOURCE
#   define _POSIX_SOURCE
#endif
/*========================== Library include files ===========================*/
#if defined(ANSI) || defined(__STDC__)
#include 	<stdarg.h>
#else
#include 	<varargs.h>
#endif
#include	<ctype.h>
#include	<math.h>
#include 	"string.h"
#include	"stddef.h"
#include	"stdlib.h"
#include	"time.h"
#include	<stdio.h>
/*================= System Library include files - unix only ================*/
#if defined(unix) || defined(__unix__)
   /*
    *  We must protect the inclusion of <sys/types.h>.  
    */
#ifndef SYS_TYPES_INCLUDED
#define SYS_TYPES_INCLUDED
#   include <sys/types.h>
#endif
#include <sys/time.h>
#include <sys/times.h>

#ifndef CLK_TCK
#   include <sys/param.h>
#   ifdef HZ
#      define CLK_TCK HZ
#   else
#      define CLK_TCK 60	/* Really old unices defined 60. May be wrong.*/
#   endif
#endif
#endif
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
/*============================================================================*/
/******************************************************************************
 *  Vector routines - special cases for cray and convex to call vectorised    *
 *  library versions.  Other machines do it the hard way.		      *
 *  CRAY versions.							      *
 ******************************************************************************/
#if defined(_CRAY)
#define VECTS_DEFINED

double vdot(n,x,ix,y,iy)
int	n;
real	*x, *y;
int	ix, iy;
{
   double SDOT();
   return(SDOT(&n, x, &ix, y, &iy));
}
double sum(n,x,ix)
int	n;
real	*x;
int	ix;
{
   double SSUM();
   return(SSUM(&n, x, &ix));
}
void	vscale(n,s,x,ix)
int	n;
double	s;
real	*x;
int	ix;
{
   void SSCAL();
   SSCAL(&n, &s, x, &ix);
}
#ifdef NOT_FOR_NOW
void saxpy(n, sa, sx, ix, sy, iy)
int	n, ix, iy;
double	sa;
double	sx[], sy[];
{
   void SAXPY();
   SAXPY(&n, &sa, sx, &ix, sy, &iy);
}
#endif
#ifndef _CRAYMPP
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   int	ISRCHFLT();
   if( n <= 0 )
      return(0);
   return( ISRCHFLT(&n, x, &ix, &s) - 1);
}
#else
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   int i, im=n*ix;
VECTORIZE
   for( i = (n-1)*ix; i >= 0; i -= ix )
      if( x[i] < s )
	 im = i;
   return(im/ix);
}
#endif
#ifndef _CRAYMPP
                            /*ARGSUSED*/
void	gather(n, a, b, ix, lim)
int	n;
real	a[], b[];
int	ix[];
{
   void GATHER();
   GATHER(&n, a, b+1, ix);
}
void    gatheri(n, a, b, ix)
int     n;
int     a[], b[];
int     ix[];
{
   void GATHER();
   GATHER(&n, a, b+1, ix);
}
#else
                            /*ARGSUSED*/
void	gather(n, a, b, ix, lim)
int	n, lim;
real	a[], b[];
int	ix[];
{
   int i;
VECTORIZE
   for( i = 0; i < n; i++)
   {
#ifdef DEBUG
      if(ix[i] < 0 || ix[i] >= lim )
	 message(0,0,2, "Gather index out of bounds %d %d\n", i, ix[i]);
#endif
      a[i] = b[ix[i]];
   }
}
void	gatheri(n, a, b, ix)
int	n;
int	a[], b[];
int	ix[];
{
   int i;
VECTORIZE
   for( i = 0; i < n; i++)
      a[i] = b[ix[i]];
}
#endif
#ifdef NOT_FOR_NOW
void	scatter(n, a, ix, b)
int	n;
real	a[], b[];
int	ix[];
{
   void SCATTER();
   SCATTER(&n, a+1, ix, b);
}
#endif
#ifndef _CRAYMPP
int wheneq(n, x, ix, s, idx)
int n, ix, idx[];
int x[], s;
{
  int nval = 0, i;
  void WHENEQ();

  WHENEQ(&n, x, &ix, &s, idx, &nval);
  for(i=0; i < nval; i++)
     idx[i]--;
  return( nval );
}
#else
int wheneq(n, x, ix, s, idx)
int n, ix, idx[];
int x[], s;
{
  int nval = 0, i;
 
  for(i = 0; i != n*ix; i += ix)
    if(x[i] == s)
    {
      idx[nval] = i;
      nval++;
   }
  return( nval );
}
#endif

#endif
/******************************************************************************
 *  Vector routines - CONVEX versions.  The calls are to VECLIB functions.    *
 ******************************************************************************/
#if defined(__convexc__)
#define VECTS_DEFINED

#define SINGLE (sizeof(real) == sizeof(float))
double vdot(n,x,ix,y,iy)
int	n;
real	*x, *y;
int	ix, iy;
{
   double ddot_(); float sdot_();
   if(SINGLE)
     return(sdot_(&n, x, &ix, y, &iy));
   else
     return(ddot_(&n, x, &ix, y, &iy));
}
double sum(n,x,ix)
int	n;
real	*x;
int	ix;
{
   double dsum_(); float ssum_();
   if(SINGLE)
      return(ssum_(&n, x, &ix));
   else
      return(dsum_(&n, x, &ix));
}
void	vscale(n,s,x,ix)
int	n;
double	s;
real	*x;
int	ix;
{
   float s2 = s;
   void dscal_(), sscal_();
   if(SINGLE)
      sscal_(&n, &s2, x, &ix);
   else
      dscal_(&n, &s, x, &ix);
}
#ifdef NOT_FOR_NOW
void saxpy(n, sa, sx, ix, sy, iy)
int	n, ix, iy;
double	sa;
double	sx[], sy[];
{
   float s2 = sa;
   void daxpy_(), saxpy_();
   if(SINGLE)
      saxpy_(&n, &s2, sx, &ix, sy, &iy);
   else
      daxpy_(&n, &sa, sx, &ix, sy, &iy);
}
#endif
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   float s2 = s;
   int	i, idsvlt_(), issvlt_();
   if(SINGLE)
      return( (i = issvlt_(&n, x, &ix, &s2)) ? i - 1 : n);
   else
      return( (i = idsvlt_(&n, x, &ix, &s)) ? i - 1 : n);
}
                            /*ARGSUSED*/
void	gather(n, a, b, ix, lim)
int	n;
real	a[], b[];
int	ix[];
{
   void dgthr_(), sgthr_();
   if(SINGLE)
      sgthr_(&n, b+1, a, ix);
   else
      dgthr_(&n, b+1, a, ix);
}
#ifdef NOT_FOR_NOW
void	scatter(n, a, ix, b)
int	n;
real	a[], b[];
int	ix[];
{
   void dsctr_(), ssctr_();
   if(SINGLE)
      ssctr_(&n, b, ix, a+1);
   else
      dsctr_(&n, b, ix, a+1);
}
#endif

void    gatheri(n, a, b, ix)
int     n;
int     a[], b[];
int     ix[];
{
   void igthr_();
   igthr_(&n, b+1, a, ix);
}

int wheneq(n, x, ix, s, idx)
int n, ix, idx[];
int x[], s;
{
   int nval = 0, i;
   void ilsteq_();
   
   ilsteq_(&n, x, &ix, &s, &nval, idx);
   for(i=0; i < nval; i++)
      idx[i]--;
   return( nval );
}
#endif
/******************************************************************************
 * alliant versions - with cray librarries				      *
 ******************************************************************************/
#if defined(alliant)
#define VECTS_DEFINED

double vdot(n,x,ix,y,iy)
int	n;
real	*x, *y;
int	ix, iy;
{
   double sdot_();
   return(sdot_(&n, x, &ix, y, &iy));
}
double sum(n,x,ix)
int	n;
real	*x;
int	ix;
{
   double ssum_();
   return(ssum_(&n, x, &ix));
}
void	vscale(n,s,x,ix)
int	n;
double	s;
real	*x;
int	ix;
{
   void sscal_();
   sscal_(&n, &s, x, &ix);
}
#ifdef NOT_FOR_NOW
void saxpy(n, sa, sx, ix, sy, iy)
int	n, ix, iy;
double	sa;
double	sx[], sy[];
{
   void saxpy_();
   saxpy_(&n, &sa, sx, &ix, sy, &iy);
}
#endif
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   int	isrchflt_();
   return( isrchflt_(&n, x, &ix, &s) - 1);
}
                            /*ARGSUSED*/
void	gather(n, a, b, ix, lim)
int	n;
real	a[], b[];
int	ix[];
{
   void gather_();
   gather_(&n, a, b+1, ix);
}
#ifdef NOT_FOR_NOW
void	scatter(n, a, ix, b)
int	n;
real	a[], b[];
int	ix[];
{
   void scatter_();
   scatter_(&n, a+1, ix, b);
}
#endif

#endif
/******************************************************************************
 *  Vector handling functions - C versions. 		     		      *
 ******************************************************************************/
#ifndef VECTS_DEFINED
double vdot(n,x,ix,y,iy)
int	n;
real	x[], y[];
int	ix, iy;
{
   register double	dot=0.0;
   register int i, j;
   if( ix == iy && ix == 1)
   {
VECTORIZE
      for(i = 0; i < n; i++)
         dot += x[i] * y[i];
   }
   else
   {
VECTORIZE
      for(i = j = 0; i != n*ix; i += ix)
      {
         dot += x[i] * y[j];
         j += iy;
      }
   }
   return(dot);
}
double sum(n,x,ix)
register int	n;
register real	x[];
register int	ix;
{
   register double	s=0.0;
   int i;

   if( ix == 1 )
   {
VECTORIZE
     for( i = 0; i < n; i++)
       s += x[i];
   } 
   else
   {
VECTORIZE
      for(i = 0; i != n*ix; i += ix)
	s += x[i];
   }
   return(s);
}
void	vscale(n,s,x,ix)
register int	n;
register double	s;
register real	x[];
register int	ix;
{
   int i;
VECTORIZE
   for(i = 0; i != n*ix; i += ix)
      x[i] *= s;
}
#ifdef NOT_FOR_NOW
void saxpy(n, sa, sx, ix, sy, iy)
int	n, ix, iy;
double	sa;
double	sx[], sy[];
{
   int i, j;
VECTORIZE
   for(i = j = 0; i != n*ix; i += ix, j += iy)
      sy[j] += sa * sx[i];
}
#endif
/*
 *  Two alternatives for search_lt.  The first can vectorize using an
 *  IDFMIN type function, the second using an isrchtl.  They don't
 *  do QUITE the same thing; the first returns the index of the 
 *  smallest value,, the second that of the first one less than s.
 */
#ifdef ardent
int	search_lt(n, x, ix, s)
int	n;
real 	x[];
int	ix;
double	s;
{
   int i, im=n*ix;
   double min = x[0];
VECTORIZE
   for( i = 0; i != n*ix; i += ix )
      if( x[i] < min )
      {
	 im = i;
	 min = x[i];
      }
   if ( min < s )
      return(im/ix);
   else
      return n;
}
#else
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   int i, im=n*ix;
VECTORIZE
   for( i = (n-1)*ix; i >= 0; i -= ix )
      if( x[i] < s )
	 im = i;
   return(im/ix);
}
#endif
                            /*ARGSUSED*/
void	gather(n, a, b, ix, lim)
int	n, lim;
real	a[], b[];
int	ix[];
{
   int i;
VECTORIZE
   for( i = 0; i < n; i++)
   {
#ifdef DEBUG
      if(ix[i] < 0 || ix[i] >= lim )
	 message(0,0,2, "Gather index out of bounds %d %d\n", i, ix[i]);
#endif
      a[i] = b[ix[i]];
   }
}
void	gatheri(n, a, b, ix)
int	n;
int	a[], b[];
int	ix[];
{
   int i;
VECTORIZE
   for( i = 0; i < n; i++)
      a[i] = b[ix[i]];
}
#ifdef NOT_FOR_NOW
void	scatter(n, a, ix, b)
int	n;
real	a[], b[];
int	ix[];
{
   int i;
VECTORIZE
   for( i = 0; i < n; i++)
      a[ix[i]] = b[i];
}
#endif
int wheneq(n, x, ix, s, idx)
int n, ix, idx[];
int x[], s;
{
  int nval = 0, i;
 
  for(i = 0; i != n*ix; i += ix)
    if(x[i] == s)
    {
      idx[nval] = i;
      nval++;
   }
  return( nval );
}
#endif
/******************************************************************************
 *  random number generators. Note they assume 'unsigned' of at least 32 bits *
 ******************************************************************************/
#define  im  1771875L
#define  ia  2416L
#define  ic  374441L
static unsigned long jran = 1;
void	smdrand(seed)
unsigned long seed;
{
   jran = seed;
}
double	mdrand()
{
   jran = (jran*ia + ic) % im;
   return (double)jran/(double)im;
}
/******************************************************************************
 *  Precision. Return smallest eps s.t. 1.0+eps > 1			      *
 ******************************************************************************/
double	precision()
{
   static	int	first=1;
   static	double	eps = 0.5;
   double VOLATILE	eps2, *eps1 = &eps2;	/* Use pointer to force store */
   double VOLATILE      junk, *ptr = &junk;
   
   if(first)
   {
      first = 0;
      do
      {
         eps = eps * 0.5;
         *eps1 = 1.0 + eps;
	 *ptr  = 1.0;		/* This prevents optimization of *eps1 */
      }   while(*eps1 > 1.0);
      eps *= 2.0;
   }
   return(eps);
}
/******************************************************************************
 *  cpu.  Return (double) the current cpu time in seconds.		      *
 *  rt_clock(). Return elapsed time in s.				      *
 ******************************************************************************/
#if defined(unix) || defined(__unix__)
double        cpu()
{
   struct tms buf;
 
   (void)times(&buf);

   return (buf.tms_utime + buf.tms_stime)/(double)CLK_TCK;
}

/*
 * The BSD version way to get the rt clock is via gettimeofday.  But
 * some systems <sys/time.h>, noticably Ultrix, DON'T include the
 * necessary struct definitions if _POSIX_SOURCE is defined!.  Aargh.
 * However in that case we can rely on POSIX behaviour of times -
 * (ie it's return value) - so the alternative rt_clock ought to
 * work.   The test is for the macro associated with the "timezone"
 * struct definition.
 */
#if defined(BSD) && defined(DST_NONE)

double rt_clock()
{
   int gettimeofday();
   struct timeval tp;
   (void)gettimeofday(&tp,(struct timezone *)0);
   return(tp.tv_sec + tp.tv_usec*0.000001);
}

#else			/* System V or POSIX.			      */

double rt_clock()
{
   struct tms buf;
 
   return times(&buf)/(double)CLK_TCK;
}

#endif			/* USG or BSG					      */
#else			/* Not Unix					      */
#if !defined(CLOCKS_PER_SEC) && defined(CLK_TCK)
#   define  CLOCKS_PER_SEC CLK_TCK
#endif
double	cpu()
{
   return((double)clock() / CLOCKS_PER_SEC);
}
double rt_clock()
{
   return time((time_t *)0);
}
#endif
/******************************************************************************
 *  cctime()  Return ctime() but without the newline			      *
 ******************************************************************************/
char	*cctime(timeloc)
time_mt	*timeloc;
{
   time_t  loc = *timeloc;
   char	*p = ctime(&loc);
   *strchr(p, '\n') = '\0';
   return(p);
}
/******************************************************************************
 *  atime.  Return a pointer to the time and date string		      *
 ******************************************************************************/
char	*atime()
{
   time_mt t = time((time_t*)0);
   return(cctime(&t));
}
/******************************************************************************
 * Zero real.  Set array to floating-point zero.			      *
 ******************************************************************************/
void	zero_real(r,n)
real	r[];
int	n;
{
   int i;
VECTORIZE
   for(i = 0; i < n; i++)
      r[i] = 0.0;
}
void	zero_double(r,n)
double	r[];
int	n;
{
   int i;
VECTORIZE
   for(i = 0; i < n; i++)
      r[i] = 0.0;
}
void	zero_dbls(r,n)
double	r[];
size_mt	n;
{
   long i;
VECTORIZE
   for(i = 0; i < n; i++)
      r[i] = 0.0;
}
/******************************************************************************
 *  replace - replace file1 by file2, renaming or overwriting		      *
 ******************************************************************************/
#if defined(unix) || defined(__unix__)
int	replace(file1, file2)
char	*file1, *file2;
{
   int f;
   char *backup = aalloc(strlen(file2)+2, char);
   (void)strcat(strcpy(backup, file2), "%");

   f = rename(file2, backup);		/* Backup old version - ignore failure*/
   f = rename(file1, file2);		/* Rename old version to new	      */
   xfree(backup);
   return(f);
}
#   else
int	replace(file1, file2)
char	*file1, *file2;
{
   int f;
   int rename();
   f = rename(file1, file2);
   if(f != 0 && remove(file2) == 0)	/* If already exists try remove target*/
      f = rename(file1, file2);
   return(f);
}
#endif
/******************************************************************************
 *  Purge Remove 1 old version of a file.  (Do nothing if O/S has no versions)*
 ******************************************************************************/
#ifdef vms 
void	purge(file)
char	*file;
{
   char *name = aalloc(strlen(file)+4,char);
   strcpy(name, file);
   if( strchr(name, ';') == NULL)
      (void)remove(strcat(name,";-1"));
   tfree(name);
}
#else				/*  VMS				*/
#   if defined(unix) || defined(__unix__)
void	purge(file)
char	*file;
{
   int unlink();
   char *name = aalloc(strlen(file)+2,char);
   (void)strcat(strcpy(name, file), "%");

   (void)unlink(name);
   xfree(name);
}
#      else
/*ARGSUSED*/
void	purge(file)
char	*file;
{
}
#   endif			
#endif				/*  VMS				*/
/******************************************************************************
 * err_fn  error function. See Abramowitz & Stegun p299.		      *
 ******************************************************************************/
#define E1	 0.254829592
#define E2	-0.284496736
#define E3	 1.421413741
#define E4	-1.453152027
#define E5	 1.061405429

#define PP	 0.3275911
#define POLY5(t)   ((t)*(E1 + (t)*(E2 + (t)*(E3 + (t)*(E4 + (t)*E5)))))

double	err_fn(x)
double	x;
{
   double t = 1.0/(1.0 + PP * x);
   if ( x >= 0.0 )
      return( 1.0 - POLY5(t) * exp(-x*x) );
   else
      return( -err_fn(-x) );
}
/******************************************************************************
 * inhibit_vectorization().  Self explanatory dummy function		      *
 ******************************************************************************/
void inhibit_vectorization()
{
}
