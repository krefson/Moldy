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
 * Alloc	Interface functions to dynamic store allocators.	      *
 * talloc()	Allocate storage and test for success or failure	      *
 * arralloc()	Allocate rectangular dope-vector (ie using pointers) array    *
 *									      *
 * N.B.         Portability.						      *
 *    These functions make some assumptions which are not guaranteed by the   *
 * ANSI standard.							      *
 * 1)   arralloc() relies on a common format for pointers to different data   *
 *      types, and assumes that the representation of a "data" pointer is the *
 *      same as of an integer pointer.  (N.B.  it can not be used to allocate *
 *      character data). It will not work (and cannot be made to work) on     *
 * 	machines for which this is false. 				      *
 ******************************************************************************
 *      Revision Log
 *       $Log$
 *       Revision 2.16.4.1  2003/08/15 09:13:10  moldydv
 *       Modified NOMEM error message to handle both singular and plural.
 *
 *       Revision 2.16  2002/09/19 09:26:27  kr
 *       Tidied up header declarations.
 *       Changed old includes of string,stdlib,stddef and time to <> form
 *
 *       Revision 2.15  2001/02/14 15:19:28  keith
 *       Rejigged definition of word_mt in alloc.c for machines where
 *       sizeof(float) < sizeof(int). (ie Cray T3E)
 *
 *       Revision 2.14  2000/12/06 17:45:28  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.13  2000/10/20 15:15:46  keith
 *       Incorporated all mods and bugfixes from Beeman branch up to Rel. 2.16
 *
 *       Revision 2.12  2000/04/27 17:57:05  keith
 *       Converted to use full ANSI function prototypes
 *
 *       Revision 2.11.2.1  2000/08/11 17:40:03  keith
 *       Workaround for Portland group C compiler bug.
 *
 *       Revision 2.11  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 *
 *       Revision 2.10  1998/01/27 15:41:42  keith
 *       Got rid of copying of arg ptr ap and replaced with a rescan of
 *       the argument list in 'arralloc()'.  This is because va_list
 *       isn't necessarily an lvalue in ANSI/ISO C (and really isn't in
 *       WATCOM C.
 *
 *       Revision 2.9  1997/11/27 16:24:49  keith
 *       Removed titan-specific THREAD_SYS stuff for cleanliness.
 *
 *       Revision 2.8  1994/06/08 13:08:08  keith
 *       New version of array allocator which breaks up requests for DOS.
 *       Now must use specific "afree()" paired with arralloc().
 *
 * Revision 2.6  1994/02/21  16:55:58  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Added sanity test for 16-bit architectures.
 *
 * Revision 2.5  1994/01/18  13:32:09  keith
 * Null update for XDR portability release
 *
 * Revision 2.3  93/10/28  10:27:37  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.0  93/03/15  14:48:54  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.22  93/03/15  14:41:32  keith
 * Added GPL copyleft notice to permit release and distribution.
 * N.B.  Previous versions were copyright (C) by author and 
 * only licensed by explicit permission.
 * 
 * Revision 1.21  93/03/12  20:11:49  keith
 * Fixed mistake of typing word_mt to double -- must be int.
 * Documented it better.
 * 
 * Revision 1.20  93/03/09  15:58:16  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.19  93/03/05  15:00:30  keith
 * Added generic type "word_t" for portability.
 * Moved include line for security in non-ansi gcc environments
 * 
 * Revision 1.18  92/06/26  17:02:42  keith
 * Got rid of assumption that memory returned by talloc() or
 * arralloc() is zeroed.  This enhances ANSI compatibility.
 * Removed memory zeroing from alloc.c() in consequence.
 * 
 * Revision 1.17  92/06/05  13:37:02  keith
 * Conditionally undefed va_dcl for ANSI, stdarg.h case --
 * just prevents warning from gcc.
 * 
 * Revision 1.16  91/10/17  14:22:21  keith
 * Added debugging code
 * 
 * Revision 1.15  91/08/19  16:44:11  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_mt and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * 
 * Revision 1.14  91/03/12  15:42:10  keith
 * Tidied up typedefs size_mt and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.13  91/03/07  17:52:32  keith
 * Macros in support of parallel version for titan added.
 * 
 * Revision 1.12  90/10/23  20:13:14  keith
 * Added dummy function call to inhibit vectorization.
 * This allows use of 'ivdep' compiler options and also
 * works round certain bugs in cray's scc compiler.
 * 
 * Revision 1.11  90/08/29  18:21:09  keith
 * Replaced calloc() call with malloc() and memset().
 * On the CRAY XMP calloc() is very inefficient.
 * 
 * 
 * Revision 1.10  90/08/23  12:46:39  keith
 * Re-implemented arralloc() using more elegant recursive algorithm.
 * 
 * Revision 1.9  90/05/16  18:39:36  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.8  90/05/16  14:19:23  keith
 * *** empty log message ***
 * 
 * Revision 1.7  90/05/02  17:51:09  keith
 * Include of stddef.h added to get size_mt (removed from defs.h)
 * 
 * Revision 1.6  90/04/25  14:20:16  keith
 * Modified to allow for machines with word ptr != char ptr.
 * 
 * Revision 1.5  90/03/26  16:54:46  keith
 * Added portability warning to header comments.
 * 
 * Revision 1.4  90/01/15  17:22:25  keith
 * New version of arralloc() orders memory so that pointers come FIRST.
 * This means you can simply free() the pointer returned (if l.b. = 0).
 * 
 * Revision 1.3  89/09/21  14:56:01  keith
 * Modified talloc() to return null rather than exit if 0 bytes requested.
 * 
 * Revision 1.2  89/05/24  13:54:26  keith
 * Changed ifdef's to select on defined(__STDC__) macro
 * 
 * Revision 1.1  89/04/27  16:52:17  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header$";
#endif
/*========================== program include files ===========================*/
#include "defs.h"
#include "messages.h"
/*========================== Library include files ===========================*/
#include <stdarg.h>
#include <stdlib.h>
#if defined(DEBUGX) || defined (DEBUGY)
#include <stdio.h>
#endif
#ifdef DEBUGZ
#include <string.h>
#endif
#ifdef DBMALLOC
#include	<stddef.h>
typedef size_mt size_t;
#include <dbmalloc.h>
#endif
/*========================== External function declarations ==================*/
void	inhibit_vectorization(void);		/* Self-explanatory dummy     */
#ifdef	DEBUG
int	malloc_verify();
int	malloc_debug();
#endif
void	note(char *, ...);		/* Write a message to the output file */
void	message(int *, ...);		/* Write a warning or error message   */
/*============================================================================*/
/******************************************************************************
 * word_mt								      *
 * This is the core of the non-ANSI conformance of "arralloc" and *must* be   *
 * right for portability.  Word_mt must be typed to the *smallest* size of    *
 * object to be allocated.  On a word-addressed machine  it ought to be the   *
 * smallest addressible object.  Char is only safe for byte-addressible       *
 * architectures (with the exception of the Cray, which works).  Moldy doesn't*
 * call it for anything smaller than an int, which being the "naturally" sized*
 * object is the optimum.  If pointer representations of any actually required*
 * object (ie NOT char) do differ then these functions will have to be        *
 * rewritten.								      *
 *									      *
 * Wide_mt is the widest type for alignment purposes.  Try double.	      *
 ******************************************************************************/
#if SIZEOF_FLOAT < SIZEOF_INT
typedef float word_mt;
#else
typedef int word_mt;
#endif
#ifdef ALLOC_SEPARATELY
typedef double wide_mt;
#endif
/******************************************************************************
 * talloc()	Call Calloc to allocate memory, test result and stop if failed*
 ******************************************************************************/
gptr	*talloc(int n, size_mt size, int line, char *file)
{
   gptr *p;

#ifdef STDC_HEADERS
   /*
    * Test for malloc arg which would overflow.  Since size_mt is long
    * and size_t may be int this could happen on 16 bit machines.
    * We can only rely on size_t as parameter to malloc if libs are
    * ANSI conformant.
    */
   if( (size_mt)(size_t)(n*size) != n*size )
      message(NULLI, NULLP, FATAL, NOMEM, line, file,
	       (int)n, (n==1)?' ':'s', (unsigned long)size);
#endif
   p = malloc(n*size);
#ifdef DEBUGX
   fprintf(stderr,"Alloc: %16s line %3d: %d x %lu bytes (%p to %p)\n", 
	   file, line, n, size, p, p+n*size);
#endif
   if(p == 0 && (n*size != 0))
     message(NULLI, NULLP, FATAL, NOMEM, line, file,
	     (int)n, (n==1)?" ":"s ", (unsigned long)size);
#ifdef DEBUGZ
   (void)memset((gptr*)p,0x10,n*size);
#endif
   return(p);
}
/******************************************************************************
 * Cfree - synonym to free()						      *
 ******************************************************************************/
void	tfree(gptr *p)
{
#ifdef DEBUG
   if( ! malloc_verify() )
      message(NULLI, NULLP, FATAL, "Internal Error: Heap corrupt");
#endif
   free((gptr*)p);
}

#ifdef ALLOC_SEPARATELY
union u {struct {int ndim, noffset, len;} b; word_mt * p; wide_mt w;};
#define bsize (sizeof(union u)/sizeof(word_mt*))
#define bwsize (sizeof(union u)/sizeof(word_mt))
void afree(pp)
gptr *pp;
{
   word_mt **p = (word_mt**)pp;
   int i; union u *up = (union u *)(p-bsize);
   if( up->b.ndim > 1 )
   {
      for( i=0; i < up->b.len; i++)
	 afree((gptr*)(p[i] + up->b.noffset));
   }
   tfree((gptr*)(p-bsize));
}
      
#else /* ALLOC_SEPARATELY*/
void 	afree(gptr *p)
{
   tfree(p);
}
#endif /* ALLOC_SEPARATELY*/
/******************************************************************************
 *  arralloc.   Allocate a psuedo array of any dimensionality and type with   *
 *  specified lower and upper bounds for each dimension.  Each dimension is   *
 *  an array of pointers, and the actual data is laid out in standard 'c'     *
 *  fashion ie last index varies most rapidly.  All storage is got in one     *
 *  block, and so can be freed in one go.  				      *
 *  array = (double*) arralloc(sizeof(double), 3, 0, 10, -10, 10, 0, 5);      *
 *  xfree(array);					     	      *
 *  (N.B. if lower bound of 1st dimension != 0 then free array+l.b.           *
 ******************************************************************************/
#ifdef ALLOC_SEPARATELY

word_mt **subarray(size_mt size, int ndim, int len, va_list ap)
{
   word_mt **p;
   word_mt  *d;
   union u *up;
   int blen, i, lb = va_arg(ap, int), ub = va_arg(ap,int);

   if( ndim > 1 )
   {
#ifdef DEBUGY
      fprintf(stderr,"[%d...%d]", lb, ub);
#endif
      blen = len+bsize;
      p = (word_mt**)talloc(blen, (size_mt)sizeof(word_mt *), 
			    __LINE__, __FILE__);
      up = (union u *)p;      p += bsize;
      up->b.ndim = ndim;	
      up->b.len = len;	
      up->b.noffset = lb*(ndim>2?sizeof(word_mt*):size)/sizeof(word_mt);
      for( i=0; i<len; i++)
	 p[i] = (word_mt*)subarray(size, ndim-1, ub-lb+1, ap) - up->b.noffset;
      return p;
   } else 
   {
      blen = len*(size/sizeof(word_mt))+bwsize;
      d = (word_mt*)talloc(blen, (size_mt)sizeof(word_mt), __LINE__, __FILE__);
      up = (union u *)d;      d += bwsize;
      up->b.ndim = ndim;
      return (word_mt **)d;
   }      
}

                /*VARARGS*/
gptr		*arralloc(size_mt size, int ndim, ...)

{
   va_list	ap;
   word_mt		*p;
   int		lb, ub;
   va_start(ap, ndim);

#ifdef DEBUGY
   fprintf(stderr,"%dD array of %lu byte elements:", ndim, size);
#endif
   if( size % sizeof(word_mt) != 0 )  /* Code only works for 'word' objects */
      message(NULLI, NULLP, FATAL, WDPTR, size);

   lb = va_arg(ap, int); ub = va_arg(ap, int);
#ifdef DEBUGY
   fprintf(stderr,"[%d...%d]", lb, ub);
#endif
   
   p=(word_mt*)subarray(size, ndim, ub-lb+1, ap) 
      - lb*(ndim>1?sizeof(word_mt*):size)/sizeof(word_mt);

#ifdef DEBUGY
   putc('\n',stderr);
#endif

   va_end(ap);   

   return (gptr*)p;
}
#else /* ALLOC_SEPARATELY*/
#define CSA(a) ((char*)(a))
#define ALIGN(a,base,b)	((word_mt*)(CSA(base)+((CSA(a)-CSA(base))+(b)-1)/(b)*(b) ))
word_mt w;
static
void 	subarray(size_mt size, int ndim, long int prdim, word_mt ***pp, word_mt **qq, word_mt *base, va_list ap)
{
   word_mt	*dd = ALIGN(qq,base,size),	**dpp = (word_mt**)pp;
   int i,	lb = va_arg(ap, int),
		dim = va_arg(ap, int) - lb + 1;

   if(ndim > 0)		/* General case - set up pointers to pointers  */
   {
      for( i = 0; i < prdim; i++)
      {
	 pp[i] = qq + i*dim - lb;
      }

      subarray(size, ndim-1, prdim*dim, (word_mt***)qq, qq+prdim*dim, base, ap);
      w=base[0];
   }
   else			/* Last recursion - set up pointers to data   */
      for( i = 0; i < prdim; i++)
	 dpp[i] = dd + (i*dim - lb)*size/sizeof(word_mt);
}
            
                /*VARARGS*/
gptr		*arralloc(size_mt size, int ndim, ...)

{
   va_list	ap;
   word_mt		**p, **start;
   int		lb, ub, idim;
   long		n_ptr = 0, n_data = 1;
   va_start(ap, ndim);

#ifdef DEBUGY
   fprintf(stderr,"%dD array of %lu byte elements:", ndim, size);
#endif
   if( size % sizeof(word_mt) != 0 )  /* Code only works for 'word' objects */
      message(NULLI, NULLP, FATAL, WDPTR, size);
   /*
    * Cycle over dims, checking bounds and accumulate # pointers & data items.
    */
   for(idim = 0; idim < ndim; idim++)
   {
      lb = va_arg(ap, int);
      ub = va_arg(ap,int);
#ifdef DEBUGY
      fprintf(stderr,"[%d...%d]", lb, ub);
#endif
      if(ub < lb)
         message(NULLI, NULLP, FATAL, INSIDE, lb, ub);
      n_data *= ub - lb + 1;
      if( idim < ndim-1 )
	 n_ptr  += n_data;
   }
#ifdef DEBUGY
   putc('\n',stderr);
#endif
   /*
    *  Allocate space  for pointers and data.
    */
#ifdef DEBUGY
   {
      size_mt mallen=(n_data+1)*size+n_ptr*sizeof(word_mt**);
      fprintf(stderr,"Calling talloc(%lu)....\n", mallen);
   }
#endif
   start = (word_mt**)talloc(1,
			     (size_mt)((n_data+1)*size+n_ptr*sizeof(word_mt**)),
		
	     __LINE__, __FILE__);
   /*
    * Rescan argument list to pass to subarray()
    */
   va_end(ap);
   va_start(ap, ndim);

   /*
    * Set up pointers to form dope-vector array.
    */
   subarray(size, ndim-1, 1L, &p, start, (word_mt*)start, ap);

   va_end(ap);   

   return (gptr*)p;
}
#endif /* ALLOC_SEPARATELY*/
