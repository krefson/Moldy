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
 */
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
	       (int)n, (n==1)?" ":"s ", (unsigned long)size);
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
