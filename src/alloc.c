/******************************************************************************
 * Alloc	Interface functions to dynamic store allocators.	      *
 * talloc()	Allocate storage and test for success or failure	      *
 * arralloc()	Allocate rectangular dope-vector (ie using pointers) array    *
 *									      *
 * N.B.         Portability.						      *
 *    These functions make some assumptions which are not guaranteed by the   *
 * ANSI standard.							      *
 * 1)   Various other modules rely on talloc() setting the store to floating- *
 *      point zero.  This may not be the same as binary zero.  In that case   *
 *      the macros "ralloc" etc in "defs.h" will have to be converted into    *
 *      functions here which do type-dependant initialisation to zero.        *
 * 2)   arralloc() relies on a common format for pointers to different data   *
 *      types, and assumes that the representation of a "data" pointer is the *
 *      same as of an integer pointer.  (N.B.  it can not be used to allocate *
 *      character data). It will not work (and cannot be made to work) on     *
 * 	machines for which this is false. 				      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	alloc.c,v $
 * Revision 1.7  90/05/02  17:51:09  keith
 * Include of stddef.h added to get size_t (removed from defs.h)
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
 * Changed ifdef's to select on __STDC__ macro
 * 
 * Revision 1.1  89/04/27  16:52:17  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/alloc.c,v 1.7 90/05/02 17:51:09 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#if ANSI || __STDC__
#   include <stdarg.h>
#else
#   include <varargs.h>
#endif
#include "stddef.h"
/*========================== program include files ===========================*/
#include "defs.h"
#include "messages.h"
/*========================== Library declarations ============================*/
char	*calloc();
int	free();
/*========================== External function declarations ==================*/
void	message();				/* Error handling routine     */
#ifdef	DEBUG
int	malloc_verify();
int	malloc_debug();
#endif
/*============================================================================*/
/******************************************************************************
 * talloc()	Call Calloc to allocate memory, test result and stop if failed*
 ******************************************************************************/
char	*talloc(n, size, line, file)
long	n;
size_t	size;
int	line;
char	*file;
{
   char *p = calloc((unsigned) n, (unsigned) size);
   if(p == NULL && (n*size != 0))
      message(NULLI, NULLP, FATAL, NOMEM, line, file,
	      (int)n, (unsigned long)size);
   return(p);
}
/******************************************************************************
 * Cfree - synonym to free()						      *
 ******************************************************************************/
void	cfree(p)
char	*p;
{
#ifdef DEBUG
   if( ! malloc_verify() )
      message(NULLI, NULLP, FATAL, "Internal Error: Heap corrupt");
#endif
   free(p);
}
/******************************************************************************
 *  arralloc.   Allocate a psuedo array of any dimensionality and type with   *
 *  specified lower and upper bounds for each dimension.  Each dimension is   *
 *  an array of pointers, and the actual data is laid out in standard 'c'     *
 *  fashion ie last index varies most rapidly.  All storage is got in one     *
 *  block, and so can be freed in one go.  				      *
 *  array = (double*) arralloc(sizeof(double), 3, 0, 10, -10, 10, 0, 5);      *
 *  cfree((char*) array);					     	      *
 *  (N.B. if lower bound of 1st dimension != 0 then free array+l.b.           *
 ******************************************************************************/
#if ANSI || __STDC__
#define	va_alist size_t size, int ndim, ...
#define va_dcl /* */
#endif

#define	MAXDIM	11
                 /*VARARGS*/
typedef char	canon_t;
canon_t		*arralloc(va_alist)
va_dcl
{
   int		lb[MAXDIM], ub[MAXDIM];
   va_list	ap;
   long		stride, idim, n_ptr, n_data;
   typedef	int* 	ptr_t;
   ptr_t	*pointer, *end, pointed, start;
#if ANSI || __STDC__
   va_start(ap, ndim);
#else
   size_t	size;			/* size of array element	      */
   int		ndim;			/* Number of dimensions		      */

   va_start(ap);
   size = va_arg(ap, size_t);
   ndim = va_arg(ap, int);
#endif
   if(ndim > MAXDIM)  message(NULLI, NULLP, FATAL, TOODIM, ndim, MAXDIM);
   for(idim = 0; idim < ndim; idim++)
   {
      lb[idim] = va_arg(ap, int);
      ub[idim] = va_arg(ap, int);
   }
   va_end(ap);
 
   n_ptr = 0; n_data = 1;
   for(idim = 0; idim < ndim; idim++)
   {
#ifdef DEBUG
      printf("[%d...%d]", lb[idim], ub[idim]);
#endif
      if(ub[idim] < lb[idim])
         message(NULLI, NULLP, FATAL, INSIDE, lb[idim], ub[idim]);
      n_data *= ub[idim] - lb[idim] + 1;
   }
#ifdef DEBUG
   putchar('\n');
#endif
   for(idim = ndim-2; idim >= 0; idim--)
      n_ptr = (1 + n_ptr) * (ub[idim] - lb[idim] + 1);
 
   if( size % sizeof(ptr_t*) != 0 )
      message(NULLI, NULLP, FATAL, WDPTR, size);

   start = (ptr_t)talloc(1L, (size_t)(n_data*size + (n_ptr+1)*sizeof(ptr_t)),
			   __LINE__, __FILE__);
   size /= sizeof(ptr_t*);		/* Size is now in terms of "words"   */


   pointer = (ptr_t *)start;
   pointed = (ptr_t)(pointer + ub[0] - lb[0] + 1);
   /*
    *  Set up pointers for all but last dimension
    */
   for(idim = 1; idim < ndim-1 ; idim++)
   {
      end = (ptr_t *)pointed;
      stride = ub[idim] - lb[idim] + 1;
      for(; pointer < end; pointer++)
      {
	 *pointer = pointed - lb[idim];
	 pointed += stride;
      }
   }
   /*
    *  Pointers to last dimension
    */
   if( idim == ndim - 1 )
   {
      end = (ptr_t *)pointed;
      pointed = start + ((pointed-start)+size-1)/size*size;
      stride = ub[idim] - lb[idim] + 1;
      for(; pointer < end; pointer++)
      {
	 *pointer = pointed - lb[idim]*size;
	 pointed += stride*size;
      }
   }
   if(ndim > 1)
      return (canon_t*)(start - lb[0]);
   else
      return (canon_t*)(start - lb[0]*size);
}
