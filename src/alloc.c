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
 * Revision 1.9  90/05/16  18:39:36  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.8  90/05/16  14:19:23  keith
 * *** empty log message ***
 * 
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
static char *RCSid = "$Header: /mnt/keith/moldy/RCS/alloc.c,v 1.10 90/08/14 14:23:22 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#if ANSI || __STDC__
#   include <stdarg.h>
#else
#   include <varargs.h>
#endif
#include "stddef.h"
#include "stdlib.h"
/*========================== program include files ===========================*/
#include "defs.h"
#include "messages.h"
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
void	tfree(p)
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
 *  tfree((char*) array);					     	      *
 *  (N.B. if lower bound of 1st dimension != 0 then free array+l.b.           *
 ******************************************************************************/
#define CSA(a) ((char*)(a))
#define ALIGN(a,base,b)	((int*)(CSA(base)+((CSA(a)-CSA(base))+(b)-1)/(b)*(b) ))

void 	subarray(size, ndim, prdim, pp, qq, base, ap)
size_t  size;
int	ndim;
long	prdim;
int	***pp, **qq, *base;
va_list	ap;
{
   int	*dd = ALIGN(qq,base,size),	**dpp = (int**)pp;
   int i,	lb = va_arg(ap, int),
		dim = va_arg(ap, int) - lb + 1;

   if(ndim > 0)		/* General case - set up pointers to pointers  */
   {
      for( i = 0; i < prdim; i++)
	 pp[i] = qq + i*dim - lb;

      subarray(size, ndim-1, prdim*dim, (int***)qq, qq+prdim*dim, base, ap);
   }
   else			/* Last recursion - set up pointers to data   */
      for( i = 0; i < prdim; i++)
	 dpp[i] = dd + (i*dim - lb)*size/sizeof(int);
}
            
#if ANSI || __STDC__
#define	va_alist size_t size, int ndim, ...
#define va_dcl /* */
#endif
                /*VARARGS*/
char		*arralloc(va_alist)
va_dcl
{
   va_list	ap, ap2;
   int		**p, **start;
   int		lb, ub, idim;
   long		n_ptr = 0, n_data = 1;
#if ANSI || __STDC__
   va_start(ap, ndim);
#else
   size_t	size;			/* size of array element	      */
   int		ndim;			/* Number of dimensions		      */

   va_start(ap);
   size = va_arg(ap, size_t);
   ndim = va_arg(ap, int);
#endif

   if( size % sizeof(int) != 0 )  /* Code only works for 'word' objects */
      message(NULLI, NULLP, FATAL, WDPTR, size);
   /*
    * Cycle over dims, checking bounds and accumulate # pointers & data items.
    */
   ap2 = ap;			/* Save ap for later use by subarray() */
   for(idim = 0; idim < ndim; idim++)
   {
      lb = va_arg(ap, int);
      ub = va_arg(ap,int);
      if(ub < lb)
         message(NULLI, NULLP, FATAL, INSIDE, lb, ub);
      n_data *= ub - lb + 1;
      if( idim < ndim-1 )
	 n_ptr  += n_data;
   }
   /*
    *  Allocate space  for pointers and data.
    */
   start = (int**)talloc(1L, (size_t)((n_data+1)*size+n_ptr*sizeof(void**)),
			  __LINE__, __FILE__);
   /*
    * Set up pointers to form dope-vector array.
    */
   subarray(size, ndim-1, 1L, &p, start, (int*)start, ap2);

   va_end(ap);   

   return (char*)p;
}
