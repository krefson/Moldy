/******************************************************************************
 * Alloc	Interface functions to dynamic store allocators.	      *
 * talloc()	Allocate storage and test for success or failure	      *
 * arralloc()	Allocate rectangular dope-vector (ie using pointers) array    *
 ******************************************************************************
 *      Revision Log
 *       $Log:	alloc.c,v $
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
static char *RCSid = "$Header: /home/tigger/keith/md/moldy/RCS/alloc.c,v 1.4 89/12/22 19:30:29 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#if ANSI || __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif
/*========================== program include files ===========================*/
#include "defs.h"
#include "messages.h"
/*========================== Library declarations ============================*/
char	*calloc();
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
      message(NULLI, NULLP, FATAL, MEMORY, line, file,
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
char	*arralloc(va_alist)
va_dcl
{
   int		lb[MAXDIM], ub[MAXDIM];
   va_list	ap;
   long		xsize, stride, idim, n_ptr, n_data;
   char		**pointer, **end, *pointed, *start;
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
 
   start = talloc(1L, (size_t)(n_data*size + (n_ptr+1)*sizeof(char*)),
		  __LINE__, __FILE__);

   pointer = (char**)start;
   pointed = (char*)(pointer + ub[0] - lb[0] + 1);
   xsize = sizeof (char*);
   for(idim = 1; idim < ndim ; idim++)
   {
      end = (char**)pointed;
      if( idim == ndim - 1 ) 		/* Last iteration - set up pointers   */
      {					/* to data this time.		      */
	 xsize = size;			/* Increment is now sizeof(data)      */
	 pointed = start + ((pointed-start)*sizeof(char)+size-1)/size*size;
      }					/* Make sure 'pointed' is data-aligned*/
      stride = ub[idim] - lb[idim] + 1;
      for(; pointer < end; pointer++)
      {
	 *pointer = pointed - lb[idim]*xsize;
	 pointed += stride*xsize;
      }
   }
   return (start - lb[0]*(ndim > 1 ? sizeof(char*) : size));
}
