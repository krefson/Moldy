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
 * matrix	Functions for manipulating 3x3 matrices			      *
 *		Contents:						      *
 * mat_vec_mul()	Multiply matrix by array of vectors		      *
 * mat_mul()		Multiply two 3x3 matrices			      *
 * mat_sca_mul()	Multiply 3x3 matrix by scalar			      *
 * mat_add()		Add two 3x3 matrices				      *
 * transpose()		Transpose a 3x3 matrix				      *
 * det()		Return determinant of 3x3 matrix		      *
 * invert()		Invert a 3x3 matrix				      *
 ******************************************************************************
 *      Revision Log
 *       $Log: matrix.c,v $
 * Revision 2.5  1994/01/18  13:32:44  keith
 * Null update for XDR portability release
 *
 * Revision 2.3  93/10/28  10:28:00  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.1  93/08/18  20:54:18  keith
 * Tidied up clashes over ABS, MIN, MAX macros.
 * 
 * Revision 2.0  93/03/15  14:49:14  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.7  93/03/09  15:59:01  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.6  91/08/15  18:12:08  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.5  90/09/28  13:29:43  keith
 * Inserted braces around VECTORIZE directives and changed include files
 * for STARDtardent 3000 series (via cond. comp symbol "ardent").
 * 
 * Revision 1.4  89/12/18  17:54:23  keith
 * Rewrote transpose() to compile correctly under -va option (for Stellar).
 * 
 * Revision 1.3  89/10/26  11:26:55  keith
 * Mat_vec_mul() vectorised.
 * 
 * Revision 1.2  89/10/24  17:16:21  keith
 * Modified transpose() so as not to vectorise under '-va' option on convex.
 * 
 * Revision 1.1  89/04/20  16:00:50  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/matrix.c,v 2.5 1994/01/18 13:32:44 keith Stab $";
#endif
/*========================== Program include files ===========================*/
#include 	"defs.h"
#include 	"messages.h"
/*========================== Library include files ===========================*/
#include	"string.h"
/*========================== External function declarations ==================*/
#if defined(ANSI) || defined(__STDC__)
void		message(int *, ...);	/* Write a warning or error message   */
#else
void		message();		/* Write a warning or error message   */
#endif
/*============================================================================*/
#define ABS(x)		((x) > 0 ? (x) : -(x))
/******************************************************************************
 * 3 x 3  Matrix - vector multiply  (of multiple vectors)                     *
 * The input and output vectors need not necessarily be distinct              *
 ******************************************************************************/
void mat_vec_mul(m, in_vec, out_vec, number)
int		number;		/* Number of vectors to be multiplied         */
mat_mt		m;		/* Matrix                                     */
vec_mp		in_vec,		/* Input vector, [number][3]          (in/out)*/
		out_vec;	/* Output vector.  CAN BE SAME AS INPUT  (out)*/
{
   int i;
   register real	a0, a1, a2;
   
   /*
    * There are two blocks of code depending on whether the input and output  *
    * arrays are distinct.  If they are not, in_vec is explicitly assigned    *
    * otherwise there could be problems with optimising compilers. Also,      *
    * temporary storage is required to avoid writing over the input vector    *
    * and trying to re-use the 'old' values.
    */
   if(in_vec == out_vec)	/* ie parameters point to the same array      */
   {
VECTORIZE
      for(i = 0; i < number; i++)
      {
         a0 = in_vec[i][0];  a1 = in_vec[i][1];  a2 = in_vec[i][2];

         in_vec[i][0] = m[0][0]*a0 + m[0][1]*a1 + m[0][2]*a2;
         in_vec[i][1] = m[1][0]*a0 + m[1][1]*a1 + m[1][2]*a2;
         in_vec[i][2] = m[2][0]*a0 + m[2][1]*a1 + m[2][2]*a2;
      }
   }
   else 				 /* parameters are distinct arrays  */
   {
VECTORIZE
      for(i = 0; i < number; i++)
      {
         a0 = in_vec[i][0];  a1 = in_vec[i][1];  a2 = in_vec[i][2];

         out_vec[i][0] = m[0][0]*a0 + m[0][1]*a1 + m[0][2]*a2;
         out_vec[i][1] = m[1][0]*a0 + m[1][1]*a1 + m[1][2]*a2;
         out_vec[i][2] = m[2][0]*a0 + m[2][1]*a1 + m[2][2]*a2;
      }
   }
#if 0
   /* Can not test for overlap portably in ANSI C.			*/
   else			/* in_vec and out_vec overlap, but are not identical  */
      message(NULLI, NULLP, FATAL, OVRLP1, in_vec, out_vec, number);
#endif
}
/******************************************************************************
 * mat_mul   Multiply two 3 x 3 matrices. Result can NOT overwrite input.     *
 ******************************************************************************/
void mat_mul(a, b, c)
mat_mt	a,					/* Input matrix 1        (in) */
	b,					/* Input matrix 2        (in) */
	c;					/* Result matrix        (out) */
{
   register int	i, j;				/* Counters		      */
   
   if(c == a || c == b)
      message(NULLI, NULLP, FATAL, OVRLAP, "mat_mul");

   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
         c[i][j] = a[i][0]*b[0][j] + a[i][1]*b[1][j] + a[i][2]*b[2][j];
}
/******************************************************************************
 *  mat_sca_mul.  Multiply a 3x3 matrix by a scalar                           *
 ******************************************************************************/
void mat_sca_mul(s, a, b)
register real	s;			/* Scalar			(in)  */
mat_mt	a, 				/* Input matrix			(in)  */
	b;				/* Result matrix	       (out)  */
{	 
   register int i, j;
   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
         b[i][j] = s * a[i][j];
}            
/******************************************************************************
 * mat_add   Add two 3 x 3 matrices.                                          *
 ******************************************************************************/
void mat_add(a, b, c)
mat_mt	a,					/* Input matrix 1        (in) */
	b,					/* Input matrix 2        (in) */
	c;					/* Result matrix        (out) */
{
   register int	i, j;				/* Counters		      */
   
   for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
         c[i][j] = a[i][j] + b[i][j];
}
/******************************************************************************
 * Transpose  Transpose a 3 x 3 matrix.  Will handle case of a = b            *
 ******************************************************************************/
void transpose(a, b)
mat_mt	a,					/* Input matrix          (in) */
	b;					/* Transposed matrix    (out) */
{
   mat_mt tmp;

   memcp(tmp, a, sizeof tmp);

   b[0][0] = tmp[0][0];	b[1][1] = tmp[1][1];	b[2][2] = tmp[2][2];
   b[0][1] = tmp[1][0];	b[1][0] = tmp[0][1];
   b[0][2] = tmp[2][0];	b[2][0] = tmp[0][2];
   b[1][2] = tmp[2][1];	b[2][1] = tmp[1][2];
}
/******************************************************************************
 *  Det.  Determinant of a 3 x 3 matrix   				      *
 ******************************************************************************/
double det(a)
mat_mt	a;					/* Matrix		 (in) */
{
   int	i, j, k;				/* Counters		      */
   register double	deter = 0.0;
   
   for(i = 0, j = 1, k = 2; i < 3; i++, j=(j+1)%3, k=(k+1)%3)
      deter += a[0][i] * (a[1][j]*a[2][k] - a[1][k]*a[2][j]);
   return(deter);
}      
/******************************************************************************
 * invert.  Calculate the inverse of a 3x3 matrix.  Adjoint method.           *
 ******************************************************************************/
void invert(a, b)
mat_mt	a,					/* Input matrix          (in) */
	b;					/* Inverse matrix	(out) */
{
   int	i, j, k, l, m, n;			/* Counters		      */
   register real	deter;			/* Reciprocal of determinant  */
	
   if(a == b)
      message(NULLI, NULLP, FATAL, OVRLAP, "invert");
   if((deter = det(a)) == 0.0)
      message(NULLI, NULLP, FATAL, SNGMAT, "invert");
   deter = 1.0 / deter;
   for(i = 0, j = 1, k = 2; i < 3; i++, j=(j+1)%3, k=(k+1)%3)
      for(l = 0, m = 1, n = 2; l < 3; l++, m=(m+1)%3, n=(n+1)%3)
         b[l][i] = deter*(a[j][m]*a[k][n] - a[j][n]*a[k][m]);
}
