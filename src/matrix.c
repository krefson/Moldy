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
 *       $Log$
 */
#ifndef lint
static char *RCSid = "$Header$";
#endif
/*========================== Program include files ===========================*/
#include "defs.h"
#include "messages.h"
/*========================== External function declarations ==================*/
void	message();		/* Error message and exit handler	      */
/*============================================================================*/
/******************************************************************************
 * 3 x 3  Matrix - vector multiply  (of multiple vectors)                     *
 * The input and output vectors need not necessarily be distinct              *
 ******************************************************************************/
void mat_vec_mul(m, in_vec, out_vec, number)
int		number;		/* Number of vectors to be multiplied         */
mat_t		m;		/* Matrix                                     */
vec_p		in_vec,		/* Input vector, [number][3]          (in/out)*/
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
      for(i = 0; i < number; i++)
      {
         a0 = in_vec[i][0];  a1 = in_vec[i][1];  a2 = in_vec[i][2];

         in_vec[i][0] = m[0][0]*a0 + m[0][1]*a1 + m[0][2]*a2;
         in_vec[i][1] = m[1][0]*a0 + m[1][1]*a1 + m[1][2]*a2;
         in_vec[i][2] = m[2][0]*a0 + m[2][1]*a1 + m[2][2]*a2;
      }
   else if(ABS(in_vec-out_vec) >= number) /* parameters are distinct arrays  */
      for(i = 0; i < number; i++)
      {
         a0 = in_vec[i][0];  a1 = in_vec[i][1];  a2 = in_vec[i][2];

         out_vec[i][0] = m[0][0]*a0 + m[0][1]*a1 + m[0][2]*a2;
         out_vec[i][1] = m[1][0]*a0 + m[1][1]*a1 + m[1][2]*a2;
         out_vec[i][2] = m[2][0]*a0 + m[2][1]*a1 + m[2][2]*a2;
      }
   else			/* in_vec and out_vec overlap, but are not identical  */
      message(NULLI, NULLP, FATAL, OVRLP1, in_vec, out_vec, number);
}
/******************************************************************************
 * mat_mul   Multiply two 3 x 3 matrices. Result can NOT overwrite input.     *
 ******************************************************************************/
void mat_mul(a, b, c)
mat_t	a,					/* Input matrix 1        (in) */
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
mat_t	a, 				/* Input matrix			(in)  */
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
mat_t	a,					/* Input matrix 1        (in) */
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
mat_t	a,					/* Input matrix          (in) */
	b;					/* Transposed matrix    (out) */
{
   register int	i, j;				/* Counters		      */
   register real tmp;
   
   for(i = 0; i < 3; i++)
      for(j = i; j < 3; j++)
      {
         tmp     = a[i][j];
         b[i][j] = a[j][i];
         b[j][i] = tmp;
      }
}
/******************************************************************************
 *  Det.  Determinant of a 3 x 3 matrix   				      *
 ******************************************************************************/
double det(a)
mat_t	a;					/* Matrix		 (in) */
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
mat_t	a,					/* Input matrix          (in) */
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
