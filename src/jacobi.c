/******************************************************************************
 * Jacobi	Find the eigenvalues and eigenvectors of a 3x3 matrix.	      *
 *		Modified from Press W.H., Flannery B.P. Teukolsky S.A. &      *
 *		Vetterling W.T. "Numerical Recipies in C" CUP, 1988 p364.     *
 ******************************************************************************
 *      Revision Log
 *       $Log$
 */
#ifndef lint
static char *RCSid = "$Header$";
#endif
/*========================== Library include files ===========================*/
#include <math.h>
/*========================== Program include files ===========================*/
#include	"messages.h"
#include	"defs.h"
/*========================== External function declarations ==================*/
void		message();
/*========================== Macros ==========================================*/
#define	NP	3
#define ROTATE(a, i, j, k, l) {g=a[i][j]; h=a[k][l]; a[i][j] = g - s*(h+g*tau);\
			       a[k][l] = h + s*(g-h*tau);}

void	jacobi(a, n, d, v, nrot)
real	a[][NP],
	d[],
	v[][NP];
int	n,
	*nrot;
{
      real	b[NP], z[NP];
      int	ip, iq, i,j;
      real 	sm, thresh, g, h, t, c, s, tau, theta;

      for(ip=0; ip < n; ip++)		/* Initialise V to identity matrix    */
        for(iq=0; iq < n; iq++)
          v[ip][iq]= (ip==iq) ? 1.0 : 0.0;

      for(ip=0; ip < n; ip++)		/* Initialise b and d to diagonal of a*/
      {
        b[ip]=a[ip][ip];
        d[ip]=b[ip];
        z[ip]=0.0;
      }
      *nrot=0;
      for(i=1; i < 50; i++)
      {
        sm=0.0;
        for(ip=0; ip < n-1; ip++)
          for(iq=ip+1; iq < n; iq++)
            sm += fabs(a[ip][iq]);

        if(sm == 0.0) return;

	if (i < 4)
	   thresh = 0.2*sm/(n * n);
	else
	   thresh = 0.0;

        for(ip=0; ip < n-1; ip++)
        {
          for(iq=ip+1; iq < n; iq++)
          {
            g = 100.0 * fabs(a[ip][iq]);
            if((i > 4) && (fabs(d[ip]) + g == fabs(d[ip]))
                       && (fabs(d[iq]) + g == fabs(d[iq])))
              a[ip][iq] = 0.0;
            else if(fabs(a[ip][iq]) > thresh)
            {
              h=d[iq]-d[ip];
              if(fabs(h) + g == fabs(h))
                t=a[ip][iq]/h;
              else
              {
                theta=0.5*h/a[ip][iq];
                t=1./(fabs(theta)+sqrt(1.+theta*theta));
                t = (theta < 0.) ? -t : t;
	      }
              c=1./sqrt(1.0+t*t);
              s=t*c;
              tau=s/(1.+c);
              h=t*a[ip][iq];
              z[ip]=z[ip]-h;
              z[iq]=z[iq]+h;
              d[ip]=d[ip]-h;
              d[iq]=d[iq]+h;
              a[ip][iq]=0.0;
              for(j=0; j <= ip-1; j++)
	         ROTATE(a, j, ip, j, iq)
              for(j=ip+1; j <= iq-1; j++)
	         ROTATE(a, ip, j, j, iq)
              for(j=iq+1; j < n; j++)
	         ROTATE(a, ip, j, iq, j)
              for(j=0; j < n; j++)
	         ROTATE(v, j, ip, j, iq)
              (*nrot)++;
            }
          }
        }
        for(ip=0; ip < n; ip++)
        {
          b[ip]=b[ip]+z[ip];
          d[ip]=b[ip];
          z[ip]=0.;
        }
      }
      message(NULLI, NULLP, FATAL, "Jacobi took 50 iterations");
}
