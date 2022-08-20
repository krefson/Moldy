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
 * rdf		Functions to accumulate and calculate radial distribution     *
 *		functions. Contents:					      *
 * init_rdf()		Prepare to collect rdf's. Must be called first	      *
 * rdf_calc() (NOT USED)Bin site-site distances and accumulate RDF's	      *
 * rdf_accum()		Bin site-site distances and accumulate RDF's	      *
 * print_rdf()		Calculate RDF from binned data and output it.	      *
 * rdf[idi][idj][ibin]	RDF database (also accessed by 'restart')	      *
 ******************************************************************************
 */
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#ifdef stellar
#include <fastmath.h>
#else
#include <math.h>
#endif
#include 	<string.h>
#include	<stdio.h>
/*========================== Program include files ===========================*/
#include	"structs.h"
/*========================== External function declarations ==================*/
gptr            *talloc(int n, size_mt size, int line, char *file);
				       /* Interface to memory allocator       */
void            tfree(gptr *p);	       /* Free allocated memory	      	      */
double	det(mat_mt);
void	invert(mat_mt, mat_mt);
void	new_line(void);
void	put_line(int c);
double	precision(void);
void	inhibit_vectorization(void);		/* Self-explanatory dummy     */
/*========================== External data references ========================*/
extern contr_mt	control;    		    /* Main simulation control parms. */
/*====================================data definitions  ======================*/
static   float	***rdf;				/* The RDF 'array'	      */
static   float	*rdf_base;			/* base of data area          */
static   int    rdf_size;
gptr *rdf_ptr(int *size)
{
   *size = rdf_size;
   return (gptr*)rdf_base;
}
/*========================== Macros ==========================================*/
#define MATMUL(i, m, r) (m[i][0]*r[0] + m[i][1]*r[1] + m[i][2]*r[2])
#define floor(x)  (double)((int)((x) + 10) - 10) /* Vectorisable macro floor()*/
/*============================================================================*/

/******************************************************************************
 *  rdf_init.  Prepare to bin rdf's.  Allocate memory and pointers	      *
 ******************************************************************************/
void	init_rdf(system_mp system)		/* System info struct	      */
{
   int		max_id = system->max_id;
   int		idi, idj;
   float	*base;

   rdf = aalloc(max_id, float ** );
   rdf_size = control.nbins * max_id * (max_id - 1) / 2;
   base = rdf_base = aalloc(rdf_size,float);
   memst(base, 0, rdf_size*sizeof(float));
   for(idi = 1; idi < max_id; idi++)
      rdf[idi] = aalloc(max_id, float * );
   for(idi = 1; idi < max_id; idi++)
      for(idj = idi; idj < max_id; idj++)
      {
         rdf[idi][idj] = rdf[idj][idi] = base;
         base += control.nbins;
      }
}
/******************************************************************************
 *  rdf_accum.  Calculate site pair distances and bin for RDF.                *
 ******************************************************************************/
void	rdf_accum(double density, int lo, int hi, real *rsq, int iid, int *id, int *nab)
{
   int  j, bin;
   int  *nj = nab + lo;
   float  **rdfi = rdf[iid];
   double rbin = control.nbins / control.limit;
   double	invrho = 1.0/density;

   for(j=lo; j<hi; j++,nj++)
   {
      bin = rbin*sqrt(rsq[j]);
      if( bin < control.nbins)
	 rdfi[id[*nj]][bin]+=invrho;      
   }
}
/******************************************************************************
 * print_rdf.  Calculate the radial distribution function from the binned pair*
 * distances in rdf.  If ib is bin index, B(ib) is # binned, n is # sites of  *
 * each type (or id), rho is number density of all sites, sum is total of all *
 * bins B(ib) for this id and b = bin width then:-			      *
 * 	g(r) = B(ib)/(4 pi r**2 rho b) * (n-1)/sum,  where r(ib) =b*(ib+1/2)  *
 ******************************************************************************/
void	print_rdf(system_mt *system, spec_mt *species, site_mt *site_info)
{
   int		idi, idj, col, ibin, is;
   int		*nfrac = ialloc(system->max_id);  /* Per site count of system*/
   spec_mt	*spec;
   double	bin = control.limit/control.nbins,
                bincb = bin*bin*bin;
   double	norm;
   char		buf[32];
   
   memst(nfrac,0,system->max_id*sizeof(*nfrac));
   for(spec = species; spec < species+system->nspecies; spec++)
      for(is = 0; is < spec->nsites; is++)
	 nfrac[spec->site_id[is]] += spec->nmols;

   put_line('_');
   (void)printf("\tRadial Distribution Functions\tBin width=%g", bin);
   new_line();

   for(idi = 1; idi < system->max_id; idi++)
      for(idj = idi; idj < system->max_id; idj++)
      {
	 if( nfrac[idi] > 0 && nfrac[idj] > 0 ) 
	 {
	    (void)printf("\t%s-%s RDF", site_info[idi].name, site_info[idj].name);
	    new_line();
	    col = 0;

	    norm = 3.0*system->nsites*control.rdf_interval /
	      (4.0*PI*bincb*nfrac[idi]*nfrac[idj]*control.rdf_out);
	    if( idi == idj )
	       norm += norm;
	    for(ibin = 0; ibin < control.nbins; ibin++)
	    {
	       sprintf(buf, " %7f",rdf[idi][idj][ibin]*norm/
		       ( 3*(SQR(ibin) + ibin) + 1));
	       col += strlen(buf);
	       if(col > control.page_width)
	       {
		  col = strlen(buf);
		  new_line();
	       }
	       fputs(buf, stdout);
	       rdf[idi][idj][ibin] = 0.0;				/* Reset      */
	    }
	    new_line();
	 }
      }
   put_line('_');
   xfree(nfrac);
}
