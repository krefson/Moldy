/******************************************************************************
 * rdf		Functions to accumulate and calculate radial distribution     *
 *		functions. Contents:					      *
 * init_rdf()		Prepare to collect rdf's. Must be called first	      *
 * rdf_calc()		Bin site-site distances and accumulate RDF's	      *
 * print_rdf()		Calculate RDF from binned data and output it.	      *
 * rdf[idi][idj][ibin]	RDF database (also accessed by 'restart')	      *
 ******************************************************************************
 *      Revision Log
 *       $Log$
 */
#ifndef lint
static char *RCSid = "$Header$";
#endif
/*========================== Library include files ===========================*/
#include	<math.h>
#include 	"string.h"
/*========================== Program include files ===========================*/
#include	"structs.h"
/*========================== Library declarations ============================*/
void    cfree();
/*========================== External function declarations ==================*/
double	det();
void	new_line();
void	put_line();
double	precision();
/*========================== External data references ========================*/
extern contr_t	control;
/*========================== External data definitions  ======================*/
int	***rdf;				/* The RDF 'array'	      */
/*============================================================================*/

/******************************************************************************
 *  rdf_init.  Prepare to bin rdf's.  Allocate memory and pointers	      *
 ******************************************************************************/
void	init_rdf(system)
system_p	system;				/* System info struct	      */
{
   int		*rdf_base;			/* base of data area          */
   int		max_id = system->max_id;
   int		idi, idj;

   rdf = aalloc(max_id, int ** );
   rdf_base = ialloc(control.nbins * max_id * (max_id - 1) / 2);
   for(idi = 1; idi < max_id; idi++)
      rdf[idi] = aalloc(max_id, int * );
   for(idi = 1; idi < max_id; idi++)
      for(idj = idi; idj < max_id; idj++)
      {
         rdf[idi][idj] = rdf[idj][idi] = rdf_base;
         rdf_base += control.nbins;
      }
   if(control.limit <= 0.0)				/* Choose a limit     */
      control.limit = 0.5*MIN3(system->h[0][0],system->h[1][1],system->h[2][2]);
}
/******************************************************************************
 *  rdf_calc.  Calculate site pair distances and bin for RDF.                 *
 ******************************************************************************/
void	rdf_calc(site, system, species)
vec_t		site[];				/* Site co-ordinate array     */
system_p	system;				/* System info struct	      */
spec_t	species[];			/* Species info struct array  */
{
   spec_p	spec;
   register double	r, rx, ry, rz;
   double	rbin;				/* 1.0/bin width	      */
   int		ispec, imol, isite, jsite;
   int		*id = ialloc(system->nsites);
   int		*id_ptr;
   double	lx   = system->h[0][0],	/* Temporaries for unit cell vectors  */
		lx_r = 2.0/lx,		/* and their reciprocals; for use     */
		ly   = system->h[1][1],	/* in applying the periodic boundary  */
		ly_r = 2.0/ly,		/* conditions.			      */
		lxy  = system->h[0][1],
		lz   = system->h[2][2],
		lz_r = 2.0/lz,
		lxz  = system->h[0][2],
		lyz  = system->h[1][2];

   rbin = control.nbins / control.limit;
   
/*  Construct and fill expanded site-identifier array, id                     */
   id_ptr = id;
   for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
      for(imol = 0; imol < spec->nmols; imol++)
      {
         (void)memcpy((char*)id_ptr, (char*)spec->site_id, 
		      (int)spec->nsites*sizeof(int));
         id_ptr += spec->nsites;
      }

    for(isite = 0; isite < system->nsites; isite++)
       for(jsite = isite+1; jsite < system->nsites; jsite++)
       {
          rx = site[jsite][0] - site[isite][0];
          ry = site[jsite][1] - site[isite][1];
          rz = site[jsite][2] - site[isite][2];
          
          rx -= lx  * (int)(rx * lx_r);		/* Apply pbc's    	      */
          rx -= lxy * (int)(ry * ly_r);		/* than usual because	      */
          ry -= ly  * (int)(ry * ly_r);		/* - More complicated 	      */
          rx -= lxz * (int)(rz * lz_r);		/* cell.  Note recip's	      */
          ry -= lyz * (int)(rz * lz_r);		/* to avoid division. 	      */
          rz -= lz  * (int)(rz * lz_r);		/* of non-cubic unit	      */
           
          r = sqrt(rx*rx + ry*ry + rz*rz);
	  if ( r < control.limit )
             rdf[id[isite]][id[jsite]][(int)(rbin * r)]++;
       }
    cfree((char*)id);
}
/******************************************************************************
 * print_rdf.  Calculate the radial distribution function from the binned pair*
 * distances in rdf.  If ib is bin index, B(ib) is # binned, n is # sites of  *
 * each type (or id), rho is number density of all sites, sum is total of all *
 * bins B(ib) for this id and b = bin width then:-			      *
 * 	g(r) = B(ib)/(4 pi r**2 rho b) * (n-1)/sum,  where r(ib) =b*(ib+1/2)  *
 ******************************************************************************/
void	print_rdf(system, site_info)
system_p	system;
site_t		site_info[];
{
   int		idi, idj, col, sum, ibin;
   double	bin = control.limit/control.nbins,
   		bincb = bin*bin*bin,
   		rho = system->nsites/det(system->h);
   double	norm;
   FILE		*f = control.out;
   
   put_line('_');
   (void)fprintf(f,"\tRadial Distribution Functions\tBin width=%g", bin);
   new_line();

   for(idi = 1; idi < system->max_id; idi++)
      for(idj = idi; idj < system->max_id; idj++)
      {
         (void)fprintf(f,"\t%s-%s RDF",
		       site_info[idi].name, site_info[idj].name);
         new_line();
         sum = 0; col = 0;
         for(ibin = 0; ibin < control.nbins; ibin++)
            sum += rdf[idi][idj][ibin];
         norm = (system->nsites - 1) / (4.0 * PI * rho * bincb * sum);
         for(ibin = 0; ibin < control.nbins; ibin++)
         {
            col += fprintf(f," %7f",rdf[idi][idj][ibin]*norm/SQR(0.5+ibin));
            if(col > control.page_width - 8 || ibin == control.nbins - 1) 
            {
               col = 0;
	       new_line();
	    }
            rdf[idi][idj][ibin] = 0;				/* Reset      */
         }
      }
      put_line('_');
}
