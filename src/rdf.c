/******************************************************************************
 * rdf		Functions to accumulate and calculate radial distribution     *
 *		functions. Contents:					      *
 * init_rdf()		Prepare to collect rdf's. Must be called first	      *
 * rdf_calc()		Bin site-site distances and accumulate RDF's	      *
 * print_rdf()		Calculate RDF from binned data and output it.	      *
 * rdf[idi][idj][ibin]	RDF database (also accessed by 'restart')	      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	rdf.c,v $
 * Revision 1.5  89/10/24  17:18:33  keith
 * Modified pbc algorithm to use floor() library function.
 * Now works with non-orthorhombic cell.
 * 
 * Revision 1.4  89/07/07  10:37:02  keith
 * Added check for uniformly zero RDF which otherwise gave divide by zero.
 * 
 * Revision 1.3  89/06/22  15:45:10  keith
 * Tidied up loops over species to use one pointer as counter.
 * 
 * Revision 1.2  89/06/01  21:25:16  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.1  89/04/20  16:00:53  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/tigger/keith/md/RCS/rdf.c,v 1.5 89/10/24 17:18:33 keith Exp $";
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
void	invert();
void	new_line();
void	put_line();
double	precision();
/*========================== External data references ========================*/
extern contr_t	control;
/*========================== External data definitions  ======================*/
int	***rdf;				/* The RDF 'array'	      */
/*========================== Macros ==========================================*/
#define MATMUL(i, m, r) (m[i][0]*r[0] + m[i][1]*r[1] + m[i][2]*r[2])
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
   register double	r, t;
   vec_t	rij;
   double	rbin;				/* 1.0/bin width	      */
   int		imol, isite, jsite;
   int		*id = ialloc(system->nsites);
   int		*id_ptr;
   mat_t	hinv;
   double	lx   = system->h[0][0], lxy  = system->h[0][1],
		ly   = system->h[1][1], lxz  = system->h[0][2],
		lz   = system->h[2][2], lyz  = system->h[1][2];
   invert(system->h,hinv);

   rbin = control.nbins / control.limit;
   
/*  Construct and fill expanded site-identifier array, id                     */
   id_ptr = id;
   for (spec = species; spec < &species[system->nspecies]; spec++)
      for(imol = 0; imol < spec->nmols; imol++)
      {
         (void)memcpy((char*)id_ptr, (char*)spec->site_id, 
		      (int)spec->nsites*sizeof(int));
         id_ptr += spec->nsites;
      }

    for(isite = 0; isite < system->nsites; isite++)
       for(jsite = isite+1; jsite < system->nsites; jsite++)
       {
          rij[0] = site[jsite][0] - site[isite][0];
          rij[1] = site[jsite][1] - site[isite][1];
          rij[2] = site[jsite][2] - site[isite][2];
          
          rij[0] -= lx  *      floor(MATMUL(0,hinv,rij) + 0.5);
          rij[0] -= lxy * (t = floor(MATMUL(1,hinv,rij) + 0.5));
          rij[1] -= ly  * t;
          rij[0] -= lxz * (t = floor(MATMUL(2,hinv,rij) + 0.5));
          rij[1] -= lyz * t;
          rij[2] -= lz  * t;
           
          r = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
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
void	print_rdf(system, species, site_info)
system_t	*system;
spec_t		species[];
site_t		site_info[];
{
   int		idi, idj, col, ibin, is;
   int		*nfrac = ialloc(system->max_id);  /* Per site count of system*/
   spec_t	*spec;
   double	bin = control.limit/control.nbins,
   		bincb = bin*bin*bin,
   		rho = system->nsites/det(system->h);
   double	norm;
   
   for(spec = species; spec < species+system->nspecies; spec++)
      for(is = 0; is < spec->nsites; is++)
	 nfrac[spec->site_id[is]] += spec->nmols;

   put_line('_');
   (void)printf("\tRadial Distribution Functions\tBin width=%g", bin);
   new_line();

   for(idi = 1; idi < system->max_id; idi++)
      for(idj = idi; idj < system->max_id; idj++)
      {
         (void)printf("\t%s-%s RDF", site_info[idi].name, site_info[idj].name);
         new_line();
	 col = 0;

	 norm = (system->nsites-1)*control.rdf_interval /
	    (4.0*PI*bincb*rho*nfrac[idi]*nfrac[idj]*control.rdf_out);
	 if( idi == idj )
	    norm += norm;
         for(ibin = 0; ibin < control.nbins; ibin++)
         {
            col += printf(" %7f",rdf[idi][idj][ibin]*norm/SQR(0.5+ibin));
            if(col > control.page_width - 8 || ibin == control.nbins - 1) 
            {
               col = 0;
	       new_line();
	    }
            rdf[idi][idj][ibin] = 0;				/* Reset      */
         }
      }
   put_line('_');
   cfree((char*)nfrac);
}
