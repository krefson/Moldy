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
 * Revision 1.15  92/09/22  14:48:13  keith
 * Tidied up calls to improve "lint" rating.
 * 
 * Revision 1.14  92/06/26  17:03:25  keith
 * Got rid of assumption that memory returned by talloc() or
 * arralloc() is zeroed.  This enhances ANSI compatibility.
 * Removed memory zeroing from alloc.c() in consequence.
 * 
 * Revision 1.13  91/08/15  18:12:12  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.12  91/03/12  15:43:14  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.11  90/10/23  20:13:18  keith
 * Added dummy function call to inhibit vectorization.
 * This allows use of 'ivdep' compiler options and also
 * works round certain bugs in cray's scc compiler.
 * 
 * Revision 1.10  90/05/16  18:40:38  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.9  90/05/01  12:51:19  keith
 * Corrected line-breaking algorithm in print_rdf() to never exceed line length.
 * 
 * Revision 1.8  90/02/02  15:30:19  keith
 * Split inner loop in rdf_calc() into one large, vectorisable loop and
 * one smaller, unvectorisable binning loop. 
 * Replaced library function floor() by (vectorisable) macro.
 * 
 * Revision 1.7  89/12/21  16:30:02  keith
 * Reversed indices in 'site' and 'site_force' to allow stride of 1 in ewald.
 * 
 * Revision 1.6  89/11/20  11:59:19  keith
 * Corrected normalisation in calculation of rdf.  Changed interface (cf main).
 * 
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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/rdf.c,v 1.15 92/09/22 14:48:13 keith Exp $";
#endif
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#if  defined(convexvc) || defined(stellar)
#include <fastmath.h>
#else
#include <math.h>
#endif
#include	<stdio.h>
#include 	"string.h"
/*========================== Program include files ===========================*/
#include	"structs.h"
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
double	det();
void	invert();
void	new_line();
void	put_line();
double	precision();
void	inhibit_vectorization();		/* Self-explanatory dummy     */
/*========================== External data references ========================*/
extern contr_t	control;
/*========================== External data definitions  ======================*/
int	***rdf;				/* The RDF 'array'	      */
/*========================== Macros ==========================================*/
#define MATMUL(i, m, r) (m[i][0]*r[0] + m[i][1]*r[1] + m[i][2]*r[2])
#define floor(x)  (double)((int)((x) + 10) - 10) /* Vectorisable macro floor()*/
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
   (void)memset((char*)rdf_base, 0, control.nbins * max_id * (max_id - 1) / 2*sizeof(int));
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
real		**site;				/* Site co-ordinate array     */
system_p	system;				/* System info struct	      */
spec_t	species[];			/* Species info struct array  */
{
   spec_p	spec;
   register double	t;
   double	r;
   real		*site0 = site[0], *site1 = site[1], *site2 = site[2];
   vec_t	rij;
   double	rbin;				/* 1.0/bin width	      */
   int		imol, isite, jsite, nsites = system->nsites;
   int		*id = ialloc(system->nsites),
                *bind = ialloc(system->nsites);
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
         memcp(id_ptr, spec->site_id, spec->nsites*sizeof(int));
         id_ptr += spec->nsites;
      }

    for(isite = 0; isite < nsites; isite++)
    {
VECTORIZE
       for(jsite = isite+1; jsite < nsites; jsite++)
       {
          rij[0] = site0[jsite] - site0[isite];
          rij[1] = site1[jsite] - site1[isite];
          rij[2] = site2[jsite] - site2[isite];
          
          rij[0] -= lx  *      floor(MATMUL(0,hinv,rij) + 0.5);
          rij[0] -= lxy * (t = floor(MATMUL(1,hinv,rij) + 0.5));
          rij[1] -= ly  * t;
          rij[0] -= lxz * (t = floor(MATMUL(2,hinv,rij) + 0.5));
          rij[1] -= lyz * t;
          rij[2] -= lz  * t;
           
          r = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
	  bind[jsite] = rbin*r;
       }
       for(jsite = isite+1; jsite < nsites; jsite++)
       {
	  inhibit_vectorization();
	  if( bind[jsite] < control.nbins )
             rdf[id[isite]][id[jsite]][bind[jsite]]++;
       }
    }
    xfree(id);    xfree(bind);
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
asite_t		site_info[];
{
   int		idi, idj, col, ibin, is;
   int		*nfrac = ialloc(system->max_id);  /* Per site count of system*/
   spec_t	*spec;
   double	bin = control.limit/control.nbins,
   		bincb = bin*bin*bin,
   		rho = system->nsites/det(system->h);
   double	norm;
   char		buf[32];
   
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
            sprintf(buf, " %7f",rdf[idi][idj][ibin]*norm/SQR(0.5+ibin));
	    col += strlen(buf);
            if(col > control.page_width)
            {
               col = strlen(buf);
	       new_line();
	    }
	    fputs(buf, stdout);
            rdf[idi][idj][ibin] = 0;				/* Reset      */
         }
	 new_line();
      }
   put_line('_');
   xfree(nfrac);
}
