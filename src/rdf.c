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
 *      Revision Log
 *       $Log: rdf.c,v $
 *       Revision 2.13  2002/09/19 09:26:30  kr
 *       Tidied up header declarations.
 *       Changed old includes of string,stdlib,stddef and time to <> form
 *
 *       Revision 2.12  2000/12/06 17:45:33  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.11  2000/10/20 15:15:48  keith
 *       Incorporated all mods and bugfixes from Beeman branch up to Rel. 2.16
 *
 *       Revision 2.9.2.1  2000/09/01 11:23:39  keith
 *       Fixed to cope with site-id's with no corresponding sites, ie gaps
 *       in list.  Used to generate NaN or divide-by-zero.
 *
 *       Revision 2.9  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 *
 *       Revision 2.8  1996/01/15 15:19:12  keith
 *       New function rdf_accum for parallel accululation or RDF data.
 *       rdf_ptr() now takes one param, *size, which returns # items.
 *
 *       Revision 2.7  1994/06/08 13:22:31  keith
 *       Null update for version compatibility
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Got rid of all global (external) data items except for
 * "control" struct and constant data objects.  The latter
 * (pot_dim, potspec, prog_unit) are declared with const
 * qualifier macro which evaluates to "const" or nil
 * depending on ANSI/K&R environment.
 * Also moved as many "write" instantiations of "control"
 * members as possible to "startup", "main" leaving just
 * "dump".
 *
 * Declared as "static"  all functions which should be.
 *
 * Added const qualifier to (re-)declarations of ANSI library
 * emulation routines to give reliable compilation even
 * without ANSI_LIBS macro. (#define's away for K&R
 * compilers)
 *
 * Revision 2.5  94/01/18  13:32:55  keith
 * Null update for XDR portability release
 * 
 * Revision 2.3  93/10/28  10:28:08  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.1  93/07/29  17:37:11  keith
 * Changed evaluation of rdf to be in exact accordance with
 * manual.  That is, replaced (nsites-1) with nsites and
 * 3(b+1/2)^2 with [(b+1)^3 - b^3].  Should make negligible
 * difference to results. 
 * 
 * Revision 2.0  93/03/15  14:49:19  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.19  93/03/12  12:25:17  keith
 * Got rid of unneccesary convex special case.
 * 
 * Revision 1.18  93/03/09  15:59:09  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.17  92/11/02  17:31:17  keith
 * Fixed bug where counter array "nfrac" was not initialised.
 * 
 * Revision 1.16  92/10/28  14:09:52  keith
 * Changed "site_[tp]" typedefs to avoid name clash on HP.
 * 
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
static char *RCSid = "$Header: /usr/users/moldydv/CVS/moldy/src/rdf.c,v 2.13 2002/09/19 09:26:30 kr Exp $";
#endif
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
