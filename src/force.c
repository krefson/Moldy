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
 * Force        This module contains functions to implement the 'link cell'   *
 *              interatomic force calculation (Hockney, R.W. & Eastwood J.W.  *
 *              "Computer Simulation Using Particles" McGraw-Hill (1981), 277)*
 *              It is vectorised and optimised for a CRAY XMP and Convex C1,  *
 *              but should run well on any vector machine with reasonably     *
 *              efficient scatter/gather. (And of course on scalar machines!) *
 *              The actual calculation of the potential is in a different     *
 *              module (kernel.c) for ease of modification.                   *
 ******************************************************************************
 *       $Log: force.c,v $
 *       Revision 2.27.2.2  2002/06/20 17:50:59  kr
 *       Patrick's mods to do 1/r**7 terms with an Ewald sum
 *       (very slightly tidied up).
 *
 *       Revision 2.27.2.1  2002/03/13 10:27:52  kr
 *       Trial version incorporating reciprocal-space summation for r^-2 and r^-6
 *       interactions.  This version implements a new potential "genpot46" to activate.
 *
 *       Revision 2.27  2002/03/04 16:08:12  kr
 *       Fixed a number of bugs in dumpext and dumpconv related to reading the
 *       sysinfo section of the dump files.
 *
 *       Revision 2.26  2001/02/22 10:30:16  keith
 *       Reinstated capability of "molecular" cutoffs.
 *
 *       Revision 2.25  2001/02/19 12:22:25  keith
 *       Fixed serious bug in allocation of "pbclookup" array which will
 *       cause heap corruption and SEGVs on constant-pressure runs.
 *
 *       Revision 2.24  2001/02/05 11:47:34  keith
 *       Improved test of close site-site approaches to catch all intermolecular
 *       ones.  The avoidance test for intramolecular approaches invalidated the
 *       previous approach of only catching the first such close contact in a list.
 *
 *       Revision 2.23  2000/12/08 12:28:21  keith
 *       Rewrote misleading intramolecular energy message
 *       Corrected (silent) bug in arglist of kernel in poteval()
 *
 *       Revision 2.22  2000/12/06 17:45:29  keith
 *       Tidied up all ANSI function prototypes.
 *       Added LINT comments and minor changes to reduce noise from lint.
 *       Removed some unneccessary inclusion of header files.
 *       Removed some old and unused functions.
 *       Fixed bug whereby mdshak.c assumed old call for make_sites().
 *
 *       Revision 2.21  2000/10/20 15:15:47  keith
 *       Incorporated all mods and bugfixes from Beeman branch up to Rel. 2.16
 *
 *       Revision 2.20.2.2  2000/10/20 13:59:32  keith
 *       Incorporated new neightbour list stuff into accel.c.
 *       Removed old "poteval" from accel.c.  Now use one in
 *       force.
 *       Other errors corrected and declarations tidied somewhat.
 *
 *       Revision 2.20.2.1  2000/10/20 11:48:52  keith
 *       Incorporated new neighbour list indexing algorithm from 2.16
 *
 *       Revision 2.19.2.1.2.7  2000/10/19 16:22:12  keith
 *       Fixed bug whereby intramolecular PE was subtracted on *every*
 *       thread in a parallel run.
 *
 *       Revision 2.19.2.1.2.6  2000/10/17 12:47:53  keith
 *       Rebased reloc_v array to index from 0 to NIMCELLS-1.  This
 *       is ANSI-C compliant and avoids a pointer reference.
 *
 *       Extracted code from force_inner()'s innermost loop into 3 new functions
 *       mk_r_sqr(), mk_forces() and scatter_forces().  Optimised these with local
 *       variables and using "loads-before-stores" trick for optimization.
 *
 *       Revision 2.19.2.1.2.5  2000/10/16 09:31:32  keith
 *       Implemented new scheme for encoding and handling relative neighbour lists.
 *       The neighbour cell list is encoded as an index into an extended 3D
 *       (N*nx x N*ny x N*nz) array with the usual 1-D index.  Absolute cell indices
 *       are then computed by simple integer addition. Finally the
 *       transformation back to a (nx x ny x nz) cell index and finding the
 *       periodic MD cell image translation are handled by a lookup table.
 *
 *       This gives a 2-3 fold speedup of site-neighrbour list.  Even the 108
 *       argon atom test with 10x10x10 partitioning isn't completely dominated
 *       by bookkeeping any more.
 *
 *       Revision 2.19.2.1.2.4  2000/10/13 13:57:41  keith
 *       Some tidying up of code and comments.
 *
 *       Revision 2.19.2.1.2.3  2000/10/13 12:59:09  keith
 *       Hybrid approach.  My neighbour relocation combined with
 *       H. Bekker's virial method.  Seems to be marginally slower
 *       than standard version.
 *
 *       Revision 2.19.2.1.2.2  2000/10/12 14:51:28  keith
 *       Updated rdf_inner to use new calling sequence.
 *       Moved generation of reloc_v tables to force_calc().
 *       Now workd for framework simulations.
 *
 *       Revision 2.19.2.1.2.2  2000/10/11 17:41:23  keith
 *       Parameterized indexing of MD cell image translation table (KMIN and
 *       KMAX).  Avoided indices with ii<0 since they would never be used.
 *       There are no cells with -ve x in the neighbour cell list.
 *
 *       Added molecule index array molout[nsites].  This is used to avoid
 *       flagging close approaches between sites belonging to the same
 *       molecule.
 *
 *       Revision 2.19.2.1.2.1  2000/10/11 16:11:11  keith
 *       First working version of H. Bekker's pbc algorithm.  This computes
 *       forces and stresses correctly without computing the virial in the
 *       inner loop.
 *
 *       It relies on atomic sites being assigned to cells rather than
 *       molecules, and should therefore be more efficient for systems
 *       containing "large" molecules.  This is because the neighbour
 *       list can be smaller.
 *
 *       It gives exactly the same energies, forces and stresses as the standard
 *       version for systems like controp.tips2 and control.quartz, but only in
 *       strict-cutoff mode.  Lazy cutoff mode generates slightly different numbers.
 *
 *       Revision 2.19.2.1  2000/10/11 09:15:40  keith
 *       Experimental testbed version which computes all interaction including
 *       intramolecular ones.
 *
 *       Revision 2.19  1998/12/07 14:44:29  keith
 *       Inlined and optimized spxpy().
 *
 *       Revision 2.18  1998/07/17 14:34:06  keith
 *       Attempt at better algorithm to guess size for neighbour list arrays,
 *       "n_nab_sites".  This determines the subcell with the highest density
 *       of sites and uses this rather than the average in the estimate.
 *       This should work for very inhomogeneous systems.
 *
 *       Revision 2.17  1998/05/07 17:06:11  keith
 *       Reworked all conditional compliation macros to be
 *       feature-specific rather than OS specific.
 *       This is for use with GNU autoconf.
 *
 *       Revision 2.16  1998/01/27 15:47:15  keith
 *       tidied up macro __STDC__ to include ANSI for consistency.
 *
 *       Revision 2.15  1997/11/26 10:22:27  keith
 *       Reordered declarations so that local structs come before
 *       function declarations. Otherwise "protoize" broke!
 *
 *       Revision 2.14  1996/11/04 17:34:30  keith
 *       Moderate rewriting and code re-organization.
 *       1. Simplified PBC relocation calculation, got rid of large
 *          arrays reloc[] etc, saving 32 MB for 8192 waters on T3D.
 *          There is now NO LIMIT to cutoff or RDF limit and macro
 *          parameter NSH is removed.
 *       2. Rewrote site_neighbour_list to be more transparent and got
 *          rid of silly vector sort calls. There are now separate versions for
 *          scalar and vector machines.  Code is now optimized for usual non-
 *          framework case, but it's faster than 2.10 even for frameworks.
 *       3. Corrected misleading comments, reorganized code in force_calc()
 *          to be more transparent.  Commented local variables MUCH better.
 *       4. It's also a bit faster.
 *
 *       Revision 2.13  1996/08/14 16:46:04  keith
 *       Workaround for T3D Cray compiler bug real*int/int ==> int division.
 *       Got rid of unnccessary par_abort() calls - rplaced with exit().
 *       (message on thread 0 calls par_abort()).
 *
 *       Revision 2.12  1996/05/03 16:14:20  keith
 *       Fixed bug whereby reloc could overflow in strict cutoff mode by
 *       tightening up test condition. Also fixed the calculation of
 *       n_nab_sites to reflect the increased cutoff in strict mode.
 *
 *       Revision 2.11  1996/01/17 17:12:47  keith
 *       Incorporated rdf accumulation into forces and parallelized.
 *       New functions rdf_inner(), calls rdf_accum() from rdf.c
 *
 * Revision 1.8.1.8  89/11/01  17:34:15  keith
 * Modified to use SPAXPY vectorised scattered add.
 * 
 * Revision 1.8.1.6  89/10/12  16:28:39  keith
 * Added conditional code to produce histogram of interaction distances
 * and calculate 'minimum image' energy.
 * Added preprocessor constant NSH to fix size of relocation arrays.
 * Fixed mistake in metric G in neighbour_list() (G=h'h not hh').
 * 
 * Revision 1.8.1.5  89/10/02  17:12:59  keith
 * New version of neighbour_list() which works for arbitrary cutoff radii.
 * site_neighbour_list checks for overflow.
 * Main loop limits modified, (in conjunction with mods in site_neighbour_list)
 * to correctly include all molecule-framework interactions.
 * 
 * Revision 1.8.1.4  89/09/12  16:14:41  keith
 * Fixed bug in fill_cells() which didn't increment spec properly in imol loop
 * 
 * Revision 1.8.1.3  89/08/31  11:58:04  keith
 * Fixed bug in 'BIN' macro to correctly handle case of rc<0.
 * 
 * Revision 1.8.1.2  89/08/30  17:00:33  keith
 * Fixed memory overlap bug in site_neighbour_list
 * 
 * Revision 1.8.1.1  89/08/25  15:24:43  keith
 * Mods to add framework structures to simulation model
 * 
 * Revision 1.7  89/08/22  14:48:39  keith
 * Created new variable 'n_nab_sites' for max size of vector arrays.
 * 
 * Revision 1.6  89/07/04  18:43:14  keith
 * Fixed error in kernel and force which led to sites being allocated the
 * wrong potential parameters.  Needed extra parameter to kernel.
 * 
 * Revision 1.5  89/06/22  15:44:23  keith
 * Tidied up loops over species to use one pointer as counter.
 * 
 * Revision 1.4  89/06/14  14:18:49  keith
 * Fixed #ifdef s for CRAY to handle case of UNICOS
 * Fix mistake in VCALLS conditional code.
 * 
 * Revision 1.3  89/06/01  18:01:46  keith
 * Moved `vadd()' from aux.c to force.c for ease of vectorisation.
 * Now no need to compile aux.c with vectorisation.
 * 
 * Revision 1.2  89/05/17  13:53:49  keith
 * Reorganised neighbour list construction in preparation for framework.
 * (Also goes slighty faster)
 * 
 * Revision 1.1  89/04/20  16:00:40  keith
 * Initial revision
 * 
 * Revision 1.3  90/03/29  15:44:51  keith
 * Merged force.c revisions 1.8.1.11-1.8.1.13
 * 
 * Revision 1.2  90/03/09  17:30:29  keith
 * Modified FKERNEL ifdefs for UNICOS.
 * 
 * Revision 1.1  90/01/31  13:19:28  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /usr/users/kr/CVS/moldy/src/force.c,v 2.27.2.2 2002/06/20 17:50:59 kr Exp $";
#endif
/*========================== Program include files ===========================*/
#include        "defs.h"
/*========================== Library include files ===========================*/
#include        <math.h>
#include        <stddef.h>
#include        <string.h>
/*========================== Program include files ===========================*/
#include        "structs.h"
#include        "messages.h"
/*========================== Structs local to module =========================*/
typedef struct cell_s                   /* Prototype element of linked list of*/
{                                       /* molecules within interaction range */
   int          isite, num, frame_type;
   struct cell_s *next;
}               cell_mt;

typedef struct                          /* Prototype of neighbour cell list   */
{                                       /* element.                           */
   real         x, y, z;
}               rvec_mt;
/*========================== External function declarations ==================*/
gptr            *talloc(int n, size_mt size, int line, char *file);
                                       /* Interface to memory allocator       */
void            tfree(gptr *p);        /* Free allocated memory               */
void            afree(gptr *p);        /* Free allocated array                */
int             search_lt(int n, real *x, int ix, double s)
                           ;            /* Search a vector for el. < scalar   */
void            gather(int n, real *a, real *b, int *ix, int lim);  
                                        /* Interface to CRAY gather routine   */
void            mat_mul(real (*a)[3], real (*b)[3], real (*c)[3]);
                                        /* Matrix multiplier                  */
double          det(real (*a)[3]);      /* Determinant of 3x3 matrix          */
void            invert(real (*a)[3], real (*b)[3]); /* 3x3 matrix inverter    */
void            mat_vec_mul(real (*m)[3], vec_mp in_vec, vec_mp out_vec, 
			    int number); /* Matrix by vector multiplier       */
void            transpose(real (*a)[3], real (*b)[3]); /* Do matrix transpose */
void            zero_real(real *r, int n);             /* Initialiser         */
double          precision(void);        /* Floating pt precision.             */
void            kernel();               /* Force kernel routine               */
double          mol_radius(spec_mt *species, int nspecies);  
                                        /* Radius of largest molecule.        */
void            rdf_accum();            /* Bin distances for rdf evaluation.  */
double          poteval(real *potpar, double r, int ptype, double chgsq);
gptr            *arralloc(size_mt,int,...); /* Array allocator                */
void            note(char *, ...);      /* Write a message to the output file */
void            message(int *, ...);    /* Write a warning or error message   */
/*========================== External data references ========================*/
extern  contr_mt control;                   /* Main simulation control parms. */
extern int              ithread, nthreads;
/*========================== Global variables ================================*/
/*========================== Macros ==========================================*/
/*
 * Multiplication factor for size of neighbour list arrays.  If you need
 * to increase this from 1, your system must be *highly* inhomogeneous
 * and may not make sense!
 */
#define         NMULT 2.0
#define         TOO_CLOSE       0.25    /* Error signalled if r**2 < this     */
#define		IMCELL_XTRA 1
#define		IMCELL_L (2*IMCELL_XTRA+1)
#define         NIMCELLS (IMCELL_L*IMCELL_L*IMCELL_L)
#define         NCELL(ix,iy,iz) ((iz)+(nz)*((iy)+(ny)*(ix)))
#define		NCELL_IDX(ix,iy,iz) ((iz)+IMCELL_L*(nz)*((iy)+IMCELL_L*(ny)*(ix)))
#define         LOCATE(r,eps)   NCELL(cellbin(r[0], nx, fnx, eps), \
                                      cellbin(r[1], ny, fny, eps), \
                                      cellbin(r[2], nz, fnz, eps))
#define         BIGINT 32768
#define         IFLOOR(i,n)     ((i+BIGINT*n)/n-BIGINT)
#define moda(hmat) sqrt(SQR(hmat[0][0]) + SQR(hmat[1][0]) + SQR(hmat[2][0]))
#define modb(hmat) sqrt(SQR(hmat[0][1]) + SQR(hmat[1][1]) + SQR(hmat[2][1]))
#define modc(hmat) sqrt(SQR(hmat[0][2]) + SQR(hmat[1][2]) + SQR(hmat[2][2]))
#ifdef __BOUNDS_CHECKING_ON /* The declaration of "array" potp has a lower */
#   define P0 0             /* bound of 1.  This is a violation of strict  */
#else                       /* ANSI, so turn it off when bounds checking.  */
#   define P0 1
#endif
/*============================================================================*/

/******************************************************************************
 *  cellbin.    Safe binning function for putting molecules/sites into cells. *
 *  Any error at the boundaries is disasterous and hard to detect.            *
 *  Results may depend on machine-dependant rounding etc.                     *
 ******************************************************************************/


extern int printf (const char *, ...);

int cellbin(double rc, double fnc, int nc, double eps)
{
   int ibin;

   if(rc < -0.5+eps || rc >= 0.5-eps)
   {
      if(rc < -0.5+eps && rc >= -0.5-eps)
         rc = -0.5;
      else if(rc >= 0.5-eps && rc <= 0.5+eps)
         rc = 0.5-eps;
      else 
         message(NULLI, NULLP, ERROR, 
                 "Co-ordinate out of range in BIN (fill_cells) %.17g\n",rc);
   }
   if( (ibin = floor((rc+0.5)*fnc)) >= nc || ibin < 0)
         message(NULLI, NULLP, ERROR, 
                 "Rounding problem in BIN (fill_cells) %.17g\n",rc);
   return ibin;
}

/******************************************************************************
 ******************************************************************************/
#ifdef DEBUG2
#define NBINS 200
static int bins[NBINS];

hist(int jmin, int jmax, real *rr)
{
   double rbin = 10.0;
   int j;

   for(j = jmin; j < jmax; j++)
      if(rr[j] < SQR(NBINS/rbin))
         bins[(int)(rbin*sqrt(rr[j]))]++;
}
void histout(void)
{
   int j;
   for(j = 0; j < NBINS; j++)
   {
      printf("%d%c",bins[j],(j+1)%10?' ':'\n');
      bins[j] = 0;
   }
}
#endif
/******************************************************************************
 *  Neighbour_list.  Build the list of cells within cutoff radius of cell 0   *
 ******************************************************************************/
static int   *neighbour_list(int *nnabor, mat_mt h, double cutoff,
			     int nx, int ny, int nz, int icheck)
{
   double               dist;
   int                  i, j, ix, iy, iz, mx, my, mz, inabor = 0, nnab;
   static int           onabor=0;
   int                  *nabor;
   vec_mt               s;
   mat_mt               G, htr, htrinv;

   transpose(h, htr);
   mat_mul(htr, h, G);
   invert(htr, htrinv);

   mx = ceil(cutoff*nx*moda(htrinv));
   my = ceil(cutoff*ny*modb(htrinv));
   mz = ceil(cutoff*nz*modc(htrinv));

   nnab = 4*mx*my*mz;
   nabor = ialloc(nnab);
#ifdef DEBUG1
   printf("  Distance    ix    iy    iz      sx        sy          sz\n");
#endif
   for(ix = 0; ix < mx; ix++)
      for(iy = (ix == 0 ? 0 : -my); iy < my; iy++)
         for(iz = (ix == 0 && iy == 0 ? 0 : -mz); iz < mz; iz++)
         {
            s[0] = (double)ix/nx;
            s[1] = (double)iy/ny;
            s[2] = (double)iz/nz;
            dist = 0.0;
            for(i = 0; i < 3; i++)
               for(j = 0; j < 3; j++)
                  dist += s[i]*G[i][j]*s[j];
            if(dist < SQR(cutoff))
            {
	       if( ix > IMCELL_XTRA*nx || ix < -IMCELL_XTRA*nx ||
		   iy > IMCELL_XTRA*ny || iy < -IMCELL_XTRA*ny ||
		   iz > IMCELL_XTRA*nz || iz < -IMCELL_XTRA*nz )
                  message(NULLI, NULLP, FATAL, CUTOFF, IMCELL_XTRA);
               if( inabor >= nnab )
                  message(NULLI, NULLP, FATAL,
                          "Internal error in neighbour_list()");
               nabor[inabor] = NCELL_IDX(ix,iy,iz);
               inabor++;
#ifdef DEBUG1
               printf("%12f %4d %4d %4d %12f %12f %12f\n",
                      dist,ix,iy,iz,s[0],s[1],s[2]);
#endif
            }
         }
   if( icheck )
   {
      if( inabor != onabor )
         note(NABORS,2 * inabor);
      onabor = inabor;
   }
   *nnabor = inabor;
   return(nabor);
}
/******************************************************************************
 *  Strict_Neighbour_list.  Build the list of cells within cutoff radius      *
 *  This is the strict version and includes every cell which has an interior  *
 *  point at a distance less than the cutoff from any interior point of the   *
 *  reference cell.  In fact the distance criterion is cutoff+2*(maximum mol- *
 *  ecular radius).  This ensures that all *sites* which might be closer to-  *
 *  gether than the cutoff are included.                                      *
 *     The method used is based on the fact that the closest interior points  *
 *  of a pair of parallelopiped cells are either at corners of both cells or  *
 *  at the ends of a line perpendicular to the faces of both parallelopipeds. *
 *  The face-face distance is always the shortest, if the perpendicular       *
 *  projection of the faces onto a common plane intersect with each other     *
 *  The goal is therefore to build a list containing all cells which have a   *
 *  corner-corner or face-face distance to the reference cell which is less   *
 *  than the cutoff.                                                          *
 *  Cells may be admitted to the list by multiple corner-corner or face-face  *
 *  contact criteria but must only be recorded in the final list once.  The   *
 *  easiest way to do this is to use a "map" of all the cells potentially     *
 *  within the cutoff radius and to flag occupancy.                           *
 *  It is easier to loop over all grid vectors within the cutoff and assign   *
 *  cells which have that vector as some corner-corner vector with the        *
 *  reference cell, rather than to loop over cells and calculate all corner   *
 *  pair distances.  This method calculates each distance only once instead   *
 *  of 27 times (the number of distinct corner-corner vectors between 2       *
 *  cells).                                                                   *
 *  To exploit Newton's third law the list should contain only the positive   *
 *  hemisphere (in the x direction).                                          *
 *                                                                            *
 *  The algorithm is as follows.                                              *
 *  1) Set up an empty "map"                                                  *
 *  2) Loop over all points on a grid with points at the link-cell corners    * 
 *     choose only points which are closer to the origin than the cutoff.     *
 *     Set the occupancy flag for all cells which have that as a corner-pair  *
 *     vector to the reference cell.   This is a 3x3x3 block of cells centred *
 *     on the cell whose index is the same as the gridpoint being considered. *
 *     Because we only want the "positive x" cells 2x3x3 will suffice.        *
 *  3) Add cells which have a perpendicular face-face separation within the   *
 *     cutoff.  Only the outermost cells need be considered since inner ones  *
 *     are already admitted by corner-pair distance.  Thus                    *
 *     3a) project the facing corner points of the reference cell onto the    *
 *         plane just within the cutoff.                                      *
 *     3b) add the cells with faces which overlap the projection.  Zero, two  *
 *         or four cells are added depending on whether the projected points  *
 *         coincide with the corner points, the edges or none of the faces.   *
 *  4) The final list is built by scanning the map.                           *
 ******************************************************************************/
static int   *strict_neighbour_list(int *nnabor, mat_mt h, double cutoff, 
				    int nx, int ny, int nz, int icheck)
{
   double               dist;
   int                  i, j, k, ix, iy, iz, mx, my, mz, inabor = 0, nnab;
   static int           onabor=0;
   int                  ***cellmap;
   int                  *nabor;
   vec_mt               s;
   mat_mt               G, htr, htrinv;
   int                  face_cells[4][3],mxyz[3], nxyz[3], ixyz, jxyz, kxyz;
   double               proj[3], modabc;

   transpose(h, htr);
   mat_mul(htr, h, G);
   invert(htr, htrinv);

   mx = ceil(cutoff*nx*moda(htrinv));
   my = ceil(cutoff*ny*modb(htrinv));
   mz = ceil(cutoff*nz*modc(htrinv));

   /*
    * Allocate and clear array for map of cells
    */
   nnab = 4*(mx+1)*(my+1)*(mz+1);
   cellmap = (int***)arralloc((size_mt)sizeof ***cellmap, 3, 
                              0, mx, -my-1, my, -mz-1, mz);
   memst(cellmap[0][-my-1]-mz-1,0, nnab*sizeof ***cellmap);

   /*
    * Add cells with corner-pair distances < cutoff
    */
#ifdef DEBUG1
   printf("  Distance    ix    iy    iz      sx        sy          sz\n");
#endif
   for(ix = 0; ix < mx; ix++)
      for(iy = (ix == 0 ? 0 : -my); iy < my; iy++)
         for(iz = (ix == 0 && iy == 0 ? 0 : -mz); iz < mz; iz++)
         {
            s[0] = (double)ix/nx;
            s[1] = (double)iy/ny;
            s[2] = (double)iz/nz;
            dist = 0.0;
            for(i = 0; i < 3; i++)
               for(j = 0; j < 3; j++)
                  dist += s[i]*G[i][j]*s[j];
            if(dist < SQR(cutoff))
            {
               for(i=0; i<=1; i++)
                  for(j = -1; j <= 1; j++)
                     for(k = -1; k <= 1; k++)
                        cellmap[ix+i][iy+j][iz+k] = 1;

#ifdef DEBUG1
               printf("%12f %4d %4d %4d %12f %12f %12f\n",
                      dist,ix,iy,iz,s[0],s[1],s[2]);
#endif
            }
         }
   /*
    * Add cells with face-face distance < cutoff.  Cells along x,y,z axes
    * are added in +/- directions, but only +ve ix indices added to map.
    */
   nxyz[0] = nx; nxyz[1] = ny; nxyz[2] = nz;
   mxyz[0] = mx; mxyz[1] = my; mxyz[2] = mz;
   for( ixyz=0; ixyz < 3; ixyz++)       /* Loop over directions */
   {
      jxyz = (ixyz+1) % 3; kxyz = (jxyz+1) % 3;
      proj[0] = proj[1] = proj[2] = 0.0;
      modabc = 0.0;
      for( i=0; i<3; i++ )
      {
         modabc += htrinv[i][ixyz];
         proj[i] += htrinv[i][ixyz]*htrinv[i][(ixyz+i) % 3];
      }
      for( i=0; i<3; i++ )
         proj[i] *= (mxyz[ixyz]-1)*nxyz[i]/(nxyz[ixyz] * modabc);
      /*
       * proj now contains projection vector.  Construct 4 candidate
       * cells.
       */
      for( i=0; i<3; i++ )
         face_cells[0][i] = face_cells[1][i] = face_cells[2][i] = 
            face_cells[3][i] = floor(proj[i]);
   
      face_cells[0][ixyz] = face_cells[1][ixyz] = face_cells[2][ixyz] = 
         face_cells[3][ixyz] = mxyz[ixyz];

      face_cells[1][jxyz] = face_cells[3][jxyz] = ceil(proj[jxyz]);
      face_cells[2][kxyz] = face_cells[3][kxyz] = ceil(proj[kxyz]);
      /*
       *  Now make sure we add only cells with +ve x index.  Add the
       *  inverse cell if ix<0.
       */
      for( i=0; i < 4; i++)
         if( face_cells[i][0] < 0 )
            for( j=0; j<3; j++ )
               face_cells[i][j] = -face_cells[i][j];
      /*
       * Now add the cells to the map.
       */
      for(  i=0; i < 4; i++ )
      {
         cellmap[face_cells[i][0]][face_cells[i][1]][face_cells[i][2]] = 1;
#ifdef DEBUG1
               printf("%12f %4d %4d %4d %12f %12f %12f\n",
                      0.5,face_cells[i][0],face_cells[i][1],face_cells[i][2],
                      0.0,0.0,0.0);
#endif
      }
   }
   /*
    * Scan map and build list.  N.B.  Loop indices are 1 greater than when
    * list built since we added cells outside original loop limits.
    */   
   nabor = ialloc(nnab);
   for(ix = 0; ix <= mx; ix++)
      for(iy = (ix == 0 ? 0 : -my-1); iy <= my; iy++)
         for(iz = (ix == 0 && iy == 0 ? 0 : -mz-1); iz <= mz; iz++)
         {
            if( cellmap[ix][iy][iz] )
            {
	       if( ix > IMCELL_XTRA*nx || ix < -IMCELL_XTRA*nx ||
		   iy > IMCELL_XTRA*ny || iy < -IMCELL_XTRA*ny ||
		   iz > IMCELL_XTRA*nz || iz < -IMCELL_XTRA*nz )
                  message(NULLI, NULLP, FATAL, CUTOFF, IMCELL_XTRA);
               if( inabor >= nnab )
                  message(NULLI, NULLP, FATAL,
                          "Internal error in neighbour_list()");
               nabor[inabor] = NCELL_IDX(ix,iy,iz);
	       
	       inabor++;
#ifdef DEBUG1
               printf("%12f %4d %4d %4d %12f %12f %12f\n",
                      1.0,ix,iy,iz,0.0,0.0,0.0);
#endif
            }
         }

   if( icheck )
   {
      if( inabor != onabor )
         note(NABORS,2 * inabor);
      onabor = inabor;
   }
   *nnabor = inabor;
   afree((gptr*)cellmap);
   return(nabor);
}
/******************************************************************************
 *  Fill_cells.  Allocate all the sites to cells depending on their centre of *
 *  mass co-ordinate by binning.                                              *
 ******************************************************************************/
static void    fill_cells(int nmols,         /* Number of molecules      (in) */
			  vec_mt (*c_of_m),  /* Centre of mass co-ords   (in) */
			  real **site, 	     /* Atomic site co-ordinates (in) */
			  spec_mp species,   /* Pointer to species array (in) */
			  mat_mt h, 	     /* Unit cell matrix         (in) */
                	  int nx, int ny, int nz,                                
			  cell_mt *lst,      /* Pile of cell structs     (in) */
			  cell_mt **cell,    /* Array of cells           (out)*/
			  int *frame_type)   /* Framework type counter   (out)*/
{
   int icell, imol, im, is, isite = 0;
   double eps = 8.0*precision();
   double	fnx = nx, fny = ny, fnz = nz;

   spec_mp spec = species;
   cell_mt *list = lst;
   vec_mt ssite;
   mat_mt hinv;

   *frame_type=1;
   invert(h, hinv);

   for(imol = 0, im = 0; imol < nmols; imol++, im++)
   {
      if(im == spec->nmols)
      {
         im = 0;
         spec++;
      }

      if( spec->framework )
         (*frame_type)++;
      if( ! control.molpbc || spec->framework )
      {
	 for( is = 0; is < spec->nsites; is++)
	 {
	    ssite[0] = site[0][isite];
	    ssite[1] = site[1][isite];
	    ssite[2] = site[2][isite];
	    mat_vec_mul(hinv, (vec_mt*)ssite, (vec_mt*)ssite, 1);
	    icell = LOCATE(ssite, eps);
	    list->isite = isite++;
	    list->num   = 1;
	    list->frame_type = *frame_type-1;
	    list->next = cell[icell];
	    cell[icell] = list++;
	 }
      }
      else
      {
         icell = LOCATE(c_of_m[imol], eps);
         list->isite = isite;
         list->num   = spec->nsites;
         list->frame_type = 0;
         list->next = cell[icell];
         cell[icell] = list++;
         list->next = NULL;
         isite += spec->nsites;
      }
   }
}
/******************************************************************************
 * Max_density()                                                              *
 *  Get the total number of sites in each cell and return the max. density    *
 ******************************************************************************/
double  max_density(cell_mt **cell, double vol, int ncells)
{
   int icell, nscell, nsmax=0;
   cell_mt *cell_p;

   for(icell = 0; icell < ncells; icell++)
   {
      nscell = 0;
      for(cell_p = cell[icell]; cell_p; cell_p = cell_p->next)
         nscell += cell_p->num;
      nsmax = MAX(nscell, nsmax);
   }
#if DEBUG8
   printf("Max sites/cell = %d\n", nsmax);
#endif
   return nsmax*ncells/vol;
}
/******************************************************************************
 *  site_neightbour list.  Build the list of sites withing interaction radius *
 *                         from the lists of sites in cells.                  *
 *	This version computes the PBC relocations of cells and translation    *
 *	vectors by computing the 1d index of the current cell assuming it is  *
 *      the central cell in a larger 3D cell extended by  IMCELL_XTRA in +ve  *
 *      and -ve x,y,z directions.  The neighbour cell in the extended 3D box  *
 *      is computed by simple addition of the 1-D indices.  Finally the index *
 *      of the neighbour cell in the original nx * ny * nz cell co-ordinates  *
 *      is looked up in a precomputed table, as is the label of which image   *
 *      of the central MD box it was located in.                              *
 *      This version is 2-3 times as fast as any previous version!            *
 ******************************************************************************/
int site_neighbour_list(int *nab,          /* Array of sites in list    (out) */
			int *pbctrans,     /* index of pbc translations (out) */ 
			int n_nab_sites,   /* Size of above arrays       (in) */
			int *nfnab,        /*  sites index by type      (out) */
			int n_frame_types, /* # of distinct frameworks   (in) */
			int n_nabors,      /* Number of neighbour cells  (in) */
			int ix, int iy, int iz,    /*  Current cell      (in) */
			int nx, int ny, int nz,    /*  # of subcells     (in) */
			int *nabor,        /* List of neighbour cells    (in) */
			cell_mt **cell,    /* Head of cell list          (in) */
			int (*pbclookup)[2]) /* 3D index lookup table    (in) */
{
   int  j0, jsite, jnab;                 /* Counters for cells etc            */
   int  nnab=0;				 /* Counter for size of nab           */
   int  ftype;
   cell_mt      *cmol;                   /* Pointer to current cell element   */
   int icell_idx, jcell_idx, ktrans;

   icell_idx = NCELL_IDX(ix+IMCELL_XTRA*nx,iy+IMCELL_XTRA*ny,iz+IMCELL_XTRA*nz);
 
   for(ftype = 0; ftype < n_frame_types; ftype++) /* Do Framework types first */
   {
      j0 = (ftype == 0) ? 0: 1;
      for(jnab = j0; jnab < n_nabors; jnab++)    /* Loop over neighbour cells  */
      {
	 jcell_idx = icell_idx + nabor[jnab];
	 ktrans  = pbclookup[jcell_idx][1];
	 jcell_idx=pbclookup[jcell_idx][0];
         /* Loop over molecules in this cell, filling 'nab' with its sites    */
         for(cmol = cell[jcell_idx]; cmol != 0; cmol = cmol->next)
         {
            if( cmol->frame_type == ftype )
            {
               for(jsite = 0; jsite < cmol->num; jsite++)
               {
		  nab[nnab] = cmol->isite + jsite;
		  pbctrans[nnab] = ktrans;
		  nnab++;
	       }
            }
            if(nnab > n_nab_sites) 
               message(NULLI,NULLP,FATAL,TONAB,nnab,n_nab_sites);
         }
      }
      nfnab[ftype] = nnab;
   }

   return nnab;
}
#ifdef DEBUG2
/******************************************************************************
 * minimage.  Calculate "minimum-image"  histogram and energy for debugging   *
 ******************************************************************************/
void minimage(system_mt *system,        /* System struct                 (in) */
	      spec_mt *species,         /* Array of species records      (in) */
	      real **site, 	        /* Site co-ordinate arrays       (in) */
	      real *chg, 	        /* Array of site charges         (in) */
	      real ***potp,	        /* Array of potential parameters (in) */
	      int *id)		        /* Site identifier array.        (in) */
{
   spec_mt *spec;
   double ppe, rr[3], ss[3];
   int i, im, is, isite, imol, jsite;
   double       norm = 2.0*control.alpha/sqrt(PI);      /* Coulombic prefactor*/
   double       cutoffsq = SQR(control.cutoff), /* Temporary copy for optim'n */
                cutoff100sq = 10000.0*cutoffsq;
   real         *r_sqr   = dalloc(system->nsites), /* Squared site-site distance */
                *forceij = dalloc(system->nsites); /* -V'(r) / r                 */
   mat_mp h = system->h;
   mat_mt hinv;

   ppe = 0;
   spec = species; isite = 0;
   invert(h, hinv);
   for(imol = 0, im = 0; imol < system->nmols; imol++, im++)
   {
      if(im == spec->nmols)
      {
         im = 0;
         spec++;
      }
      for(is = isite; is < isite+spec->nsites; is++)
      {
         for(jsite = 0; jsite < isite; jsite++)
         {
            for( i=0; i<3; i++)
               rr[i] = site[i][jsite] - site[i][is];
            mat_vec_mul(hinv, rr, ss, 1);
            for( i=0; i<3; i++)
               ss[i] -= floor(ss[i]+0.5);
            mat_vec_mul(h, ss, rr, 1);

            r_sqr[jsite] = SUMSQ(rr);
            if( control.strict_cutoff && r_sqr[jsite] > cutoffsq )
                  r_sqr[jsite] = cutoff100sq;
         }

         hist(0,isite,r_sqr);
         kernel(0,isite,forceij,&ppe,r_sqr,chg,chg[is],
                control.alpha,control.alpha46,system->ptype,potp[id[is]]);
      }
      isite += spec->nsites;
   }
   histout();
   note("Real-space intermolecular potential energy = %g",ppe*CONV_E);
   xfree(r_sqr); xfree(forceij);
   }
#endif
/******************************************************************************
 * Poteval            Return potential evaluated at a single point.           *
 ******************************************************************************/
double poteval(real *potpar,            /* Array of potential parameters      */
	       double r,                /* Cutoff distance                    */
	       int ptype,               /* Potential type selector            */
	       double chgsq)            /* Product of site charges            */
{
   double pe = 0.0;
   real   chgsq_r = chgsq;
   real f,rr;
   real *pp[NPOTP];
   int  i;

   for(i=0; i<NPOTP; i++)
      pp[i] = potpar+i;

   rr = SQR(r);
   kernel(0,1,&f,&pe,&rr,&chgsq_r,1.0,control.alpha,control.alpha46,ptype, pp);
   return pe;
}
/******************************************************************************
 * pe_intra.  Caculate intra-molecular energy of a molecule
 ******************************************************************************/
double pe_intra(spec_mt *spec, real chg[], int ptype, 
                pot_mt *potpar, int max_id)
{
   double eintra=0.0;
   int  jsite, isite;

   for(jsite = 0; jsite < spec->nsites; jsite++)
   {
      for(isite = jsite+1; isite < spec->nsites; isite++)
	{
         eintra += poteval(potpar[spec->site_id[jsite]*max_id+spec->site_id[isite]].p,
                           DISTANCE(spec->p_f_sites[isite],spec->p_f_sites[jsite]),
                           ptype,chg[isite]*chg[jsite]);
        }
   }
   return eintra;
}
/******************************************************************************
 * Debug routines
 ******************************************************************************/
#ifdef DEBUG3
void dump_neighbour_list(int n, int jmin, int isite, int nab[], int pbctrans[], double xi, double yi, double zi,
                         double x[], double y[], double z[], double r_sqr[], int id[])
{
   int jnab;
   if( n <= 0)
      return;
   printf("      i(type)    j(type)     k      (x,y,z)                dist\n");
   printf("Ref (%6.3f,%6.3f,%6.3f)\n--------------------------------------------------\n",xi,yi,zi);
   for(jnab = 0; jnab < n; jnab++)
      printf("%c%4d %4d(%3d) %4d(%3d) %3d (%6.3f,%6.3f,%6.3f)  %12.5f\n", 
             jnab<jmin?'*':' ',
             jnab,isite,id[isite],nab[jnab],id[nab[jnab]],pbctrans[jnab],
             x[jnab],y[jnab],z[jnab], sqrt(r_sqr[jnab]));
   printf("\n");
}
#endif
/******************************************************************************
 * mk_r_sqr().  Compute site-site vectors and squared distance between sites  *
 ******************************************************************************/
void mk_r_sqr(int jmax, int pbctrans[], rvec_mt reloc_v[], 
	      real nab_sx[], real nab_sy[], real nab_sz[], 
	      real site0, real site1, real site2, 
	      real r_sqr[], real rx[], real ry[], real rz[])
{
   int jsite, k;
   real rrx, rry, rrz;

   for(jsite=0; jsite < jmax; jsite++)
   {
      k = pbctrans[jsite];
      rrx = nab_sx[jsite] - site0 + reloc_v[k].x;
      rry = nab_sy[jsite] - site1 + reloc_v[k].y;
      rrz = nab_sz[jsite] - site2 + reloc_v[k].z;
      r_sqr[jsite] = rrx*rrx+rry*rry+rrz*rrz;
      rx[jsite] = rrx;
      ry[jsite] = rry;
      rz[jsite] = rrz;
   }
}

/******************************************************************************
 * mk_forces().  Compute vector forces on site i and neighbour sites          *
 *     Expanded using temporary vars to permit optimization.		      *
 ******************************************************************************/
void mk_forces(int jmin, int jmax,  real rx[], real ry[], real rz[], real forceij[],
	       real forcejx[], real forcejy[], real forcejz[],
	       real sf0[], real sf1[], real sf2[])
{
   int jsite;
   real force_cptx, force_cpty, force_cptz, site0, site1, site2;

   site0 = site1 = site2 = 0.0;
   for(jsite=jmin; jsite < jmax; jsite++)
   {
      force_cptx = forceij[jsite]*rx[jsite];
      force_cpty = forceij[jsite]*ry[jsite];
      force_cptz = forceij[jsite]*rz[jsite];
      site0           -= force_cptx;
      site1           -= force_cpty;
      site2           -= force_cptz;
      force_cptx      += forcejx[jsite];
      force_cpty      += forcejy[jsite];
      force_cptz      += forcejz[jsite];
      forcejx[jsite]   = force_cptx;
      forcejy[jsite]   = force_cpty;
      forcejz[jsite]   = force_cptz;
   }
   sf0[0] += site0;
   sf1[0] += site1;
   sf2[0] += site2;
}
/******************************************************************************
 * scatter_forces().  Add forces on sites in neighbour list to main force     *
 *     array.  Also compute Bekker's "g" forces, ie forces on each pbc image  *
 *     surrounding the main MD cell.                                          *
 *     Expanded using temporary vars to permit optimization.		      *
 ******************************************************************************/
void scatter_forces(int nnab, int nab[], int pbctrans[], real gforce[NIMCELLS][3],
		    real forcejx[], real forcejy[], real forcejz[], 
		    real sf0[], real sf1[], real sf2[])
{
   int inab, it, k;
   real *gfk, fx, fy, fz, sfx, sfy, sfz;

   for(inab=0; inab < nnab; inab++)
   {
      fx = forcejx[inab];
      fy = forcejy[inab];
      fz = forcejz[inab];
      it = nab[inab];
      k  = pbctrans[inab];
      gfk = gforce[k];
      sfx       = fx + sf0[it];
      sfy       = fy + sf1[it];
      sfz       = fz + sf2[it];
      fx       += gfk[0];
      fy       += gfk[1];
      fz       += gfk[2];
      sf0[it]  = sfx;
      sf1[it]  = sfy;
      sf2[it]  = sfz;
      gfk[0]    = fx;
      gfk[1]    = fy;
      gfk[2]    = fz;
   }
}
/******************************************************************************
 *  Force_inner() Paralellised inner loops of force_calc.  Loops over cells   *
 *  in MD cell with stride = nomber of processors available.  Should be       *
 *  called once for each parallel thread.                                     *
 ******************************************************************************/
void force_inner(int ithread,            
		 int nthreads,          /* Parallel node variables.      (in) */ 
		 real **site, 	        /* Site co-ordinate arrays       (in) */
		 real *chg, 	        /* Array of site charges         (in) */
		 real ***potp, 	        /* Expanded potential parameter array */
		 int *id, 	        /* Array of site_id[nsites]      (in) */
		 int n_nab_sites,       /* Dimension of site n'bor list arrays*/
		 int n_nabors,	        /* Number of elements in lists.   (in)*/
		 int *nabor, 	        /* Lists of neighbour cells       (in)*/
		 int nx, int ny, int nz,/* Number of subcells in MD cell  (in)*/
		 cell_mt **cell,        /* Array of list heads of subcells(in)*/
		 int n_frame_types,     /* ==1 for no fw, 2 if fw present (in)*/
		 system_mt *system,     /* System struct                 (in) */
		 mat_mt stress,         /* Stress virial                (out) */
		 double *pe, 	        /* Potential energy             (out) */
		 real **site_force,     /* Site force arrays            (out) */
		 int *molmap,           /* Which molecule site belongs to(in) */
		 rvec_mt *reloc_v,      /* Table of PBC reloc. vectors.  (in) */
		 int (*pbclookup)[2])   /* Extended 3D cell lookup table (in) */
{
                /*
                 * The following arrays are for 'neighbour site list'
                 * quantities and should be dimensioned to the max value of
                 * 'nnab'.  A rough approx is the ratio of the volume of
                 * the "cutoff sphere" to that of the MD cell times nsites.
                 * This may be too small for inhomogeneous systems, but at
                 * least it scales with the cutoff radius.
                 */
   int          *nab  = ialloc(n_nab_sites),    /* Neigbour site gather vector*/
                *pbctrans= ialloc(n_nab_sites);
   real         *nab_sx  = dalloc(n_nab_sites), /* 'Gathered' list of         */
                *nab_sy  = dalloc(n_nab_sites), /*   neighbour site co-ords   */
                *nab_sz  = dalloc(n_nab_sites), /*   - x,y,z components.      */
                *forcejx = dalloc(n_nab_sites), /* List of neighbour site     */
                *forcejy = dalloc(n_nab_sites), /*  forces in gathered form   */
                *forcejz = dalloc(n_nab_sites), /*  - xyz components.         */
                *rx      = dalloc(n_nab_sites), /* Reference to neigbour site */
                *ry      = dalloc(n_nab_sites), /* - site vector adjusted for */
                *rz      = dalloc(n_nab_sites), /*  periodic boundaries. xyz. */
                *r_sqr   = dalloc(n_nab_sites), /* Squared site-site distance */
                *nab_chg = dalloc(n_nab_sites), /* Gathered neig. site charges*/
                *forceij = dalloc(n_nab_sites); /* -V'(r) / r                 */
   real         **nab_pot                       /* Gathered pot par array     */
                = (real**)arralloc((size_mt)sizeof(real), 2,
                                   0, system->n_potpar-1, 0, n_nab_sites-1);
   real         **pp, **ppp;            /* Loop pointer variables for potp.   */
   real         s00, s01, s02, s11, s12, s22;
                                   /* Accumulators for forces and stresses.   */
   real         gforce[NIMCELLS][3];
   double       norm = 2.0*control.alpha/sqrt(PI);      /* Coulombic prefactor*/
   double       cutoffsq = SQR(control.cutoff), /* Temporary copy for optim'n */
                cutoff100sq = 10000.0*cutoffsq;
   int          ix, iy, iz;             /* 3-d cell indices for ref and neig. */
   int          icell,                  /* Index for cells of molecule pair   */
                nnab, jmin, jmax,       /* Number of sites in neighbour list  */
                isite, jsite, ipot, lim;/* Counters.                          */
   int          nsites = system -> nsites;      /* Temporary copy for optim'n */
   int          nfnab[2];               /* Number of non-fw and fw neighbours */
   int          k;
   cell_mt      *cmol;                  /* Loop counter for link cells.       */

   for(k = 0; k < NIMCELLS; k++)
      gforce[k][0] = gforce[k][1] = gforce[k][2] = 0.0;
/******************************************************************************
 *  Start of main loop over all subcells.                                     *
 *  First build "site neighbour list" containing all sites belonging to       *
 *  molecules in this subcell and all others in the cell neighbour list.      *
 *  Use "gather" to construct corresponding arrays of co-ordinates, charges   *
 *  and potential parameters.                                                 *
 *  Then loop over all sites in THIS cell and calculate pair distances,       *
 *  potential forces and stress.                                              *
 ******************************************************************************/
   for( icell = ithread; icell < nx*ny*nz; icell += nthreads)
   {
      if(cell[icell] == NULL) continue;       /* Empty cell - go on to next */
      ix = icell/ (ny*nz);
      iy = icell/nz - ny*ix;
      iz = icell - nz*(iy + ny*ix);
#ifdef DEBUG3
      printf("Working on cell %4d (%d,%d,%d) (sites %4d to %4d)\n", icell,
             ix,iy,iz,cell[icell]->isite,cell[icell]->isite+cell[icell]->num-1);
      printf("\n jcell\tjx jy jz\tNsites\n");
#endif
      /*
       * Build site neighbour list 'nab' from cell list.
       */ 
      nnab = site_neighbour_list(nab, pbctrans, n_nab_sites, nfnab, n_frame_types, 
                                 n_nabors, ix, iy, iz, nx, ny, nz, 
				 nabor, cell, pbclookup);
#ifdef DEBUG3
      if( nnab > 0 )
         printf(" %d entries in neighbour list.\n",nnab);
#endif
      gather(nnab, nab_sx, site[0], nab, nsites); /* Construct list of site  */
      gather(nnab, nab_sy, site[1], nab, nsites); /* co-ordinates from nabor */
      gather(nnab, nab_sz, site[2], nab, nsites); /* list.                   */
         
      gather(nnab, nab_chg, chg, nab, nsites); /* Gather site charges as well*/
      zero_real(forcejx,nnab);
      zero_real(forcejy,nnab);
      zero_real(forcejz,nnab);
         
      /*
       * Main loop over sites - generate "isite" index from contents of cell.
       */
      jmin = 0;
      for(cmol = cell[icell]; cmol != NULL; cmol = cmol->next)
      {
         /*
          * The limits of the inner loops (jmin,jmax) are complicated to handle
          * frameworks.  For a *non* framework site, isite, the lower limit
          * jmin is incremented by one for each site in this reference cell
          * (and only in the central box) to avoid the self-term and double
          * counting.  That is it implements sum i=1,N; j=i+1,N.
          * Framework sites are handled differently.  site_neighbour_list 
          * orders all framework sites at the end of the list and sets nfnab[0]
          * to delimit their beginning.  Therefore all interactions between
          * non-framework and framework sites are included and intra-framework
          * one excluded by looping from 0 to nfnab[0].
          */
         if( cmol->frame_type )
         {
            jmin = 0;
            jmax = nfnab[0];
         }
         else
         {
            jmax = nnab;
         }
         lim = cmol->isite + cmol->num;
         for(isite = cmol->isite; isite < lim; isite++)
         {                                   /* Loop over sites in molecule */
            if( ! cmol->frame_type )        /* Avoid double-counting interactions */
               jmin++;                  /* of particles in central cell       */
#ifdef DEBUG1
            printf("icell %d:\tisite=%d\tjmin=%d\n",icell,isite,jmin);
#endif
            /*
             * Construct pot'l param arrays corresponding to neighbour sites.
             */
            pp = potp[id[isite]];
            ppp = nab_pot;
            for(ipot = 0; ipot < system->n_potpar; ipot++)
               gather(jmax, *ppp++, *pp++, nab, nsites);
            
	    mk_r_sqr(jmax, pbctrans, reloc_v, nab_sx, nab_sy,  nab_sz, 
		     site[0][isite], site[1][isite], site[2][isite],
		     r_sqr,  rx,  ry,  rz);
#ifdef DEBUG3X
            if(isite == 15 || isite==16)
#endif
#if defined(DEBUG3) 
               dump_neighbour_list(jmax,jmin,isite,nab,pbctrans,site0, site1, site2,
                                   rx, ry, rz, r_sqr,id);
#endif
	    jsite = jmin;
	    do {
	       jsite += search_lt(jmax-jsite, r_sqr+jsite, 1, TOO_CLOSE);
	       if( jsite < jmax )
	       {
		  if( molmap[isite] != molmap[nab[jsite]])
		     message(NULLI, NULLP, WARNING, TOOCLS,
			     isite, nab[jsite], sqrt(TOO_CLOSE));
		  jsite++;
	       }
	    } while (jsite < jmax );
               
            if( control.strict_cutoff && ! control.molpbc )
               for(jsite = jmin; jsite < jmax; jsite++)
                  if( r_sqr[jsite] > cutoffsq )
                     r_sqr[jsite] = cutoff100sq;
               
#ifdef DEBUG2
            hist(jmin, jmax, r_sqr);
#endif
	    /*  Call the potential function kernel                            */

            kernel(jmin, jmax, forceij, pe, r_sqr, nab_chg, chg[isite],
                   control.alpha, control.alpha46, system->ptype, nab_pot);

	    mk_forces(jmin, jmax,  rx, ry, rz, forceij, forcejx, forcejy, forcejz,
		      site_force[0]+isite, site_force[1]+isite, site_force[2]+isite);
#ifdef DEBUG5
	    printf("PE = %f\n",pe[0]);
#endif
         }
      }
      scatter_forces(nnab, nab, pbctrans, gforce, forcejx, forcejy, forcejz,
		     site_force[0], site_force[1], site_force[2]);
   }

   s00 = s01 = s02 = s11 = s12 = s22 = 0.0;     /* Accumulators for stress    */
   for(isite = 0; isite < nsites; isite++)
   {
      s00 += site[0][isite]*site_force[0][isite];
      s01 += site[1][isite]*site_force[0][isite];
      s02 += site[2][isite]*site_force[0][isite];
      s11 += site[1][isite]*site_force[1][isite];
      s12 += site[2][isite]*site_force[1][isite];
      s22 += site[2][isite]*site_force[2][isite];
   }
   for(k = 0; k < NIMCELLS; k++)
   {
      s00 += reloc_v[k].x*gforce[k][0];
      s01 += reloc_v[k].y*gforce[k][0];
      s02 += reloc_v[k].z*gforce[k][0];
      s11 += reloc_v[k].y*gforce[k][1];
      s12 += reloc_v[k].z*gforce[k][1];
      s22 += reloc_v[k].z*gforce[k][2];
   }
   stress[0][0]  += s00;
   stress[0][1]  += s01;
   stress[0][2]  += s02;
   stress[1][1]  += s11;
   stress[1][2]  += s12;
   stress[2][2]  += s22;

   afree((gptr*)nab_pot);
   xfree(nab);     xfree(pbctrans);xfree(nab_chg);
   xfree(r_sqr);   xfree(forceij);
   xfree(rx);      xfree(ry);      xfree(rz);
   xfree(forcejx); xfree(forcejy); xfree(forcejz);
   xfree(nab_sx);  xfree(nab_sy);  xfree(nab_sz);
}
/******************************************************************************
 *  Rdf_inner() Paralellised inner loops of force_calc.  Based on force_inner *
 *     but only calls rdf_accum().                                            *
 ******************************************************************************/
void rdf_inner(int ithread, 
	       int nthreads,            /* Parallel node variables.      (in) */
	       real **site, 	        /* Site co-ordinate arrays       (in) */
	       int *id, 	        /* Array of site_id[nsites]      (in) */
	       int n_nab_sites,         /* Dimension of site n'bor list arrays*/
	       int n_nabors, 	        /* Number of elements in lists.   (in)*/
	       int *nabor, 	        /* Lists of neighbour cells       (in)*/
	       int nx, int ny, int nz,  /* Number of subcells in MD cell  (in)*/
	       cell_mt **cell, 	        /* Array of list heads of subcells(in)*/
	       int n_frame_types,       /* ==1 for no fw, 2 if fw present (in)*/
	       system_mt *system,       /* System struct                  (in)*/
	       rvec_mt *reloc_v,        /* Table of PBC reloc. vectors.  (in) */
	       int (*pbclookup)[2])     /* Extended 3D cell lookup table (in) */
{
   int          *nab  = ialloc(n_nab_sites),    /* Neigbour site gather vector*/
                *pbctrans= ialloc(n_nab_sites);
   real         *nab_sx  = dalloc(n_nab_sites), /* 'Gathered' list of         */
                *nab_sy  = dalloc(n_nab_sites), /*   neighbour site co-ords   */
                *nab_sz  = dalloc(n_nab_sites), /*   - x,y,z components.      */
                *r_sqr   = dalloc(n_nab_sites); /* Squared site-site distance */
   real         site0, site1, site2;
   real         rrx, rry, rrz;                  /* Scalar loop temporaries    */
   int          ix, iy, iz;             /* 3-d cell indices for ref and neig. */
   int          icell,                  /* Index for cells of molecule pair   */
                nnab, jmin, jmax,       /* Number of sites in neighbour list  */
                isite, jsite, lim;      /* Counters.                          */
   int          nsites = system -> nsites;      /* Temporary copy for optim'n */
   int          nfnab[2];               /* Number of non-fw and fw neighbours */
   int          k;
   cell_mt      *cmol;                  /* Loop counter for link cells.       */
/******************************************************************************
 *  Start of main loop over subcells.                                         *
 ******************************************************************************/
   for( icell = ithread; icell < nx*ny*nz; icell += nthreads)
   {
      if(cell[icell] == NULL) continue;       /* Empty cell - go on to next */
      ix = icell/ (ny*nz);
      iy = icell/nz - ny*ix;
      iz = icell - nz*(iy + ny*ix);

      /*
       * Build site neighbour list 'nab' from cell list.
       */ 
      nnab = site_neighbour_list(nab, pbctrans, n_nab_sites, nfnab, n_frame_types, 
                                 n_nabors, ix, iy, iz, nx, ny, nz, 
				 nabor, cell, pbclookup);
      gather(nnab, nab_sx, site[0], nab, nsites); /* Construct list of site  */
      gather(nnab, nab_sy, site[1], nab, nsites); /* co-ordinates from nabo  */
      gather(nnab, nab_sz, site[2], nab, nsites); /* list.                   */

      /*
       * Main loop over sites - generate "isite" index from contents of cell.
       */
      jmin = 0;
      for(cmol = cell[icell]; cmol != NULL; cmol = cmol->next)
      {
         if( cmol->frame_type )
         {
            jmin = 0;
            jmax = nfnab[0];
         }
         else
         {
            jmax = nnab;
         }
         lim = cmol->isite + cmol->num;
         for(isite = cmol->isite; isite < lim; isite++)
         {                                   /* Loop over sites in molecule */
            if( ! cmol->frame_type )        /* Avoid double-counting interactions */
               jmin++;                  /* of particles in central cell       */
            
            site0=site[0][isite]; 
            site1=site[1][isite]; 
            site2=site[2][isite];
            for(jsite=0; jsite < jmax; jsite++)
            {
               k = pbctrans[jsite];
               rrx = nab_sx[jsite] - site0 + reloc_v[k].x;
               rry = nab_sy[jsite] - site1 + reloc_v[k].y;
               rrz = nab_sz[jsite] - site2 + reloc_v[k].z;
               r_sqr[jsite] = rrx*rrx+rry*rry+rrz*rrz;
            }
            /*
             * Accumulate radial distribution functions
             */
            rdf_accum(jmin, jmax, r_sqr, id[isite], id, nab);
         }
      }
   }
   xfree(nab);   xfree(pbctrans);
   xfree(r_sqr);
   xfree(nab_sx);  xfree(nab_sy);  xfree(nab_sz);
}
/******************************************************************************
 * Force_calc.   This is the main intermolecular site force calculation       *
 * routine                                                                    *
 ******************************************************************************/
void force_calc(real **site,            /* Site co-ordinate arrays       (in) */
		real **site_force,      /* Site force arrays            (out) */
		system_mt *system,      /* System struct                 (in) */
		spec_mt *species,       /* Array of species records      (in) */
		real *chg, 	        /* Array of site charges         (in) */
		pot_mt *potpar,         /* Array of potential parameters (in) */
		double *pe, 	        /* Potential energy             (out) */
		mat_mt stress)	        /* Stress virial                (out) */
{
   int          isite, imol,            /* Site counter i,j                   */
                i_id, ipot;             /* Miscellaneous                      */
   int          n_frame_types;          /* ==1 for no fw, 2 if fw present.    */
   int          nsites = system->nsites,/* Local copy to keep optimiser happy */
                n_potpar = system->n_potpar,
                max_id = system->max_id;
   int          *molmap = ialloc(nsites);
   double       cutoff = control.cutoff; /* + (control.strict_cutoff?mol_diam:0);*/
   rvec_mt      reloc_v[NIMCELLS];
   int          ii,jj,kk, k;
   double       subcell = control.subcell; /* Local copy. May change it.      */
   int          *id   = ialloc(nsites), /* Array of site_id[nsites]           */
                *id_ptr;                /* Pointer to 'id' array              */
   int          n_nab_sites;            /* Dimension of site n'bor list arrays*/
   int          *nabor, *rdf_nabor;     /* Lists of neighbour cells           */
   int          n_nabors, n_rdf_nabors; /* Number of elements in lists.       */
   int          icell, ncells;          /* Subcell counter and total = nxnynz */
   int          n_cell_list;            /* Size of link-cells "heap"          */
   int          nx, ny, nz;             /* Number of subcells in MD cell      */
   static int   onx=0, ony=0, onz=0;    /* Saved values of nx, ny, nz.        */
   real         ***potp                 /* Expanded potential parameter array */
                      = (real***)arralloc((size_mt)sizeof(real), 3, P0,max_id-1,
                                          0, n_potpar-1, 0, nsites-1);
   cell_mt      *c_ptr;                 /* Heap of link cell entries for list */
   cell_mt      **cell;                 /* Array of list heads for subcells   */
   spec_mt      *spec;                  /* Temp. loop pointer to species.     */
   int          jsite, jmol;
   int		icell4d, ii0, ix, iy, iz;
   double       vol = det(system->h);
   static int init =1;
   static       double eintra;
   static       int   (*pbclookup)[2];
   static	int   imcell_offset;
   /*
    * Choose a partition into subcells if none specified.
    */
   if(subcell <= 0.0) subcell = control.cutoff/5.0;
   nx = system->h[0][0]/subcell+0.5;
   ny = system->h[1][1]/subcell+0.5;
   nz = system->h[2][2]/subcell+0.5;
   ncells = nx*ny*nz;
   if( init ) {
      int isite=0;
      for(spec = species; spec < species+system->nspecies; spec++)
      {
         if( ! spec->framework) 
            eintra+= spec->nmols
               *pe_intra(spec, chg+isite, system->ptype, potpar,max_id);
         isite += spec->nmols*spec->nsites;
      }
      note("Intramolecular potential energy correction = %g",eintra*CONV_E);
      init=0;
   }
   if( ithread == 0 )
      *pe -= eintra;

   jsite = 0; 
   jmol = 0;
   for (spec = species; spec < species+system->nspecies; spec++)
      for(imol = 0; imol < spec->nmols; imol++)
      {
         for(isite=0; isite < spec->nsites; isite++)
            molmap[jsite++] = jmol;
         jmol++;
      }
   /*  
    * Construct and fill expanded site-identifier array, id   
    */
   id_ptr = id;
   for (spec = species; spec < species+system->nspecies; spec++)
      for(imol = 0; imol < spec->nmols; imol++)
      {
         memcp(id_ptr, spec->site_id, spec->nsites*sizeof(int));
         id_ptr += spec->nsites;
      }
   /*   
    * Build arrays of pot. pars [max_id][nsites] for access in vector loops
    */
   for(ipot = 0; ipot < n_potpar; ipot++)
      for(i_id = 1; i_id < max_id; i_id++)
      {
#ifdef titan
NOVECTOR
#endif
         for(isite = 0; isite < nsites; isite++)
            potp[i_id][ipot][isite] = potpar[i_id*max_id+id[isite]].p[ipot];
      }
   /*
    * Allocate "heap" of list entries to build linked lists from.
    */
   n_cell_list = nsites;              /* This is always safe */
   c_ptr = aalloc(n_cell_list, cell_mt); 
   /*
    * Build a linked list of molecules/sites for each subcell.
    * "cell" is the array of list heads NX x NY x NZ.
    */
   cell = aalloc(ncells, cell_mt *);     
   for( icell=0; icell < ncells; icell++)
      cell[icell] = NULL;
   fill_cells(system->nmols, system->c_of_m, site, species, system->h,
              nx, ny, nz, c_ptr, cell, &n_frame_types);
   if( n_frame_types > 2 )
      message(NULLI, NULLP, FATAL,
              "Multiple framework molecules are not supported");
   /*
    * Build lists of cells within cutoff of reference cell.
    */
   if( control.strict_cutoff )
      nabor = strict_neighbour_list(&n_nabors, system->h, cutoff, nx, ny, nz, 1);
   else
      nabor = neighbour_list(&n_nabors, system->h, cutoff, nx, ny, nz, 1);
   
   /*
    * Calculate size needed for site neighbour list arrays
    */
   if(control.strict_cutoff)
      n_nab_sites=NMULT*4.19*CUBE(cutoff+subcell)*max_density(cell, vol, ncells);
   else
      n_nab_sites=NMULT*4.19*CUBE(cutoff)*max_density(cell, vol, ncells);
#ifdef DEBUG2 
   n_nab_sites = NMULT*nsites;
   minimage(system, species, site, chg, potp, id);
#endif
   if( nx != onx || ny != ony || nz != onz )
   {
      note("MD cell divided into %d subcells (%dx%dx%d)",ncells,nx,ny,nz);
      onx = nx; ony = ny; onz = nz;

      if( pbclookup )
	 xfree(pbclookup+imcell_offset);

      /*
       * Allocate and build extended 3d ->4d image cell relocation table.
       * To save 1/3 of the space unused because neighbour lists always have
       * ix > 0, we don't start at zero.  
       * Revert to ANSI-compliant llim=0 if bounds checking.
       */
#if __BOUNDS_CHECKING_ON
      imcell_offset = 0;
      ii0 = 0;
#else
      imcell_offset=IMCELL_XTRA*IMCELL_L*IMCELL_L*ncells;
      ii0 = IMCELL_XTRA;
#endif
      /*
       * N.B.  Use alloc of int[2] rather than 2-D dope-vector array
       *       for efficieny.
       */
      pbclookup = (int (*)[2])arralloc(2*sizeof(int), 1, 
				       imcell_offset, NIMCELLS*ncells-1);
      
      icell4d=imcell_offset;
      for(ii = ii0; ii < IMCELL_L; ii++)
	 for(ix = 0; ix < nx; ix++)
	    for(jj = 0; jj < IMCELL_L; jj++)
	       for(iy = 0; iy < ny; iy++)
		  for(kk = 0; kk < IMCELL_L; kk++)
		     for(iz = 0; iz < nz; iz++)
		     {
			pbclookup[icell4d][0] = NCELL(ix,iy,iz);
			pbclookup[icell4d][1] = IMCELL_L*(IMCELL_L*ii+jj)+kk;
			icell4d++;
		     }
   }
   /*
    * Build the next image cell translation table. 
    */
   k=0;
   for(ii = -IMCELL_XTRA; ii <= IMCELL_XTRA; ii++)
      for(jj = -IMCELL_XTRA; jj <= IMCELL_XTRA; jj++)
         for(kk = -IMCELL_XTRA; kk <= IMCELL_XTRA; kk++)
         {
            reloc_v[k].x = system->h[0][0]*ii+system->h[0][1]*jj+system->h[0][2]*kk;
            reloc_v[k].y = system->h[1][0]*ii+system->h[1][1]*jj+system->h[1][2]*kk;
            reloc_v[k].z = system->h[2][0]*ii+system->h[2][1]*jj+system->h[2][2]*kk;
            k++;
         }

   force_inner(ithread, nthreads, site, chg, potp, id, n_nab_sites, 
               n_nabors, nabor, nx, ny, nz, cell, n_frame_types, system,
               stress, pe, site_force, molmap, reloc_v, pbclookup);

   /*
    * Accumulate radial distribution functions
    */
   if (control.rdf_interval > 0 && 
       control.istep >= control.begin_rdf &&
       control.istep % control.rdf_interval == 0)
   {
      n_nab_sites = NMULT*4.19*CUBE(control.limit)
                                               *max_density(cell, vol, ncells);
      rdf_nabor = strict_neighbour_list(&n_rdf_nabors, system->h, 
                                        control.limit, nx, ny, nz, 0);
      rdf_inner(ithread, nthreads, site, id, n_nab_sites, n_rdf_nabors, 
                rdf_nabor, nx, ny, nz, cell, n_frame_types, system, reloc_v, pbclookup);
      xfree(rdf_nabor);
   }
#ifdef DEBUG2
   histout();
#endif
   afree((gptr*)(potp+P0));  xfree(c_ptr); 
   xfree(cell);        xfree(id); 
   xfree(nabor);       xfree(molmap);
}

