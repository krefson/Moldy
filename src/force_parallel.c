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
 *       $Log: force_parallel.c,v $
 *       Revision 2.9  1997/07/09 14:45:19  keith
 *       Brought up to daye with main-line developments in force.c.
 *       (Actually, re-introduced compiler parallelism into force.c 2.14)
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
static char *RCSid = "$Header: /home/eeyore_data/keith/md/moldy/RCS/force_parallel.c,v 2.9 1997/07/09 14:45:19 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include        "defs.h"
/*========================== Library include files ===========================*/
#ifdef  stellar
#include        <fastmath.h>
#else
#include        <math.h>
#endif
#include        "stddef.h"
#include        "string.h"
#include        <assert.h>
/*========================== Program include files ===========================*/
#include        "structs.h"
#include        "messages.h"
/*========================== External function declarations ==================*/
gptr            *talloc();             /* Interface to memory allocator       */
void            tfree();               /* Free allocated memory               */
void            afree();               /* Free allocated array                */
int             search_lt();            /* Search a vector for el. < scalar   */
void            gather();               /* Interface to CRAY gather routine   */
void            mat_mul();              /* Matrix multiplier                  */
double          det();                  /* Determinant of 3x3 matrix          */
void            invert();               /* 3x3 matrix inverter                */
void            mat_vec_mul();          /* Matrix by vector multiplier        */
void            transpose();            /* Generate 3x3 matrix transpose      */
void            zero_real();            /* Initialiser                        */
void    	zero_double();          /* Initialiser                        */
void            force_inner();          /* Inner loop forward reference       */
void            rdf_inner();            /* RDF calc forward reference         */
double          precision();            /* Floating pt precision.             */
void            kernel();               /* Force kernel routine               */
double          mol_radius();           /* Radius of largest molecule.        */
void            rdf_accum();            /* Bin distances for rdf evaluation.  */
#ifdef HAVE_STDARG_H
gptr            *arralloc(size_mt,int,...); /* Array allocator                */
void            note(char *, ...);      /* Write a message to the output file */
void            message(int *, ...);    /* Write a warning or error message   */
#else
gptr            *arralloc();            /* Array allocator                    */
void            note();                 /* Write a message to the output file */
void            message();              /* Write a warning or error message   */
#endif
int     	nprocessors();          /* Return no. of procs to execute on. */
/*========================== External data references ========================*/
extern  contr_mt control;                   /* Main simulation control parms. */
/*========================== Structs local to module =========================*/
typedef struct cell_s                   /* Prototype element of linked list of*/
{                                       /* molecules within interaction range */
   int          isite, num, frame_type;
   struct cell_s *next;
}               cell_mt;

typedef struct                          /* Prototype of neighbour cell list   */
{                                       /* element.                           */
   int          i, j, k;
}               ivec_mt;

typedef struct                          /* Prototype of neighbour cell list   */
{                                       /* element.                           */
   real         i, j, k;
}               rvec_mt;

typedef struct                          /* Prototype of neighbour cell list   */
{                                       /* element.                           */
   real         x, y, z;
   int          i, j, k;
}               irvec_mt;

/*========================== Global variables ================================*/
static irvec_mt *ifloor; /*Lookup tables for int "floor()"    */
/*========================== Macros ==========================================*/
/*
 * Multiplication factor for size of neighbour list arrays.  If you need
 * to increase this from 1, your system must be *highly* inhomogeneous
 * and may not make sense!
 */
#define         NMULT 3.0
#define         TOO_CLOSE       0.25    /* Error signalled if r**2 < this     */
#define         NCELL(ix,iy,iz) ((iz)+(nz)*((iy)+(ny)*(ix)))
#define         LOCATE(r,eps)   NCELL(cellbin(r[0], nx, eps), \
                                      cellbin(r[1], ny, eps), \
                                      cellbin(r[2], nz, eps))
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
 * spxpy    Sparse add for force.c.  N.B.  MUST NOT BE VECTORIZED as ix may   *
 *          contain duplicate entries.  This occurs if a site interacts with  *
 *          more than one periodic copy of another site.                      *
 ******************************************************************************/
static void spxpy(n, sx, sy, ix)
int     n, ix[];
real    sx[], sy[];
{
   int i;
NOVECTOR
#ifdef __STDC__
#pragma novector
#endif
   for( i = 0; i < n; i++)
   {
      sy[ix[i]] += sx[i];
   }
}

/******************************************************************************
 *  cellbin.    Safe binning function for putting molecules/sites into cells. *
 *  Any error at the boundaries is disasterous and hard to detect.            *
 *  Results may depend on machine-dependant rounding etc.                     *
 ******************************************************************************/
int cellbin(rc, nc, eps)
double rc, eps;
int nc;
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
   if( (ibin = ((rc+0.5)*nc)) >= nc || ibin < 0)
         message(NULLI, NULLP, ERROR, 
                 "Rounding problem in BIN (fill_cells) %.17g\n",rc);
   return ibin;
}

/******************************************************************************
 ******************************************************************************/
#ifdef DEBUG2
#define NBINS 200
static int bins[NBINS];

hist(jmin, jmax, rr)
int jmin, jmax;
real rr[];
{
   double rbin = 10.0;
   int j;

   for(j = jmin; j < jmax; j++)
      if(rr[j] < SQR(NBINS/rbin))
         bins[(int)(rbin*sqrt(rr[j]))]++;
}
void histout()
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
static ivec_mt   *neighbour_list(nnabor, h, cutoff, nx, ny, nz, icheck)
int     *nnabor;
mat_mt  h;
double  cutoff;
int     nx, ny, nz;
int     icheck;
{
   double               dist;
   int                  i, j, ix, iy, iz, mx, my, mz, inabor = 0, nnab;
   static int           onabor=0;
   ivec_mt              *nabor;
   vec_mt               s;
   mat_mt               G, htr, htrinv;

   transpose(h, htr);
   mat_mul(htr, h, G);
   invert(htr, htrinv);

   mx = ceil(cutoff*nx*moda(htrinv));
   my = ceil(cutoff*ny*modb(htrinv));
   mz = ceil(cutoff*nz*modc(htrinv));

   nnab = 4*mx*my*mz;
   nabor = aalloc(nnab, ivec_mt);
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
               if( inabor >= nnab )
                  message(NULLI, NULLP, FATAL,
                          "Internal error in neighbour_list()");
               nabor[inabor].i = ix;
               nabor[inabor].j = iy;
               nabor[inabor].k = iz;
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
static ivec_mt   *strict_neighbour_list(nnabor, h, cutoff, nx, ny, nz, icheck)
int     *nnabor;
mat_mt  h;
double  cutoff;
int     nx, ny, nz;
int     icheck;
{
   double               dist;
   int                  i, j, k, ix, iy, iz, mx, my, mz, inabor = 0, nnab;
   static int           onabor=0;
   int                  ***cellmap;
   ivec_mt              *nabor;
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
   nabor = aalloc(nnab, ivec_mt);
   for(ix = 0; ix <= mx; ix++)
      for(iy = (ix == 0 ? 0 : -my-1); iy <= my; iy++)
         for(iz = (ix == 0 && iy == 0 ? 0 : -mz-1); iz <= mz; iz++)
         {
            if( cellmap[ix][iy][iz] )
            {
               if( inabor >= nnab )
                  message(NULLI, NULLP, FATAL,
                          "Internal error in neighbour_list()");
               nabor[inabor].i = ix;
               nabor[inabor].j = iy;
               nabor[inabor].k = iz;
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
static void    fill_cells(c_of_m, nmols, site, species, h, nx, ny, nz, 
                          lst, cell, frame_type)
vec_mt  c_of_m[];                       /* Centre of mass co-ords        (in) */
int     nmols;                          /* Number of molecules           (in) */
real    **site;                         /* Atomic site co-ordinates      (in) */
spec_mp species;                        /* Pointer to species array      (in) */
mat_mt  h;                              /* Unit cell matrix              (in) */
int     nx, ny, nz;
cell_mt *lst;                           /* Pile of cell structs          (in) */
cell_mt *cell[];                        /* Array of cells (assume zeroed)(out)*/
int     *frame_type;                    /* Framework type counter        (out)*/
{
   int icell, imol, im=0, is, isite = 0;
   double eps = 8.0*precision();
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
            list->frame_type = *frame_type;
            list->next = cell[icell];
            cell[icell] = list++;
            list->next = NULL;
         }
         (*frame_type)++;
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
 *  site_neightbour list.  Build the list of sites withing interaction radius *
 *                         from the lists of sites in cells.                  *
 ******************************************************************************/
int     site_neighbour_list(nab, reloc, n_nab_sites, nfnab, n_frame_types,
                            n_nabors, ix, iy, iz, nx, ny, nz, nabor, cell)
int     *nab;                           /* Array of sites in list      (out) */
rvec_mt reloc[];                        /* Relocation indices for list (out) */
int     n_nab_sites;                    /* Size of above arrays         (in) */
int     nfnab[];                        /* N frame sites index by type (out) */
int     n_frame_types;                  /* Number of distinct frameworks(in) */
int     n_nabors;                       /* Number of neighbour cells    (in) */
int     ix, iy, iz, nx, ny, nz;         /* Labels of current cell       (in) */
ivec_mt *nabor;                         /* List of neighbour cells      (in) */
cell_mt **cell;                         /* Head of cell list            (in) */
{
   int  jx, jy, jz;                     /* Labels of cell in neighbour list  */
   int  j0, jnab, jsite;                /* Counters for cells etc            */
   int  nnab = 0;                       /* Counter for size of nab           */
   int  ftype;
   cell_mt      *cmol;                  /* Pointer to current cell element   */

#ifndef VECTOR
   /*
    * Usual version for scalar machines. 
    */
   irvec_mt *ifl = ifloor;
   real ri, rj, rk;
   int jcell;

   nnab = 0;
   for(ftype = 0; ftype < n_frame_types; ftype++) /* Do Framework types first */
   {
      j0 = (ftype == 0) ? 0: 1;
      for(jnab = j0; jnab < n_nabors; jnab++)    /* Loop over neighbour cells  */
      {
         jx = ix + nabor[jnab].i;       /* jx-jz are indices of neighbour cell*/
         jy = iy + nabor[jnab].j;
         jz = iz + nabor[jnab].k;
         ri  = ifl[jx].x;       /* Compute PBC relocation vector -            */
         jx -= ifl[jx].i;       /* Actually we use a lookup table for speed.  */
         rj  = ifl[jy].y;       /* Ifl[] contains integer and real versions   */
         jy -= ifl[jy].j;       /* of value since conversions are expensive.  */
         rk  = ifl[jz].z;       /* Float value is floor(jx/nx).               */
         jz -= ifl[jz].k;       /* Int value is nx*floor(jx/nx).              */
         jcell = NCELL(jx,jy,jz);
#ifdef DEBUG1
         if(jx<0 || jx>=nx || jy<0 || jy>=ny || jz<0 || jz>=nz)
            message(NULLI,NULLP,FATAL,"Bounds error on reloc (%d,%d,%d)",jx,jy,jz);
#endif
         /* Loop over molecules in this cell, filling 'nab' with its sites    */
         for(cmol = cell[jcell]; cmol != 0; cmol = cmol->next)
         {
            if( cmol->frame_type == ftype )
            {
               for(jsite = 0; jsite < cmol->num; jsite++)
               {
                  nab[nnab] = cmol->isite + jsite;
                  reloc[nnab].i = ri;
                  reloc[nnab].j = rj;
                  reloc[nnab].k = rk;
                  nnab++;
               }
            }
            if(nnab > n_nab_sites) 
               message(NULLI,NULLP,FATAL,TONAB,nnab,n_nab_sites);
         }
      }
      nfnab[ftype] = nnab;
   }
#else
   /* 
    * Version optimized for vector machines, particularly Cray PVP
    */
#ifdef PARALLEL
#ifdef CRAY
#pragma _CRI taskcommon(nnarray, ri,rj, rk, jcell)
#endif
#endif
   int ti,tj,tk;
   static int nnarray;
   static real *ri, *rj, *rk;
   static int  *jcell;

   if( n_nabors > nnarray )
   {
      /*
       * Malloc workspace arrays. Keep them around and only deallocate if
       * required size changes.
       */
      if (ri) 
      {
	 xfree(ri); xfree(jcell);
      }
      ri = dalloc(3*n_nabors);
      rj = ri + n_nabors;
      rk = rj + n_nabors;
      jcell = ialloc(n_nabors);
      nnarray = n_nabors;
   }
   
   nnab = 0;
   for(jnab = 0; jnab < n_nabors; jnab++)    /* Loop over neighbour cells  */
   {
      jx = ix + nabor[jnab].i;       /* jx-jz are indices of neighbour cell*/
      jy = iy + nabor[jnab].j;
      jz = iz + nabor[jnab].k;
      ri[jnab] = ti = IFLOOR(jx,nx);
      rj[jnab] = tj = IFLOOR(jy,ny);
      rk[jnab] = tk = IFLOOR(jz,nz);
      jcell[jnab] = NCELL(jx-ti*nx,jy-tj*ny,jz-tk*nz);
#ifdef DEBUG1
      if(jx<0 || jx>=nx || jy<0 || jy>=ny || jz<0 || jz>=nz)
	 message(NULLI,NULLP,FATAL,"Bounds error on reloc (%d,%d,%d)",jx,jy,jz);
#endif
   }
   for(ftype = 0; ftype < n_frame_types; ftype++) /* Do Framework types first */
   {
      j0 = (ftype == 0) ? 0: 1;
      for(jnab = j0 ; jnab < n_nabors; jnab++)    /* Loop over neighbour cells  */
      {
         /* Loop over molecules in this cell, filling 'nab' with its sites    */
         for(cmol = cell[jcell[jnab]]; cmol != 0; cmol = cmol->next)
         {
            if( cmol->frame_type == ftype)
            {
               for(jsite = 0; jsite < cmol->num; jsite++)
               {
                  nab[nnab] = cmol->isite + jsite;
                  reloc[nnab].i = ri[jnab];
                  reloc[nnab].j = rj[jnab];
                  reloc[nnab].k = rk[jnab];
                  nnab++;
               }
            }
         }
      }
      if(nnab > n_nab_sites) 
	 message(NULLI,NULLP,FATAL,TONAB,nnab,n_nab_sites);
      nfnab[ftype] = nnab;
   }
#endif

   return nnab;
}
/******************************************************************************
 * Force_calc.   This is the main intermolecular site force calculation       *
 * routine                                                                    *
 ******************************************************************************/
void force_calc(site, site_force, system, species, chg, potpar, pe, stress)
real            **site,                 /* Site co-ordinate arrays       (in) */
                **site_force;           /* Site force arrays            (out) */
system_mt       *system;                /* System struct                 (in) */
spec_mt         species[];              /* Array of species records      (in) */
real            chg[];                  /* Array of site charges         (in) */
pot_mt          potpar[];               /* Array of potential parameters (in) */
double          *pe;                    /* Potential energy             (out) */
mat_mt          stress;                 /* Stress virial                (out) */
{
   int          isite, imol, i,         /* Site counter i,j                   */
                i_id, ipot;             /* Miscellaneous                      */
   int          n_frame_types;          /* ==1 for no fw, 2 if fw present.    */
   int          nsites = system->nsites,/* Local copy to keep optimiser happy */
                n_potpar = system->n_potpar,
                max_id = system->max_id;
   double       mol_diam = 2.0*mol_radius(species, system->nspecies),
                cutoff = control.cutoff + (control.strict_cutoff?mol_diam:0);
   double       reloc_lim = MAX(cutoff, control.limit+mol_diam);
   double       subcell = control.subcell; /* Local copy. May change it.      */
   int          *id   = ialloc(nsites), /* Array of site_id[nsites]           */
                *id_ptr;                /* Pointer to 'id' array              */
   int          n_nab_sites = nsites*   /* Dimension of site n'bor list arrays*/
#ifdef DEBUG2 
                                MAX(1.0,NMULT*4.19*CUBE(cutoff)/det(system->h));
#else
                                NMULT*4.19*CUBE(cutoff)/det(system->h);
#endif
   ivec_mt      *nabor, *rdf_nabor;     /* Lists of neighbour cells           */
   int          n_nabors, n_rdf_nabors; /* Number of elements in lists.       */
   int          icell, ncells;          /* Subcell counter and total = nxnynz */
   int          n_cell_list;            /* Size of link-cells "heap"          */
   int          nx, ny, nz, nmax;       /* Number of subcells in MD cell      */
   static int   onx=0, ony=0, onz=0;    /* Saved values of nx, ny, nz.        */
   static int   mmax;                   /* Saved offset of ifloor for free(). */
   int          mx, my, mz;             /* Limits for ifloor array.           */
   real         ***potp                 /* Expanded potential parameter array */
                      = (real***)arralloc((size_mt)sizeof(real), 3, P0,max_id-1,
                                          0, n_potpar-1, 0, nsites-1);
   cell_mt      *c_ptr;                 /* Heap of link cell entries for list */
   cell_mt      **cell;                 /* Array of list heads for subcells   */
   spec_mt      *spec;                  /* Temp. loop pointer to species.     */
   mat_mt       htr, htrinv;            /* Transpos and inverse of h matrix   */
   int		nthreads = nprocessors(),
   		ithread;
   /*
    * Thread-local force and stress arrays For compiler-parallel version
    */
   int		j;
   real		*sf0, *sf1, *sf2, *ssf0, *ssf1, *ssf2;
   double	*pe_n = (double *)aalloc(nthreads, double);
   mat_mt	*stress_n = (mat_mt *)aalloc(nthreads, mat_mt);
   real		***s_f_n;
   
#ifdef DEBUG2
   double ppe, rr[3], ss[3];
   int im, is, i;
   mat_mp h = system->h;
   mat_mt hinv;
#endif

   /*
    * Initialise thread-local arrays
    */
   s_f_n = aalloc(nthreads, real**);
   s_f_n[0] = site_force;
   for(ithread = 1; ithread < nthreads; ithread++)
   {
      s_f_n[ithread] = (real**)arralloc((size_mt)sizeof(real),2,0,2,0,nsites-1);
      zero_real(s_f_n[ithread][0],3*nsites);
   }
   zero_real(stress_n[0][0],9*nthreads);
   zero_double(pe_n, nthreads);

   /*
    * Choose a partition into subcells if none specified.
    */
   if(subcell <= 0.0) subcell = control.cutoff/5.0;
   nx = system->h[0][0]/subcell+0.5;
   ny = system->h[1][1]/subcell+0.5;
   nz = system->h[2][2]/subcell+0.5;
   ncells = nx*ny*nz;
   if( nx != onx || ny != ony || nz != onz )
   {
      note("MD cell divided into %d subcells (%dx%dx%d)",ncells,nx,ny,nz);
      onx = nx; ony = ny; onz = nz;
      /*
       * Allocate and fill lookup tables for floor(jx/nx) etc.
       */
      if( ifloor )
         xfree(ifloor-mmax);
      
      transpose(system->h, htr);
      invert(htr, htrinv);
      mx = (int)ceil(reloc_lim*nx*moda(htrinv))+1;
      my = (int)ceil(reloc_lim*ny*modb(htrinv))+1;
      mz = (int)ceil(reloc_lim*nz*modc(htrinv))+1;
      mmax = MAX3(mx, my, mz);
      nmax = MAX3(nx, ny, nz);

      ifloor = aalloc(2*mmax+nmax, irvec_mt)+mmax;
      for(i = -mmax; i < mmax+nmax; i++)
      {
         ifloor[i].x = IFLOOR(i,nx);
         ifloor[i].y = IFLOOR(i,ny);
         ifloor[i].z = IFLOOR(i,nz);
         ifloor[i].i = nx*IFLOOR(i,nx);
         ifloor[i].j = ny*IFLOOR(i,ny);
         ifloor[i].k = nz*IFLOOR(i,nz);
      }
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
   n_cell_list = 1;
   for(spec = species; spec < species+system->nspecies; spec++)
      if( spec->framework )
         n_cell_list += spec->nmols*spec->nsites;
      else 
         n_cell_list += spec->nmols;
   c_ptr = aalloc(n_cell_list, cell_mt); 
   /*
    * Build a linked list of molecules/sites for each subcell.
    * "cell" is the array of list heads NX x NY x NZ.
    */
   cell = aalloc(ncells, cell_mt *);     
#pragma _CRI novector
   for( icell=0; icell < ncells; icell++)
      cell[icell] = NULL;
   fill_cells(system->c_of_m, system->nmols, site, species, system->h,
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
   
   
#ifdef DEBUG2
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
                norm,control.alpha,system->ptype,potp[id[is]]);
      }
      isite += spec->nsites;
   }
   histout();
   note("Direct pot. energy = %g",ppe*CONV_E);
#endif

#ifdef PARALLEL
#ifdef stellar
/*$dir parallel*/
#endif /* stellar */
#ifdef titan
#pragma ipdep
#endif /* titan */
#ifdef __convexc__
#pragma _CNX force_parallel
#endif /* --convexc__ */
#ifdef CRAY
#pragma _CRI taskloop private(ithread) shared(nthreads,site,chg, potp, id,\
		  n_nab_sites, n_nabors, nabor, nx, ny, nz, cell, n_frame_types, \
		  system, stress_n, pe_n, s_f_n)
#endif /* CRAY */
#endif /*PARALLEL */
   for(ithread = 0; ithread < nthreads; ithread++)
      force_inner(ithread, nthreads, site, chg, potp, id, n_nab_sites, 
		  n_nabors, nabor, nx, ny, nz, cell, n_frame_types, system,
		  stress_n[ithread], pe_n+ithread, s_f_n[ithread]);
/*
 *  Sum Pot, energies, forces and stress from each parallel invocation
 */
   for(ithread = 0; ithread < nthreads; ithread++)
   {
      *pe += pe_n[ithread];
NOVECTOR
      for(i = 0; i < 3; i++)
	 for(j = 0; j < 3; j++)
	    stress[i][j] += stress_n[ithread][i][j];
   }
   /*
    * Sum thread's copies of the site forces.  s_f_n[ithread] points
    * to the force arrays for each thread.  s_f_n[0] is just
    * site_force and the others are independant arrays.
    */
   sf0 = site_force[0]; sf1 = site_force[1]; sf2 = site_force[2];
   for(ithread = 1; ithread < nthreads; ithread++)
   {
      ssf0 = s_f_n[ithread][0];
      ssf1 = s_f_n[ithread][1];
      ssf2 = s_f_n[ithread][2];
VECTORIZE
#ifdef titan
#pragma ipdep
#endif
#ifdef __convexc__
#pragma _CNX vstrip (32)
#pragma _CNX force_vector
#pragma _CNX force_parallel_ext
#endif
#ifdef CRAY
#pragma _CRI ivdep
#endif
      for(isite = 0; isite < nsites; isite++)
      {
	 sf0[isite] += ssf0[isite];
	 sf1[isite] += ssf1[isite];
	 sf2[isite] += ssf2[isite];
      }
   }
   /*
    * Accumulate radial distribution functions
    */
   if (control.rdf_interval > 0 && 
       control.istep >= control.begin_rdf &&
       control.istep % control.rdf_interval == 0)
   {
      n_nab_sites = nsites*
                    NMULT*4.19*CUBE(control.limit+mol_diam)/det(system->h);
      rdf_nabor = strict_neighbour_list(&n_rdf_nabors, system->h, 
                                        control.limit+mol_diam, nx, ny, nz, 0);
      rdf_inner(0, 1, site, id, n_nab_sites, n_rdf_nabors, 
                rdf_nabor, nx, ny, nz, cell, n_frame_types, system);
      xfree(rdf_nabor);
   }
#ifdef DEBUG2
   histout();
#endif
   afree((gptr*)(potp+P0));  xfree(c_ptr); 
   xfree(cell);        xfree(id); 
   xfree(nabor);
   xfree(pe_n);   xfree(stress_n);
   for( ithread = 1; ithread < nthreads; ithread++)
      afree((gptr*)s_f_n[ithread]);
   xfree(s_f_n);
}
#ifdef titan
#ifdef PARALLEL
#pragma opt_level 2
#endif
#endif
/******************************************************************************
 *  Force_inner() Paralellised inner loops of force_calc.  Loops over cells   *
 *  in MD cell with stride = nomber of processors available.  Should be       *
 *  called once for each parallel thread.                                     *
 ******************************************************************************/
void
force_inner(ithread, nthreads, site, chg, potp, id, n_nab_sites, n_nabors, 
            nabor, nx, ny, nz, cell, n_frame_types, system,
            stress, pe, site_force)
int             ithread, nthreads;      /* Parallel node variables.      (in) */
real            **site;                 /* Site co-ordinate arrays       (in) */
real            chg[];                  /* Array of site charges         (in) */
real            ***potp;                /* Expanded potential parameter array */
int             id[];                   /* Array of site_id[nsites]      (in) */
int             n_nab_sites;            /* Dimension of site n'bor list arrays*/
int             n_nabors;               /* Number of elements in lists.   (in)*/
ivec_mt         *nabor;                 /* Lists of neighbour cells       (in)*/
int             nx, ny, nz;             /* Number of subcells in MD cell  (in)*/
cell_mt         **cell;                 /* Array of list heads of subcells(in)*/
int             n_frame_types;          /* ==1 for no fw, 2 if fw present (in)*/
system_mt       *system;                /* System struct                 (in) */
mat_mt          stress;                 /* Stress virial                (out) */
double          *pe;                    /* Potential energy             (out) */
real            **site_force;           /* Site force arrays            (out) */
{
                /*
                 * The following arrays are for 'neighbour site list'
                 * quantities and should be dimensioned to the max value of
                 * 'nnab'.  A rough approx is the ratio of the volume of
                 * the "cutoff sphere" to that of the MD cell times nsites.
                 * This may be too small for inhomogeneous systems, but at
                 * least it scales with the cutoff radius.
                 */
   int          *nab  = ialloc(n_nab_sites);    /* Neigbour site gather vector*/
   rvec_mt      *reloc = aalloc(n_nab_sites, rvec_mt); /* Site PBC shifts     */
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
   real         force_cpt, site0, site1, site2, s00, s01, s02, s11, s12, s22;
                                   /* Accumulators for forces and stresses.   */
   real         rrx, rry, rrz;                  /* Scalar loop temporaries    */
   real         h00, h01, h02, h11, h12, h22;   /* Temp copies of system->h   */
   double       norm = 2.0*control.alpha/sqrt(PI);      /* Coulombic prefactor*/
   double       cutoffsq = SQR(control.cutoff), /* Temporary copy for optim'n */
                cutoff100sq = 10000.0*cutoffsq;
   int          ix, iy, iz;             /* 3-d cell indices for ref and neig. */
   int          icell,                  /* Index for cells of molecule pair   */
                nnab, jbeg, jmin, jmax, /* Number of sites in neighbour list  */
                isite, jsite, ipot, lim;/* Counters.                          */
   int          nsites = system -> nsites;      /* Temporary copy for optim'n */
   int          nfnab[2];               /* Number of non-fw and fw neighbours */
   cell_mt      *cmol;                  /* Loop counter for link cells.       */

   s00 = s01 = s02 = s11 = s12 = s22 = 0.0;     /* Accumulators for stress    */

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
      nnab = 0;
#ifdef DEBUG1
      printf("Working on cell %4d (%d,%d,%d) (sites %4d to %4d)\n", icell,
             ix,iy,iz,cell[icell]->isite,cell[icell]->isite+cell[icell]->num-1);
      printf("\n jcell\tjx jy jz\tNsites\n");
#endif
      /*
       * Build site neighbour list 'nab' from cell list.
       */ 
      nnab = site_neighbour_list(nab, reloc, n_nab_sites, nfnab, n_frame_types, 
                                 n_nabors, ix, iy, iz, nx, ny, nz, nabor, cell);
#ifdef DEBUG4
      for(jsite=0; jsite<nnab; jsite++)
         printf("%d %d %d %d\n",nab[jsite], 
                reloc[jsite].i,reloc[jsite].j,reloc[jsite].j);
#endif
      gather(nnab, nab_sx, site[0], nab, nsites); /* Construct list of site  */
      gather(nnab, nab_sy, site[1], nab, nsites); /* co-ordinates from nabor */
      gather(nnab, nab_sz, site[2], nab, nsites); /* list.                   */

      /*
       * Apply periodic boundary conditions to neighbour site co-ords.
       * Assume h matrix is upper triangular.
       */
      h00 = system->h[0][0];   h01 = system->h[0][1];   h02 = system->h[0][2];
      h11 = system->h[1][1];   h12 = system->h[1][2];   h22 = system->h[2][2];
VECTORIZE
      for(jsite=0; jsite<nnab; jsite++)
      {
         rrx = nab_sx[jsite] 
                 + h00*reloc[jsite].i + h01*reloc[jsite].j + h02*reloc[jsite].k;
         rry = nab_sy[jsite] +          h11*reloc[jsite].j + h12*reloc[jsite].k;
         rrz = nab_sz[jsite] +                               h22*reloc[jsite].k;
         nab_sx[jsite] = rrx;
         nab_sy[jsite] = rry;
         nab_sz[jsite] = rrz;
      }
#ifdef DEBUG7
      for(jsite = 0; jsite < nnab; jsite++)
         printf("%f %f %f\n",nab_sx[jsite],nab_sy[jsite],nab_sz[jsite]);
#endif
      gather(nnab, nab_chg, chg, nab, nsites); /* Gather site charges as well*/
      zero_real(forcejx,nnab);
      zero_real(forcejy,nnab);
      zero_real(forcejz,nnab);

      jbeg = 0;                 /* Extra element alllocated makes [1] safe*/
                                /* Loop over all molecules in cell icell. */
      for(cmol = cell[icell]; cmol != NULL; cmol = cmol->next)
      {
         if( cmol->frame_type )
         {
            jmin = 0;
            jmax = nfnab[0];
         }
         else
         {
            jmin = jbeg += cmol->num;
            jmax = nnab;
         }
         lim = cmol->isite + cmol->num;
         for(isite = cmol->isite; isite < lim; isite++)
         {                                   /* Loop over sites in molecule */
            /*
             * Construct pot'l param arrays corresponding to neighbour sites.
             */
            pp = potp[id[isite]];
            ppp = nab_pot;
            for(ipot = 0; ipot < system->n_potpar; ipot++)
            {
               gather(jmax, *ppp++, *pp++, nab, nsites);
            }
#ifdef DEBUG1
            if(isite == 100)
#endif
#if defined(DEBUG1) || defined(DEBUG5)
            { int jnab;
            for(jnab = jmin; jnab < jmax; jnab++)
               printf("%4d %4d\n", jnab,nab[jnab]);
           }
#endif
            site0=site[0][isite]; site1=site[1][isite]; site2=site[2][isite];
VECTORIZE
            for(jsite=jmin; jsite < jmax; jsite++)
            {
               rrx = nab_sx[jsite] - site0;
               rry = nab_sy[jsite] - site1;
               rrz = nab_sz[jsite] - site2;
               r_sqr[jsite] = rrx*rrx+rry*rry+rrz*rrz;
               rx[jsite] = rrx;
               ry[jsite] = rry;
               rz[jsite] = rrz;
            }
            if( (jsite = jmin+search_lt(jmax-jmin, r_sqr+jmin, 1, TOO_CLOSE))
               < jmax )
               message(NULLI, NULLP, WARNING, TOOCLS,
                       isite, nab[jsite], sqrt(TOO_CLOSE));

            if( control.strict_cutoff )
               for(jsite = jmin; jsite < jmax; jsite++)
                  if( r_sqr[jsite] > cutoffsq )
                     r_sqr[jsite] = cutoff100sq;

#ifdef DEBUG2
            hist(jmin, jmax, r_sqr);
#endif
               
            /*  Call the potential function kernel                            */
            kernel(jmin, jmax, forceij, pe, r_sqr, nab_chg, chg[isite],
                   norm, control.alpha, system->ptype, nab_pot);
            site0 = site1 = site2 = 0.0;
VECTORIZE
            for(jsite=jmin; jsite < jmax; jsite++)
            {
               force_cpt  = forceij[jsite]*rx[jsite];
               s00         += force_cpt * rx[jsite];
               s02         += force_cpt * rz[jsite];
               s01         += force_cpt * ry[jsite];
               site0           -= force_cpt;
               forcejx[jsite]  += force_cpt;
            }
VECTORIZE
            for(jsite=jmin; jsite < jmax; jsite++)
            {
               force_cpt = forceij[jsite]*ry[jsite];
               s11         += force_cpt * ry[jsite];
               s12         += force_cpt * rz[jsite];
               site1           -= force_cpt;
               forcejy[jsite]  += force_cpt;
            }
VECTORIZE
            for(jsite=jmin; jsite < jmax; jsite++)
            {
               force_cpt = forceij[jsite]*rz[jsite];
               s22         += force_cpt * rz[jsite];
               site2           -= force_cpt;
               forcejz[jsite]  += force_cpt;
            }
            site_force[0][isite] += site0;
            site_force[1][isite] += site1;
            site_force[2][isite] += site2;
#ifdef DEBUG3
            printf("PE = %f\n",pe[0]);
#endif
         }
      }
      spxpy(nnab, forcejx, site_force[0], nab);
      spxpy(nnab, forcejy, site_force[1], nab);
      spxpy(nnab, forcejz, site_force[2], nab);
   }
   stress[0][0]  += s00;
   stress[0][1]  += s01;
   stress[0][2]  += s02;
   stress[1][1]  += s11;
   stress[1][2]  += s12;
   stress[2][2]  += s22;

   afree((gptr*)nab_pot);
   xfree(nab);     xfree(reloc);  xfree(nab_chg);
   xfree(r_sqr);   xfree(forceij);
   xfree(rx);      xfree(ry);      xfree(rz);
   xfree(forcejx); xfree(forcejy); xfree(forcejz);
   xfree(nab_sx);  xfree(nab_sy);  xfree(nab_sz);
}
/******************************************************************************
 *  Rdf_inner() Paralellised inner loops of force_calc.  Based on force_inner *
 *     but only calls rdf_accum().                                            *
 ******************************************************************************/
void
rdf_inner(ithread, nthreads, site, id, n_nab_sites, n_nabors, 
            nabor, nx, ny, nz, cell, n_frame_types, system)
int             ithread, nthreads;      /* Parallel node variables.      (in) */
real            **site;                 /* Site co-ordinate arrays       (in) */
int             id[];                   /* Array of site_id[nsites]      (in) */
int             n_nab_sites;            /* Dimension of site n'bor list arrays*/
int             n_nabors;               /* Number of elements in lists.   (in)*/
ivec_mt         *nabor;                 /* Lists of neighbour cells       (in)*/
int             nx, ny, nz;             /* Number of subcells in MD cell  (in)*/
cell_mt         **cell;                 /* Array of list heads of subcells(in)*/
int             n_frame_types;          /* ==1 for no fw, 2 if fw present (in)*/
system_mt       *system;                /* System struct                  (in)*/
{
   int          *nab  = ialloc(n_nab_sites);    /* Neigbour site gather vector*/
   rvec_mt      *reloc = aalloc(n_nab_sites, rvec_mt); /* Site PBC shifts     */
   real         *nab_sx  = dalloc(n_nab_sites), /* 'Gathered' list of         */
                *nab_sy  = dalloc(n_nab_sites), /*   neighbour site co-ords   */
                *nab_sz  = dalloc(n_nab_sites), /*   - x,y,z components.      */
                *r_sqr   = dalloc(n_nab_sites); /* Squared site-site distance */
   real         site0, site1, site2;
   real         rrx, rry, rrz;                  /* Scalar loop temporaries    */
   real         h00, h01, h02, h11, h12, h22;   /* Temp copies of system-> h  */
   int          ix, iy, iz;             /* 3-d cell indices for ref and neig. */
   int          icell,                  /* Index for cells of molecule pair   */
                nnab, jbeg, jmin, jmax, /* Number of sites in neighbour list  */
                isite, jsite, lim;      /* Counters.                          */
   int          nsites = system -> nsites;      /* Temporary copy for optim'n */
   int          nfnab[2];               /* Number of non-fw and fw neighbours */
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
      nnab = site_neighbour_list(nab, reloc, n_nab_sites, nfnab, n_frame_types, 
                                 n_nabors, ix, iy, iz, nx, ny, nz, nabor, cell);
      gather(nnab, nab_sx, site[0], nab, nsites); /* Construct list of site  */
      gather(nnab, nab_sy, site[1], nab, nsites); /* co-ordinates from nabo  */
      gather(nnab, nab_sz, site[2], nab, nsites); /* list.                   */

      /*
       * Apply periodic boundary conditions to neighbour site co-ords.
       * Assume h matrix is upper triangular.
       */
      h00 = system->h[0][0];   h01 = system->h[0][1];   h02 = system->h[0][2];
      h11 = system->h[1][1];   h12 = system->h[1][2];   h22 = system->h[2][2];
VECTORIZE
      for(jsite=0; jsite<nnab; jsite++)
      {
         rrx = nab_sx[jsite] 
                 + h00*reloc[jsite].i + h01*reloc[jsite].j + h02*reloc[jsite].k;
         rry = nab_sy[jsite] +          h11*reloc[jsite].j + h12*reloc[jsite].k;
         rrz = nab_sz[jsite] +                               h22*reloc[jsite].k;
         nab_sx[jsite] = rrx;
         nab_sy[jsite] = rry;
         nab_sz[jsite] = rrz;
      }

      jbeg = 0;                 /* Extra element alllocated makes [1] safe*/
                                /* Loop over all molecules in cell icell. */
      for(cmol = cell[icell]; cmol != NULL; cmol = cmol->next)
      {
         if( cmol->frame_type )
         {
            jmin = 0;
            jmax = nfnab[0];
         }
         else
         {
            jmin = jbeg += cmol->num;
            jmax = nnab;
         }
         lim = cmol->isite + cmol->num;
         for(isite = cmol->isite; isite < lim; isite++)
         {                                   /* Loop over sites in molecule */
            site0=site[0][isite]; site1=site[1][isite]; site2=site[2][isite];
VECTORIZE
            for(jsite=jmin; jsite < jmax; jsite++)
            {
               rrx = nab_sx[jsite] - site0;
               rry = nab_sy[jsite] - site1;
               rrz = nab_sz[jsite] - site2;
               r_sqr[jsite] = rrx*rrx+rry*rry+rrz*rrz;
            }
            /*
             * Accumulate radial distribution functions
             */
            rdf_accum(jmin, jmax, r_sqr, id[isite], id, nab);
         }
      }
   }
   xfree(nab);     xfree(reloc);
   xfree(r_sqr);
   xfree(nab_sx);  xfree(nab_sy);  xfree(nab_sz);
}
