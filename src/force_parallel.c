/******************************************************************************
 * Force	This module contains functions to implement the 'link cell'   *
 *		interatomic force calculation (Hockney, R.W. & Eastwood J.W.  *
 *		"Computer Simulation Using Particles" McGraw-Hill (1981), 277)*
 *		It is vectorised and optimised for a CRAY XMP and Convex C1,  *
 *		but should run well on any vector machine with reasonably     *
 *		efficient scatter/gather. (And of course on scalar machines!) *
 *		The actual calculation of the potential is in a different     *
 *		module (kernel.c) for ease of modification.		      *
 ******************************************************************************
 *      Revision Log
 * Revision 1.3  90/03/29  15:44:51  keith
 * Merged force.c revisions 1.8.1.11-1.8.1.13
 * 
 * Revision 1.2  90/03/09  17:30:29  keith
 * Modified FKERNEL ifdefs for UNICOS.
 * 
 * Revision 1.1  90/01/31  13:19:28  keith
 * Initial revision
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
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/force_parallel.c,v 1.20 91/11/26 12:47:25 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#ifdef  convexvc
#include 	<fastmath.h>
#else
#include 	<math.h>
#endif
#include	"stddef.h"
#include 	"string.h"
/*========================== Program include files ===========================*/
#include 	"structs.h"
#include 	"messages.h"
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
void    	note();                 /* Make a note in output file         */
gptr    	*arralloc();            /* General purpose array allocator    */
int     	search_lt();            /* Search a vector for el. < scalar   */
double  	vdot();                 /* Vector dot product                 */
double  	sum();                  /* Sum a vector                       */
void    	gather();               /* Interface to CRAY gather routine   */
void    	gatheri();              /* Integer gather                     */
int     	wheneq();               /* Array indexer                      */
void    	mat_mul();              /* Matrix multiplier                  */
double		det(); 			/* Determinant of 3x3 matrix	      */
void    	message();              /* Send message to stderr             */
void    	invert();               /* 3x3 matrix inverter                */
void    	mat_vec_mul();          /* Matrix by vector multiplier        */
void    	spaxpy();               /* Scattered vector add               */
void    	transpose();            /* Generate 3x3 matrix transpose      */
void    	zero_real();            /* Initialiser                        */
void    	force_inner();          /* Inner loop forward reference       */
int     	nprocessors();          /* Return no. of procs to execute on. */
double  	precision();            /* Floating pt precision.             */
void    	kernel();               /* Force kernel routine               */
/*========================== External data references ========================*/
extern  contr_t control;
/*========================== Structs local to module =========================*/
typedef struct cell_s			/* Prototype element of linked list of*/
{					/* molecules within interaction range */
   int 		isite, num, frame_type;
   struct cell_s *next;
}		cell_t;

typedef struct				/* Prototype of neighbour cell list   */
{					/* element.			      */
   int	 	i, j, k;
}		ivec_t;

typedef struct				/* Prototype for a 3Nx x 3Ny x 3Nz    */
{					/* pbc relocation array.	      */
   int 		ncell, rel;
}		reloc_t;
/*========================== Global variables ================================*/
static          int nx = 0, ny = 0, nz = 0;
/*========================== Macros ==========================================*/
#define         NCELL(ix,iy,iz) ((iz)+nz*((iy)+ny*(ix)))
/*
 * Maximum cutoff radius relative to MD cell dimension.
 */
#define		NSH 1
/*
 * Multiplication factor for size of neighbour list arrays.  If you need
 * to increase this from 1, your system must be *highly* inhomogeneous
 * and may not make sense!
 */
#define         NMULT 1.0
#define		NSHELL (2*NSH+1)
#define         NREL(ix,iy,iz) ((iz)+NSH+NSHELL*((iy)+NSH+NSHELL*((ix)+NSH)))
#define         CELLMAX 5
#define         TOO_CLOSE       0.25    /* Error signalled if r**2 < this     */
#define		LOCATE(r,eps)	NCELL(cellbin(r[0], nx, eps), \
				      cellbin(r[1], ny, eps), \
				      cellbin(r[2], nz, eps))
#define moda(hmat) (hmat[0][0])
#define modb(hmat) sqrt(SQR(hmat[0][1]) + SQR(hmat[1][1]))
#define modc(hmat) sqrt(SQR(hmat[0][2]) + SQR(hmat[1][2]) + SQR(hmat[2][2]))
/*============================================================================*/
/******************************************************************************
 *  cellbin.    Safe binning function for putting molecules/sites into cells. *
 *  Any error at the boundaries is disasterous and hard to detect.            *
 *  Results may depend on machine-dependant rounding etc.		      *
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
static ivec_t    *neighbour_list(nnabor, h, cutoff)
int  	*nnabor;
mat_t   h;
double  cutoff;
{
   double               dist;
   int                  i, j, ix, iy, iz, mx, my, mz, inabor = 0, nnab;
   ivec_t		*nabor;
   vec_t                s;
   mat_t                G, htr, htrinv;

   transpose(h, htr);
   mat_mul(htr, h, G);
   invert(htr, htrinv);

   mx = ceil(cutoff*nx*moda(htrinv));
   my = ceil(cutoff*ny*modb(htrinv));
   mz = ceil(cutoff*nz*modc(htrinv));

   nnab = 4*mx*my*mz;
   nabor = aalloc(nnab, ivec_t);
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
   note(NABORS,2 * inabor);
   *nnabor = inabor;
   return(nabor);
}
/******************************************************************************
 *  Fill_cells.  Allocate all the sites to cells depending on their centre of *
 *  mass co-ordinate by binning.                                              *
 ******************************************************************************/
static void    fill_cells(c_of_m, nmols, site, species, h, lst, cell, frame_type)
vec_t   c_of_m[];                       /* Centre of mass co-ords        (in) */
int     nmols;                          /* Number of molecules           (in) */
real	**site;				/* Atomic site co-ordinates      (in) */
spec_p  species;                        /* Pointer to species array      (in) */
mat_t	h;				/* Unit cell matrix              (in) */
cell_t  *lst;				/* Pile of cell structs          (in) */
cell_t  *cell[];                        /* Array of cells (assume zeroed)(out)*/
int	*frame_type;			/* Framework type counter	 (out)*/
{
   int icell, imol, im=0, is, isite = 0;
   double eps = 8.0*precision();
   spec_p spec = species;
   cell_t *list = lst;
   vec_t ssite;
   mat_t hinv;

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
	    mat_vec_mul(hinv, (vec_t*)ssite, (vec_t*)ssite, 1);
	    icell = LOCATE(ssite, eps);
	    list->isite = isite++;
	    list->num   = 1;
	    list->frame_type = *frame_type;
	    list->next = cell[icell];
	    cell[icell] = list++;
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
	 isite += spec->nsites;
      }
   }
}
/******************************************************************************
 *  site_neightbour list.  Build the list of sites withing interaction radius *
 *                         from the lists of sites in cells.		      *
 ******************************************************************************/
int	site_neighbour_list(nab, reloc_i, n_nab_sites, nfnab, n_frame_types,
			    n_nabors, ix, iy, iz, nabor, cell, reloc, work )
int	*nab;				/* Array of sites in list      (out) */
int	*reloc_i;			/* Relocation indices for list (out) */
int	n_nab_sites;			/* Size of above arrays		(in) */
int	nfnab[];			/* N frame sites index by type (out) */
int	n_frame_types;			/* Number of distinct frameworks(in) */
int	n_nabors;			/* Number of neighbour cells    (in) */
int	ix, iy, iz;			/* Labels of current cell	(in) */
ivec_t	*nabor;				/* List of neighbour cells	(in) */
cell_t	**cell;				/* Head of cell list		(in) */
reloc_t	***reloc;			/* Relocation index array	(in) */
int   	*work;                          /* Workspace                         */
{
   int	jx, jy, jz;			/* Labels of cell in neighbour list  */
   int	jcell, jnab, jsite, rl;		/* Counters for cells etc	     */
   int	nnab = 0;			/* Counter for size of nab	     */
   int	ftype, c, nnabf = 0;
   cell_t	*cmol;			/* Pointer to current cell element   */
   int	*reloc_if = work,
        *frame_key = reloc_if + n_nab_sites,
   	*nabf = frame_key + n_nab_sites,
        *idx = nabf + n_nab_sites;
   
   for(jnab = 0; jnab < n_nabors; jnab++)    /* Loop over neighbour cells  */
   {
      jx = ix + nabor[jnab].i;
      jy = iy + nabor[jnab].j;
      jz = iz + nabor[jnab].k;
#ifdef DEBUG1
      if(jx<-nx || jx>=2*nx || jy<-ny || jy>=2*ny || jz<-nz || jz>=2*nz)
	 message(NULLI,NULLP,FATAL,"Bounds error on reloc (%d,%d,%d)",jx,jy,jz);
#endif
      jcell = reloc[jx][jy][jz].ncell;
      rl    = reloc[jx][jy][jz].rel;

      /* Loop over molecules in this cell, filling 'nab' with its sites    */
      for(cmol = cell[jcell]; cmol != NULL; cmol = cmol->next)
      {
         if( jnab > 0 || cmol->frame_type == 0 )
	    for(jsite = 0; jsite < cmol->num; jsite++)
	    {
	       nabf[nnabf] = cmol->isite + jsite;
	       reloc_if[nnabf] = rl;
	       frame_key[nnabf] = cmol->frame_type;
	       nnabf++;
	    }
	 if(nnabf > n_nab_sites) 
 	    message(NULLI,NULLP,FATAL,TONAB,nnabf,n_nab_sites);
      }
   }

   nnab = 0;
   for( ftype = 0; ftype < n_frame_types; ftype++ )
   {
      c = wheneq(nnabf, frame_key, 1, ftype, idx);
      gatheri(c, nab+nnab, nabf, idx);
      gatheri(c, reloc_i+nnab, reloc_if, idx);

      nfnab[ftype] = nnab += c;
   }

   return nnab;
}
/******************************************************************************
 *  reloc_alloc.  Allocate and fill relocation array                          *
 ******************************************************************************/
static void    reloc_alloc(h, reloc)
mat_t   h;
real    reloc[][CUBE(NSHELL)];
{
   int  tx, ty, tz;
   for(tx = -NSH; tx <= NSH; tx++)
      for(ty = -NSH; ty <= NSH; ty++)
         for(tz = -NSH; tz <= NSH; tz++)
         {
            reloc[0][NREL(tx,ty,tz)] = tx*h[0][0] + ty*h[0][1] + tz*h[0][2];
            reloc[1][NREL(tx,ty,tz)] = ty*h[1][1] + tz*h[1][2];
            reloc[2][NREL(tx,ty,tz)] = tz*h[2][2];
         }
}
#ifdef titan
#ifdef PARALLEL
#pragma opt_level 3
#pragma pproc force_inner
#endif
#endif
/******************************************************************************
 * Force_calc.   This is the main intermolecular site force calculation       *
 * routine                                                                    *
 ******************************************************************************/
void force_calc(site, site_force, system, species, chg, potpar, pe, stress)
real            **site,                 /* Site co-ordinate arrays       (in) */
                **site_force;           /* Site force arrays            (out) */
system_p        system;                 /* System struct                 (in) */
spec_t          species[];              /* Array of species records      (in) */
real            chg[];                  /* Array of site charges         (in) */
pot_t           potpar[];               /* Array of potential parameters (in) */
double          *pe;                    /* Potential energy             (out) */
mat_t           stress;                 /* Stress virial                (out) */
{
   int		isite, imol,            /* Site counter i,j                   */
   		i_id, ipot, i, j;	/* Miscellaneous		      */
   int		n_frame_types;
   int          nsites = system->nsites,/* Local copy to keep optimiser happy */
                n_potpar = system->n_potpar,
   		max_id = system->max_id;
   int          ncells = nx*ny*nz;	/* Total number of cells	      */
   int          tx, ty, tz;             /* Temporaries for # unit cell shifts */
   int		ix, iy, iz;
   int		*id      = ialloc(nsites),   	/* Array of site_id[nsites]   */
   		*id_ptr;                /* Pointer to 'id' array              */
   real         ***potp				/* Expanded pot'l parameters  */
   		= (real***)arralloc(sizeof(real), 3,
                                    1, max_id-1, 0, n_potpar-1, 0, nsites-1);
   		/*
		 * The following arrays are for 'neighbour site list'
		 * quantities and should be dimensioned to the max value of
		 * 'nnab'.  A rough approx is the ratio of the volume of
		 * the "cutoff sphere" to that of the MD cell times nsites.
		 * This may be too small for inhomogeneous systems, but at
		 * least it scales with the cutoff radius.
		 */
   int		n_nab_sites = nsites*		/* Max # sites in n'bor list  */
                    NMULT*4.19*CUBE(control.cutoff)/det(system->h);
   static	int n_cell_list = 1;
   cell_t       *c_ptr = aalloc(n_cell_list, cell_t );
   spec_p       spec;
   cell_t       **cell;
   real		*sf0, *sf1, *sf2, *ssf0, *ssf1, *ssf2;
   
   static boolean       init = true;
   static int           n_nabors;
   static ivec_t        *nabor;
   static reloc_t       ***reloc;
   int		nthreads = nprocessors(),
   		ithread;
   double	*pe_n = (double *)aalloc(nthreads, double);
   mat_t	*stress_n = (mat_t *)aalloc(nthreads, mat_t);
   real		***s_f_n;

#ifdef DEBUG2
   double ppe, rrx, rry, rrz;
   int im, is;
   mat_p h = system->h;
#endif
  
   s_f_n = aalloc(nthreads, real**);
   s_f_n[0] = site_force;
   for(ithread = 1; ithread < nthreads; ithread++)
      s_f_n[ithread] = (real**)arralloc(sizeof(real), 2, 0, 2, 0, nsites-1);

   if(init)
   {
      if(control.subcell <= 0.0) control.subcell = control.cutoff/5.0;
      nx = system->h[0][0]/control.subcell+0.5;
      ny = system->h[1][1]/control.subcell+0.5;
      nz = system->h[2][2]/control.subcell+0.5;
      ncells = nx*ny*nz;
      note("MD cell divided into %d subcells (%dx%dx%d)",ncells,nx,ny,nz);

      if(control.cutoff > NSH*MIN3(system->h[0][0],system->h[1][1],system->h[2][2]))
	 message(NULLI, NULLP, FATAL, CUTOFF, NSH);
      nabor = neighbour_list(&n_nabors, system->h, control.cutoff);
   
      reloc = (reloc_t***)arralloc(sizeof(reloc_t),3,
            -NSH*nx, (NSH+1)*nx-1, -NSH*ny, (NSH+1)*ny-1, -NSH*nz, (NSH+1)*nz-1);
 
      for(ix = 0; ix < nx; ix++)
         for(tx = -NSH*nx; tx <= NSH*nx; tx += nx)
            for(iy = 0; iy < ny; iy++)
               for(ty = -NSH*ny; ty <= NSH*ny; ty += ny)
                  for(iz = 0; iz < nz; iz++)
                     for(tz = -NSH*nz; tz <= NSH*nz; tz += nz)
                     {
                        reloc[ix+tx][iy+ty][iz+tz].ncell = NCELL(ix, iy, iz);
                        reloc[ix+tx][iy+ty][iz+tz].rel= NREL(tx/nx,ty/ny,tz/nz);
                     }
#ifdef titan
#pragma asis
#endif
      for(spec = species; spec < species+system->nspecies; spec++)
	 if( spec->framework )
	    n_cell_list += spec->nmols*spec->nsites;
	 else 
	    n_cell_list += spec->nmols;
      xfree(c_ptr);
      c_ptr = aalloc(n_cell_list, cell_t);
      init = false;
   }

   cell = aalloc(ncells, cell_t *);

/*  Construct and fill expanded site-identifier array, id                     */
   id_ptr = id;
   for (spec = species; spec < species+system->nspecies; spec++)
      for(imol = 0; imol < spec->nmols; imol++)
      {
         memcp(id_ptr, spec->site_id, spec->nsites*sizeof(int));
         id_ptr += spec->nsites;
      }
/*   Build arrays of pot. pars [max_id][nsites] for access in vector loops    */
   for(ipot = 0; ipot < n_potpar; ipot++)
      for(i_id = 1; i_id < max_id; i_id++)
      {
#ifdef titan
NOVECTOR
#endif
         for(isite = 0; isite < nsites; isite++)
            potp[i_id][ipot][isite] = potpar[i_id*max_id+id[isite]].p[ipot];
      }

   fill_cells(system->c_of_m, system->nmols, site, species, system->h,
	      c_ptr, cell, &n_frame_types);
   if( n_frame_types > 2 )
      message(NULLI, NULLP, FATAL,
	      "Multiple framework molecules are not supported");
   

#ifdef DEBUG2
   ppe = 0;
   spec = species; isite = 0;
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
	    rrx = site[0][jsite] - site[0][is];
	    rry = site[1][jsite] - site[1][is];
	    rrz = site[2][jsite] - site[2][is];
	    rrx -= h[0][0]*floor(rrx/h[0][0]+0.5);
	    rry -= h[1][1]*floor(rry/h[1][1]+0.5);
	    rrz -= h[2][2]*floor(rrz/h[2][2]+0.5);
	    r_sqr[jsite] = rrx*rrx+rry*rry+rrz*rrz;
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
/*
 *   Start of main loop over processors
 */
#ifdef stellar
/*$dir parallel*/
#endif
#ifdef titan
#pragma ipdep
#endif
#ifdef __convexc__
#pragma _CNX force_parallel
#endif
   for(ithread = 0; ithread < nthreads; ithread++)
      force_inner(ithread, nthreads, site, chg, potp, id,
		  n_nab_sites, n_nabors, nabor, cell, reloc, n_frame_types, 
		  system,
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
      for(isite = 0; isite < nsites; isite++)
      {
	 sf0[isite] += ssf0[isite];
	 sf1[isite] += ssf1[isite];
	 sf2[isite] += ssf2[isite];
      }
   }

#ifdef DEBUG2
   histout();
#endif
   xfree((potp+1));  xfree(c_ptr); 
   xfree(cell);        xfree(id); 
   xfree(pe_n);   xfree(stress_n);
   for( ithread = 1; ithread < nthreads; ithread++)
      xfree(s_f_n[ithread]);
}
#ifdef titan
#ifdef PARALLEL
#pragma opt_level 2
#endif
#endif
/******************************************************************************
 *  Force_inner() Paralellised inner loops of force_calc.  Loops over cells   *
 *  in MD cell with stride = nomber of processors available.  Should be       *
 *  called once for each parallel thread.				      *
 ******************************************************************************/
void
force_inner(ithread, nthreads, site, chg, potp, id, n_nab_sites, n_nabors, 
	    nabor, cell, reloc, n_frame_types, system,
	    stress, pe, site_force)
int	ithread, nthreads;
real	**site, **site_force;
real	chg[];
real	***potp;
int	id[];
int	n_nab_sites;
int	n_nabors;
int	n_frame_types;
system_t *system;
cell_t	**cell;
reloc_t	***reloc;
ivec_t	*nabor;
real	*pe;
mat_t	stress;
{
   int          *nab  = ialloc(n_nab_sites),	/* Neigbour site gather vector*/
                *reloc_i = ialloc(n_nab_sites),	/* Vector of pbc relocations  */
   		*work = ialloc(4*n_nab_sites);  /* Workspace for s_n_list     */
   real         *nab_sx  = dalloc(n_nab_sites),	/* 'Gathered' list of         */
   		*nab_sy  = dalloc(n_nab_sites),	/*   neighbour site co-ords   */
                *nab_sz  = dalloc(n_nab_sites),	/*   - x,y,z components.      */
                *forcejx = dalloc(n_nab_sites),	/* List of neighbour site     */
                *forcejy = dalloc(n_nab_sites),	/*  forces in gathered form   */
                *forcejz = dalloc(n_nab_sites),	/*  - xyz components.	      */
                *rx      = dalloc(n_nab_sites),	/* Reference to neigbour site */
                *ry      = dalloc(n_nab_sites),	/* - site vector adjusted for */
                *rz      = dalloc(n_nab_sites),	/*  periodic boundaries. xyz. */
                *r_sqr   = dalloc(n_nab_sites),	/* Squared site-site distance */
                *nab_chg = dalloc(n_nab_sites),	/* Gathered neig. site charges*/
                *forceij = dalloc(n_nab_sites),	/* -V'(r) / r		      */
                *R = dalloc(n_nab_sites);    	/* pbc site relocation cpt    */
   real         **nab_pot			/* Gathere'd pot par array    */
   		= (real**)arralloc(sizeof(real), 2,
				   0, system->n_potpar-1, 0, n_nab_sites-1);
   real         force_cpt, site0, site1, site2, s00, s01, s02, s11, s12, s22;
   real         reloc_v[3][CUBE(NSHELL)];	/* PBC relocation vectors     */
   real 	**pp, **ppp;
   int          ix, iy, iz;		/* 3-d cell indices for ref and neig. */
   int          icell,			/* Index for cells of molecule pair   */
                nnab, jbeg, jmin, jmax,	/* Number of sites in neighbour list  */
   		isite, jsite, ipot, lim;
   int		nsites = system -> nsites;
   int		nfnab[2];
   cell_t       *cmol;
   double       norm = 2.0*control.alpha/sqrt(PI);	/* Coulombic prefactor*/
   s00 = s01 = s02 = s11 = s12 = s22 = 0.0;	/* Accumulators for stress    */
   reloc_alloc(system->h, reloc_v);
#ifdef DEBUG6
   { int i;
     for(i = 0; i < CUBE(NSHELL); i++)
	printf("%f %f %f\n", reloc_v[0][i], reloc_v[1][i], reloc_v[2][i]);
  }
#endif

/******************************************************************************
 *  Start of main loops.  Loop over species, ispec, molecules, imol, and      *
 *  sites, isite_mol on imol.  Isite is index into 'site' of site specified   *
 *  by isite_mol, imol and ispec.  Jbase is set to first site of imol+1, ispec*
 *  so jsite, counts from jbase to end.  Thus isite, jsite run over all site  *
 *  site pairs on distinct molecules.                                         *
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
      nnab = site_neighbour_list(nab, reloc_i,n_nab_sites,nfnab,n_frame_types, 
				 n_nabors, ix, iy, iz, nabor, cell, reloc, work);
#ifdef DEBUG4
      for(jsite=0; jsite<nnab; jsite++)
	 printf("%d %d\n",nab[jsite], reloc_i[jsite]);
#endif
      gather(nnab, nab_sx, site[0], nab, nsites);     /* Construct list of site     */
      gather(nnab, nab_sy, site[1], nab, nsites);     /* co-ordinated from nabor    */
      gather(nnab, nab_sz, site[2], nab, nsites);     /* list.                      */

      gather(nnab, R, reloc_v[0], reloc_i, CUBE(NSHELL));
VECTORIZE
      for(jsite=0; jsite<nnab; jsite++)
	 nab_sx[jsite] += R[jsite];
      gather(nnab, R, reloc_v[1], reloc_i, CUBE(NSHELL));
VECTORIZE
      for(jsite=0; jsite<nnab; jsite++)
	 nab_sy[jsite] += R[jsite];
      gather(nnab, R, reloc_v[2], reloc_i, CUBE(NSHELL));
VECTORIZE
      for(jsite=0; jsite<nnab; jsite++)
	 nab_sz[jsite] += R[jsite];
#ifdef DEBUG7
      for(jsite = 0; jsite < nnab; jsite++)
	 printf("%f %f %f\n",nab_sx[jsite],nab_sy[jsite],nab_sz[jsite]);
#endif
      gather(nnab, nab_chg, chg, nab, nsites);       /* Gather site charges as well*/
      zero_real(forcejx,nnab);
      zero_real(forcejy,nnab);
      zero_real(forcejz,nnab);

      jbeg = 0;			/* Extra element alllocated makes [1] safe*/
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
               rx[jsite] = nab_sx[jsite] - site0;
               ry[jsite] = nab_sy[jsite] - site1;
               rz[jsite] = nab_sz[jsite] - site2;
               r_sqr[jsite] = rx[jsite]*rx[jsite] + ry[jsite]*ry[jsite]
                                                  + rz[jsite]*rz[jsite];
            }
#ifdef DEBUG2
	    hist(jmin, jmax, r_sqr);
#endif
            if( (jsite = jmin+search_lt(jmax-jmin, r_sqr+jmin, 1, TOO_CLOSE))
	       < jmax )
               message(NULLI, NULLP, WARNING, TOOCLS,
		       isite, nab[jsite], sqrt(TOO_CLOSE));

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

   xfree(nab_pot); xfree(work);  
   xfree(nab);     xfree(reloc_i);  xfree(nab_chg);
   xfree(r_sqr);   xfree(R);       xfree(forceij);
   xfree(rx);      xfree(ry);      xfree(rz);
   xfree(forcejx); xfree(forcejy); xfree(forcejz);
   xfree(nab_sx);  xfree(nab_sy);  xfree(nab_sz);
}
