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
 *       $Log:	force.c,v $
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
static char *RCSid = "$Header: /home/tigger/keith/md/RCS/force.c,v 1.8.1.8 89/11/01 17:34:15 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#ifdef  convexvc
#include 	<fastmath.h>
#else
#include 	<math.h>
#endif
#include 	"string.h"
/*========================== Program include files ===========================*/
#include "structs.h"
#include "messages.h"
/*========================== Library declarations ============================*/
void    cfree();
char	*realloc();
/*========================== External function declarations ==================*/
void    note();                         /* Make a note in output file         */
void    *arralloc();                    /* General purpose array allocator    */
int     search_lt();			/* Search a vector for el. < scalar   */
double  vdot();                         /* Vector dot product                 */
double  sum();                          /* Sum a vector                       */
void    scatter();                      /* Interface to CRAY scatter routine  */
void    gather();                       /* Interface to CRAY gather routine   */
void	gatheri();			/* Integer gather		      */
int	wheneq();			/* Array indexer		      */
void    mat_mul();                      /* Matrix multiplier                  */
void    message();                      /* Send message to stderr             */
void	invert();			/* 3x3 matrix inverter		      */
void	mat_vec_mul();			/* Matrix by vector multiplier        */
void	spaxpy();			/* Scattered vector add		      */
void	transpose();			/* Generate 3x3 matrix transpose      */
void    zero_real();                    /* Initialiser                        */
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
   short 	i, j, k;
}		ivec_t;

typedef struct				/* Prototype for a 3Nx x 3Ny x 3Nz    */
{					/* pbc relocation array.	      */
   int 		ncell, rel;
}		reloc_t;
/*========================== Global variables ================================*/
static          int nx = 0, ny = 0, nz = 0;
/*========================== Macros ==========================================*/
#define         NCELL(ix,iy,iz) ((iz)+nz*((iy)+ny*(ix)))
#define		NSH 1
#define		NSHELL (2*NSH+1)
#define         NREL(ix,iy,iz) ((iz)+NSH+NSHELL*((iy)+NSH+NSHELL*((ix)+NSH)))
#define         CELLMAX 5
#define         TOO_CLOSE       0.25    /* Error signalled if r**2 < this     */
#define		BIN(rc, nc)	((int)((rc+2.5)*nc)%nc)
#define		LOCATE(r)	NCELL(BIN(r[0], nx), BIN(r[1],ny), BIN(r[2],nz))
#define moda(hmat) (hmat[0][0])
#define modb(hmat) sqrt(SQR(hmat[0][1]) + SQR(hmat[1][1]))
#define modc(hmat) sqrt(SQR(hmat[0][2]) + SQR(hmat[1][2]) + SQR(hmat[2][2]))
/*============================================================================*/
/******************************************************************************
 * As of version 4.0, the CRAY C compiler can not vectorise the 'kernel'      *
 * function.  A replacement in FORTRAN is provided which will vectorise. To   *
 * use it, compile this module with -DFKERNEL to define the symbol FKERNEL.   *
 ******************************************************************************/
#ifdef FKERNEL
#ifdef unix
#define KERNEL kernel_
#endif
void    KERNEL();                       /* Force kernel routine               */
#else
void    kernel();                       /* Force kernel routine               */
#endif
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
#ifdef DEBUG
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
#ifdef DEBUG
               printf("%12f %4d %4d %4d %12f %12f %12f\n",
                      dist,ix,iy,iz,s[0],s[1],s[2]);
#endif
            }
         }
   			/* Return unused storage to heap		*/
   if( (nabor=(ivec_t *)realloc((char*)nabor, (size_t)inabor*sizeof(ivec_t)))
      == NULL)
      message(NULLI,NULLP, FATAL,"Realloc fails in neighbour_list()");
   note(NABORS,2 * inabor);
   *nnabor = inabor;
   return(nabor);
}
/******************************************************************************
 *  Fill_cells.  Allocate all the sites to cells depending on their centre of *
 *  mass co-ordinate by binning.                                              *
 ******************************************************************************/
static void    fill_cells(c_of_m, nmols, site, spec, h, list, cell, frame_type)
vec_t   c_of_m[];                       /* Centre of mass co-ords        (in) */
int     nmols;                          /* Number of molecules           (in) */
real	**site;				/* Atomic site co-ordinates      (in) */
spec_p  spec;                           /* Pointer to species array      (in) */
mat_t	h;				/* Unit cell matrix              (in) */
cell_t  *list;				/* Pile of cell structs          (in) */
cell_t  *cell[];                        /* Array of cells (assume zeroed)(out)*/
int	*frame_type;			/* Framework type counter	 (out)*/
{
   int icell, imol, im=0, is, isite = 0;
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
	    icell = LOCATE(ssite);
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
	 icell = LOCATE(c_of_m[imol]);
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
			    n_nabors, ix, iy, iz, nabor, cell, reloc )
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
{
   int	jx, jy, jz;			/* Labels of cell in neighbour list  */
   int	jcell, jnab, jsite, rl;		/* Counters for cells etc	     */
   int	nnab = 0;			/* Counter for size of nab	     */
   int	ftype, c, nnabf = 0;
   cell_t	*cmol;			/* Pointer to current cell element   */
   int	*reloc_if = ialloc(4*n_nab_sites),
        *frame_key = reloc_if + n_nab_sites,
   	*nabf = frame_key + n_nab_sites,
        *idx = nabf + n_nab_sites;
   
   for(jnab = 0; jnab < n_nabors; jnab++)    /* Loop over neighbour cells  */
   {
      jx = ix + nabor[jnab].i;
      jy = iy + nabor[jnab].j;
      jz = iz + nabor[jnab].k;
#ifdef DEBUG
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
	 if(nnabf > n_nab_sites) message(NULLI,NULLP,FATAL,TONAB,nnabf);
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

   cfree((char*)reloc_if);
   
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
/******************************************************************************
 * Vadd vector addition							      *
 ******************************************************************************/
static void vadd(n, a, b)
int	n;
real	a[], b[];
{
   int i;
VECTORIZE
   for(i=0; i<n; i++)
      a[i] += b[i];
}
/******************************************************************************
 *  vaadd,calcdist. These are functions because of a compiler bug in cray C4.0*
 ******************************************************************************/
#ifdef VCALLS
static vaadd(nnab, forcej, force_comp, forceij, r)
int nnab;
real forcej[], force_comp[], forceij[],r[];
{
   int jsite;
VECTORIZE
   for(jsite=0; jsite < nnab; jsite++)
      forcej[jsite] += force_comp[jsite] = forceij[jsite]*r[jsite];
}
#endif
static void calcdist(j0, nnab, r_sqr, rx, ry, rz, nab_sx, nab_sy, nab_sz,
		     sitex, sitey, sitez)
int j0, nnab;
real r_sqr[], rx[], ry[], rz[], nab_sx[], nab_sy[], nab_sz[];
double sitex, sitey, sitez;
{
   int jsite;
VECTORIZE
            for(jsite=j0; jsite < nnab; jsite++)
            {
               rx[jsite] = nab_sx[jsite] - sitex;
               ry[jsite] = nab_sy[jsite] - sitey;
               rz[jsite] = nab_sz[jsite] - sitez;
               r_sqr[jsite] = rx[jsite]*rx[jsite] + ry[jsite]*ry[jsite]
                                                  + rz[jsite]*rz[jsite];
            }
}
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
   int          *id_ptr;                /* Pointer to 'id' array              */
   		/*
		 * The following arrays are for 'neighbour site list'
		 * quantities and should be dimensioned to the max value of
		 * 'nnab'. Nsites is certainly too big, but as 'ewald' uses
		 * more store, it shouldn't increase total program use.
		 */
   int		n_nab_sites = nsites;		/* Max # sites in n'bor list  */
   int		*id = ialloc(n_nab_sites);   	/* Array of site_id[nsites]   */
   real         ***potp				/* Expanded pot'l parameters  */
   		= (real***)arralloc(sizeof(real), 3,
                                    1, max_id-1, 0, n_potpar-1, 0, nsites-1);
   static	int n_cell_list = 1;
   cell_t       *c_ptr = aalloc(n_cell_list, cell_t );
   spec_p       spec;
   cell_t       **cell;
   
   static boolean       init = true;
   static int           n_nabors;
   static ivec_t        *nabor;
   static reloc_t       ***reloc;
   int		nthreads = nprocessors(),
   		ithread;
   double	*pe_n = (double *)aalloc(nthreads, double);
   mat_t	*stress_n = (mat_t *)aalloc(nthreads, mat_t);
   real		***s_f_n;

#ifdef DEBUG
   int jnab;
#endif
#ifdef DEBUG2
   double ppe, rrx, rry, rrz;
   int im, is;
   mat_p h = system->h;
#endif
  
   if (nthreads > 1)
       s_f_n = (real***)arralloc(sizeof(real), 3,
				 0, nthreads-2, 0, 2, 0, nsites-1);

   if(init)
   {
      if(control.subcell <= 0.0) control.subcell = control.cutoff/5.0;
      nx = system->h[0][0]/control.subcell+0.5;
      ny = system->h[1][1]/control.subcell+0.5;
      nz = system->h[2][2]/control.subcell+0.5;
      ncells = nx*ny*nz;
      note("MD cell divided into %d subcells (%dx%dx%d)",ncells,nx,ny,nz);

      if(control.cutoff > NSH*MIN3(system->h[0][0],system->h[1][1],system->h[2][2]))
	 message(NULLI, NULLP, FATAL,
		 "Cutoff radius > cell dimension not supported at present");
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
      for(spec = species; spec < species+system->nspecies; spec++)
	 if( spec->framework )
	    n_cell_list += spec->nmols*spec->nsites;
	 else 
	    n_cell_list += spec->nmols;
      (void)cfree((char*)c_ptr);
      c_ptr = aalloc(n_cell_list, cell_t);
      init = false;
   }

   cell = aalloc(ncells, cell_t *);

/*  Construct and fill expanded site-identifier array, id                     */
   id_ptr = id;
   for (spec = species; spec < species+system->nspecies; spec++)
      for(imol = 0; imol < spec->nmols; imol++)
      {
         (void)memcpy((char*)id_ptr, (char*)spec->site_id, 
                      (int)spec->nsites*sizeof(int));
         id_ptr += spec->nsites;
      }
/*   Build arrays of pot. pars [max_id][nsites] for access in vector loops    */
   for(ipot = 0; ipot < n_potpar; ipot++)
      for(i_id = 1; i_id < max_id; i_id++)
VECTORIZE
         for(isite = 0; isite < nsites; isite++)
            potp[i_id][ipot][isite] = potpar[i_id*max_id+id[isite]].p[ipot];

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
/*$dir parallel*/
   for(ithread = 0; ithread < nthreads; ithread++)
      force_inner(ithread, nthreads, site, chg, potp, id,
		  n_nab_sites, n_nabors, nabor, cell, reloc, n_frame_types, 
		  system,
		  stress_n[ithread], pe_n+ithread,
		  ithread ? s_f_n[ithread-1] : site_force);
/*
 *  Sum Pot, energies, forces and stress from each parallel invocation
 */
   for(ithread = 0; ithread < nthreads; ithread++)
   {
      *pe += pe_n[ithread];
      for(i = 0; i < 3; i++)
	 for(j = 0; j < 3; j++)
	    stress[i][j] += stress_n[ithread][i][j];
   }
   for(ithread = 0; ithread < nthreads-1; ithread++)
VECTORIZE
      for(isite = 0; isite < nsites; isite++)
      {
	 site_force[0][isite] += s_f_n[ithread][0][isite];
	 site_force[1][isite] += s_f_n[ithread][1][isite];
	 site_force[2][isite] += s_f_n[ithread][2][isite];
      }

#ifdef DEBUG2
   histout();
#endif
   cfree((char*)(potp+1));  cfree((char*)c_ptr); 
   cfree((char*)cell);        cfree((char*)id); 
   cfree((char*)pe_n);   cfree((char*)stress_n);
   if( nthreads > 1)
      cfree((char*)s_f_n);
}
/******************************************************************************
 *  Force_inner() Paralellised inner loops of force_calc.  Loops over cells   *
 *  in MD cell with stride = nomber of processors available.  Should be       *
 *  called once for each parallel thread.				      *
 ******************************************************************************/
force_inner(ithread, nthreads, site, chg, potp, id, n_nab_sites, n_nabors, 
	    nabor, cell, reloc, n_frame_types, system,
	    stress, pe, site_force)
int	ithread, nthreads;
real	**site, **site_force;
real	chg[];
real	***potp;
int	id[];
int	n_nab_sites;
system_t *system;
cell_t	**cell;
reloc_t	***reloc;
ivec_t	*nabor;
real	*pe;
mat_t	stress;
{
   int          *nab  = ialloc(n_nab_sites),	/* Neigbour site gather vector*/
                *reloc_i = ialloc(n_nab_sites);	/* Vector of pbc relocations  */
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
#ifdef VCALLS
   real         *force_comp = dalloc(n_nab_sites);
#endif
   real         **nab_pot			/* Gathere'd pot par array    */
   		= (real**)arralloc(sizeof(real), 2,
				   0, system->n_potpar-1, 0, n_nab_sites-1);
   real         force_cpt, site0, site1, site2, s00, s01, s02, s11, s12, s22;
   real         reloc_v[3][CUBE(NSHELL)];	/* PBC relocation vectors     */
   short        ix, iy, iz;		/* 3-d cell indices for ref and neig. */
   int          icell,			/* Index for cells of molecule pair   */
                nnab,	j0, jmin, jmax,	/* Number of sites in neighbour list  */
   		isite, jsite, ipot, lim;
   int		nfnab[2];
   cell_t       *cmol;
   double       norm = 2.0*control.alpha/sqrt(PI);	/* Coulombic prefactor*/
   s00 = s01 = s02 = s11 = s12 = s22 = 0.0;	/* Accumulators for stress    */
   reloc_alloc(system->h, reloc_v);

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
#ifdef DEBUG
      printf("Working on cell %4d (%d,%d,%d) (sites %4d to %4d)\n", icell,
             ix,iy,iz,cell[icell]->isite,cell[icell]->isite+cell[icell]->num-1);
      printf("\n jcell\tjx jy jz\tNsites\n");
#endif
      /*
       * Build site neighbour list 'nab' from cell list.
       */ 
      nnab = site_neighbour_list(nab, reloc_i,n_nab_sites,nfnab,n_frame_types, 
				 n_nabors, ix, iy, iz, nabor, cell, reloc);

      gather(nnab, nab_sx, site[0], nab);     /* Construct list of site     */
      gather(nnab, nab_sy, site[1], nab);     /* co-ordinated from nabor    */
      gather(nnab, nab_sz, site[2], nab);     /* list.                      */

      gather(nnab, R, reloc_v[0], reloc_i);
      vadd(nnab, nab_sx, R);
      gather(nnab, R, reloc_v[1], reloc_i);
      vadd(nnab, nab_sy, R);
      gather(nnab, R, reloc_v[2], reloc_i);
      vadd(nnab, nab_sz, R);

      gather(nnab, nab_chg, chg, nab);       /* Gather site charges as well*/
      zero_real(forcejx,nnab);
      zero_real(forcejy,nnab);
      zero_real(forcejz,nnab);

      j0 = 0;			/* Extra element alllocated makes [1] safe*/
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
	    jmin = j0 += cmol->num;
	    jmax = nnab;
	 }
         lim = cmol->isite + cmol->num;
         for(isite = cmol->isite; isite < lim; isite++)
         {                                   /* Loop over sites in molecule */
            for(ipot = 0; ipot < system->n_potpar; ipot++)
               gather(jmax, nab_pot[ipot], potp[id[isite]][ipot], nab);
#ifdef DEBUG
            if(isite == 100)
            for(jsite = jmin; jsite < jmax; jsite++)
               printf("%4d %4d\n", jsite,nab[jsite]);
#endif

            calcdist(jmin, jmax, r_sqr, rx, ry, rz,nab_sx,nab_sy,nab_sz,
		     site[0][isite],site[1][isite],site[2][isite]);
#ifdef DEBUG2
	    hist(jmin, jmax, r_sqr);
#endif
            if( (jsite = jmin+search_lt(jmax-jmin, r_sqr+jmin, 1, TOO_CLOSE))
	       < jmax )
               message(NULLI, NULLP, WARNING, TOOCLS,
		       isite, nab[jsite], sqrt(TOO_CLOSE));

            /*  Call the potential function kernel                            */
#ifdef FKERNEL
            KERNEL(&jmin, &jmax, forceij, pe, r_sqr, nab_chg, chg+isite, &norm,
		   &control.alpha, &system->ptype, &n_nab_sites, nab_pot[0]);
#else
            kernel(jmin, jmax, forceij, pe, r_sqr, nab_chg, chg[isite],
		   norm, control.alpha, system->ptype, nab_pot);
#endif
#ifdef VCALLS
            vaadd(jmax-jmin,forcejx+jmin,force_comp+jmin,forceij+jmin,rx+jmin);
            site_force[0][isite] -= sum(jmax-jmin, force_comp+jmin,1);
            s00                  += vdot(jmax-jmin,force_comp+jmin,1,rx+jmin,1);
            s01                  += vdot(jmax-jmin,force_comp+jmin,1,ry+jmin,1);
            s02                  += vdot(jmax-jmin,force_comp+jmin,1,rz+jmin,1);

            vaadd(jmax-jmin,forcejy+jmin,force_comp+jmin,forceij+jmin,ry+jmin);
            site_force[1][isite] -= sum(jmax-jmin, force_comp+jmin,1);
            s11                  += vdot(jmax-jmin,force_comp+jmin,1,ry+jmin,1);
            s12                  += vdot(jmax-jmin,force_comp+jmin,1,rz+jmin,1);

            vaadd(jmax-jmin,forcejz+jmin,force_comp+jmin,forceij+jmin,rz+jmin);
            site_force[2][isite] -= sum(jmax-jmin, force_comp+jmin,1);
            s22                  += vdot(jmax-jmin,force_comp+jmin,1,rz+jmin,1);
#else
            site0 = site1 = site2 = 0.0;
VECTORIZE
            for(jsite=jmin; jsite < jmax; jsite++)
            {
               force_cpt  = forceij[jsite]*rx[jsite];
               site0           -= force_cpt;
               forcejx[jsite]  += force_cpt;
               s00         += force_cpt * rx[jsite];
               s02         += force_cpt * rz[jsite];
               s01         += force_cpt * ry[jsite];
            }
VECTORIZE
            for(jsite=jmin; jsite < jmax; jsite++)
            {
               force_cpt = forceij[jsite]*ry[jsite];
               site1           -= force_cpt;
               forcejy[jsite]  += force_cpt;
               s11         += force_cpt * ry[jsite];
               s12         += force_cpt * rz[jsite];
            }
VECTORIZE
            for(jsite=jmin; jsite < jmax; jsite++)
            {
               force_cpt = forceij[jsite]*rz[jsite];
               site2           -= force_cpt;
               forcejz[jsite]  += force_cpt;
               s22         += force_cpt * rz[jsite];
            }
            site_force[0][isite] += site0;
            site_force[1][isite] += site1;
            site_force[2][isite] += site2;
#endif
         }
      }
      spaxpy(nnab, 1.0, forcejx, site_force[0], nab);
      spaxpy(nnab, 1.0, forcejy, site_force[1], nab);
      spaxpy(nnab, 1.0, forcejz, site_force[2], nab);
   }
   stress[0][0]  += s00;
   stress[0][1]  += s01;
   stress[0][2]  += s02;
   stress[1][1]  += s11;
   stress[1][2]  += s12;
   stress[2][2]  += s22;

   cfree((char*)nab_pot);   
   cfree((char*)nab);     cfree((char*)reloc_i);  cfree((char*)nab_chg);
   cfree((char*)r_sqr);   cfree((char*)R);       cfree((char*)forceij);
   cfree((char*)rx);      cfree((char*)ry);      cfree((char*)rz);
   cfree((char*)forcejx); cfree((char*)forcejy); cfree((char*)forcejz);
   cfree((char*)nab_sx);  cfree((char*)nab_sy);  cfree((char*)nab_sz);
#ifdef VCALLS
   cfree((char*)force_comp);
#endif
}
