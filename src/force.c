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
 * Revision 1.1  89/04/20  16:00:40  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: force.c,v 1.1.1.1 89/05/17 12:45:50 keith Exp $";
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
/*========================== External function declarations ==================*/
void    note();                         /* Make a note in output file         */
void    *arralloc();                    /* General purpose array allocator    */
int     search_lt();			/* Search a vector for el. < scalar   */
void    vadd();                         /* Vector add                         */
double  vdot();                         /* Vector dot product                 */
double  sum();                          /* Sum a vector                       */
void    scatter();                      /* Interface to CRAY scatter routine  */
void    gather();                       /* Interface to CRAY gather routine   */
void    mat_mul();                      /* Matrix multiplier                  */
void    message();                      /* Send message to stderr             */
void    zero_real();                    /* Initialiser                        */
/*========================== External data references ========================*/
extern  contr_t control;
/*========================== Structs local to module =========================*/
typedef struct cell_s			/* Prototype element of linked list of*/
{					/* molecules within interaction range */
   int 		isite, num;
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
#define         NREL(ix,iy,iz) ((iz)+1+3*((iy)+1+3*((ix)+1)))
#define         CELLMAX 5
#define         TOO_CLOSE       0.25    /* Error signalled if r**2 < this     */
/*============================================================================*/
/******************************************************************************
 * As of version 4.0, the CRAY C compiler can not vectorise the 'kernel'      *
 * function.  A replacement in FORTRAN is provided which will vectorise. To   *
 * use it, compile this module with -DFKERNEL to define the symbol FKERNEL.   *
 ******************************************************************************/
#ifdef FKERNEL
#ifndef CRAY
#define KERNEL kernel_
#endif
void    KERNEL();                       /* Force kernel routine               */
#else
void    kernel();                       /* Force kernel routine               */
#endif
/******************************************************************************
 *  Neighbour_list.  Build the list of cells within cutoff radius of cell 0   *
 ******************************************************************************/
static void    neighbour_list(nabor, n_nabors, h, cutoff)
ivec_t  nabor[];
int     *n_nabors;
mat_t   h;
double  cutoff;
{
   double               dist;
   int                  i, j, ix, iy, iz;
   vec_t                s;
   mat_t                G;
   
   mat_mul(h, h, G);
#ifdef DEBUG
   printf("  Distance    ix    iy    iz      sx        sy          sz\n");
#endif
   for(ix = nx/2; ix < nx+1; ix++)
      for(iy = (ix==nx/2 ? ny/2 : 0); iy < ny+1; iy++)
         for(iz = (ix==nx/2 && iy==ny/2 ? nz/2 : 0); iz < nz+1; iz++)
         {
            s[0] = (ix + 0.5*(nx%2))/nx - 0.5;
            s[1] = (iy + 0.5*(ny%2))/ny - 0.5;
            s[2] = (iz + 0.5*(nz%2))/nz - 0.5;
            dist = 0.0;
            for(i = 0; i < 3; i++)
               for(j = 0; j < 3; j++)
                  dist += s[i]*G[i][j]*s[j];
            if(dist < SQR(cutoff))
            {
               nabor[*n_nabors].i = ix - nx / 2;
               nabor[*n_nabors].j = iy - ny / 2;
               nabor[*n_nabors].k = iz - nz / 2;
               (*n_nabors)++;
#ifdef DEBUG
               printf("%12f %4d %4d %4d %12f %12f %12f\n",
                      dist,ix,iy,iz,s[0],s[1],s[2]);
#endif
            }
         }
   note(NABORS,2 * *n_nabors);
}
/******************************************************************************
 *  Fill_cells.  Allocate all the sites to cells depending on their centre of *
 *  mass co-ordinate by binning.                                              *
 ******************************************************************************/
static void    fill_cells(c_of_m, nmols, spec, cell, mol)
vec_t   c_of_m[];                       /* Centre of mass co-ords        (in) */
int     nmols;                          /* Number of molecules           (in) */
spec_p  spec;                           /* Pointer to species array      (in) */
cell_t  *cell[];                        /* Array of cells (assume zeroed)(out)*/
cell_t  mol[];                          /* What 'cell' points to         (out)*/
{
   int ix, iy, iz, icell, imol, isite = 0;
#ifdef DEBUG
   printf("  ix    iy    iz  icell isite\n");
#endif
   for(imol = 0; imol < nmols; imol++)
   {
      if(imol == spec->nmols) spec++;

      mol[imol].isite = isite;
      mol[imol].num  = spec->nsites;

      ix = (int)((c_of_m[imol][0] + 0.5)*nx)%nx;/* Get 'cell' co-ords ix,iy,iz*/
      iy = (int)((c_of_m[imol][1] + 0.5)*ny)%ny;/* in [0,n*-1) from co-ords   */
      iz = (int)((c_of_m[imol][2] + 0.5)*nz)%nz;/* c_of_m in [-0.5,0.5]       */
      icell = NCELL(ix, iy, iz);                /* 1-d cell index             */
      mol[imol].next = cell[icell];             /* Link in 'mol[imol]'        */
      cell[icell]   = &(mol[imol]);             /* at head of list.           */
#ifdef DEBUG
      printf(" %5d %5d %5d %5d %5d\n",ix, iy, iz, icell, isite);
#endif
      isite += spec->nsites;
   }
}
/******************************************************************************
 *  site_neightbour list.  Build the list of sites withing interaction radius *
 *                         from the lists of sites in cells.		      *
 ******************************************************************************/
int	site_neighbour_list(nab, reloc_i,
			    n_nabors, ix, iy, iz, nabor, cell, reloc )
int	*nab;				/* Array of sites in list      (out) */
int	*reloc_i;			/* Relocation indices for list (out) */
int	n_nabors;			/* Number of neighbour cells    (in) */
int	ix, iy, iz;			/* Labels of current cell	(in) */
ivec_t	*nabor;				/* List of neighbour cells	(in) */
cell_t	**cell;				/* Head of cell list		(in) */
reloc_t	***reloc;			/* Relocation index array	(in) */
{
   int	jx, jy, jz;			/* Labels of cell in neighbour list  */
   int	jcell, jnab, jsite, rl;		/* Counters for cells etc	     */
   int	nnab = 0;			/* Counter for size of nab	     */
   cell_t	*cmol;			/* Pointer to current cell element   */
   for(jnab = 0; jnab < n_nabors; jnab++)    /* Loop over neighbour cells  */
   {
      jx = ix + nabor[jnab].i;
      jy = iy + nabor[jnab].j;
      jz = iz + nabor[jnab].k;
      jcell = reloc[jx][jy][jz].ncell;
      rl    = reloc[jx][jy][jz].rel;

      /* Loop over molecules in this cell, filling 'nab' with its sites    */
      for(cmol = cell[jcell]; cmol != NULL; cmol = cmol->next)
      {
	 for(jsite = 0; jsite < cmol->num; jsite++)
	 {
	    nab[nnab] = cmol->isite + jsite;
	    reloc_i[nnab] = rl;
	    nnab++;
	 }
      }
   }
   return nnab;
}
/******************************************************************************
 *  reloc_alloc.  Allocate and fill relocation array                          *
 ******************************************************************************/
static void    reloc_alloc(h, reloc)
mat_t   h;
real    reloc[][27];
{
   int  tx, ty, tz;
   for(tx = -1; tx <= 1; tx++)
      for(ty = -1; ty <= 1; ty++)
         for(tz = -1; tz <= 1; tz++)
         {
            reloc[0][NREL(tx,ty,tz)] = tx*h[0][0] + ty*h[0][1] + tz*h[0][2];
            reloc[1][NREL(tx,ty,tz)] = ty*h[1][1] + tz*h[1][2];
            reloc[2][NREL(tx,ty,tz)] = tz*h[2][2];
         }
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
static void calcdist(j0, nnab, r_sqr, rx, ry, rz, nab_sx, nab_sy, nab_sz, sitei)
int j0, nnab;
real r_sqr[], rx[], ry[], rz[], nab_sx[], nab_sy[], nab_sz[];
vec_t  sitei;
{
   int jsite;
VECTORIZE
            for(jsite=j0; jsite < nnab; jsite++)
            {
               rx[jsite] = nab_sx[jsite] - sitei[0];
               ry[jsite] = nab_sy[jsite] - sitei[1];
               rz[jsite] = nab_sz[jsite] - sitei[2];
               r_sqr[jsite] = rx[jsite]*rx[jsite] + ry[jsite]*ry[jsite]
                                                  + rz[jsite]*rz[jsite];
            }
}
/******************************************************************************
 * Force_calc.   This is the main intermolecular site force calculation       *
 * routine                                                                    *
 ******************************************************************************/
void force_calc(site, site_force, system, species, chg, potpar, pe, stress)
vec_t           site[],                 /* Site co-ordinate arrays       (in) */
                site_force[];           /* Site force arrays            (out) */
system_p        system;                 /* System struct                 (in) */
spec_t          species[];              /* Array of species records      (in) */
real            chg[];                  /* Array of site charges         (in) */
pot_t           potpar[];               /* Array of potential parameters (in) */
double          *pe;                    /* Potential energy             (out) */
mat_t           stress;                 /* Stress virial                (out) */
{
   int          ispec,                  /* Counter for species                */
                imol,                   /* Counter for molecules i            */
                icell,			/* Index for cells of molecule pair   */
                nnab,	j0,		/* Number of sites in neighbour list  */
                isite, jsite,           /* Site counter i,j                   */
   		i_id, jnab, lim, ipot;	/* Miscellaneous		      */
#ifdef FKERNEL
   int          nnab2;
#endif
   short        ix, iy, iz;		/* 3-d cell indices for ref and neig. */
   int          nsites = system->nsites,/* Local copy to keep optimiser happy */
                n_potpar = system->n_potpar,
   		max_id = system->max_id;
   double       norm = 2.0*control.alpha/sqrt(PI);	/* Coulombic prefactor*/
   int          ncells = nx*ny*nz;	/* Total number of cells	      */
   int          tx, ty, tz;             /* Temporaries for # unit cell shifts */
   int          *id_ptr;                /* Pointer to 'id' array              */
   		/*
		 * The following arrays are for 'neighbour site list'
		 * quantities and should be dimensioned to the max value of
		 * 'nnab'. Nsites is certainly too big, but as 'ewald' uses
		 * more store, it shouldn't increase total program use.
		 */
   int          *nab  = ialloc(nsites),		/* Neigbour site gather vector*/
                *nab3 = ialloc(nsites),		/* Stride 3 version of above  */
   		*id = ialloc(nsites),   	/* Array of site_id[nsites]   */
                *reloc_i = ialloc(nsites);	/* Vector of pbc relocations  */
   real         *R = dalloc(nsites);    	/* pbc site relocation cpt    */
   real         *nab_sx  = dalloc(nsites),	/* 'Gathered' list of         */
   		*nab_sy  = dalloc(nsites),	/*   neighbour site co-ords   */
                *nab_sz  = dalloc(nsites),	/*   - x,y,z components.      */
                *forcejx = dalloc(nsites),	/* List of neighbour site     */
                *forcejy = dalloc(nsites),	/*  forces in gathered form   */
                *forcejz = dalloc(nsites),	/*  - xyz components.	      */
                *rx      = dalloc(nsites),	/* Reference to neigbour site */
                *ry      = dalloc(nsites),	/* - site vector adjusted for */
                *rz      = dalloc(nsites),	/*  periodic boundaries. xyz. */
                *r_sqr   = dalloc(nsites),	/* Squared site-site distance */
                *nab_chg = dalloc(nsites),	/* Gathered neig. site charges*/
                *forceij = dalloc(nsites);	/* -V'(r) / r		      */
#ifdef VCALLS
   real         *force_comp = dalloc(nsites),
#else
   real         force_cpt, site0, site1, site2, s00, s01, s02, s11, s12, s22;
#endif
   real         reloc_v[3][27];			/* PBC relocation vectors     */
   real         ***potp				/* Expanded pot'l parameters  */
   		= (real***)arralloc(sizeof(real), 3,
                                    0,n_potpar-1, 1, max_id-1, 0, nsites-1);
   real         **nab_pot			/* Gathere'd pot par array    */
   		= (real**)arralloc(sizeof(real), 2,
				   0, n_potpar-1, 0, nsites-1);
   cell_t       *mol = aalloc(system->nmols, cell_t );
   spec_p       spec;
   cell_t       *cmol;
   cell_t       **cell;
   
   static boolean       init = true;
   static int           n_nabors;
   static ivec_t        *nabor;
   static reloc_t       ***reloc;
  

   if(init)
   {
      if(control.subcell <= 0.0) control.subcell = control.cutoff/5.0;
      nx = system->h[0][0]/control.subcell+0.5;
      ny = system->h[1][1]/control.subcell+0.5;
      nz = system->h[2][2]/control.subcell+0.5;
      ncells = nx*ny*nz;
      note("MD cell divided into %d subcells (%dx%dx%d)",ncells,nx,ny,nz);
      nabor = aalloc(CELLMAX*ncells, ivec_t );
      neighbour_list(nabor, &n_nabors, system->h, control.cutoff);
   
      reloc = (reloc_t***)arralloc(sizeof(reloc_t),3,
                                   -nx, 2*nx-1, -ny, 2*ny-1, -nz, 2*nz-1);
 
      for(ix = 0; ix < nx; ix++)
         for(tx = -nx; tx <= nx; tx += nx)
            for(iy = 0; iy < ny; iy++)
               for(ty = -ny; ty <= ny; ty += ny)
                  for(iz = 0; iz < nz; iz++)
                     for(tz = -nz; tz <= nz; tz += nz)
                     {
                        reloc[ix+tx][iy+ty][iz+tz].ncell = NCELL(ix, iy, iz);
                        reloc[ix+tx][iy+ty][iz+tz].rel= NREL(tx/nx,ty/ny,tz/nz);
                     }
 
      init = false;
   }

   cell = aalloc(ncells, cell_t *);

/*  Construct and fill expanded site-identifier array, id                     */
   id_ptr = id;
   for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
      for(imol = 0; imol < spec->nmols; imol++)
      {
         (void)memcpy((char*)id_ptr, (char*)spec->site_id, 
                      (int)spec->nsites*sizeof(int));
         id_ptr += spec->nsites;
      }
/*   Build arrays of pot. pars [max_id][nsites] for access in vector loops    */
   for(ipot = 0; ipot < n_potpar; ipot++)
      for(i_id = 1; i_id < max_id; i_id++)
         for(isite = 0; isite < nsites; isite++)
            potp[ipot][i_id][isite] = potpar[i_id*max_id+id[isite]].p[ipot];

   fill_cells(system->c_of_m, system->nmols, species, cell, mol);
   reloc_alloc(system->h, reloc_v);

   s00 = s01 = s02 = s11 = s12 = s22 = 0.0;	/* Accumulators for stress    */
/******************************************************************************
 *  Start of main loops.  Loop over species, ispec, molecules, imol, and      *
 *  sites, isite_mol on imol.  Isite is index into 'site' of site specified   *
 *  by isite_mol, imol and ispec.  Jbase is set to first site of imol+1, ispec*
 *  so jsite, counts from jbase to end.  Thus isite, jsite run over all site  *
 *  site pairs on distinct molecules.                                         *
 ******************************************************************************/
   icell = -1;
   for(ix = 0; ix < nx; ix++)
      for(iy = 0; iy < ny; iy++)
         for(iz = 0; iz < nz; iz++)
   {
      if(cell[++icell] == NULL) continue;       /* Empty cell - go on to next */
      nnab = 0;
#ifdef DEBUG
      printf("Working on cell %4d (%d,%d,%d) (sites %4d to %4d)\n", icell,
             ix,iy,iz,cell[icell].isite,cell[icell].isite+cell[icell].num-1);
      printf("\n jcell\tjx jy jz\tNsites\n");
#endif
      /*
       * Build site neighbour list 'nab' from cell list.
       */ 
      nnab = site_neighbour_list(nab, reloc_i,
				 n_nabors, ix, iy, iz, nabor, cell, reloc);
      if(nnab > nsites) message(NULLI,NULLP,FATAL,TONAB,nnab);

VECTORIZE
      for(jnab = 0; jnab < nnab; jnab++)     /* Build stride 3 nabor list  */
         nab3[jnab] = 3*nab[jnab];

      gather(nnab, nab_sx, site[0]  , nab3); /* Construct list of site     */
      gather(nnab, nab_sy, site[0]+1, nab3); /* co-ordinated from nabor    */
      gather(nnab, nab_sz, site[0]+2, nab3); /* list.                      */

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

      j0 = 0;
            			/* Loop over all molecules in cell icell. */
      for(cmol = cell[icell]; cmol != NULL; cmol = cmol->next)
      {
	 j0 += cmol->num;
         lim = cmol->isite + cmol->num;
         for(isite = cmol->isite; isite < lim; isite++)
         {                                   /* Loop over sites in cell icell */
            for(ipot = 0; ipot < n_potpar; ipot++)
               gather(nnab, nab_pot[ipot], potp[ipot][id[isite]], nab);
#ifdef DEBUG
            if(isite == 100)
            for(jnab = j0; jnab < nnab; jnab++)
               printf("%4d %4d\n", jnab,nab[jnab]);
#endif

            calcdist(j0, nnab, r_sqr, rx, ry, rz,
		     nab_sx,nab_sy,nab_sz,site[isite]);

            if( (jsite = j0+search_lt(nnab-j0, r_sqr+j0, 1, TOO_CLOSE)) < nnab )
               message(NULLI, NULLP, WARNING, TOOCLS,
		       isite, nab[jsite], sqrt(TOO_CLOSE));

            /*  Call the potential function kernel                            */
#ifdef FKERNEL
            nnab2 = nnab-j0;
            KERNEL(&nnab2, forceij+j0, pe, r_sqr+j0, nab_chg+j0, chg+isite,
		   &norm, &control.alpha, &system->ptype, &nsites, nab_pot[0]);
#else
            kernel(nnab-j0, forceij+j0, pe, r_sqr+j0, nab_chg+j0, chg[isite],
		   norm, control.alpha, system->ptype, nab_pot);
#endif
#ifdef VCALLS
            vaadd(nnab-j0, forcejx+j0, force_comp+j0, forceij+j0, rx+j0);
            site_force[isite][0] -= sum(nnab-j0, force_comp+j0,1);
            s00                  += vdot(nnab-j0, force_comp+j0, 1, rx+j0, 1);
            s01                  += vdot(nnab-j0, force_comp+j0, 1, ry+j0, 1);
            s02                  += vdot(nnab-j0, force_comp+j0, 1, rz+j0, 1);

            vaadd(nnab-j0, forcejy+j0, force_comp+j0, forceij+j0, ry+j0);
            site_force[isite][1] -= sum(nnab-j0, force_comp+j0,1);
            s11                  += vdot(nnab-j0, force_comp+j0, 1, ry+j0, 1);
            s12                  += vdot(nnab-j0, force_comp+j0, 1, rz+j0, 1);

            vaadd(nnab-j0, forcejz+j0, force_comp+j0, forceij+j0, rz+j0);
            site_force[isite][2] -= sum(nnab-j0, force_comp+j0,1);
            s22                  += vdot(nnab-j0, force_comp+j0, 1, rz+j0, 1);
#else
            site0 = site1 = site2 = 0.0;
VECTORIZE
            for(jsite=j0; jsite < nnab; jsite++)
            {
               force_cpt  = forceij[jsite]*rx[jsite];
               site0           -= force_cpt;
               forcejx[jsite]  += force_cpt;
               s00         += force_cpt * rx[jsite];
               s02         += force_cpt * rz[jsite];
               s01         += force_cpt * ry[jsite];
            }
VECTORIZE
            for(jsite=j0; jsite < nnab; jsite++)
            {
               force_cpt = forceij[jsite]*ry[jsite];
               site1           -= force_cpt;
               forcejy[jsite]  += force_cpt;
               s11         += force_cpt * ry[jsite];
               s12         += force_cpt * rz[jsite];
            }
VECTORIZE
            for(jsite=j0; jsite < nnab; jsite++)
            {
               force_cpt = forceij[jsite]*rz[jsite];
               site2           -= force_cpt;
               forcejz[jsite]  += force_cpt;
               s22         += force_cpt * rz[jsite];
            }
            site_force[isite][0] += site0;
            site_force[isite][1] += site1;
            site_force[isite][2] += site2;
#endif
         }
      }
      gather(nnab, forceij, site_force[0], nab3);
      vadd(nnab, forceij, forcejx);
      scatter(nnab, site_force[0], nab3, forceij);
      gather(nnab, forceij, site_force[0]+1, nab3);
      vadd(nnab, forceij, forcejy);
      scatter(nnab, site_force[0]+1, nab3, forceij);
      gather(nnab, forceij, site_force[0]+2, nab3);
      vadd(nnab, forceij, forcejz);
      scatter(nnab, site_force[0]+2, nab3, forceij);
   }
   stress[0][0]  += s00;
   stress[0][1]  += s01;
   stress[0][2]  += s02;
   stress[1][1]  += s11;
   stress[1][2]  += s12;
   stress[2][2]  += s22;

   cfree((char*)potp[0][1]);     cfree((char*)nab_pot[0]);
   cfree((char*)r_sqr);   cfree((char*)mol); 
   cfree((char*)R);       cfree((char*)reloc_i);  cfree((char*)nab_chg);
   cfree((char*)rx);      cfree((char*)ry);      cfree((char*)rz);
   cfree((char*)forcejx); cfree((char*)forcejy); cfree((char*)forcejz);
   cfree((char*)nab);     cfree((char*)nab3);    cfree((char*)cell);
   cfree((char*)nab_sx);  cfree((char*)nab_sy);  cfree((char*)nab_sz);
   cfree((char*)id);      cfree((char*)forceij);
#ifdef VCALLS
   cfree((char*)force_comp);
#endif
}
