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
 * Parallel - support and interface routines to parallel MP libraries.	      *
 ******************************************************************************
 *       $Log: parallel.c,v $
 * Revision 1.1  1994/07/07  16:59:44  keith
 * Initial revision
 *
 */
/*========================== program include files ===========================*/
#include	"defs.h"
#include	"structs.h"
#include	"messages.h"
/*========================== system  include files ===========================*/
#include	<signal.h>
#include	"string.h"
#ifdef TCGMSG
#include	<sndrcv.h>
#endif
#ifdef MPI
#include	<mpi.h>
#endif
static long	lval;
#define		ADDR(expr) (lval=(expr),&lval)
#define M_REAL (sizeof(real)==sizeof(double)?MPI_DOUBLE:MPI_FLOAT)
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
gptr		*av_ptr();
void		init_averages();
void		allocate_dynamics();
extern int 	ithread, nthreads;
/*====================== Utilities for interface functions ===================*/
#ifdef TCGMSG
/*
 * This gets called a lot.  A "prettier" C interface to BRDCST_.
 */
cbrdcst(type, buf, lenbuf, ifrom)
long 	type;
gptr	*buf;
long	lenbuf,
	ifrom;
{
   switch (type & ~(long)0x7FFF)
   {
   case MSGINT: lenbuf *= sizeof(int);
      break;
   case MSGCHR: lenbuf *= sizeof(char);
      break;
   case MSGDBL: lenbuf *= sizeof(double);
      break;
   default:
      break;
   }
   BRDCST_(&type, buf, &lenbuf, &ifrom);
}
#endif
#ifdef BSP
/*
 * BSP currently doesn't handle auto & heap vars.  This interface
 * copies to static storage.
 */
#define NBUFMAX 1048576
static char tmpbuf[NBUFMAX];

cbrdcst(ifrom, inbuf, outbuf, size)
int	ifrom;
char	*inbuf, *outbuf;
int	size;
{
   if( size > NBUFMAX)
      message(NULLI, NULLP, FATAL, 
	      "cbrdcst called with too large n (%d) - maximum is %d",
	      size, NBUFMAX);
   memcp(tmpbuf, inbuf, size);
   bspbroadcast(ifrom, tmpbuf, tmpbuf, size);
   memcp(outbuf, tmpbuf, size);
}
#endif
/*====================== Parallel lib interface functions ====================*/
/******************************************************************************
 * par_sigintreset().  Reset signal handler to parallel lib default upon trap *
 *		       of SIGINT.					      *
 ******************************************************************************/
#ifdef TCGMSG
extern	void    SigintHandler();

par_sigintreset()
{
   signal(SIGINT, SigintHandler);
}
#endif
#ifdef BSP
par_sigintreset()
{
   signal(SIGINT, SIG_DFL);
}
#endif
#ifdef MPI
par_sigintreset()
{
   signal(SIGINT, SIG_DFL);
}
#endif
/******************************************************************************
 * par_imax().  Calculate global maximum over all processors.		      *
 ******************************************************************************/
#ifdef TCGMSG
par_imax(idat)
int *idat;
{
       IGOP_(ADDR(10+MSGINT), idat, ADDR(1), "max");
}
#endif
#ifdef BSP
void imax(i1, i2, i3, size)
int *i1, *i2, *i3;
int	size;
{
  *i1 = MAX(*i2,*i3);
}

par_imax(idat)
int *idat;
{
   bspreduce(imax, idat, idat, sizeof(int));
}
#endif
#ifdef MPI
par_imax(idat)
int *idat;
{
   int result;
   MPI_Allreduce(idat, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   *idat = result;
}
#endif
/******************************************************************************
 * par_rsum()/dsum.  Calculate sum of "reals"/doubles  over all processors.   *
 ******************************************************************************/
#ifdef TCGMSG
par_rsum(buf, n)
real *buf;
int  n;
{
   DGOP_(ADDR(MSGDBL), buf, &n, "+");
}
par_dsum(buf, n)
real *buf;
int  n;
{
   DGOP_(ADDR(MSGDBL), buf, &n, "+");
}
#endif
#ifdef BSP
void vradd(res, x, y, nb)
real res[], x[], y[];
int nb;
{
   int i, n=nb/sizeof(real);
   for(i = 0; i < n; i++)
      res[i] = x[i] + y[i];
}
void vdadd(res, x, y, nb)
double res[], x[], y[];
int nb;
{
   int i, n=nb/sizeof(double);
   for(i = 0; i < n; i++)
      res[i] = x[i] + y[i];
}

par_rsum(buf, n)
real *buf;
int  n;
{
   if( n*sizeof(real) > NBUFMAX)
      message(NULLI, NULLP, FATAL, 
	      "par_rsum called with too large n (%d) - maximum is %d",
	      n, NBUFMAX/sizeof(real));
   memcp(tmpbuf, buf, n*sizeof(real));
   bspreduce(vradd, tmpbuf, tmpbuf, n*sizeof(real));
   memcp(buf, tmpbuf, n*sizeof(real));
}
par_dsum(buf, n)
double *buf;
int  n;
{
   if( n*sizeof(double) > NBUFMAX)
      message(NULLI, NULLP, FATAL, 
	      "par_rsum called with too large n (%d) - maximum is %d",
	      n, NBUFMAX/sizeof(double));
   memcp(tmpbuf, buf, n*sizeof(double));
   bspreduce(vdadd, tmpbuf, tmpbuf, n*sizeof(double));
   memcp(buf, tmpbuf, n*sizeof(double));
}
#endif
#ifdef MPI
/*
 * MPI demands seperate send and receive buffers.  Malloc one and keep
 * it around.  Extend if necessary.
 */
par_rsum(buf, n)
real *buf;
int  n;
{
   static real *tmpbuf = 0;
   static int  tmpsize = 0;
   if(n > tmpsize)
   {
      if( tmpbuf )
	 free(tmpbuf);
      tmpbuf = dalloc(n);
      tmpsize = n;
   }
   MPI_Allreduce(buf, tmpbuf, n, M_REAL, MPI_SUM, MPI_COMM_WORLD);
   memcp(buf, tmpbuf, n*sizeof(real));
}
par_dsum(buf, n)
double *buf;
int  n;
{
   static double *tmpbuf = 0;
   static int  tmpsize = 0;
   if(n > tmpsize)
   {
      if( tmpbuf )
	 free(tmpbuf);
      tmpbuf = aalloc(n, double);
      tmpsize = n;
   }
   MPI_Allreduce(buf, tmpbuf, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   memcp(buf, tmpbuf, n*sizeof(double));
}
#endif
/******************************************************************************
 * par_begin().  Initialize parallel libs.				      *
 ******************************************************************************/
#ifdef TCGMSG
par_begin(argc, argv, ithread, nthreads)
int	*argc;
char	***argv;
int	*ithread;
int	*nthreads;
{
   int i;
   PBEGIN_(*argc, *argv);
   *nthreads = NNODES_();
   *ithread  = NODEID_();
   for(i = 1; i < *argc && strcmp((*argv)[i],"-master"); i++)
      ;
   *argc = i;
}
#endif
#ifdef BSP
par_begin(argc, argv, ithread, nthreads)
int	*argc;
char	***argv;
int	*ithread;
int	*nthreads;
{
   bspstart(0, nthreads, ithread);
}
#endif
#ifdef MPI
par_begin(argc, argv, ithread, nthreads)
int	*argc;
char	***argv;
int	*ithread;
int	*nthreads;
{
   MPI_Init(argc, argv);
   MPI_Comm_size(MPI_COMM_WORLD, nthreads);
   MPI_Comm_rank(MPI_COMM_WORLD, ithread);
}
#endif
/******************************************************************************
 * par_finish().  Parallel lib wind-up function.			      *
 ******************************************************************************/
#ifdef TCGMSG
par_finish()
{
   PEND_();
}
#endif
#ifdef BSP
par_finish()
{
   bspfinish();
}
#endif
#ifdef MPI
par_finish()
{
   MPI_Finalize();
}
#endif
/******************************************************************************
 * par_abort().  Parallel lib abort function.				      *
 ******************************************************************************/
#ifdef TCGMSG
par_abort(code)
int code;
{
   Error("",code);
}
#endif
#ifdef BSP
par_abort(code)
int code;
{
#ifdef NOTYET
   bspabort(code);
#endif
}
#endif
#ifdef MPI
par_abort(code)
int code;
{
   MPI_Abort(MPI_COMM_WORLD, code);
}
#endif
/******************************************************************************
 *  copy_sysdef                                                            *
 ******************************************************************************/
void	copy_sysdef(system, spec_ptr, site_info, pot_ptr)
system_mp	system;			/* Pointer to system array (in main)  */
spec_mp		*spec_ptr;		/* Pointer to be set to species array */
site_mp		*site_info;		/* To be pointed at site_info array   */
pot_mp		*pot_ptr;		/* To be pointed at potpar array      */
{
   spec_mp	spec;
   int		n_pot_recs;

   /*
    * Fetch "system" struct
    */
#ifdef BSP
   cbrdcst(0, (gptr*)system, (gptr*)system, lsizeof(system_mt));
#endif
   /*
    * N.B.  TCGMSG does not support transfer of structs using xdr.
    * Therefore this code will only work between homogeneous nodes.
    * (But that's all BSP does anyway)
    */
#ifdef TCGMSG
   cbrdcst(0, (gptr*)system, lsizeof(system_mt), 0);
#endif
   /*
    * MPI has facilities to describe struct types and transfer.
    * Until someone gets around to actually doing it, this will
    * work on homogeneous processors.
    */
#ifdef MPI
   MPI_Bcast((gptr*)system, lsizeof(system_mt), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
   /* Allocate space for species, site_info and potpar arrays and set pointers*/
   if( ithread > 0 )
   {
      *spec_ptr  = aalloc(system->nspecies,                spec_mt );
      *site_info = aalloc(system->max_id,                  site_mt );
      *pot_ptr   = aalloc(system->max_id * system->max_id, pot_mt );
   }
   /*  read species array into allocated space				      */
#ifdef BSP
   cbrdcst(0,(gptr*)*spec_ptr,(gptr*)*spec_ptr, 
	      system->nspecies*lsizeof(spec_mt));
#endif
#ifdef TCGMSG
   cbrdcst(0, (gptr*)*spec_ptr, system->nspecies*lsizeof(spec_mt),0);
#endif
#ifdef MPI
   MPI_Bcast((gptr*)*spec_ptr, system->nspecies*lsizeof(spec_mt), 
	     MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
   /* Fill p_f_sites and site_id arrays for each species */
   for (spec = *spec_ptr; spec < &(*spec_ptr)[system->nspecies]; spec++)
   {
      if( ithread > 0 )
      {      
	 spec->p_f_sites = ralloc(spec->nsites);  /* Allocate the species -     */
	 spec->site_id   = ialloc(spec->nsites);  /* specific arrays	      */
      }
#ifdef BSP
      cbrdcst(0,(gptr*)spec->p_f_sites,(gptr*)spec->p_f_sites, 
		 3*spec->nsites*lsizeof(real));
      cbrdcst(0,(gptr*)spec->site_id,(gptr*)spec->site_id,   
		 spec->nsites*lsizeof(int));
#endif
#ifdef TCGMSG
      cbrdcst(MSGDBL, (gptr*)spec->p_f_sites, 3*spec->nsites, 0);
      cbrdcst(MSGINT, (gptr*)spec->site_id,     spec->nsites, 0);
#endif
#ifdef MPI
      MPI_Bcast((gptr*)spec->p_f_sites, 3*spec->nsites, M_REAL, 
		0, MPI_COMM_WORLD);
      MPI_Bcast((gptr*)spec->site_id, spec->nsites, MPI_INT, 
		0, MPI_COMM_WORLD);
#endif
   }
   /* Fill site_info array */
#ifdef BSP
   cbrdcst(0, (gptr*)*site_info, (gptr*)*site_info, 
	      system->max_id*lsizeof(site_mt));
#endif
#ifdef TCGMSG
   cbrdcst(0, (gptr*)*site_info, system->max_id*lsizeof(site_mt), 0);
#endif
#ifdef MPI
   MPI_Bcast((gptr*)*site_info, system->max_id*lsizeof(site_mt), 
	     MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
   /*
    * Potential Parameters.
    */
   n_pot_recs = SQR(system->max_id);
   if( ithread > 0 )
      *pot_ptr = (pot_mt*)aalloc(n_pot_recs*sizeof(pot_mt), char);
#ifdef BSP
   cbrdcst(0, (gptr*)*pot_ptr, (gptr*)*pot_ptr, n_pot_recs*sizeof(pot_mt));
#endif
#ifdef TCGMSG
   cbrdcst(0, (gptr*)*pot_ptr, n_pot_recs*sizeof(pot_mt), 0);
#endif
#ifdef MPI
   MPI_Bcast((gptr*)*pot_ptr, n_pot_recs*sizeof(pot_mt), 
	     MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
}
/******************************************************************************
 *  copy_dynamics()							      *
 ******************************************************************************/
void	copy_dynamics(system)
system_mp	system;
{
   gptr		*ap;			/* Pointer to averages database       */
   size_mt	asize;			/* Size of averages database	      */

#ifdef BSP
   cbrdcst(0, (gptr*)system->c_of_m, (gptr*)system->c_of_m, 
	      3*system->nmols*lsizeof(real));
   cbrdcst(0, (gptr*)system->vel, (gptr*)system->vel,    
	      3*system->nmols*lsizeof(real));
   cbrdcst(0, (gptr*)system->velp, (gptr*)system->velp,   
	      3*system->nmols*lsizeof(real));
   cbrdcst(0, (gptr*)system->acc, (gptr*)system->acc,    
	      3*system->nmols*lsizeof(real));
   cbrdcst(0, (gptr*)system->acco, (gptr*)system->acco,   
	      3*system->nmols*lsizeof(real));
   cbrdcst(0, (gptr*)system->accvo, (gptr*)system->accvo,  
	      3*system->nmols*lsizeof(real));
   if(system->nmols_r > 0)
   {
      cbrdcst(0, (gptr*)system->quat, (gptr*)system->quat,    
		 4*system->nmols_r*lsizeof(real));
      cbrdcst(0, (gptr*)system->qdot, (gptr*)system->qdot,    
		 4*system->nmols_r*lsizeof(real));
      cbrdcst(0, (gptr*)system->qdotp, (gptr*)system->qdotp,   
		 4*system->nmols_r*lsizeof(real));
      cbrdcst(0, (gptr*)system->qddot, (gptr*)system->qddot,   
		 4*system->nmols_r*lsizeof(real));
      cbrdcst(0, (gptr*)system->qddoto, (gptr*)system->qddoto,  
		 4*system->nmols_r*lsizeof(real));
      cbrdcst(0, (gptr*)system->qddotvo, (gptr*)system->qddotvo, 
		 4*system->nmols_r*lsizeof(real));
   }
   cbrdcst(0, (gptr*)system->h,      (gptr*)system->h,       9*lsizeof(real));
   cbrdcst(0, (gptr*)system->hdot,   (gptr*)system->hdot,    9*lsizeof(real));
   cbrdcst(0, (gptr*)system->hdotp,  (gptr*)system->hdotp,   9*lsizeof(real));
   cbrdcst(0, (gptr*)system->hddot,  (gptr*)system->hddot,   9*lsizeof(real));
   cbrdcst(0, (gptr*)system->hddoto, (gptr*)system->hddoto,  9*lsizeof(real));
   cbrdcst(0, (gptr*)system->hddotvo,(gptr*)system->hddotvo, 9*lsizeof(real));

   ap = av_ptr(&asize,0);	      /* get addr, size of database   */
   cbrdcst(0,ap,ap, asize);
#endif
#ifdef TCGMSG
   cbrdcst(MSGDBL, (gptr*)system->c_of_m,3*system->nmols, 0);
   cbrdcst(MSGDBL, (gptr*)system->vel,   3*system->nmols, 0);
   cbrdcst(MSGDBL, (gptr*)system->velp,  3*system->nmols, 0);
   cbrdcst(MSGDBL, (gptr*)system->acc,   3*system->nmols, 0);
   cbrdcst(MSGDBL, (gptr*)system->acco,  3*system->nmols, 0);
   cbrdcst(MSGDBL, (gptr*)system->accvo, 3*system->nmols, 0);
   if(system->nmols_r > 0)
   {
      cbrdcst(MSGDBL, (gptr*)system->quat,    4*system->nmols_r, 0);
      cbrdcst(MSGDBL, (gptr*)system->qdot,    4*system->nmols_r, 0);
      cbrdcst(MSGDBL, (gptr*)system->qdotp,   4*system->nmols_r, 0);
      cbrdcst(MSGDBL, (gptr*)system->qddot,   4*system->nmols_r, 0);
      cbrdcst(MSGDBL, (gptr*)system->qddoto,  4*system->nmols_r, 0);
      cbrdcst(MSGDBL, (gptr*)system->qddotvo, 4*system->nmols_r, 0);
   }
   cbrdcst(MSGDBL, (gptr*)system->h,       9, 0); 
   cbrdcst(MSGDBL, (gptr*)system->hdot,    9, 0);
   cbrdcst(MSGDBL, (gptr*)system->hdotp,   9, 0); 
   cbrdcst(MSGDBL, (gptr*)system->hddot,   9, 0);
   cbrdcst(MSGDBL, (gptr*)system->hddoto,  9, 0); 
   cbrdcst(MSGDBL, (gptr*)system->hddotvo, 9, 0);

   ap = av_ptr(&asize,0);	      /* get addr, size of database   */
   cbrdcst(0, ap, asize,0);
#endif
#ifdef MPI
   MPI_Bcast((gptr*)system->c_of_m,3*system->nmols, M_REAL, 0, MPI_COMM_WORLD);
   MPI_Bcast((gptr*)system->vel,   3*system->nmols, M_REAL, 0, MPI_COMM_WORLD);
   MPI_Bcast((gptr*)system->velp,  3*system->nmols, M_REAL, 0, MPI_COMM_WORLD);
   MPI_Bcast((gptr*)system->acc,   3*system->nmols, M_REAL, 0, MPI_COMM_WORLD);
   MPI_Bcast((gptr*)system->acco,  3*system->nmols, M_REAL, 0, MPI_COMM_WORLD);
   MPI_Bcast((gptr*)system->accvo, 3*system->nmols, M_REAL, 0, MPI_COMM_WORLD);
   if(system->nmols_r > 0)
   {
      MPI_Bcast((gptr*)system->quat,    4*system->nmols_r, 
		M_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast((gptr*)system->qdot,    4*system->nmols_r, 
		M_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast((gptr*)system->qdotp,   4*system->nmols_r, 
		M_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast((gptr*)system->qddot,   4*system->nmols_r, 
		M_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast((gptr*)system->qddoto,  4*system->nmols_r, 
		M_REAL, 0, MPI_COMM_WORLD);
      MPI_Bcast((gptr*)system->qddotvo, 4*system->nmols_r, 
		M_REAL, 0, MPI_COMM_WORLD);
   }
   MPI_Bcast((gptr*)system->h,       9, M_REAL, 0, MPI_COMM_WORLD); 
   MPI_Bcast((gptr*)system->hdot,    9, M_REAL, 0, MPI_COMM_WORLD);
   MPI_Bcast((gptr*)system->hdotp,   9, M_REAL, 0, MPI_COMM_WORLD); 
   MPI_Bcast((gptr*)system->hddot,   9, M_REAL, 0, MPI_COMM_WORLD);
   MPI_Bcast((gptr*)system->hddoto,  9, M_REAL, 0, MPI_COMM_WORLD); 
   MPI_Bcast((gptr*)system->hddotvo, 9, M_REAL, 0, MPI_COMM_WORLD);

   ap = av_ptr(&asize,0);	      /* get addr, size of database   */
   MPI_Bcast(ap, asize, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
}
/******************************************************************************
 * replicate().  Make a copy of all Moldy's constant and dynamic data.        *
 *               This is for parallel implementations and allows "start_up"   *
 *		 to be called on one processor.  This version calls the       *
 *		 Oxford BSP library.					      *
 ******************************************************************************/
void replicate(control, system, spec_ptr, site_info, pot_ptr, restart_header)
contr_mt  *control;
system_mt *system;
spec_mt   **spec_ptr;
site_mt	  **site_info;
pot_mt    **pot_ptr;
restrt_mt *restart_header;
{
   int av_convert;
   /*
    *  Fetch the top-level structs
    */
#ifdef BSP
   cbrdcst(0, control, control, sizeof *control);
   cbrdcst(0, restart_header, restart_header, sizeof *restart_header);
#endif
#ifdef TCGMSG
   cbrdcst(0, (gptr*)control, sizeof *control, 0);
   cbrdcst(0, (gptr*)restart_header, sizeof *restart_header, 0);
#endif
#ifdef MPI
   MPI_Bcast((gptr*)control, sizeof *control, MPI_BYTE, 0, MPI_COMM_WORLD);
   MPI_Bcast((gptr*)restart_header, sizeof *restart_header, 
	     MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
   /*
    * Now get "species" struct array and internal dynamic arrays
    */
   copy_sysdef(system, spec_ptr, site_info, pot_ptr);
   if( ithread > 0 )
   {
      /*
       * Now we have all the information to create the dynamic variables
       */
      allocate_dynamics(system, *spec_ptr);
      /*
       * Initialise averages database.
       */
      init_averages(system->nspecies, (char*)0,
		    control->roll_interval, control->roll_interval,&av_convert);
   }
   /*
    * Copy the dynamic vars from the other processes.
    */
   copy_dynamics(system);
}

