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
 * Revision 2.10  1994/10/17  10:49:41  keith
 * Changed arg list of bspstart to match changed library version.
 *
 * Revision 2.9  1994/07/11  11:15:30  keith
 * Tidied up startup routine with par_broadcast() function.
 * Documented parallel routine interface calls for porting.
 *
 * Revision 2.8  1994/07/07  17:00:26  keith
 * Interface to BSP, TCGMSG and MPI message-passing libraries.
 *
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
#ifdef BSP
/*
 * BSP currently doesn't handle auto & heap vars.  This interface
 * copies to static storage.
 */
#define NBUFMAX 1048576
static char tmpbuf[NBUFMAX];
#endif
/*====================== Parallel lib interface functions ====================*
 *  The following set of functions define the interface between moldy and     *
 *  a message-passing parallel library.  To port to a new library it suffices *
 *  to implement these in terms of the library primitives.  Thus all of the   *
 *  implementation-dependent parallel code is confined to this file.          *
 *									      *
 *  The functions required are:						      *
 *  par_sigintreset()	:  Moldy sets a handler for SIGINT.  This fn is called*
 *			   from the signal handler to restore the default.    *
 *  par_begin(int *argc, char ***argv, int *ithread, int *nthreads) 	      *
 *                      :  Initialize the library and return the number of    *
 *			   processes and the ID of this process.	      *
 *  par_finish()	:  Terminate the parallel run normally.		      *
 *  par_abort(int code) :  Terminate the run abnormally.  Return code if poss.*
 *  par_broadcast(void *buf, int n, size_mt size, int ifrom)		      *
 *			:  Broadcast the specified buffer from ifrom to all.  *
 *  par_{r,d}sum(void *buf, int n) :  Perform a global parallel sum reduction *
 *			   on the buffer containing n {reals,doubles}.	      *
 *  par_imax(int *idat) :  Perform a global "maximum" reduction on the single *
 *       		   int argument.				      *
 *									      *
 *  Note that there is no provision for heterogeneous execution by way of     *
 *  type identification.  Though some MP libraries (eg MPI) do provide the    *
 *  hooks it is too hard to implement for the control and other structs.      *
 *  It's also hard to see why this might ever be useful for a MD run.	      *
 *============================================================================*/
/******************************************************************************
 * par_sigintreset().  Reset signal handler to parallel lib default upon trap *
 *		       of SIGINT.					      *
 ******************************************************************************/
#ifdef TCGMSG
extern	void    SigintHandler();

void
par_sigintreset()
{
   signal(SIGINT, SigintHandler);
}
#endif
#ifdef BSP
void
par_sigintreset()
{
   signal(SIGINT, SIG_DFL);
}
#endif
#ifdef MPI
void
par_sigintreset()
{
   signal(SIGINT, SIG_DFL);
}
#endif
/******************************************************************************
 * par_imax().  Calculate global maximum over all processors.		      *
 ******************************************************************************/
#ifdef TCGMSG
void
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

void
par_imax(idat)
int *idat;
{
   bspreduce(imax, idat, idat, sizeof(int));
}
#endif
#ifdef MPI
void
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
void
par_rsum(buf, n)
real *buf;
int  n;
{
   DGOP_(ADDR(MSGDBL), buf, &n, "+");
}
void
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

void
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
void
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
void
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
void
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
 * par_broadcast(). Broadcast data to all processors.			      *
 ******************************************************************************/
#ifdef TCGMSG
void
par_broadcast(buf, n, size, ifrom)
gptr	*buf;
int	n;
size_mt	size;
int	ifrom;
{
   long type = 0;
   long	lenbuf = n*size;
   long ifrm = ifrom;
   BRDCST_(&type, buf, &lenbuf, &ifrm);
}
#endif
#ifdef BSP
void
par_broadcast(buf, n, size, ifrom)
gptr	*buf;
int	n;
size_mt	size;
int	ifrom;
{
   if( n*size > NBUFMAX)
      message(NULLI, NULLP, FATAL, 
	      "cbrdcst called with too large n (%d) - maximum is %d",
	      n*size, NBUFMAX);
   memcp(tmpbuf, buf, n*size);
   bspbroadcast(ifrom, tmpbuf, tmpbuf, n*size);
   memcp(buf, tmpbuf, n*size);
}
#endif
#ifdef MPI
void
par_broadcast(buf, n, size, ifrom)
gptr	*buf;
int	n;
size_mt	size;
int	ifrom;
{
   MPI_Bcast(buf, n*size, MPI_BYTE, ifrom, MPI_COMM_WORLD);
}
#endif
/******************************************************************************
 * par_begin().  Initialize parallel libs.				      *
 ******************************************************************************/
#ifdef TCGMSG
void
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
   for(i = 1; i < *argc; i++ )
      if( !strcmp((*argv)[i],"-master"))
      {
	 *argc = i-1;
	 (*argv)++;
	 break;
      }
}
#endif
#ifdef BSP
void
par_begin(argc, argv, ithread, nthreads)
int	*argc;
char	***argv;
int	*ithread;
int	*nthreads;
{
   bspstart(*argc, *argv, 0, nthreads, ithread);
}
#endif
#ifdef MPI
void
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
void
par_finish()
{
   PEND_();
}
#endif
#ifdef BSP
void
par_finish()
{
   bspfinish();
}
#endif
#ifdef MPI
void
par_finish()
{
   MPI_Finalize();
}
#endif
/******************************************************************************
 * par_abort().  Parallel lib abort function.				      *
 ******************************************************************************/
#ifdef TCGMSG
void
par_abort(code)
int code;
{
   Error("",code);
}
#endif
#ifdef BSP
void
par_abort(code)
int code;
{
#ifdef NOTYET
   bspabort(code);
#endif
}
#endif
#ifdef MPI
void
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
#ifdef SPMD
   /*
    * Fetch "system" struct
    */
   par_broadcast((gptr*)system, 1, lsizeof(system_mt), 0);
   /* Allocate space for species, site_info and potpar arrays and set pointers*/
   if( ithread > 0 )
   {
      *spec_ptr  = aalloc(system->nspecies,                spec_mt );
      *site_info = aalloc(system->max_id,                  site_mt );
      *pot_ptr   = aalloc(system->max_id * system->max_id, pot_mt );
   }
   /*  read species array into allocated space				      */
   par_broadcast((gptr*)*spec_ptr, system->nspecies, lsizeof(spec_mt),0);
   /* Fill p_f_sites and site_id arrays for each species */
   for (spec = *spec_ptr; spec < &(*spec_ptr)[system->nspecies]; spec++)
   {
      if( ithread > 0 )
      {      
	 spec->p_f_sites = ralloc(spec->nsites);  /* Allocate the species -     */
	 spec->site_id   = ialloc(spec->nsites);  /* specific arrays	      */
      }
      par_broadcast((gptr*)spec->p_f_sites, 3*spec->nsites, sizeof(real), 0);
      par_broadcast((gptr*)spec->site_id,     spec->nsites, sizeof(int), 0);
   }
   /* Fill site_info array */
   par_broadcast((gptr*)*site_info, system->max_id, sizeof(site_mt), 0);
   /*
    * Potential Parameters.
    */
   n_pot_recs = SQR(system->max_id);
   if( ithread > 0 )
      *pot_ptr = (pot_mt*)aalloc(n_pot_recs*sizeof(pot_mt), char);
   par_broadcast((gptr*)*pot_ptr, n_pot_recs, sizeof(pot_mt), 0);
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
#ifdef SPMD
   par_broadcast((gptr*)system->c_of_m,3*system->nmols, sizeof(real), 0);
   par_broadcast((gptr*)system->vel,   3*system->nmols, sizeof(real), 0);
   par_broadcast((gptr*)system->velp,  3*system->nmols, sizeof(real), 0);
   par_broadcast((gptr*)system->acc,   3*system->nmols, sizeof(real), 0);
   par_broadcast((gptr*)system->acco,  3*system->nmols, sizeof(real), 0);
   par_broadcast((gptr*)system->accvo, 3*system->nmols, sizeof(real), 0);
   if(system->nmols_r > 0)
   {
      par_broadcast((gptr*)system->quat,    4*system->nmols_r, sizeof(real), 0);
      par_broadcast((gptr*)system->qdot,    4*system->nmols_r, sizeof(real), 0);
      par_broadcast((gptr*)system->qdotp,   4*system->nmols_r, sizeof(real), 0);
      par_broadcast((gptr*)system->qddot,   4*system->nmols_r, sizeof(real), 0);
      par_broadcast((gptr*)system->qddoto,  4*system->nmols_r, sizeof(real), 0);
      par_broadcast((gptr*)system->qddotvo, 4*system->nmols_r, sizeof(real), 0);
   }
   par_broadcast((gptr*)system->h,       9, sizeof(real), 0); 
   par_broadcast((gptr*)system->hdot,    9, sizeof(real), 0);
   par_broadcast((gptr*)system->hdotp,   9, sizeof(real), 0); 
   par_broadcast((gptr*)system->hddot,   9, sizeof(real), 0);
   par_broadcast((gptr*)system->hddoto,  9, sizeof(real), 0); 
   par_broadcast((gptr*)system->hddotvo, 9, sizeof(real), 0);

   ap = av_ptr(&asize,0);	      /* get addr, size of database   */
   par_broadcast(ap, 1, asize,0);
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
#ifdef SPMD
   par_broadcast((gptr*)control, 1, sizeof *control, 0);
   par_broadcast((gptr*)restart_header, 1, sizeof *restart_header, 0);
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
#endif
}

