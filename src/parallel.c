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
 *       Revision 2.21.2.1.2.1  2000/12/07 15:58:28  keith
 *       Mainly cosmetic minor modifications and added special comments to
 *       shut lint up.
 *
 *       Revision 2.21.2.1  2000/08/29 16:49:19  keith
 *       Fixed RNG to be synchronous on multiprocessor -- needed for
 *       scale-options=8.
 *
 *       Revision 2.21  1998/07/17 14:54:06  keith
 *       Ported SHMEM version to IRIX 6/ SGI Origin 2000
 *
 *       Revision 2.21  1998/07/17 12:06:53  keith
 *       Ported SHMEM parallel interface to Irix for Origin systems
 *
 *       Revision 2.20  1997/08/26 16:40:59  keith
 *       Updated with interface to new BSPlib standard and tested with Oxford
 *       BSPlib.
 *       Entries for old BSP remain - define preprocessor symbol BSP0 to use.
 *
 *       Revision 2.19  1996/10/15 13:50:45  keith
 *       Corrections to Cray SHMEM library interface.
 *
 *       Revision 2.18  1996/09/03 15:04:51  keith
 *       Added optional code for MPI global sums to guarantee identical
 *       results on all processors. Compile with -DUNSYMM to activate.
 *
 *       Revision 2.17  1996/03/19 12:33:50  keith
 *       Optimised T3D par_collect_all.  Now does not call barrier()
 *       within loop, but uses pt-to-pt synchronization.
 *       (b) Works "in-place" in memory rather than in static buffers.
 *       This may well not be robust because of sbreak() call which
 *       may break malloc().
 *
 *       Revision 2.14  1996/01/17 17:08:42  keith
 *       Added "par_isum()" for rdf calculation.
 *       Added security "exit()" call to par_abort().
 *       Corrected bug where par_dsum called vradd.
 *       Added "init_rdf" call for all threads.
 *
 *       Revision 2.13  1995/12/22 11:42:04  keith
 *       Modified buffer handling for BSP interface.  It used to complain and
 *       stop if buffer was too small. Now it divides data into chunks smaller
 *       than the buffer and transfers them one at a time.
 *
 *       Revision 2.12  1995/12/06 10:44:50  keith
 *       Nose-Hoover and Gaussian (Hoover constrained) thermostats added.
 *
 *       Revision 2.11  1994/11/24 14:48:17  keith
 *       Fixed problem with arg lists for TCGMSG
 *
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
#ifndef lint
static char *RCSid = "$Header: /home/minphys2/keith/CVS/moldy/src/parallel.c,v 2.21.2.1.2.1 2000/12/07 15:58:28 keith Exp $";
#endif
/*========================== program include files ===========================*/
#include	"defs.h"
#include	"structs.h"
#include	"messages.h"
/*========================== system  include files ===========================*/
#include	<signal.h>
#include	"string.h"
#ifdef BSP
#include	<bsp.h>
#include	<bsp_level1.h>
#endif
#ifdef DEBUG
#include        <stdio.h>
#endif
#ifdef TCGMSG
#include	<sndrcv.h>
#endif
#ifdef MPI
#include	<mpi.h>
#endif
#ifdef SHMEM
#include	<malloc.h>
#include	<mpp/shmem.h>
#include        <unistd.h>
#ifndef HAVE_SBREAK
#define sbreak(x) sbrk(x)
#endif
#ifndef HAVE_SHMEM_LONG_GET
#define shmem_long_get shmem_get
#define shmem_long_put shmem_put
#endif
#endif
/*========================== External function declarations ==================*/
gptr            *talloc();	       /* Interface to memory allocator       */
void		tfree();
gptr		*av_ptr();
gptr            *rdf_ptr();
void		init_averages();
void		allocate_dynamics();
void		init_rdf();
unsigned long	getseed();
void		smdrand();
extern int 	ithread, nthreads;
/*====================== Utilities for interface functions ===================*/
#ifdef TCGMSG
static long	lval;
#define		ADDR(expr) (lval=(expr),&lval)
#endif
#ifdef BSP0
/*
 * BSP version 0 doesn't handle auto & heap vars.  This interface
 * copies to static storage.  A 1MB buffer is large enough to
 * handle up to 44000 atomic sites in one go, and data for larger 
 * systems is parcelled up appropriately and sent in chunks.
 */
#define NBUFMXW 131072
#define NBUFMAX (NBUFMXW*sizeof(double))
static double tmpbuf[NBUFMXW];
#endif
#ifdef SHMEM
/*
 * SHMEM  doesn't handle auto & heap vars either. Handle as for BSP.
 * Make sure buffer is a good power-of 2 words and cache aligned.
 */
#pragma _CRI cache_align tmpbuf, pWrk, pSync, pSyncb
#define NBUFMXW 131072
#define NBUFMAX (NBUFMXW*sizeof(double))
static double tmpbuf[NBUFMXW];
/* 
 * SHMEM synchronization and work arrays. pWrk must be bigger than
 * _SHMEM_REDUCE_MIN_WRKDATA_SIZE.  Since we have no compile-time MAX
 * function we declare pSync inelegantly as the sum of the sizes.
 */
static double pWrk[2][NBUFMXW/2+1];
static long pSync[2][_SHMEM_REDUCE_SYNC_SIZE];
static long pSyncb[_SHMEM_BCAST_SYNC_SIZE];
static int  psi = 0;
#endif
#ifdef MPI
#define M_REAL (sizeof(real)==sizeof(double)?MPI_DOUBLE:MPI_FLOAT)
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
#ifdef BSP0
void
par_sigintreset()
{
   signal(SIGINT, SIG_DFL);
}
#endif
#ifdef BSP
void
par_sigintreset()
{
   signal(SIGINT, SIG_DFL);
}
#endif
#ifdef SHMEM
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
void
par_imax(idat)
int *idat;
{
   *idat = bsp_maxI(*idat);
}
#endif
#ifdef BSP0
static void imax(i1, i2, i3, size)
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
#ifdef SHMEM
void
par_imax(idat)
int *idat;
{
   shmem_int_max_to_all(idat, idat, 1, 0, 0, nthreads, (int*)pWrk[psi], pSync[psi]);
   psi = ! psi;
#ifdef SHMEM_BARRIER
   shmem_barrier_all();
#endif
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
 * par_isum().  Calculate sum of int array over all processors.               *
 ******************************************************************************/
#ifdef TCGMSG
void
par_isum(buf, n)
int *buf;
int  n;
{
   IGOP_(ADDR(MSGINT), buf, &n, "+");
}
#endif
#ifdef BSP
static void viadd(res, x, y, nb)
int res[], x[], y[];
int *nb;
{
   int i, n=*nb/sizeof(int);
   for(i = 0; i < n; i++)
      res[i] = x[i] + y[i];
}

void
par_isum(buf, n)
int *buf;
int  n;
{
   bsp_fold(viadd, buf, buf, n*sizeof(int));
}
#endif
#ifdef BSP0
static void viadd(res, x, y, nb)
int res[], x[], y[];
int nb;
{
   int i, n=nb/sizeof(int);
   for(i = 0; i < n; i++)
      res[i] = x[i] + y[i];
}

void
par_isum(buf, n)
int *buf;
int  n;
{
   int m;
   
   /*
    * BSP only allows operations on statically allocated buffers. *sigh*
    * Use loop to perform general operation copying in and out of a
    * fixed-size, static buffer.
    */
   while( n > 0 )
   {
      m = MIN(n, NBUFMAX/sizeof(int));
      memcp(tmpbuf, buf, m*sizeof(int));
      bspreduce(viadd, tmpbuf, tmpbuf, m*sizeof(int));
      memcp(buf, tmpbuf, m*sizeof(int));
      buf += m;
      n -= m;
   }
}
#endif
#ifdef SHMEM
void
par_isum(buf, n)
int *buf;
int  n;
{
   int m;
   
   /*
    * SHMEM only allows operations on statically allocated buffers. *sigh*
    * Use loop to perform general operation copying in and out of a
    * fixed-size, static buffer.
    */
   while( n > 0 )
   {
      m = MIN(n, NBUFMAX/sizeof(int));
      memcp(tmpbuf, buf, m*sizeof(int));
      shmem_int_sum_to_all((int*)tmpbuf, (int*)tmpbuf, m, 0, 0, nthreads, 
			   (int*)pWrk[psi], pSync[psi]);
      memcp(buf, tmpbuf, m*sizeof(int));
      psi = ! psi;
      buf += m;
      n -= m;
#ifdef SHMEM_BARRIER
      shmem_barrier_all();
#endif
   }
}
#endif
#ifdef MPI
/*
 * MPI demands seperate send and receive buffers.  Malloc one and keep
 * it around.  Extend if necessary.
 */
void
par_isum(buf, n)
int *buf;
int  n;
{
   static int *tmpbuf = 0;
   static int  tmpsize = 0;

   if( n <= 0 )
      return;
   if(n > tmpsize)
   {
      if( tmpbuf )
	 xfree(tmpbuf);
      tmpbuf = aalloc(n, int);
      tmpsize = n;
   }
   MPI_Allreduce(buf, tmpbuf, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   memcp(buf, tmpbuf, n*sizeof(int));
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
static void vradd(res, x, y, nb)
real res[], x[], y[];
int *nb;
{
   int i, n=*nb/sizeof(real);
   for(i = 0; i < n; i++)
      res[i] = x[i] + y[i];
}
static void vdadd(res, x, y, nb)
double res[], x[], y[];
int *nb;
{
   int i, n=*nb/sizeof(double);
   for(i = 0; i < n; i++)
      res[i] = x[i] + y[i];
}

void
par_rsum(buf, n)
real *buf;
int  n;
{
   bsp_fold(vradd, buf, buf, n*sizeof(real));
}
void
par_dsum(buf, n)
double *buf;
int  n;
{
   bsp_fold(vradd, buf, buf, n*sizeof(double));
}
#endif
#ifdef BSP0
static void vradd(res, x, y, nb)
real res[], x[], y[];
int nb;
{
   int i, n=nb/sizeof(real);
   for(i = 0; i < n; i++)
      res[i] = x[i] + y[i];
}
static void vdadd(res, x, y, nb)
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
   int m;
   
   /*
    * BSP only allows operations on statically allocated buffers. *sigh*
    * Use loop to perform general operation copying in and out of a
    * fixed-size, static buffer.
    */
   while( n > 0 )
   {
      m = MIN(n, NBUFMAX/sizeof(real));
      memcp(tmpbuf, buf, m*sizeof(real));
      bspreduce(vradd, tmpbuf, tmpbuf, m*sizeof(real));
      memcp(buf, tmpbuf, m*sizeof(real));
      buf += m;
      n -= m;
   }
}
void
par_dsum(buf, n)
double *buf;
int  n;
{
   int m;
   
   /*
    * BSP only allows operations on statically allocated buffers. *sigh*
    * Use loop to perform general operation copying in and out of a
    * fixed-size, static buffer.
    */
   while( n > 0 )
   {
      m = MIN(n, NBUFMAX/sizeof(double));
      memcp(tmpbuf, buf, m*sizeof(double));
      bspreduce(vdadd, tmpbuf, tmpbuf, m*sizeof(double));
      memcp(buf, tmpbuf, m*sizeof(double));
      buf += m;
      n -= m;
   }
}
#endif
#ifdef SHMEM
void
par_rsum(buf, n)
real *buf;
int  n;
{
   int m;
   
   if( sizeof(real) == sizeof(float))
   {
      /*
       * Use loop to perform general operation copying in and out of a
       * fixed-size, static buffer.
       */
      while( n > 0 )
      {
         m = MIN(n, NBUFMAX/sizeof(real));
         memcp(tmpbuf, buf, m*sizeof(real));
         shmem_float_sum_to_all((float*)tmpbuf, (float*)tmpbuf, m,0,0, 
				nthreads, (float*)pWrk[psi], pSync[psi]);
         memcp(buf, tmpbuf, m*sizeof(real));
	 psi = ! psi;
         buf += m;
         n -= m;
#ifdef SHMEM_BARRIER
	 shmem_barrier_all();
#endif
      }
    }
    else if ( sizeof(real) == sizeof(double))
    {
      while( n > 0 )
      {
         m = MIN(n, NBUFMAX/sizeof(real));
         memcp(tmpbuf, buf, m*sizeof(real));
         shmem_double_sum_to_all((double*)tmpbuf, (double*)tmpbuf, m,0,0, 
				 nthreads, (double*)pWrk[psi], pSync[psi]);
         memcp(buf, tmpbuf, m*sizeof(real));
	 psi = ! psi;
         buf += m;
         n -= m;
#ifdef SHMEM_BARRIER
	 shmem_barrier_all();
#endif
      }
    }      
}
void
par_dsum(buf, n)
double *buf;
int  n;
{
   int m;
   
   /*
    * Use loop to perform general operation copying in and out of a
    * fixed-size, static buffer.
    */
   while( n > 0 )
   {
      m = MIN(n, NBUFMAX/sizeof(double));
      memcp(tmpbuf, buf, m*sizeof(double));
      shmem_double_sum_to_all((double*)tmpbuf, (double*)tmpbuf, m, 0, 0, 
			      nthreads, (double*)pWrk[psi], pSync[psi]);
      memcp(buf, tmpbuf, m*sizeof(double));
      psi = ! psi;
      buf += m;
      n -= m;
#ifdef SHMEM_BARRIER
      shmem_barrier_all();
#endif
   }
}
#endif
#ifdef MPI
/*
 * MPI demands seperate send and receive buffers.  Malloc one and keep
 * it around.  Extend if necessary.
 *
 * The MPI standard recommends but does not guarantee that the results 
 * of an Allreduce operation are identical on all threads.  For
 * efficiency we assume that it is here, but provide conditional code
 * in case it isn't. The test on accel.c will trap if results differ,
 * in which case recompile with -DUNSYMM.
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
	 xfree(tmpbuf);
      tmpbuf = dalloc(n);
      tmpsize = n;
   }
#ifdef UNSYMM
   MPI_Reduce(buf, tmpbuf, n, M_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Bcast(tmpbuf, n, M_REAL, 0, MPI_COMM_WORLD);
#else
   MPI_Allreduce(buf, tmpbuf, n, M_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
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
	 xfree(tmpbuf);
      tmpbuf = aalloc(n, double);
      tmpsize = n;
   }
#ifdef UNSYMM
   MPI_Reduce(buf, tmpbuf, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Bcast(tmpbuf, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
   MPI_Allreduce(buf, tmpbuf, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
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
   bsp_push_reg(buf, n*size);
   bsp_bcast(ifrom, buf, buf, n*size);
   bsp_pop_reg(buf);
}
#endif
#ifdef BSP0
void
par_broadcast(buf, n, size, ifrom)
gptr	*buf;
int	n;
size_mt	size;
int	ifrom;
{
   int m;
   long nbyt = n*size;	/* Must have a signed type for loop test */
   
   /*
    * BSP only allows operations on statically allocated buffers. *sigh*
    * Use loop to perform general operation copying in and out of a
    * fixed-size, static buffer.
    */
   while( nbyt > 0 )
   {
      m = MIN(nbyt, NBUFMAX);
      memcp(tmpbuf, buf, m);
      bspbroadcast(ifrom, tmpbuf, tmpbuf, m);
      memcp(buf, tmpbuf, m);
      buf = (char*)buf + m;
      nbyt -= m;
   }
}
#endif
#ifdef SHMEM
void
par_broadcast(buf, n, size, ifrom)
gptr	*buf;
int	n;
size_mt	size;
int	ifrom;
{
   int m;
   long nbyt = n*size;	/* Must have a signed type for loop test */
   
   /*
    * Usual comments about fixed-size buffers apply.  The shmem
    * broadcast routine works in 8 byte word units, but that's OK
    * since we copy the exact length in and out of the real arrays.
    */
   while( nbyt > 0 )
   {
      m = MIN(nbyt, NBUFMAX);
      if( ithread == ifrom )
	 memcp(tmpbuf, buf, m);
      shmem_broadcast64((long*)tmpbuf, (long*)tmpbuf, 
		      (m+sizeof(long)-1)/sizeof(long), 
		      ifrom, 0, 0, nthreads, pSyncb);
      if( ithread != ifrom )
	 memcp(buf, tmpbuf, m);
      buf = (char*)buf + m;
      nbyt -= m;
      shmem_barrier_all();
   }
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
 * par_collect_all().  Global gather to all.                                  *
 ******************************************************************************/
#ifdef BSP
void
par_collect_all(send, recv, n, stride, nblk)
real	*send, *recv;
int	n, nblk, stride;
{
   int  i, right, left, iblk, ibeg;
   
   bsp_push_reg(recv, nblk*stride*sizeof(real));
   /*
    * Copy send data to proper place in receive buffer unless it is
    * already there.
    */
   if( recv+ithread*n != send )
   {
      for(iblk = 0; iblk < nblk; iblk++)
	 memcp(recv+ithread*n+iblk*stride, send+iblk*stride, n*sizeof(real));
   }

   for (i=1; i<nthreads; i*=2) 
   {
      left  = (nthreads + ithread - i) % nthreads;
      right = (ithread + i) % nthreads;
      ibeg = ithread + 1 - i;

      if( ibeg >= 0 )
	 for(iblk = 0; iblk < nblk; iblk++)
	    bsp_hpput(right, recv+ibeg*n+iblk*stride, recv,
                            (ibeg*n+iblk*stride)*sizeof(real), i*n*sizeof(real));
      else
	 for(iblk = 0; iblk < nblk; iblk++)
	 {
	    bsp_hpput(right, recv+(ibeg+nthreads)*n+iblk*stride, recv,
                            ((ibeg+nthreads)*n+iblk*stride)*sizeof(real), 
		            -ibeg*n*sizeof(real));
	    bsp_hpput(right, recv+iblk*stride, recv,
                            (iblk*stride)*sizeof(real), 
                            (ithread+1)*n*sizeof(real));
	 }
      bsp_sync();
   }
   bsp_pop_reg(recv);
}
#endif
#ifdef BSP0
void
par_collect_all(send, recv, n, stride, nblk)
real	*send, *recv;
int	n, nblk, stride;
{
   int  i, right, left, iblk, ibeg;
   real *recvbuf = (real*)tmpbuf;
   int  nbuf = NBUFMAX/sizeof(real);

   if(nblk*stride > nbuf)
      message(NULLI, NULLP, FATAL, "Par_Collect_All: Buffer too small %d,%d\n",
	      nblk*stride,nbuf);
   
   for(iblk = 0; iblk < nblk; iblk++)
      memcp(recvbuf+ithread*n+iblk*stride, send+iblk*stride, n*sizeof(real));
  
   for (i=1; i<nthreads; i*=2) 
   {
      left  = (nthreads + ithread - i) % nthreads;
      right = (ithread + i) % nthreads;
      ibeg = ithread + 1 - i;
      bspsstep(111);

      if( ibeg >= 0 )
	 for(iblk = 0; iblk < nblk; iblk++)
	    bspstore(right, recvbuf+ibeg*n+iblk*stride, 
		            recvbuf+ibeg*n+iblk*stride, i*n*sizeof(real));
      else
	 for(iblk = 0; iblk < nblk; iblk++)
	 {
	    bspstore(right, recvbuf+(ibeg+nthreads)*n+iblk*stride, 
		            recvbuf+(ibeg+nthreads)*n+iblk*stride, 
		            -ibeg*n*sizeof(real));
	    bspstore(right, recvbuf+iblk*stride, 
		            recvbuf+iblk*stride, (ithread+1)*n*sizeof(real));
	 }
      bspsstep_end(111);
   }
   memcp(recv, recvbuf, nblk*stride*sizeof(real));
}
#endif
#ifdef SHMEM
#define NSYNC 16
void
par_collect_all(send, recv, n, stride, nblk)
real	*send, *recv;
int	n, nblk, stride;
{
   int  i, right, left, iblk, ibeg, iSync;
   static long *recvs;
   long        *recvbuf;
   static volatile long SyncB[NSYNC];
                   long Sync0 = ithread;
   long        *cur_break;
   
   for(iSync=0; iSync < NSYNC; iSync++)
      SyncB[iSync] =  _SHMEM_SYNC_VALUE;
   iSync = 0;

   recvs = (long*)recv;
   /*
    * Copy send data to proper place in receive buffer unless it is
    * already there.
    */
   if( recv+ithread*n != send )
   {
      for(iblk = 0; iblk < nblk; iblk++)
	 memcp(recv+ithread*n+iblk*stride, send+iblk*stride, n*sizeof(real));
   }

   shmem_barrier_all();		/* Ensure all PEs ready to receive data */
   shmem_set_cache_inv();
   
   /*
    * This version lives dangerously and works "in-place" with
    * the existing receive buffers.  
    */
   for (i=1; i<nthreads; i*=2) 
   {
      left  = (nthreads + ithread - i) % nthreads;
      right = (ithread + i) % nthreads;
      ibeg = ithread + 1 - i;
      /* 
       * Get Addr of receive buffer on remote PE 
       */
      shmem_long_get((long*)&recvbuf, (long*)&recvs, 1, right);
#ifdef _T3D
      /*
       * The documentation warns that the REMOTE address must not
       * extend beyond allocated memory on the LOCAL PE. Therefore
       * we extend the break if necessary and pray that malloc() can
       * cope. (Actually, we _could_ store it and reset on exit.)
       */
      cur_break = sbreak(0); 
      if( recvbuf+nblk*stride > cur_break )
      {
	 if( sbreak(recvbuf+nblk*stride - cur_break) == 0)
	    message(NULLI, NULLP, FATAL, 
		    "Par_collect_all (%d): sbreak() failed\n");
      }
#endif
      if( ibeg >= 0 )
	 for(iblk = 0; iblk < nblk; iblk++)
	    shmem_long_put(recvbuf+ibeg*n+iblk*stride, 
		      recvs+ibeg*n+iblk*stride, i*n, right);
      else
	 for(iblk = 0; iblk < nblk; iblk++)
	 {
	    shmem_long_put(recvbuf+(ibeg+nthreads)*n+iblk*stride, 
		      recvs+(ibeg+nthreads)*n+iblk*stride, -ibeg*n, right);
	    shmem_long_put(recvbuf+iblk*stride, 
		      recvs+iblk*stride, (ithread+1)*n, right);
	 }
      shmem_quiet();	                   /* Wait for puts to complete  */
      shmem_long_put((long*)&SyncB[iSync], &Sync0, 1, right); 
                                           /* Signal to "right" its done */
      while( SyncB[iSync] == _SHMEM_SYNC_VALUE )
	 ;			/* Spin wait until data has arrived here */
#ifdef DEBUG
      if( SyncB[iSync] != left )
	 message(NULLI, NULLP, FATAL, 
		 "Par Collect: (%d) Sync out of order. Expected %d, got %d\n",
		 ithread, left, SyncB[iSync]);
#endif
      SyncB[iSync] = _SHMEM_SYNC_VALUE;
      iSync++;
   }
   shmem_clear_cache_inv();
}
#endif
#ifdef MPI
void
par_collect_all(send, recv, n, stride, nblk)
real	*send, *recv;
int	n, nblk, stride;
{
   int i;
   int  blens[2];
   MPI_Datatype vtype, block, types[2];
   MPI_Aint displs[2];

   /*
    * Use the defined datatypes of MPI to collect the whole array from
    * distributed slices across processors.  The "vtype" vector defines
    * the actual data, "nblk" blocks of "n" elements with a stride of "stride".
    * "block" is a struct with identical data to "vtype" BUT with a upper
    * bound of exactly 1 block.  This fools MPI_Allgather into assembling
    * the scattered data into an interleaved array.
    */
   MPI_Type_vector(nblk, n, stride, M_REAL, &vtype);
   blens[0]  = 1;     blens[1]  = 1;
   types[0]  = vtype; types[1]  = MPI_UB;
   displs[0] = 0;     displs[1] = n*sizeof(real);
   MPI_Type_struct(2, blens, displs, types, &block);
   MPI_Type_commit(&block);
   
   MPI_Allgather(send, 1, block, recv, 1, block, MPI_COMM_WORLD); 
   MPI_Type_free(&vtype);
   MPI_Type_free(&block);   
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
   bsp_begin(bsp_nprocs());
   *nthreads = bsp_nprocs();
   *ithread  = bsp_pid();
}
#endif
#ifdef BSP0
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
#ifdef SHMEM
void
par_begin(argc, argv, ithread, nthreads)
int	*argc;
char	***argv;
int	*ithread;
int	*nthreads;
{
   int i;
   start_pes(0);
   *nthreads = _num_pes();
   *ithread  = _my_pe();
   
   for(i=0; i < _SHMEM_BCAST_SYNC_SIZE; i++)
     pSyncb[i] = _SHMEM_SYNC_VALUE;
   for(i=0; i < _SHMEM_REDUCE_SYNC_SIZE; i++)
     pSync[0][i] = pSync[1][i] =  _SHMEM_SYNC_VALUE;
   shmem_barrier_all();
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
   bsp_end();
}
#endif
#ifdef BSP0
void
par_finish()
{
   bspfinish();
}
#endif
#ifdef SHMEM
void
par_finish()
{
  shmem_barrier_all();
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
   exit(code);
}
#endif
#ifdef BSP
void
par_abort(code)
int code;
{
   bsp_abort("");
   exit(code);
}
#endif
#ifdef BSP0
void
par_abort(code)
int code;
{
#ifdef NOTYET
   bspabort(code);
#endif
   bspfinish();
   exit(code);
}
#endif
#ifdef SHMEM
#  ifdef HAVE_GLOBALEXIT
void
par_abort(code)
int code;
{
   globalexit(code);
}
#  else
void
par_abort(code)
int code;
{
   kill(getpid(),SIGKILL);
}

#  endif
#endif
#ifdef MPI
void
par_abort(code)
int code;
{
   MPI_Abort(MPI_COMM_WORLD, code);
   exit(code);
}
#endif
/******************************************************************************
 *  copy_sysdef                                                            *
 ******************************************************************************/
/*ARGSUSED*/
void	copy_sysdef(system, spec_ptr, site_info, pot_ptr)
system_mp	system;			/* Pointer to system array (in main)  */
spec_mp		*spec_ptr;		/* Pointer to be set to species array */
site_mp		*site_info;		/* To be pointed at site_info array   */
pot_mp		*pot_ptr;		/* To be pointed at potpar array      */
{
#ifdef SPMD
   spec_mp	spec;
   int		n_pot_recs;
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
      *pot_ptr = aalloc(n_pot_recs, pot_mt);
   par_broadcast((gptr*)*pot_ptr, n_pot_recs, sizeof(pot_mt), 0);
#endif
}
/******************************************************************************
 *  copy_dynamics()							      *
 ******************************************************************************/
/*ARGSUSED*/
void	copy_dynamics(system)
system_mp	system;
{
#ifdef SPMD
   gptr		*ap;			/* Pointer to averages database       */
   size_mt	asize;			/* Size of averages database	      */

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

   par_broadcast((gptr*)system->ta,      system->nspecies, sizeof(real), 0);
   par_broadcast((gptr*)system->tap,     system->nspecies, sizeof(real), 0);
   par_broadcast((gptr*)system->tadot,   system->nspecies, sizeof(real), 0);
   par_broadcast((gptr*)system->tadoto,  system->nspecies, sizeof(real), 0);
   par_broadcast((gptr*)system->tadotvo, system->nspecies, sizeof(real), 0);
					                   	      
   par_broadcast((gptr*)system->ra,      system->nspecies, sizeof(real), 0);
   par_broadcast((gptr*)system->rap,     system->nspecies, sizeof(real), 0);
   par_broadcast((gptr*)system->radot,   system->nspecies, sizeof(real), 0);
   par_broadcast((gptr*)system->radoto,  system->nspecies, sizeof(real), 0);
   par_broadcast((gptr*)system->radotvo, system->nspecies, sizeof(real), 0);

   ap = av_ptr(&asize,0);	      /* get addr, size of database   */
   par_broadcast(ap, 1, asize,0);

   /*
    * N.B. We do NOT broadcast the accumulated RDF info.  That would
    * be incorrect since it is later summed. Leave on thread zero and
    * zeros on other threads.  When globally summed it will be correct.
    */
#endif
}
/******************************************************************************
 * replicate().  Make a copy of all Moldy's constant and dynamic data.        *
 *               This is for parallel implementations and allows "start_up"   *
 *		 to be called on one processor.                               *
 ******************************************************************************/
/*ARGSUSED*/
void replicate(control, system, spec_ptr, site_info, pot_ptr, restart_header)
contr_mt  *control;
system_mt *system;
spec_mt   **spec_ptr;
site_mt	  **site_info;
pot_mt    **pot_ptr;
restrt_mt *restart_header;
{
#ifdef SPMD
   int av_convert;
   int jran;
   /*
    *  Fetch the top-level structs
    */
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
      /*
       * Initialise radial distribution function database.
       */
      if(control->rdf_interval > 0)
         init_rdf(system);
   }
   /*
    * Copy the dynamic vars from the other processes.
    */
   copy_dynamics(system);
   /*
    * Synchronize pseudo-random number generator.
    */
   jran = getseed();
   par_broadcast((gptr*)&jran, 1, sizeof jran, 0);
   smdrand(jran);
#endif
}
