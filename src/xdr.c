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
 * xdr  Moldy-specific xdr routines for storing binary data in machine-       *
 *      independent format.	For compatibility with existing binary 	      *
 *      formats, strings are stored as fixed-length opaque data.	      *
 ******************************************************************************
 *      Revision Log
 *       $Log: xdr.c,v $
 *       Revision 2.10  1996/03/06 18:20:45  keith
 *       Added cast in xdr_vector() call to supress spurious warning message.
 *
 *       Revision 2.9  1994/10/17 10:54:06  keith
 *       Got rid of dummy xdr_array which really screwed things up!
 *
 * Revision 2.8  1994/07/07  17:03:39  keith
 * Fixed up missing xdr_vector to be compiled in only if NEED_XDR_VECTOR defined.
 *
 * Revision 2.7  1994/06/08  13:17:38  keith
 * Changed all timestep-related parameters to type "long". This means
 * that 16-bit DOS compilers can do more than 32767 timesteps.
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Got rid of all global (external) data items except for
 * "control" struct and constant data objects.  The latter
 * (pot_dim, potspec, prog_unit) are declared with CONST
 * qualifier macro which evaluates to "const" or nil
 * depending on ANSI/K&R environment.
 * Also moved as many "write" instantiations of "control"
 * members as possible to "startup", "main" leaving just
 * "dump".
 *
 * Declared as "static"  all functions which should be.
 *
 * Added CONST qualifier to (re-)declarations of ANSI library
 * emulation routines to give reliable compilation even
 * without ANSI_LIBS macro. (#define's away for K&R
 * compilers)
 *
 * Revision 2.5  94/01/18  13:33:07  keith
 * Null update for XDR portability release
 * 
 * Revision 2.4  94/01/18  13:15:18  keith
 * Put casts in function calls to satisfy picky-picky-picky SGI compiler.
 * Added return values to dummy xdr functions to get rid of VMS warnings.
 * 
 * Revision 2.4  93/12/20  16:42:04  keith
 * Put casts in function calls to satisfy picky-picky-picky SGI compiler.
 * 
 * Revision 2.3  93/10/28  10:28:17  keith
 * Corrected declarations of stdargs functions to be standard-conforming
 * 
 * Revision 2.2  93/09/06  14:42:46  keith
 * Fixed portability problems/bugs in XDR code.
 * 
 * Revision 2.1  93/07/19  13:29:08  keith
 * Support for XDR backup/dump routines.
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/users/keith/data/md/moldy/RCS/xdr.c,v 2.11 1996/10/19 11:53:43 keith Exp $";
#endif
/*========================== program include files ===========================*/
#include	"structs.h"
/*========================== Library include files ===========================*/
#include	"stddef.h"
#include 	"xdr.h"
/*============================================================================*/

#ifdef USE_XDR

bool_t xdr_real(xdrs, rp)
XDR     *xdrs;
real    *rp;
{
   if( sizeof(real) == sizeof(double) )
      return xdr_double(xdrs, (double*)rp);
   else if( sizeof(real) == sizeof(float) )
      return xdr_float(xdrs, (float*)rp);
   else
      return FALSE;
}

bool_t xdr_contr(xdrs, cp)
XDR      *xdrs;
contr_mt *cp;
{
   return
      xdr_opaque(xdrs, cp->title, L_name) &&
      xdr_long(xdrs, &cp->istep) &&
      xdr_long(xdrs, &cp->nsteps) &&
      xdr_double(xdrs, &cp->step) &&
      xdr_vector(xdrs, (gptr*)&cp->print_sysdef, 7, sizeof(boolean), (xdrproc_t)xdr_bool) &&
      xdr_opaque(xdrs, cp->sysdef, 6*L_name) &&
      xdr_vector(xdrs, (gptr*)cp->spare, 23, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_double(xdrs, &cp->ttmass) &&
      xdr_double(xdrs, &cp->rtmass) &&
      xdr_int(xdrs, &cp->pad) &&
      xdr_int(xdrs, &cp->const_temp) &&
      xdr_bool(xdrs, &cp->xdr_write) &&
      xdr_bool(xdrs, &cp->strict_cutoff) &&
      xdr_int(xdrs, &cp->strain_mask) &&
      xdr_int(xdrs, &cp->nbins) &&
      xdr_u_long(xdrs, &cp->seed) &&
      xdr_vector(xdrs, (gptr*)&cp->page_width, 2, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_vector(xdrs, (gptr*)&cp->scale_interval, 7, sizeof(long), (xdrproc_t)xdr_long) &&
      xdr_vector(xdrs, (gptr*)&cp->dump_level, 2, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_vector(xdrs, (gptr*)&cp->backup_interval, 6, sizeof(long), (xdrproc_t)xdr_long) &&
      xdr_vector(xdrs, (gptr*)&cp->temp, 10, sizeof(double), (xdrproc_t)xdr_double);
}

bool_t xdr_system(xdrs, sp)
XDR      *xdrs;
system_mt *sp;
{
   return
      xdr_vector(xdrs, (gptr*)&sp->nsites, 8, sizeof(int), (xdrproc_t)xdr_int) &&
      /*
       * This is an awful hack.  There are 28 real[3]* pointers
       * next.  Their stored values are NEVER re-used so we just
       * output a placeholder.  For compatibility of XDR/non-XDR
       * files on 4 byte big-endian ieee architectures we emit
       * 4 bytes each.  DON'T use sizeof as that would make XDR
       * file M/C dependent.
       */
      xdr_opaque(xdrs, (gptr*)&sp->c_of_m, 28*XDR_4PTR_SIZE);
}

/*
 * This version for reading restart files written by 2.10 or before.
 */
bool_t xdr_system_2(xdrs, sp)
XDR      *xdrs;
system_mt *sp;
{
   return
      xdr_vector(xdrs, (gptr*)&sp->nsites, 8, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_opaque(xdrs, (gptr*)&sp->c_of_m, 18*XDR_4PTR_SIZE);
}

bool_t xdr_species(xdrs, sp)
XDR      *xdrs;
spec_mt *sp;
{
   return
      xdr_vector(xdrs, (gptr*)sp->inertia, 6, sizeof(real), (xdrproc_t)xdr_real) &&
      xdr_vector(xdrs, (gptr*)&sp->nsites, 4, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_opaque(xdrs, sp->name, 32) &&
      /*
       * This is an awful hack.  There are 14 real[3]* pointers
       * next.  Their stored values are NEVER re-used so we just
       * output a placeholder.  For compatibility of XDR/non-XDR
       * files on 4 byte big-endian ieee architectures we emit
       * 4 bytes each.  DON'T use sizeof as that would make XDR
       * file M/C dependent.
       */
      xdr_opaque(xdrs, (gptr*)&sp->site_id, 14*XDR_4PTR_SIZE) &&
      xdr_int(xdrs,(int*)&sp->pad[0]) &&
      xdr_int(xdrs,(int*)&sp->pad[1]);
}

bool_t xdr_site(xdrs, sp)
XDR      *xdrs;
site_mt  *sp;
{
   return
      xdr_double(xdrs, &sp->mass) &&
      xdr_double(xdrs, &sp->charge) &&
      xdr_opaque(xdrs, sp->name, 8) &&
      xdr_int(xdrs, &sp->flag) &&
      xdr_int(xdrs, &sp->pad);
}

static unsigned int xdr_npotpar;

void xdr_set_npotpar(npotpar)
int	npotpar;
{
   xdr_npotpar = npotpar;
}

bool_t xdr_pot(xdrs, sp)
XDR      *xdrs;
pot_mt   *sp;
{
   return
      xdr_int(xdrs, &sp->flag) &&
      xdr_int(xdrs, &sp->pad) &&
      xdr_vector(xdrs, (gptr*)sp->p, xdr_npotpar, sizeof(real), (xdrproc_t)xdr_real);
}


bool_t xdr_restrt(xdrs, sp)
XDR         *xdrs;
restrt_mt   *sp;
{
   return
      xdr_u_long(xdrs, &sp->timestamp) &&
      xdr_u_long(xdrs, &sp->prev_timestamp) &&
      xdr_opaque(xdrs, sp->init_date, DLEN+L_name+16) &&
      xdr_int(xdrs, &sp->seq);
}

bool_t xdr_dump(xdrs, sp)
XDR         *xdrs;
dump_mt   *sp;
{
   return
      xdr_opaque(xdrs, sp->title, L_name+16) &&
      xdr_vector(xdrs, (gptr*)&sp->istep, 2, sizeof(long), (xdrproc_t)xdr_long) &&
      xdr_vector(xdrs, (gptr*)&sp->dump_level, 4, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_vector(xdrs, (gptr*)&sp->timestamp, 3, sizeof(unsigned long), (xdrproc_t)xdr_u_long);
}

static
bool_t xdr_av_head_t(xdrs,ap)
XDR	  *xdrs;
av_head_mt *ap;
{
   return
      xdr_vector(xdrs, (gptr*)&ap->nav, 4, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_double(xdrs, &ap->align);
}

static
bool_t xdr_old_av_u_t(xdrs,ap)
XDR	  *xdrs;
old_av_u_mt *ap;
{
   return
      xdr_vector(xdrs, (gptr*)&ap->cnt.av, 2, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_opaque(xdrs, (gptr*)(&ap->cnt.av+2*sizeof(int)), 
		 (7+MAX_ROLL_INTERVAL)*XDR_DOUBLE_SIZE-2*XDR_INT_SIZE);
}

static size_mt  xdr_av_size;
static int    av_convert;

void   xdr_set_av_size_conv(size, av_conv)
size_mt	   size;
int	   av_conv;
{
   xdr_av_size = size;
   av_convert = av_conv;
}

bool_t xdr_averages(xdrs, ap)
XDR	   *xdrs;
gptr	   *ap;
{
   /*
    * The global flag av_convert requires explanation.  It is
    * set by init_averages() which is always called before
    * this function.  av_convert=1 if we are reading data in
    * the old (fixed length = MAX_ROLL_INTERVAL) format and
    * 0 or 2 otherwise. 
    */
   unsigned    navst;		/* Total # of average data items */
   old_av_u_mt *apo = (old_av_u_mt *)ap;
   av_head_mt  *aph = (av_head_mt  *)ap;

   if( xdrs->x_op == XDR_DECODE && av_convert == 1 )
   {
      navst = (xdr_av_size-sizeof(old_av_u_mt)) / sizeof(double);
      return
	 xdr_old_av_u_t(xdrs,apo) &&
	 xdr_vector(xdrs, (gptr*)(apo+1), navst, sizeof(double), (xdrproc_t)xdr_double);
      }
   else
   {
      navst = (xdr_av_size-sizeof(av_head_mt)) / sizeof(double);
      return
	 xdr_av_head_t(xdrs,aph) &&
	 xdr_vector(xdrs, (gptr*)(aph+1), navst, sizeof(double), (xdrproc_t)xdr_double);
   }
}

#ifdef NEED_XDR_VECTOR
#define LASTUNSIGNED	((u_int)0-1)
/*
 * xdr_vector():
 *
 * XDR a fixed length array. Unlike variable-length arrays,
 * the storage of fixed length arrays is static and unfreeable.
 * > basep: base of the array
 * > size: size of the array
 * > elemsize: size of each element
 * > xdr_elem: routine to XDR each element
 */
bool_t
xdr_vector(xdrs, basep, nelem, elemsize, xdr_elem)
	register XDR *xdrs;
	register char *basep;
	register u_int nelem;
	register u_int elemsize;
	register xdrproc_t xdr_elem;	
{
	register u_int i;
	register char *elptr;

	elptr = basep;
	for (i = 0; i < nelem; i++) {
		if (! (*xdr_elem)(xdrs, elptr, LASTUNSIGNED)) {
			return(FALSE);
		}
		elptr += elemsize;
	}
	return(TRUE);	
}
#endif
#else
void	xdr_set_npotpar (npotpar) int npotpar; {}
void	xdr_set_av_size_conv (size, av_conv) size_mt size; int av_conv; {}
bool_t	xdr_site () {return 0;}
bool_t	xdr_restrt () {return 0;}
bool_t	xdr_averages () {return 0;}
bool_t	xdr_real () {return 0;}
bool_t	xdr_contr () {return 0;}
bool_t	xdr_system () {return 0;}
bool_t	xdr_system_2 () {return 0;}
bool_t	xdr_species () {return 0;}
bool_t	xdr_pot () {return 0;}
bool_t	xdr_int () {return 0;}
bool_t	xdr_bool () {return 0;}
#endif
