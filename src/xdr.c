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
 *       $Log:	xdr.c,v $
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
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/xdr.c,v 2.4 93/12/20 16:42:04 keith Exp $";
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
      xdr_int(xdrs, &cp->istep) &&
      xdr_int(xdrs, &cp->nsteps) &&
      xdr_double(xdrs, &cp->step) &&
      xdr_vector(xdrs, (gptr*)&cp->print_sysdef, 7, sizeof(boolean), xdr_bool) &&
      xdr_opaque(xdrs, cp->sysdef, 6*L_name) &&
      xdr_vector(xdrs, (gptr*)cp->spare, 29, sizeof(int), xdr_int) &&
      xdr_bool(xdrs, &cp->xdr_write) &&
      xdr_bool(xdrs, &cp->strict_cutoff) &&
      xdr_int(xdrs, &cp->strain_mask) &&
      xdr_int(xdrs, &cp->nbins) &&
      xdr_u_long(xdrs, &cp->seed) &&
      xdr_vector(xdrs, (gptr*)&cp->page_width, 17, sizeof(int), xdr_int) &&
      xdr_vector(xdrs, (gptr*)&cp->temp, 10, sizeof(double), xdr_double);
}

bool_t xdr_system(xdrs, sp)
XDR      *xdrs;
system_mt *sp;
{
   return
      xdr_vector(xdrs, (gptr*)&sp->nsites, 8, sizeof(int), xdr_int) &&
      /*
       * This is an awful hack.  There are 18 real[3]* pointers
       * next.  Their stored values are NEVER re-used so we just
       * output a placeholder.  For compatibility of XDR/non-XDR
       * files on 4 byte big-endian ieee architectures we emit
       * 4 bytes each.  DON'T use sizeof as that would make XDR
       * file M/C dependent.
       */
      xdr_opaque(xdrs, (gptr*)&sp->c_of_m, 18*XDR_4PTR_SIZE);
}

bool_t xdr_species(xdrs, sp)
XDR      *xdrs;
spec_mt *sp;
{
   return
      xdr_vector(xdrs, (gptr*)sp->inertia, 6, sizeof(real), xdr_real) &&
      xdr_vector(xdrs, (gptr*)&sp->nsites, 4, sizeof(int), xdr_int) &&
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
      xdr_vector(xdrs, (gptr*)sp->p, xdr_npotpar, sizeof(real), xdr_real);
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
      xdr_vector(xdrs, (gptr*)&sp->istep, 6, sizeof(int), xdr_int) &&
      xdr_vector(xdrs, (gptr*)&sp->timestamp, 3, sizeof(unsigned long), xdr_u_long);
}

bool_t xdr_av_head_t(xdrs,ap)
XDR	  *xdrs;
av_head_mt *ap;
{
   return
      xdr_vector(xdrs, (gptr*)&ap->nav, 4, sizeof(int), xdr_int) &&
      xdr_double(xdrs, &ap->align);
}

bool_t xdr_old_av_u_t(xdrs,ap)
XDR	  *xdrs;
old_av_u_mt *ap;
{
   return
      xdr_vector(xdrs, (gptr*)&ap->cnt.av, 2, sizeof(int), xdr_int) &&
      xdr_opaque(xdrs, (gptr*)(&ap->cnt.av+2*sizeof(int)), 
		 (7+MAX_ROLL_INTERVAL)*XDR_DOUBLE_SIZE-2*XDR_INT_SIZE);
}

static size_t xdr_av_size;

void   xdr_set_av_size(size)
size_t	   size;
{
   xdr_av_size = size;
}

bool_t xdr_averages(xdrs, ap)
XDR	   *xdrs;
gptr	   *ap;
{
   /*
    * The following global requires explanation.  It is
    * set by init_averages() which is always called before
    * this function.  av_convert=1 if we are reading data in
    * the old (fixed length = MAX_ROLL_INTERVAL) format and
    * 0 or 2 otherwise. 
    */
   extern int av_convert;
   unsigned    navst;		/* Total # of average data items */
   old_av_u_mt *apo = (old_av_u_mt *)ap;
   av_head_mt  *aph = (av_head_mt  *)ap;

   if( xdrs->x_op == XDR_DECODE && av_convert == 1 )
   {
      navst = (xdr_av_size-sizeof(old_av_u_mt)) / sizeof(double);
      return
	 xdr_old_av_u_t(xdrs,apo) &&
	 xdr_vector(xdrs, (gptr*)(apo+1), navst, sizeof(double), xdr_double);
      }
   else
   {
      navst = (xdr_av_size-sizeof(av_head_mt)) / sizeof(double);
      return
	 xdr_av_head_t(xdrs,aph) &&
	 xdr_vector(xdrs, (gptr*)(aph+1), navst, sizeof(double), xdr_double);
   }
}

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

#else
void	xdr_set_npotpar () {}
void	xdr_set_av_size () {}
bool_t	xdr_site () {return 0;}
bool_t	xdr_restrt () {return 0;}
bool_t	xdr_averages () {return 0;}
bool_t	xdr_real () {return 0;}
bool_t	xdr_contr () {return 0;}
bool_t	xdr_system () {return 0;}
bool_t	xdr_species () {return 0;}
bool_t	xdr_pot () {return 0;}
bool_t	xdr_int () {return 0;}
bool_t	xdr_bool () {return 0;}
#endif
