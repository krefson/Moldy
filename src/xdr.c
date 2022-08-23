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
 */
/*========================== program include files ===========================*/
#include	"structs.h"
/*========================== Library include files ===========================*/
#include	<stddef.h>
#include 	"xdr.h"
/*============================================================================*/

#ifdef USE_XDR

bool_t xdr_real(XDR *xdrs, real *rp)
{
   /*CONSTCOND*/
   if( sizeof(real) == sizeof(double) )
      return xdr_double(xdrs, (double*)rp);
   /*CONSTCOND*/
   else if( sizeof(real) == sizeof(float) )
      return xdr_float(xdrs, (float*)rp);
   else
      return FALSE;
}

bool_t xdr_contr(XDR *xdrs, contr_mt *cp)
{
   return
      xdr_opaque(xdrs, cp->title, L_name) &&
      xdr_long(xdrs, &cp->istep) &&
      xdr_long(xdrs, &cp->nsteps) &&
      xdr_double(xdrs, &cp->step) &&
      xdr_vector(xdrs, (gptr*)&cp->print_sysdef, 4, sizeof(boolean), 
		 (xdrproc_t)xdr_bool) &&
      xdr_int(xdrs, &cp->scale_options) &&
      xdr_bool(xdrs, &cp->surface_dipole) &&
      xdr_bool(xdrs, &cp->lattice_start) &&
      xdr_opaque(xdrs, cp->sysdef, 6*L_name) &&
      xdr_vector(xdrs, (gptr*)cp->spare, 20, sizeof(int), 
		 (xdrproc_t)xdr_int) &&
      xdr_bool(xdrs, &cp->nosymmetric_rot) &&
      xdr_double(xdrs, &cp->ewald_accuracy) &&
      xdr_double(xdrs, &cp->ttmass) &&
      xdr_double(xdrs, &cp->rtmass) &&
      xdr_int(xdrs, &cp->const_pressure) &&
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

bool_t xdr_system(XDR *xdrs, system_mt *sp)
{
   return
      xdr_vector(xdrs, (gptr*)&sp->nsites, 8, sizeof(int), (xdrproc_t)xdr_int) &&
      /*
       * This is an awful hack.  There are 9 real[3]* pointers
       * next.  Their stored values are NEVER re-used so we just
       * output a placeholder.  For compatibility of XDR/non-XDR
       * files on 4 byte big-endian ieee architectures we emit
       * 4 bytes each.  DON'T use sizeof as that would make XDR
       * file M/C dependent.
       */
      xdr_opaque(xdrs, (gptr*)&sp->c_of_m, 9*XDR_4PTR_SIZE) &&
      xdr_vector(xdrs, (gptr*)&sp->ts, 5, sizeof(real),(xdrproc_t)xdr_real);
}

bool_t xdr_species(XDR *xdrs, spec_mt *sp)
{
   return
      xdr_vector(xdrs, (gptr*)sp->inertia, 6, sizeof(real), (xdrproc_t)xdr_real) &&
      xdr_vector(xdrs, (gptr*)&sp->nsites, 4, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_opaque(xdrs, sp->name, L_spec) &&
      /*
       * This is an awful hack.  There are 8 real[3]* pointers
       * next.  Their stored values are NEVER re-used so we just
       * output a placeholder.  For compatibility of XDR/non-XDR
       * files on 4 byte big-endian ieee architectures we emit
       * 4 bytes each.  DON'T use sizeof as that would make XDR
       * file M/C dependent.
       */
      xdr_opaque(xdrs, (gptr*)&sp->site_id, 8*XDR_4PTR_SIZE);
}

bool_t xdr_site(XDR *xdrs, site_mt *sp)
{
   return
      xdr_double(xdrs, &sp->mass) &&
      xdr_double(xdrs, &sp->charge) &&
      xdr_opaque(xdrs, sp->name, L_site) &&
      xdr_int(xdrs, &sp->flag) &&
      xdr_int(xdrs, &sp->pad);
}

static unsigned int xdr_npotpar;

void xdr_set_npotpar(int npotpar)
{
   xdr_npotpar = npotpar;
}

bool_t xdr_pot(XDR *xdrs, pot_mt *sp)
{
   return
      xdr_int(xdrs, &sp->flag) &&
      xdr_int(xdrs, &sp->pad) &&
      xdr_vector(xdrs, (gptr*)sp->p, xdr_npotpar, sizeof(real), (xdrproc_t)xdr_real);
}


bool_t xdr_restrt(XDR *xdrs, restrt_mt *sp)
{
   return
      xdr_u_long(xdrs, &sp->timestamp) &&
      xdr_u_long(xdrs, &sp->prev_timestamp) &&
      xdr_opaque(xdrs, sp->init_date, DLEN+L_name+L_vsn) &&
      xdr_int(xdrs, &sp->seq);
}

bool_t xdr_dump(XDR *xdrs, dump_mt *sp)
{
   return
      xdr_opaque(xdrs, sp->title, L_name+L_vsn) &&
      xdr_vector(xdrs, (gptr*)&sp->istep, 2, sizeof(long), (xdrproc_t)xdr_long) &&
      xdr_vector(xdrs, (gptr*)&sp->dump_level, 4, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_vector(xdrs, (gptr*)&sp->timestamp, 4, sizeof(unsigned long), (xdrproc_t)xdr_u_long);
}

bool_t xdr_dump_sysinfo_hdr(XDR *xdrs, dump_sysinfo_mt *sp)
{
   int i, ispec;
   i = 
      xdr_float(xdrs, &sp->deltat) &&
      xdr_vector(xdrs, (gptr*)&sp->nmols, 3, sizeof(int), (xdrproc_t)xdr_int);
   return i;
}

bool_t xdr_dump_sysinfo_2_23(XDR *xdrs, dump_sysinfo_mt *sp)
{
   int i, ispec;
   i = 
      xdr_float(xdrs, &sp->deltat) &&
      xdr_vector(xdrs, (gptr*)&sp->nmols, 3, sizeof(int), (xdrproc_t)xdr_int);
   for(ispec = 0; ispec < sp->nspecies; ispec++)
   {
      i = i && xdr_vector(xdrs, (gptr*)&sp->mol[ispec].inertia, 6, sizeof(float), (xdrproc_t)xdr_float);
      i = i && xdr_vector(xdrs, (gptr*)&sp->mol[ispec].nmols, 3, sizeof(int), (xdrproc_t)xdr_int);
      i = i && xdr_opaque(xdrs, (gptr*)&sp->mol[ispec].name, L_spec);
   }
   return i;
}


bool_t xdr_dump_sysinfo_2_18(XDR *xdrs, dump_sysinfo_mt *sp)
{
   int i, ispec;
   i = 
      xdr_float(xdrs, &sp->deltat) &&
      xdr_vector(xdrs, (gptr*)&sp->nmols, 3, sizeof(int), (xdrproc_t)xdr_int);
   for(ispec = 0; ispec < sp->nspecies; ispec++)
   {
      i = i && xdr_opaque(xdrs, (gptr*)&sp->mol[ispec].name, L_spec);
      i = i && xdr_vector(xdrs, (gptr*)&sp->mol[ispec].nmols, 2, sizeof(int), (xdrproc_t)xdr_int);
   }
   return i;
}

bool_t xdr_dump_sysinfo(XDR *xdrs, dump_sysinfo_mt *sp, 
			     int vmajor, int vminor)
{
   if( vmajor >= 2 )
   {
      if(vminor >= 23) 
	 return xdr_dump_sysinfo_2_23(xdrs, sp);
      else if(vminor >= 18) 
	 return xdr_dump_sysinfo_2_18(xdrs, sp);
   }
   return false;
}

static
bool_t xdr_av_head_t(XDR *xdrs, av_head_mt *ap)
{
   return
      xdr_vector(xdrs, (gptr*)&ap->nav, 4, sizeof(int), (xdrproc_t)xdr_int) &&
      xdr_double(xdrs, &ap->align);
}

static size_mt  xdr_av_size;
static int    av_convert;

void   xdr_set_av_size_conv(size_mt size, int av_conv)
{
   xdr_av_size = size;
   av_convert = av_conv;
}

bool_t xdr_averages(XDR *xdrs, gptr *ap)
{
   unsigned    navst;		/* Total # of average data items */
   av_head_mt  *aph = (av_head_mt  *)ap;

   navst = (xdr_av_size-sizeof(av_head_mt)) / sizeof(double);
   return
      xdr_av_head_t(xdrs,aph) &&
      xdr_vector(xdrs, (gptr*)(aph+1), navst, sizeof(double), (xdrproc_t)xdr_double);
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

/*ARGSUSED*/
bool_t xdr_averages(XDR *xdrs, gptr *ap){return 0;}
/*ARGSUSED*/
bool_t xdr_contr(XDR *xdrs, contr_mt *cp){return 0;}
/*ARGSUSED*/
bool_t xdr_dump(XDR *xdrs, dump_mt *sp){return 0;}
/*ARGSUSED*/
bool_t xdr_dump_sysinfo(XDR *xdrs, dump_sysinfo_mt *sp, int vmajor, int vminor)
{return 0;}
/*ARGSUSED*/
bool_t xdr_pot(XDR *xdrs, pot_mt *sp){return 0;}
/*ARGSUSED*/
bool_t xdr_real(XDR *xdrs, real *rp){return 0;}
/*ARGSUSED*/
bool_t xdr_restrt(XDR *xdrs, restrt_mt *sp){return 0;}
/*ARGSUSED*/
bool_t xdr_site(XDR *xdrs, site_mt *sp){return 0;}
/*ARGSUSED*/
bool_t xdr_species(XDR *xdrs, spec_mt *sp){return 0;}
/*ARGSUSED*/
bool_t xdr_system(XDR *xdrs, system_mt *sp){return 0;}
/*ARGSUSED*/
bool_t xdr_system_2(XDR *xdrs, system_mt *sp){return 0;}
/*ARGSUSED*/
void   xdr_set_av_size_conv(size_mt size, int av_conv){return ;}
/*ARGSUSED*/
void   xdr_set_npotpar(int npotpar){return ;}

/*ARGSUSED*/
bool_t	xdr_bool (void) {return 0;}
/*ARGSUSED*/
bool_t	xdr_int (void) {return 0;}
#endif
