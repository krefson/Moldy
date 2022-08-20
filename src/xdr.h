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
 */
/*========================== Library include files ===========================*/
#ifdef USE_XDR

/*
 * Some <rpc/types.h> descended from Sun's original include a declaration
 * of "malloc" in an unprotected fashion.  Try to define it out of the
 * way -- include "stdlib.h" if necessary to put it back.
 * In case an implementation (eg SGI) does it right by including <stdlib.h>
 * ensure that any Moldy module includes "stdlib.h" *before* "xdr.h".
*/
#ifndef STDC_HEADERS
#define free xxfree
#define exit xxexit
#define malloc xxmalloc
#define calloc xxcalloc
#define realloc xxrealloc
#endif

/*
 * A Horrible hack.  defs.h declares MIN and MAX macros, but so does
 * <rpc/types.h>.  Undefine and redefine them here.
 */

#undef MIN
#undef MAX

#ifdef vms
#include	"rpc_types.h"
#include	"rpc_xdr.h"
#else
#include	"time.h"
#include	<rpc/types.h>
#include	<rpc/xdr.h>
#endif

#ifdef MIN
#undef  MIN
#endif
#ifdef MAX
#undef  MAX
#endif
#define MIN(x,y)	((x) < (y) ? (x) : (y))
#define MAX(x,y)	((x) > (y) ? (x) : (y))


#ifndef STDC_HEADERS
#undef free
#undef exit
#undef malloc
#undef calloc
#undef realloc
#endif

#else
typedef	char XDR;
typedef int bool_t;
typedef bool_t (*xdrproc_t)();
#endif
/*============================================================================*/

bool_t xdr_real(XDR *xdrs, real *rp);
bool_t xdr_contr(XDR *xdrs, contr_mt *cp);
bool_t xdr_system(XDR *xdrs, system_mt *sp);
bool_t xdr_system_2(XDR *xdrs, system_mt *sp);
bool_t xdr_species(XDR *xdrs, spec_mt *sp);
bool_t xdr_site(XDR *xdrs, site_mt *sp);
void   xdr_set_npotpar(int npotpar);
bool_t xdr_pot(XDR *xdrs, pot_mt *sp);
bool_t xdr_restrt(XDR *xdrs, restrt_mt *sp);
bool_t xdr_dump(XDR *xdrs, dump_mt *sp);
bool_t xdr_dump_sysinfo_hdr(XDR *xdrs, dump_sysinfo_mt *sp);
bool_t xdr_dump_sysinfo(XDR *xdrs, dump_sysinfo_mt *sp, int vmajor, int vminor);
void   xdr_set_av_size_conv(size_mt size, int av_conv);
bool_t xdr_averages(XDR *xdrs, gptr *ap);

#ifndef USE_XDR
bool_t	xdr_int(void);
bool_t  xdr_bool(void);
#endif

#define XDR_INT_SIZE 4
#define XDR_4PTR_SIZE 4
#define XDR_ULONG_SIZE 4
#define XDR_FLOAT_SIZE 4
#define XDR_DOUBLE_SIZE 8
#define XDR_REAL_SIZE ( (sizeof(real)==sizeof(double))?XDR_DOUBLE_SIZE:XDR_FLOAT_SIZE)

#define XDR_RESTRT_SIZE  (2*XDR_ULONG_SIZE+(DLEN)+(L_name)+L_vsn+XDR_INT_SIZE)
#define XDR_DUMP_SIZE    ((L_name)+L_vsn+6*XDR_INT_SIZE+4*XDR_ULONG_SIZE)
#define XDR_SYSINFO_SIZE_PRE22(nspecies) (XDR_FLOAT_SIZE+(3+2*nspecies)*XDR_INT_SIZE + 32*nspecies)
#define XDR_SYSINFO_SIZE(nspecies) (XDR_FLOAT_SIZE+3*XDR_INT_SIZE+\
				    nspecies*( 6*XDR_FLOAT_SIZE\
					      +3*XDR_INT_SIZE\
					      +L_spec))
