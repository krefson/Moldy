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
 *       $Log: xdr.h,v $
 *       Revision 2.9  1996/09/25 16:29:12  keith
 *       Fixed restart structure correctly - broken in prev version.
 *       Thermostat parameters may not be properly read.
 *
 *       Revision 2.8  1996/03/07 15:01:57  keith
 *       Made "malloc in types.h" macro protection conditional on "ANSI_LIBS"
 *
 *       Revision 2.7  1994/06/08 13:22:31  keith
 *       Null update for version compatibility
 *
 * Revision 2.6  1994/02/17  16:38:16  keith
 * Significant restructuring for better portability and
 * data modularity.
 *
 * Revision 2.5.1.1  1994/02/03  18:36:12  keith
 * Tidied up and got rid of most of the global data items.
 *
 * Revision 2.5  94/01/18  17:35:45  keith
 * Incorporated all portability experience to multiple platforms since 2.2.
 * Including ports to VAX/VMS and Open VMS on Alpha AXP and Solaris.
 * 
 * Revision 2.5  94/01/18  13:33:08  keith
 * Null update for XDR portability release
 * 
 * Revision 2.4  94/01/18  13:23:19  keith
 * Incorporated all portability experience to multiple platforms since 2.2.
 * Including ports to VAX/VMS and Open VMS on Alpha AXP and Solaris.
 * 
 * Revision 2.3  93/10/14  18:18:16  keith
 * Fixed prortability problems to IBM RS6000
 * 
 * Revision 2.2  93/09/06  14:42:47  keith
 * Fixed portability problems/bugs in XDR code.
 * 
 * Revision 2.1  93/07/19  13:29:47  keith
 * Support for XDR backup/dump routines.
 * 
 */
#ifndef lint
static char *RCSidh = "$Header: /home/eeyore_data/keith/md/moldy/RCS/xdr.h,v 2.9 1996/09/25 16:29:12 keith Exp $";
#endif
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

#ifdef vms
#include	"rpc_types.h"
#include	"rpc_xdr.h"
#else
#include	"time.h"
#include	<rpc/types.h>
#include	<rpc/xdr.h>
#endif


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

bool_t xdr_real();
bool_t xdr_contr();
bool_t xdr_system();
bool_t xdr_system_2();
bool_t xdr_species();
bool_t xdr_site();
void   xdr_set_npotpar();
bool_t xdr_pot();
bool_t xdr_restrt();
bool_t xdr_dump();
void   xdr_set_av_size_conv();
bool_t xdr_averages();

#ifndef USE_XDR
bool_t	xdr_int();
bool_t  xdr_bool();
#endif

#define XDR_INT_SIZE 4
#define XDR_4PTR_SIZE 4
#define XDR_ULONG_SIZE 4
#define XDR_FLOAT_SIZE 4
#define XDR_DOUBLE_SIZE 8
#define XDR_REAL_SIZE ( (sizeof(real)==sizeof(double))?XDR_DOUBLE_SIZE:XDR_FLOAT_SIZE)

#define XDR_RESTRT_SIZE  (2*XDR_ULONG_SIZE+(DLEN)+(L_name)+16+XDR_INT_SIZE)
#define XDR_DUMP_SIZE    ((L_name)+16+6*XDR_INT_SIZE+3*XDR_ULONG_SIZE)
