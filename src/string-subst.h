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
/*
 * $Header: /home/eeyore/keith/md/moldy/RCS/string.h,v 2.1 93/07/19 13:28:19 keith Exp $
 *
 * $Log:	string.h,v $
 * Revision 2.1  93/07/19  13:28:19  keith
 * Added XDR capability for backup and dump files.
 * 
 * Revision 2.0  93/03/15  14:49:29  keith
 * Added copyright notice and disclaimer to apply GPL
 * to all modules. (Previous versions licensed by explicit 
 * consent only).
 * 
 * Revision 1.10  93/03/09  15:59:26  keith
 * Changed all *_t types to *_mt for portability.
 * Reordered header files for GNU CC compatibility.
 * 
 * Revision 1.9  91/08/16  15:26:17  keith
 * Checked error returns from fread, fwrite, fseek and fclose more
 * rigourously.   Called strerror() to report errors.
 * 
 * Revision 1.8  91/08/15  18:12:24  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * Revision 1.7  91/03/12  15:43:40  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.6  90/04/14  16:02:06  keith
 * Replaced calls to system headers with own definitions for portability.
 * 
 * 
 */
#ifndef STRING_ALREADY
#define STRING_ALREADY

#ifdef ANSI_LIBS
#   include <string.h>
#else

#include "stddef.h"

extern char
	*strcat(),
        *strcpy(),
        *strncpy(),
	*strchr(),
	*strtok(),
	*strdup(),
	*strstr(),
        *strerror();

extern int
	strcmp();
#ifdef __GNUC__
extern size_t
#else
extern int
#endif
	strlen();

extern gptr
	*memcpy(),
	*memset();

#endif

#endif
