/*
 * $Header: /home/eeyore/keith/md/moldy/RCS/string.h,v 1.5 90/03/26 18:05:18 keith Exp $
 *
 * $Log:	string.h,v $
 * Revision 1.4  90/03/09  17:33:52  keith
 * Modified conditionals for unicos.
 * 
 * Revision 1.3  89/11/20  12:05:01  keith
 * Added 'strdup'.
 * 
 * Revision 1.2  89/06/09  12:18:51  keith
 * Recognised sysV and SysV macros
 * 
 * Revision 1.1  89/04/25  17:38:05  keith
 * Initial revision
 * 
 * 
 */
#ifndef STRING_ALREADY
#define STRING_ALREADY

#include "defs.h"		/* Defs.h sets USG symbol for sysV	      */

#ifdef unix

#ifdef USG
#	include	<string.h>
#	include	<memory.h>
#else  /* BSD */
#	define strchr(c,p) index(c,p)
#	include	<strings.h>
	char		*memcpy();
	char		*strdup();
#endif

#else	/* Not unix */

#ifndef CRAY
#	include	<string.h>
#endif

#endif /* unix      */

#endif
