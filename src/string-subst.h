/*
 * $Header: /home/eeyore/keith/md/moldy/RCS/string.h,v 1.8 91/08/14 14:24:00 keith Exp $
 *
 * $Log:	string.h,v $
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
	*strdup();
extern int
	strcmp();
extern int
	strlen();

extern gptr
	*memcpy(),
	*memset();

#endif

#endif
