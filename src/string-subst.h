/*
 * $Header: /home/eeyore/keith/md/moldy/RCS/string.h,v 1.6 90/04/14 16:02:06 keith Exp $
 *
 * $Log:	string.h,v $
 * Revision 1.6  90/04/14  16:02:06  keith
 * Replaced calls to system headers with own definitions for portability.
 * 
 * 
 */
#ifndef STRING_ALREADY
#define STRING_ALREADY

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

extern char
	*memcpy(),
	*memset();

#endif
