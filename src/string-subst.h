/*
 * $Header: /home/eeyore/keith/md/moldy/RCS/string.h,v 1.5 90/03/27 17:38:29 keith Exp $
 *
 * $Log:	string.h,v $
 * 
 */
#ifndef STRING_ALREADY
#define STRING_ALREADY

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
