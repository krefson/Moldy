/*
 * $Header: string.h,v 1.1 89/04/25 12:32:16 keith Exp $
 *
 * $Log:	string.h,v $
 * 
 */
#ifndef STRING_ALREADY

#ifdef unix

#ifdef USG
#include	<string.h>
#include	<memory.h>
#else  /* BSD */
#define strchr(c,p) index(c,p)
#include	<strings.h>
char		*memcpy();
#endif

#else	/* Not unix */

#ifndef CRAY
#include	<string.h>
#endif

#endif /* unix      */

#endif
