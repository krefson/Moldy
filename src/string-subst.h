/*
 * $Header: string.h,v 1.1 89/04/25 17:38:05 keith Exp $
 *
 * $Log:	string.h,v $
 * Revision 1.1  89/04/25  17:38:05  keith
 * Initial revision
 * 
 * 
 */
#ifndef STRING_ALREADY

#ifdef unix

#if USG || sysV || SysV
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
