/*
 * $Header: /home/tigger/keith/md/RCS/string.h,v 1.3 89/11/20 12:05:01 keith Exp $
 *
 * $Log:	string.h,v $
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

#ifdef unix

#if USG || sysV || SysV || cray
#include	<string.h>
#include	<memory.h>
#else  /* BSD */
#define strchr(c,p) index(c,p)
#include	<strings.h>
char		*memcpy();
char		*strdup();
#endif

#else	/* Not unix */

#ifndef CRAY
#include	<string.h>
#endif

#endif /* unix      */

#endif
