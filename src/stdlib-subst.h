#ifndef _MOLDY_STDLIB_H_
#define _MOLDY_STDLIB_H_

#ifdef __STDC__
#   include <stdlib.h>
#else
   extern int atoi ();
   extern char * calloc ();
   extern void free ();
   extern char * malloc ();
   extern char * realloc ();
   extern void exit ();
   extern int abs ();
#endif

#endif /* __STDLIB_H_ */





