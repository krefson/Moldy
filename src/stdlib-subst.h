#ifndef _MOLDY_STDLIB_H_
#define _MOLDY_STDLIB_H_

#ifdef ANSI_LIBS
#   include <stdlib.h>
#else

#  include "stddef.h"

   extern int atoi ();
   extern char * calloc ();
   extern void free ();
   extern char * malloc ();
   extern char * realloc ();
   extern void exit ();
   extern int abs ();
#endif

#endif /* __STDLIB_H_ */





