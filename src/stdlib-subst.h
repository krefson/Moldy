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
#ifndef _STDLIB_H_
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
   extern long strtol();
#endif

#endif /* _MOLDY_STDLIB_H_ */
#endif /* __STDLIB_H_ */





