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
#ifndef MOLDY_STDDEF_H
#define MOLDY_STDDEF_H

#ifdef ANSI_LIBS
#   include <stddef.h>
#else
#   ifndef __SIZE_T
#      ifndef __SIZE_T__
#         ifndef _SIZE_T
#            ifndef _SIZE_T_
#		ifndef _GCC_SIZE_T
#                  ifdef sun
#                     include <sys/stdtypes.h>
#                  else
#		      define  _GCC_SIZE_T
#                     define  __SIZE_T
#                     define  __SIZE_T__
#                     define  _SIZE_T
#                     define  _SIZE_T_
#                     if defined(unix) && (!defined(__GNUC__))
                         typedef unsigned size_t;
#                     else
                         typedef unsigned long size_t;
#                     endif
#                  endif
#		endif
#            endif
#         endif
#      endif
#   endif
#   ifndef NULL
#      define NULL 0
#   endif
#endif

#endif /* _STDDEF_H */
