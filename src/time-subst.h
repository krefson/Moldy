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
/*
 *  Time.h   replacement for ANSI one
 */
#ifndef __time_h
#ifndef _TIME_H
#ifndef _TIME_H_
#ifndef __TIME_H
#ifndef __TIME_H__
#ifndef _TIME_INCLUDED

#include <time.h>

#ifndef ANSI_LIBS
#   if defined(unix) || defined(__unix__)
/*
 *  We must protect the inclusion of <sys/types.h>.  But size_t may already
 *  be defined in "stddef.h" , so we define it out of the way.  The twist
 *  comes for a GNU CC compilation.  That *may* have fixed includes, but
 *  we can't rely on that, so define size_t to the same symbol as the GNU
 *  header file does to stop it complaining.   Roll on ANSI.
 */
#      ifndef SYS_TYPES_INCLUDED
#         define SYS_TYPES_INCLUDED
#         define size_t ___size_t
#         include <sys/types.h>
#         undef  size_t
#      endif
       extern time_t	time();
#   else
#      ifdef CRAY
          typedef long	time_t;
          extern time_t	time();
#      endif
#   endif
#endif
#define __time_h

#endif
#endif
#endif
#endif
#endif
#endif
