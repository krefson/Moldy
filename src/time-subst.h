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
#ifndef _H_TIME
#ifndef __TIME_H__
#ifndef _TIME_INCLUDED

#include <time.h>
#ifndef ANSI_LIBS
#   if defined(unix) || defined(__unix__)
/*
 *  We must protect the inclusion of <sys/types.h>.  
 */
#      ifndef SYS_TYPES_INCLUDED
#         define SYS_TYPES_INCLUDED
#         include <sys/types.h>
#      endif
       extern time_t	time();
       extern char * ctime();
#   endif
#endif
#define __time_h
#define _TIME_H
#define _TIME_H_
#define __TIME_H
#define _H_TIME
#define __TIME_H__
#define _TIME_INCLUDED

#endif
#endif
#endif
#endif
#endif
#endif
#endif
