/*
 *  Time.h   replacement for ANSI one
 */
#ifndef MOLDY_TIME_H_INCLUDED
#define MOLDY_TIME_H_INCLUDED
#ifndef __STDC__
#   ifdef unix
#      ifndef SYS_TYPES_INCLUDED
#         define SYS_TYPES_INCLUDED
#         include <sys/types.h>
#      endif
       extern time_t	time();
#   else
#      ifdef CRAY
          typedef long	time_t;
          extern time_t	time();
#      endif
#   endif
#endif
#include <time.h>
#endif
