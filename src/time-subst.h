/*
 *  Time.h   replacement for ANSI one
 */
#ifndef MOLDY_TIME_H_INCLUDED
#define MOLDY_TIME_H_INCLUDED
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
#include <time.h>
#endif
