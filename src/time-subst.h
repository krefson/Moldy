/*
 *  Time.h   replacement for ANSI one
 */
#ifndef MOLDY_TIME_H_INCLUDED
#define MOLDY_TIME_H_INCLUDED
#ifndef ANSI_LIBS
#   ifdef unix
#      ifdef __GNUC__			/* Gcc already protects size_t      */
#         ifndef SYS_TYPES_INCLUDED
#            define SYS_TYPES_INCLUDED
#            include <sys/types.h>
#         endif
#      else				/* Protect size_t by hand.	    */
#         ifndef SYS_TYPES_INCLUDED
#            define SYS_TYPES_INCLUDED
#            define size_t XXSIZE_T
#            include <sys/types.h>
#            undef  size_t
#         endif
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
