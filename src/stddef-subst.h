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
