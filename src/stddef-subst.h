#ifndef MOLDY_STDDEF_H
#define MOLDY_STDDEF_H

#ifdef ANSI_LIBS
#   include <stddef.h>
#else
#   ifdef unix
#      ifndef _SIZE_T
#         define  _SIZE_T
          typedef unsigned size_t;
#      endif
#   else
       typedef unsigned long size_t;
#   endif
#endif

#endif /* _STDDEF_H */
