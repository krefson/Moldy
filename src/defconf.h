/******************************************************************************
 *  Non-unix - specific Configurational information.  	                      *
 ******************************************************************************/
#if defined(__vms) && !defined(vms)
#   define vms
#endif
#ifdef vms
#   define ALLOC_ALIGN
#   define HAVE_POPEN
#   define HAVE_GETOPT
#   define STDC_HEADERS
#   define BACKUP_FILE		"MDBACKUP.DAT"
#   define TEMP_FILE		"MDTEMPXXXX.DAT"
#   define LOCKEX		"$LCK"
#endif

#if defined(__MSDOS__)
#   define ALLOC_SEPARATELY
#   define STDC_HEADERS
#   define BACKUP_FILE		"MDBCK"
#   define TEMP_FILE		"MDTEMPX"
#   define LOCKEX		"$LK"
#endif

#if defined(_WIN32) || defined(__WIN32) || defined(__WIN32__)
#   ifndef __WIN32__
#      define __WIN32__
#   endif
#   define STDC_HEADERS
#   define LOCKEX		"_lck"
#   define ALLOC_ALIGN
#endif

#ifdef CMS
#   define BACKUP_FILE		"MDBACKUP MOLDY A1"
#   define TEMP_FILE		"MDTEMP XXXXXXXX A1"
#endif

