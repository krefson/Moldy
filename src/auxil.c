/******************************************************************************
 * Aux		This file contains auxilliary service functions and (almost)  *
 *		all machine-dependant stuff.  The ANSI standard interface is  *
 *              used as much as possible for library functions; if a machine  *
 *		is missing something an ANSI interface routine is supplied.   *
 ******************************************************************************
 *      Revision Log
 *       $Log:	aux.c,v $
 * Revision 1.3  89/06/01  21:23:04  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.2  89/06/01  18:00:42  keith
 * Moved `vadd()' from aux.c to force.c for ease of vectorisation.
 * Now no need to compile aux.c with vectorisation.
 * 
 * Revision 1.1  89/04/27  15:01:14  keith
 * Initial revision
 * 
 * 
 */
#ifndef lint
static char *RCSid = "$Header: aux.c,v 1.7 89/06/09 15:45:11 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#include	<stdio.h>
#include	<ctype.h>
#include	<math.h>
#include 	"string.h"
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library declarations ============================*/
#ifdef unix
time_t	time();
#endif
/*============================================================================*/
/******************************************************************************
 *  Vector routines - special cases for cray and convex to call vectorised    *
 *  library versions.  Other machines do it the hard way.		      *
 *  CRAY versions.							      *
 ******************************************************************************/
#if defined(CRAY)
double vdot(n,x,ix,y,iy)
int	n;
real	*x, *y;
int	ix, iy;
{
   double SDOT();
   return(SDOT(&n, x, &ix, y, &iy));
}
double sum(n,x,ix)
int	n;
real	*x;
int	ix;
{
   double SSUM();
   return(SSUM(&n, x, &ix));
}
void	vscale(n,s,x,ix)
int	n;
double	s;
real	*x;
int	ix;
{
   void SSCAL();
   SSCAL(&n, &s, x, &ix);
}
void saxpy(n, sa, sx, ix, sy, iy)
int	n, ix, iy;
double	sa;
double	sx[], sy[];
{
   void SAXPY();
   SAXPY(&n, &sa, sx, &ix, sy, &iy);
}
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   int	ISRCHFLT();
   return( ISRCHFLT(&n, x, &ix, &s) - 1);
}
void	gather(n, a, b, ix)
int	n;
real	a[], b[];
int	ix[];
{
   void GATHER();
   GATHER(&n, a, b+1, ix);
}
void	scatter(n, a, ix, b)
int	n;
real	a[], b[];
int	ix[];
{
   void SCATTER();
   SCATTER(&n, a+1, ix, b);
}

#endif
/******************************************************************************
 *  Vector routines - CONVEX versions.  The calls are to VECLIB functions.    *
 ******************************************************************************/
#if defined(convexvc)
#define SINGLE (sizeof(real) == sizeof(float))
double vdot(n,x,ix,y,iy)
int	n;
real	*x, *y;
int	ix, iy;
{
   double ddot_(); float sdot_();
   if(SINGLE)
     return(sdot_(&n, x, &ix, y, &iy));
   else
     return(ddot_(&n, x, &ix, y, &iy));
}
double sum(n,x,ix)
int	n;
real	*x;
int	ix;
{
   double dsum_(); float ssum_();
   if(SINGLE)
      return(ssum_(&n, x, &ix));
   else
      return(dsum_(&n, x, &ix));
}
void	vscale(n,s,x,ix)
int	n;
double	s;
real	*x;
int	ix;
{
   float s2 = s;
   void dscal_(), sscal_();
   if(SINGLE)
      sscal_(&n, &s2, x, &ix);
   else
      dscal_(&n, &s, x, &ix);
}
void saxpy(n, sa, sx, ix, sy, iy)
int	n, ix, iy;
double	sa;
double	sx[], sy[];
{
   float s2 = sa;
   void daxpy_(), saxpy_();
   if(SINGLE)
      saxpy_(&n, &s2, sx, &ix, sy, &iy);
   else
      daxpy_(&n, &sa, sx, &ix, sy, &iy);
}
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   float s2 = s;
   int	i, idsvlt_(), issvlt_();
   if(SINGLE)
      return( (i = issvlt_(&n, x, &ix, &s2)) ? i - 1 : n);
   else
      return( (i = idsvlt_(&n, x, &ix, &s)) ? i - 1 : n);
}
void	gather(n, a, b, ix)
int	n;
real	a[], b[];
int	ix[];
{
   void dgthr_(), sgthr_();
   if(SINGLE)
      sgthr_(&n, b+1, a, ix);
   else
      dgthr_(&n, b+1, a, ix);
}
void	scatter(n, a, ix, b)
int	n;
real	a[], b[];
int	ix[];
{
   void dsctr_(), ssctr_();
   if(SINGLE)
      ssctr_(&n, b, ix, a+1);
   else
      dsctr_(&n, b, ix, a+1);
}
#endif
/******************************************************************************
 * alliant versions - with cray librarries				      *
 ******************************************************************************/
#if defined(alliant)
double vdot(n,x,ix,y,iy)
int	n;
real	*x, *y;
int	ix, iy;
{
   double sdot_();
   return(sdot_(&n, x, &ix, y, &iy));
}
double sum(n,x,ix)
int	n;
real	*x;
int	ix;
{
   double ssum_();
   return(ssum_(&n, x, &ix));
}
void	vscale(n,s,x,ix)
int	n;
double	s;
real	*x;
int	ix;
{
   void sscal_();
   sscal_(&n, &s, x, &ix);
}
void saxpy(n, sa, sx, ix, sy, iy)
int	n, ix, iy;
double	sa;
double	sx[], sy[];
{
   void saxpy_();
   saxpy_(&n, &sa, sx, &ix, sy, &iy);
}
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   int	isrchflt_();
   return( isrchflt_(&n, x, &ix, &s) - 1);
}
void	gather(n, a, b, ix)
int	n;
real	a[], b[];
int	ix[];
{
   void gather_();
   gather_(&n, a, b+1, ix);
}
void	scatter(n, a, ix, b)
int	n;
real	a[], b[];
int	ix[];
{
   void scatter_();
   scatter_(&n, a+1, ix, b);
}

#endif
/******************************************************************************
 *  Vector handling functions - versions for scalar machines.		      *
 ******************************************************************************/
#if ! defined(CRAY) && ! defined(convexvc) && ! defined(alliant)
double vdot(n,x,ix,y,iy)
int	n;
real	x[], y[];
int	ix, iy;
{
   register double	dot=0.0;
   int i, j;
   for(i = (n-1)*ix, j = (n-1)*iy; i >= 0; i -= ix, j -= iy)
      dot += x[i] * y[j];
   return(dot);
}
double sum(n,x,ix)
register int	n;
register real	x[];
register int	ix;
{
   register double	s=0.0;
   int i;
   for(i = (n-1)*ix; i >= 0; i -= ix)
      s += x[i];
   return(s);
}
void	vscale(n,s,x,ix)
register int	n;
register double	s;
register real	x[];
register int	ix;
{
   int i;
   for(i = (n-1)*ix; i >= 0; i -= ix)
      x[i] *= s;
}
void saxpy(n, sa, sx, ix, sy, iy)
int	n, ix, iy;
double	sa;
double	sx[], sy[];
{
   int i, j;
   for(i = (n-1)*ix, j = (n-1)*iy; i >= 0; i -= ix, j -= iy)
      sy[j] += sa * sx[i];
}
int	search_lt(n, x, ix, s)
int	n;
real	x[];
int	ix;
double	s;
{
   int i;
   for( i = 0; i < n*ix; i += ix )
      if( x[i] < s )
         break;
   return(i);
}
void	gather(n, a, b, ix)
int	n;
real	a[], b[];
int	ix[];
{
   int i;
   for(i = n-1; i >= 0; i--)
      a[i] = b[ix[i]];
}
void	scatter(n, a, ix, b)
int	n;
real	a[], b[];
int	ix[];
{
   int i;
   for(i = n-1; i >= 0; i--)
      a[ix[i]] = b[i];
}
#endif
/******************************************************************************
 *  random number generators. Note they assume 'unsigned' of at least 32 bits *
 ******************************************************************************/
#define  im  1771875L
#define  ia  2416L
#define  ic  374441L
static unsigned long jran = 1;
void	smdrand(seed)
unsigned long seed;
{
   jran = seed;
}
double	mdrand()
{
   jran = (jran*ia + ic) % im;
   return (double)jran/(double)im;
}
/******************************************************************************
 *  Precision. Return smallest eps s.t. 1.0+eps > 1			      *
 ******************************************************************************/
double	precision()
{
   static	int	first=1;
   double		eps = 0.5, eps2;
   double  		*eps1 = &eps2;		/* Use pointer to force store */
   
   if(first)
   {
      first = 0;
      do
      {
         eps = eps * 0.5;
         *eps1 = 1.0 + eps;
      }   while(*eps1 > 1.0);
      eps *= 2.0;
   }
   return(eps);
}
/******************************************************************************
 *  cctime()  Return ctime() but without the newline			      *
 ******************************************************************************/
char	*cctime(timeloc)
time_t	*timeloc;
{
   char	*p = ctime(timeloc);
   *strchr(p, '\n') = '\0';
   return(p);
}
/******************************************************************************
 *  atime.  Return a pointer to the time and date string		      *
 ******************************************************************************/
char	*atime()
{
   time_t t = time((time_t*)0);
   return(cctime(&t));
}
/******************************************************************************
 *  cpu.  Return (double) the current cpu time in seconds.		      *
 ******************************************************************************/
#ifdef unix
#ifdef USG

#define TICK 60.0
#define size_t NOTHING
#define time_t NOTHING
#include <sys/types.h>
#include <sys/times.h>
#undef size_t
#undef time_t

double        cpu()
{
   struct tms buf;
 
   (void)times(&buf);

   return((buf.tms_utime + buf.tms_stime)/TICK);
}
#else

#define KERNEL          /* This stops sys/time.h from including time.h again  */
#define time_t NOTHING
#include <sys/time.h>
#include <sys/resource.h>
#undef time_t

double	cpu()	/* The standard unix 'clock' wraps after 36 mins.	      */
{
   struct rusage ru;
 
   (void)getrusage(RUSAGE_SELF, &ru);

   return(ru.ru_utime.tv_sec  + ru.ru_stime.tv_sec
	  + 1.0e-6 * (ru.ru_utime.tv_usec + ru.ru_stime.tv_usec));
}
#endif
#else
#ifdef CRAY
#define CLK_TCK 1000000
#endif
double	cpu()
{
   return((double)clock() / CLK_TCK);
}
#endif
/******************************************************************************
 * Zero real.  assumes memset or bzero.					      *
 ******************************************************************************/
#if defined(unix) && !defined(USG)   	/* BSD Unix - use bzero		      */
void	zero_real(r,n)
real	*r;
int	n;
{
   int bzero();
   (void)bzero((char *)r, n * sizeof(real));
}
#else					/* SYSV unix or ANSI - use memset()   */
void	zero_real(r,n)
real	*r;
int	n;
{
   (void)memset((char *)r, 0, n * sizeof(real));
}
#endif
/******************************************************************************
 * memcpy   replacement for BSD machines which don't have it		      *
 ******************************************************************************/
#if defined(unix) && !defined(USG)
char *memcpy(s1, s2, n)
char *s1, *s2;
int n;
{
   int bcopy();
   (void)bcopy((char *)s2, (char *)s1, n);
   return(s1);
}
#endif
/******************************************************************************
 * remove.  delete (unlink) a file.   (ANSI replacement)		      *
 ******************************************************************************/
#if defined(unix) || defined(CRAY)
remove(file)
char	*file;
{
   return (unlink(file));
}
#endif
/******************************************************************************
 *  replace - replace file1 by file2, renaming or overwriting		      *
 ******************************************************************************/
#ifdef vms
int	replace(file1, file2)
char	*file1, *file2;
{
   int f;  char name[L_name];
   f = rename(file1, file2);
   if(f == 0) (void)remove(strcat(strcpy(name,file2),";-1"));
   return(f);
}
#else

#if defined(CRAY) && !defined(unix)	/* COS only - not UNICOS	      */

void cos_parse(name, own, id, dn)
char	*name, *own, *id, *dn;
{
   char	*tok[3], *t = NULL, *name2 = strdup(name);
   tok[0] = tok[1] = tok[2] = "";

   while ( (t = strtok(t ? 0 : name2, "/")) != NULL )
   {
      tok[0] = tok[1]; tok[1] = tok[2]; tok[2] = t;
   }
   (void)strcpy(dn, tok[2]);
   (void)strcpy(id, tok[1]);
   (void)strcpy(own, tok[0]);
   cfree(name2);
}

#define PDNLEN 16

int     replace(file1, file2)   /* Actually saves a DN as a PDN and purges    */
char    *file1, *file2;
{
   void SAVE(), RELEASE(), ACCESS(), DELETE();
   char from[PDNLEN], to[PDNLEN], tmp[PDNLEN], id[PDNLEN], own[PDNLEN];
   int flag; int i, oldflag = -1;
   /*
    * The CRAY JCL routines expect nulls up to the end of each word (8 bytes)
    */
   for( i = 0; i < PDNLEN; i++ )
      from[i] = to[i] = tmp[i] = id[i] = '\0';
   strncpy(from, file1, PDNLEN-1); 
   (void)tmpnam(tmp);

   cos_parse(file2, own, id, to);

   ACCESS(&oldflag,"DN",tmp,"PDN",to,"ID", id, "UQ");
   SAVE(&flag,"DN",from,"PDN",to, "ID", id);
   if( flag )     return(flag);
   RELEASE(&flag,"DN",from);            
   if( flag )   return(flag);
   if( ! oldflag )
      DELETE(&flag,"DN",tmp); 
   return(flag);
}
#else

int	replace(file1, file2)
char	*file1, *file2;
{
   int f;
   f = rename(file1, file2);
   if(f != 0 && remove(file2) )		/* If already exists try remove target*/
      f = rename(file1, file2);
   return(f);
}
#endif
#endif
/******************************************************************************
 * err_fn  error function. See Abramowitz & Stegun p299.		      *
 ******************************************************************************/
#define E1	 0.254829592
#define E2	-0.284496736
#define E3	 1.421413741
#define E4	-1.453152027
#define E5	 1.061405429

#define PP	 0.3275911
#define POLY5(t)   ((t)*(E1 + (t)*(E2 + (t)*(E3 + (t)*(E4 + (t)*E5)))))

double	err_fn(x)
double	x;
{
   double t = 1.0/(1.0 + PP * x);
   if ( x >= 0.0 )
      return( 1.0 - POLY5(t) * exp(-x*x) );
   else
      return( -err_fn(-x) );
}
/******************************************************************************
 *  vprintf for those machines which don't have it			      *
 ******************************************************************************/
#if ! defined(sun) && ! defined(vms) && ! defined(ANSI) && ! defined(__STDC__)

#include <varargs.h>

#ifdef CRAY

int	vprintf (fmt, args)
char   	*fmt;
va_list args;
{
    return(XPRINTF(&fmt, &args, &stdout));
}

#else
#if defined(convex) || defined(sequent)

int	vprintf (fmt, args)
char	*fmt;
va_list	args;
{
   return(_doprnt(fmt, args, stdout));
}


#else		/* No _doprnt or equivalent */

#include <ctype.h>
int	vprintf (format, ap)
char	*format;
va_list	ap;
{
    int     pos, charsout, fpos, error;
    char    fmt[1024], temps[1024];
    char    ch;

    pos = charsout = error = 0;
    while (format[pos] != 0)
    {
        if (format[pos] != '%')
        {
            putchar (format[pos++]);
        }
        else
        {
            fmt[0] = '%';
            pos++;
            fpos = 1;
            while (format[pos] != 0 && !islower (format[pos]) &&
                format[pos] != '%')
                fmt[fpos++] = format[pos++];
            if (format[pos] == 0)
            {
            /* error in format string */
                error++;
            }
            ch = fmt[fpos++] = format[pos++];
            fmt[fpos] = 0;
#ifdef DEBUG
            printf ("Format is \"%s\"\n", fmt);
#endif
        /* Now fmt contains the format for this part of the output and ch
           contains the conversion code. */
            temps[0] = 0;
            switch (ch)
            {
              case 'n':
                *va_arg (ap, int *) = charsout;
                break;
              case 's':
                sprintf (temps, fmt, va_arg (ap, char *));
                break;
              case 'p':
                sprintf (temps, fmt, va_arg (ap, void *));
                break;
              case 'c':
              case 'd':
              case 'i':
              case 'o':
              case 'u':
              case 'x':
              case 'X':
                sprintf (temps, fmt, va_arg (ap, int));
                break;
              case 'e':
              case 'E':
              case 'f':
              case 'g':
              case 'G':
                sprintf (temps, fmt, va_arg (ap, double));
                break;
              default:
                (void)strcpy (temps, fmt);
                error++;
            }
#ifdef DEBUG
            printf ("temps is \"%s\"\n", temps);
#endif
            puts (temps);
            charsout += strlen (temps);
#ifdef DEBUG
            printf ("still to interpret \"%s\"\n", (char *) &format[pos]);
#endif
        }
    }
    if (error)
        return (-1);
    return (charsout);
}
#endif		/* convex or sequent */
#endif		/* Cray		     */
#endif		/* Vprintf needed    */

