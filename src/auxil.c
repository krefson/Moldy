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
/******************************************************************************
 * Aux		This file contains auxilliary service functions and (almost)  *
 *		all machine-dependant stuff.  The ANSI standard interface is  *
 *              used as much as possible for library functions; if a machine  *
 *		is missing something an ANSI interface routine is supplied.   *
 ******************************************************************************
 */
/*========================== program include files ===========================*/
#include	"defs.h"
/*========================== Library include files ===========================*/
#include	<math.h>
#include 	<string.h>
#include	<stdio.h>
#include        <ctype.h>
#include	<stdlib.h>
/*================= System Library include files - unix only ================*/
#if defined(HAVE_TIMES) && defined(HAVE_SYS_TIMES_H)
#   include <sys/times.h>
#else
#   if defined HAVE_GETRUSAGE
#      include <sys/resource.h>
#   endif
#endif

#if defined(HAVE_GETTIMEOFDAY) && \
     (defined(TIMES_RETURNS_STATUS) || \
      !defined(HAVE_TIMES) || !defined(HAVE_SYS_TIMES_H))
#   ifdef HAVE_SYS_TIME_H
#      include <sys/time.h>
#      ifndef TIME_WITH_SYS_TIME
#         define EXCLUDE_TIME_H
#      endif 
#   endif
#endif

#ifndef EXCLUDE_TIME_H
#   include <time.h>
#endif

#ifndef CLK_TCK
#   ifdef HAVE_SYS_PARAM_H
#      include <sys/param.h>
#   endif
#   ifdef HZ
#      define CLK_TCK HZ
#   else
#      define CLK_TCK 60	/* Really old unices defined 60. May be wrong.*/
#   endif
#endif
#if !defined(CLOCKS_PER_SEC) && defined(CLK_TCK)
#   define  CLOCKS_PER_SEC CLK_TCK
#endif
/*========================== External function declarations ==================*/
gptr            *talloc(int n, size_mt size, int line, char *file);	       /* Interface to memory allocator       */
void            tfree(gptr *p);	       /* Free allocated memory	      	      */
/*============================================================================*/
/******************************************************************************
 *  Vector routines - special cases for cray and convex to call vectorised    *
 *  library versions.  Other machines do it the hard way.		      *
 *  CRAY versions.							      *
 ******************************************************************************/
#if defined(_CRAY1)
#define VECTS_DEFINED

double vdot(int n, real *x, int ix, real *y, int iy)
{
   double SDOT();
   return(SDOT(&n, x, &ix, y, &iy));
}
double sum(int n, real *x, int ix)
{
   double SSUM();
   return(SSUM(&n, x, &ix));
}
void	vscale(int n, double s, real *x, int ix)
{
   void SSCAL();
   SSCAL(&n, &s, x, &ix);
}
int	search_lt(int n, real *x, int ix, double s)
{
   int	ISRCHFLT();
   if( n <= 0 )
      return(0);
   return( ISRCHFLT(&n, x, &ix, &s) - 1);
}
void	gather(int n, real *a, real *b, int *ix, int lim)
{
   void GATHER();
   GATHER(&n, a, b+1, ix);
}
#endif
/******************************************************************************
 *  Vector routines - CONVEX versions.  The calls are to VECLIB functions.    *
 ******************************************************************************/
#if defined(__convexc__)
#define VECTS_DEFINED

#define SINGLE (sizeof(real) == sizeof(float))
double vdot(int n, real *x, int ix, real *y, int iy)
{
   double ddot_(); float sdot_();
   if(SINGLE)
     return(sdot_(&n, x, &ix, y, &iy));
   else
     return(ddot_(&n, x, &ix, y, &iy));
}
double sum(int n, real *x, int ix)
{
   double dsum_(); float ssum_();
   if(SINGLE)
      return(ssum_(&n, x, &ix));
   else
      return(dsum_(&n, x, &ix));
}
void	vscale(int n, double s, real *x, int ix)
{
   float s2 = s;
   void dscal_(), sscal_();
   if(SINGLE)
      sscal_(&n, &s2, x, &ix);
   else
      dscal_(&n, &s, x, &ix);
}
int	search_lt(int n, real *x, int ix, double s)
{
   float s2 = s;
   int	i, idsvlt_(), issvlt_();
   if(SINGLE)
      return( (i = issvlt_(&n, x, &ix, &s2)) ? i - 1 : n);
   else
      return( (i = idsvlt_(&n, x, &ix, &s)) ? i - 1 : n);
}
                            /*ARGSUSED*/
void	gather(int n, real *a, real *b, int *ix, int lim)
{
   void dgthr_(), sgthr_();
   if(SINGLE)
      sgthr_(&n, b+1, a, ix);
   else
      dgthr_(&n, b+1, a, ix);
}
#endif
/******************************************************************************
 * alliant versions - with cray librarries				      *
 ******************************************************************************/
#if defined(alliant)
#define VECTS_DEFINED

double vdot(int n, real *x, int ix, real *y, int iy)
{
   double sdot_();
   return(sdot_(&n, x, &ix, y, &iy));
}
double sum(int n, real *x, int ix)
{
   double ssum_();
   return(ssum_(&n, x, &ix));
}
void	vscale(int n, double s, real *x, int ix)
{
   void sscal_();
   sscal_(&n, &s, x, &ix);
}
int	search_lt(int n, real *x, int ix, double s)
{
   int	isrchflt_();
   return( isrchflt_(&n, x, &ix, &s) - 1);
}
                            /*ARGSUSED*/
void	gather(int n, real *a, real *b, int *ix, int lim)
{
   void gather_();
   gather_(&n, a, b+1, ix);
}
#endif
/******************************************************************************
 *  Vector handling functions - C versions. 		     		      *
 ******************************************************************************/
#ifndef VECTS_DEFINED
double vdot(int n, real *x, int ix, real *y, int iy)
{
   register double	dot=0.0;
   register int i, j;
   if( ix == iy && ix == 1)
   {
VECTORIZE
      for(i = 0; i < n; i++)
         dot += x[i] * y[i];
   }
   else
   {
VECTORIZE
      for(i = j = 0; i != n*ix; i += ix)
      {
         dot += x[i] * y[j];
         j += iy;
      }
   }
   return(dot);
}
double sum(register int n, register real *x, register int ix)
{
   register double	sa=0.0,sb=0.0,sc=0.0,sd=0.0,s1=0.0,s2=0.0,s3=0.0,s4=0.0;
   int i;

   if( ix == 1 )
   {

     for( i = 0; i < n-3; i+=4)
     {
	s1 += sa;
	sa = x[i];
	s2 += sb;
	sb = x[i+1];
	s3 += sc;
	sc = x[i+2];
	s4 += sd;
	sd = x[i+3];
	
     }
     s1 += sa;
     s2 += sb;
     s3 += sc;
     s4 += sd;
     for( ; i < n; i++)
       s1 += x[i];
   } 
   else
   {

      for(i = 0; i != n*ix; i += ix)
	s1 += x[i];
   }
   return(s1+s2+s3+s4);
}

void	vscale(register int n, register double s, register real *x, register int ix)
{
   int i;
VECTORIZE
   for(i = 0; i != n*ix; i += ix)
      x[i] *= s;
}

int	search_lt(int n, real *x, int ix, double s)
{
   int i,j;

   if( n < 0 )
      n = 0;

   if ( ix == 1 ) 
   {
      for( i = 0; i < n; i ++ )
	 if( x[i] < s )
	    return i;
   }
   else
   {
VECTORIZE
      for( i = 0, j=0; i < n*ix; i += ix, j++ )
         if( x[i] < s )
	    return j;
   }
   return n;
}
                            /*ARGSUSED*/
void	gather(int n, real *a, real *b, int *ix, int lim)
{
   int i;
VECTORIZE
#ifdef _CRAY_T3E
#pragma _CRI cache_bypass  b
#endif
   for( i = 0; i < n; i++)
   {
#ifdef DEBUG
      if(ix[i] < 0 || ix[i] >= lim )
	 message(0,0,2, "Gather index out of bounds %d %d\n", i, ix[i]);
#endif
      a[i] = b[ix[i]];
   }
}

#endif
/******************************************************************************
 *  random number generators. Note they assume 'unsigned' of at least 32 bits *
 ******************************************************************************/
#define  im  1771875L
#define  ia  2416L
#define  ic  374441L
static unsigned long jran = 1, jran1 = 1;
void	smdrand(long unsigned int seed)
{
   jran = seed;
}
unsigned long getseed(void)
{
   return jran;
}
double	mdrand(void)
{
   jran = (jran*ia + ic) % im;
   return (double)jran/(double)im;
}
void	smdrand1(long unsigned int seed)
{
   jran = seed;
}
double	mdrand1(void)
{
   jran1 = (jran1*ia + ic) % im;
   return (double)jran1/(double)im;
}
/******************************************************************************
 *  Precision. Return smallest eps s.t. 1.0+eps > 1			      *
 ******************************************************************************/
double	precision(void)
{
   static	int	first=1;
   static	double	eps = 0.5;
   double volatile	eps2, *eps1 = &eps2;	/* Use pointer to force store */
   double volatile      junk, *ptr = &junk;
   
   if(first)
   {
      first = 0;
      do
      {
         eps = eps * 0.5;
         *eps1 = 1.0 + eps;
	 *ptr  = 1.0;		/* This prevents optimization of *eps1 */
      }   while(*eps1 > 1.0);
      eps *= 2.0;
   }
   return(eps);
}
/******************************************************************************
 *  cpu.  Return (double) the current cpu time in seconds.		      *
 ******************************************************************************/
#if defined(HAVE_TIMES) && defined(HAVE_SYS_TIMES_H)
double        cpu(void)
{
   struct tms buf;
 
   (void)times(&buf);

   return (buf.tms_utime + buf.tms_stime)/(double)CLK_TCK;
}
#else
#   if defined HAVE_GETRUSAGE
double  cpu(void)   /* The standard unix 'clock' wraps after 36 mins.            */
{
   struct rusage ru;
   int getrusage(int, struct rusage *);
 
   (void)getrusage(RUSAGE_SELF, &ru);

   return(ru.ru_utime.tv_sec  + ru.ru_stime.tv_sec
          + 1.0e-6 * (ru.ru_utime.tv_usec + ru.ru_stime.tv_usec));
}
#   else
double	cpu(void)
{
   return((double)clock() / CLOCKS_PER_SEC);
}
#   endif
#endif
/******************************************************************************
 *  rt_clock(). Return elapsed time in s.				      *
 ******************************************************************************/
#if defined(HAVE_TIMES) && defined(HAVE_SYS_TIMES_H) && !defined(TIMES_RETURNS_STATUS)
double rt_clock(void)
{
   struct tms buf;
   static int use_time=0;
   clock_t t;
 
   if( ! use_time )
   {
      t=times(&buf);
      if( t <= 0 )   /*Times failed.  Old BSD system? */
	 use_time = 1;
      else
	 return t/(double)CLK_TCK;
   }
   return time((time_t *)0);
}
#else
#   if defined(HAVE_GETTIMEOFDAY)

double rt_clock(void)
{
   int gettimeofday(struct timeval *, void *);
   struct timeval tp;
   (void)gettimeofday(&tp,(struct timezone *)0);
   return(tp.tv_sec + tp.tv_usec*0.000001);
}

#   else 

double rt_clock(void)
{
   return (double)time((time_t *)0);
}
#   endif
#endif

/******************************************************************************
 *  cctime()  Return ctime() but without the newline			      *
 ******************************************************************************/
char	*cctime(time_mt *timeloc)
{
   time_t  loc = (time_t)*timeloc;
   char	*p = ctime(&loc);
   *strchr(p, '\n') = '\0';
   return(p);
}
/******************************************************************************
 *  atime.  Return a pointer to the time and date string		      *
 ******************************************************************************/
char	*atime(void)
{
   time_mt t = (time_mt)time((time_t*)0);
   return(cctime(&t));
}
/******************************************************************************
 * Zero real.  Set array to floating-point zero.			      *
 ******************************************************************************/
void	zero_real(real *r, int n)
{
   int i;
VECTORIZE
   for(i = 0; i < n; i++)
      r[i] = 0.0;
}
void	zero_double(double *r, int n)
{
   int i;
VECTORIZE
   for(i = 0; i < n; i++)
      r[i] = 0.0;
}
void	zero_dbls(double *r, size_mt n)
{
   long i;
VECTORIZE
   for(i = 0; i < n; i++)
      r[i] = 0.0;
}
/******************************************************************************
 *  replace - replace file1 by file2, renaming or overwriting		      *
 ******************************************************************************/
#ifdef UNRESTRICTED_FILE_NAMES
int	replace(char *file1, char *file2)
{
   int f;
   char *backup = aalloc(strlen(file2)+2, char);
   (void)strcat(strcpy(backup, file2), "%");

   (void)rename(file2, backup);		/* Backup old version - ignore failure*/
   f = rename(file1, file2);		/* Rename old version to new	      */
   xfree(backup);
   return(f);
}
#   else
int	replace(char *file1, char *file2)
{
   int f;
   int rename(const char *, const char *);
   f = rename(file1, file2);
   if(f != 0 && remove(file2) == 0)	/* If already exists try remove target*/
      f = rename(file1, file2);
   return(f);
}
#endif
/******************************************************************************
 *  Purge Remove 1 old version of a file.  (Do nothing if O/S has no versions)*
 ******************************************************************************/
#ifdef vms 
void	purge(file)
char	*file;
{
   char *name = aalloc(strlen(file)+4,char);
   strcpy(name, file);
   if( strchr(name, ';') == 0)
      (void)remove(strcat(name,";-1"));
   tfree(name);
}
#else				/*  VMS				*/
#   ifdef UNRESTRICTED_FILE_NAMES
void	purge(char *file)
{
   int unlink(const char *);
   char *name = aalloc(strlen(file)+2,char);
   (void)strcat(strcpy(name, file), "%");

   (void)unlink(name);
   xfree(name);
}
#      else
/*ARGSUSED*/
void	purge(char *file)
{
}
#   endif			
#endif				/*  VMS				*/

/******************************************************************************
 * save_version_inc. Add a version number to save_file name                   *
 ******************************************************************************/
#ifdef UNRESTRICTED_FILE_NAMES
void save_version_inc(char *save_name, int namesize)
{
   int len = strlen(save_name), i, vsn=0;
   FILE *rst;
   
   for(i=len-1; isdigit(save_name[i]); i--)
      ;
   if( save_name[i] != '.' )
   {
      if( namesize - len >= 2 ) 
      {
	 (void)strcat(save_name,".0");
	 len +=2;
	 i++;
      }
      else
      {
	 /* No room to add extension.  Return and give up? */
	 return;
      }
   }
   if( i < len-1)
      vsn=atoi(save_name+i+1)+1;

   (void)sprintf(save_name+i+1,"%d", vsn);
   /*
    * Better not overwrite a file which exists.
    */
   while( rst = fopen(save_name,"rb") )
   {
      vsn++;
      (void)sprintf(save_name+i+1,"%d", vsn);
      fclose(rst);
   }
}
#else
void save_version_inc(char *save_name, int namelen)
{
}
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

double	err_fn(double x)
{
   double t = 1.0/(1.0 + PP * x);
   if ( x >= 0.0 )
      return( 1.0 - POLY5(t) * exp(-x*x) );
   else
      return( -err_fn(-x) );
}
/******************************************************************************
 * inhibit_vectorization().  Self explanatory dummy function		      *
 ******************************************************************************/
void inhibit_vectorization(void)
{
}
