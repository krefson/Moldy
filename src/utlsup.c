#ifndef lint
static char *RCSid = "$Header: /home/eeyore_data/keith/moldy/src/RCS/utlsup.c,v 1.2 1999/10/29 16:44:28 keith Exp $";
#endif

#include "defs.h"
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#include <errno.h>
#include <math.h>
#include "stdlib.h"
#include "stddef.h"
#include "string.h"
#include <stdio.h>
#include "structs.h"
#include "messages.h"

void mat_vec_mul();
void invert();

char	*comm;
/******************************************************************************
 * Dummies of moldy routines so that utils may be linked with moldy library   *
 ******************************************************************************/
void 	init_rdf()
{}
gptr *rdf_ptr()
{return 0;}
void new_lins()
{}
int lines_left()
{return 0;}
void new_page()
{}
void	new_line()
{
   (void)putchar('\n');
}
void	banner_page()
{}
void	note()
{}

/******************************************************************************
 *  message.   Deliver error message to possibly exiting.  It can be called   *
 *             BEFORE output file is opened, in which case output to stderr.  *
 ******************************************************************************/
#ifdef HAVE_STDARG_H
#   undef  va_alist
#   define      va_alist int *nerrs, ...
#   ifdef va_dcl
#      undef va_dcl
#   endif
#   define va_dcl /* */
#endif
/*VARARGS*/
void    message(va_alist)
va_dcl
{
   va_list      ap;
   char         *buff;
   int          sev;
   char         *format;
   static char  *sev_txt[] = {" *I* "," *W* "," *E* "," *F* "};
#ifdef HAVE_STDARG_H
   va_start(ap, nerrs);
#else
   int          *nerrs;

   va_start(ap);
   nerrs = va_arg(ap, int *);
#endif
   buff  = va_arg(ap, char *);
   sev   = va_arg(ap, int);
   format= va_arg(ap, char *);

   (void)fprintf(stderr,"%s: ",comm);
   (void)vfprintf(stderr, format, ap);
   va_end(ap);
   fputc('\n',stderr);

   if(buff != NULL)                     /* null ptr means don't print buffer  */
   {
      (void)fprintf(stderr,"     buffer contents=\"%s\"",buff);
      fputc('\n',stderr);
   }
   if(sev >= ERROR && nerrs != NULL)
      (*nerrs)++;
   if(sev == FATAL)
      exit(3);
}

/******************************************************************************
 *  error.   Deliver error message to possibly exiting.                       *
 ******************************************************************************/
#ifdef HAVE_STDARG_H
#undef  va_alist
#define	va_alist char *format, ...
#ifdef  va_dcl
#   undef  va_dcl
#endif
#define va_dcl /* */
#endif
/*VARARGS*/
void	error(va_alist)
va_dcl
{
   va_list	ap;
#ifdef HAVE_STDARG_H
   va_start(ap, format);
#else
   char		*format;

   va_start(ap);
   format= va_arg(ap, char *);
#endif

   (void)fprintf(stderr, "%s: ",comm);
   (void)vfprintf(stderr, format, ap);
   fputc('\n',stderr);
   va_end(ap);

   exit(3);
}
/******************************************************************************
 * mystrdup().  Routine for copying one string to another.                    *
 ******************************************************************************/
char * mystrdup(s)
char *s;
{
   char * t = NULL;
   if(s) t=malloc(strlen(s)+1);
   return t?strcpy(t,s):0;
}
/******************************************************************************
 * get_int().  Read an integer from stdin, issuing a prompt and checking      *
 * validity and range.  Loop until satisfied, returning EOF if appropriate.   *
 ******************************************************************************/
int get_int(prompt, lo, hi)
char	*prompt;
int	lo, hi;
{
   char		ans_str[80];
   int		ans_i, ans_flag;

   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, "%d", &ans_i) == 1 && ans_i >= lo && ans_i <= hi)
	 ans_flag++;
   }
   if( ans_flag )
      return(ans_i);
   else
      return(EOF);
}
/******************************************************************************
 * get_real().  Read a real from stdin, issuing a prompt and checking         *
 * validity and range.  Loop until satisfied, returning EOF if appropriate.   *
 ******************************************************************************/
real get_real(prompt, lo, hi)
char	*prompt;
real 	lo, hi;
{
   char		ans_str[80];
   real		ans_r; 
   int		ans_flag;
   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, "%lf", &ans_r) == 1 && ans_r >= lo && ans_r <= hi)
	 ans_flag++;
   }
   if( ans_flag )
      return(ans_r);
   else
      return(EOF);
}
/******************************************************************************
 * get_sym().  Read a character from stdin and match to supplied set          *
 ******************************************************************************/
int get_sym(prompt, cset)
char	*prompt;
char	*cset;
{
   char		ans_c, ans_str[80];
   int		ans_flag;

   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, " %c", &ans_c) == 1 && strchr(cset, ans_c))
	 ans_flag++;
   }
   if( ans_flag )
      return(ans_c);
   else
      return(EOF);
}
/******************************************************************************
 * get_str().  Read a string from stdin, issuing a prompt                     *
 ******************************************************************************/
char	*get_str(prompt)
char	*prompt;
{
   char		ans_str[80];
   char		*str = malloc(80);
   int		ans_flag;

   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, "%s", str) == 1)
	 ans_flag++;
   }
   if( ans_flag )
      return(str);
   else
      return(NULL);
}
/******************************************************************************
 * forstr.  Parse string str of format s-f:n (final 2 parts optional),        *
 *          returning integer values of s,f,n. f defaults to s and n to 1     *
 ******************************************************************************/
int
forstr(instr, start, finish, inc)
char	*instr;
int	*start, *finish, *inc;
{
   char	*p, *pp, *str = mystrdup(instr);
   
   if( (p = strchr(str,':')) != NULL)
   {
      *inc = strtol(p+1, &pp, 0);
      if( pp == p+1 )
	 goto limerr;
      if(*pp)
	 goto unrec;
      *p = 0;
   }
   else
      *inc = 1;
   if( (p = strchr(str,'-')) != NULL)
   {
      *p = 0;
      *start = strtol(str, &pp, 0);
      if( pp == str )
	 goto limerr;
      if(*pp)
	 goto unrec;
      *finish = strtol(p+1, &pp, 0);
      if( pp == p+1 )
	 goto limerr;
      if(*pp)
	 goto unrec;
   }
   else
   {
      *start = *finish = strtol(str, &pp, 0);
      if( pp == str )
	 goto limerr;
      if(*pp)
	 goto unrec;
   }
   if( *start > *finish || *start < 0 || *inc <= 0 )
   {
      fputs("Limits must satisfy", stderr);
      fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
      goto limerr;
   }
   return 0;
 unrec:
   fprintf(stderr, "Unrecognised characters in range string \"%s\"\n", pp);
 limerr:
   return -1;
}
/******************************************************************************
 * dump_to_moldy.  Fill the 'system' arrays with the dump data in 'buf' (see  *
 * dump.c for format), expanding floats to doubles if necessary.              *
 ******************************************************************************/
#define DUMP_SIZE(level)  (( (level & 1) + (level>>1 & 1) + (level>>2 & 1) ) * \
           (3*sys.nmols + 4*sys.nmols_r + 9)+ (level>>3 & 1) * \
           (3*sys.nmols + 3*sys.nmols_r + 9) + (level & 1))
void
dump_to_moldy(buf, system)
float	*buf;
system_mt *system;
{
   int i;
   float	*c_of_m = buf;
   float	*quat   = buf+3*system->nmols;
   float	*h      = buf+3*system->nmols + 4*system->nmols_r;
   mat_mt	hinv;

/* $dir no_recurrence */
   for(i = 0; i < system->nmols; i++)
   {
      system->c_of_m[i][0] = c_of_m[3*i];
      system->c_of_m[i][1] = c_of_m[3*i+1];
      system->c_of_m[i][2] = c_of_m[3*i+2]; 
   }
/* $dir no_recurrence */
   for(i = 0; i < system->nmols_r; i++)
   {
      system->quat[i][0] = quat[4*i];
      system->quat[i][1] = quat[4*i+1];
      system->quat[i][2] = quat[4*i+2];
      system->quat[i][3] = quat[4*i+3];
   }
/* $dir no_recurrence */
   for(i = 0; i < 3; i++)
   {
      system->h[i][0] = h[3*i];
      system->h[i][1] = h[3*i+1];
      system->h[i][2] = h[3*i+2];
   }
   invert(system->h, hinv);
   mat_vec_mul(hinv, system->c_of_m, system->c_of_m, system->nmols);
}
/******************************************************************************
 * traj_con().  Connect molecular c_of_m`s into continuous trajectories      *
 ******************************************************************************/
void
traj_con(system, prev_cofm, n)
system_mt       *system;
vec_mt          *prev_cofm;
int		n;
{
    int		i, imol;

    for( imol = 0; imol < system->nmols; imol++ )
       for( i = 0; i < 3; i++)
       {
          if( n > 0 ) 
              system->c_of_m[imol][i] = system->c_of_m[imol][i]
                     - floor(system->c_of_m[imol][i]-prev_cofm[imol][i]+0.5);
          prev_cofm[imol][i] = system->c_of_m[imol][i];
       }
}
/******************************************************************************
 * traj_con2().  Connect molecular c_of_m`s into continuous trajectories      * 
 ******************************************************************************/
void
traj_con2(species, prev_cofm, traj_cofm, sp_range)
spec_mt		species[];
vec_mt		*prev_cofm;
vec_mt		*traj_cofm;
int		sp_range[3];
{
   spec_mt	*spec;
   int		i, imol, totmol=0;
 
   for( spec = species+sp_range[0]; spec <= species+sp_range[1]; spec+=sp_range[2])
     for( imol = 0; imol < spec->nmols; totmol++, imol++)
	if( prev_cofm == 0 ) 
	   for( i = 0; i < 3; i++)
	      traj_cofm[totmol][i] = spec->c_of_m[imol][i];
	else 
	   for( i = 0; i < 3; i++)
	      traj_cofm[totmol][i] = spec->c_of_m[imol][i] - 
		 floor(spec->c_of_m[imol][i]-prev_cofm[totmol][i]+0.5);
}
/******************************************************************************
 * range_in. Function for determining ranges of atom positions to include     *
 ******************************************************************************/
int
range_in(system, range, range_flag)
system_mt	*system;
real		range[3][2];
int		*range_flag;
{
   double	box[3], blimit;
   mat_mp	h = system->h;
   int		i;

   for( i=0; i<3; i++)
   {
       box[i] = sqrt(SQR(h[0][i]) + SQR(h[1][i]) + SQR(h[2][i]));
       blimit = MAX3(box[0], box[1], box[2]);

       if( range_flag[i])
          switch(i)
          {
           case 0:
             range[i][0] = get_real("Enter x minimum: ", -1.0*blimit,blimit);
             range[i][1] = get_real("Enter x maximum: ", range[i][0],2.0*blimit);
             break;
           case 1:
             range[i][0] = get_real("Enter y minimum: ", -1.0*blimit,blimit);
             range[i][1] = get_real("Enter y maximum: ", range[i][0],2.0*blimit);
             break;
           case 2:
             range[i][0] = get_real("Enter z minimum: ", -1.0*blimit,blimit);
             range[i][1] = get_real("Enter z maximum: ", range[i][0],2.0*blimit);
             break;
          }    
       else
       {
          range[i][0] = -1.0*blimit;
          range[i][1] = blimit;
       }
   }
   return 0;
}
