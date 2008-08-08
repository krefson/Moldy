#ifndef lint
static char *RCSid = "$Header: /home/moldy/CVS/moldy/src/utlsup.c,v 1.20 2005/03/06 18:27:13 cf Exp $";
#endif
#include "defs.h"
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include "structs.h"
#include "messages.h"
#include "specdata.h"
#ifdef USE_XDR
#   include     "xdr.h"
   XDR          xdrs;
#endif

char	*comm;

#define STRLEN 80

void invert (real (*a)[3], real (*b)[3]);
void mat_vec_mul (real (*m)[3], vec_mp in_vec, vec_mp out_vec, int number);
/******************************************************************************
 * Dummies of moldy routines so that utils may be linked with moldy library   *
 ******************************************************************************//*VARARGS*/
void    init_rdf(void)
{}
/*VARARGS*/
gptr *rdf_ptr(void)
{return 0;}
/*VARARGS*/
void new_lins(void)
{}
int lines_left(void)
{return 0;}
void new_page(void)
{}
void    new_line(void)
{
   (void)putchar('\n');
}
/*VARARGS*/
void    banner_page(void)
{}
/*VARARGS*/
void    note(void)
{}
/******************************************************************************
 *  message.   Deliver error message to possibly exiting.  It can be called   *
 *             BEFORE output file is opened, in which case output to stderr.  *
 ******************************************************************************/
/*VARARGS*/
void    message(int *nerrs, ...)
{
   va_list      ap;
   char         *buff;
   int          sev;
   char         *format;

   va_start(ap, nerrs);
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
/*VARARGS*/
void    error(char *format, ...)
{
   va_list	ap;
   va_start(ap, format);

   (void)fprintf(stderr, "%s: ",comm);
   (void)vfprintf(stderr, format, ap);
   fputc('\n',stderr);
   va_end(ap);

   exit(3);
}
/******************************************************************************
 * trim().  Trim leading and trailing spaces from string variables            *
 ******************************************************************************/
char    *trim(char *s)
{
char    *t = s, *p;

   /* Remove leading spaces */
   while( *t == ' ' && *t != '\0')
      t++;

   strcpy(s, t);

   /* Now remove trailing spaces */
   p = s + strlen(s);
   while( p > s && !isalnum(*p))
      *p-- = '\0';

   return (s);
}
/******************************************************************************
 * mystrdup().  Routine for copying one string to another.                    *
 ******************************************************************************/
char * mystrdup(const char *s)
{
   char * t = NULL;
   if(s) t=malloc(strlen(s)+1);
   return t?strcpy(t,s):0;
}
/******************************************************************************
 *  Tokenise().  Parse the string of fields to be returned and return a mask  *
 *  in a char array.  Format is 1,3,6-9,3 . . . ie comma-separated with cont- *
 *  iguous range specified with hyphen.  Numbering starts at 1.               *
 ******************************************************************************/
int     tokenise(char *fields, char *mask, int len)
{
   char *s;
   int  lo, hi, i, n;

   for(i = 0; i < len; i++)
      mask[i] = 0;

   while( ( s = strtok(fields,",") ) != NULL )
   {
      n = sscanf(s, "%d-%d", &lo, &hi);
      if( n == 0 )
         return 0;

      if( n == 1 )
         hi = lo;

      if( lo < 1 || hi < lo || hi > len)
         return 0;

      for( i = lo-1; i < hi; i++)
         mask[i] = 1;
      fields = NULL;
   }
   return 1;
}
/******************************************************************************
 * get_tokens(). Routine for breaking down string into substrings.            *
 *               Returns no of substrings (not including comments).           *
 ******************************************************************************/
int   get_tokens(char *buf, char **linev, char *sep)
{
   char *p;
   int linec = 0;

   while( (p = strtok(buf,sep) ) != NULL)
   {
      if( strncmp(p,"#",1) == 0) /* Ignore comments */
         break;
      *(linev++) = p;
      linec++;
      buf = NULL;
   }

   return linec;
}
/******************************************************************************
 * get_int().  Read an integer from stdin, issuing a prompt and checking      *
 * validity and range.  Loop until satisfied, returning EOF if appropriate.   *
 ******************************************************************************/
int get_int(char *prompt, int lo, int hi)
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
real get_real(char *prompt, real lo, real hi)
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
      return((real)EOF);
}
/******************************************************************************
 * get_sym().  Read a character from stdin and match to supplied set          *
 ******************************************************************************/
int get_sym(char *prompt, char *cset)
{
   char		ans_c, ans_str[STRLEN];
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
char	*get_str(char *prompt)
{
   char		ans_str[STRLEN];
   char		*str = malloc(STRLEN);
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
forstr(char *instr, int *start, int *finish, int *inc)
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
 *  get_line  read an input line with or without skipping blank/comment lines *
 ******************************************************************************/
char    *get_line(char *line, int len, FILE *file, int skip)
{
   char *s;
   int  i;

   if( skip )
   {
      do
      {
         s = fgets(line, len, file);            /* Read one line of input     */
         if(s == NULL) break;                   /* exit if end of file        */
         i = strlen(s) - 1;
         while(i >= 0 && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n'))
            s[i--] = '\0';                      /* Strip trailing white space */
      }
      while(*s == '\0' || *s == '#');           /* Repeat if blank or comment */
   }
   else
   {
      s = fgets(line, len, file);               /* Read one line of input     */
      if(s != NULL)                             /* ignore if end of file      */
      {
         i = strlen(s) - 1;
         while(i >= 0 && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n'))
            s[i--] = '\0';                      /* Strip trailing white space */
      }
   }

   if(s == NULL)
      *line = '\0';                             /* Return null at eof         */
   return(line);
}

/* vdotf.  Dot product between two vector floats */
double vdotf(int n, float *x, int ix, float *y, int iy)
{
   register double      dot=0.0;
   register int i, j;
   if( ix == iy && ix == 1)
   {
      for(i = 0; i < n; i++)
         dot += x[i] * y[i];
   }
   else
   {
      for(i = j = 0; i < n*ix; i += ix)
      {
         dot += x[i] * y[j];
         j += iy;
      }
   }
   return(dot);
}

/******************************************************************************
 * dump_to_moldy.  Fill the 'system' arrays with the dump data in 'buf' (see  *
 * dump.c for format), expanding floats to doubles if necessary.              *
 ******************************************************************************/
void
dump_to_moldy(float *buf, system_mt *system)
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
 *  dump_info. Extract no of slices and dump level from header info file.     *
 ******************************************************************************/
int dump_info(FILE *Fp, int *level)
{
   int      numslice=0, idata;
   char     line[LLEN];

   while( !feof(Fp) )
   {
     get_line(line,LLEN,Fp,1);

     if( sscanf(line, "Dump level\t\t\t= %d", &idata) > 0)
       *level = idata;
     if( sscanf(line, "Number of dumps\t\t\t= %d", &idata) > 0)
       numslice += idata;
   }

   numslice--;

   return numslice;
}
/******************************************************************************
 * traj_con().  Connect molecular c_of_m`s into continuous trajectories      *
 ******************************************************************************/
void
traj_con(system_mt *system, vec_mt (*prev_cofm), int n)
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
traj_con2(spec_mt *species, vec_mt (*prev_cofm), vec_mt (*traj_cofm), int nspecies)
{
   spec_mt	*spec;
   int		i, imol, totmol=0;
 
   for( spec = species; spec < species+nspecies; spec++)
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
range_in(system_mt *system, real (*range)[3])
{
   double	box[3], blimit;
   mat_mp	h = system->h;
   int		i;

   for( i=0; i<3; i++)
   {
       box[i] = sqrt(SQR(h[0][i]) + SQR(h[1][i]) + SQR(h[2][i]));
       blimit = MAX3(box[0], box[1], box[2]);

       if( range[i][2] )
       {
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
          range[i][2]--;
       }
       else
       {
          range[i][0] = -1.0*blimit;
          range[i][1] = blimit;
       }
   }
   return 0;
}
/******************************************************************************
 * in_region. Function for checking if atom lies within/out region            *
 ******************************************************************************/
int
in_region(real *pos, real (*range)[3])
{
int	xi, yi, zi;

    xi = ( pos[0] >= range[0][0] && pos[0] <= range[0][1] ) ? 0 : 1;
    yi = ( pos[1] >= range[1][0] && pos[1] <= range[1][1] ) ? 0 : 1;        
    zi = ( pos[2] >= range[2][0] && pos[2] <= range[2][1] ) ? 0 : 1;        

    if( (range[0][2] == xi) && (range[1][2] == yi) && (range[2][2] == zi) )
       return 1;   /* Molecule c_of_m lies within selected region */
    else
       return 0;   /* Molecule c_of_m isn't within selected region */ 
}
/******************************************************************************
 *  open_dump(). Open a moldy dump file for read or write.                    *
 ******************************************************************************/
FILE  *open_dump(char *fname, char *mode)
{
   FILE *dumpf;

   dumpf = fopen(fname, mode);

#ifdef USE_XDR
   if( dumpf )
   {
      if( mode[0] == 'w' || (mode[0] && mode[1] == '+') ||  (mode[1] && mode[2] == '+'))
         xdrstdio_create(&xdrs, dumpf, XDR_ENCODE);
      else
         xdrstdio_create(&xdrs, dumpf, XDR_DECODE);
   }
#endif
    return dumpf;
}

int close_dump(FILE *dumpf)
{
#ifdef USE_XDR
   xdr_destroy(&xdrs);
#endif
   return fclose(dumpf);
}
int rewind_dump(FILE *dumpf, int xdr)
{
#ifdef USE_XDR
   if( xdr )
      xdr_setpos(&xdrs, 0);
#endif
   return fseek(dumpf, 0L, SEEK_SET);
}
/******************************************************************************
 *  check_control. Check if file contains control info.                       *
 ******************************************************************************/
int check_control(FILE *file)
{
  int num_items;
  char  line[LLEN], name[LLEN], value[LLEN];

  get_line(line, LLEN, file, 1);
  num_items = sscanf(line, " %[^= ] = %127[^#]", name, value);
  rewind(file);
  if( num_items == 2 )
    return 1; /* Control info present */
  else
    return 0; /* No control info present */
}
/******************************************************************************
 * str_cut().  Separates alphabetical and numeric parts of a string.          *
 ******************************************************************************/
int     str_cut(char *in, char *out) /* Input and output strings must be different */
{
   int          i,j,k=strlen(in);
   int          value;

   for( i=0; i < k; i++)
     if( isdigit(in[i]) || in[i] == '-' )
        break;

   for( j=0; j < k; j++)
     if( isalpha(in[j]) )
        break;

   if( i <= j )
   {
      strncpy(out,in+j,k-j);   /* Copy alphabetical part */
      in[j] ='\0';             /* Truncate alphabetical part */
      value=atoi(in);          /* Convert numeric part */
   }
   else
   {
      strncpy(out,in,i);       /* Copy alphabetical part */
      value=atoi(in+i);        /* Convert numeric part */
   }
   if( strrchr(in,'-') != NULL && value > 0 )
        value *= -1;

   return value;
}
