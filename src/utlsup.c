#ifndef lint
static char *RCSid = "$Header: /usr/users/kr/CVS/moldy/src/utlsup.c,v 1.15.8.2 2003/07/30 09:31:30 moldydv Exp $";
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
#endif

#define DATAREC "%d record%ssuccessfully read from %s"

void mat_vec_mul(real (*m)[3], vec_mp in_vec, vec_mp out_vec, int number);
void invert(real (*a)[3], real (*b)[3]);

char	*comm;
#ifdef USE_XDR
   XDR          xdrs;
#endif
/******************************************************************************
 * Dummies of moldy routines so that utils may be linked with moldy library   *
 ******************************************************************************/
/*VARARGS*/
void 	init_rdf(void)
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
void	new_line(void)
{
   (void)putchar('\n');
}
/*VARARGS*/
void	banner_page(void)
{}
/*VARARGS*/
void	note(void)
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
 * mystrdup().  Routine for copying one string to another.                    *
 ******************************************************************************/
char * mystrdup(char *s)
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
#define STRLEN 80
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
         s = fgets(line, len, file);               /* Read one line of input     */
         if(s == NULL) break;                      /* exit if end of file        */
         i = strlen(s) - 1;
         while(i >= 0 && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n'))
            s[i--] = '\0';                         /* Strip trailing white space */
      }
      while(*s == '\0' || *s == '#');              /* Repeat if blank or comment */
   }
   else
   {
      s = fgets(line, len, file);               /* Read one line of input     */
      if(s != NULL)                             /* ignore if end of file        */
      {
         i = strlen(s) - 1;
         while(i >= 0 && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n'))
            s[i--] = '\0';                         /* Strip trailing white space */
      }
   }

   if(s == NULL)
      *line = '\0';                             /* Return null at eof         */
   return(line);
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
 *  read_dump_header. Read the header of a moldy dump file.                   *
 ******************************************************************************/
int read_dump_header(char *fname, FILE *dumpf, dump_mt *hdr_p, boolean *xdr_write,
		size_mt sysinfo_size, dump_sysinfo_mt *dump_sysinfo)
{
   int      errflg = true;      /* Provisionally !!   */
   char     vbuf[sizeof hdr_p->vsn + 1];
   int      vmajor,vminor;

   *xdr_write = false;
#ifdef USE_XDR
   /*
       * Attempt to read dump header in XDR format
       */
   if( xdr_dump(&xdrs, hdr_p) )
   {
      strncpy(vbuf,hdr_p->vsn,sizeof hdr_p->vsn);
      vbuf[sizeof hdr_p->vsn] = '\0';
      if( strstr(vbuf,"(XDR)") )
      {
         errflg = false;
         *xdr_write = true;
      }
   }
#endif
   /*
    * If we failed, try to read header as native struct image.
    */
   if( ! *xdr_write )
   {
      if( fseek(dumpf, 0L, 0) )
         message(NULLI, NULLP, WARNING, SEFAIL, fname, strerror(errno));
      else if( fread((gptr*)&*hdr_p, sizeof(dump_mt), 1, dumpf) == 0 )
         message(NULLI, NULLP, WARNING, DRERR, fname, strerror(errno));
      else
         errflg = false;
   }
   if( ! errflg )
   {
      /*
       * Parse header version
       */
      errflg = true;
      if( sscanf(hdr_p->vsn, "%d.%d", &vmajor, &vminor) < 2 )
         message(NULLI, NULLP, WARNING, INDVSN, hdr_p->vsn);
      if( vmajor < 2 || vminor <= 17)
         message(NULLI, NULLP, WARNING, OLDVSN, hdr_p->vsn);
      else
         errflg = false;
   }
   if( errflg ) return errflg;

   if( dump_sysinfo == 0)
      return errflg;
   else if ( sysinfo_size == sizeof(dump_sysinfo_mt) )
   {
      /*
       * Now check for sysinfo and read fixed part of it.  This is needed to
       * determine species count to read the whole thing.
       */
#ifdef USE_XDR
      if( *xdr_write ) {
         if( ! xdr_dump_sysinfo_hdr(&xdrs, dump_sysinfo) )
            message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
         errflg = false;
      } else
#endif
      {
         if( fread((gptr*)dump_sysinfo,sizeof(dump_sysinfo_mt), 1, dumpf) == 0)
            message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
         errflg = false;
      }
   }
   else
   {
      /*
       * Now check for sysinfo and read it all.  N.B.  Buffer must be
       * allocated to full expected size by prior call to read_dump_header.
       */
#ifdef USE_XDR
      if( *xdr_write ) {
         if( ! xdr_dump_sysinfo(&xdrs, dump_sysinfo, vmajor, vminor) )
            message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
      if (sizeof(dump_sysinfo_mt)
          + sizeof(mol_mt) * (dump_sysinfo->nspecies-1) > sysinfo_size)
      {
         /*
          * We have already overrun the end of the "dump_sysinfo" buffer.
          * Perhaps we can exit gracefully before crashing?
          */
         message(NULLI, NULLP, FATAL, RDHERR,  sizeof(dump_sysinfo_mt)
                 + sizeof(mol_mt) * (dump_sysinfo->nspecies-1), sysinfo_size);
      }
         errflg = false;
      } else
#endif
      {
         if( fread((gptr*)dump_sysinfo, sysinfo_size, 1, dumpf) == 0)
            message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
         errflg = false;
      }
   }
   return errflg;
}
