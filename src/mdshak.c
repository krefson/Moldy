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
#ifndef lint
static char *RCSid = "$Header: /home/eeyore_data/keith/md/moldy/RCS/mdshak.c,v 2.15 1998/01/09 11:35:11 keith Exp $";
#endif

#include "defs.h"
#if defined(ANSI) || defined(__STDC__)
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
#include "ReadDCD.h"
#if defined(ANSI) || defined(__STDC__)
gptr	*arralloc(size_mt,int,...); 	/* Array allocator		      */
#else
gptr	*arralloc();	        	/* Array allocator		      */
#endif

void	invert();
void	mat_vec_mul();
void	make_sites();
char	*strlower();
void	read_sysdef();
void	initialise_sysdef();
void	re_re_header();
void	re_re_sysdef();
void	allocate_dynamics();
void	lattice_start();
void	read_restart();
void	init_averages();
int	getopt();
gptr	*talloc();
FILE	*popen();
/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
static   char		*comm;
contr_mt		control;

#define OUTBIN 2
#define SHAK   0
#define XYZ 1
#define DCD 3
/******************************************************************************
 * Dummies of 'moldy' routines so that mdshak may be linked with moldy library*
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
void	conv_potentials()
{}
void	conv_control()
{}
/******************************************************************************
 *  message.   Deliver error message to possibly exiting.  It can be called   *
 *             BEFORE output file is opened, in which case outt to stderr.    *
 ******************************************************************************/
#if defined(ANSI) || defined(__STDC__)
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
#if defined(ANSI) || defined(__STDC__)
   va_start(ap, nerrs);
#else
   int          *nerrs;

   va_start(ap);
   nerrs = va_arg(ap, int *);
#endif
   buff  = va_arg(ap, char *);
   sev   = va_arg(ap, int);
   format= va_arg(ap, char *);

   (void)fprintf(stderr, "%s: ",comm);
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
 *  message.   Deliver error message to possibly exiting. 		      *
 ******************************************************************************/
#if defined(ANSI) || defined(__STDC__)
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
#if defined(ANSI) || defined(__STDC__)
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
static char * mystrdup(s)
char *s;
{
   char * t=malloc(strlen(s)+1);
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
 * get_sym().  Read a character from stdin and match to suplied set	      *
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
 * get_str().  Read an string from stdin, issuing a prompt.		      *
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
 ******************************************************************************/
/******************************************************************************
 * forstr.  Parse string str of format s-f:n  (final 2 parts optional),       *
 *          returning integer values of s,f,n. f defaults to s and n to 1     *
 ******************************************************************************/
int
forstr(instr, start, finish, inc)
char	*instr;
int	*start, *finish, *inc;
{
   char	*p, *pp, *str = mystrdup(instr);
   long strtol();
   
   if( (p = strchr(str,':')) != NULL)
   {
      *inc = strtol(p+1, &pp, 0);
      if( pp == p+1 )
	 goto limerr;
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
      *finish = strtol(p+1, &pp, 0);
      if( pp == p+1 )
	 goto limerr;
   }
   else
   {
      *start = *finish = strtol(str, &pp, 0);
      if( pp == str )
	 goto limerr;
   }
   return 0;
 limerr:
   return -1;
}
/******************************************************************************
 * dump_to_moldy.  Fill the 'system' arrays with the dump data in 'buf' (see  *
 * dump.c for format), expanding floats to doubles if necessary.              *
 ******************************************************************************/
#define DUMP_SIZE(level)  (( (level & 1) + (level>>1 & 1) + (level>>2 & 1) ) * \
			            (3*sys.nmols + 4*sys.nmols_r + 9)+ \
			     (level>>3 & 1) * \
			            (3*sys.nmols + 3*sys.nmols_r + 9) +\
			     (level & 1))
void
dump_to_moldy(buf, system)
float	*buf;
system_mt *system;
{
   int i;
   float	*c_of_m = buf;
   float	*quat   = buf+3*system->nmols;
   float	*h      = buf+3*system->nmols + 4*system->nmols_r;
   mat_mt hinv;

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
 ******************************************************************************/
void mat_vec_mul3(m, vec, number)
int             number;         /* Number of vectors to be multiplied         */
real            m[3][3];        /* Matrix                                     */
real            **vec;          /* Output vector.  CAN BE SAME AS INPUT  (out)*/
{
   int i;
   register double        a0, a1, a2;
   
   for(i = 0; i < number; i++)
   {
      a0 = vec[0][i];  a1 = vec[1][i];  a2 = vec[2][i];
      
      vec[0][i] = m[0][0]*a0 + m[0][1]*a1 + m[0][2]*a2;
      vec[1][i] = m[1][0]*a0 + m[1][1]*a1 + m[1][2]*a2;
      vec[2][i] = m[2][0]*a0 + m[2][1]*a1 + m[2][2]*a2;
   }
}
/******************************************************************************
 * shakal_out().  Write a system configuration to stdout in the form of an    *
 * input data file for the graphics program SCHAKAL88.			      *
 ******************************************************************************/
void
schakal_out(n, system, species, site_info, insert)
int	n;
system_mt	*system;
spec_mt		species[];
site_mt		site_info[];
char		*insert;
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   spec_mt	*spec;
   double	a, b, c, alpha, beta, gamma;
   mat_mp	h = system->h;
   mat_mt	hinv;
   int		imol, isite, is;

   invert(h,hinv);

   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

   printf("CELL %f %f %f %f %f %f\n", a, b, c, alpha, beta, gamma);
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      make_sites(system->h, spec->c_of_m, spec->quat, spec->p_f_sites,
		 spec->framework, site, spec->nmols, spec->nsites);

      mat_vec_mul3(hinv, site, spec->nsites*spec->nmols);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 puts("MOL");
	 for(is = 0; is < spec->nsites; is++)
	 {
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
	       (void)printf("ATOM %-8s %7.4f %7.4f %7.4f\n",
			    site_info[spec->site_id[is]].name,
			    site[0][isite], site[1][isite], site[2][isite]);
	    isite++;
	 }
      }
   }

   if( insert != NULL)
      (void)printf("%s\n", insert);

   (void)printf("END %d\n", n);
   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*) site);
}
/******************************************************************************
 * xyz_out().  Write a system configuration to stdout in the form of an    *
 * input data file for the graphics program XYZ (rasmol -xyz file)	      *
 ******************************************************************************/
void
xyz_out(n, system, species, site_info, insert)
int	n;
system_mt	*system;
spec_mt		species[];
site_mt		site_info[];
char		*insert;
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   spec_mt	*spec;
   int		imol, isite, is;

/* We count the number of atoms */
   isite=0;
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
         {
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
              isite++;
         }
      }
   }
/* Now we write the Xyz header */
   (void)printf("%d\n",isite);
/* It would be nice to have here the real title */
   (void)printf("%s\n",control.title);
   
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      make_sites(system->h, spec->c_of_m, spec->quat, spec->p_f_sites,
		 spec->framework, site, spec->nmols, spec->nsites);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
	 {
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
	       (void)printf("%-8s %7.4f %7.4f %7.4f\n",
			    site_info[spec->site_id[is]].name,
			    site[0][isite], site[1][isite], site[2][isite]);
	    isite++;
	 }
      }
   }

   if( insert != NULL)
      (void)printf("%s\n", insert);

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*) site);
}
/******************************************************************************
 * dcd_out().  Write a system configuration to stdout in the form of an    *
 * DCD data file for the graphics program VMD                              *
 ******************************************************************************/
void
dcd_out(n, irec, inc, system, species, site_info, insert)
int	n, irec, inc;
system_mt	*system;
spec_mt		species[];
site_mt		site_info[];
char		*insert;
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   float	**sitef = (float**)arralloc(sizeof(float),2,
					    0,2,0,system->nsites-1);
   spec_mt	*spec;
   int		isite, is, i, imol, isitem;

   isitem=0; 
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      make_sites(system->h, spec->c_of_m, spec->quat, spec->p_f_sites,
		 spec->framework, site, spec->nmols, spec->nsites);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
	 {
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
	    {
	       for(i=0; i<3; i++)
		  sitef[i][isitem] = site[i][isite];
	       isitem++;
	    }
	    isite++;
	 }
      }
   }
/* On first call  write the DCD header.  Always write to stdout. */
   if( n == 0 )
      write_dcdheader(stdout, control.title, isitem, irec, 0, inc, control.step);
   
   write_dcdstep(stdout, isitem, sitef[0], sitef[1], sitef[2]);

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*)site);
   afree((gptr*)sitef);
}
/******************************************************************************
 * atoms_out().  Write a system configuration to stdout in the form of an     *
 * binary atomic co-ordinates.						      *
 ******************************************************************************/
void
atoms_out(n, system, species, site_info, atom_sel)
int	n;
system_mt	*system;
spec_mt		species[];
site_mt		site_info[];
char		*atom_sel;
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   spec_mt	*spec;
   float	fsite[3];
   mat_mp	h = system->h;
   mat_mt	hinv;
   int		imol, isite, is,i;

   invert(h,hinv);

   for(spec = species; spec < species+system->nspecies; spec++)
   {
      make_sites(system->h, spec->c_of_m, spec->quat, spec->p_f_sites,
		 spec->framework, site, spec->nmols, spec->nsites);

      mat_vec_mul3(hinv, site, spec->nsites*spec->nmols);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
	 {
	    for(i=0; i<3; i++)
	       fsite[i]=site[i][isite];
            fwrite((char*)fsite, sizeof fsite, 1, stdout);
	    isite++;
	 }
      }
   }
   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*) site);
}
 /******************************************************************************
 * Centre_mass.  Shift system centre of mass to origin (in discrete steps),   *
 ******************************************************************************/
void
centre_mass(species, nspecies, c_of_m)
spec_mt		species[];
int		nspecies;
vec_mt		c_of_m;
{
   double	mass;
   spec_mt	*spec;
   int		imol;
   vec_mt	*s_c_of_m;

   mass = c_of_m[0] = c_of_m[1] = c_of_m[2] = 0.0;
   for(spec = species; spec < species + nspecies; spec++ )
   {
      s_c_of_m = spec->c_of_m;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 c_of_m[0] += spec->mass*s_c_of_m[imol][0];
	 c_of_m[1] += spec->mass*s_c_of_m[imol][1];
	 c_of_m[2] += spec->mass*s_c_of_m[imol][2];
      }
      mass += spec->nmols*spec->mass;
   }

   c_of_m[0] /= mass;
   c_of_m[1] /= mass;
   c_of_m[2] /= mass;
   c_of_m[0] = floor(c_of_m[0]+0.5);
   c_of_m[1] = floor(c_of_m[1]+0.5);
   c_of_m[2] = floor(c_of_m[2]+0.5);
}
/******************************************************************************
 * Shift.  Translate all co-ordinates.					      *
 ******************************************************************************/
void	shift(r, nmols, s)
vec_mt	r[];
int	nmols;
vec_mt	s;
{
   int imol;
   for(imol = 0; imol < nmols; imol++)
   {
      r[imol][0] -= s[0];
      r[imol][1] -= s[1];
      r[imol][2] -= s[2];
   }
}
/******************************************************************************
 * moldy_out.  Select output routine and handle file open/close		      *
 * Translate system relative to either centre of mass of posn of framework.   *
 ******************************************************************************/
void
moldy_out(n, irec, inc, system, species, site_info, atom_sel, outsw, insert)
int	n, irec, inc;
system_mt	*system;
spec_mt		species[];
site_mt		site_info[];
char		*atom_sel;
int		outsw;
char		*insert;
{
   spec_mp	spec, frame_spec  = NULL;
   vec_mt	c_of_m;
   
   for(spec = species; spec < species+system->nspecies; spec++)
      if( spec->framework )
	 frame_spec = spec;

   if( frame_spec != NULL )
      shift(system->c_of_m, system->nmols, frame_spec->c_of_m[0]);
   else
   {
      centre_mass(species, system->nspecies, c_of_m);
      shift(system->c_of_m, system->nmols, c_of_m);
   }
   switch (outsw)
   {
   case DCD:
      dcd_out(n, irec, inc, system, species, site_info, insert);
      break;
    case SHAK:
      schakal_out(n, system, species, site_info, insert);
      break;
    case XYZ:
      xyz_out(n, system, species, site_info, insert);
      break;
    case OUTBIN:
      atoms_out(n, system, species, site_info, atom_sel);
      break;
   }
}
/******************************************************************************
 * main().   Driver program for generating SCHAKAL input files from MOLDY     *
 * files.    Acceptable inputs are sys-spec files, or restart files. Actual   *
 * configrational info can be read from dump files, lattice-start files or    *
 * restart files.  Call: mdshak [-s sys-spec-file] [-r restart-file].   If    *
 * neither specified on command line, user is interrogated.		      *
 ******************************************************************************/
int
main(argc, argv)
int	argc;
char	*argv[];
{
   int	c, cflg = 0, ans_i, sym, data_source = 0;
   char 	line[80];
   extern char	*optarg;
   int		errflg = 0;
   int		intyp = 0;
   int		start, finish, inc;
   int		rflag;
   int		irec;
   int		iout = 0;
   int		outsw=0;
   char		*filename = NULL, *dump_name = NULL;
   char		*dumplims = NULL, *atom_sel = NULL;
   char		*insert = NULL;
   char		*tempname;
   char		dumpcommand[256];
   int		dump_size;
   float	*dump_buf;
   FILE		*Fp, *Dp;
   restrt_mt	restart_header;
   system_mt	sys;
   spec_mt	*species;
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf;
   int		av_convert;
   
#define MAXTRY 100
   control.page_length=1000000;

   comm = argv[0];
   if( strstr(comm, "mdshak") )
     outsw = SHAK;
   else if (strstr(comm, "mdxyz") )
     outsw = XYZ;
   else if (strstr(comm, "mddcd") || strstr(comm, "mdvmd") )
     outsw = DCD;
   else
     outsw = OUTBIN;

   while( (c = getopt(argc, argv, "a:bo:cr:s:d:t:i:xhv") ) != EOF )
      switch(c)
      {
       case 'a': 
	 atom_sel = optarg;
	 break;
       case 'b':
	 outsw=OUTBIN;
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error("failed to open file \"%s\" for output", optarg);
	 break;
       case 'c':
	 cflg++;
	 break;
       case 'r':
	 if( intyp )
	    errflg++;
	 intyp = data_source = c;
	 filename = optarg;
	 break;
       case 's':
	 if( intyp )
	    errflg++;
	 intyp = data_source = c;
	 filename = optarg;
	 break;
       case 'd':
	 dump_name = optarg;
	 break;
       case 't':
	 dumplims = optarg;
	 break;
       case 'i':
	 insert = optarg;
	 break;
       case 'x':
	 outsw = XYZ;
	 break;
       case 'h':
	 outsw = SHAK;
	 break;
       case 'v':
	 outsw = DCD;
	 break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fprintf(stderr,
	      "Usage: %s [-x|-v] [-h] [-c] [-s sys-spec-file|-r restart-file] ",
	      comm);
      fputs("[-d dump-files] [-t s[-f[:n]]] [-o output-file]\n", stderr);
      exit(2);
   }

   if( dump_name )
      data_source = 'd';

   if(intyp == 0)
   {
      fputs("How do you want to  specify the simulated system?\n", stderr);
      fputs("Do you want to use a system specification file (1)", stderr);
      fputs(" or a restart file (2)\n", stderr);
      if( (ans_i = get_int("? ", 1, 2)) == EOF )
	 exit(2);
      intyp = ans_i-1 ? 'r': 's';
      if( intyp == 's' )
      {
	 fputs( "Do you need to skip 'control' information?\n", stderr);
	 if( (sym = get_sym("y or n? ","yYnN")) == 'y' || sym == 'Y')
	    cflg++;
      }

      if( (filename = get_str("File name? ")) == NULL )
	 exit(2);
   }

   switch(intyp)
   {
    case 's':
      if( (Fp = fopen(filename,"r")) == NULL)
	 error("Couldn't open sys-spec file \"%s\" for reading", filename);
      if( cflg )
      {
	 do
	 {
	    fscanf(Fp, "%s",line);
	    (void)strlower(line);
	 }
	 while(! feof(stdin) && strcmp(line,"end"));
      }
      read_sysdef(Fp, &sys, &species, &site_info, &potpar);
      qpf = qalloc(sys.nspecies);
      initialise_sysdef(&sys, species, site_info, qpf);
      break;
    case 'r':
      if( (Fp = fopen(filename,"rb")) == NULL)
	 error("Couldn't open restart file \"%s\" for reading -\n%s\n", 
	       filename, strerror(errno));
      re_re_header(Fp, &restart_header, &control);
      control.rdf_interval = 0;       /* Don't attempt to read RDF data */
      re_re_sysdef(Fp, restart_header.vsn, &sys, &species, &site_info, &potpar);
      break;
    default:
      error("Internal error - invalid input type", "");
   }
   allocate_dynamics(&sys, species);

   if( data_source == 0 )		/* If called interactively	      */
   {
      fputs( "Where is the configurational information kept?\n", stderr);
      if( intyp == 's' )
      {
	 fputs( "In a lattice start file(1) or a dump dataset(2)?\n", stderr);
	 if( (ans_i = get_int("? ", 1, 2)) == EOF)
	    exit(2);
	 data_source = ans_i-1 ? 'd' : 's';
      }
      else if( intyp == 'r' )
      {
	 fputs( "In a restart file(1) or a dump dataset(2)?\n", stderr);
	 if( (ans_i = get_int("? ", 1, 2)) == EOF)
	    exit(2);
	 data_source = ans_i-1 ? 'd' : 'r';
      }
   }

   switch(data_source)			/* To read configurational data	      */
   {
    case 's':				/* Lattice_start file		      */
	lattice_start(Fp, &sys, species, qpf);
	moldy_out(1, 1, 1, &sys, species, site_info, atom_sel, outsw, insert);
      break;
    case 'r':				/* Restart file			      */
	init_averages(sys.nspecies, restart_header.vsn,
		      control.roll_interval, control.roll_interval,
		      &av_convert);
	read_restart(Fp, restart_header.vsn, &sys, av_convert);
	moldy_out(1, 1, 1, &sys, species, site_info, atom_sel, outsw, insert);
      break;
    case 'd':				/* Dump dataset			      */
	if( dump_name == 0 )
	{
	   fputs("Enter canonical name of dump files (as in control)\n",stderr);
	   if( (dump_name = get_str("Dumps? ")) == NULL)
	      exit(2);
	}

	/*
	 *  Ensure that the dump limits start, finish, inc are set up,
	 *  either on command line or by user interaction.
	 */
	do
	{
	   rflag = 0;
	   if( dumplims == NULL )
	   {
	      fputs("Please specify range of dump records in form", stderr);
	      fputs(" in form start-finish:increment\n", stderr);
	      dumplims = get_str("s-f:n? ");
	   }
	   if( forstr(dumplims, &start, &finish, &inc) )
	   {
	      rflag++;
	      fputs("Invalid range for dump records \"", stderr);
	      fputs(dumplims, stderr);
	      fputs("\"\n", stderr);
	   }
	   if( start > finish || start < 0 || inc <= 0 )
	   {
	      rflag++;
	      fputs("Dump record limits must satisfy", stderr);
	      fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
	   }
	   if( rflag)
	   {
	      (void)free(dumplims);
	      dumplims = NULL;
	   }
	} while(rflag);   
		
	/*
	 * Allocate buffer for data
         */
	dump_size = DUMP_SIZE(~0)*sizeof(float);
	if( (dump_buf = (float*)malloc(dump_size)) == 0)
	   error("malloc failed to allocate dump record buffer (%d bytes)",
		 dump_size);
	/*
	 * Loop over dump records, ascertaining which file they are in
	 * and opening it if necessary.  Call output routine.
	 */
#if defined (HAVE_POPEN) 
	sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %s %s",
		sys.nmols,sys.nmols_r, dumplims, dump_name);
	if( (Dp = popen(dumpcommand,"r")) == 0)
	   error("Failed to execute \'dumpext\" command - \n%s",
		 strerror(errno));
#else
	tempname = tmpnam((char*)0);
	sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %s -o %s %s",
		sys.nmols,sys.nmols_r, dumplims, tempname, dump_name);
	system(dumpcommand);
	if( (Dp = fopen(tempname,"rb")) == 0)
	   error("Failed to open \"%s\"",tempname);
#endif
	for(irec = start; irec <= finish; irec+=inc)
	{
	   
	   if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
              error("Error reading record %d in dump file - \n%s\n",
		    irec, strerror(errno));

	   dump_to_moldy(dump_buf, &sys);

	   moldy_out(iout++, irec, inc, &sys, species, site_info, atom_sel, outsw, insert);
#ifdef DEBUG
	   fprintf(stderr,"Sucessfully read dump record %d from file  \"%s\"\n",
		   irec%header.maxdumps, dump_name);
#endif
	}
#if defined (HAVE_POPEN) 
	pclose(Dp);
#else
	fclose(Dp);
	remove(tempname);
#endif
	break;
      default:
	break;
     }
   return 0;    
}
      

		   
			     
