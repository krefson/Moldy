/* MSD - Mean Square Displacement code for use with MOLDY by Keith Refson
Copyright (C) 1995 Craig Fisher
 
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
/*======================== Structure vars ====================================*/
typedef struct
{
int	start;
int	finish;
int	inc;
} limits;
/******************************************************************************
 * Dummies of 'moldy' routines so that msd may be linked with moldy library   *
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

   (void)fprintf(stderr, "msd: ");
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

   (void)fprintf(stderr, "msd: ");
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
 * get_sym().  Read a character from stdin and match to supplied set	      *
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
 * get_str().  Read a string from stdin, issuing a prompt.		      *
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
   if( *start > *finish || *start < 0 || *inc <= 0 )
      goto limerr;
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
 * msd_calc. Routine for calculating mean square displacements of selected    *
 *   	     species over given timeslices.				      *
 ******************************************************************************/
void
msd_calc(system, species, Tp, nav, iavd, mlims, slims, totmols)
system_mt       *system;
spec_mt		species[];
FILE		*Tp;
int 		nav, iavd;
limits		*mlims, *slims;
int		totmols;
{
int 		x, i, iav, it;
int		nspec, nmol;
double		**msd = (double**)arralloc(sizeof(double),2,0,mlims->finish-
			mlims->start,0,slims->finish-slims->start);
spec_mt		*spec;
vec_mp		tr, prev_tr;

  tr = ralloc(totmols);
  prev_tr = ralloc(totmols);

  zero_real(*msd, (slims->finish-slims->start+1)*(mlims->finish-mlims->start+1));
  zero_real(tr, 3*totmols);
  zero_real(prev_tr, 3*totmols);

  for(iav=0; iav<nav; iav+=iavd)
  {
     fseek(Tp,(long)(iav*totmols*sizeof(vec_mt)),SEEK_SET);
     if( fread(prev_tr,sizeof(vec_mt),totmols,Tp) < totmols || ferror(Tp))
 	 error("Error reading t(init) from temp trajectory file - \n",strerror(errno));

     for( it = 0; it < mlims->finish-mlims->start+1; it+=mlims->inc)
     {
        fseek(Tp,(long)((it+iav+1)*totmols*sizeof(vec_mt)),SEEK_SET);

        if( fread(tr,3*sizeof(double),totmols,Tp) < totmols || ferror(Tp))
 	  error("Error reading t(final) from temp trajectory file - \n",strerror(errno));
 	nspec = 0;
 	nmol=0;
        for(spec=species+slims->start-1; spec<species+slims->finish;
        		 spec+=slims->inc)
        {
           for(i=0; i<spec->nmols; i++)
           {
             for (x=0; x<3; x++)
	       msd[it][nspec]=msd[it][nspec]+SQR(tr[x][nmol]-prev_tr[x][nmol]);
	     nmol++; 
	   }
	   nspec++;
	}
     }
  }
  for( it = 0; it < mlims->finish - mlims->start + 1; it+=mlims->inc)
  {
    nspec = 0;
    for (spec=species+slims->start-1; spec<species+slims->finish;
    		 spec+=slims->inc)
    { 
       msd[it][nspec] = msd[it][nspec]*iavd/(nav*species->nmols);
       printf("%12.10f  ", msd[it][nspec]);
       nspec++;
    }
    printf("\n");
  }
  if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
  (void)free(tr, prev_tr);
  (void)xfree(msd);
}
/******************************************************************************
 * main().   Driver program for calculating mean square displacements from    *
 * 	     dump files.  						      *
 * Call: msd [-r restart-file] [-s sys-spec-file] dump-file  	              *
 * If dumpfile not specified on command line, user is interrogated.	      *
 ******************************************************************************/
contr_mt	control;

int
main(argc, argv)
int	argc;
char	*argv[];
{
   int	c, cflg =0, ans_i, sym, i, x, iav, it, nmol, nspec;
   char 	line[80];
   extern char	*optarg;
   int		errflg = 0, intyp = 0;
   int		dflag, irec, nav;
   long int	frac = 0, iavd = 0, totmols;
   char		*filename = NULL, *dump_name = NULL;
   char		*msdlims = NULL, *speclims = NULL, *dumplims = NULL;
   char		*tempname, *trname;
   char		dumpcommand[256];
   int		dump_size, maxslice;
   float	*dump_buf;
   FILE		*Fp, *Dp, *Tp, *dump_file;
   restrt_mt	restart_header;
   pot_mt	*potpar;
   quat_mp	*qpf;
   system_mt	sys;
   spec_mt	*species;
   contr_mt	control_junk;
   site_mt	*site_info;
   dump_mt      header;
   double	***ttr, **msd;
   spec_mt	*spec;
   limits	mlims, slims, dlims;
	 
#define MAXTRY 100

   while( (c = getopt(argc, argv, "r:s:o:m:z:d:i:f:") ) != EOF )
     switch(c)
     {
     case 'r':
     case 's':
	 if( intyp )
	    errflg++;
	 intyp = c;
	 filename = optarg;
	 break;
     case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error("failed to open file \"%s\" for output", optarg);
	 break;
     case 'd':
         dumplims = optarg;
         break;
     case 'z':
         speclims = optarg;
         break;
     case 'm':
         msdlims = optarg;
         break;
     case 'f':
         frac = strtol(optarg,(char**)0,0);
         break;
     case 'i':
         iavd = strtol(optarg,(char**)0,0);
         break;
     default:
     case 'h':
     case '?':
	 errflg++;
     }

   if( errflg )
   {
      fputs("Usage: msd [-s sys-spec-file] [-r restart-file] [-m s[-f[:n]]]\n",stderr);
      fputs("[-z s[-f[:n]]] [-d s[-f[:n]]] [-o output-file] dump-file\\s\n", stderr);
      exit(2);
   }

   if(intyp == 0)
   {
      fputs("Is the system identification information in:\n",stderr);
      fputs("(1) a System Specification file or (2) a Restart file", stderr);
      if( (ans_i = get_int("? ", 1, 2)) == EOF )
	 exit(2);
      intyp = ans_i-1 ? 'r': 's';
      if( intyp == 's' )
      {
	 fputs( "Do you need to skip 'control' information ", stderr);
	 if( (sym = get_sym("(Y/N)? ","yYnN")) == 'y' || sym == 'Y')
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
      re_re_header(Fp, &restart_header, &control_junk);
      re_re_sysdef(Fp, &sys, &species, &site_info, &potpar);
      break;
      
    default:
      error("Internal error - invalid input type", "");
   }
  
   allocate_dynamics(&sys, species);

   /* Dump dataset */
   dump_name = argv[optind];
   if( dump_name == 0 )
   {
      fputs("Enter canonical name of dump files (as in control)\n",stderr);
      if( (dump_name = get_str("Dump files? ")) == NULL)
	 exit(2);
   }
/*
 *  Ensure that the species limits start, finish, inc are set up,
 *  either on command line or with default values.
 */
   if( speclims == NULL )
   {
      slims.start = 1;
      slims.finish = sys.nspecies;
      slims.inc = 1;
   }
   else
   {
      if( forstr(speclims, &(slims.start), &(slims.finish), &(slims.inc)) ||
      		 slims.start == 0)
 	 error("Invalid range for species: \"%s\"\n",speclims);
      if( slims.finish > sys.nspecies )
            error("There are only %d species in this system\n", sys.nspecies);
   } 
/*
 * Allocate buffer for data
 */
/*   dump_size = DUMP_SIZE(~0)*sizeof(float);
   if( (dump_buf = (float*)malloc(dump_size)) == 0)
      error("Failed to allocate dump record buffer (%d bytes)",dump_size);*/
		 
/* Calculate total number of timeslices */
/*   if( (dump_file = fopen(dump_name, "rb")) == 0 )
      error("Failed to open dump file \"%s\"\n",dump_name);

   fseek(dump_file, 0L, SEEK_SET);
   if( fread((char*)&header, sizeof(dump_mt), 1, dump_file) == 0 ||
    ferror(dump_file) ) 
      error("Error reading header of dump file %s\n",dump_name); 
   maxslice = (header.istep-1)/header.dump_interval + header.ndumps - 1;
   fclose(dump_file); */
/*
 *  Ensure that the dump record limits start, finish, inc are set up,
 *  either on command line or with default values.
 */
   if( dumplims == NULL )
   {
      maxslice = 100;
      dlims.start = 0;
      dlims.finish = maxslice;
      dlims.inc = 1;
      dumplims = "0-99";
   }
   else
   {
     if( forstr(dumplims, &(dlims.start), &(dlims.finish), &(dlims.inc)) )
       error("Invalid range for dump records: \"%s\"\n",dumplims);
     if( dlims.finish > maxslice )
       error("The last record in this dump file is %d\n", maxslice);
   }
/*
 * Loop over dump records, ascertaining which file they are in
 * and opening it if necessary.
 */
/*#if defined (unix) || defined (__unix__)
   sprintf(dumpcommand,"dumpext -R %d -Q %d -b -c 0 -t %s %s", sys.nmols,
      sys.nmols_r, dumplims, dump_name);
   if( (Dp = popen(dumpcommand,"r")) == 0 )
      error("Failed to execute \'dumpext\" command - \n%s", strerror(errno));
#else
   tempname = tmpnam((char*)0);
   sprintf(dumpcommand,"dumpext -R %d -Q %d -b -c 0 -t %s -o %s %s",
	sys.nmols, sys.nmols_r, dumplims, tempname, dumplims, dump_name);
   system(dumpcommand);
   if( (Dp = fopen(tempname,"rb")) == 0 )
	error("Failed to open \"%s\"",tempname);
#endif */	   
/*
 *  Ensure that the msd limits start, finish, inc are set up,
 *  either on command line or by user interaction.
 */
   do
   {
      dflag = 0;
      if( msdlims == NULL )
      {
         fputs("Please specify RANGE for MSD calculations ", stderr);
         fputs("in form start-finish:increment\n", stderr);
         msdlims = get_str("s-f:n? ");
      }
      if( forstr(msdlims, &(mlims.start), &(mlims.finish), &(mlims.inc)) ||
            mlims.start == 0 )
      {
         dflag++;
	 fputs("Invalid range for MSD calculations \"", stderr);
	 fputs(msdlims, stderr);
	 fputs("\"\n", stderr);
	 (void)free(msdlims);
	 msdlims = NULL;
      }
   } while(dflag);
           
/* Default values for frac and iavd */    
   if ( frac == 0 )
        frac = dlims.finish - dlims.start + 1;
   if ( iavd == 0 )
        iavd = 1;         
/* Calculate number of timeslices to average msd over */          
   nav = (dlims.finish - dlims.start + 1)/frac;
   if( mlims.finish > (dlims.finish - dlims.start - nav + 1))
	mlims.finish = dlims.finish - dlims.start - nav + 1;

/* Read centres of mass of all molecules in system */
   ttr = arralloc(sizeof(double),3,0,mlims.finish+nav-1,0,sys.nmols-1,0,2);
   
   if( (Tp = fopen("trajj","r")) == NULL)
    	fputs("Failed to open trajj",stderr);
/*
   for (irec=0;irec<mlims.finish+nav;irec++)
   {
      if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
         error("Error reading record %d in dump file - \n%s\n", irec,
            strerror(errno));
      dump_to_moldy(dump_buf, &sys);

      for ( nmol=0,spec=species+slims.start-1; spec<species+slims.finish;
      		 spec+=slims.inc )
      {
       Join trajectories 
	for( i=0; i<spec->nmols; nmol++,i++ )
  	  for(x=0; x<3; x++)
   	  {
    	    if( irec > 0 )
   	       spec->c_of_m[i][x] = spec->c_of_m[i][x] - 
   	           floor(spec->c_of_m[i][x] - ttr[irec-1][nmol][x]+0.5);
   	    ttr[irec][nmol][x] = spec->c_of_m[i][x];
   	  }
   	 Convert to cartesian coordinates 
   	mat_vec_mul(sys.h, spec->c_of_m, spec->c_of_m, spec->nmols);
      } 
      mat_vec_mul(sys.h, ttr[irec], ttr[irec], sys.nmols);                 
   	
   } */

totmols = 2304;
for(irec=0;irec<mlims.finish+nav;irec++)
  for(i=0;i<totmols;i++) 
    if(fscanf(Tp,"     %8f     %8f     %8f",ttr[irec][i][0],ttr[irec][i][1],ttr[irec][i][2])==0)
         error("Error reading record %d in dump file - \n%s\n", irec,
            strerror(errno));
/* MSD calculation */
 /*  msd_calc(&sys, species, ttr, nav, iavd, &mlims, &slims, nmol); */
fclose(Tp);
msd = arralloc(sizeof(double),2,0,mlims.finish-
			mlims.start,0,slims.finish-slims.start);

  zero_real(*msd, (slims.finish-slims.start+1)*(mlims.finish-mlims.start+1));

  for(iav=0; iav<nav; iav+=iavd)
  {
     for( it = 0; it < mlims.finish-mlims.start+1; it+=mlims.inc)
     {
 	nspec = 0;
 	nmol=0;
        for(spec=species+slims.start-1; spec<species+slims.finish;
        		 spec+=slims.inc)
        {
           for(i=0; i<spec->nmols; i++)
           {
             for (x=0; x<3; x++)
	       msd[it][nspec]=msd[it][nspec]+SQR(ttr[it+iav+1][nmol][x]-ttr[iav][nmol][x]);
	     nmol++; 
	   }
	   nspec++;
	}
     }
  }
  for( it = 0; it < mlims.finish - mlims.start + 1; it+=mlims.inc)
  {
    nspec = 0;
    for (spec=species+slims.start-1; spec<species+slims.finish;
    		 spec+=slims.inc)
    { 
       msd[it][nspec] = msd[it][nspec]*iavd/(nav*spec->nmols);
       printf("%12.10f  ", msd[it][nspec]);
       nspec++;
    }
    printf("\n");
  }
  if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
  (void)xfree(msd);
  
/*#if defined (unix) || defined (__unix__)
   pclose(Dp);
#else
   fclose(Dp);
   remove(tempname);
#endif*/
   
   return 0;    
}
      

		   
			     
