/* MSD - Program for calculating mean square displacements
Copyright (C) 1996 Craig Fisher
For use with dump files from MOLecular DYnamics simulation code, Moldy,
by Keith Refson
 
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
what you give them.   Help stamp out software-hoarding! */
/***************************************************************************
***********  *
 * msd		code for calculating mean square displacements of centres of mass of  * 
 *		molecules from MolDy dump files					      *
 *		nb. msd intervals taken relative to extract dump slices		      *
 ***************************************************************************
*********** 
 *  Revision Log
 *  $Log: msd.c,v $
 *  Revision 1.1  1996/10/24 16:50:12  craig 
 *  Initial revision
 *
 */
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
/******************************************************************************
 * Dummies of 'moldy' routines so that msd may be linked with moldy library*
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
 * minimum.	return the lower of two integer values			      *
 ******************************************************************************/
int
minimum(a,b)
int	a,b;
{
    return (a < b) ? a: b;
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
   {
      fputs("Limits must satisfy", stderr);
      fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
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
}
/******************************************************************************
 * msd_calc().  Calculate mean square displacements of each species for one   * 
 * time slice                           		        	      *	
 ******************************************************************************/
void
msd_calc(system, species, init_spec, msd, imsd)
system_mt	*system;
spec_mt		species[];
spec_mt		init_spec[];
real		**msd;
int		imsd;
{
   spec_mt	*spec;
   int		i, imol, ispec=0;
   real		*sum = dalloc(system->nspecies);

   zero_real(sum,system->nspecies);

   for(spec = species; spec < species+system->nspecies; ispec++, init_spec++,
    				spec++)
   {   
      for(imol=0;imol<spec->nmols;imol++)
      {             
        /* Summate squared differences */
         for( i=0; i<3; i++)
            sum[ispec] +=
SQR(spec->c_of_m[imol][i]-init_spec->c_of_m[imol][i]);           
      } 
      msd[imsd][ispec] += sum[ispec]/spec->nmols;
   }

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));      
}
/******************************************************************************
 * traj_con().  Connect molecular c_of_m`s into continuous trajectories       * 
 ******************************************************************************/
void
traj_con(system, species, prev_slice)
system_mt	*system;
spec_mt		species[];
spec_mt		prev_slice[];

{
   spec_mt	*spec;
   double	box[3];
   mat_mt	hinv;
   mat_mp	h = system->h;
   int		i, imol;
 
   invert(h, hinv);
   box[0] = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   box[1] = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   box[2] = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));

   for(spec = species; spec < species+system->nspecies; prev_slice++, spec++)
   {                                   
      for( imol=0; imol<spec->nmols; imol++)
      {
        for (i = 0; i < 3; i++)
           spec->c_of_m[imol][i] = spec->c_of_m[imol][i] - box[i]*floor(
              (spec->c_of_m[imol][i]-prev_slice->c_of_m[imol][i])/box[i]+0.5);
      }
   }
}
/******************************************************************************
 * alloc_init().  Create arrays for each species	        	      * 
 ******************************************************************************/
void
alloc_init(nspecies, spec_pp)
int		nspecies;
spec_mp		*spec_pp;
{
   *spec_pp    = aalloc(nspecies, spec_mt);
}   
  
/******************************************************************************
 * msd_out().  Output routine for displaying msd results		      * 
 ******************************************************************************/
void
msd_out(nspecies, species, msd, nav, imsd)
int		nspecies;
spec_mt		*species;
real		**msd;
int		nav[];
int		imsd;
{
   int		ispec, i;
   spec_mp	spec;

   for( spec = species, ispec = 0; ispec < nspecies; spec++, ispec++ )
   {
       fputs(spec->name, stderr);
       for( i = 0; i <= imsd; i++)
       {
         if( nav[i] != 0)
       	    msd[i][ispec]/= nav[i];
       	 (void)printf("%9.7f\n", msd[i][ispec]);
       }
   }
}
/******************************************************************************
 * copy_cofm().  Duplicate values of c_of_m`s to another array                *
 ******************************************************************************/
void
copy_cofm(system, species, dupl_spec)
system_mt	*system;
spec_mt		species[];
spec_mt		dupl_spec[];
{
   spec_mt	*spec;
   int		imol;

   for(spec = species; spec < species+system->nspecies; dupl_spec++,spec++)
   {
        for(imol=0; imol < spec->nmols; imol++)
        {
        dupl_spec->c_of_m[imol][0] = spec->c_of_m[imol][0];
        dupl_spec->c_of_m[imol][1] = spec->c_of_m[imol][1];
        dupl_spec->c_of_m[imol][2] = spec->c_of_m[imol][2];
        }
   }
}
/******************************************************************************
 * init_species().  Create arrays of c_of_m`s for each molecule of species    *
 ******************************************************************************/
void
init_species(system, species, init_spec)
system_mt	*system;
spec_mt		species[];
spec_mt		init_spec[];
{
   spec_mt	*spec;

   for(spec = species; spec < species+system->nspecies; init_spec++,spec++)
         init_spec->c_of_m = ralloc(spec->nmols);
}
/******************************************************************************
 * main().   Driver program for calculating MSDs from MOLDY dump files        *
 * Acceptable inputs are sys-spec files or restart files. Actual 	      *
 * configurational info must be read from dump files.			      *
 * Call: msd [-s sys-spec-file] [-r restart-file]. 			      *
 * If neither specified on command line, user is interrogated.		      *
 ******************************************************************************/
contr_mt		control;

int
main(argc, argv)
int	argc;
char	*argv[];
{
   int	c, cflg = 0, ans_i, sym = 0;
   char 	line[80];
   extern char	*optarg;
   int		errflg = 0;
   int		intyp = 0;
   int		start, finish, inc;
   int		mstart, mfinish, minc;
   int		rflag, mflag;
   int		irec, it;
   char		*filename = NULL, *dump_name = NULL;
   char		*dumplims = NULL, *atom_sel = NULL;
   char		*insert = NULL, *msdlims = NULL;
   char		*tempname;
   char		dumpcommand[256];
   int		dump_size;
   float	*dump_buf;
   FILE		*Fp, *Dp;
   restrt_mt	restart_header;
   system_mt	sys;
   spec_mt	*species;
   spec_mt	*init_spec, *prev_slice;
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf;
   contr_mt	control_junk;
   int		*nav;
   real		**msd;

#define MAXTRY 100
   control.page_length=1000000;

   while( (c = getopt(argc, argv, "r:s:d:t:x:o:") ) != EOF )
      switch(c)
      {
       case 'r':
       case 's':
	 if( intyp )
	    errflg++;
	 intyp = c;
	 filename = optarg;
	 break;
       case 'd':
	 dump_name = optarg;
	 break;
       case 't':
	 dumplims = optarg;
         break;
       case 'x':
	 msdlims = optarg;
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error("failed to open file \"%s\" for output", optarg);
	 break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fputs("Usage: msd [-s sys-spec-file] [-r restart-file] ",stderr);
      fputs("[-d dump-files] [-t s[-f[:n]]] [-x s[-f[:n]]] [-o output-file]\n",
            stderr);
      exit(2);
   }

   if(intyp == 0)
   {
      fputs("How do you want to specify the simulated system?\n", stderr);
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
      re_re_header(Fp, &restart_header, &control_junk);
      re_re_sysdef(Fp, restart_header.vsn, &sys, &species, &site_info, &potpar);
      break;
    default:
      error("Internal error - invalid input type", "");
   }
   allocate_dynamics(&sys, species);

  /* Dump dataset			      */
   if( dump_name == 0 )
   {
	fputs("Enter canonical name of dump files (as in control)\n",stderr);
	if( (dump_name = get_str("Dump file name? ")) == NULL)
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
          fputs(" start-finish:increment\n", stderr);
          dumplims = get_str("s-f:n? ");
       }
       if( forstr(dumplims, &start, &finish, &inc) )
       {
          rflag++;
          fputs("Invalid range for dump records \"", stderr);
          fputs(dumplims, stderr);
          fputs("\"\n", stderr);
       }
       if( rflag)
       {
          (void)free(dumplims);
          dumplims = NULL;
       } 
   } while(rflag);

   if( msdlims != NULL)
   {
      do
      {
         mflag = 0;
         if( forstr(msdlims, &mstart, &mfinish, &minc) )
         {  
	   mflag++;
           fputs("Invalid range for msd intervals \"", stderr);
	   fputs(msdlims, stderr);
	   fputs("\"\n", stderr);
         }
         if( mstart*inc > finish-start || mfinish*inc > finish-start)
         {
            fputs("MSD interval exceeds dump range\n",stderr);
            mflag++;
         }
         if( mflag )
         {
            (void)free(msdlims);
            msdlims = NULL;
            fputs("Please specify msd intervals in form", stderr);
            fputs(" start-finish:increment\n", stderr);
            msdlims = get_str("s-f:n? ");
         }
       } while(mflag);         
    }
    else
    {
      /* Use default values for msd interval limits */
       mstart = 1;
       mfinish = finish-start;
       minc = inc;
     } 
      
  /*
   * Allocate buffer for data
   */
     dump_size = DUMP_SIZE(~0)*sizeof(float);
   
   /* create array of msds ie. one for each species and time interval */
   msd = arralloc(sizeof(real),2,0,(mfinish-mstart)/minc,0,sys.nspecies-1); 
   nav = ialloc((mfinish-mstart)/minc);  /* number of calculations of each
msd */

   for( it=0; it<=mfinish-mstart; nav[it/minc]=0,it+=minc) /* initialize
array of nav to 0 */

  /* create arrays for initial c_of_m`s and previous c_of_m`s for each
species */
   alloc_init(sys.nspecies, &init_spec);  
   alloc_init(sys.nspecies, &prev_slice);

 /* Outer loop for setting starting time slice to calculate
                  displacement difference */ 
   for(it = start; it <= finish-mstart*inc; it+=inc)
   {
      if( (dump_buf = (float*)malloc(dump_size)) == 0)
         error("malloc failed to allocate dump record buffer (%d bytes)",
       dump_size);

     /*
      * Loop over dump records, ascertaining which file they are in
      * and opening it if necessary.  Calculate msd`s and call output routine.
      */
#if defined (unix) || defined (__unix__)
      sprintf(dumpcommand,"dumpext -R %d -Q %d -b -c 0 -t %d-%d:%d %s",
	sys.nmols, sys.nmols_r, it, minimum(it+mfinish*inc,finish), inc, dump_name);
   
      if( (Dp = popen(dumpcommand,"r")) == 0)
        error("Failed to execute \'dumpext\" command - \n%s",
            strerror(errno));
#else
      tempname = tmpnam((char*)0);
      sprintf(dumpcommand,"dumpext -R %d -Q %d -b -c 0 -t %d-%d:%d -o %s %s",
	   sys.nmols,sys.nmols_r, it, minimum(it+mfinish*inc,finish), inc, tempname, 
	        dump_name);
      system(dumpcommand);
      if( (Dp = fopen(tempname,"rb")) == 0)
         error("Failed to open \"%s\"",tempname);
#endif 
     /* Inner loop for calculating displacement at time slice irec from
that at it */
      for( irec = it; irec <= minimum(it+mfinish*inc,finish); irec+=inc)
      {
         if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
           error("Error reading record %d in dump file - \n%s\n",
              irec, strerror(errno));
         dump_to_moldy(dump_buf, &sys);  /*read dump data */

         if( irec == start)
         {
            init_species(&sys, species, init_spec); /* create arrays for
each c_of_m */
            init_species(&sys, species, prev_slice);
         }
         if( irec == it)
            copy_cofm(&sys, species, init_spec);
	 else
	    traj_con(&sys, species, prev_slice);
	        
        /* Only perform calculation if msd interval is in range mstart to
mfinish */ 
         if( irec-it >= mstart*inc && irec != it)
         {
            if( fmod((real)(irec-it-mstart*inc),(real)(minc*inc)) == 0.0)
            {
               msd_calc(&sys, species, init_spec, msd,
(irec-it-mstart*inc)/(minc*inc));
               nav[(irec-it-mstart*inc)/(minc*inc)]++;
            } 
         } 
        /* Make copy of c_of_m`s for joining trajectory of next time slice */
            copy_cofm(&sys, species, prev_slice); 
#ifdef DEBUG
      fprintf(stderr,"Sucessfully read dump record %d from file  \"%s\"\n",
	   irec%header.maxdumps, dump_name);
#endif
      }
      xfree(dump_buf);
#if defined (unix) || defined (__unix__)
   pclose(Dp);
#else
   fclose(Dp);
   remove(tempname);
#endif
   }
   /* Display species, msd and number averaged over */
       msd_out(sys.nspecies, species, msd, nav, (mfinish-mstart)/minc);

   return 0;    
}
      

		   
			     

