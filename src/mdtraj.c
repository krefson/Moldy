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
what you give them.   Help stamp out software-hoarding! */
/**************************************************************************************
 * mdtraj    	code for calculating and printing trajectories of 		      *       
 *              centres of mass of molecules from MolDy dump files		      *
 *		nb. output in form columns: x1 y1 z1 against timeslice	 	      *
 ************************************************************************************** 
 *  Revision Log
 *  $Log: mdtraj.c,v $
 *  Revision 1.3  1997/07/10 18:05:26  craig
 *  Previous slice spec_mt converted to cofm array
 *  Traj_cofm zeroed after call to arralloc
 *  Removed error in upper limit of species range
 *
 *  Revision 1.2  1997/07/04 14:55:21  craig
 *  Trajectory connection now performed in fractional coordinates
 *  Removed alloc_init function
 *
 *  Revision 1.1  1997/04/10 19:01:12  craig 
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
 * Dummies of 'moldy' routines so that mdtraj may be linked with moldy library*
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

   (void)fprintf(stderr, "mdtraj: ");
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

   (void)fprintf(stderr, "mdtraj: ");
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
 * get_real().  Read an int from stdin, issuing a prompt and checking         *
 * validity and range.  Loop until satisfied, returning EOF if appropriate.   *
 ******************************************************************************/
int get_real(prompt, lo, hi)
char	*prompt;
real 	lo, hi;
{
   char		ans_str[80];
   int	ans_r; 
   int		ans_flag;
   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, "%d", &ans_r) == 1 && ans_r >= lo && ans_r <= hi)
	 ans_flag++;
   }
   if( ans_flag )
      return(ans_r);
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
 * forstr.  Parse string str of format s-f:n (final 2 parts optional),        *
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
   if( *start+*inc > *finish+1 || *start < 0 || *inc <= 0 )
   {
      fputs("Limits must satisfy", stderr);
      fputs(" finish >= start + increment, start >= 0 and increment > 0\n", stderr);
      goto limerr;
   }
   return 0;
 limerr:
   return -1;
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
   double	box[3];
   mat_mp	h = system->h;
   int		i;

   for( i=0; i<3; i++)
   {
       box[i] = sqrt(SQR(h[0][i]) + SQR(h[1][i]) + SQR(h[2][i]));

       if( range_flag[i])
          switch(i)
          {
           case 0:
             range[i][0] = get_real("Enter x minimum: ", -1*box[i],box[i]);
             range[i][1] = get_real("Enter x maximum: ", range[i][0],2*box[i]);
             break;
           case 1:
             range[i][0] = get_real("Enter y minimum: ", -1*box[i],box[i]);
             range[i][1] = get_real("Enter y maximum: ", range[i][0],2*box[i]);
             break;
           case 2:
             range[i][0] = get_real("Enter z minimum: ", -1*box[i],box[i]);
             range[i][1] = get_real("Enter z maximum: ", range[i][0],2*box[i]);
             break;
          }    
       else
       {
          range[i][0] = -1*box[i];
          range[i][1] = box[i];
       }
   }
   return 0;
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
   invert(system->h,hinv);
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
      a0 = vec[i][0];  a1 = vec[i][1];  a2 = vec[i][2];

      vec[i][0] = m[0][0]*a0 + m[0][1]*a1 + m[0][2]*a2;
      vec[i][1] = m[1][0]*a0 + m[1][1]*a1 + m[1][2]*a2;
      vec[i][2] = m[2][0]*a0 + m[2][1]*a1 + m[2][2]*a2;
   }
}
/******************************************************************************
 * traj_con().  Connect molecular c_of_m`s into continuous trajectories       * 
 ******************************************************************************/
void
traj_con(system, species, prev_cofm, traj_cofm, sp_range)
system_mt	*system;
spec_mt		species[];
real		**prev_cofm,
		**traj_cofm;
int		sp_range[3];
{
   spec_mt	*spec;
   int		i, imol, ispec, totmol=0;
 
   for( ispec = sp_range[0], spec=species+sp_range[0]; spec <= species+sp_range[1]; 
                                              ispec+=sp_range[2], spec+=sp_range[2])
     for( imol=0; imol<spec->nmols; totmol++, imol++)
	if( prev_cofm == 0 ) 
	{
	   for( i = 0; i < 3; i++)
	      traj_cofm[totmol][i] = spec->c_of_m[imol][i];
	}
	else 
	{
	   for( i = 0; i < 3; i++)
	      traj_cofm[totmol][i] = spec->c_of_m[imol][i] - 
		 floor(spec->c_of_m[imol][i]-prev_cofm[totmol][i]+0.5);
	}
}
/******************************************************************************
 * traj_out().  Output routine for displaying trajectories                    * 
 ******************************************************************************/
void
traj_out(system, species, traj_cofm, nslices, range, sp_range)
system_mt	*system;
spec_mt		species[];
real		***traj_cofm;
real		range[3][2];
int		nslices, sp_range[3];
{
   int		totmol, imol, i, itime;
   spec_mp	spec;

   totmol = 0;
   for( spec = species+sp_range[0]; spec <= species+sp_range[1]; spec+=sp_range[2])
   {
     (void)printf("# %s\n",spec->name);
     for( imol = 0; imol < spec->nmols; totmol++, imol++)
     {
       if( traj_cofm[0][totmol][0] >= range[0][0] && traj_cofm[0][totmol][0] <= range[0][1] &&
          traj_cofm[0][totmol][1] >= range[1][0] && traj_cofm[0][totmol][1] <= range[1][1] &&
             traj_cofm[0][totmol][2] >= range[2][0] && traj_cofm[0][totmol][2] <= range[2][1])
          for( itime = 0; itime <= nslices; itime++)
          {
             if( totmol == 0)
                 mat_vec_mul3(system->h, traj_cofm[itime], system->nmols); 
             for(i= 0; i < 3; i++)
                 (void)printf("%f ",traj_cofm[itime][totmol][i]);
             (void)printf("\n");
           }
           (void)printf("\n\n");
     }
   }
}
/******************************************************************************
 * create_array(). Link pointers and base array to make larger pseudo-array   *
 ******************************************************************************/
void
create_array(system, species, sp_cofm, base_cofm)
system_mt	*system;
spec_mt		species[];
vec_mp		*sp_cofm, base_cofm;
{
   spec_mt	*spec;
   int		ispec, imol=0;

   for( ispec=0, spec=species; ispec < system->nspecies; ispec++, spec++)
   {  
      sp_cofm[ispec] = base_cofm + imol;
      imol+=spec->nmols;
   }
}
/******************************************************************************
 * main().   Driver program for calculating trajectories from MOLDY dumps     *
 * Acceptable inputs are sys-spec files or restart files. Actual 	      *
 * configurational info must be read from dump files.			      *
 * Call: mdtraj [-s sys-spec-file] [-r restart-file]. 			      *
 * If neither specified on command line, user is interrogated.		      *
 * Options [-x][-y][-z] prompt for inclusion limits in selected direction     *
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
   int		nslices, islice=0;
   int		sp_range[3];
   int		rflag, sflag;
   int		irec;
   char		*filename = NULL, *dump_name = NULL;
   char		*dumplims = NULL, *atom_sel = NULL;
   char		*insert = NULL, *speclims = NULL;
   char		*tempname;
   char		dumpcommand[256];
   int		dump_size;
   float	*dump_buf;
   FILE		*Fp, *Dp;
   restrt_mt	restart_header;
   system_mt	sys;
   spec_mt	*species;
 /*  spec_mt	*init_spec, *prev_slice; */
   real 	***traj_cofm;
   real		range[3][2];
   int		range_flag[3]={0,0,0};
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf;
   contr_mt	control_junk;
   int		i,j,it;

#define MAXTRY 100
   control.page_length=1000000;

   while( (c = getopt(argc, argv, "r:s:d:t:m:o:xyz") ) != EOF )
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
       case 'm':
	 speclims = optarg;
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error("failed to open file \"%s\" for output", optarg);
	 break;
       case 'x':
         range_flag[0] = 1;
         break;
       case 'y':
	 range_flag[1] = 1;
         break;
       case 'z':
         range_flag[2] = 1;
         break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fputs("Usage: mdtraj [-s sys-spec-file] [-r restart-file] ",stderr);
      fputs("[-d dump-files] [-t s[-f[:n]]] [-m s[-f[:n]]] ",stderr);
      fputs("[-x] [-y] [-z] [-o output-file]\n",stderr);
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

   if( speclims != NULL)
   {
      do
      {
         sflag = 0;
         if( forstr(speclims, &(sp_range[0]), &(sp_range[1]), &(sp_range[2])))
         {  
	   sflag++;
           fputs("Invalid range for molecule selection \"", stderr);
	   fputs(speclims, stderr);
	   fputs("\"\n", stderr);
         }
         if( sp_range[1] > sys.nspecies)
         {
            fputs("Molecule selection exceeds no. of species\n",stderr);
            sflag++;
         }
         if( sflag )
         {
            (void)free(speclims);
            speclims = NULL;
            fputs("Please specify molecule selection in form", stderr);
            fputs(" start-finish:increment\n", stderr);
            speclims = get_str("s-f:n? ");
         }
       } while(sflag);         
    }
    else
    {
      /* Use default values for molecule selection limits */
       sp_range[0] = 0;
       sp_range[1] = sys.nspecies-1;
       sp_range[2] = 1;
     } 
  /*
   * Allocate buffer for data
   */
     dump_size = DUMP_SIZE(~0)*sizeof(float);

   
  /* create arrays for previous c_of_m`s for each species */
   /*  prev_slice = aalloc(sys.nspecies, spec_mt); */
 
     nslices = (finish-start)/inc; /* no. of time slices */

  if( (dump_buf = (float*)malloc(dump_size)) == 0)
      error("malloc failed to allocate dump record buffer (%d bytes)",
          dump_size);
#if defined (HAS_POPEN) 
     sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d %s",
        sys.nmols, sys.nmols_r, start, finish, inc, dump_name);
   
     if( (Dp = popen(dumpcommand,"r")) == 0)
        error("Failed to execute \'dumpext\" command - \n%s",
            strerror(errno));
#else
     tempname = tmpnam((char*)0);
     sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d -o %s %s",
         sys.nmols,sys.nmols_r, start, finish, inc, tempname, dump_name);
     system(dumpcommand);
     if( (Dp = fopen(tempname,"rb")) == 0)
        error("Failed to open \"%s\"",tempname);
#endif
     /* Allocate memory for trajectory data and zero */
     traj_cofm = arralloc(sizeof(real),3,0,nslices,0,sys.nmols-1,0,2);
     zero_real(traj_cofm[0][0], (nslices+1)*sys.nmols*3);

/* Loop for calculating trajectories from current and previous time slices */ 
     for(irec = start; irec <= finish; irec+=inc)
     {
        if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
           error("Error reading record %d in dump file - \n%s\n",
              irec, strerror(errno));
        dump_to_moldy(dump_buf, &sys);  /*read dump data */

        if( irec == start )
	{
          range_in(&sys, range, range_flag);
          traj_con(&sys, species, 0, traj_cofm[islice], sp_range);
	}
	else
          traj_con(&sys, species, traj_cofm[islice-1], traj_cofm[islice], sp_range);

        islice++;

#ifdef DEBUG
        fprintf(stderr,"Sucessfully read dump record %d from file  \"%s\"\n",
	   irec%header.maxdumps, dump_name);
#endif
        xfree(dump_buf);
     }
     traj_out(&sys, species, traj_cofm, nslices, range, sp_range);
#if defined (HAS_POPEN) 
   pclose(Dp);
#else
   fclose(Dp);
   remove(tempname);
#endif
   return 0;    
}