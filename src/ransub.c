/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1999 Craig Fisher
Copyright (C) 1988, 1992, 1993, 1997 Keith Refson
 
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
 * ransub    	Code for randomly substituting species                                *
 *              in Moldy configuration files                                          *
 *		Randomly replaces species "m" with n molecules of species "u"	      *
 *		Output written in Moldy system specification format		      *
 ************************************************************************************** 
 *  Revision Log
 *  $Log: ransub.c,v $
 *  Revision 2.1  1999/06/24 16:05:44  craig
 *  Improved randomization of random number reseeder.
 *
 *  Revision 2.0  1999/06/03 15:39:34  craig
 *  Corrected memory freeing of dump limits.
 *  Tidied up use of structure variable 'dop'.
 *  Added loop to check if species being replaced exists.
 *
 *  Revision 1.4  1999/04/08 17:57:29  craig
 *  Options to specify dopant mass, charge and symbol added.
 *
 *  Revision 1.3  1999/03/23 15:09:45  craig
 *  Removed unnecessary variable 'is' from sys_spec_out.
 *
 *  Revision 1.2  1999/03/12 15:15:39  craig
 *  Altered energy conversion units to be consistent with defs.h values
 *
 *  Revision 1.1  1999/03/10 18:10:19  craig 
 *  Corrected bug limiting max no of substitutions to no of first species
 *  Upper limit to no of substituted species set to total no in system 
 *
 *  Revision 1.0  1999/03/05 17:41:12  craig 
 *  Initial revision
 *
 */
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
#ifdef HAVE_STDARG_H
gptr	*arralloc(size_mt,int,...); 	/* Array allocator */
#else
gptr	*arralloc();	        	/* Array allocator */
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
void    conv_potentials();
int	getopt();
gptr	*talloc();
char	*atime();
FILE	*popen();
/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
static  unit_mt prog_unit = {MUNIT, LUNIT, TUNIT, QUNIT};
static  char            *comm;
static  unit_mt input_unit = {MUNIT, LUNIT, TUNIT, _ELCHG};
contr_mt                control;

/* Time units for different energy units */
#define EV 1.018050697e-14; 		/* electron volts */
#define KJMOL 1.0e-13;			/* kilojoules per mole */
#define KCALS 4.88882131e-14;		/* kilocalories per mole */
#define E2A 2.682811715e-15;		/* electron charge squared per angstrom */
/*========================== External data references ========================*/

extern  const pots_mt   potspec[];          /* Potential type specification */

/******************************************************************************
 * Dummies of moldy routines so that ransub may be linked with moldy library  *
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
 * Structure declarations                                                     *
 ******************************************************************************/
typedef struct {
   int		*pos;       /* Array of substituting species' position nos */
   int		mols;       /* No of molecules of substituting species */
   double	mass;       /* Mass of substituting species */
   char         *charge;    /* Charge of substituting species */
   char		*name;      /* Name of substituting species */
   char		*sym;       /* Symbol of substituting species */
} dopant;
/******************************************************************************
 *  message.   Deliver error message to possibly exiting.  It can be called   *
 *             BEFORE output file is opened, in which case outt to stderr.    *
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

   (void)fprintf(stderr, "ransub: ");
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
 *  message.   Deliver error message to possibly exiting.                     *
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

   (void)fprintf(stderr, "ransub: ");
   (void)vfprintf(stderr, format, ap);
   fputc('\n',stderr);
   va_end(ap);

   exit(3);
}
static char * mystrdup(s)
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
 * sys_spec_out().  Write a system configuration to stdout in the form of a   *
 * system specification file for MOLDY                                        *
 ******************************************************************************/
void
sys_spec_out(system, species, molname, dop, site_info, potpar, intyp)
system_mt       *system;
spec_mt         species[];
char            *molname;
dopant          *dop;
site_mt         site_info[];
pot_mt          *potpar;
int             intyp;
{
   double       **site = (double**)arralloc(sizeof(double),2,
                                            0,2,0,system->nsites-1);
   spec_mt      *spec;
   double       a, b, c, alpha, beta, gamma;
   mat_mp       h = system->h;
   int          i, imol, isite, itot=1;
   int		specmol;
   int          pot_type, idi, idj, idij, ip;
   int          n_potpar = system->n_potpar;
   int		nunits;
   char         *specname;

   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

/* If c_of_m's read from sys-spec file, charges not converted */

   if( intyp != 's')
   {
      fputs("What energy units would you like the potentials in:\n",stderr);
      fputs("(1) eV, (2) kJ/mol, (3) kcal/mol, or (4) e**2/A",stderr);
      nunits = get_int(" ? ", 1, 4);
      switch(nunits)
      {
      case 1:
         input_unit.t = EV;
         break;
      case 2:
         input_unit.t = KJMOL;
         break;
      case 3:
         input_unit.t = KCALS;
         break;
      case 4:
         input_unit.t = E2A;
         break;
      }
      conv_potentials(&prog_unit, &input_unit, potpar, system->n_potpar,
           system->ptype, site_info, system->max_id);
   }

/* Write header for sys_spec file */
   (void)printf("# System specification file written by RANSUB on %s\n",atime());

/* Write site data for each molecule */
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      if( strcmp(strlower(spec->name), molname) )
             specmol = spec->nmols;
      else
             specmol = spec->nmols - dop->mols; /* Subtract number of substituting species */

      (void)printf("%s  %d  %s\n", spec->name, specmol,
                    spec->framework ? "framework" : "");
      for(isite=0; isite < spec->nsites; isite++)
         (void)printf("%d %9g %9g %9g %9g %9g %s\n",
                        spec->site_id[isite],
                        spec->p_f_sites[isite][0],
                        spec->p_f_sites[isite][1],
                        spec->p_f_sites[isite][2],
                        site_info[spec->site_id[isite]].mass,
                        site_info[spec->site_id[isite]].charge,
                        site_info[spec->site_id[isite]].name);

      if( !strcmp(strlower(spec->name), molname) && dop->name != NULL && dop->mols > 0 )
      {
         (void)printf("%s  %d  %s\n", dop->name, dop->mols,
                    spec->framework ? "framework" : "");
         (void)printf("%d %9g %9g %9g %9g %9g %s\n",
                        system->max_id,
                        spec->p_f_sites[0][0],
                        spec->p_f_sites[0][1],
                        spec->p_f_sites[0][2],
             dop->mass < 0 ? site_info[spec->site_id[0]].mass:dop->mass,
             dop->charge == NULL ? site_info[spec->site_id[0]].charge:atof(dop->charge),
             dop->sym == NULL ? site_info[spec->site_id[0]].name:dop->sym);
      }
   }
   (void)printf("end\n");

/* Write potential parameters for pairs of site_ids */

   (void)printf("%s\n",potspec[system->ptype].name);
   for(idi = 1; idi < system->max_id; idi++)
     for(idj = idi; idj < system->max_id; idj++)
     {
        idij = idj + idi*system->max_id;
        (void)printf("%5d %5d", idi, idj);
        for(ip = 0; ip < n_potpar; ip++)
           (void)printf("   %g",potpar[idij].p[ip]);
        (void)putchar('\n');
    }
    (void)printf("end\n");

/* Now we write the box dimensions */
   (void)printf("%g  %g  %g  %g  %g  %g  1  1  1\n",
          a,b,c,alpha,beta,gamma);

/* Followed by the molecules' centre of mass positions */
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      for(imol = 0; imol < spec->nmols; imol++)
      {
        specname = spec->name;
        if( !strcmp(strlower(spec->name), molname))
          for( i = 0; i < dop->mols; i++)
             if( dop->pos[i] == imol )
                specname = dop->name;
        (void)printf("%s ", specname);
        for( i = 0; i < 3; i++)
          (void)printf("%9g ",
             spec->c_of_m[imol][i]+0.5 - floor(spec->c_of_m[imol][i]+0.5));
        if(spec->quat != NULL)
           (void)printf("%9g %9g %9g %9g",spec->quat[imol][0],spec->quat[imol][1],
                 spec->quat[imol][2],spec->quat[imol][3]);
        (void)putchar('\n');
      }
   }
   (void)printf("end\n");

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
}
/******************************************************************************
 * random_pos.  Choose positions to be replaced randomly                      *
 ******************************************************************************/
void
random_pos(totmol, submol, subpos)
int             totmol, submol, subpos[];
{
   int		ranpos, subflag;
   int		i, j;

   srand(time(NULL)*rand());

   for( i = 0; i < submol; i++ )
   {   
       do
       {
          subflag = 0;
          ranpos = rand() % totmol;
          for ( j = 0; j < i; j++ ) 
             if( ranpos == subpos[j] )
                 subflag = 1;
       }
       while (subflag);

       subpos[i] = ranpos;
   }
}
/******************************************************************************
 * main().   Driver program for substituting species in MOLDY sys_spec files  *
 * Acceptable inputs are sys-spec files, restart files or dump files.         *
 * Call: ransub [-s sys-spec-file] [-r restart-file].                         *
 * If neither specified on command line, user is interrogated.                *
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
   int		rflag, mflag;
   int		irec;
   char		*filename = NULL, *dump_name = NULL;
   char		*dumplims = NULL;
   char		*molname = NULL;
   char		*tempname;
   char		dumpcommand[256];
   int		dump_size;
   float	*dump_buf;
   FILE		*Fp, *Dp;
   restrt_mt	restart_header;
   system_mt	sys;
   spec_mt	*species, *spec;
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf;
   contr_mt	control_junk;
   int          av_convert;
   int          maxmol;
   dopant       dop = {0, -1, -1.0, NULL, NULL, NULL};
   
#define MAXTRY 100
   control.page_length=1000000;

   while( (c = getopt(argc, argv, "cr:s:d:t:m:n:u:o:w:q:z:") ) != EOF )
      switch(c)
      {
       case 'c':
         cflg++;
         break;
       case 'r':
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
	 dumplims = mystrdup(optarg);
	 break;
       case 'm':
	 molname = strlower(mystrdup(optarg));
	 break;
       case 'n':
	 dop.mols = atoi(optarg);
	 break;
       case 'u':
	 dop.name = strlower(optarg);
	 break;
       case 'w':
         dop.mass = atof(optarg);
	 break;
       case 'q':
         dop.charge = optarg;
	 break;
       case 'z':
         dop.sym = optarg;
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
      fputs("Usage: ransub [-r restart-file | -s sys-spec-file] ",stderr);
      fputs("[-c] [-d dump-files] [-t s[-f[:n]]] [-m replaced species] ",stderr);
      fputs("[-u substituting species] [-n no of substitutions] ",stderr);
      fputs("[-w mass] [-q charge] [-z symbol] [-o output-file]\n",stderr);
      exit(2);
   }

   if( dump_name )
      data_source = 'd';

   if(intyp == 0)
   {
      fputs("How do you want to specify the simulated system?\n", stderr);
      fputs("Do you want to use a system specification file (1)", stderr);
      fputs(" or a restart file (2)", stderr);
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
   maxmol = sys.nmols;

  /*
   * Request species to be replaced if not already provided
   */
   do
   {
      mflag = 0;
      if( molname == NULL)
         if( dop.name != NULL || dop.mols > 0)
         {
            fputs("What is the name of the species to be replaced",stderr);
            molname = get_str("? ");
         }
         else
            mflag++;

      for(spec = species; spec < species+sys.nspecies; spec++)
         if( !strcmp(strlower(spec->name), molname) )
      {
          maxmol = spec->nmols;
          mflag++;
      }
      if(!mflag)
      {
          fprintf(stderr,"Species \"%s\" cannot be found\n", molname);
          (void)free(molname);
          molname = NULL;
      }
   } while (!mflag);

   if( dop.name == NULL && dop.mols > 0 )
   {
        fputs("What is the name of the substituting species ",stderr);
        dop.name = get_str("? ");
   }

   if( dop.name != NULL && dop.mols < 0 )
   {
        fprintf(stderr, "How many %s species do you want to replace", molname);
	dop.mols = get_int("? ",0,maxmol);
   }

   if( dop.mols < 0 )
        dop.mols = 0;

   if( dop.mols > maxmol )
      dop.mols = maxmol;

   (&dop)->pos = ialloc(dop.mols);         

   if( data_source == 0 )               /* If called interactively            */
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

   switch(data_source)                  /* To read configurational data       */
   {
    case 's':                           /* Lattice_start file                 */
        lattice_start(Fp, &sys, species, qpf);
        random_pos(maxmol, dop.mols, (&dop)->pos);
        sys_spec_out(&sys, species, molname, &dop, site_info, potpar, intyp);
      break;
    case 'r':                           /* Restart file                       */
        init_averages(sys.nspecies, restart_header.vsn,
                      control_junk.roll_interval, control_junk.roll_interval,
                      &av_convert);
        read_restart(Fp, restart_header.vsn, &sys, av_convert);
        random_pos(maxmol, dop.mols, (&dop)->pos);
        sys_spec_out(&sys, species, molname, &dop, site_info, potpar, intyp);
      break;
    case 'd':
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
      
    /*
     * Allocate buffer for data
     */
     dump_size = DUMP_SIZE(~0)*sizeof(float);

     if( (dump_buf = (float*)malloc(dump_size)) == 0)
       error("malloc failed to allocate dump record buffer (%d bytes)",
           dump_size);
#if defined (HAVE_POPEN) 
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

     for(irec = start; irec <= finish; irec+=inc)
     {
       if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
           error("Error reading record %d in dump file - \n%s\n",
              irec, strerror(errno));
        dump_to_moldy(dump_buf, &sys);  /*read dump data */

#ifdef DEBUG
        fprintf(stderr,"Sucessfully read dump record %d from file  \"%s\"\n",
	   irec, dump_name);
#endif
     /* Perform random substitution */
        random_pos(maxmol, dop.mols, (&dop)->pos);
        sys_spec_out(&sys, species, molname, &dop, site_info, potpar, intyp);
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