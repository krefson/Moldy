#if ANSI || __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#include <math.h>
#include <stdio.h>
#include "string.h"
#define error(str, args) message(NULLI, NULLP, FATAL, str, args)
#include "structs.h"
#include "messages.h"

char 	*malloc();
void	invert();
void	mat_vec_mul();
char	*arralloc();
void	make_sites();
char	*strlower();
void	read_sysdef();
void	initialise_sysdef();
void	re_re_header();
void	re_re_sysdef();
void	allocate_dynamics();
void	lattice_start();
void	read_restart();
/******************************************************************************
 * Dummies of 'moldy' routines so that mdshak may be linked with moldy library*
 ******************************************************************************/
void 	init_rdf()
{}
int	***rdf;
void	init_averages()
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
 *  message.   Deliver error message to possibly exiting. 		      *
 ******************************************************************************/
#if ANSI || __STDC__
#undef  va_alist
#define	va_alist int *nerrs, ...
#define va_dcl /* */
#endif
/*VARARGS*/
void	message(va_alist)
va_dcl
{
   va_list	ap;
   char		*buff;
   int		sev;
   char		*format;
#if ANSI || __STDC__
   va_start(ap, nerrs);
#else
   int		*nerrs;

   va_start(ap);
   nerrs = va_arg(ap, int *);
#endif

   buff  = va_arg(ap, char *);
   sev   = va_arg(ap, int);
   format= va_arg(ap, char *);

   (void)fprintf(stderr, "mdshak: ");
   (void)vfprintf(stderr, format, ap);
   fputc('\n',stderr);
   va_end(ap);

   if(buff != NULL)                     /* null ptr means don't print buffer  */
   {
      (void)printf("     buffer contents=\"%s\"",buff);
      new_line();
   }
   if(sev >= ERROR && nerrs != NULL)
      (*nerrs)++;
   if(sev == FATAL)
      exit(3);
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
char *av_ptr(size_p)
int *size_p;
{
   *size_p = 0;
   return NULL;
}
/******************************************************************************
 * forstr.  Parse string str of format s-f:n  (final 2 parts optional),       *
 *          returning integer values of s,f,n. f defaults to s and n to 1     *
 ******************************************************************************/
int
forstr(str, start, finish, inc)
char	*str;
int	*start, *finish, *inc;
{
   char	*p, *pp;
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
void
dump_to_moldy(buf, system)
float	*buf;
system_t *system;
{
   int i;
   float	*c_of_m = buf;
   float	*quat   = buf+3*system->nmols;
   float	*h      = buf+3*system->nmols + 4*system->nmols_r;
   mat_t hinv;

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
schakal_out(n, system, species, site_info)
int	n;
system_t	*system;
spec_t		species[];
site_t		site_info[];
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   spec_t	*spec;
   double	a, b, c, alpha, beta, gamma;
   mat_p	h = system->h;
   mat_t	hinv;
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

   (void)printf("END %d\n", n);
}
/******************************************************************************
 * Centre_mass.  Shift system centre of mass to origin (in discrete steps),   *
 ******************************************************************************/
void
centre_mass(system, species)
system_t	*system;
spec_t		species[];
{
   double	mass;
   vec_t	c_of_m;
   spec_t	*spec;
   int		imol;

   mass = c_of_m[0] = c_of_m[1] = c_of_m[2] = 0.0;
   for(spec = species; spec < species + system->nspecies; spec++ )
   {
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 c_of_m[0] += spec->mass*spec->c_of_m[imol][0];
	 c_of_m[1] += spec->mass*spec->c_of_m[imol][1];
	 c_of_m[2] += spec->mass*spec->c_of_m[imol][2];
      }
      mass += spec->nmols*spec->mass;
   }

   c_of_m[0] /= mass;
   c_of_m[1] /= mass;
   c_of_m[2] /= mass;
   c_of_m[0] = floor(c_of_m[0]+0.5);
   c_of_m[1] = floor(c_of_m[1]+0.5);
   c_of_m[2] = floor(c_of_m[2]+0.5);
   for(imol = 0; imol < system->nmols; imol++)
   {
      system->c_of_m[imol][0] -= c_of_m[0];
      system->c_of_m[imol][1] -= c_of_m[1];
      system->c_of_m[imol][2] -= c_of_m[2];
   }
}
/******************************************************************************
 * moldy_out.  Select output routine and handle file open/close		      *
 ******************************************************************************/
void
moldy_out(n, system, species, site_info)
int	n;
system_t	*system;
spec_t		species[];
site_t		site_info[];
{
   centre_mass(system, species);
   schakal_out(n, system, species, site_info);
}
/******************************************************************************
 * main().   Driver program for generating SCHAKAL input files from MOLDY     *
 * files.    Acceptable inputs are sys-spec files, or restart files. Actual   *
 * configrational info can be read from dump files, lattice-start files or    *
 * restart files.  Call: mdshak [-s sys-spec-file] [-r restart-file].   If    *
 * neither specified on command line, user is interrogated.		      *
 ******************************************************************************/
contr_t	control;
unit_t	input_unit;
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
   int		idump, idump0, jdump, irec;
   int		iout = 0;
   char		*filename = NULL, *dump_name = NULL;
   char		*dumplims = NULL;
   char		cur_dump[256];
   int		dump_size;
   float	*dump_buf;
   FILE		*Fp, *Dp;
   dump_t	header;
   restrt_t	restart_header;
   system_t	system;
   spec_t	*species;
   site_t	*site_info;
   pot_t	*potpar;
   quat_t	*qpf;
   contr_t	control_junk;
   control.page_length=1000000;
#define MAXTRY 100

   while( (c = getopt(argc, argv, "o:cr:s:d:n:") ) != EOF )
      switch(c)
      {
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error("failed to open file \"%s\" for output", optarg);
	 break;
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
       case 'n':
	 dumplims = optarg;
	 break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fputs("Usage: mdshak [-s sys-spec-file] [-r restart-file] ",stderr);
      fputs("[-d dump-files] [-n s[-f[:n]]]\n", stderr);
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
      {
	 error("Couldn't open sys-spec file \"%s\" for reading\n", filename);
	 exit(2);
      }
      if( cflg )
      {
	 do
	 {
	    fscanf(Fp, "%s",line);
	    (void)strlower(line);
	 }
	 while(! feof(stdin) && strcmp(line,"end"));
      }
      read_sysdef(Fp, &system, &species, &site_info, &potpar);
      qpf = qalloc(system.nspecies);
      initialise_sysdef(&system, species, site_info, qpf);
      break;
    case 'r':
      if( (Fp = fopen(filename,"rb")) == NULL)
      {
	 error("Couldn't open restart file \"%s\" for reading\n", filename);
	 exit(2);
      }
      re_re_header(Fp, &restart_header, &control_junk);
      re_re_sysdef(Fp, &system, &species, &site_info, &potpar);
      break;
    default:
      error("Internal error - invalid input type", "");
   }
   allocate_dynamics(&system, species);

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
	lattice_start(Fp, &system, species, qpf);
	moldy_out(1, &system, species, site_info);
      break;
    case 'r':				/* Restart file			      */
	read_restart(Fp, &system);
	moldy_out(1, &system, species, site_info);
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
		
	idump0 = -1;
	do			/* Search for a dump file matching pattern */
	   sprintf(cur_dump, dump_name, ++idump0);
	while( (Dp = fopen(cur_dump, "rb")) == NULL && idump0 < MAXTRY);
	if( Dp == NULL )	/* If we didn't find one . . 		   */
	   error("I can't find any dump files to match \"%s\".",dump_name);
	/*
	 * At this stage we should have the first dump file open.
	 * Read the header.
	 */
        fread(&header, sizeof header, 1, Dp);
	if( ferror(Dp) )
	   error("Couldn't read dump header","");
	if( (header.dump_level & 1) == 0 )
	   error("Dump at level %d doesn't contain co-ordinate information",
		 header.dump_level);
	/*
	 * Allocate buffer for data
         */
	dump_size = header.dump_size*sizeof(float);
	dump_buf = (float*)malloc(dump_size);
	idump = idump0;
	/*
	 * Loop over dump records, ascertaining which file they are in
	 * and opening it if necessary.  Call output routine.
	 */
	for(irec = start; irec <= finish; irec+=inc)
	{
	   jdump = irec/header.maxdumps + idump0;   /* Which file is rec. in? */
	   if( jdump != idump )			    /*   currently open file? */
	   {
	      (void)fclose(Dp);
	      sprintf(cur_dump, dump_name, jdump);  /* Make new name          */
	      if( freopen(cur_dump, "rb", Dp) == NULL )
		 error("Failed to open dump file \"%s\"", cur_dump);
	      idump = jdump;			    /* Mark as current dump   */
	   }
	   /*
	    *  Now go and get the data
	    */
	   fseek(Dp, sizeof header + irec%header.maxdumps*dump_size, 0);
	   fread(dump_buf, dump_size, 1, Dp);
	   if(ferror(Dp))
	      error("Error reading record %d in dump file",
		    irec%header.maxdumps);

	   dump_to_moldy(dump_buf, &system);

	   moldy_out(iout++, &system, species, site_info);
#ifdef DEBUG
	   fprintf(stderr,"Sucessfully read dump record %d from file  \"%s\"\n",
		   irec%header.maxdumps, cur_dump);
#endif
	}
	break;
      default:
	break;
     }
   return 0;    
}
      

		   
			     
