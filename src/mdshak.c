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
static char *RCSid = "$Header: /home/eeyore_data/keith/moldy/src/RCS/mdshak.c,v 2.20 1999/10/29 16:44:28 keith Exp $";
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
#include "ReadDCD.h"
#include "utlsup.h"
#ifdef HAVE_STDARG_H
gptr	*arralloc(size_mt,int,...); 	/* Array allocator		      */
#else
gptr	*arralloc();	        	/* Array allocator		      */
#endif

void    moldy_out();
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
void    zero_real();
/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
contr_mt		control;

#define SHAK   0
#define XYZ 1
#define OUTBIN 2
#define DCD 3
#define PDB 4
#define CSSR 5
#define DUMP_SIZE(level)  (( (level & 1) + (level>>1 & 1) + (level>>2 & 1) ) * \
           (3*sys.nmols + 4*sys.nmols_r + 9)+ (level>>3 & 1) * \
           (3*sys.nmols + 3*sys.nmols_r + 9) + (level & 1))
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
   int		trajsw = 0;
   vec_mt       *prev_cofm;
   
#define MAXTRY 100
   control.page_length=1000000;

   comm = argv[0];
   if( strstr(comm, "mdshak") )
     outsw = SHAK;
   else if (strstr(comm, "mdpdb") )
     outsw = PDB;
   else if (strstr(comm, "mdxyz") )
     outsw = XYZ;
   else if (strstr(comm, "mddcd") || strstr(comm, "mdvmd") )
     outsw = DCD;
   else if (strstr(comm, "mdcssr") )
     outsw = CSSR;
   else
     outsw = OUTBIN;

   while( (c = getopt(argc, argv, "a:d:o:cr:s:d:t:i:hxbvpyg") ) != EOF )
      switch(c)
      {
       case 'a': 
	 atom_sel = optarg;
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
	 dumplims = mystrdup(optarg);
	 break;
       case 'i':
	 insert = optarg;
	 break;
       case 'h':
	 outsw = SHAK;
	 break;
       case 'x':
	 outsw = XYZ;
	 break;
       case 'b':
	 outsw = OUTBIN;
	 break;
       case 'v':
	 outsw = DCD;
	 break;
       case 'p':
         outsw = PDB;
	 break;
       case 'y':
         trajsw = 1;
	 break;
       case 'g':
	 outsw = CSSR;
	 break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fprintf(stderr,
	      "Usage: %s [-x|-v|-h|-p|-g] [-y] [-c] [-s sys-spec-file|-r restart-file] ",
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
	moldy_out(0, 0, 1, &sys, sys.h, species, site_info, outsw, intyp, insert);
      break;
    case 'r':				/* Restart file			      */
	init_averages(sys.nspecies, restart_header.vsn,
		      control.roll_interval, control.roll_interval,
		      &av_convert);
	read_restart(Fp, restart_header.vsn, &sys, av_convert);
	moldy_out(0, 0, 1, &sys, sys.h, species, site_info, outsw, intyp, insert);
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

           if( trajsw )  /* Write coordinates for continuous trajectory */
           {
               if( irec == start )
               {
                   prev_cofm = aalloc(sys.nmols, vec_mt);
                   zero_real(prev_cofm, 3*sys.nmols);
               }
               traj_con(&sys, prev_cofm, irec-start);
           }

	   moldy_out(iout++, irec, inc, &sys, sys.h, species, site_info, outsw, intyp, insert);
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
