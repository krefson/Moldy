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

#include "defs.h"
#include <stdarg.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "structs.h"
#include "utlsup.h"
#include "messages.h"

/*======================== Global variables ==================================*/
int ithread=0, nthreads=1;
/******************************************************************************
 * main().   Driver program for generating SCHAKAL input files from MOLDY     *
 * files.    Acceptable inputs are sys-spec files, or restart files. Actual   *
 * configurational info can be read from dump files, lattice-start files or   *
 * restart files.  Call: mdshak [-s sys-spec-file] [-r restart-file].   If    *
 * neither specified on command line, user is interrogated.		      *
 ******************************************************************************/
int
main(int argc, char **argv)
{
   int	c, cflg, ans_i, data_source = 0;
   char 	line[80];
   extern char	*optarg;
   int		errflg = 0;
   int		intyp = 0;
   int		start = 0, finish = 0, inc = 1;
   int		tflag = 0;
   int		irec;
   int		iout = 0;
   int		outsw=0;
   char		*filename = NULL, *dump_base = NULL;
   char		*dump_names = NULL;
   char		*dumplims = NULL;
   char		*insert = NULL;
   char		*tempname = NULL;
   char		*dumpcommand;
   int		dump_size;
   float	*dump_buf;
   FILE		*Fp = NULL, *Hp, *Dp;
   restrt_mt	restart_header;
   system_mt	sys;
   spec_mt	*species;
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf = NULL;
   int		av_convert;
   int		trajsw = 0;
   vec_mt       *prev_cofm = NULL;
   int          arglen, ind, genflg=0;
   int		verbose = 0;
   int		dump_level = 0;
   
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
   else if (strstr(comm, "mdarc") )
     outsw = ARC;
   else if (strstr(comm, "mdxtl") )
     outsw = XTL;
   else if (strstr(comm, "mdins") )
     outsw = SHELX;
   else
     outsw = OUTBIN;

   while( (c = getopt(argc, argv, "o:r:s:d:t:i:yf:v?") ) != EOF )
      switch(c)
      {
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
            error(NOOUTF, optarg);
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
	 dump_base = optarg;
	 break;
       case 't':
         if( tflag++ == 0)
           dumplims = mystrdup(optarg);
         else
           errflg++;
	 break;
       case 'i':
	 insert = optarg;
	 break;
       case 'y':
         trajsw = 1;
	 break;
       case 'f':
	 if( !strcasecmp(optarg, "shak") )
	    outsw = SHAK;
	 else if (!strcasecmp(optarg, "pdb") )
	    outsw = PDB;
	 else if (!strcasecmp(optarg, "xyz") )
	    outsw = XYZ;
	 else if (!strcasecmp(optarg, "dcd") || !strcasecmp(optarg, "vmd") )
	    outsw = DCD;
	 else if (!strcasecmp(optarg, "cssr") )
	    outsw = CSSR;
	 else if (!strcasecmp(optarg, "arc") )
	    outsw = ARC;
	 else if (!strcasecmp(optarg, "xtl") )
	    outsw = XTL;
         else if (!strcasecmp(optarg, "ins") )
            outsw = SHELX;
	 else if (!strcasecmp(optarg, "bin") )
	    outsw = OUTBIN;
	 break;
       case 'v':
         verbose++;
         break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fprintf(stderr,
	      "Usage: %s -s sys-spec-file|-r restart-file [-f out-type] [-y] ",
	      comm);
      fputs("[-d dump-files] [-t s[-f[:n]]] [-v] [-o output-file]\n", stderr);
      exit(2);
   }

   if( dump_base )
      data_source = 'd';

   if(intyp == 0)
   {
      fputs("How do you want to specify the simulated system?\n", stderr);
      fputs("Do you want to use a system specification file (1)", stderr);
      fputs(" or a restart file (2)\n", stderr);
      if( (ans_i = get_int("? ", 1, 2)) == EOF )
	 exit(2);
      intyp = ans_i-1 ? 'r': 's';

      if( (filename = get_str("File name? ")) == NULL )
	 exit(2);
   }

   switch(intyp)
   {
    case 's':
      if( (Fp = fopen(filename,"r")) == NULL)
	 error(NOSYSSPEC, filename);
      cflg = check_control(Fp);
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
	 error(NORESTART, filename, strerror(errno)); 
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
	if( dump_base == 0 )
	{
	   fputs("Enter canonical name of dump files (as in control)\n",stderr);
	   if( (dump_base = get_str("Dump file name? ")) == NULL)
             exit(2);
	}

        if( strstr(dump_base,"%d") )
           genflg++;

        /* Prepare dump file name for reading */

        if( genflg == 0 && optind < argc ) {
           arglen = strlen(dump_base);
           for(ind=optind; ind < argc; ind++) {
	      arglen += strlen(argv[ind]) + 1;
           }
           dump_names=malloc(arglen);
           dump_names[0] = 0;
           strcat(dump_names, dump_base);
           if(optind < argc) strcat(dump_names," ");
           for(ind=optind; ind < argc; ind++) {
              strcat(dump_names,argv[ind]);
              if(ind < argc-1) strcat(dump_names," ");
           }
        }
        else
           dump_names = dump_base;

	if( (dumpcommand = malloc(256+strlen(dump_names))) == 0)
	   error(COMMEM, 256+strlen(dump_names));

#if defined (HAVE_POPEN)
	sprintf(dumpcommand,"dumpext -c -1 %s%s", verbose?"-v ":"", dump_names);
	if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
	if( (Hp = popen(dumpcommand,"r")) == 0)
	   error(DUMPCOMM, dumpcommand);
#else
	tempname = tmpnam((char*)0);
	sprintf(dumpcommand,"dumpext -c -1 -o %s %s%s", tempname, verbose?"-v ":"", dump_names);
	if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
	system(dumpcommand);
	if( (Hp = fopen(tempname,"rb")) == 0)
	   error(FILEOPEN,tempname);
#endif
	finish = dump_info(Hp, &dump_level);

	(void)free(dumpcommand);

        if( verbose ) message(NULLI, NULLP, INFO, HEADER);

#if defined (HAS_POPEN)
	pclose(Hp);
#else
	fclose(Hp);
	remove(tempname);
#endif

   /* Check dump file contains necessary data */
        if( !(dump_level & 1) )
	  error(NOCOMP, "C of M positions", dump_level);

   /*
    *  Ensure that the dump limits start, finish, inc are set up.
    */
	if( tflag )
	  if( forstr(dumplims, &start, &finish, &inc) )
	    error(INVSLICES, dumplims);

	/*
	 * Allocate buffer for data
         */
	dump_size = DUMP_SIZE(~0,sys.nmols,sys.nmols_r)*sizeof(float);
	if( (dump_buf = (float*)malloc(dump_size)) == 0)
	   error(BUFFMEM, dump_size);
	if( (dumpcommand = malloc(256+strlen(dump_names))) == 0)
	   error(COMMEM, 256+strlen(dump_names));

	/*
	 * Loop over dump records, ascertaining which file they are in
	 * and opening it if necessary.  Call output routine.
	 */
#if defined (HAVE_POPEN) 
	sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d %s %s",
		sys.nmols,sys.nmols_r, start, finish, inc, dump_names, verbose?"-v ":"");
	if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
	if( (Dp = popen(dumpcommand,"r")) == 0)
	   error(DUMPCOMM, dumpcommand);
#else
	tempname = tmpnam((char*)0);
	sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d -o %s %s %s",
		sys.nmols,sys.nmols_r, start, finish, inc, tempname, dump_names, verbose?"-v ":"");
	if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
	system(dumpcommand);
	if( (Dp = fopen(tempname,"rb")) == 0)
	   error(FILEOPEN, tempname);
#endif
	for(irec = start; irec <= finish; irec+=inc)
	{
           if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
           {
              if( !strcmp(strerror(errno),"Success") )
                 error(DUMPREC, irec, dump_base, strerror(errno));
              else
                 error(DUMPREC0, irec, dump_base, strerror(errno));
           }

	   dump_to_moldy(dump_buf, &sys);

           if( trajsw )  /* Write coordinates for continuous trajectory */
           {
               if( irec == start )
               {
                   prev_cofm = aalloc(sys.nmols, vec_mt);
                   zero_real(prev_cofm[0], 3*sys.nmols);
               }
               traj_con(&sys, prev_cofm, irec-start);
           }

	   moldy_out(iout++, irec, inc, &sys, sys.h, species, site_info, outsw, intyp, insert);
#ifdef DEBUG
	   fprintf(stderr,"Successfully read dump record %d from file \"%s\"\n",
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
   if( verbose ) message(NULLI, NULLP, INFO, COMPLETE);
   return 0;    
}
