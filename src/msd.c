/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1997 Craig Fisher
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
#ifndef lint
static char *RCSid = "$Header$";
#endif
/**************************************************************************************
 * msd    	Code for calculating mean square displacements of centres of mass     *
 *              of molecules from MolDy dump files.			              *
 *		Output in columnar form "x y z total" for successive time intervals.  *
 *		Selection of species using -g: 1 = species 1, 2 = species 2, etc.     *
 *		Default msd time intervals:			     	              *
 *                             0 to (total no. of dump slices-1)/2, step size 1       *
 *		Option -w outputs trajectory coordinates in specified format          *
 *		"x y z" against time for each particle of selected species.           *
 *		nb. msd time intervals taken relative to extracted dump slices.       *
 ************************************************************************************** 
 *  Revision Log
 *  $Log$
 *  Revision 2.4.10.4  2003/10/21 10:27:55  kr
 *  Fixed a couple of bugs in the trajectory output
 *
 *  Revision 2.4.10.3  2003/07/31 02:55:58  moldydv
 *  Removed obsolete variables.
 *  Updated function descriptions to reflect changes in input syntax.
 *  Improved error message handling when problem with dump files.
 *  More commenting.
 *
 *  Revision 2.4.10.2  2003/07/30 09:30:53  moldydv
 *  Incorporated Keith's changes to read sysinfo from dump header.
 *  Removed ispec increment causing error in msd_out.
 *  Corrected absolute time step for msd output.
 *  Added more comments for clarity.
 *
 *  Revision 2.4.10.1  2003/07/29 09:38:20  moldydv
 *  Trajectory output options combined into one option, '-w'.
 *  Species now selected with -g in 'true' selector format.
 *
 *  Revision 2.4  2002/09/19 09:26:29  kr
 *  Tidied up header declarations.
 *  Changed old includes of string,stdlib,stddef and time to <> form
 *
 *  Revision 2.3  2000/12/06 10:47:33  keith
 *  Fixed call of make_sites() in utlsup.c to be compatible with new version.
 *  Tidied up declarations and added lint flags to reduce lint noise.
 *
 *  Revision 2.2  2000/11/09 16:54:13  keith
 *  Updated utility progs to be consistent with new dump format
 *
 *  Revision 2.1  2000/04/27 17:57:10  keith
 *  Converted to use full ANSI function prototypes
 *
 *  Revision 2.0  1999/11/18 09:58:33  keith
 *  checked in with -k by keith at 1999/11/25 14:27:58
 *
 *  Revision 2.0  1999/11/25  14:05:33  craig
 *  Selection of positions outside of limits now possible with -X,-Y,-Z
 *  Added new function "in_region" to utlsup.c
 *  Removed variable range_flag and replaced with extra dimension in range array.
 *  Corrected gnuplot output to only print blank lines when particle within region.
 *  Rewrote msd output header as one printf statement.
 *
 *  Revision 1.25  1999/11/18 09:58:33  keith
 *  Made headers into comments for gnuplot reading.
 *
 *  Revision 1.24  1999/11/15 11:50:39  keith
 *  Extended msd with spatial selection abilities of mdtraj.
 *
 *  Revision 1.23  1999/11/15  17:16:13  craig
 *  Extended position limits to apply to msd as well as traj calcs.
 *
 *  Revision 1.22  1999/11/12  11:37:51  craig
 *  Fixed bug in species number iteration.
 *
 *  Revision 1.22  1999/11/01 17:25:43  keith
 *  Got rid of initialization of array for pre-ANSI compilers.
 *
 *  Revision 1.21  1999/10/29 16:44:28  keith
 *  Updated usage message.
 *  Now uses mat_vec_mul() not mat_vec_mul3().
 *
 *  Revision 1.21  1999/10/25  15:34:32  craig
 *  Re-added zeroing of range_flag for non-ANSI machines.
 *  Tidied up usage message.
 *  Changed to generic matrix multiplier for converting coords.
 *
 *  Revision 1.20  1999/10/11 14:05:19  keith
 *  Removed common functions to "utlsup.c".
 *
 *  Revision 1.19  1999/10/08 10:54:51  keith
 *  Corrected it_inc limits for when only one or two time slices are selected.
 *  Reduced memory alloc to msd to account for selected species.
 *
 *  Revision 1.19  1999/10/01  10:59:09  craig
 *  Corrected it_inc limits for when only one or two time slices are selected.
 *  Reduced memory alloc to msd to account for selected species.
 *
 *  Revision 1.18  1999/07/22 13:35:59  keith
 *  Various fixes from Craig Fisher.
 *
 *  Revision 1.18  1999/06/03  18:00:09  craig
 *  Corrected memory freeing of dump, msd and species limits.
 *  Corrected for case when only one time slice selected.
 *  Modified init_inc input to read integer between limits.
 *
 *  Revision 1.17  1999/05/11  15:27:49  craig
 *  Corrected species iteration error in msd routine.
 *
 *  Revision 1.16  1998/10/23  10:31:23  craig
 *  Removed initialization of range_flag[3] since ANSI feature
 *  Removed unnecessary linking with lattice_start
 *  Corrected bug in get_real subroutine
 *
 *  Revision 1.15  1998/05/07 17:06:11  keith
 *  Reworked all conditional compliation macros to be
 *  feature-specific rather than OS specific.
 *  This is for use with GNU autoconf.
 *
 *  Revision 1.14  1998/01/28 09:55:05  keith
 *  Changed HAS_POPEN to more natural HAVE_POPEN.
 *  Fixed minor portability problem with struct initialization.
 *
 *  Revision 1.13  1998/01/27 17:46:58  keith
 *  Fixed minor portability problem with struct initialization.
 *
 *  Revision 1.12  1998/01/09 11:34:14  keith
 *  Added casts to arralloc() calls for portability.
 *  Changed to "HAVE_POPEN" macro from system-specifics
 *
 *  Revision 1.11  1997/11/27 16:00:38  keith
 *  Lintified it to avoid compiler warnings and complaints from
 *  nervous users.
 *
 *  Revision 1.10  1997/11/26 10:06:00  keith
 *  Corrected usage message.
 *  Made -r and -s options mutually exclusive.
 *
 *  Revision 1.9  1997/10/15 13:13:09  keith
 *  Minor tidying up by CF.
 *
 *  Revision 1.8  1997/10/13  11:16:10  craig
 *  Removed unused variable declarations
 *
 *  Revision 1.8  1997/10/09  11:21:40  craig
 *  Option for renaming program "mdtraj" added for default trajectory calculation
 *  Option 'c' added to parameter list to skip control information
 *  Changed hmat allocation from arralloc to aalloc
 *  Removed freeing of dumplims which was causing crash
 *  Msd limits modified to cope with dump limits interval = 1 (special case)
 *
 *  Revision 1.7  1997/10/08 13:30:55  keith
 *  Fixed dump_buf mem free bug.
 *
 *  Revision 1.6  1997/08/12 14:03:05  keith
 *  Combined version to produce output for GNUPLOT or IDL
 *
 *  Revision 1.2  1997/07/16 14:20:47  craig
 *  Option for various output formats for trajectory coords
 *
 *  Revision 1.1  1997/07/14 15:36:24  keith
 *  Modified by KR.  MSD calc put into separate function and optimised
 *  for roughly 4x speedup.
 *
 *  Revision 1.0  1997/07/11 16:55:26  craig
 *  Initial revision
 *
 */
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
#include "utlsup.h"
int     in_region(real *pos, real (*range)[3]);
int	getopt(int, char *const *, const char *);
extern  int optind;
/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
#define MSD  0
#define TRAJ 1
#define GNU  0
#define GEN  1
#define INNER 1
#define OUTER 2
/******************************************************************************
 * traj_gnu().  Output routine for displaying trajectories                    *
 *		- coords vs time for each species/atom for GNUplot            * 
 ******************************************************************************/
void
traj_gnu(spec_mt *species, vec_mt (**traj_cofm), int nslices, real (*range)[3],
	       	char *spec_mask, int nspecies)
{
   int          totmol=0, imol, i, itime;
   int		ispec=0;
   spec_mp      spec;

   for( spec = species, ispec=0; spec < species+nspecies; spec++, ispec++)
   {
      if( spec_mask[ispec] )
         for( imol = 0; imol < spec->nmols; imol++)
         {
            (void)printf("# %s\n",spec->name);
            if( in_region( traj_cofm[0][totmol], range) )
            {
               for( itime = 0; itime < nslices; itime++)
               {
                  for( i = 0; i < 3; i++)
                     (void)printf("%f ",traj_cofm[itime][totmol+imol][i]);
                  (void)printf("\n");
               }
               (void)printf("\n\n");
	    }
         }
      totmol += spec->nmols;
   }
   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
}
/******************************************************************************
 * traj_idl().  Output routine for displaying trajectories                    *
 *		- simple columnar format e.g. for IDL		              * 
 ******************************************************************************/
void
traj_idl(spec_mt *species, vec_mt (**traj_cofm), int nslices, real (*range)[3],
	       	char *spec_mask, int nspecies)
{
   int		totmol, imol, i, itime;
   int		ispec=0;
   spec_mp	spec;

   for( itime = 0; itime < nslices; itime++)
   {
      totmol = 0;
      for( spec = species, ispec=0; spec < species+nspecies; spec++, ispec++)
      {
	 if( spec_mask[ispec] )
            for( imol = 0; imol < spec->nmols; imol++)
               if( in_region( traj_cofm[0][totmol], range) )
                  for( i = 0; i < 3; i++)
                     (void)printf("%f ",traj_cofm[itime][totmol+imol][i]);
	 totmol += spec->nmols;
      }
      (void)printf("\n");
   }
   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
}
/***********************************************************************
 * msd_calc. Calculate msds from trajectory array		       *
 ***********************************************************************/    
void
msd_calc(spec_mt *species,	/* Species data */
         char *spec_mask,	/* Species selection mask */
         int nspecies,		/* No of species in system */
         int mstart,		/* Initial msd time slice */
  	 int mfinish,		/* Last msd time slice */
         int minc,		/* Msd time increment */
         int max_av,		/* Maximum no slices in average */
         int it_inc,		/* Initial time slice increment */
         real (*range)[3],	/* Spatial range included in msd's */
         vec_mt (**traj_cofm),	/* Continuous trajectory data */
         real ***msd)		/* Msd data for each time slice */
{
   int it, irec, totmol, imsd, ispec, imol, nmols, cmols, i;
   int nspec, tmol;
   spec_mp      spec;
   double       msdtmp, stmp;
   vec_mt	*tct0, *tct1;

   /* Outer loop for selecting initial time slice */
   for(it = 0; it <= (max_av-1)*it_inc; it+=it_inc)

      /* Inner loop for calculating displacement from initial slice */ 
      for(irec = mstart; irec <= mfinish; irec+=minc)
      {
	 imsd = (irec-mstart)/minc;
	 tct0 = traj_cofm[it];
	 tct1 = traj_cofm[it+irec];

	 for(i=0; i<3; i++)
	 {
	    totmol=0; ispec=0;
            for(spec = species; spec < species+nspecies; spec++, ispec++)
	    {
	       nmols = spec->nmols;
	       if( spec_mask[ispec] )
	       {
	          msdtmp = 0.0;
                  cmols= 0;
	          for( imol = 0; imol < nmols; imol++)
	          {
		     tmol = totmol+imol;
                     if( in_region( traj_cofm[0][tmol], range) )
                     {
		        stmp = tct1[tmol][i] - tct0[tmol][i] ;
		        msdtmp += SQR(stmp);
                        cmols++;
                     }
	          }
	          msd[imsd][ispec][i] += (cmols == 0 ? 0 : msdtmp / cmols);
	       }
	       totmol += nmols;
	    }
         }
      }
}
/******************************************************************************
 * msd_out().  Output routine for displaying msd results                      *
 ******************************************************************************/
void
msd_out(spec_mt *species,	/* Species data */
        real ***msd,            /* Msd data for each time slice */
        int max_av,             /* No of time slices averaged over */
        int nmsd,               /* No of msd time intervals */
        char *spec_mask,        /* Species selection mask */
        int nspecies,           /* No of species in system */
        double tstep)           /* Absolute time step in ps */
{
   int          ispec=0, imsd, i;
   real         totmsd;
   spec_mp      spec;

   /* Write species labels */
   (void)printf("#          ");
   for(spec = species; spec < species+nspecies; spec++, ispec++)
      if( spec_mask[ispec] )
         (void)printf("   %40s",spec->name);

   /* Write column labels */
   (void)printf("\n#         t ");
   for(spec = species, ispec=0; spec < species+nspecies; spec++, ispec++)
      if( spec_mask[ispec] )
         (void)printf("  dX(t)**2   dY(t)**2   dZ(t)**2   dR(t)**2    ");
   (void)printf("\n");

   /* Write msd data for each time slice */
   for( imsd = 0; imsd < nmsd; imsd++)
   {
      printf("%10.3f  ", imsd*tstep);
      for(spec = species, ispec=0; spec < species+nspecies; spec++, ispec++)
         if( spec_mask[ispec] )
         {
            totmsd = 0;
            for( i=0; i<3; i++)
            {
               msd[imsd][ispec][i] /= max_av;
               totmsd += msd[imsd][ispec][i];
               (void)printf(" %10.6f", msd[imsd][ispec][i]);
            }
            (void)printf(" %10.6f   ",totmsd);
         }
      (void)printf("\n");
   }

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
}
/******************************************************************************
 * main().  Driver program for calculating trajectories/msds from MOLDY dumps *
 * System and configurational info must be read from dump files. 	      *
 * Call: msd [-d dump-file] [-t dump-time-slices] [-m msd-time-slices         *
 * If not specified on command line, user is interrogated.		      *
 * Options [-x][-y][-z] prompt for limits in given direction to be applied    *
 *        when outputting trajectories. Molecules selected based on initial   *
 *	  positions.							      *
 ******************************************************************************/
contr_mt                control;

int
main(int argc, char **argv)
{
   int	c, cflg = 0;
   extern char	*optarg;
   int		errflg = 0;
   int		outsw = MSD, trajsw = GNU;
   int		start, finish, inc;			/* Range specifiers for dump files */
   int		mstart, mfinish, minc;			/* Range specifiers for msds */ 
   int		it_inc = 1;				/* Default increment for initial time slice */
   int		dflag, iflag, sflag, mflag;
   int		i, irec, ispec;				/* Counters */
   char		*filename = NULL, *dump_base = NULL;
   char		*dump_name = NULL, *dump_names = NULL;
   char		*dumplims = NULL;
   char		*msdlims = NULL;
   char		*tempname = NULL;
   char		dumpcommand[256];
   int		dump_size;
   float	*dump_buf;
   FILE         *Fp, *Dp, *dump_file;
   dump_mt	dump_header;
   dump_sysinfo_mt *dump_sysinfo;
   size_mt	sysinfo_size;
   system_mt	sys;
   spec_mt	*species, *spec;
   vec_mt 	**traj_cofm;		/* Cofm data for continuous trajectories */
   mat_mt	*hmat;			/* h matrix for each slice */
   real		range[3][3];		/* Spatial range to include in msd calcs */
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf;
   int          nmsd, max_av, nspecies, nslices;
   real         ***msd;			/* Calculated MSD data for each time slice and species */
   char         *spec_list = "1-50";    /* List of species to calculate for (with default) */
   char         spec_mask[MAX_SPECIES]; /* Mask for selecting species */
   int          arglen, ind, genflg;
   int          xdr = 0;
   real		msd_step; 		/* Absolute msd time interval (in ps) */

   zero_real(range[0],9);

#define MAXTRY 100
   control.page_length=1000000;

   comm = argv[0];
   if( strstr(comm, "msd") )
      outsw = MSD;
   else if( strstr(comm, "mdtraj") )
      outsw = TRAJ;

   while( (c = getopt(argc, argv, "cd:t:m:i:g:o:w:xXyYzZ") ) != EOF )
      switch(c)
      {
       case 'c':
         cflg++;
         break;
       case 'd':
	 dump_base = optarg;
	 break;
       case 't':
	 dumplims = mystrdup(optarg);
         break;
       case 'm':
	 msdlims = mystrdup(optarg);
         break;
       case 'g':
	 spec_list = optarg;
	 break;
       case 'i':
	 it_inc = atoi(optarg);
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error("failed to open file \"%s\" for output", optarg);
	 break;
       case 'x':
         range[0][2] = INNER;
         break;
       case 'y':
         range[1][2] = INNER;
         break;
       case 'z':
         range[2][2] = INNER;
         break;
       case 'X':
         range[0][2] = OUTER;
         break;
       case 'Y':
         range[1][2] = OUTER;
         break;
       case 'Z':
         range[2][2] = OUTER;
         break;
       case 'w':
	 if( !strcasecmp(optarg, "gnu") )
	    trajsw = GNU;
	 else if (!strcasecmp(optarg, "gen") )
	    trajsw = GEN;
	 outsw = TRAJ;
         break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fprintf(stderr,
         "Usage: %s [-c] [-d dump-files] [-t s[-f[:n]]] ",comm);
      fputs("[-m s[-f[:n]]] [-g s[-f[:n]]] ",stderr);
      fputs("[-i initial-time-increment] [-w trajectory-format] ",stderr);
      fputs("[-x|-X] [-y|-Y] [-z|-Z] [-o output-file]\n",stderr);
      exit(2);
   }

   if( tokenise(mystrdup(spec_list), spec_mask, MAX_SPECIES) == 0 )
      error("Invalid species specification \"%s\": usage eg 1,3,5-9,4",spec_list);

   /* Prepare dump file name for reading */
   genflg = 0;
   if( dump_base != 0 ) genflg++;

   if( genflg == 0 && optind < argc ) {
      arglen = 0;
      for(ind=optind; ind < argc; ind++) {
	 arglen += strlen(argv[ind]) + 1;
      }
      dump_names=malloc(arglen);
      dump_names[0] = 0;
      for(ind=optind; ind < argc; ind++) {
	 strcat(dump_names,argv[ind]);
	 if(ind < argc-1) strcat(dump_names," ");
      }
   }

   /* Dump dataset			      */
   if( dump_base == 0 && dump_names == 0)
   {
	fputs("Enter canonical name of dump files (as in control)\n",stderr);
	if( (dump_base = get_str("Dump file name? ")) == NULL) exit(2);
	genflg++;
   }

   if( dump_names == 0 ) dump_names = dump_base;

   /* Prepare to read header info */
   if( genflg )
   {
      dump_name = malloc(strlen(dump_base) + 2);
      sprintf(dump_name, dump_base, 0);
   }
   else
   {
      dump_name = argv[optind];
   }

   if( (dump_file = open_dump(dump_name, "rb")) == NULL)
   {
      fprintf(stderr, "Failed to open dump file \"%s\"\n", dump_name);
      exit(2);
   }

   /*
    * Read dump header. On first call we need to read only the first
    * part of sysinfo to determine nspecies and consequently the size
    * of the buffer needed to hold all of it.
    */
   sysinfo_size = sizeof(dump_sysinfo_mt);
   dump_sysinfo = (dump_sysinfo_mt*)malloc(sysinfo_size);
   if( read_dump_header(dump_name, dump_file, &dump_header, &xdr, 
			sysinfo_size, dump_sysinfo) 
       || ferror(dump_file) || feof(dump_file) )
   {
      fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
      exit(2);
   }
   sysinfo_size = sizeof(dump_sysinfo_mt) 
      + sizeof(mol_mt) * (dump_sysinfo->nspecies-1);
   (void)free(dump_sysinfo);

   /*
    * Allocate space for and read dump sysinfo.
    */
   dump_sysinfo = (dump_sysinfo_mt*)malloc(sysinfo_size);
   /*
    * Rewind and reread header, this time including sysinfo.
    */
   (void)rewind_dump(dump_file, xdr);
   if( read_dump_header(dump_name, dump_file, &dump_header, &xdr, 
			sysinfo_size, dump_sysinfo) 
       || ferror(dump_file) || feof(dump_file) )
   {
      fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
      exit(2);
   }

   sys.nspecies = dump_sysinfo->nspecies;
   sys.nmols    = dump_sysinfo->nmols;
   sys.nmols_r  = dump_sysinfo->nmols_r;

   species = aalloc(dump_sysinfo->nspecies, spec_mt );
   for( ispec = 0, spec = species; ispec < dump_sysinfo->nspecies; ispec++, spec++)
   {
      spec->nmols = dump_sysinfo->mol[ispec].nmols;
      spec->rdof = dump_sysinfo->mol[ispec].rdof;
      spec->framework = dump_sysinfo->mol[ispec].framework;
      strncpy(spec->name, dump_sysinfo->mol[ispec].name, L_spec);
   }

   allocate_dynamics(&sys, species);

  /*
   *  Ensure that the dump limits start, finish, inc are set up,
   *  either on command line or by user interaction.
   */
   do
   {
      dflag = 0;
      if( dumplims == NULL )
      {
          fputs("Please specify range of dump records in form", stderr);
          fputs(" start-finish:increment\n", stderr);
          dumplims = get_str("s-f:n? ");
      }
      if( forstr(dumplims, &start, &finish, &inc) )
      {
          dflag++;
          fputs("Invalid range for dump records \"", stderr);
          fputs(dumplims, stderr);
          fputs("\"\n", stderr);
      }
      if( dflag)
      {
          (void)free(dumplims);
          dumplims = NULL;
      } 
   } while(dflag);

   nslices = (finish-start)/inc+1;  /* no. of time slices in traj_cofm */

   /* Ensure initial time slice increment is valid for given dump range */
   do
   {
      iflag = 0;
      if( (it_inc <= 0) || it_inc > (nslices>2?nslices-2:1) )
      {
         fputs("Invalid initial time slice increment\n",stderr);
         fputs("Please specify initial time slice increment between",stderr);
         fprintf(stderr," 1 and %d\n",nslices>2?nslices-2:1);
         it_inc=get_int("Increment? ",1,nslices>2?nslices-2:1);
         iflag++;
      }
   } while(iflag);

  /*
   * Ensure that the msd limits mstart, mfinish, minc are set up,
   * either on command line or by user interaction.
   */
   if( msdlims != NULL)
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
         if( (mstart*inc > finish-start) || (mfinish*inc > finish-start))
         {
            mflag++;
            fputs("MSD interval exceeds dump range\n",stderr);
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
   else
   {
     /* Use default values for msd interval limits */
      mstart = 0;
      if(nslices == 2)
         mfinish = 1;
      else
         mfinish = floor((nslices-1)/2); /* Midpoint of longest time span */
      minc = 1;
   }

   nspecies = sys.nspecies;

  /*
   * Allocate buffer for data
   */
   dump_size = DUMP_SIZE(~0, sys.nmols, sys.nmols_r)*sizeof(float);

  /* Allocate memory for trajectory data and zero */
   traj_cofm = (vec_mt**)arralloc(sizeof(vec_mt),2,0,nslices-1,0,sys.nmols-1);
   zero_real(traj_cofm[0][0], nslices*sys.nmols*3);

  /* Allocate array to store unit cell matrices */
   hmat = aalloc(nslices, mat_mt);

   if( (dump_buf = (float*)malloc(dump_size)) == 0)
      error("malloc failed to allocate dump record buffer (%d bytes)",
          dump_size);
#if defined (HAVE_POPEN) 
   sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d %s",
        sys.nmols, sys.nmols_r, start, finish, inc, dump_names);
   
   if( (Dp = popen(dumpcommand,"r")) == 0)
        error("Failed to execute \'dumpext\" command - \n%s",
            strerror(errno));
#else
   tempname = tmpnam((char*)0);
   sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d -o %s %s",
         sys.nmols, sys.nmols_r, start, finish, inc, tempname, dump_names);
   system(dumpcommand);
   if( (Dp = fopen(tempname,"rb")) == 0)
        error("Failed to open \"%s\"",tempname);
#endif

/* Loop for calculating trajectories from current and previous time slices */ 
   for(irec = 0; irec <= finish-start; irec+=inc)
   {
        if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
           if( !strcmp(strerror(errno),"Success") )
              error("Error reading record %d in dump file \"%s\"\n",irec, dump_name);
           else
              error("Error reading record %d in dump file \"%s\" - \n%s\n",
                 irec, dump_name, strerror(errno));
        dump_to_moldy(dump_buf, &sys);  /* read dump data */

	memcpy(hmat[irec/inc], sys.h, sizeof(mat_mt));

        if( irec == 0)
	{
          range_in(&sys, range);
          traj_con2(species, (vec_mt*)0, traj_cofm[irec/inc], nspecies);
	}
	else
          traj_con2(species, traj_cofm[irec/inc-1], traj_cofm[irec/inc], nspecies);

#ifdef DEBUG
        fprintf(stderr,"Sucessfully read dump record %d from file  \"%s\"\n",
	   irec, dump_name);
#endif
   }
   xfree(dump_buf);

#if defined (HAS_POPEN) 
   pclose(Dp);
#else
   fclose(Dp);
   remove(tempname);
#endif

/* Convert trajectories from frac coords to Cartesian coords */
   for( i = 0; i < nslices; i++)
      mat_vec_mul(hmat[i], traj_cofm[i], traj_cofm[i], sys.nmols);

/*
 * Output either msd values or trajectory coords
 */
   if( outsw == MSD)
   {
  /* Calculate msd parameters */
     nmsd = (mfinish-mstart)/minc+1; /* No of msd time intervals */
     max_av = (nslices - mfinish)/it_inc; /* Max no of msd calcs to average over */
     msd_step = dump_sysinfo->deltat*inc*minc;

     if (max_av < 1)
          max_av = 1;

  /* Allocate memory for msd array and zero */
     msd = (real***)arralloc(sizeof(real),3,0,nmsd-1,0,nspecies-1,0,2);
     zero_real(msd[0][0],nmsd*nspecies*3);

  /* Calculate and print msd values */
     msd_calc(species, spec_mask, nspecies,  mstart, mfinish, minc, max_av, it_inc,
		     range, traj_cofm, msd);
     msd_out(species, msd, max_av, nmsd, spec_mask, nspecies, msd_step);
   }
   else /* Otherwise output trajectories in selected format */
     switch(trajsw)
     {
       case GEN:
          traj_idl(species, traj_cofm, nslices, range, spec_mask, nspecies);
          break;
       case GNU:
       default:
	  traj_gnu(species, traj_cofm, nslices, range, spec_mask, nspecies);
     }
   return 0;    
}
