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
static char *RCSid = "$Header: /home/moldy/CVS/moldy/src/msd.c,v 2.9 2005/01/13 11:29:49 cf Exp $";
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
 *  $Log: msd.c,v $
 *  Revision 2.9  2005/01/13 11:29:49  cf
 *  Fixed formatting error in dumpext command line.
 *
 *  Revision 2.8  2005/01/11 17:56:46  kr
 *  Fix for dumpcommand buffer overrun.
 *
 *  Revision 2.7  2004/12/07 13:00:01  cf
 *  Merged with latest utilities.
 *
 *  Revision 2.4.10.7  2004/12/07 11:03:37  cf
 *  Added read_dump_header and made static.
 *  Added verbose option for dumpext.
 *
 *  Revision 2.4.10.6  2004/12/07 10:35:56  cf
 *  Incorporated Keith's corrections and additions.
 *
 *  Revision 2.4.10.5  2004/12/06 19:10:07  cf
 *  System data read from dump header file.
 *  Removed unused variables.
 *
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

static int verbose;
int     in_region(real *pos, real (*range)[3]);
/*======================== Global variables ==================================*/
int ithread=0, nthreads=1;

#define MSD  0
#define TRAJ 1
#define GNU  0
#define GEN  1
#define INNER 1
#define OUTER 2

typedef struct list_mt
{
   struct list_mt       *next;
   int                  i;
   char *p;
   int num;
} list_mt;
/******************************************************************************
 * List manipulation procedures                                               *
 ******************************************************************************/void
insert(list_mt *entry, list_mt *head)
{
   while( head->next != NULL && entry->i > head->next->i)
      head = head->next;
                                                                                
   entry->next = head->next;
   head->next  = entry;
}
/******************************************************************************
 *  header_to_moldy. Extract data from header info file created by dumpext.   *
 ******************************************************************************/
static
void header_to_moldy(FILE *Fp, system_mt *sys, spec_mp *spec_pp,
                     double *delta_t, int *nslices, int *level)
{
   int      n=0, nspec=-1;
   int      idata;
   char     line[132], sdata[64];
   double   ddata;
   spec_mt  *species=NULL;

   *nslices = 0;

   while( !feof(Fp) )
   {
     get_line(line,132,Fp,1);

     if( sscanf(line, "File name\t\t\t= %s", sdata) > 0)
       n++;
     if( sscanf(line, "Time between dumps\t\t= %lf %s", &ddata, sdata) > 1)
         *delta_t = ddata;
     if( n == 1)
     {
       if( sscanf(line, "Dump level\t\t\t= %d", &idata) > 0)
         *level = idata;
       if( sscanf(line, "Number of particles\t\t= %d\n", &idata) > 0)
         sys->nmols = idata;
       if( sscanf(line, "Number of polyatomics\t\t= %d\n", &idata) > 0)
         sys->nmols_r = idata;
       if( sscanf(line, "Number of species\t\t= %d", &idata) > 0)
       {
         sys->nspecies = idata;
         *spec_pp = aalloc(sys->nspecies, spec_mt );
       }
 
       if( sscanf(line, "Species %d name\t\t\t= %s\n", &idata, sdata) > 1)
       {
         nspec++;
         species = *spec_pp+nspec;
         strncpy(species->name, sdata, L_spec);
       }
       if( sscanf(line, "  Number of ions\t\t= %d\n", &idata) > 0 ||
           sscanf(line, "  Number of molecules\t\t= %d\n", &idata) > 0)
             species->nmols = idata;

       if( sscanf(line, "  Ion is a framework\n") ||
           sscanf(line, "  Molecule is a framework\n"))
             species->framework = true;
       if( sscanf(line, "  Rotational deg. of freedom\t= %d", &idata) > 0)
          species->rdof = idata;
     }
     if( sscanf(line, "Number of dumps\t\t\t= %d", &idata) > 0)
       *nslices += idata;
   }
   *nslices-=1;
}
/******************************************************************************
 * traj_gnu().  Output routine for displaying trajectories                    *
 *		- coords vs time for each species/atom for GNUplot            * 
 ******************************************************************************/
void
traj_gnu(spec_mt *species, vec_mt (**traj_cofm), int nslices, real (*range)[3],
	       	char *spec_mask, int nspecies)
{
   register int totmol=0, imol, i, itime, ispec;
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
      error(WRITERR, strerror(errno));
}
/******************************************************************************
 * traj_idl().  Output routine for displaying trajectories                    *
 *		- simple columnar format e.g. for IDL		              * 
 ******************************************************************************/
void
traj_idl(spec_mt *species, vec_mt (**traj_cofm), int nslices, real (*range)[3],
	       	char *spec_mask, int nspecies)
{
   register int	totmol, imol, i, itime, ispec;
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
      error(WRITERR, strerror(errno));
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
   register int it, irec, totmol, imsd, ispec, imol, nmols, cmols;
   register int i, tmol;
   register double       msdtmp, stmp;
   spec_mp      spec;
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
   register int ispec=0, imsd, i;
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
      error(WRITERR, strerror(errno));
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
int
main(int argc, char **argv)
{
   int	c;
   extern char	*optarg;
   int		errflg = 0;
   int		outsw = MSD, trajsw = GNU;
   int		start = 0, finish = 0, inc = 1;	/* Range specifiers for dump files */
   int		mstart, mfinish, minc;		/* Range specifiers for msds */ 
   int		it_inc = 1;			/* Default increment for initial time slice */
   int		iflag, mflag, tflag = 0;
   int		i, irec;			/* Counters */
   char		*dump_base = NULL;
   char		*dump_names = NULL;
   char		*dumplims = NULL;
   char		*msdlims = NULL;
   char		*tempname = NULL;
   char		*dumpcommand;
   int		dump_size;
   float	*dump_buf;
   FILE         *Dp, *Hp;
   system_mt	sys;
   spec_mt	*species;
   vec_mt 	**traj_cofm;		/* Cofm data for continuous trajectories */
   mat_mt	*hmat;			/* h matrix for each slice */
   real		range[3][3];		/* Spatial range to include in msd calcs */
   int          nmsd, max_av, nspecies, nslices=0;
   real         ***msd;			/* Calculated MSD data for each time slice and species */
   char         *spec_list = NULL;      /* List of species to calculate for (with default) */
   char         *spec_mask = NULL;      /* Mask for selecting species */
   int          arglen, ind, genflg=0;
   real		msd_step; 		/* Absolute msd time interval (in ps) */
   double	delta_t;
   int		dump_level = 1;

   zero_real(range[0],9);

#define MAXTRY 100

   comm = argv[0];
   if( strstr(comm, "msd") )
      outsw = MSD;
   else if( strstr(comm, "mdtraj") )
      outsw = TRAJ;

   while( (c = getopt(argc, argv, "t:m:i:g:o:w:xXyYzZv") ) != EOF )
      switch(c)
      {
       case 't':
         if( tflag++ == 0)
	   dumplims = mystrdup(optarg);
         else
           errflg++;
         break;
       case 'm':
	 msdlims = mystrdup(optarg);
         break;
       case 'g':
	 spec_list = mystrdup(optarg);
	 break;
       case 'i':
	 it_inc = atoi(optarg);
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
            error(NOOUTF, optarg);
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
         "Usage: %s [-t s[-f[:n]]] [-m s[-f[:n]]] ",comm);
      fputs("[-g s[-f[:n]]] [-i initial-time-increment] ",stderr);
      fputs("[-w trajectory-format] [-x|-X] [-y|-Y] [-z|-Z] ",stderr);
      fputs("[-v] [-o output-file] dump-files\n",stderr);
      exit(2);
   }

   if( optind <= argc)
      dump_base = argv[optind];

   /* Dump dataset			      */
   if( dump_base == 0 )
   {
	fputs("Enter canonical name of dump files (as in control)\n",stderr);
	if( (dump_base = get_str("Dump file name? ")) == NULL)
           exit(2);
   }

   /* Prepare dump file name for reading */
   if( strstr(dump_base,"%d") )
      genflg++;

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
   header_to_moldy(Hp, &sys, &species, &delta_t, &finish, &dump_level);

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

   allocate_dynamics(&sys, species);

/* Check species selection list */
   spec_mask = (char*)calloc(sys.nspecies+1,sizeof(char));

   if( spec_list == NULL)
     {
     spec_list = malloc((int)log10(sys.nspecies) + 4);
     sprintf(spec_list,"1-%d",sys.nspecies);
     }

   if( tokenise(mystrdup(spec_list), spec_mask, sys.nspecies) == 0 )
      error(INVSPECIES, spec_list, sys.nspecies);

  /*
   *  Ensure that the dump limits start, finish, inc are set up.
   */
   if( tflag )
     if( forstr(dumplims, &start, &finish, &inc) )
        error(INVSLICES, dumplims);

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
     error(BUFFMEM, dump_size);
   if( (dumpcommand = malloc(256+strlen(dump_names))) == 0)
     error(COMMEM, 256+strlen(dump_names));

#if defined (HAVE_POPEN) 
   sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d %s%s",
     sys.nmols, sys.nmols_r, start, finish, inc, verbose?"-v ":"", dump_names);
   
   if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
   if( (Dp = popen(dumpcommand,"r")) == 0)
     error(DUMPCOMM, dumpcommand);
#else
   tempname = tmpnam((char*)0);
   sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d -o %s %s%s",
       sys.nmols, sys.nmols_r, start, finish, inc, tempname, verbose?"-v ":"", dump_names);
   if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
   system(dumpcommand);
   if( (Dp = fopen(tempname,"rb")) == 0)
     error(FILEOPEN, tempname);
#endif

/* Loop for calculating trajectories from current and previous time slices */ 
   for(irec = 0; irec <= finish-start; irec+=inc)
   {
     if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
     {
        if( !strcmp(strerror(errno),"Success") )
           error(DUMPREC, irec, dump_base, strerror(errno));
        else
           error(DUMPREC0, irec, dump_base, strerror(errno));
     }

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
     fprintf(stderr,"Successfully read dump record %d from file \"%s\"\n",
        irec, dump_base);
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
     msd_step = delta_t*inc*minc;

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
   if( verbose ) message(NULLI, NULLP, INFO, COMPLETE);
   return 0;    
}
