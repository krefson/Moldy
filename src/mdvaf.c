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
#ifndef lint
static char *RCSid = "$Header: /home/moldy/CVS/moldy/src/mdvaf.c,v 1.22 2005/03/06 18:26:18 cf Exp $";
#endif
/**************************************************************************************
 * mdvaf    	Code for calculating velocity autocorrelation functions (vaf) and     *
 *              velocity total correlation functions (vtf) from MolDy dump files.     *
 *		Output in columnar form "x y z total" for successive time intervals.  *
 *		Selection of species using -g: 1 = species 1, 2 = species 2, etc.     *
 *		Default mdvaf time intervals:			     	              *
 *                             1 to (total no. of dump slices-1)/2, step size 1       *
 *		nb. mdvaf time intervals taken relative to extracted dump slices.     *
 ************************************************************************************** 
 *  Revision Log
 *  $Log: mdvaf.c,v $
 *  Revision 1.22  2005/03/06 18:26:18  cf
 *  Returned velocity variables to floats.
 *  Moved vdotf to utlsup.c
 *  Removed unnecessary zeroing of vel variable.
 *
 *  Revision 1.21  2005/02/07 11:25:21  cf
 *  Changed floats to reals.
 *
 *  Revision 1.20  2005/02/04 14:52:15  cf
 *  Reads header info with dumpext to determine system info and maximum time slice range.
 *  Common utility messages/errors moved to utlsup.h.
 *
 *  Revision 1.19  2005/01/13 11:53:42  cf
 *  Removed unused -3 and -c options.
 *
 *  Revision 1.18  2005/01/13 11:30:15  cf
 *  Fixed formatting error in dumpext command line.
 *
 *  Revision 1.17  2005/01/11 17:06:27  kr
 *  Fixed error which failed to pass correct file args to "dumpext".
 *  Fixed real stinker of a buffer overflow on dumpcommand by making it dynamic.
 *  Added some helpful messages on the "-v" output.
 *
 *  Revision 1.16  2004/12/07 13:00:02  cf
 *  Merged with latest utilities.
 *
 *  Revision 1.13.10.5  2004/12/07 11:03:37  cf
 *  Added read_dump_header and made static.
 *  Added verbose option for dumpext.
 *
 *  Revision 1.13.10.4  2004/12/07 10:35:56  cf
 *  Incorporated Keith's corrections and additions.
 *
 *  Revision 1.13.10.3  2004/12/06 19:08:50  cf
 *  Removed unused variables.
 *  Removed option -c for skipping control info.
 *  Formatted output to be similar to msd.
 *  Choice of angular or linear velocities clearer.
 *  Automatically calculates x,y,z components.
 *
 *  Revision 1.13.10.2  2003/07/31 02:52:45  moldydv
 *  System info now read from dump header, not sys-spec or restart files.
 *  Updated function descriptions to reflect this.
 *  Improved error message handling for when dump level too low.
 *
 *  Revision 1.13.10.1  2003/07/29 09:40:53  moldydv
 *  Species now specified with -g in 'true' selector format.
 *
 *  Revision 1.13  2002/09/19 09:26:29  kr
 *  Tidied up header declarations.
 *  Changed old includes of string,stdlib,stddef and time to <> form
 *
 *  Revision 1.12  2001/05/18 11:02:00  keith
 *  Added "-3" option to calculate and print XYZ components separately
 *
 *  Revision 1.11  2001/02/19 18:06:04  keith
 *  Fix to work with a mixture of polyatomic and monatomic species.
 *
 *  Revision 1.10  2000/11/13 16:01:24  keith
 *  Changed dump format to contain principle-frame angular velocities.
 *  Adapted mdvaf.c to calculate angular acf's too - added "-a" flag.
 *
 *  Revision 1.9  2000/11/09 16:54:13  keith
 *  Updated utility progs to be consistent with new dump format
 *
 *  Revision 1.8  2000/04/27 17:57:10  keith
 *  Converted to use full ANSI function prototypes
 *
 *  Revision 1.7  1999/11/18 09:31:47  keith
 *  Fixed bug in zeroing vaf array.
 *
 *  Revision 1.7  1999/11/18  10:45:24  craig
 *  Fixed bug in zeroing vaf array.
 *
 *  Revision 1.6  1999/11/12  11:15:48  craig
 *  Fixed bug in species number iteration.
 *
 *  Revision 1.5  1999/11/01 17:24:16  keith
 *  Got rid of declaration of "comm" as it's already in "utlsup.h".
 *
 *  Revision 1.4  1999/10/29 16:44:28  keith
 *  Added function to read velocities from dump files.
 *
 *  Revision 1.4  1999/10/12 17:41:24  craig
 *  Added new function 'vel_to_moldy' to read velocity data
 *
 *  Revision 1.3  1999/10/11 14:13:29  keith
 *  Removed common functions to "utlsup.c".
 *  Removed some unused variables.
 *
 *  Revision 1.1  1999/10/11 10:50:07  keith
 *  Initial revision
 *
 *  Revision 1.0  1999/07/07 17:10:55  craig
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

double vdotf(int n, float *x, int ix, float *y, int iy);
/*======================== Global variables ==================================*/
extern int optind;
static int verbose;
int ithread=0, nthreads=1;

#define VAF  0
#define VTF  1

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
/***********************************************************************
 * vaf_calc. Calculate vaf from velocity array		               *
 ***********************************************************************/    
void
vaf_calc(spec_mt *species, char *spec_mask, int nspecies, int vstart, int vfinish,
	 int vinc, int max_av, int it_inc, float (**vel)[3], real **vaf, int aflg)
{
   register int i, it, irec, totmol, ivaf, ispec;
   spec_mp      spec;
   float	(*vel0)[3], (*vel1)[3];

   /* Outer loop for selecting initial time slice */
   for(it = 0; it <= (max_av-1)*it_inc; it+=it_inc)

      /* Inner loop for multiplying by current slice */ 
      for(irec = vstart; irec <= vfinish; irec+=vinc)
      {
	 ivaf = (irec-vstart)/vinc;
	 vel0 = vel[it];
	 vel1 = vel[it+irec];
         totmol=0; ispec=0;
	 for( spec = species; spec < species+nspecies; spec++, ispec++)
           {
           if( spec_mask[ispec] )
              if( aflg == 0 || spec->rdof > 0)
                 for(i=0; i < 3; i++)
	            vaf[ivaf][3*ispec+i] += vdotf(spec->nmols, vel0[totmol]+i, 
				   3, vel1[totmol]+i, 3)/spec->nmols;
           if( aflg == 0 || spec->rdof > 0)
              totmol+=spec->nmols;
           }
      }
}
/***********************************************************************
 * vtf_calc. Calculate vtf from velocity array		               *
 ***********************************************************************/    
void
vtf_calc(spec_mt *species, char *spec_mask, int nspecies, int vstart, int vfinish, int vinc, 
	 int max_av, int it_inc, float (**vel)[3], real **vtf, int aflg)
{
   register int it, irec, totmol, ivtf, ispec, imol, i;
   spec_mp      spec;
   float	(*vel0)[3], (*vel1)[3];
   real		vtftmp0[3], vtftmp1[3];

   /* Outer loop for selecting initial time slice */
   for(it = 0; it <= (max_av-1)*it_inc; it+=it_inc)

      /* Inner loop for calculating displacement from initial slice */ 
      for(irec = vstart; irec <= vfinish; irec+=vinc)
      {
	 ivtf = (irec-vstart)/vinc;
	 vel0 = vel[it];
	 vel1 = vel[it+irec];
         totmol=0; ispec=0;
	 for( spec = species; spec < species+nspecies; spec++, ispec++ )
            if( spec_mask[ispec] )
               if( aflg == 0 || spec->rdof > 0)
               {
                  for( i = 0; i < 3; i++)
                     vtftmp0[i] = vtftmp1[i] = 0.0;
	          for( imol = 0; imol < spec->nmols; totmol++, imol++)
                     for( i = 0; i < 3; i++)
                     {
                        vtftmp0[i] += spec->charge*vel0[totmol][i];
                        vtftmp1[i] += spec->charge*vel1[totmol][i];
                     }
	          vtf[ivtf][ispec] += DOTPROD(vtftmp0,vtftmp1);
               }
     }
}
/******************************************************************************
 * vaf_out().  Output routine for displaying vaf results                      *
 ******************************************************************************/
void
vaf_out(spec_mt *species, real **vaf, int max_av, int nvaf, char *spec_mask, int nspecies, int aflg, float tstep, int outsw)
{
   register int i, ispec=0, ivaf;
   real		total;
   spec_mp      spec;

   for( spec=species; spec < species+nspecies; spec++, ispec++)
      if( spec_mask[ispec] )
      {
         if( aflg == 0 || spec->rdof > 0)
         {
            (void)printf("#\n# %s\n",spec->name); 
            (void)printf("#     t    ");
            if( outsw )
              (void)printf("        VTF\n");
            else
              (void)printf("        X(t)           Y(t)           Z(t)          VAF(t)\n");
            for( ivaf = 0; ivaf < nvaf; ivaf++)
            {
              if( outsw )
              {
                 vaf[ivaf][ispec] /= max_av;
                 (void)printf("%10.3f %14.7f\n", ivaf*tstep, vaf[ivaf][ispec]);
              }
              else
              {
                 total = 0.0;
                 for( i=0; i<3; i++)
                 {
	            vaf[ivaf][3*ispec+i] /= max_av;
                    total += vaf[ivaf][3*ispec+i];
                 }
	            (void)printf("%10.3f %14.7f %14.7f %14.7f %14.7f\n", ivaf*tstep, vaf[ivaf][3*ispec],
	 		vaf[ivaf][3*ispec+1],vaf[ivaf][3*ispec+2], total);
              }
            }
         }
      }
   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
}
/******************************************************************************
 * main().  Driver program for calculating vaf/vtf from MOLDY dump files      *
 * Acceptable inputs are sys-spec files or restart files. Actual 	      *
 * configurational info must be read from dump files.			      *
 * Call: mdvaf [-s sys-spec-file] [-r restart-file] [-d dump-file] 	      *
 * If not specified on command line, user is interrogated.		      *
 ******************************************************************************/
int
main(int argc, char **argv)
{
   int	aflg = 0,c;
   extern char	*optarg;
   int		errflg = 0;
   int		start = 0, finish = 0, inc = 1;
   int		vstart, vfinish, vinc;
   int		nslices;
   int		iflag, vflag, tflag = 0;
   int		outsw;
   int		irec, it_inc = 1;
   char         *dump_base = NULL;
   char		*dump_names = NULL;
   char         *dumplims = NULL;
   char		*tempname = NULL;
   char		*vaflims = NULL;
   char		*dumpcommand;
   int		dump_size;
   float	*dump_buf;
   FILE         *Dp, *Hp;
   system_mt    sys;
   spec_mt      *species;
   float        (**vel)[3];
   int          nvaf, max_av, nspecies;
   real         **vaf;
   char         *spec_list = NULL;
   char         *spec_mask = NULL;
   int          arglen, ind, genflg=0;
   float	vaf_step = 0.0;
   double	delta_t;
   int          dump_level = 1;

#define MAXTRY 100

   comm = argv[0];
   if( strstr(comm, "mdvtf") )
     outsw = VTF;
   else
     outsw = VAF;

   verbose = 0;

   while( (c = getopt(argc, argv, "at:l:i:g:o:qv") ) != EOF )
      switch(c)
      {
       case 'a':  /* Calculate angular velocity function */
	 aflg++;
	 break;
       case 't':  /* Dump file limits */
         if( tflag++ == 0) 
	   dumplims = mystrdup(optarg);
         else
           errflg++;
         break;
       case 'l':  /* Limits for vtf/vaf calculation */
	 vaflims = mystrdup(optarg);
         break;
       case 'g':  /* Species selection */
         spec_list = optarg;
	 break;
       case 'i':  /* Time increment */
	 it_inc = atoi(optarg);
	 break;
       case 'q':  /* Calculate vtf rather than vaf */
         if( outsw == VAF)
	   outsw = VTF;
         else
           outsw = VAF;
	 break;
       case 'o':  /* Name of output file */
	 if( freopen(optarg, "w", stdout) == NULL )
            error(NOOUTF, optarg);
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
      fprintf(stderr,"Usage: %s [-a] [-q] [-t s[-f[:n]]] ",comm);
      fputs("[-l s[-f[:n]]] [-g s[-f[:n]] [-i init-inc] ",stderr);
      fputs("[-v] [-o output-file] dump-files\n",stderr);
      exit(2);
   }

   if( optind <= argc)
      dump_base = argv[optind];

   /* Dump dataset                            */
   if( dump_base == 0 && dump_names == 0)
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
   if( !(dump_level & 2) )
     error(NOCOMP, "C of M velocities", dump_level);

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

   nslices = (finish-start)/inc+1;  /* no. of time slices in vel */

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
   * Ensure that the vaf limits vstart, vfinish, vinc are set up,
   * either on command line or by user interaction.
   */
   if( vaflims != NULL)
      do
      {
         vflag = 0;
         if( forstr(vaflims, &vstart, &vfinish, &vinc) )
         {
           vflag++;
           fputs("Invalid range for vaf intervals \"", stderr);
           fputs(vaflims, stderr);
           fputs("\"\n", stderr);
         }
         if( (vstart*inc > finish-start) || (vfinish*inc > finish-start))
         {
            vflag++;
            fprintf(stderr,"%s interval exceeds dump range\n",(outsw?"VTF":"VAF"));
         }
         if( vflag )
         {
            (void)free(vaflims);
            vaflims = NULL;
            fprintf(stderr,"Please specify %s intervals in form", (outsw?"VTF":"VAF"));
            fputs(" start-finish:increment\n", stderr);
            vaflims = get_str("s-f:n? ");
         }
      } while(vflag);
   else
   {
     /* Use default values for vaf interval limits */
      vstart = 0;
      if(nslices == 2)
         vfinish = 1;
      else
         vfinish = floor((nslices-1)/2); /* Midpoint of longest time span */
      vinc = 1;
   }

   nspecies = sys.nspecies;

  /*
   * Allocate buffer for data
   */
   dump_size = 3*(aflg?sys.nmols_r:sys.nmols)*sizeof(float);

  /* Allocate memory for velocity data and zero */
   vel = (float (**)[3])arralloc(sizeof(float[3]),2,0,nslices-1,0,(aflg?sys.nmols_r:sys.nmols)-1);

   if( (dump_buf = (float*)malloc(dump_size)) == 0)
      error(BUFFMEM, dump_size);
   if( (dumpcommand = malloc(256+strlen(dump_names))) == 0)
      error(COMMEM, 256+strlen(dump_names));

#if defined (HAVE_POPEN) 
   sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c %d -t %d-%d:%d %s%s",
      sys.nmols, sys.nmols_r, aflg?7:6, start, finish, inc, verbose?"-v ":"", dump_names);
   if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
   
   if( (Dp = popen(dumpcommand,"r")) == 0)
     error(DUMPCOMM, dumpcommand);
#else
   tempname = tmpnam((char*)0);
   sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c %d -t %d-%d:%d -o %s %s%s",
       sys.nmols, sys.nmols_r, aflg?7:6, start, finish, inc, tempname, verbose?"-v ":"", dump_names);
   if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
   system(dumpcommand);
   if( (Dp = fopen(tempname,"rb")) == 0)
        error("Failed to open \"%s\"",tempname);
#endif

/* Loop for calculating trajectories from current and previous time slices */ 
   for(irec = 0; irec <= finish-start; irec+=inc)
   {
     if( verbose ) fprintf(stderr,"Reading dump record %d from dumpext output\n",irec);
     if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
     {
        if( !strcmp(strerror(errno),"Success") )
           error(DUMPREC, irec, dump_base, strerror(errno));
        else
           error(DUMPREC0, irec, dump_base, strerror(errno));
     }
     memcpy(vel[irec/inc], dump_buf, dump_size);  /* read dump data */

#ifdef DEBUG
     fprintf(stderr,"Sucessfully read dump record %d from file \"%s\"\n",
         irec, dump_base);
#endif
   }
   xfree(dump_buf);

#if defined (HAVE_POPEN) 
   pclose(Dp);
#else
   fclose(Dp);
   remove(tempname);
#endif

  /* Calculate vaf parameters */
     nvaf = (vfinish-vstart)/vinc+1; /* No of vaf time intervals */
     max_av = (nslices - vfinish)/it_inc; /* Max no of vaf calcs to average over */
     vaf_step = delta_t*inc*it_inc;

     if (max_av < 1)
        max_av = 1;

     vaf = (real**)arralloc(sizeof(real),2,0,nvaf-1,0,
                             outsw?nspecies-1:3*nspecies-1);
     zero_real(vaf[0],nvaf*(outsw?nspecies:3*nspecies));

  /* Calculate and print vaf/vtf values */
     if( outsw )
        vtf_calc(species, spec_mask, nspecies, vstart, vfinish, vinc, max_av, it_inc, vel, vaf, aflg);
     else
        vaf_calc(species, spec_mask, nspecies, vstart, vfinish, vinc, max_av, it_inc, vel, vaf, aflg);
     vaf_out(species, vaf, max_av, nvaf, spec_mask, nspecies, aflg, vaf_step, outsw);
   if( verbose ) message(NULLI, NULLP, INFO, COMPLETE);
   return 0;    
}
