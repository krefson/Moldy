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
static char *RCSid = "$Header: /usr/users/moldy/CVS/moldy/src/mdvaf.c,v 1.13.10.1 2003/07/29 09:40:53 moldydv Exp $";
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
gptr	*arralloc(size_mt,int,...); 	/* Array allocator		      */

char	*strlower(char *s);
void	read_sysdef(FILE *file, system_mp system, spec_mp *spec_pp, 
                    site_mp *site_info, pot_mp *pot_ptr);
void	initialise_sysdef(system_mp system, spec_mt *species, site_mt *site_info, 
                          quat_mt (*qpf));
void	re_re_header(FILE *restart, restrt_mt *header, contr_mt *contr);
void	re_re_sysdef(FILE *restart, char *vsn, system_mp system, spec_mp *spec_ptr, 
                     site_mp *site_info, pot_mp *pot_ptr);
void	allocate_dynamics(system_mp system, spec_mt *species);
void	read_restart(FILE *restart, char *vsn, system_mp system, int av_convert);
void	init_averages(int nspecies, char *vsn, long int roll_interval, 
   long int old_roll_interval, int *av_convert);
int	getopt(int, char *const *, const char *);
void	tfree(gptr *p);
extern  int optind;
/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
contr_mt                control;

#define VAF  0
#define VTF  1
void	zero_float(float *r, int n)
{
   int i;
   for(i=0; i < n; i++)
      r[i] = 0.0;
}
double vdotf(int n, float *x, int ix, float *y, int iy)
{
   register double      dot=0.0;
   register int i, j;
   if( ix == iy && ix == 1)
   {
      for(i = 0; i < n; i++)
         dot += x[i] * y[i];
   }
   else
   {
      for(i = j = 0; i < n*ix; i += ix)
      {
         dot += x[i] * y[j];
         j += iy;
      }
   }
   return(dot);
}
/***********************************************************************
 * vaf_calc. Calculate vaf from velocity array		               *
 ***********************************************************************/    
void
vaf_calc(spec_mt *species, char *spec_mask, int nspecies, int vstart, int vfinish,
	 int vinc, int max_av, int it_inc, float (**vel)[3], float **vaf, int aflg, 
	 int cptflg)
{
   int it, irec, totmol, ivaf, ispec, icpt, nspec;
   spec_mp      spec;
   float	(*vel0)[3], (*vel1)[3];

   /* Outer loop for selecting initial time slice */
   for(it = 0; it <= (max_av-1)*it_inc; it+=it_inc)

      /* Inner loop for calculating displacement from initial slice */ 
      for(irec = vstart; irec <= vfinish; irec+=vinc)
      {
	 ivaf = (irec-vstart)/vinc;
	 vel0 = vel[it];
	 vel1 = vel[it+irec];
         totmol=0; ispec=0;
	 for( spec = species; spec < species+nspecies; spec++, ispec++)
            if( spec_mask[ispec] )
            {
               if( cptflg ) 
	       {
	          for(icpt=0; icpt < 3; icpt++)
	          {
	             if( aflg == 0 || spec->rdof > 0)
	             {
		        vaf[ivaf][3*ispec+icpt] += vdotf(spec->nmols, vel0[0]+icpt, 
					   3, vel1[0]+icpt, 3)/ spec->nmols;
   	             }
	          }
	       } else {
	           if( aflg == 0 || spec->rdof > 0)
		       vaf[ivaf][ispec] += vdotf(3*spec->nmols, vel0[0], 1, 
					  vel1[0], 1) / spec->nmols;
	       }
            }
      }
}
/***********************************************************************
 * vtf_calc. Calculate vtf from velocity array		               *
 ***********************************************************************/    
void
vtf_calc(spec_mt *species, char *spec_mask, int nspecies, int vstart, int vfinish, int vinc, 
	 int max_av, int it_inc, float (**vel)[3], float **vtf)
{
   int it, irec, totmol, ivtf, ispec, imol, i, nspec;
   spec_mp      spec;
   float	(*vel0)[3], (*vel1)[3];
   float	vtftmp0[3], vtftmp1[3];

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
            {
               for( i = 0; i < 3; i++)
                  vtftmp0[i] = vtftmp1[i] = 0;
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
vaf_out(spec_mt *species, float **vaf, int max_av, int nvaf, char *spec_mask, int nspecies, int aflg, int cptflg)
{
   int          ispec=0, ivaf;
   spec_mp      spec;

   for( spec=species; spec < species+nspecies; spec++, ispec++)
      if( spec_mask[ispec] )
      {
         if( aflg == 0 || spec->rdof > 0)
         {
            (void)printf("# %s\n",spec->name); 
            for( ivaf = 0; ivaf < nvaf; ivaf++)
            {
	       if( cptflg ) 
	       {
	         vaf[ivaf][3*ispec] /= max_av;
	         vaf[ivaf][3*ispec+1] /= max_av;
	         vaf[ivaf][3*ispec+2] /= max_av;
	         (void)printf("%6d %10.7f %10.7f %10.7f\n",ivaf, vaf[ivaf][3*ispec],
			vaf[ivaf][3*ispec+1],vaf[ivaf][3*ispec+2]);
      	       } else {
	         vaf[ivaf][ispec] /= max_av;
	         (void)printf("%10.7f\n",vaf[ivaf][ispec]);
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
   int	aflg = 0,c, cflg = 0, cptflg = 0;
   extern char	*optarg;
   int		errflg = 0;
   int		intyp = 0;
   int		start, finish, inc;
   int		vstart, vfinish, vinc;
   int		nslices;
   int		dflag, iflag, sflag, vflag;
   int		outsw;
   int		irec, ispec, it_inc = 1;
   char         *filename = NULL, *dump_base = NULL;
   char         *dump_name = NULL, *dump_names = NULL;
   char         *dumplims = NULL;
   char		*vaflims = NULL;
   char		*tempname;
   char		dumpcommand[256];
   int		dump_size;
   float	*dump_buf;
   FILE         *Fp, *Dp, *dump_file;
   dump_mt      dump_header;
   dump_sysinfo_mt *dump_sysinfo;
   size_mt      sysinfo_size;
   system_mt    sys;
   spec_mt      *species, *spec;
   site_mt      *site_info;
   pot_mt       *potpar;
   float        (**vel)[3];
   quat_mt	*qpf;
   contr_mt	control, control_junk;
   int          nvaf, max_av, nspecies;
   float        **vaf;
   char         *spec_list = "1-50";
   char         spec_mask[MAX_SPECIES];
   int          arglen, ind, genflg;
   int          xdr = 0;

#define MAXTRY 100
   control.page_length=1000000;

   comm = argv[0];
   if( strstr(comm, "mdvtf") )
     outsw = VTF;
   else
     outsw = VAF;


   while( (c = getopt(argc, argv, "3acd:t:v:i:g:o:q") ) != EOF )
      switch(c)
      {
       case '3':
	 cptflg++;
	 break;
       case 'a':
	 aflg++;
	 break;
       case 'c':
         cflg++;
         break;
       case 'd':
	 dump_base = optarg;
	 break;
       case 't':
	 dumplims = mystrdup(optarg);
         break;
       case 'v':
	 vaflims = mystrdup(optarg);
         break;
       case 'g':
         spec_list = optarg;
	 break;
       case 'i':
	 it_inc = atoi(optarg);
	 break;
       case 'q':
	 outsw = VTF;
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
      fprintf(stderr,"Usage: %s [-a] [-c] ",comm);
      fputs("-d dump-files -t s[-f[:n]] [-v s[-f[:n]]] ",stderr);
      fputs("[-g s[-f[:n]]] [-i init-inc] [-q] [-o output-file]\n",stderr);
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

   /* Dump dataset                            */
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
            fputs("VAF interval exceeds dump range\n",stderr);
         }
         if( vflag )
         {
            (void)free(vaflims);
            vaflims = NULL;
            fputs("Please specify VAF intervals in form", stderr);
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
   vel = (float (**)[3])arralloc(sizeof(float[3]),2,0,nslices-1,0,sys.nmols-1);
   zero_float(vel[0][0], nslices*sys.nmols*3);

   if( (dump_buf = (float*)malloc(dump_size)) == 0)
      error("malloc failed to allocate dump record buffer (%d bytes)",
          dump_size);
#if defined (HAVE_POPEN) 
   sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c %d -t %d-%d:%d %s",
        sys.nmols, sys.nmols_r, aflg?7:6, start, finish, inc, dump_name);
   
   if( (Dp = popen(dumpcommand,"r")) == 0)
        error("Failed to execute \'dumpext\" command - \n%s",
            strerror(errno));
#else
   tempname = tmpnam((char*)0);
   sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c %d -t %d-%d:%d -o %s %s",
         sys.nmols, sys.nmols_r, aflg?7:6, start, finish, inc, tempname, dump_name);
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
	memcpy(vel[irec/inc], dump_buf, dump_size);
#ifdef DEBUG
        fprintf(stderr,"Sucessfully read dump record %d from file  \"%s\"\n",
	   irec, dump_name);
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

     if (max_av < 1)
          max_av = 1;

     vaf = (float**)arralloc(sizeof(float),2,0,nvaf-1,0,
			     cptflg?3*nspecies-1:nspecies-1);
     zero_float(vaf[0],nvaf*(cptflg?3*nspecies:nspecies));

  /* Calculate and print vaf/vtf values */
     if( outsw )
        vtf_calc(species, spec_mask, nspecies, vstart, vfinish, vinc, max_av, it_inc, vel, vaf);
     else
        vaf_calc(species, spec_mask, nspecies, vstart, vfinish, vinc, max_av, it_inc, vel, vaf, aflg, cptflg);
     vaf_out(species, vaf, max_av, nvaf, spec_mask, nspecies, aflg, cptflg);
   return 0;    
}
