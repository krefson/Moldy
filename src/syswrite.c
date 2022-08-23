/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 2001 Craig Fisher
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
 * syswrite code for converting, SHAKAL, CSSR or PDB coords into Moldy sys-spec files *
 *              -y selects .pot file from which potentials are read                   *
 **************************************************************************************
 */
#include "defs.h"
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include "structs.h"
#include "messages.h"
#include "utlsup.h"
#include "specdata.h"
#include "sginfo.h"
#include "readers.h"
#include "elements.h"

/*========================== External data references ========================*/
extern  const pots_mt   potspec[];           /* Potential type specification  */

/*======================== Global variables ==================================*/
int ithread=0, nthreads=1;

/******************************************************************************
 * add_suffix().  Add numerical suffix to string.                             *
 ******************************************************************************/
int     add_suffix(char *string, int number)
{
   int     n_digit;
   int     length = strlen(string);

   n_digit = floor(log10(number))+1;
   if( length + n_digit > NLEN-1)       /* Longer than allowable name length? */
      string[NLEN-n_digit-1] = '\0';

   sprintf(string, "%s%d", string, number);

   return 0;
}
/******************************************************************************
 * copy_atom_data()  Make copy of atom/molecule names and positions           *
 ******************************************************************************/
int copy_atom_data(site_mp site, double *charge, char (*label)[NLEN], int natoms)
{
   int iatom;
   for(iatom = 0; iatom < natoms; iatom++)
   {
       strncpy(site[iatom].name, label[iatom], L_site);
       site[iatom].charge = charge[iatom];
       site[iatom].mass = 0.0;
       site[iatom].pad = -1;
   }
   return 0;
}
/******************************************************************************
 * main().   Driver program for converting files                              *
 ******************************************************************************/
int
main(int argc, char **argv)
{
   int		u, i, j, k;
   extern char	*optarg;
   int		errflg = 0;
   char		*filename = NULL, *specname = NULL;
   char         *filetype = "";
   char		*elefile = NULL;
   char		*potfile = NULL;
   char		title[TITLE_SIZE] = "";
   char		*buffer[NLEN];
   double	cell[6];                      /* Cell parameters */
   mat_mt	h;                            /* Hessian */
   int		insw=-1;                      /* Switch for input format */
   int		typetot=0, atomtot=0;         /* No of species, atoms */
   int		idij, ip;
   int          cflag = 0;                    /* Flag for periodic system */
   int		nflag;                        /* Flag for site type matching */
   int          n_elem;                       /* No of records read from element data file */
   spec_data	elem_data[NELEM];             /* Element info */
   const ele_data   *elem;                    /* Pointer to element info */
   spec_mt 	spec[MAX_SPECIES];            /* Species info */
   site_mt      *site, *st, specsite[MAX_SPECIES];  /* Site info */ 
   pot_mt       *pot_par;                     /* Potential parameters */
   int		ptype = -1, n_potpar=NPOTP;
   char         spgr[16];                     /* Space Group in Herman Maugain form */
   int          repeat[3] = {0,0,0};          /* Simulation box repeat factors from command line*/
   double       simbox[3] = {1,1,1};          /* Simulation box repeat factors */
   vec_mt	x[MAX_ATOMS]; 	              /* C of M coordinates */
   double       charge[MAX_ATOMS];            /* Site charge array */
   char         label[MAX_ATOMS][NLEN];       /* Site name array */
   int		num_mols = 1;		      /* No of molecules for non-lattice start species */

   comm = argv[0];
   while( (u = getopt(argc, argv, "i:o:y:e:n:l:a:b:c:?") ) != EOF )
      switch(u)
      {
       case 'i':
	 filename = optarg;
         filetype = read_ftype(filename);
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error(NOOUTF, optarg);
	 break;
       case 'y':
         potfile = optarg;
	 break;
       case 'e':
	 elefile = optarg;
         break;
       case 'n':
	 num_mols = atoi(optarg);
         if( num_mols < 1 )
           num_mols = 1;
         break;
       case 'l':
	 specname = optarg;
         break;
       case 'a':
	 repeat[0] = atoi(optarg);
         break;
       case 'b':
	 repeat[1] = atoi(optarg);
         break;
       case 'c':
	 repeat[2] = atoi(optarg);
         break;
       default:
       case '?': 
	 errflg++;
      }

   if( errflg )
   {
      fprintf(stderr,"Usage: %s [-i input-file] [-o output-file] ",comm);
      fputs("[-n no-of-molecules] [-l species-label] [-e element-data-file] ",stderr);
      fputs("[-y potential-parameter-file] [-a a-direction-cell-repeat-factor] ",stderr);
      fputs("[-b a-direction-cell-repeat-factor] [-c c-direction-cell-repeat-factor]\n", stderr);

      exit(2);
   }

   /* If no filename given on command line, request from user */
   while( insw < 0)
   {
      if( !strcasecmp(filetype, "cssr") )
           insw = CSSR;
      else if( !strcasecmp(filetype, "pdb") )
           insw = PDB;
      else if( !strcasecmp(filetype, "shak") )
           insw = SHAK;
      else if( !strcasecmp(filetype, "ins") || !strcasecmp(filetype, "res"))
           insw = SHELX;
      else if( !strcasecmp(filetype, "xtl") )
           insw = XTL;
      else if( !strcasecmp(filetype, "xyz") )
           insw = XYZ;

      if( insw < 0)
      {
         if( (filename = get_str("Structure file name? (PDB, CSSR, SCHAKAL, SHELX, XTL, XYZ) ")) == NULL )
            exit(2);
         filetype = read_ftype(filename);
      }
   }

   zero_real(h[0],9);
   zero_double(charge,MAX_ATOMS);
   strcpy(spgr,"P 1");

   switch(insw)  /* Read in data according to format selected */
   {
      case CSSR:
         atomtot = read_cssr(filename, h, label, x, charge, title, spgr);
	 break;
      case PDB:
         atomtot = read_pdb(filename, h, label, x, charge, title, spgr);
	 break;
      case SHAK:
         atomtot = read_shak(filename, h, label, x, charge, title, simbox);
	 break;
      case SHELX:
         atomtot = read_shelx(filename, h, label, x, title, spgr);
	 break;
      case XTL:
         atomtot = read_xtl(filename, h, label, x, charge, title, spgr);
	 break;
      case XYZ:
         atomtot = read_xyz(filename, h, label, x, title);
	 break;
      default:
         error(UNKSTRUCT, filename);
   }

   if( atomtot <= 0  )
      error("No known atoms found");

   if( det(h) != 0)
   {
      cflag++;  /* Periodic system */
   }
   else
   {            /* Non-periodic system */
      spec[0].nmols = num_mols;
      if( specname != NULL )
        j = get_tokens(specname, buffer, "\\/");
      else
        j = get_tokens(filename, buffer, "\\/. ");

      if( j > 1 )
        specname = mystrdup(*(buffer+j-2));
      else
        specname = mystrdup(*buffer);
      strncpy(spec[0].name, specname, (strlen(specname) > NLEN ? NLEN : strlen(specname)));
   }

   if( strcmp(spgr,"P 1") && atomtot > 0 && cflag )
      atomtot = sgexpand(MAX_ATOMS, atomtot, x, label, charge, spgr);

   /* Check if element data file is available */
   if( elefile ) 
     if( n_elem = read_ele(elem_data, elefile) < 0 )
       message(NULLI,NULLP,WARNING,NOELEM,elefile);

   site = (site_mt*)arralloc(sizeof(site_mt),1,0,atomtot-1);
   copy_atom_data(site, charge, label, atomtot);

   /* Loop to determine no of different species */
   for(st = site; st < site+atomtot; st++)
   {
      strlower(st->name);                    /* Remainder of symbol lowercase */
      st->name[0] = xtoupper(st->name[0]);   /* First character of symbol uppercase */

      nflag = 0;   /*  Flag for whether site assigned (1) or not assigned (0) to species */
      if( cflag )  /* Lattice start: Treat each atom as single species */
      {
        if( typetot > 0 )
           for( i = 0; i < typetot; i++)
              /* Does site match any species identified so far? */
              if( !strcmp(st->name, specsite[i].name) &&
                 ( insw == SHAK || insw == XYZ || st->charge == specsite[i].charge ) )
              {
                   st->pad = i;                 /* Assign this site to species i */
                   spec[i].nmols++;
                   nflag++;
              }
        if (!nflag)                          /* If no matches, create new species type */
        {
           if( typetot > MAX_SPECIES )
             error("Too many species found - current limit is %d (MAX_SPECIES in specdata.h)", MAX_SPECIES);

           /* First check user-specified element data file */
           if( elefile )
             for( j=0; j < n_elem; j++)        /* Search element/species data for match with site */
               if( !strcmp(st->name, elem_data[j].symbol) )
               {
                 strncpy(spec[typetot].name,elem_data[j].name,NLEN);
                 st->mass = elem_data[j].mass;
                 if( insw == SHAK || insw == XYZ)
                   st->charge = elem_data[j].charge;
                 nflag++;
                 break;
               }

           /* If not found, check default element data */
           if( !nflag)
             for (elem = ElementData; elem->name; elem++)
               if(  (!strcasecmp(st->name, elem->symbol)) )
               {
                 strncpy(spec[typetot].name,elem->name, NLEN);
                 st->mass = elem->mass;
                 if( insw == SHAK || insw == XYZ)
                   st->charge = elem->charge;
                 nflag++;
                 break;
               }

           if( !nflag)
             strcpy(spec[typetot].name, st->name);

           strcpy(specsite[typetot].name, st->name);
           specsite[typetot].mass = st->mass;
           specsite[typetot].charge = st->charge;
           spec[typetot].nmols=1;
           st->pad = typetot;
           typetot++;
        }
      }
      else /* Treat as single species */
      {
         for( i = 0; i < typetot; i++)  /* typetot now refers to number of different site types */
           if( !strcmp(st->name, specsite[i].name) &&        /* Does site match any species identified so far? */
              ( insw == XYZ || st->charge == specsite[i].charge ) )
             {
               st->pad = i;
               nflag++;
               break;
             }
         if (!nflag)                          /* If no matches, create new site type */
         {
           /* First check user-specified element data file */
           if( elefile )
             for( j=0; j < n_elem; j++)        /* Search element/species data for match with site */
               if( !strcmp(st->name, elem_data[j].symbol) )
               {
                 st->mass = elem_data[j].mass;
                 nflag++;
                 break;
               }

           /* If not found, check default element data */
           if (!nflag)
             for (elem = ElementData; elem->name; elem++)
               if(  (!strcasecmp(st->name, elem->symbol)) )
                 if( !strcmp(st->name, ElementData[j].symbol) )
                 {
                   st->mass = ElementData[j].mass;
                   nflag++;
                   break;
                 }

           /* Create new site entry */
           strcpy(specsite[typetot].name, st->name);
           specsite[typetot].mass = st->mass;
           specsite[typetot].charge = st->charge;
           st->pad = typetot;
           typetot++;
         }
      }
   }

   /* Modify species name if same name but different charge to another species */
   if( cflag )  /* Lattice start: Treat each atom as single species */
     for( i = 0; i < typetot-1; i++)
     {
       k = 0;
       for(j = i+1; j < typetot; j++)
         if( !strcmp(spec[i].name, spec[j].name) )
         {
            k++;
            add_suffix(spec[j].name, k);
         }
     }
   /* Initialize potential parameter arrays */
   pot_par = aalloc(SQR(typetot), pot_mt);

   for( i = 0; i < SQR(typetot); i++)
   {
      pot_par[i].flag = 0;
      zero_real(pot_par[i].p, NPOTP);
   }

   /* Read potential parameter file if available */
   if( (potfile != NULL) && ((ptype = read_pot(potfile, &pot_par, specsite, 0, typetot)) < 0) )
       message(NULLI,NULLP,WARNING,NOPOTL,potfile);

   /* nb. Command line repeat factors override those from Schakal file*/
   for( i=0; i<3; i++)
      if( repeat[i] > 0)
         simbox[i] = (double)repeat[i];

   /* Write data to file in Moldy input form */
   if (strchr(title,'\n') != NULL )
      printf("# %s",title);
   else
      if(strlen(title) > 0)
         printf("# %s\n",title);

   printf("# System specification file written by SYSWRITE on %s\n",atime());

   if( cflag )
   {
     for( i=0; i < typetot; i++)
     {
       printf("%-14s %d\n",spec[i].name, (int)(spec[i].nmols*simbox[0]*simbox[1]*simbox[2]));
       printf("%-3d    0    0    0  %12.12g %12.12g  %s\n",i+1,
                     specsite[i].mass,specsite[i].charge,specsite[i].name);
     }
   }
   else
   {
     printf("%-14s %d\n",spec[0].name, spec[0].nmols);
     st = site;
     for(i = 0; i < atomtot; i++)
     {
        printf("%-3d  %10g  %10g  %10g  %12.12g %12.12g  %s\n", st->pad+1,
               x[i][0], x[i][1], x[i][2],
               specsite[st->pad].mass, specsite[st->pad].charge, st->name);
        st++;
     }
   }
   puts("end");

   if( ptype >= 0)
   {
      n_potpar = potspec[ptype].npar;
      (void)printf("%s\n",potspec[ptype].name);
   }
   else
      (void)puts("Potential parameters");

   /* Write potential parameters for pairs of site_ids */
   for( i = 0; i < typetot; i++)
      for( j = i; j < typetot; j++)
      {
        idij = j+i*typetot;
        (void)printf(" %6d %6d", i+1, j+1);
        for(ip = 0; ip < n_potpar; ip++)
           (void)printf(" %10g",pot_par[idij].p[ip]);
        (void)putchar('\n');
      }

   if( cflag )  /* Only write if crystalline lattice */
   {
     cell[0] = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
     cell[1] = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
     cell[2] = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
     cell[3] = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/cell[1]/cell[2]);
     cell[4] = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/cell[0]/cell[2]);
     cell[5] = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/cell[0]/cell[1]);
  
     (void)printf("end\n");

     for( i = 0; i<6; i++)
         printf("%g ",cell[i]);
     for( i = 0; i<3; i++)
         printf("%g ",simbox[i]);
     puts("");

     for( i=0; i < atomtot; i++)
     {
       /* Ensure that coords are between [0,1) */
        for( j = 0; j < 3; j++)
        {
           if ( fabs(fmod(x[i][j],1)) < 5e-7 )
              x[i][j] = 0.0;
           else
              x[i][j] -= floor(x[i][j]);
        }
  
        printf("%-14s %10g %10g %10g\n",spec[site[i].pad].name,x[i][0],x[i][1],x[i][2]); 
     }
     puts("end");
   }
   return 0;
}
