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
 *  Revision Log
 *  $Log: syswrite.c,v $
 *  Revision 2.5  2002/09/18 09:59:19  kr
 *  Rolled in several changes by Craig Fisher:
 *  Ransub can now read polyatomic species
 *  Syswrite can handle polyatomics from CSSR PDB or SCHACKAL files
 *
 *  Revision 2.4  2002/06/21 11:29:07  kr
 *  Got rid of K&R varargs-compatibility stuff.
 *
 *  Revision 2.3  2001/08/09 16:41:08  keith
 *  Fixed some bugs.
 *  Converted some expressions to use array syntax.
 *
 *  Revision 2.2  2001/08/09 11:46:56  keith
 *  Tidied up against some compiler warnings.
 *  Added license file for SgInfo routines with permission of
 *  Ralf W. Grosse-Kunstleve
 *
 *  Revision 2.1  2001/08/09 09:36:36  keith
 *  Incorporated Craig's new "Syswrite" utility.
 *
 * Revision 1.1  2001/04/25  18:27:41  fisher
 * Initial revision
 *
 */
#ifndef lint
static char *RCSid = "$Header: /usr/users/kr/CVS/moldy/src/syswrite.c,v 2.5 2002/09/18 09:59:19 kr Exp $";
#endif
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
#include "sginfo.h"
#include "specdata.h"

/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
contr_mt                control;

#define xtoupper(c) (islower(c) ? toupper(c) : c)

/*========================== External data references ========================*/
extern  const pots_mt   potspec[];           /* Potential type specification  */

double  det(mat_mt a);
char	*read_ftype(char *filename);
int     read_ele(spec_data *element, char *filename);
int     read_pdb(char *, mat_mp, char (*)[32], vec_mp, double *, char *, char *);
int     read_cssr(char *, mat_mp, char (*)[32], vec_mp, double *, char *, char *);
int     read_shak(char *, mat_mp, char (*)[32], vec_mp, double *, char *, double *);
int     sgexpand(int , int , vec_mt *, char (*)[NLEN], double *, char *);
/******************************************************************************
 * add_suffix().  Add numerical suffix to string.                             *
 ******************************************************************************/
int     add_suffix(char *string, int number)
{
   char    suffix;
   int     i, digit, n_digit;
   int     length = strlen(string);

   n_digit = floor(log10(number))+1;
   if( length + n_digit > NLEN-1)       /* Longer than allowable name length */
      string[NLEN-n_digit-1] = '\0';

   for ( i = n_digit; i > 0; i--)
   {
      digit = floor(number/pow(10,i-1));
      suffix = digit + '0';             /* Calculate digit's ascii code      */
      strcat(string, &suffix);
      number = number - digit * pow(10,i-1);
   }

   return 0;
}
/******************************************************************************
 * read_pot().  Read potential data from file.                                *
 ******************************************************************************/
static
int       read_pot(char *potfile, pot_mp *pot_ptr, site_mt *spec, int max_id)
{
char       atom1[4], atom2[4];
double     chg1, chg2;
double     p_tmp;
pot_mt     pot;
int        i, n_items;
char       name[LLEN],             /* Temporary potential type name      */
           line[LLEN],             /* Store for input line from file     */
           pline[LLEN];            /* Used in pot'l paramater parsing    */
int        ptype=-1;               /* Potential type index               */
int        nerrs = 0;              /* Accumulated error count            */
int        idi, idj;
site_mt	   *spi,*spj;              /* Temporary species pointers         */

FILE       *Fpot;

    if( (Fpot = fopen(potfile,"r")) == NULL)
        return ptype;

    n_items = sscanf(get_line(line,LLEN,Fpot,1), "%s", name);
    if( n_items <= 0 )
       message(NULLI,NULLP,FATAL,SYSEOF,"potential type specification");

    for(i = 0; potspec[i].name; i++)             /* Is 'name' a known type? */
       if(strcmp(strlower(name), potspec[i].name) == 0)
          break;

    if(! potspec[i].name)                        /* Did the loop find 'name'? */
       message(&nerrs,line,FATAL,UNKPOT,name);   /* no                        */
    ptype = i;                                   /* yes                       */

    while(sscanf(get_line(line,LLEN,Fpot,1),"%s",name) > 0
                    && strcmp(strlower(name), "end") != 0)
    {
        n_items = 0;
        if(sscanf(line,"%4s %lf %4s %lf %[^#]",atom1,&chg1,atom2,&chg2,pline) <= 2)
           message(&nerrs,line,ERROR,NOPAIR);
        else
        {
                                              /* Now read in parameters */
           (void)strcat(pline, "$");              /* Add marker to end      */
           while(n_items < NPOTP && sscanf(pline,"%lf %[^#]", &p_tmp, pline) > 1 )
              pot.p[n_items++] = p_tmp;
        }

       for( idi = 0; idi < max_id; idi++)
         for( idj = idi; idj < max_id; idj++)
         {
           spi = spec+idi;
           spj = spec+idj;

           if( ((!strcmp(atom1,spi->name)) && (chg1 == spi->charge) &&
             (!strcmp(atom2,spj->name)) && (chg2 == spj->charge)) ||
                ((!strcmp(atom1,spj->name)) && (chg1 == spj->charge) &&
                   (!strcmp(atom2,spi->name)) && (chg2 == spi->charge)) )
           {
               (*pot_ptr)[idj + idi * max_id] = pot;
           }
         }
    }

   fclose(Fpot);

   if(nerrs > 0)                        /* if any errors have been detected */
      message(&nerrs,NULLP,FATAL,ERRS,nerrs,(nerrs>1)?'s':' ');

    return ptype;
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
   int		u;
   extern char	*optarg;
   int		errflg = 0;
   char		*filename = NULL;
   char         *filetype = NULL;
   char		*elename = "elements.dat";
   char		*potname = NULL;
   char		elefile[PATHLN] = "";
   char		potfile[PATHLN] = "";
   char		title[TITLE_SIZE] = "";
   double	cell[6];                      /* Cell parameters */
   mat_mt	h;                            /* Hessian */
   int		insw=-1;                      /* Switch for input format */
   int		spectot=0, atomtot=0;         /* No of species, atoms */
   int		i,j,k,idij,ip;
   int		nflag;                        /* Flag for site type matching */
   int          n_elem;                       /* No of records read from element data file */
   spec_data    element[NELEM];               /* Element info */
   spec_mt 	spec[NSPEC];                  /* Species info */
   site_mt      *site, *st, specsite[NSPEC];  /* Site info */ 
   pot_mt       *pot_par;                     /* Potential parameters */
   int		ptype = -1, n_potpar=NPOTP;
   char         spgr[16];                     /* Space Group in Herman Maugain form */
   double       simbox[3] = {1,1,1};          /* Simulation box repeat factors */
   vec_mt	x[MAX_ATOMS]; 	              /* C of M coordinates */
   double       charge[MAX_ATOMS];            /* Site charge array */
   char         label[MAX_ATOMS][NLEN];       /* Site name array */

   comm = argv[0];

   while( (u = getopt(argc, argv, "o:y:p:g:h:e:") ) != EOF )
      switch(u)
      {
       case 'g':
         if( insw > -1)
            errflg++;
	 filename = optarg;
	 insw = CSSR;
	 break;
       case 'p':
         if( insw > -1 )
            errflg++;
	 filename = optarg;
	 insw = PDB;
	 break;
       case 'h':
         if( insw > -1 )
            errflg++;
	 filename = optarg;
	 insw = SHAK;
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error("Failed to open file \"%s\" for output",optarg);
	 break;
       case 'y':
         potname = optarg;
         /* Create full path name, but don't exceed max length of string */
         strncat(strncat(potfile,POTPATH, PATHLN-strlen(potfile)), potname, PATHLN-strlen(potfile));
	 break;
       case 'e':
	 elename = optarg;
         break;
       default:
       case '?': 
	 errflg++;
      }

   if( errflg )
   {
      fputs("Usage: syswrite [-p pdb-input-file | -g cssr-input-file | -h schakal-input-file ] ",stderr);
      fputs("[-o output-file] [-y potential-parameter-file]\n",stderr);

      exit(2);
   }

   /* If no filename given on command line, request from user */
   if( insw < 0)
   {
      if( (filename = get_str("Structure file name? (PDB, CSSR or SCHAKAL)")) == NULL )
         exit(2);

      filetype = strlower(read_ftype(filename));

      if( !strcmp(filetype, "cssr") )
           insw = CSSR;
      if( !strcmp(filetype, "pdb") )
           insw = PDB;
      if( !strcmp(filetype, "shak") )
           insw = SHAK;
   }

   zero_real(h[0],9);
   zero_real(charge,MAX_ATOMS);
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
      default:
         error("Structure file \"%s\" of unknown format", filename);
   }

   if( det(h) == 0 )
      error("Invalid lattice parameters in %s", filename);

   if( strcmp(spgr,"P 1")  && atomtot > 0 )
      atomtot = sgexpand(MAX_ATOMS, atomtot, x, label, charge, spgr);

   /* Create full path name for element data, but don't exceed max length of string */
   strncat(strncat(elefile, ELEPATH, PATHLN-strlen(elefile)), elename, PATHLN-strlen(elefile));

   /* Check if element data file is available */
   if( (n_elem = read_ele(element, elefile)) < 0 )
      message(NULLI,NULLP,WARNING,NOELEM,elefile);

   site = (site_mt*)arralloc(sizeof(site_mt),1,0,atomtot-1);
   copy_atom_data(site, charge, label, atomtot);

   /* Loop to determine no of different species */
   for(st = site; st < site+atomtot; st++)
   {
      toupper(st->name[0]);   /* First character of symbol uppercase */
      strlower(st->name+1);   /* Remainder of symbol lowercase */

      nflag = 0;    /*  Flag for whether site (1) assigned or (0) not assigned to species */
      for( i = 0; i < spectot; i++)
         if( !strcmp(st->name,specsite[i].name)          /* Does site match any species identified so far? */
                 && (st->charge == specsite[i].charge ))
         {
              st->pad = i;                 /* Assign this site to species i */
              spec[i].nmols++;
              nflag++;
         }
      if (!nflag)                          /* If no matches, create new species type */
      {
         if( spectot > NSPEC )
            error("Too many species found - current limit is %d (NSPEC in specdata.h)", NSPEC);
         for( j=0; j < n_elem; j++)        /* Search element/species data for match with site */
           if( !strcmp(st->name,element[j].symbol) )
           {
              strncpy(spec[spectot].name,element[j].name,NLEN);
              st->mass = element[j].mass;
              nflag++;
              break;
           }

         if( !nflag)
             strcpy(spec[spectot].name,st->name);

         strcpy(specsite[spectot].name,st->name);
         specsite[spectot].mass = st->mass;
         specsite[spectot].charge = st->charge;
         spec[spectot].nmols=1;
         st->pad = spectot;
         spectot++;
      }
   }

   /* Modify species name if same name but different charge to another species */
   for( i = 0; i < spectot-1; i++)
   {
      k = 0;
      for(j = i+1; j < spectot; j++)
         if( !strcmp(spec[i].name,spec[j].name) )
         {
            k++;
            add_suffix(spec[j].name, k);
         }
   }

   /* Initialize potential parameter arrays */
   pot_par = aalloc(SQR(spectot), pot_mt);

   for( i = 0; i < SQR(spectot); i++)
   {
      pot_par[i].flag = 0;
      zero_real(pot_par[i].p, NPOTP);
   }

   /* Read potential parameter file if available */
   if( (potname != NULL) && ((ptype = read_pot(potfile, &pot_par, specsite, spectot)) < 0) )
       message(NULLI,NULLP,WARNING,NOPOTL,potfile);

   /* Write data to file in Moldy input form */
   if (strchr(title,'\n') != NULL )
      printf("# %s",title);
   else
      if(strlen(title) > 0)
         printf("# %s\n",title);

   printf("# System specification file written by SYSWRITE on %s\n",atime());

   for( i=0; i < spectot; i++)
   {
       printf("%-14s %d\n",spec[i].name,spec[i].nmols);
       printf("%-3d    0    0    0  %12.12g %12.12g  %s\n",i+1,
                     specsite[i].mass,specsite[i].charge,specsite[i].name);
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
   for( i = 0; i < spectot; i++)
      for( j = i; j < spectot; j++)
      {
        idij = j+i*spectot;
        (void)printf(" %6d %6d", i+1, j+1);
        for(ip = 0; ip < n_potpar; ip++)
           (void)printf(" %10g",pot_par[idij].p[ip]);
        (void)putchar('\n');
      }

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

   return 0;
}
