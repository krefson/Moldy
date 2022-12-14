/* Copyright (C) 2001 Craig Fisher
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
/************************************************************************************
 * Readers   Routines for reading SHAKAL, CSSR and PDB structure files              *
 *           Contents:                                                              *
 * is_symbol()       Return element no if label in element list                     *
 * multi_char()      Return value of run-time string as if multichar constant       *
 * read_ftype        Read filename suffix to guess file type                        *
 * read_pot          Read potential file                                            *
 * cell_to_prim      Calculate primitive cell matrix from unit cell parms           *
 * read_cssr()       Read data from Cambridge Search and Structure Retrieval file   *
 * read_pdb()        Read data from Brookhaven Protein Databank (PDB) file          *
 * read_shak()       Read data from Schakal file                                    *
 * read_shelx()      Read data from shelx file                                      *
 * read_xtl()        Read data from Biosym xtl file                                 *
 * read_xyz()        Read data from generic xyz file                                *
 * read_ele()        Read element data from file                                    *
 *****************************************************************+******************
 */
#define _XOPEN_SOURCE_EXTENDED
#include "defs.h"
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "structs.h"
#include "messages.h"
#include "utlsup.h"
#include "specdata.h"
#include "list.h"

/* For space group info */
#define SGCOREDEF__
#include "sginfo.h"

#define TITLE_SIZE  80
#define DATAREC "%d record%ssuccessfully read from %s"
#define XSATOMS "Only reading first %d atoms from \"%s\". Increase MAX_ATOMS to read more."

/*========================== External data references ========================*/
extern  const ele_data ElementData[];
extern  const pots_mt   potspec[];           /* Potential type specification  */
/* extern  const T_TabSgName TabSgName[]; */      /* Space group names, etc */

extern int transformation_matrix (char *buf, T_RTMx *trans_matrix);
extern int symm_gen (T_RTMx matrix, mat_mp apos, char (*atype)[NLEN], double *charge, int max, int natoms, int abegin, int aend);
extern void sgtransform (T_RTMx m, mat_mp x, mat_mp xp, int natoms);
extern int sgexpand (int maxnatoms, int natoms, vec_mt (*a_lst), char (*label)[NLEN], double *charge, char *spgr);

#define DEBUG_IS_ELEM 0
/****************************************************/
/* does a given string represent an element symbol? */
/****************************************************/
int is_symbol(char *input)
{
int i, j, m, n;
int num_elements;
char *label;

/* checks */
if (input == NULL)
  return(0);
/* zero length string matches something, so force it to return 0 */
if (!strlen(input))
  return(0);

/* duplicate for manipulation */
label = mystrdup(input);

/* remove anything but alphabetic chars */
for (i=0 ; i<strlen(label) ; i++)
  if (!isalpha(*(label+i)))
    *(label+i) = ' ';
trim(label);

m = strlen(label);

/* catch Deuterium */
if (strncasecmp(label, "D", strlen(label)) == 0)
  *label = 'H';

#if DEBUG_IS_ELEM
printf("Looking for [%s]...", label);
#endif

/* attempt to match atom type with database */
j=0;
num_elements = sizeof *ElementData;

/* FIXME - elim dependence on const */
for(i=1 ; i<num_elements ; i++)
  {
/* only compare if lengths match */
  n = strlen(ElementData[i].symbol);
  if (n == m)
    {
    if (strcasecmp(label, ElementData[i].symbol) == 0)
      {
      j = i;
      break;
      }
    }
  }

#if DEBUG_IS_ELEM
if (j)
  fprintf(stderr,"found.\n");
else
  fprintf(stderr,"not found.\n");
#endif
return(j);
}
/******************************************************************************
 * multi_char().  Return the value of a (run-time) string as if it had been   *
 *                declared as a multi-character constant. Up to 4-byte ints.  *
 ******************************************************************************/
int multi_char(char *s)
{
   static int unity2[] = {'\1\0','\0\1'},
              unity3[] = {'\1\0\0','\0\1\0','\0\0\1'},
              unity4[] = {'\1\0\0\0','\0\1\0\0','\0\0\1\0','\0\0\0\1'};
   int *unity;
   int m=0, l=MIN(strlen(s),sizeof *unity);
   switch (l)
   {
   default:
      return 0;
   case 1:
      return *s;
   case 2:
      unity = unity2;
      break;
   case 3:
      unity = unity3;
      break;
   case 4:
      unity = unity4;
      break;
   }
   while( l-- > 0 )
      m |= *unity++ * *s++;
   return m;
}
/******************************************************************************
 * cell_to_prim().  Convert unit cell params to primitive vector matrix.      *
 ******************************************************************************/
void cell_to_prim(double *cell, mat_mt h)
{
double   ca, cb, cg, sg;               /* Cos and sin of cell angles */

   ca = cos(cell[3]*DTOR); cb = cos(cell[4]*DTOR);
   cg = cos(cell[5]*DTOR); sg = sin(cell[5]*DTOR);

   h[0][0] = cell[0];              /* Set up primitive cell matrix   */
   h[0][1] = cell[1] * cg;         /* from lengths and angles.       */
   h[1][1] = cell[1] * sg;
   h[0][2] = cell[2] * cb;
   h[1][2] = cell[2] / sg * (ca - cb*cg);
   h[2][2] = cell[2] / sg * sqrt(1 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg);
   h[1][0] = h[2][0] = h[2][1] = 0.0;
}
/******************************************************************************
 * read_ftype().  Read suffix from filename.                                  *
 ******************************************************************************/
char    *read_ftype(char *filename)
{
char    *s = filename+strlen(filename);

   if( strstr(filename,".") != NULL)
   {
     while( s != NULL && *s != '.' )
       s--;

     if( *s == '.')
       s++;
     return (s);
   }
   else
     return ("");
}
/******************************************************************************
 * read_pot().  Read potential data from file.                                *
 ******************************************************************************/
int       read_pot(char *potfile, pot_mp *pot_ptr, site_mt *spec, int add_sites, int max_id)
{
char       atom1[4], atom2[4];
double     chg1, chg2;
double     p_tmp;
pot_mt     pot;
int        i, m, n, n_items;
char       name[LLEN],             /* Temporary potential type name      */
           line[LLEN],             /* Store for input line from file     */
           pline[LLEN];            /* Used in pot'l paramater parsing    */
int        ptype=-1;               /* Potential type index               */
int        nerrs = 0;              /* Accumulated error count            */
int        idi, idj;
site_mt    *spi,*spj;              /* Temporary species pointers         */
boolean    flag1, flag2;           /* Atom labels = element symbol flags */
int        match1, match2;         /* Atoms match sites */
int        code1, code2;           /* Element number */

FILE       *Fpot;

    if( (Fpot = fopen(potfile,"r")) == NULL)
        return ptype;

    n_items = sscanf(get_line(line,LLEN,Fpot,1), "%s", name);
    if( n_items <= 0 )
       message(NULLI,NULLP,FATAL,SYSEOF,"potential type specification");

    for(i = 0; potspec[i].name; i++)             /* Is 'name' a known type? */
       if(strcasecmp(name, potspec[i].name) == 0)
          break;

    if(! potspec[i].name)                        /* Did the loop find 'name'? */       message(&nerrs,line,FATAL,UNKPOT,name);   /* no                        */    ptype = i;                                   /* yes                       */
    while(sscanf(get_line(line,LLEN,Fpot,1),"%s",name) > 0
                    && strcasecmp(name, "end") != 0)
    {
        n_items = 0;
        flag1 = flag2 = false;
        match1 = match2 = 0;

        if(sscanf(line,"%4s %lf %4s %lf %[^#]",atom1,&chg1,atom2,&chg2,pline) <= 2)
           message(&nerrs,line,ERROR,NOPAIR);
        else
        {
                                              /* Now read in parameters */
           (void)strcat(pline, "$");              /* Add marker to end      */
           while(n_items < NPOTP && sscanf(pline,"%lf %[^#]", &p_tmp, pline) > 1 )
              pot.p[n_items++] = p_tmp;
        }

       /* Is atom1 an element symbol or other label? */
       m = is_symbol(atom1);
       if (m)
         if( strncasecmp(atom1, ElementData[m].symbol, strlen(ElementData[m].symbol)))
           flag1 = true;

       /* Is atom2 an element symbol or other label? */
       n = is_symbol(atom2);
       if (n)
         if( strncasecmp(atom2, ElementData[n].symbol, strlen(ElementData[n].symbol)))
           flag2 = true;

       for( idi = 0; idi < max_id; idi++)
       {
         spi = spec+idi;
         for( idj = idi; idj < max_id+add_sites; idj++)
         {
           spj = spec+idj;

           code1 = is_symbol(spi->name);
           code2 = is_symbol(spj->name);

           if(flag1)
             {
             if( (m == code1) && (chg1 == spi->charge))
               match1 = 1;
             else
               if( (m == code2) && (chg1 == spj->charge))
                 match1 = 2;
               else
                 match1 = 0;
             }
           else
             {
             if( (!strcmp(atom1,spi->name)) && (chg1 == spi->charge))
               match1 = 1;
             else
               if( (!strcmp(atom1,spj->name)) && (chg1 == spj->charge))
                 match1 = 2;
               else
                 match1 = 0;
             }

           if(flag2)
             {
             if( (n == code2) && (chg2 == spj->charge))
               match2 = 2;
             else
               if( (n == code1) && (chg2 == spi->charge))
                 match2 = 1;
               else
                 match2 = 0;
             }
           else
             {
             if( (!strcmp(atom2,spj->name)) && (chg2 == spj->charge))
               match2 = 2;
             else
               if( (!strcmp(atom2,spi->name)) && (chg2 == spi->charge))
                 match2 = 1;
               else
                 match2 = 0;
             }

           if( (match1 == 1 && match2 == 2) ||
               (match1 == 2 && match2 == 1) )
           {
               /* Assign additional potential parameters */
               (*pot_ptr)[idi + idj * max_id] = pot;
               (*pot_ptr)[idj + idi * max_id] = pot;
           }
         }
       }
    }

   fclose(Fpot);

   if(nerrs > 0)                        /* if any errors have been detected */
      message(&nerrs,NULLP,FATAL,ERRS,nerrs,(nerrs>1)?'s':' ');

    return ptype;
}
/******************************************************************************
 * read_cssr().  Read structural data from cssr file.                         *
 ******************************************************************************/
int      read_cssr(char *filename, mat_mp h, char (*label)[NLEN], vec_mp x, double *charge, char *title, char *spgr)
{
int      i;		  /* Counter */
int	 coord_sys=-1;	  /* -1 = unknown, 0 = frac coords, 1 = Cartesian coords */
int	 nocell=0;        /* 0 = no cell, 1 = cell parms given */
int      natoms;          /* Number of atoms in structure */
double   cell[6];         /* Cell parameters */
mat_mt   hinv;            /* Inverse h matrix */
char     temp_name[6];    /* Temporary storage for atom label */
char     line[LLEN];      /* Storage for line read from file */
FILE     *Fp;             /* cssr file pointer */
int	 sgno=0, sgopt=0; /* Space group number and associated option */

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   /* Read cssr file header */
   if( sscanf(get_line(line,LLEN,Fp,0),"%*38c%8lf%8lf%8lf",&cell[0],&cell[1],&cell[2]) < 3)
      nocell++;

   if( sscanf(get_line(line,LLEN,Fp,0),"%*21c%8lf%8lf%8lf    SPGR =%3d %11[a-zA-Z0-9-/ ] OPT =%2d",
          &cell[3],&cell[4],&cell[5],&sgno,spgr,&sgopt) < 5)
      nocell++;

   if( sscanf(get_line(line,LLEN,Fp,0),"%4d%4d %60c", &natoms, &coord_sys, title) < 2)
            error("EOF or unexpected format on line 3 in \"%s\"", filename);

   if( natoms > MAX_CSSR )
      error("\"%s\" contains too many atoms! (max: %d)\n", filename, MAX_CSSR);

   if( !coord_sys )
   {
      if( nocell )
         error("No unit cell parameters supplied for fractional coordinates");

      if( ! (cell[0] > 0 && cell[1] > 0 && cell[2] > 0 &&
             cell[3] > 0 && cell[3] < 180.0 &&
             cell[4] > 0 && cell[4] < 180.0 &&
             cell[5] > 0 && cell[5] < 180.0) )
                 error("Invalid unit cell in \"%s\"", filename);
   }

   trim(spgr);
   if( sgno > 1 && (strcmp(spgr,"") == 0 || strcmp(spgr,"P 1") == 0))
      sprintf(spgr,"%d",sgno);

   /* Match cssr options to axis choice */
   if(sgopt)
   {
     if( sgno == 5 || sgno == 7 || sgno == 8 ||
         sgno == 12 || sgno == 13 || sgno == 14 )
     {
        if( sgopt < 4)
          strcat(spgr,":b");
        else if( sgopt > 3 && sgopt < 7)
          strcat(spgr,":c");
        else if( sgopt > 6 && sgopt < 10)
          strcat(spgr,":a");
        sgopt = sgopt % 3;
        if(!sgopt) sgopt = 3;
          sprintf(spgr,"%s%d",spgr,sgopt);
     }
     /* Match cssr options to sginfo list order */
     else if( sgno == 9 || sgno == 15 )
     {
        if( sgopt < 4)
          strcat(spgr,":b");
        else if( sgopt > 3 && sgopt < 7)
          strcat(spgr,":-b");
        else if( sgopt > 6 && sgopt < 10)
          strcat(spgr,":c");
        else if( sgopt > 9 && sgopt < 13)
          strcat(spgr,":-c");
        else if( sgopt > 12 && sgopt < 16)
          strcat(spgr,":a");
        else if( sgopt > 15 && sgopt < 19)
          strcat(spgr,":-a");
        sgopt = sgopt % 3;
        if(!sgopt) sgopt = 3;
          sprintf(spgr,"%s%d",spgr,sgopt);
     }
     else if( sgno == 50 || sgno == 59 || sgno == 68)
     {
        sgopt = sgopt % 2;
        if(!sgopt)
          strcat(spgr,":2");
     }
     else if( sgno == 146 || sgno == 148 || sgno == 155 ||
              sgno == 160 || sgno == 161 || sgno == 166 ||
              sgno == 167)
     {
         if( sgopt == 2)
           strcat(spgr,":r");
         else
           strcat(spgr,":h");
     }
     /* Range of space groups with possible alternative origins */
     else if( sgno == 70 || (sgno > 84 && sgno < 143) ||
              sgno > 201)
        sprintf(spgr,"%s:%d",spgr,sgopt);
  }

   if( fgets(line,sizeof(line),Fp) == NULL )
      error("Unexpected end of input on line 4 of \"%s\"", filename);

   for( i = 0; i < natoms; i++)
   {
      if( sscanf(get_line(line,LLEN,Fp,0),"%*4d %4s  %9lf %9lf %9lf %*4d%*4d%*4d%*4d%*4d%*4d%*4d%*4d %7lf",
               temp_name, &x[i][0], &x[i][1], &x[i][2], &charge[i]) < 5 )
          error("EOF or unexpected format in atom definitions in \"%s\"", filename);

      str_cut(temp_name, label[i]); /* Convert cssr atom name to atom symbol */
      trim(label[i]);
   }
   fclose(Fp);

   /*
    * Convert to fractional co-ordinates if in real coords
    */
   if( !nocell )
   {
      cell_to_prim(cell,h);
      if( coord_sys > 0 )
      {
         invert(h,hinv);
         mat_vec_mul(hinv, x, x, natoms);   /* Convert to fractional coords */
      }
   }

   return natoms;
}
/******************************************************************************
 * read_pdb().  Read structural data from pdb file.                           *
 ******************************************************************************/
int      read_pdb(char *filename, mat_mp h, char (*label)[NLEN], vec_mp x, double *charge, char *title, char *spgr)
{
int      natoms = 0;
double   cell[6];                                  /* Cell parameters */
char     dummy[4], chg[2];
mat_mt   hinv;
double   tr[3];
char     line[LLEN], keyword[LLEN];
int      irow;
int  scaleflg = 0, crystflag = 0;
FILE     *Fp;

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   while( strcasecmp(keyword,"end") !=0 )
   {
       if( sscanf(get_line(line,LLEN,Fp,0),"%s",keyword) < 1)
           error("File \"%s\" has incorrect format", filename);

       if(strcasecmp(keyword,"title") == 0 )         /* Title */
       {
          if(isspace(line[9]))
             strncpy(title,line+10,TITLE_SIZE);
          else {
             strncat(title," ",TITLE_SIZE);
             strncat(title,line+10,TITLE_SIZE);
          }
       }

       if(strcasecmp(keyword,"cryst1") == 0)   /* CRYST1 - Specify unit cell */
       {
          if( sscanf(line,"CRYST1%9lf%9lf%9lf%7lf%7lf%7lf%11c",
             &cell[0],&cell[1],&cell[2],&cell[3],&cell[4],&cell[5],spgr) < 6 )
                error("Error in CRYST1 line of \"%s\" -- should have at least 6 parameters", filename);

          trim(spgr);
          crystflag++; /* Unit cell defined for this structure */
       }

       if(strncasecmp(keyword,"scale",5) == 0 )         /* SCALE */
       {
          sscanf(line, "SCALE%1d", &irow);
          if( irow < 1 || irow > 3 )
             error("Error in \"%s\" - Error in SCALE[1-3] - %6s unknown", filename, line);

          sscanf(line+10,"%10lf%10lf%10lf     %10lf",
                 &hinv[irow-1][0], &hinv[irow-1][1], &hinv[irow-1][2], &tr[irow-1]);

          scaleflg+=irow;
       }

       if(strcasecmp(keyword,"hetatm") == 0 || strcasecmp(keyword,"atom") == 0)
       {
          if( natoms >= MAX_ATOMS )
          {
             message(NULLI, NULLP, WARNING, XSATOMS, MAX_ATOMS, filename);
             natoms = MAX_ATOMS;
             break;
          }

          sscanf(line+12,"%4s%*1c%*3s%*2c%*4d%*4c%8lf%8lf%8lf", dummy,
               &x[natoms][0], &x[natoms][1], &x[natoms][2]);
          str_cut(dummy, label[natoms]);

          trim(label[natoms]);
          sprintf(chg,"  ");
          sscanf(line+78,"%2s", chg);
          charge[natoms] = str_cut(chg, dummy);

          natoms++;
       }
   }
   fclose(Fp);

   if( crystflag ) /* Only transform coords if unit cell defined */
   {
       if( ! (cell[0] > 0 && cell[1] > 0 && cell[2] > 0 &&
          cell[3] > 0 && cell[3] < 180.0 &&
          cell[4] > 0 && cell[4] < 180.0 &&
          cell[5] > 0 && cell[5] < 180.0))
      {
         error("Error in \"%s\" - Invalid unit cell %f %f %f %f %f %f", filename,
                       cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
      }
      /*
       * Convert to fractional co-ordinates
       */
      cell_to_prim(cell,h);
      if( scaleflg != 6 )
         invert(h,hinv);

      mat_vec_mul(hinv, x, x, natoms);   /* Convert to fractional coords */
   }

   return natoms;
}
/******************************************************************************
 * read_shak().  Read structural data from SCHAKAL file.                      *
 ******************************************************************************/
int      read_shak(char *filename, mat_mt h, char (*label)[NLEN], vec_mp x, double *charge, char *title, double *simbox)
{
int      i, natoms = 0;
double   cell[6];                                    /* Cell parameters */
double   box[6];
char     line[LLEN], *buf, *bufend, last;
int      linec, linep; char *linev[32];
int      blankflg = 0;                    /* Blanks in atom labels?            */
int      multi_flag = 0;                  /* Flag for multiple units in file   */
int      dupl_flag, symm_unit_reset=1;
int      symm_unit_begin=0, symm_unit_end=0; /* Initialization for lint only */
T_RTMx   trans_matrix;
double   sfact;
int      keywd = 0;
FILE     *Fp;

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   while( !feof(Fp) )
   {
      if( fgets(line, sizeof line, Fp) == NULL)
         break;              /* Force end on EOF     */

      bufend = strchr(line,'\n');
      if( bufend == NULL)
         bufend = line+strlen(line);
      *bufend = '\0';
      buf = line;
      last = '\0';

      while(buf < bufend )
      {
         dupl_flag=0;
         linec = get_tokens(buf, linev, " ");
         if( linec < 1 )
            break;
         linep = 0;
         if( last != ';' )
         {
            (void) strlower(linev[0]);
            linev[0][2] = '\0';
            keywd = multi_char(linev[0]);
            if( linec == 1)
               buf += strlen(buf);
            else
               buf += linev[1]-linev[0];
            linec--;
            linep++;
         }
         switch( keywd )
         {
          case 0:               /* Blank line                   */
          case ':':
             break;
          case 'ti':    /* Title                                */
             strcpy(title,buf);
             break;
          case 'li':    /* List  - No action                    */
             break;
          case 'ce':    /* Cell - Specify unit cell             */
             if( multi_flag )   /* Only first dataset can specify cell  */
                break;
             i = 0;
             cell[0] = cell[1] = cell[2] = 1.0;
             cell[3] = cell[4] = cell[5] = 90.0;
             switch(linec)
             {
              case 0:
                break;
              case 1:                    /* Cubic - only "a" specified   */
                cell[0] = cell[1] = cell[2] = atof(linev[linep]);
                break;
              case 3:                    /* Orthorhombic -- a,b,c given  */
              case 6:                    /* Full 6 parameters given      */
                while(linep <= linec)
                   cell[i++] = atof(linev[linep++]);
                break;
              default:
                error( "Error in \'%s\'- CELL should have 0, 1, 3 or 6 parameters",filename);
                return 0;
             }

             if( ! (cell[0] > 0 && cell[1] > 0 && cell[2] > 0 &&
                    cell[3] > 0 && cell[3] < 180.0 &&
                    cell[4] > 0 && cell[4] < 180.0 &&
                    cell[5] > 0 && cell[5] < 180.0))
             {
                error("Invalid unit cell %f %f %f %f %f %f in \'%s\'",
                        cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], filename);
                return 0;
             }
             cell_to_prim(cell,h);
             break;
          case 'at':    /* Atom - add atom to list              */
             if( natoms >= MAX_ATOMS )
             {
                message(NULLI, NULLP, WARNING, XSATOMS, MAX_ATOMS, filename);
                natoms = MAX_ATOMS;
                break;
             }

             if( linec < 4 + blankflg)
             {
                error("Insufficient items on ATOM card in %s - %s", filename, buf);
                return 0;
             }
             if( blankflg )
             {
                strncpy(label[natoms],strcat(linev[linep],linev[linep+1]),NLEN);
                linep += 2;
             }
             else
                strncpy(label[natoms], linev[linep++],NLEN);

             label[natoms][NLEN-1] = '\0';

            for(i = 0; i < 3; i++)
               x[natoms][i] = atof(linev[linep++]);

            if( symm_unit_reset )
            {
               symm_unit_reset = 0;
               symm_unit_begin = natoms;
            }
            trim(label[natoms]);
            natoms++;
            symm_unit_end = natoms;
            break;
          case 'as':    /* Assume Blanks - modifies parsing     */
             if( blankflg == 0 )
                blankflg++;
            break;
          case 'du':    /* Duplicate symmetry                   */
             dupl_flag++;
          case 'sy':    /* Duplicate explicit atoms only        */
             if( ! transformation_matrix(buf, &trans_matrix))
                 error("Syntax Error on SYMM or DUPL card in %s\n%s\n", filename, buf);

             natoms = symm_gen(trans_matrix, x, label, charge,
                               MAX_ATOMS,natoms,symm_unit_begin,
                               dupl_flag?natoms:symm_unit_end);
             symm_unit_reset++;
             break;
          case 'tr':    /* Translate atoms.                     */
             if( ! transformation_matrix(buf, &trans_matrix))
                 error("Syntax Error on TRANS card in %s\n%s\n", filename, buf);
             sgtransform(trans_matrix, x, x, natoms);
             break;
          case 'ad':    /* For adding corners of unit cell      */
             break;
          case 'bo':    /* Box size specifiers                  */
             switch(linec)       /* Various numbers of parameters allowed  */
             {
              case 6:            /* All 6 lower and upper bounds         */
                linep = 1;
                for(i = 0; i < 3; i++)
                {
                   sfact = atof(linev[linep++]);
                   box[i] = sfact;
                   sfact = atof(linev[linep++]);
                   box[i+3] = sfact;
                }
                break;
              case 3:            /* Expansion in all 3 directions        */
                linep = 1;
                for(i = 0; i < 3; i++)
                {
                   sfact = atof(linev[linep++]);
                   box[i] = -sfact;
                   box[i+3] = 1+sfact;
                }

                break;
              case 1:            /* General expansion factor.            */
                sfact = atof(linev[1]);

                box[0] = box[1] = box[2] = -sfact;
                box[3] = box[4] = box[5] = 1+sfact;
                break;
              case 0:            /* No expansion factor                  */
                box[0] = box[1] = box[2] = 0.0;
                box[3] = box[4] = box[5] = 1.0;
                break;
              default:
                error("Incorrect number of parameters (%d) for BOX label in %s", linec-1, filename);
             }
             for( i = 0; i < 3; i++)
                simbox[i] = box[i+3] - box[i];
             symm_unit_reset++;
             break;
          case 'pa':    /* Pack -                                  */
             symm_unit_reset++;
             error("Unsupported label %s in %s", linev[0], filename);
          case 'mo':    /* Mol separate molecules for PACK.        */
             break;
          case 'gr':    /* Group - sets group index for atoms      */
             break;
          case 'co':    /* CONN - no action                        */
             break;
          case 'zp':    /* ZPAR - no action                        */
             break;
          case 'en':    /* End - Add unit cell, etc                */
             if( natoms == 0  )
                return 0;

             if( linec > 0 )     /* End of independent molecule -          */
             {                   /*     more molecules follow so continue  */
                symm_unit_reset++;
                multi_flag++;
                break;
             }

             return natoms;      /* Finished reading file   */
          default:
             if(! ispunct(linev[0][0] ) )
                error("Unknown label %s\nbuffer=\"%s\"in %s",linev[0], buf, filename);
             break;
         }

        /*
         *  Is there another command on this line?
         */
         while( buf < bufend && *buf != ':' && *buf != ';' )
            buf++;

         if( buf < bufend )             /* Store separator as flag.          */
            last = *buf++;
      }
   }
   return natoms;                       /* Finished reading file   */
}
/******************************************************************************
 * read_xtl().  Read structural data from xtl file.                           *
 ******************************************************************************/
int      read_xtl(char *filename, mat_mp h, char (*label)[NLEN], vec_mp x, double *charge, char *title, char *spgr)
{
int      i, num_tokens, natoms = 0;
int	 shift, num_columns;
double   cell[6];                                  /* Cell parameters */
mat_mt   hinv;
char     line[LLEN], *buffer[LLEN];
char	 *temp_name=NULL;
int      sgno = 0;
int	 dimension = 3;
boolean	 crystflag = false, cart = false;
int	 column[9];
FILE     *Fp;

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

  /* initialize columns */
   for(i=9;i--;)
     column[i] = -1;

   while( !feof(Fp) )
   {
      get_line(line,LLEN,Fp,0);
      if(strncasecmp(line,"tit",3) == 0 )         /* Title */
      {
         if( strlen(line+5) == 0)
            strncat(title, filename, TITLE_SIZE);
         else
            strncat(title, line+5, TITLE_SIZE);
         trim(title);
      }

      if (strncasecmp(line,"dim ",3) == 0)  /* Periodicity */
      {
         if( sscanf(strlower(line),"dimension %d", &dimension) != 1)
            error("Error in DIMENSION line of \"%s\" -- should have 1 parameter", filename);
      }
      if( dimension != 3)
            error("Error - \"%s\" does not contain a 3D structure", filename);

      if(strncasecmp(line,"cel",3) == 0)   /* CELL - Specify unit cell */
      {
         if( sscanf(get_line(line,LLEN,Fp,0),"%10lf %10lf %10lf %10lf %10lf %10lf",
            &cell[0],&cell[1],&cell[2],&cell[3],&cell[4],&cell[5]) < 6 )
               error("Error in CELL line of \"%s\" -- should have 6 parameters", filename);
         crystflag = true;
      }

      if(strncasecmp(line,"sym ",3) == 0 )         /* Space group */
      {
        num_tokens = get_tokens(line, buffer, " ");
        for (i=1 ; i<num_tokens-1 ; i++)
        {
          /* seek space group number */
          if (strncasecmp("num", *(buffer+i),3) == 0)
            sscanf(*(buffer+i+1),"%d", &sgno);
          /* seek space group label */
          if (strncasecmp("lab", *(buffer+i),3) == 0)
            strcpy(spgr, *(buffer+i+1));

          if( sgno > 1 && (strcmp(spgr,"") == 0 || strcmp(spgr,"P 1") == 0))
            sprintf(spgr,"%d",sgno);

          if (strncasecmp("qua", *(buffer+i),3) == 0)
          {
            strcat(spgr,":");

            if(strcasecmp("origin", *(buffer+i+1)) == 0)
              strcat(spgr,*(buffer+i+2));

            if( strcasecmp(*(buffer+i),"hexagonal") == 0 )
                 strcat(spgr,"h");

            if( strcasecmp(*(buffer+i),"rhombohedral") == 0 )
                 strcat(spgr,"r");
          }
        }
      }
      /* Read in atom data */
      /* nb. columns could be in any order */
      if(strncasecmp(line,"ato",3) == 0)
      {
         get_line(line,LLEN,Fp,0);
         num_tokens = get_tokens(line, buffer, " ");
         /* Determine order of columns */
         for(i=0; i<num_tokens; i++)
         {
           /* NAME secondary keyword */
           if (strncasecmp("nam", *(buffer+i), 3) == 0)
              column[0] = i;
           /* X secondary keyword */
           if (strncasecmp("x", *(buffer+i), 1) == 0)
              column[1] = i;
           /* Y secondary keyword */
           if (strncasecmp("y", *(buffer+i), 1) == 0)
              column[2] = i;
           /* Z secondary keyword */
           if (strncasecmp("z", *(buffer+i), 1) == 0)
              column[3] = i;
           /* SCATTER secondary keyword */
           if (strncasecmp("sca", *(buffer+i), 3) == 0)
              column[4] = i;
           /* CHARGE secondary keyword */
           if (strncasecmp("cha", *(buffer+i), 3) == 0)
              column[5] = i;
           /* CARX secondary keyword */
           if (strncasecmp("carx", *(buffer+i), 4) == 0)
           {
              column[6] = i;
              cart = true;
           }
           /* CARY secondary keyword */
           if (strncasecmp("cary", *(buffer+i), 4) == 0)
           {
              column[7] = i;
              cart = true;
           }
           /* CARZ secondary keyword */
           if (strncasecmp("carz", *(buffer+i), 4) == 0)
           {
              column[8] = i;
              cart = true;
           }
         }
         for (;;)
         {
           get_line(line,LLEN,Fp,0);
           num_columns = get_tokens(line, buffer, " ");
           if(!buffer)
             break;

           if (strncasecmp("eof", *buffer, 3) != 0)
           {
             if( natoms >= MAX_ATOMS )
             {
                message(NULLI, NULLP, WARNING, XSATOMS, MAX_ATOMS, filename);
                natoms = MAX_ATOMS;
                break;
             }

/* enough tokens */
/* Use scattering element in lookup, as label not always element symbol */
/* otherwise check if first item is a valid atom type */
/* C. Fisher 2004 */
             if (column[4] != -1)
             {
               temp_name = mystrdup(*(buffer+column[4]));
               if( num_columns > num_tokens)
                 strcat(temp_name, *(buffer+column[4]+1));
             }
             else
               if(column[0] != -1)
                 temp_name = mystrdup(*(buffer+column[0]));
               else
                 break;
             if( temp_name != NULL)
             {
               str_cut(temp_name, label[natoms]); /* Extract atom symbol from symbol+number */
               trim(label[natoms]);
             }

             if( num_columns > num_tokens) /* SCATTER factor can exist as two items */
               shift=1;      /* Shift column by one to allow for extra scatter term */
             else
               shift=0;

             if( !cart )
             {
               if( column[1] > column[4] && shift==1 )
                 x[natoms][0] = atof(*(buffer+column[1]+shift));
               else
                 x[natoms][0] = atof(*(buffer+column[1]));
               if( column[2] > column[4] && shift==1 )
                 x[natoms][1] = atof(*(buffer+column[2]+shift));
               else
                 x[natoms][1] = atof(*(buffer+column[2]));
               if( column[3] > column[4] && shift==1 )
                 x[natoms][2] = atof(*(buffer+column[3]+shift));
               else
                 x[natoms][2] = atof(*(buffer+column[3]));
             }
             else
             {
               if( column[6] > column[4] && shift==1 )
                 x[natoms][0] = atof(*(buffer+column[6]+shift));
               else
                 x[natoms][0] = atof(*(buffer+column[6]));
               if( column[7] > column[4] && shift==1 )
                 x[natoms][1] = atof(*(buffer+column[7]+shift));
               else
                 x[natoms][1] = atof(*(buffer+column[7]));
               if( column[8] > column[4] && shift==1 )
                 x[natoms][2] = atof(*(buffer+column[8]+shift));
               else
                 x[natoms][2] = atof(*(buffer+column[8]));
             }
             if( column[5] != -1)
             {
               if( column[5] > column[4] && shift==1 )
                 charge[natoms] = atof(*(buffer+column[5]+shift));
               else
                 charge[natoms] = atof(*(buffer+column[5]));
             }
             natoms++;
           } 
           else
             break;
        }
      }
   }
   fclose(Fp);

   if( crystflag ) /* Only transform coords if unit cell defined */
   {
       if( ! (cell[0] > 0 && cell[1] > 0 && cell[2] > 0 &&
          cell[3] > 0 && cell[3] < 180.0 &&
          cell[4] > 0 && cell[4] < 180.0 &&
          cell[5] > 0 && cell[5] < 180.0))
      {
         error("Error in \"%s\" - Invalid unit cell %f %f %f %f %f %f", filename,
                       cell[0], cell[1], cell[2], cell[3], cell[4], cell[5]);
      }

      cell_to_prim(cell,h);
      if( cart)
      {
        invert(h,hinv);
        mat_vec_mul(hinv, x, x, natoms);   /* Convert to fractional coords */
      }
   }

   return natoms;
}
/******************************************************************************
 * read_xyz().  Read structural data from xyz file.                           *
 ******************************************************************************/
int      read_xyz(char *filename, mat_mp h, char (*label)[NLEN], vec_mp x, char *title)
{
int      natoms = 0;
char     line[LLEN];
char     temp_name[8];
int	 i;
FILE     *Fp;

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   if( sscanf(get_line(line,LLEN,Fp,0),"%d", &natoms) < 1)
      error("EOF or unexpected format on line 1 in \"%s\"", filename);

   if( natoms >= MAX_ATOMS )
   {
      message(NULLI, NULLP, WARNING, XSATOMS, MAX_ATOMS, filename);
      natoms = MAX_ATOMS;
   }

   if( sscanf(get_line(line,LLEN,Fp,0),"%s", title) < 1)
      error("EOF or unexpected format on line 2 in \"%s\"", filename);

   strncpy(title,line,TITLE_SIZE);
   trim(title);
   if( strlen(title) == 0)
      strncat(title,filename,TITLE_SIZE);

   for( i = 0; i < natoms; i++)
   {
      if( sscanf(get_line(line,LLEN,Fp,0),"%s %lf %lf %lf",
            temp_name, &x[i][0], &x[i][1], &x[i][2]) < 4 )
         error("EOF or unexpected format in atom definitions in \"%s\"", filename);

      str_cut(temp_name, label[i]); /* Convert atom name to atom symbol */
      trim(label[i]);
   }

   fclose(Fp);

   return natoms;
}
/******************************************************************************
 * read_shelx(). Read structural data from shelx .ins or .res file.           *
 ******************************************************************************/
int      read_shelx(char *filename, mat_mp h, char (*label)[NLEN], vec_mp x, char *title, char *spgr)
{
int       i, num_tokens, natoms = 0, n_elem;
const ele_data *elem;
double    cell[6];                                  /* Cell parameters */
char      line[LLEN], *buffer[LLEN];
FILE      *Fp;
T_SgInfo  SgInfo;
const T_LatticeInfo  *LatticeInfo = LI_P; /* Defaults to primitive */
int       Latt_N = 1;   /* Defaults to centrosymmetric primitive */
char      *xyz = NULL; /* Symmetry operators */
ROOT      *xyz_list = NULL; /* Linked list for xyz symmetry operator list */
NODE      *node;    /* Node for xyz symmetry operator list */
T_RTMx    SeitzMx;
char      *sglabel;
char      *temp_name = NULL;
char      dummy[4];

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   while ( get_line(line,LLEN,Fp,0))
   {
       num_tokens = get_tokens(line, buffer, " ");
       if( num_tokens < 1)
           error("File \"%s\" has incorrect format", filename);

       if(strncasecmp(*buffer,"title",5) == 0 )         /* Title */
          for( i=0; i<num_tokens; i++)
             strncat(title,*(buffer+i),TITLE_SIZE);

       if(strncasecmp(*buffer,"cell",4) == 0)   /* CELL - Specify unit cell */
       {
         if (num_tokens > 2)
         {
           cell[0] = atof(*(buffer+2));
           cell[1] = atof(*(buffer+3));
           cell[2] = atof(*(buffer+4));
           cell[3] = atof(*(buffer+5));
           cell[4] = atof(*(buffer+6));
           cell[5] = atof(*(buffer+7));
         }
         else
           error("Error in CELL line of \"%s\" -- should have at least 6 parameters", filename);

         cell_to_prim(cell,h);
       }

       if(strncasecmp(*buffer,"latt",4) == 0 )
       {
          Latt_N = atof(*(buffer+1));

          switch (abs(Latt_N))
          {
          case 1: LatticeInfo = LI_P; break;
          case 5: LatticeInfo = LI_A; break;
          case 6: LatticeInfo = LI_B; break;
          case 7: LatticeInfo = LI_C; break;
          case 2: LatticeInfo = LI_I; break;
          case 3: LatticeInfo = LI_R; break;
          case 4: LatticeInfo = LI_F; break;
          }
       }

     if(strncasecmp(*buffer,"symm",4) == 0 )
     {
       /* Read xyz symmetry operators */
       xyz = mystrdup(*(buffer+1));
       for(i=2; i<num_tokens; i++)
         xyz = strcat(xyz, strcat(" ",*(buffer+i)));

       if( insert_data(&xyz_list, strlower(xyz), 1) < 0) /* Add to operator list */
          error("Error creating first node in sym list - \n%s\n",strerror(errno));
     }

     /* End of useful header info */
     if( strncasecmp("sfac", *buffer, 4) == 0 ||
        strncasecmp("unit", *buffer, 4) == 0)
     break;
   }

/* atom coordinate search */
   for (;;)
   {
     get_line(line,LLEN,Fp,0);
     num_tokens = get_tokens(line, buffer, " ");

     if (!buffer)
       break;

     if(strncasecmp(*buffer,"end",3) == 0 ||
        strncasecmp(*buffer,"hklf",4) == 0)
       break;

   /* Check atom is real element */

     sprintf(dummy,"    ");
     str_cut(mystrdup(*buffer), dummy); /* Convert atom label to atom symbol */
     trim(dummy);
     for( elem=ElementData; elem; elem++)
       if( !strcmp(dummy,elem->symbol) &&
           n_elem > 0 && num_tokens > 4 &&
           strcmp(".",*(buffer+2)) != 0.0 &&
           strcmp(".",*(buffer+3)) != 0.0 &&
           strcmp(".",*(buffer+4)) != 0.0 )
         {
           strcpy(label[natoms],elem->symbol);
           trim(label[natoms]);
           x[natoms][0] = atof(*(buffer+2));
           x[natoms][1] = atof(*(buffer+3));
           x[natoms][2] = atof(*(buffer+4));
           natoms++;
           break;
         }
   }
   InitSeitzMx(&SeitzMx, 1);
   InitSgInfo(&SgInfo);

   SgInfo.MaxList = 192;  /* absolute maximum number of symops */

   SgInfo.ListSeitzMx
      = malloc(SgInfo.MaxList * sizeof (*SgInfo.ListSeitzMx));

   SgInfo.ListRotMxInfo
      = malloc(SgInfo.MaxList * sizeof (*SgInfo.ListRotMxInfo));

   if( Latt_N > 0 )
     if (AddInversion2ListSeitzMx(&SgInfo) < 0)
       goto ShelxReadError;

   if (AddLatticeTr2ListSeitzMx(&SgInfo, LatticeInfo) < 0)
     goto ShelxReadError;

  /* Loop through sym operators */
   if( VALID(&xyz_list))
   {
      node = xyz_list->head;
      do
      {
      xyz = node->data;

      if (ParseSymXYZ(xyz, &SeitzMx, STBF) < 0)
         goto ShelxReadError;

      if (Add2ListSeitzMx(&SgInfo, &SeitzMx) < 0)
         goto ShelxReadError;

      node = node->next;
      } while(node != NULL);
   }

  if (CompleteSgInfo(&SgInfo) < 0)
    goto ShelxReadError;

  /* process space group name */
  strcpy(sglabel, SgInfo.TabSgName->SgLabels);
  for (i=0 ; i<strlen(sglabel) ; i++)
    if (*(sglabel+i) == '_')
      *(sglabel+i) = ' ';

  i=strlen(sglabel);
  while(i>0)
    {
    if (*(sglabel+i) == '=')
      {
      i++;
      break;
      }
    i--;
    }

  /* fill in space group name */
  if(sglabel) /* Modified by C.Fisher 2005 */
    {
    strcpy(spgr, sglabel+i);
    trim(spgr);
    free(sglabel);
    }

  fclose(Fp);

  return natoms;

ShelxReadError:

  fprintf(stderr,"Error found in Shelx file %s\n", filename);
  fclose(Fp);
  return(1);
}
/******************************************************************************
 * read_ele().  Read element data from file.                                  *
 ******************************************************************************/
int          read_ele(spec_data *element, char *filename)
{
int        n=0;                  /* No of records read */
char       name[NLEN];
char       symbol[4];
double     mass, chg;
spec_data  *ele = element;
char       line[LLEN];
FILE       *Fe;

     if( (Fe = fopen(filename,"r")) == NULL)
        return -1;
     
     while(!feof(Fe))
     {
        if( (sscanf(get_line(line,LLEN,Fe,1),"%32s %4s %lf %lf", name, symbol, &mass, &chg) < 4) && !feof(Fe) )
            error("Unexpected format in \"%s\"", filename, strerror(errno));
        strcpy(ele->name, name);
        strcpy(ele->symbol, symbol);
        ele->mass = mass;
        ele->charge = chg;
        ele++; n++;
     }
     fclose(Fe);

     message(NULLI,NULLP,INFO,DATAREC,n,n==1?" ":"s ",filename);
     return n;
}
