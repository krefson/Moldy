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
 * tokenize()        Break down string into substrings and return no of substrings  *
 * space_minus()     Insert spaces before every '-' in s string except exponent     *
 * multi_char()      Return value of run-time string as if multichar constant       *
 * trim()            Remove leading and trailing spaces from string                 *
 * str_cut()         Divide string into alphabetical and numerical parts            *
 * read_ftype        Read filename suffix to guess file type                        *
 * cell_to_prim      Calculate primitive cell matrix from unit cell parms           *
 * read_cssr()       Read data from Cambridge Search and Structure Retrieval file   *
 * read_pdb()        Read data from Brookhaven Protein Databank (PDB) file          *
 * read_shak()       Read data from Schakal file                                    *
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
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "structs.h"
#include "messages.h"
#include "utlsup.h"
#include "sginfo.h"
#include "specdata.h"

#define TITLE_SIZE  80
#define DATAREC "%d record%ssuccessfully read from %s"

#define xtoupper(c) (islower(c) ? toupper(c) : c)

/*========================== External data references ========================*/
extern  const pots_mt   potspec[];           /* Potential type specification  */

extern int transformation_matrix (char *buf, T_RTMx *trans_matrix);
extern int symm_gen (T_RTMx matrix, mat_mp apos, char (*atype)[32], double *charge, int max, int natoms, int abegin, int aend);
extern void sgtransform (T_RTMx m, mat_mp x, mat_mp xp, int natoms);
extern int sgexpand (int maxnatoms, int natoms, vec_mt (*a_lst), char (*label)[32], double *charge, char *spgr);

/******************************************************************************
 * tokenize().                                                                *
 ******************************************************************************/
int   tokenize(char *buf, char **linev)
{
   char *p;
   char *sep = " +*/=:;,'\"<>()[]!?@\\^_`{}~";
   int linec = 0;

   while( (p = strtok(buf,sep) ) != NULL)
   {
      *(linev++) = p;
      linec++;
      buf = NULL;
   }

   return linec;
}
/******************************************************************************
 * space_minus().  Insert spaces before every '-' in s string except exponent.*
 ******************************************************************************/
void space_minus(char *str, int n)
{
   char *in = strdup(str), *savein = in, *strbegin = str;

   do
   {
      if( *in == '-' &&(
        (in > savein && in[-1] != 'E' && in[-1] != 'e') || in == savein))
         *str++ = ' ';
      *str++ = *in++;
   } while( *(in-1) && str-strbegin < n-1);
   free(savein);
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
 * trim().  Trim leading and trailing spaces from string variables            *
 ******************************************************************************/
char    *trim(char *s)
{
char    *t = s, *p;

   /* Remove leading spaces */
   while( *t == ' ' && *t != '\0')
      t++;

   strcpy(s, t);

   /* Now remove trailing spaces */
   p = s + strlen(s);
   while( p > s && !isalnum(*p))
      *p-- = '\0';

   return (s);
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
 * str_cut().  Separates alphabetical and numeric parts of a string.         *
 ******************************************************************************/
int     str_cut(char *in, char *out)
{
   int          i,j,k=strlen(in);
   int          value;

   for( i=0; i < k; i++)
     if( isdigit(in[i]) )
        break;

   for( j=0; j < k; j++)
     if( isalpha(in[j]) )
        break;

   if( i <= j )
   {
      strncat(out,in+j,k-j);    /* Copy alphabetical part */
      in[j] ='\0';              /* Truncate alphabetical part */
      value=atoi(in);           /* Convert numeric part */
   }
   else
   {
      strncat(out,in,i);        /* Copy alphabetical part */
      value=atoi(in+i);         /* Convert numeric part */
   }
   if( strrchr(out,'-') != NULL && value > 0 )
        value *= -1;

   return value;
}
/******************************************************************************
 * read_ftype().  Read suffix from filename.                                  *
 ******************************************************************************/
char    *read_ftype(char *filename)
{
char      *s = filename;

   while( *s != '.' && *s != '\0' )
      s++;

   if( *s == '.')
        s++;
   return (s);
}
/******************************************************************************
 * read_cssr().  Read structural data from cssr file.                         *
 ******************************************************************************/
int      read_cssr(char *filename, mat_mp h, char (*label)[32], vec_mp x, double *charge, char *title, char *spgr)
{
int      i, coord_sys=-1, nocell=0;
int      natoms;
double   cell[6];                                       /* Cell parameters */
mat_mt   hinv;
char     temp_name[6];
double   chg;
char     line[LLEN];
FILE     *Fp;
int	answer;

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   /* Read cssr file header */
   if( sscanf(get_line(line,LLEN,Fp,0),"%*38c%8lf%8lf%8lf",&cell[0],&cell[1],&cell[2]) < 3)
      nocell++;

   if( sscanf(get_line(line,LLEN,Fp,0),"%*21c%8lf%8lf%8lf    SPGR =%*3d %11c",
          &cell[3],&cell[4],&cell[5],spgr) < 4)
      nocell++;

   if( sscanf(get_line(line,LLEN,Fp,0),"%4d%4d %60c", &natoms, &coord_sys, title) < 2)
            error("EOF or unexpected format on line 3 in \"%s\"", filename);

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
int      read_pdb(char *filename, mat_mp h, char (*label)[32], vec_mp x, double *charge, char *title, char *spgr)

                                    /* Species' C of M coordinates */
{
int      natoms = 0;
double   cell[6];                                  /* Cell parameters */
char     temp_name[6], chg[2];
mat_mt   hinv;
double   tr[3];
char     line[LLEN], keyword[LLEN];
int      irow;
int      scaleflg = 0, crystflag = 0;
FILE     *Fp;

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   while( strcmp(strlower(keyword),"end") !=0 )
   {
       if( sscanf(get_line(line,LLEN,Fp,0),"%s",keyword) < 1)
           error("File \"%s\" has incorrect format", filename);

       if(strcmp(strlower(keyword),"title") == 0 )         /* Title */
       {
          if(isspace(line[9]))
             strncpy(title,line+10,TITLE_SIZE);
          else {
             strncat(title," ",TITLE_SIZE);
             strncat(title,line+10,TITLE_SIZE);
          }
       }

       if(strcmp(strlower(keyword),"cryst1") == 0)   /* CRYST1 - Specify unit cell */
       {
          if( sscanf(line,"CRYST1%9lf%9lf%9lf%9lf%9lf%9lf%11c",
             &cell[0],&cell[1],&cell[2],&cell[3],&cell[4],&cell[5],spgr) < 6 )
                error("Error in CRYST1 line of \"%s\" -- should have at least 6 parameters", filename);

          trim(spgr);
          crystflag++; /* Unit cell defined for this structure */
       }

       if(strncmp(strlower(keyword),"scale",5) == 0 )         /* SCALE */
       {
          sscanf(line, "SCALE%1d", &irow);
          if( irow < 1 || irow > 3 )
             error("Error in \"%s\" - Error in SCALE[1-3] - %6s unknown", filename, line);

          sscanf(line+10,"%10lf%10lf%10lf     %10lf",
                 &hinv[irow-1][0], &hinv[irow-1][1], &hinv[irow-1][2], &tr[irow-1]);

          scaleflg+=irow;
       }

       if(strcmp(strlower(keyword),"hetatm") == 0 || strcmp(strlower(keyword),"atom") == 0)
       {
          if( natoms > MAX_ATOMS )
              fprintf(stderr,"Warning: \"%s\" contains too many atoms! (%d)\n", filename, natoms);

          sscanf(line+12,"%4s%*1c%*3s%*2c%*4d%*4c%8lf%8lf%8lf", temp_name,
               &x[natoms][0], &x[natoms][1], &x[natoms][2]);
          str_cut(temp_name, label[natoms]);

          trim(label[natoms]);
          strcpy(chg,"  ");
          sscanf(line+78,"%2s", chg);
          charge[natoms] = str_cut(chg, chg);

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
int      read_shak(char *filename, mat_mt h, char (*label)[32], vec_mp x, double *charge, char *title, double *simbox)
{
int      i, natoms = 0;
double   cell[6];                                    /* Cell parameters */
double   box[6];
char     line[LLEN], combuf[LLEN], *buf, *bufend, last;
int      linec, linep; char *linev[32];
int      blankflg = 0;                    /* Blanks in atom labels?            */
int      multi_flag = 0;                  /* Flag for multiple units in file   */
int      dupl_flag, symm_unit_reset=1;
int      symm_unit_begin=0, symm_unit_end=0; /* Initialiation for lint only */
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

      space_minus(line, sizeof line);
      bufend = strchr(line,'\n');
      if( bufend == NULL)
         bufend = line+strlen(line);
      *bufend = '\0';
      buf = line;
      last = '\0';

      while(buf < bufend )
      {
         dupl_flag=0;
         linec = tokenize(strcpy(combuf, buf), linev);
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
 * read_ele().  Read elemental data from file.                                *
 ******************************************************************************/
int          read_ele(spec_data *element, char *filename)
{
int        n=0;                  /* No of records read */
char       name[NLEN];
char       symbol[4];
double     mass, chg;
spec_data  *ele = element;
char       line[LLEN];
char       *success_read;
FILE       *Fe;

     if( (Fe = fopen(filename,"r")) == NULL)
        return -1;

     for( ele = element; ele < element+NELEM; ele++ )
     {
        if( (sscanf(get_line(line,LLEN,Fe,1),"%32s %4s %lf %lf", name, symbol, &mass, &chg) < 4) && !feof(Fe) )
            error("Unexpected format in \"%s\"", filename);
        if( feof(Fe) )
            break;
        strcpy(ele->name, name);
        strcpy(ele->symbol, symbol);
        ele->mass = mass;
        ele->charge = chg;
        n++;
     }
     fclose(Fe);

     message(NULLI,NULLP,INFO,DATAREC,n,n==1?" ":"s ",filename);
     return n;
}
