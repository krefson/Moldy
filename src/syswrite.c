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
 *  Revision 2.1  2001/08/09 09:36:36  keith
 *  Incorporated Craig's new "Syswrite" utility.
 *
 * Revision 1.1  2001/04/25  18:27:41  fisher
 * Initial revision
 *
 */
#ifndef lint
static char *RCSid = "$Header: /home/minphys2/keith/CVS/moldy/src/syswrite.c,v 2.1 2001/08/09 09:36:36 keith Exp $";
#endif
#include "defs.h"
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include "stdlib.h"
#include "stddef.h"
#include "string.h"
#include <ctype.h>
#include "structs.h"
#include "messages.h"
#include "utlsup.h"
#include "sginfo.h"
#include "specdata.h"

/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
contr_mt                control;

#define  PDB              0
#define  CSSR             1
#define  SHAK             2

#define TITLE_SIZE  80

#define xtoupper(c) (islower(c) ? toupper(c) : c)

/*========================== External data references ========================*/
extern  const pots_mt   potspec[];           /* Potential type specification  */

/******************************************************************************
 * tokenize().                                                                *
 ******************************************************************************/

int str_cut (char *in, char *out);
extern int transformation_matrix (char *buf, T_RTMx *trans_matrix);
extern int symm_gen (T_RTMx matrix, mat_mp apos, char (*atype)[32], double *charge, int max, int natoms, int abegin, int aend);
extern void sgtransform (T_RTMx m, mat_mp x, mat_mp xp, int natoms);
extern int sgexpand (int maxnatoms, int natoms, vec_mt (*a_lst), char (*label)[32], double *charge, char *spgr);

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
void
cell_to_prim(double *cell, mat_mt h)
{
double	 ca, cb, cg, sg;               /* Cos and sin of cell angles */

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
 * read_suffix().  Read suffix from filename.                                 *
 ******************************************************************************/
char    *read_suffix(char *filename)
{
char      *s = filename;

   while( *s != '.' && *s != '\0' )
      s++;
  
   if( *s == '.')
        s++;
   return (s);
}
/******************************************************************************
 * read_pot().  Read potential data from file.                                *
 ******************************************************************************/
int       read_pot(char *potfile, pot_mp *pot_ptr, spec_data *spec, int max_id)
                   
                                  /* To be pointed at potpar array      */
                
                        /* Maximum no of site ids (nb. 1 less than in moldy) */
{
char       atom1[4], atom2[4];
double	   chg1, chg2;
double     p_tmp;
pot_mt     pot;
int	   i, n_items;
char       name[LLEN],             /* Temporary potential type name      */
           line[LLEN],             /* Store for input line from file     */
           pline[LLEN];            /* Used in pot'l paramater parsing    */
int        ptype=-1;               /* Potential type index               */
int        nerrs = 0;              /* Accumulated error count            */
int        idi, idj;
spec_data  *spi,*spj;              /* Temporary species pointers         */

FILE       *Fpot;

    if( (Fpot = fopen(potfile,"r")) == NULL)
        return ptype;

    n_items = sscanf(get_line(line,LLEN,Fpot), "%s", name);

    if( n_items <= 0 )
       message(NULLI,NULLP,FATAL,SYSEOF,"potential type specification");

    for(i = 0; potspec[i].name; i++)             /* Is 'name' a known type? */
       if(strcmp(strlower(name), potspec[i].name) == 0)
          break;

    if(! potspec[i].name)                        /* Did the loop find 'name'? */
       message(&nerrs,line,FATAL,UNKPOT,name);   /* no                        */
    ptype = i;                                   /* yes                       */

    while(sscanf(get_line(line,LLEN,Fpot),"%s",name) > 0
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

           if( ((!strcmp(atom1,spi->symbol)) && (chg1 == spi->charge) &&
             (!strcmp(atom2,spj->symbol)) && (chg2 == spj->charge)) || 
                ((!strcmp(atom1,spj->symbol)) && (chg1 == spj->charge) &&
                   (!strcmp(atom2,spi->symbol)) && (chg2 == spi->charge)) )
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

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   /* Read cssr file header */
   if( sscanf(get_line(line,LLEN,Fp),"%*38c%8lf%8lf%8lf",&cell[0],&cell[1],&cell[2]) < 3)
      nocell++;

   if( sscanf(get_line(line,LLEN,Fp),"%*21c%8lf%8lf%8lf    SPGR =%*3d %11c",
          &cell[3],&cell[4],&cell[5],spgr) < 4)
      nocell++;

   if( sscanf(get_line(line,LLEN,Fp),"%4d%4d %60c", &natoms, &coord_sys, title) < 2)
            error("EOF or unexpected format on line 3 in \"%s\"", filename);

   trim(spgr);

   if( ! (cell[0] > 0 && cell[1] > 0 && cell[2] > 0 &&
          cell[3] > 0 && cell[3] < 180.0 &&
          cell[4] > 0 && cell[4] < 180.0 &&
          cell[5] > 0 && cell[5] < 180.0) )
            error("Invalid unit cell in \"%s\"", filename);

   if( nocell && ! coord_sys )
      error("No unit cell parameters supplied for fractional coordinates");

   if( fgets(line,sizeof(line),Fp) == NULL )
      error("Unexpected end of input on line 4 of \"%s\"", filename);


   if( ! coord_sys && cell[0] >= 0.0 )
   {
      cell_to_prim(cell, h);
   }
   else
   {
      h[0][0] = h[1][1] = h[2][2] = 1.0;
      h[0][1] = h[1][0] = h[0][2] = h[2][0] = h[1][2] = h[2][1] = 0.0;
   }
   invert(h, hinv);
   mat_vec_mul(hinv, x, x, natoms);

      
   for( i = 0; i < natoms; i++)
   {
      if( sscanf(get_line(line,LLEN,Fp),"%*4d %4s  %9lf %9lf %9lf %*4d%*4d%*4d%*4d%*4d%*4d%*4d%*4d %7lf",
               temp_name, &x[i][0], &x[i][1], &x[i][2], &chg) < 5 )
         error("EOF or unexpected format in atom definitions in \"%s\"", filename);

      charge[i] = chg; 

      str_cut(temp_name, label[i]); /* Convert cssr atom name to atom symbol */
      trim(label[i]);
   }
   fclose(Fp);

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
int      scaleflg = 0;
FILE     *Fp;

   if( (Fp = fopen(filename,"r")) == NULL)
      error("Failed to open configuration file \"%s\" for reading", filename);

   while( strcmp(strlower(keyword),"end") !=0 )
   {
       if( sscanf(get_line(line,LLEN,Fp),"%s",keyword) < 1)
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
   if( scaleflg != 6 )
   {
      cell_to_prim(cell,h);
      invert(h,hinv);
   }
   else
      cell_to_prim(cell,h);

   mat_vec_mul(hinv, x, x, natoms);   /* Convert to fractional coords */

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
 * copy_atom_data()  Make copy of atom/molecule names and positions           *
 ******************************************************************************/
int copy_atom_data(site_mp site, double *charge, char (*label)[32], int natoms)
{
   int iatom;
   for(iatom = 0; iatom < natoms; iatom++)
   {
       strncpy((site+iatom)->name, label[iatom], NLEN);
       (site+iatom)->charge = charge[iatom];
       (site+iatom)->mass = 0.0;
       (site+iatom)->pad = -1;
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
   char		elefile[PATHLN];
   char		potfile[PATHLN];
   char		title[TITLE_SIZE];
   double	cell[6];                      /* Cell parameters */
   mat_mt	h;                            /* Hessian */
   int		insw=-1;                      /* Switch for input format */
   int		spectot=0, atomtot=0;         /* No of species, atoms */
   int		i,j,k,idij,ip;
   int		nflag;                        /* Flag for site type matching */
   int          n_elem;                       /* No of records read from element data file */
   spec_data    element[NELEM];               /* Element info */
   spec_data    spec[NSPEC];                  /* Species info */
   site_mt      *site, *st;                   /* Site info */ 
   pot_mt       *pot_par;                     /* Potential parameters */
   int		ptype = -1, n_potpar=NPOTP;
   char         spgr[16];                     /* Space Group in Herman Maugain form */
   double       simbox[3] = {1,1,1};          /* Simulation box repeat factors */
   double	x[MAX_ATOMS][3];              /* C of M coordinates */
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

      filetype = strlower(read_suffix(filename));

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

   if( strcmp(spgr,"P 1")  && atomtot > 0 )
      atomtot = sgexpand(MAX_ATOMS, atomtot, x, label, charge, spgr);

   site = calloc(sizeof(site_mt),atomtot);
   copy_atom_data(site, charge, label, atomtot);

   /* Create full path name, but don't exceed max length of string */
   strncat(strncat(elefile, ELEPATH, PATHLN-strlen(elefile)), elename, PATHLN-strlen(elefile));

   /* Check if element data file is available */
   if( (n_elem = read_ele(element, elefile)) < 0 )
      message(NULLI,NULLP,WARNING,NOELEM,elefile);

   /* Loop to determine no of different species */
   for(st = site; st < site+atomtot; st++)
   {
      toupper(st->name[0]);   /* First character of symbol uppercase */
      strlower(st->name+1);   /* Remainder of symbol lowercase */

      nflag = 0;    /*  Flag for whether site (1) assigned or (0) not assigned to species */
      for( i = 0; i < spectot; i++)
         if( !strcmp(st->name,(spec+i)->symbol)          /* Does site match any species identified so far? */
                 && (st->charge == (spec+i)->charge) )
         {
              st->pad = i;                 /* Assign this site to species i */
              (spec+i)->nmols++;
              nflag++;
         }
      if (!nflag)			   /* If no matches, create new species type */
      {
         for( j=0; j < n_elem; j++)        /* Search element/species data for match with site */
           if( !strcmp(st->name,(element+j)->symbol) )
           {
              strncpy((spec+spectot)->name,(element+j)->name,NLEN);
              st->mass = (element+j)->mass;
              if( insw == SHAK)
                st->charge = (element+j)->charge;
              nflag++;
              break;
           }

         if( !nflag)
             strcpy((spec+spectot)->name,st->name);

         strcpy((spec+spectot)->symbol,st->name);
         (spec+spectot)->mass = st->mass;
         (spec+spectot)->charge = st->charge;
         (spec+spectot)->nmols=1;
         st->pad = spectot;
         spectot++;
      }
   }

   /* Modify species name if same name but different charge to another species */
   for( i = 0; i < spectot-1; i++)
   {
      k = 0;
      for(j = i+1; j < spectot; j++)
         if( !strcmp((spec+i)->name,(spec+j)->name) )
         {
            k++;
            add_suffix((spec+j)->name, k);
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
   if( (potname != NULL) && ((ptype = read_pot(potfile, &pot_par, spec, spectot)) < 0) )
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
       printf("%-14s %d\n",(spec+i)->name,(spec+i)->nmols);
       printf("%-3d    0    0    0  %12.12g %12.12g  %s\n",i+1,
                     (spec+i)->mass,(spec+i)->charge,(spec+i)->symbol);
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

      printf("%-14s %10g %10g %10g\n",(spec+((site+i)->pad))->name,x[i][0],x[i][1],x[i][2]); 
   }

   puts("end");

   return 0;
}
