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
 * spgrexp code for expanding unit cells from primitive coordinates and space group   *
 **************************************************************************************
 *  Revision Log
 *  $Log$
 *  Revision 2.4.8.1  2003/07/29 09:45:20  moldydv
 *  Improvements, added functions and bug fixes for several utilities.
 *
 *  Revision 2.4  2002/09/19 09:26:30  kr
 *  Tidied up header declarations.
 *  Changed old includes of string,stdlib,stddef and time to <> form
 *
 *  Revision 2.3  2002/06/21 11:29:07  kr
 *  Got rid of K&R varargs-compatibility stuff.
 *
 *  Revision 2.2  2001/08/09 11:46:56  keith
 *  Tidied up against some compiler warnings.
 *  Added license file for SgInfo routines with permission of
 *  Ralf W. Grosse-Kunstleve
 *
 *  Revision 2.1  2001/08/09 09:36:35  keith
 *  Incorporated Craig's new "Syswrite" utility.
 *
 * Revision 1.1  2001/03/23  06:32:03  fisher
 * Initial revision
 *
 */

#include <stdio.h>
#include <ctype.h>
#include "sginfo.h"
#include "specdata.h"

#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "messages.h"
#include "defs.h"

#define SPEC_TOL 0.005    /* Special Position Tolerance */

void error(char *, ...);
/******************************************************************************
 * Structure declarations                                                     *
 ******************************************************************************/
struct symm_s {
struct symm_s *next;
T_RTMx trans_matrix;
};

typedef struct symm_s symm_t;

/******************************************************************************
 * sgtransform() Apply transformation matrix to a single set of co-ordinates. *
 * Allow "in-situ" copying, ie x == xp.                                       *
 ******************************************************************************/
void    sgtransform(T_RTMx m, mat_mp x, mat_mp xp, int natoms)
{
   int i, iatom;
   float y[3];
   for(iatom = 0; iatom < natoms; iatom++)
   {
      for(i = 0; i < 3; i++)
         y[i] = m.s.R[3*i + 0]*x[iatom][0] + m.s.R[3*i + 1]*x[iatom][1]
              + m.s.R[3*i + 2]*x[iatom][2] + m.s.T[i]/(1.0*STBF);
      for(i = 0; i < 3; i++)
         xp[iatom][i] = y[i];
   }
}
/******************************************************************************
 * symm_gen() Generate atoms by symmetry operation in matrix form.            *
 ******************************************************************************/
int     symm_gen(T_RTMx matrix, mat_mp apos, char (*atype)[NLEN], double *charge, int max, int natoms, int abegin, int aend)
{
   int          iatom, i, iatt, idupl;
   double       trial_pos[3], da, db, dc;

   for(iatom = abegin; iatom < aend; iatom++)
   {
      sgtransform(matrix, apos+iatom, &trial_pos, 1);

      idupl = 0;
      for(iatt=0; iatt<natoms; iatt++)
      {
         da = fmod(16.5+trial_pos[0]-apos[iatt][0],1.0)-0.5;
         db = fmod(16.5+trial_pos[1]-apos[iatt][1],1.0)-0.5;
         dc = fmod(16.5+trial_pos[2]-apos[iatt][2],1.0)-0.5;

#ifdef DEBUG
         fprintf(stderr,"%d %d %3s: %f   %f %f %f\n",natoms, iatt, atype[iatom],
               sqrt(SQR(da)+SQR(db)+SQR(dc)), da, db, dc);
#endif
         if( SQR(da) + SQR(db) + SQR(dc) < SQR(SPEC_TOL))
         {
            if(strcmp(atype[iatt],atype[iatom]))
            {
               error("Generated %s atom overlaps %s atom (%d,%d) ",atype[iatom],atype[iatt],iatom,iatt);
#ifdef DEBUG
               fprintf(stderr,"Generated atom (%f,%f,%f); orig atom (%f,%f,%f); dist %f\n",
                     trial_pos[0],trial_pos[1],trial_pos[2],apos[iatt][0],apos[iatt][1],apos[iatt][2],sqrt(SQR(da)+SQR(db)+SQR(dc)));
#endif
            }
            else
            {
#ifdef DEBUG
               fprintf(stderr,"Generated %s atom overlaps %s atom (%d,%d)\n",atype[iatt],atype[iatom],iatt,iatom);
#endif
               idupl++;
               break;
            }
         }
      }
#ifdef DEBUG
      fprintf(stderr,"Generated atom %d (%f,%f,%f) idupl=%d\n",natoms,trial_pos[0],trial_pos[1],trial_pos[2],idupl);
#endif
      if( idupl == 0 )
      {
         if( natoms >= max )
         {
            fprintf(stderr,"No of atoms >= %d after applying symmetry operations.\n",max);
            return max;
         }
         for(i = 0; i < 3; i++)
            apos[natoms][i] = fmod(16.0+trial_pos[i],1.0);
         strncpy(atype[natoms],atype[iatom],NLEN);
         charge[natoms] = charge[iatom];
         natoms++;
      }
   }
   return natoms;
}
/******************************************************************************
 *  buildsginfo  build space group information                                *
 ******************************************************************************/
int BuildSgInfo(T_SgInfo *SgInfo, const char *SgName)
{
  int   VolLetter = 'A';
  const T_TabSgName  *tsgn;


  tsgn = FindTabSgNameEntry(SgName, VolLetter);
  if (tsgn == NULL) return -1; /* no matching table entry */
  SgName = tsgn->HallSymbol;

  /* Allocate memory for the list of Seitz matrices and
     a supporting list which holds the characteristics of
     the rotation parts of the Seitz matrices.
   */

  SgInfo->MaxList = 192; /* absolute maximum number of symmetry ops */

  SgInfo->ListSeitzMx
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));

  if (SgInfo->ListSeitzMx == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  SgInfo->ListRotMxInfo
    = malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));

  if (SgInfo->ListRotMxInfo == NULL) {
    SetSgError("Not enough core");
    return -1;
  }

  /* Initialize the SgInfo structure
   */

  InitSgInfo(SgInfo);
  SgInfo->TabSgName = tsgn; /* in case we know the table entry */

  /* Translate the Hall symbol and generate the whole group
   */

  ParseHallSymbol(SgName, SgInfo);
  if (SgError != NULL) return -1;

  return CompleteSgInfo(SgInfo);
}
/******************************************************************************
 *  sgexpand    calculate total number of atoms from space group symmetry     *
 ******************************************************************************/
int sgexpand(int maxnatoms, int natoms, vec_mt *a_lst, char (*label)[NLEN], double *charge, char *spgr)
{
  T_SgInfo  SgInfo;
  T_RTMx    inversion = {{{-1,0,0, 0,-1,0, 0,0,-1},{0,0,0}}};
  T_RTMx    centre    = {{{1, 0, 0, 0, 1, 0, 0, 0, 1}, {0, 0, 0}}};
  int isymm;

  if(BuildSgInfo(&SgInfo, spgr) < 0)
     error("sgexpand: BuildSgInfo failed");

#ifdef DEBUG
  fprintf(stderr,"Space Group %s %s\n",spgr,SgInfo.HallSymbol);
#endif
  /*
   * Apply space group symmetry ops.  Start from one to avoid unity!
   */
  for(isymm = 1; isymm < SgInfo.nList; isymm++)
  {

#ifdef DEBUG
     fprintf(stderr,"Applying matrix [%d %d %d; %d %d %d; %d %d %d] and translation (%d,%d,%d)\n",
	   SgInfo.ListSeitzMx[isymm].s.R[0],SgInfo.ListSeitzMx[isymm].s.R[1],SgInfo.ListSeitzMx[isymm].s.R[2],
	   SgInfo.ListSeitzMx[isymm].s.R[3],SgInfo.ListSeitzMx[isymm].s.R[4],SgInfo.ListSeitzMx[isymm].s.R[5],
	   SgInfo.ListSeitzMx[isymm].s.R[6],SgInfo.ListSeitzMx[isymm].s.R[7],SgInfo.ListSeitzMx[isymm].s.R[8],
	   SgInfo.ListSeitzMx[isymm].s.T[0],SgInfo.ListSeitzMx[isymm].s.T[1],SgInfo.ListSeitzMx[isymm].s.T[2]);
#endif
     natoms=symm_gen(SgInfo.ListSeitzMx[isymm], a_lst, label, charge, maxnatoms, natoms, 0, natoms);
  }
  /*
   * Add inversion centre if necessary
   */
  if( SgInfo.Centric == -1 )
  {
#ifdef DEBUG
     fprintf(stderr,"Applying matrix [%d %d %d; %d %d %d; %d %d %d] and translation (%d,%d,%d)\n",
	   inversion.s.R[0],inversion.s.R[1],inversion.s.R[2],
	   inversion.s.R[3],inversion.s.R[4],inversion.s.R[5],
	   inversion.s.R[6],inversion.s.R[7],inversion.s.R[8],
	   inversion.s.T[0],inversion.s.T[1],inversion.s.T[2]);
#endif
     natoms=symm_gen(inversion, a_lst, label, charge, maxnatoms, natoms, 0, natoms);
  }
  /*
   * Add lattice centre2 if necessary
   */
  if( SgInfo.LatticeInfo->Code != 'P' )
     for(isymm = 0; isymm < SgInfo.LatticeInfo->nTrVector; isymm++)
     {
	centre.s.T[0] = SgInfo.LatticeInfo->TrVector[3*isymm];
	centre.s.T[1] = SgInfo.LatticeInfo->TrVector[3*isymm+1];
	centre.s.T[2] = SgInfo.LatticeInfo->TrVector[3*isymm+2];
#ifdef DEBUG
	fprintf(stderr,"Applying matrix [%d %d %d; %d %d %d; %d %d %d] and translation (%d,%d,%d)\n",
	      centre.s.R[0],centre.s.R[1],centre.s.R[2],
	      centre.s.R[3],centre.s.R[4],centre.s.R[5],
	      centre.s.R[6],centre.s.R[7],centre.s.R[8],
	      centre.s.T[0],centre.s.T[1],centre.s.T[2]);
#endif
	natoms=symm_gen(centre, a_lst, label, charge, maxnatoms, natoms, 0, natoms);
     }
  
  return natoms;
}
