/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1999, 2001 Craig Fisher
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
 * ransub    	Code for randomly substituting species                                *
 *              in Moldy configuration files                                          *
 *		Randomly replaces species "m" with n molecules of species "u"	      *
 *		Output written in Moldy system specification format		      *
 ************************************************************************************** 
 *  Revision Log
 *  $Log: ransub.c,v $
 *  Revision 1.17.10.2  2003/08/01 00:45:02  moldydv
 *  Added afree and sgexpand to function declarations.
 *
 *  Revision 1.17.10.1  2003/07/29 09:36:04  moldydv
 *  Polyatomic dopants can now be read in xtl and xyz formats.
 *  Options -f,-t,-p added for specifying Euler angles phi, theta and psi of all polyatomic dopant molecules.
 *  Corrected file pathname specifier.
 *  Quaternions now written correctly when replacing monatomic with polyatomic.
 *  Changed 'abs' to 'fabs'.
 *  Replaced explicit lengths of name variables with NLEN.
 *  Fixed bug when counting substituted positions.
 *
 *  Revision 1.17  2002/09/20 15:45:08  kr
 *  Corrected some C errors/removed dependence on gcc extensions.
 *
 *  Revision 1.16  2002/09/19 09:26:30  kr
 *  Tidied up header declarations.
 *  Changed old includes of string,stdlib,stddef and time to <> form
 *
 *  Revision 1.15  2002/09/18 09:59:18  kr
 *  Rolled in several changes by Craig Fisher:
 *  Ransub can now read polyatomic species
 *  Syswrite can handle polyatomics from CSSR PDB or SCHACKAL files
 *
 *  Revision 1.14  2002/06/21 11:18:10  kr
 *  Got rid of K&R varargs-compatibility stuff.
 *
 *  Revision 1.13   2002/05/02 14:57:19  fisher
 *  Modified sys-spec output routine to add dopant (solute) species to end of list.
 *  Check added to see if dopant species already present in system.
 *  Sorting of substituted positions list for improved efficiency.
 *  If no potentials specified, assume same pots as first site of substituted species.
 *  Substitution aborted if solute and solvent species have same name.
 *  External data files now only read if no. of solute particles > 0
 *  Tidied up error messages.
 *
 *  Revision 1.12  2001/08/09 11:46:55  keith
 *  Tidied up against some compiler warnings.
 *  Added license file for SgInfo routines with permission of
 *  Ralf W. Grosse-Kunstleve
 *
 *  Revision 1.11  2001/08/08 16:31:52  keith
 *  Incorporated Craig's modifications.
 *  Compiles but not properly tested.
 *
 *  Revision 1.7  2001/04/24 16:17:21  fisher
 *  elem.h renamed specdata.h to avoid confusion with other software packages.
 *  Modifications and improvements regarding options -y and -e.
 *
 *  Revision 1.6.2.1  2001/03/27 17:42:42  keith
 *  New version from Craig:
 *
 *  Removed relevant definitions to header file elem.h.
 *  Option -y added for reading potential parameters from text file.
 *  Option -e added for reading species data from text file.
 *  Minor modifications to program structure and variable names for clarity.
 *
 * Revision 2.1  2001/03/23  01:29:11  fisher
 * Removed shared definitions to header file elem.h.
 *
 * Revision 2.0  2001/02/19  06:18:32  fisher
 * Option -y added for reading potential parameters from text file.
 * Option -e added for reading species data from text file.
 * Minor modifications to program structure and variable names for clarity.
 *
 *  Revision 1.6  2000/02/16 11:46:09  craig
 *  Incorporated site-pbc branch "bekker" into main "Beeman" branch.
 *
 *  Revision 1.6  2000/02/16 11:46:09  craig
 *  Corrected memory leak when performing strcmp of NULL value.
 *
 *  Revision 1.5  1999/10/29 16:44:28  keith
 *  Bugfixes.
 *
 *  Revision 1.5  1999/10/25 10:41:46  craig
 *  Corrected memory leak when performing strcmp of NULL value.
 *  Added check for correct entry of substituting species' name.
 *
 *  Revision 1.4  1999/10/11 14:07:08  keith
 *  Removed common utility functions to "utlsup.c".
 *
 *  Revision 1.3  1999/09/24 11:02:28  keith
 *  Minor changes to random seeder and terminology.
 *
 *  Revision 1.3  1999/09/24 16:47:36  craig
 *  Minor changes to random seeder and terminology.
 *
 *  Revision 1.2  1999/09/21 11:16:29  keith
 *  Fixed compile problem on pre-ANSI compilers
 *
 *  Revision 1.1  1999/07/22 14:02:26  keith
 *  Initial revision
 *
 *  Revision 1.6  1999/06/24 16:05:44  craig
 *  Improved randomization of random number reseeder.
 *
 *  Revision 1.5  1999/06/03 15:39:34  craig
 *  Corrected memory freeing of dump limits.
 *  Tidied up use of structure variable 'dop'.
 *  Added loop to check if species being replaced exists.
 *
 *  Revision 1.4  1999/04/08 17:57:29  craig
 *  Options to specify dopant mass, charge and symbol added.
 *
 *  Revision 1.3  1999/03/23 15:09:45  craig
 *  Removed unnecessary variable 'is' from sys_spec_out.
 *
 *  Revision 1.2  1999/03/12 15:15:39  craig
 *  Altered energy conversion units to be consistent with defs.h values
 *
 *  Revision 1.1  1999/03/10 18:10:19  craig 
 *  Corrected bug limiting max no of substitutions to no of first species
 *  Upper limit to no of substituted species set to total no in system 
 *
 *  Revision 1.0  1999/03/05 17:41:12  craig 
 *  Initial revision
 *
 */
#ifndef lint
static char *RCSid = "$Header: /usr/users/moldy/CVS/moldy/src/ransub.c,v 1.17.10.2 2003/08/01 00:45:02 moldydv Exp $";
#endif  

#include "defs.h"
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include "structs.h"
#include "messages.h"
#include "utlsup.h"
#include "specdata.h"
#include "sginfo.h"

void	read_sysdef(FILE *file, system_mp system, spec_mp *spec_pp, site_mp *site_info, pot_mp *pot_ptr);
void	initialise_sysdef(system_mp system, spec_mt *species, site_mt *site_info, quat_mt (*qpf));
void	re_re_header(FILE *restart, restrt_mt *header, contr_mt *contr);
void	re_re_sysdef(FILE *restart, char *vsn, system_mp system, spec_mp *spec_ptr, site_mp *site_info, pot_mp *pot_ptr);
void	allocate_dynamics(system_mp system, spec_mt *species);
void	lattice_start(FILE *file, system_mp system, spec_mp species, quat_mt (*qpf));
void	read_restart(FILE *restart, char *vsn, system_mp system, int av_convert);
void	init_averages(int nspecies, char *vsn, long int roll_interval, long int old_roll_interval, int *av_convert);
void    conv_potentials(const unit_mt *unit_from, const unit_mt *unit_to, pot_mt *potpar, int npotpar, int ptype, site_mt *site_info, int max_id);
void    q_to_rot(real *quat,mat_mt rot);
int	getopt(int, char *const *, const char *);
void    afree(gptr *p);
char    *atime(void);
double  det(mat_mt a);
char	*read_ftype(char *filename);
int     read_ele(spec_data *element, char *filename);
int     read_pdb(char *, mat_mp, char (*)[NLEN], vec_mp, double *, char *, char *);
int     read_cssr(char *, mat_mp, char (*)[NLEN], vec_mp, double *, char *, char *);
int     read_shak(char *, mat_mp, char (*)[NLEN], vec_mp, double *, char *, double *);
int     read_xtl(char *, mat_mp, char (*)[NLEN], vec_mp, double *, char *, char *);
int     read_xyz(char *, mat_mp, char (*)[NLEN], vec_mp, char *);
int     sgexpand(int , int , vec_mt *, char (*)[NLEN], double *, char *);
/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
extern const  unit_mt prog_unit;
static  unit_mt input_unit = {MUNIT, LUNIT, TUNIT, _ELCHG};
contr_mt               control;

/* Time units for different energy units */
#define EV 1.018050697e-14  		/* electron volts */
#define KJMOL 1.0e-13 			/* kilojoules per mole */
#define KCALS 4.88882131e-14 		/* kilocalories per mole */
#define E2A 2.682811715e-15 		/* electron charge squared per angstrom */

#define OFF               0
#define ON                1

/*========================== External data references ========================*/

extern  const pots_mt   potspec[];          /* Potential type specification */

/******************************************************************************
 * prep_pot().  Convert units from system to user specified units             *
 ******************************************************************************/
int     prep_pot(system_mt *system, site_mt *site_info, pot_mt *potpar)
{
int   nunits;

      fputs("What energy units would you like the potentials in:\n",stderr);
      fputs("(1) eV, (2) kJ/mol, (3) kcal/mol, or (4) e**2/A",stderr);
      nunits = get_int(" ? ", 1, 4);
      switch(nunits)
      {
      case 1:
         input_unit.t = EV;
         break;
      case 2:
         input_unit.t = KJMOL;
         break;
      case 3:
         input_unit.t = KCALS;
         break;
      case 4:
         input_unit.t = E2A;
         break;
      }
      conv_potentials(&prog_unit, &input_unit, potpar, system->n_potpar,
          system->ptype, site_info, system->max_id);

      return 0;
}
/******************************************************************************
 * read_pot().  Read potential data from file.                                *
 ******************************************************************************/
static
int        read_pot(char *potfile, pot_mp *pot_ptr, int max_id, int new_sites, site_mt site_info[])
{
char       atom1[4], atom2[4];
double     chg1, chg2;
double     p_tmp;
pot_mt     pot;
int        i,idi,idj, n_items;
char       name[LLEN],             /* Temporary species name             */
           line[LLEN],             /* Store for input line from file     */
           pline[LLEN];            /* Used in pot'l paramater parsing    */
int        ptype = -1;             /* Potential type index               */
int        nerrs = 0;              /* Accumulated error count            */
site_mt    spi, spj;

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
       if(sscanf(line,"%s %lf %s %lf %[^#]",&atom1,&chg1,&atom2,&chg2,pline) <= 2)
            message(&nerrs,line,ERROR,NOPAIR);
       else
       {
                                             /* Now read in parameters   */
          (void)strcat(pline, "$");          /*   Add marker to end      */
          while(n_items < NPOTP && sscanf(pline,"%lf %[^#]", &p_tmp, pline) > 1 )
              pot.p[n_items++] = p_tmp;
       }

       for( idi = 1; idi < max_id; idi++)
       {
          spi = site_info[idi];
          for( idj = max_id; idj < max_id+new_sites; idj++)
          {
             spj = site_info[idj];
             if( ((!strcmp(atom1,spi.name)) && (chg1 == spi.charge) &&
              (!strncmp(atom2,spj.name,4)) && (chg2 == spj.charge)) ||
                 ((!strcmp(atom1,spj.name)) && (chg1 == spj.charge) &&
                    (!strcmp(atom2,spi.name)) && (chg2 == spi.charge)) )
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
 * sort_pos. Sort positions to be replaced into ascending order               *
 ******************************************************************************/
void
sort_pos(int *pos, int n)   /* A simple bubble sort algorithm */
{
   int  i,j,temp;

   for( i = 0; i < n; i++)
     for(j = i + 1; j < n; j++)
        if( pos[i] > pos[j])
        {
           temp = pos[i];
           pos[i] = pos[j];
           pos[j] = temp;
        }
   return;
}
/******************************************************************************
 * quat_gen. Generate random set of quaternions.                              *
 ******************************************************************************/
void
quat_gen(quat_mt quaternion)  /* Randomly generates Euler angles and converts them to quaternions */
{
   double   euler[3];
   int i;

   srand(time(NULL)+rand());   /* Randomly re-seed random generator */

   for(i = 0; i < 3; i++)
   {
      euler[i] = rand() % 360;    /* Generate random angle in degrees */
      euler[i] *= PI/180.0;       /* Convert to radians */
   }

   quaternion[0] = cos(0.5*euler[1])*cos(0.5*(euler[0]+euler[2]));
   quaternion[1] = sin(0.5*euler[1])*cos(0.5*(euler[0]-euler[2]));
   quaternion[2] = sin(0.5*euler[1])*sin(0.5*(euler[0]-euler[2]));
   quaternion[3] = cos(0.5*euler[1])*sin(0.5*(euler[0]+euler[2]));

   return;
}
/******************************************************************************
 * ran_quat. Assign new random quaternions to polyatomic species.             *
 ******************************************************************************/
void
ran_quat(system_mt *system, spec_mt *species, char *molname)
{
    spec_mp	spec;
    int		i, imol;

    for(spec = species; spec < species+system->nspecies; spec++)
      if( !strcmp(strlower(spec->name), molname) )
         if( spec->quat != NULL )   /* If polyatomic, generate quaternions */
            for(imol = 0; imol < spec->nmols; imol++)
               quat_gen(spec->quat[imol]);

   return;
}
/******************************************************************************
 * sys_spec_out().  Write a system configuration to stdout in the form of a   *
 * system specification file for MOLDY                                        *
 ******************************************************************************/
void
sys_spec_out(system_mt *system, spec_mt *species, spec_mt *dopant, char *molname,
                  int *positions, site_mt *site_info, double *euler, pot_mt *potpar)
{
   spec_mt      *spec;
   double       a, b, c, alpha, beta, gamma;
   mat_mp       h = system->h;
   mat_mt       hinv;
   int          i, imol, isite, ipos;
   int          idi, idj, idij, ip;
   int          specmol, id=-1;
   int          n_potpar = system->n_potpar;
   char         *specname;
   int          max_id = system->max_id;
   int          dflag = 0;
   int		namelength = 0;
   quat_mt      quaternion;
   boolean	quat_valid = true;
   vec_mt       *site = ralloc(dopant->nsites);
   vec_mt       pf_cofm, pf_origin;
   vec_mt       spec_origin;
   mat_mt       rot_mat;

   zero_real(pf_cofm,3);
   invert(h,hinv);

   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

   if( site_info[0].pad )                               /* Calculate dopant centre of mass */
   {
      if( dopant->nsites > 1 )
      {
         for( isite = 0; isite < dopant->nsites; isite++)
         {
            if( site_info[dopant->site_id[isite]].mass < 0)
               site_info[dopant->site_id[isite]].mass = 0;
            if( site_info[dopant->site_id[isite]].charge == 1e6)
               site_info[dopant->site_id[isite]].charge = 0;

            for( i=0; i < 3; i++)
               pf_cofm[i] += dopant->p_f_sites[isite][i]*site_info[dopant->site_id[isite]].mass;
            dopant->mass += site_info[dopant->site_id[isite]].mass;
         }
         if(dopant->mass < 1.0)              /* Lighter than 1 amu ?              */
            message(NULLI,NULLP,FATAL,ZMASS,dopant->name,dopant->mass);

         for( i=0; i < 3; i++)
         {
            pf_cofm[i] /= dopant->mass;     /* Dopant centre of mass in principal frame */
            pf_cofm[i] -= dopant->p_f_sites[0][i];  /* COFM relative to first site in molecule */
         }
      }
      else
         for( i=0; i < 3; i++)
            pf_cofm[i] = dopant->p_f_sites[0][i];
   }

/* Write header for sys_spec file */
   (void)printf("# System specification file written by RANSUB on %s\n",atime());

/* Write site data for each molecule */
   for(spec = species, i = 0; spec < species+system->nspecies; spec++)
   {
      /* Check if species matches species to be substituted */
      if( (molname != NULL) && !strcmp(strlower(spec->name), molname) )
      {
         specmol = spec->nmols - dopant->nmols; /* Subtract number of substituting species */
         id = i;    /* Label to identify species being substituted */
      }
      else
         specmol = spec->nmols;

      /* Check if dopant species already present in system */
      if( (dopant->name != NULL) && !strcmp(strlower(spec->name),dopant->name) )
      {
         specmol += dopant->nmols;
         message(NULLI,NULLP,WARNING,"Species `%s' already present in system - properties left unchanged",dopant->name);
         dflag++;
      }

      i++;

      if( specmol > 0 )      /* Write data for original species */
      {
         (void)printf("%s  %d  %s\n", spec->name, specmol,
                    spec->framework ? "framework" : "");
         for(isite=0; isite < spec->nsites; isite++)
            (void)printf("%d %9g %9g %9g %9g %9g %s\n",
                        spec->site_id[isite],
                        fabs(spec->p_f_sites[isite][0]-spec->p_f_sites[0][0]) > 1e-6 ?
                             spec->p_f_sites[isite][0]-spec->p_f_sites[0][0] : 0.0,
                        fabs(spec->p_f_sites[isite][1]-spec->p_f_sites[0][1]) > 1e-6 ?
                             spec->p_f_sites[isite][1]-spec->p_f_sites[0][1] : 0.0,
                        fabs(spec->p_f_sites[isite][2]-spec->p_f_sites[0][2]) > 1e-6 ?
                             spec->p_f_sites[isite][2]-spec->p_f_sites[0][2] :0.0,
                        site_info[spec->site_id[isite]].mass,
                        site_info[spec->site_id[isite]].charge,
                        site_info[spec->site_id[isite]].name);
      }
      namelength = MAX(namelength,strlen(spec->name)+1);
   }
   if( !dflag && id >= 0 && dopant->nmols > 0 && strncmp("#",dopant->name,1) )  /* Write data for species added */
   {
      (void)printf("%s  %d  %s\n", dopant->name, dopant->nmols,
           (species+id)->framework ? "framework" : "");

      for(isite=0; isite < dopant->nsites; isite++)
      {
         (void)printf("%d %9g %9g %9g %9g %9g %s\n",
               dopant->site_id[isite],
               fabs(dopant->p_f_sites[isite][0]-dopant->p_f_sites[0][0]) > 1e-6 ?
                    dopant->p_f_sites[isite][0]-dopant->p_f_sites[0][0] : 0.0,
               fabs(dopant->p_f_sites[isite][1]-dopant->p_f_sites[0][1]) > 1e-6 ?
                    dopant->p_f_sites[isite][1]-dopant->p_f_sites[0][1] : 0.0,
               fabs(dopant->p_f_sites[isite][2]-dopant->p_f_sites[0][2]) > 1e-6 ?
                    dopant->p_f_sites[isite][2]-dopant->p_f_sites[0][2] : 0.0,
         site_info[dopant->site_id[isite]].mass < 0 ? site_info[(species+id)->site_id[0]].mass:
                    site_info[dopant->site_id[isite]].mass,
         site_info[dopant->site_id[isite]].charge == 1e6 ? site_info[(species+id)->site_id[0]].charge:
                     site_info[dopant->site_id[isite]].charge,
         !strcmp(site_info[dopant->site_id[isite]].name,"") ? site_info[(species+id)->site_id[0]].name:
                     site_info[dopant->site_id[isite]].name);
      }
   }
   (void)printf("end\n");

/* Write potential parameters for pairs of site_ids */

   (void)printf("%s\n",potspec[system->ptype].name);
   for(idi = 1; idi < max_id; idi++)
   {
      for(idj = idi; idj < max_id; idj++)
      {
         idij = idj + idi*max_id;
         (void)printf("%5d %5d", idi, idj);
         for(ip = 0; ip < n_potpar; ip++)
            (void)printf(" %10g",potpar[idij].p[ip]);
         (void)putchar('\n');
      }
   }
   (void)printf("end\n");

/* Now we write the box dimensions */
   (void)printf("%g  %g  %g  %g  %g  %g  1  1  1\n",
          a,b,c,alpha,beta,gamma);

/* Followed by the molecules' centres of mass */
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      ipos = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
         specname = spec->name;
         if( spec->quat != NULL )   /* If polyatomic, record quaternions */
            for( i=0; i < 4; i++)
               quaternion[i] = spec->quat[imol][i];
         else
 	    quat_valid = false;

         if( molname != NULL && !strcmp(strlower(spec->name), molname) && dopant->nmols > 0 ) /* Species being replaced */
         {
            if( ipos < dopant->nmols && positions[ipos] == imol )  /* Match species position with list of substituted positions */
            {
               specname = dopant->name;
               
               if( quat_valid)
               {
                  q_to_rot(quaternion, rot_mat);
                  mat_vec_mul(rot_mat, (vec_mt*)spec->p_f_sites[0], (vec_mt*)spec_origin, 1);
               }
               else
               {
                  for(i = 0; i < 3; i++)
                     spec_origin[i] = spec->p_f_sites[0][i];
               }
               mat_vec_mul(hinv, (vec_mt*)spec_origin, (vec_mt*)spec_origin, 1);    /* Convert to fractional coords */

               if( dopant->nsites > 1 && ! spec->framework )
               {
                  if( euler[0]+euler[1]+euler[2] != 3e6)   /* Use euler angles if specified */
                  {
                     quaternion[0] = cos(0.5*euler[1])*cos(0.5*(euler[0]+euler[2]));
                     quaternion[1] = sin(0.5*euler[1])*cos(0.5*(euler[0]-euler[2]));
                     quaternion[2] = sin(0.5*euler[1])*sin(0.5*(euler[0]-euler[2]));
                     quaternion[3] = cos(0.5*euler[1])*sin(0.5*(euler[0]+euler[2]));
                  }
                  else if( spec->quat == NULL || dopant->rdof > 0)
                     quat_gen(quaternion); /* Generate quaternions for poly replacing monoatomic species */
		  quat_valid = true;
               }

               if( site_info[0].pad )  /* Place dopant pf origin at solvent pf origin */
               {
                  if( quat_valid )
                  {
                     q_to_rot(quaternion, rot_mat);
                     mat_vec_mul(rot_mat, (vec_mt*)pf_cofm, (vec_mt*)pf_origin, 1);
                  }
                  else
                  {
                     for(i = 0; i < 3; i++)
                        pf_origin[i] = pf_cofm[i];
                  }
               
                  mat_vec_mul(hinv, (vec_mt*)pf_origin, (vec_mt*)pf_origin, 1);    /* Convert to fractional coords */

                  for(i = 0; i < 3; i++)
                     spec->c_of_m[imol][i] += spec_origin[i]+pf_origin[i];  /* Shift c_of_m to dopant's pf origin */
               }
               if( dopant->nsites == 1 )
		  quat_valid = false;

               ipos++;
            }
         }

         (void)printf("%-*s  ", namelength,specname);          /* Write species name and coords */
         for( i = 0; i < 3; i++)
           (void)printf("%9g ",
              spec->c_of_m[imol][i]+0.5 - floor(spec->c_of_m[imol][i]+0.5));

         if( quat_valid )                       /* Write quaternions if polyatomic */
            (void)printf("%9g %9g %9g %9g",quaternion[0],quaternion[1],
                        quaternion[2],quaternion[3]);
         (void)putchar('\n');
      }
   }
   (void)printf("end\n");

   afree((gptr*) site);
   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
}
/******************************************************************************
 * random_pos.  Choose positions to be replaced randomly                      *
 ******************************************************************************/
void
random_pos(int totmol, int submol, int *subpos)
{
   int		ranpos, subflag;
   int		i, j;

   srand(time(NULL)+rand());   /* Generate random number */

   for( i = 0; i < submol; i++ )
   {   
       do
       {
          subflag = 0;
          ranpos = rand() % totmol;
          for ( j = 0; j < i; j++ ) 
             if( ranpos == subpos[j] )
                 subflag++;
       }
       while (subflag);

       subpos[i] = ranpos;
   }
   sort_pos(subpos, submol);
   return;
}
/******************************************************************************
 * create_total_sites. Combine site arrays into single array.                 *
 ******************************************************************************/
int
create_total_sites(int max_id, int new_sites, site_mt site_info[], site_mt dopant_sites[], site_mt totsites[])
{
   int       i;

   for( i=0; i < max_id; i++)
   {
      strcpy(totsites[i].name, site_info[i].name);
      totsites[i].mass = site_info[i].mass;
      totsites[i].charge = site_info[i].charge;
      totsites[i].pad = site_info[i].pad;
   }

   if( new_sites > 0 )
      for( i=0; i < new_sites; i++)
      {
         strcpy((totsites+max_id+i)->name, dopant_sites[i].name);
         totsites[max_id+i].mass = dopant_sites[i].mass;
         totsites[max_id+i].charge = dopant_sites[i].charge;
         totsites[max_id+i].pad = dopant_sites[i].pad;
      }

   return 0;
}
/******************************************************************************
 * copy_pot. Copy potential array to new array (possibly of greater dimension)*
 ******************************************************************************/
void
copy_pot(pot_mt *new_pot, pot_mt *old_pot, int max_id, int new_sites)
{
   int     i, j, k;
   int     idi, idj, idij_new, idij_old;

   for( i = 0; i < max_id; i++)
      for( j = 0; j < max_id; j++)
      {
         idij_old = j + i*max_id;
         idij_new = j + i*(max_id+new_sites);
         new_pot[idij_new].flag = old_pot[idij_old].flag;

         for( k=0; k < NPOTP; k++)
            new_pot[idij_new].p[k] = old_pot[idij_old].p[k];
      }
   return;
}
/******************************************************************************
 * main().   Driver program for substituting species in MOLDY sys_spec files  *
 * Acceptable inputs are sys-spec files, restart files or dump files.         *
 * Call: ransub [-s sys-spec-file] [-r restart-file].                         *
 * If neither specified on command line, user is interrogated.                *
 ******************************************************************************/
int
main(int argc, char **argv)
{
   int  c, cflg = 0, ans_i, sym, data_source = 0;
   char         line[80];
   extern char  *optarg;
   int          errflg = 0;
   int          intyp = 0;
   int          start, finish, inc;
   int          rflag, mflag, uflag;
   int          i,j,k,irec;
   int          n_elem = -1;       /* No of records read from element data file */
   char         *filename = NULL, *dump_name = NULL;
   char         *dumplims = NULL;
   char         *molname = NULL;
   char         *elename = "elements.dat";
   char         *potname = NULL;
   char         *dopfile = NULL;
   char         elefile[PATHLN] = "";
   char         potfile[PATHLN] = "";
   char         *tempname;
   char         dumpcommand[256];
   int          dump_size;
   float        *dump_buf;
   FILE         *Fp, *Dp;
   restrt_mt    restart_header;
   system_mt    sys;
   spec_mt      *species, *spec, *new_species;
   site_mt      *site_info;
   pot_mt       *potpar;
   pot_mp       new_pot;
   quat_mt      *qpf;
   contr_mt     control_junk;
   int          av_convert;
   int          maxmol;
   spec_data    element[NELEM];
   spec_data    dopant = {"","",-1.0, 1e6};
   spec_mt      dopspec;
   site_mt      *dopsite, *totsite;
   int          ndopsites = -1, *pos;
   vec_mt       p_f_coords[MAX_ATOMS];
   double       charge[MAX_ATOMS];
   mat_mt       h;
   char         spgr[16];                     /* Space Group in Herman Maugain form */
   char         label[MAX_ATOMS][NLEN];       /* Site name array */
   char         title[TITLE_SIZE];
   int          newsites=0;
   int          insw = -1;
   double       simbox[3];
   real		euler[3];
   boolean      strict_match = OFF;
   boolean      cofm_shift = OFF;
   boolean      new_quat = OFF;

#define MAXTRY 100
   control.page_length=1000000;
   dopspec.nmols = dopspec.rdof = 0;
   zero_real(h[0],9);
   zero_real(charge,MAX_ATOMS);
   euler[0] = euler[1] = euler[2] = 1e6;
   strcpy(spgr,"P 1");

   comm = argv[0];
   if( strstr(comm, "ranquat") )
     new_quat = ON;

   while( (c = getopt(argc, argv, "cr:s:m:n:u:o:w:q:z:e:y:a:hxkjf:t:p:") ) != EOF )
      switch(c)
      {
       case 'c':
         cflg++;
         break;
       case 'r':
       case 's':
	 if( intyp )
	    errflg++;
	 intyp = data_source = c;
	 filename = optarg;
	 break;
       case 'm':
	 molname = strlower(mystrdup(optarg));
	 break;
       case 'n':
         dopspec.nmols = atoi(optarg);
	 break;
       case 'u':
	 strncpy(dopant.name, strlower(optarg),NLEN);
	 break;
       case 'w':
         dopant.mass = atof(optarg);
	 break;
       case 'q':
         dopant.charge = atof(optarg);
	 break;
       case 'z':
         strncpy(dopant.symbol, optarg,4);
	 break;
       case 'e':
         elename = optarg;
         break;
       case 'y':
         potname = optarg;
         /* Create full path name, but don't exceed max length of string */
         strncat(strncat(potfile, POTPATH, PATHLN-strlen(potfile)), potname, PATHLN-strlen(potfile));
         break;
       case 'a':
         dopfile = optarg;   /* Structure file for (polyatomic) dopant species */
         break;
       case 'x':
         strict_match = ON;  /* Match atoms/ions by name and charge rather than name only */
         break;
       case 'h':
         cofm_shift = ON;    /* Position molecules using first sites rather than COFMs */
         break;
       case 'k':
         dopspec.rdof = 1;   /* Generate new quaternions for all dopant molecules */
         break;
       case 'j':
         new_quat = ON;   /* Generate new quaternions for all polyatomic molecules */
         dopspec.rdof = 0;    /* New quaternions for dopant molecules redundant */
         break;
       case 'f':
         euler[0] = atof(optarg);   /* Euler angle 'phi' */
         break;
       case 't':
         euler[1] = atof(optarg);   /* Euler angle 'theta' */
         break;
       case 'p':
         euler[2] = atof(optarg);   /* Euler angle 'psi' */
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
      fprintf(stderr,
             "Usage: %s [-r restart-file | -s sys-spec-file] ",comm);
      fputs("[-c] [-m solvent-species] ",stderr);
      fputs("[-u solute-species] [-n no-of-substitutions] [-w mass] [-q charge] [-z symbol] ",stderr);
      fputs("[-a solute-structure-file] [-f] [-t] [-p] ",stderr);
      fputs("[-x] [-h] [-k] [-o output-file]\n",stderr);
      exit(2);
   }

   if(intyp == 0)
   {
      fputs("How do you want to specify the simulated system?\n", stderr);
      fputs("Do you want to use a system specification file (1)", stderr);
      fputs(" or a restart file (2)", stderr);
      if( (ans_i = get_int("? ", 1, 2)) == EOF )
	 exit(2);
      intyp = ans_i-1 ? 'r': 's';
      if( intyp == 's' )
      {
	 fputs( "Do you need to skip 'control' information?\n", stderr);
	 if( (sym = get_sym("y or n? ","yYnN")) == 'y' || sym == 'Y')
	    cflg++;
      }

      if( (filename = get_str("File name? ")) == NULL )
	 exit(2);
   }

   switch(intyp)
   {
    case 's':
      if( (Fp = fopen(filename,"r")) == NULL)
	 error("Couldn't open sys-spec file \"%s\" for reading", filename);
      if( cflg )
      {
	 do
	 {
	    fscanf(Fp, "%132s",line);
	    (void)strlower(line);
	 }
	 while(! feof(stdin) && strcmp(line,"end"));
      }
      read_sysdef(Fp, &sys, &species, &site_info, &potpar);
      qpf = qalloc(sys.nspecies);
      initialise_sysdef(&sys, species, site_info, qpf);
      break;
    case 'r':
      if( (Fp = fopen(filename,"rb")) == NULL)
	 error("Couldn't open restart file \"%s\" for reading -\n%s\n", 
	       filename, strerror(errno));
      re_re_header(Fp, &restart_header, &control_junk);
      re_re_sysdef(Fp, restart_header.vsn, &sys, &species, &site_info, &potpar);
      prep_pot(&sys, site_info, potpar);
      break;
    default:
      error("Internal error - invalid input type", "");
   }
   allocate_dynamics(&sys, species);
   maxmol = sys.nmols;

  /*
   * Request substituting species if not already provided
   */
   do
   {
      uflag = 0;
      if( !strcmp(dopant.name,"") && (dopspec.nmols > 0 || ( molname != NULL && !new_quat)) )
      {
           fputs("What is the name of the substituting species ",stderr);
           strncpy(dopant.name, get_str("? "),NLEN);
      }
      else
         uflag++;
   } while (!uflag);

  /*
   * Request species to be replaced if not already provided
   */
   do
   {
      mflag = 0;
      if( molname == NULL)
         if( strcmp(dopant.name,"") || dopspec.nmols > 0)
         {
            fputs("What is the name of the species to be replaced",stderr);
            molname = get_str("? ");
         }
         else
            if( new_quat )
            {
               fputs("What is the name of the species to be rotated",stderr);
               molname = get_str("? ");
            }
            else
               mflag++;

      if( molname != NULL)
      {
         for(spec = species; spec < species+sys.nspecies; spec++)
            if( !strcmp(strlower(spec->name), molname) )
         {
             maxmol = spec->nmols;
             mflag++;
         }
         if(!mflag)
         {
             fprintf(stderr,"Species \"%s\" cannot be found\n", molname);
             (void)free(molname);
             molname = NULL;
         }
      }
   } while (!mflag);

   if( (molname != NULL || strcmp(dopant.name,"")) )
   {
       /* Halt program if names of solute and solvent species identical */
       if( !strcmp(molname,dopant.name) )
          error("Solute and solvent species identical - substitution aborted\n");

       while( dopspec.nmols < 0 )
       {
          fprintf(stderr, "How many %s species do you want to replace", molname);
          dopspec.nmols = get_int("? ",0,maxmol);
       }
   }

   if( dopspec.nmols < 0 )
      dopspec.nmols = 0;

   if( dopspec.nmols > maxmol )
      dopspec.nmols = maxmol;

   if( dopspec.nmols == 0)
      pos = ialloc(1);
   else
      pos = ialloc(dopspec.nmols);

   data_source = intyp;

   /* Create full path name, but don't exceed max length of string */
   strncat(strncat(elefile, ELEPATH, PATHLN-strlen(elefile)), elename, PATHLN-strlen(elefile));

   /* Read dopant species data from element data file if available */
   n_elem = read_ele(element, elefile);
   if(n_elem < 0 )
      message(NULLI,NULLP,WARNING,NOELEM,elefile);
  
   if( dopspec.nmols > 0 )
   {
      /* If dopant molecular structure file specified, read in data */
      if( dopfile != NULL )
      {
         if( !strncmp(strlower(read_ftype(dopfile)),"pdb",3) )
            insw = PDB;
         else
            if( !strncmp(strlower(read_ftype(dopfile)),"cssr",4) )
               insw = CSSR;
            else
               if( !strncmp(strlower(read_ftype(dopfile)),"shak",4) )
                  insw = SHAK;
               else
                  if( !strncmp(strlower(read_ftype(dopfile)),"xtl",3) )
                     insw = XTL;
                  else
                     if( !strncmp(strlower(read_ftype(dopfile)),"xyz",3) )
                        insw = XYZ;

         switch(insw)  /* Read in data according to format selected */
         {
            case PDB:
               ndopsites = read_pdb(dopfile, h, label, p_f_coords, charge, title, spgr);
               break;
            case CSSR:
               ndopsites = read_cssr(dopfile, h, label, p_f_coords, charge, title, spgr);
               break;
            case SHAK:
               ndopsites = read_shak(dopfile, h, label, p_f_coords, charge, title, simbox);
               break;
            case XTL:
               ndopsites = read_xtl(dopfile, h, label, p_f_coords, charge, title, spgr);
               break;
            case XYZ:
               ndopsites = read_xyz(dopfile, h, label, p_f_coords, title);
               break;
            default:
               error("Configuration file \"%s\" of unknown format", dopfile);
         }

         if( ndopsites < 0)
            message(NULLI,NULLP,WARNING,NOSUB,dopfile);

         if( strcmp(spgr,"P 1")  && ndopsites > 0 )
            ndopsites = sgexpand(MAX_ATOMS, ndopsites, p_f_coords, label, charge, spgr);

         if( det(h) != 0 )
             mat_vec_mul(h, p_f_coords, p_f_coords, ndopsites);   /* Convert to Cartesian coords */
      }
      else
      {
         ndopsites = 1;
         if( !strcmp(dopant.symbol,"" ) )
            strcpy(label[0],dopant.name);
         else
            strcpy(label[0],dopant.symbol);
         p_f_coords[0][0] = p_f_coords[0][1] = p_f_coords[0][2] = 0.0;
      }

      dopsite = (site_mt*)arralloc(sizeof(site_mt),1,0,ndopsites-1);
      if( cofm_shift )
         site_info[0].pad = ON;
      else
         site_info[0].pad = OFF;
      strcpy(dopspec.name, dopant.name);
      dopspec.nsites = ndopsites;
      dopspec.site_id = ialloc(ndopsites+1);

      dopspec.p_f_sites = ralloc(ndopsites);
      dopspec.c_of_m = NULL;
      dopspec.quat = NULL;
      dopspec.mass = 0.0;

      /* Check Euler angles for polyatomic dopants */
      if( ndopsites > 1)
        if( euler[0] != 1e6 || euler[1] != 1e6 || euler[2] != 1e6)
           for(i = 0; i < 3; i++)
           {
              if( euler[i] == 1e6)
                 euler[i] = 0.0;
	      else
                 euler[i] = fmod(euler[i],360.0);
              euler[i] *= PI/180.0;       /* Convert to radians */
           }

      for( i=0; i < ndopsites; i++) /* Add polyatomic info to dopant species and site_info */
      {
         dopspec.site_id[i] = 0;   /* Initialize site id */

         for( k=0; k < 3; k++)
            if( fabs(p_f_coords[i][k]) >= 1e-6)
               dopspec.p_f_sites[i][k] = p_f_coords[i][k];  /* Site coords */
            else
               dopspec.p_f_sites[i][k] = 0.0;

         /* Compare sites in new molecule with existing sites */
         for( j=0; j < sys.max_id; j++)
         {
            if( strict_match )   /* Match both atom name and charge if strict match selected */
            {
               if( !strcmp(label[i],(site_info+j)->name) && (fabs((site_info+j)->charge-charge[i]) < 1e-10) )
               {
                  dopspec.site_id[i] = j;
                  break;
               }
            }
            else
            {
               if( !strcmp(label[i],(site_info+j)->name) ) /* Match atoms to name only */
               {
                  dopspec.site_id[i] = j;
                  break;
               }
            }
         }

         /* Compare sites in new molecule with other new sites */
         for( j=0; j < newsites; j++)
         {
            if( strict_match )   /* Match both atom name and charge if strict match selected */
            {
               if( !strcmp(label[i],(dopsite+j)->name) && (fabs(charge[j]-charge[i]) < 1e-10) )
               {
                  dopspec.site_id[i] = sys.max_id+j;
                  break;
               }
            }
            else
            {
               if( !strcmp(label[i],(dopsite+j)->name) )
               {
                  dopspec.site_id[i] = sys.max_id+j;
                  break;
               }
            }
         }
         if( ndopsites == 1 )  /* If monatomic, compare with existing species */
           for( j=0; j < sys.nspecies; j++ )
              if( !strcmp(label[i],(species+j)->name) )
                  dopspec.site_id[i] = (species+j)->site_id[0];

         /* If not found, add new site to site array */
         if( (dopspec.site_id[i] == 0) && strncmp("#", dopspec.name,1) )
         {
            dopspec.site_id[i] = sys.max_id+newsites;
            if( ndopsites == 1 ) /* If monatomic species, assign user specified properties */
            {
               dopsite->mass = dopant.mass;
               dopsite->charge = dopant.charge;
               strncpy(dopsite->name, dopant.symbol, 4);
            }
            else                /* Otherwise set to dummy values */
            {
               (dopsite+newsites)->mass = -1;
               (dopsite+newsites)->charge = charge[i];
               strcpy((dopsite+newsites)->name, "");
            }

            if( label[i] != NULL )
               strncpy((dopsite+newsites)->name, label[i], 4);

            for( irec=0; irec < n_elem; irec++)
            {
               if( (ndopsites == 1) && (!strcmp(strlower(dopspec.name),strlower((element+irec)->name))) )
                  strncpy(dopsite->name, (element+irec)->symbol, 4);
               if( !strcmp((dopsite+newsites)->name,(element+irec)->symbol) )
               {
                  if( (dopsite+newsites)->mass < 0)
                     (dopsite+newsites)->mass = (element+irec)->mass;
                  if( (dopsite+newsites)->charge == 1e6 || !strict_match )
                     (dopsite+newsites)->charge = (element+irec)->charge;
               }
            }
            newsites++;
         }
      }
   }

   new_pot = aalloc(SQR(sys.max_id+newsites), pot_mt );

   copy_pot(new_pot,potpar,sys.max_id, newsites);  /* Create enlarged potential parm array to accommodate new species */

   totsite = (site_mt*)arralloc(sizeof(site_mt),1,0,sys.max_id+newsites-1);

   create_total_sites(sys.max_id, newsites, site_info, dopsite, totsite); /* Combine dopant and system site arrays */

   /* If potential parameter file exists, read in data */
   if( (potname != NULL) && (read_pot(potfile, &new_pot, sys.max_id, newsites, totsite) < 0) )
        message(NULLI,NULLP,WARNING,NOPOTL,potfile);

   sys.max_id += newsites;

   switch(data_source)                  /* To read configurational data       */
   {
    case 's':                           /* Lattice_start file                 */
        lattice_start(Fp, &sys, species, qpf);
        if( new_quat && molname != NULL)
           ran_quat(&sys, species, molname);
        random_pos(maxmol, dopspec.nmols, pos);
        sys_spec_out(&sys, species, &dopspec, molname, pos, totsite, euler,  new_pot);
      break;
    case 'r':                           /* Restart file                       */
        init_averages(sys.nspecies, restart_header.vsn,
                      control_junk.roll_interval, control_junk.roll_interval,
                      &av_convert);
        read_restart(Fp, restart_header.vsn, &sys, av_convert);
        if( new_quat && molname != NULL)
           ran_quat(&sys, species, molname);
        random_pos(maxmol, dopspec.nmols, pos);
        sys_spec_out(&sys, species, &dopspec, molname, pos, totsite, euler, new_pot);
      break;
    default:
      break;
    }
   return 0;
}
