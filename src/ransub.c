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
 *  Revision 1.13  2002/06/17 08:58:12  kr
 *   Modified sys-spec output routine to add dopant (solute) species to end of list.
 *   Check added to see if dopant species already present in system.
 *   Sorting of substituted positions list for improved efficiency.
 *   If no potentials specified, assume same pots as first site of substituted species.
 *   Substitution aborted if solute and solvent species have same name.
 *   External data files now only read if no. of solute particles > 0
 *    Tidied up error messages.
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
static char *RCSid = "$Header: /usr/users/kr/CVS/moldy/src/ransub.c,v 1.13 2002/06/17 08:58:12 kr Exp $";
#endif  

#include "defs.h"
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include "stdlib.h"
#include "stddef.h"
#include "string.h"
#include "time.h"
#include <stdio.h>
#include "structs.h"
#include "messages.h"
#include "utlsup.h"
#include "specdata.h"

void	read_sysdef(FILE *file, system_mp system, spec_mp *spec_pp, site_mp *site_info, pot_mp *pot_ptr);
void	initialise_sysdef(system_mp system, spec_mt *species, site_mt *site_info, quat_mt (*qpf));
void	re_re_header(FILE *restart, restrt_mt *header, contr_mt *contr);
void	re_re_sysdef(FILE *restart, char *vsn, system_mp system, spec_mp *spec_ptr, site_mp *site_info, pot_mp *pot_ptr);
void	allocate_dynamics(system_mp system, spec_mt *species);
void	lattice_start(FILE *file, system_mp system, spec_mp species, quat_mt (*qpf));
void	read_restart(FILE *restart, char *vsn, system_mp system, int av_convert);
void	init_averages(int nspecies, char *vsn, long int roll_interval, long int old_roll_interval, int *av_convert);
void    conv_potentials(const unit_mt *unit_from, const unit_mt *unit_to, pot_mt *potpar, int npotpar, int ptype, site_mt *site_info, int max_id);
int	getopt(int, char *const *, const char *);
gptr	*talloc(int n, size_mt size, int line, char *file);
char    *atime(void);
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
/*========================== External data references ========================*/

extern  const pots_mt   potspec[];          /* Potential type specification */

/******************************************************************************
 * read_pot2().  Read potential data from file.                                *
 ******************************************************************************/

extern int read_ele (spec_data *element, char *filename);

int        read_pot2(char *potfile, pot_mp *pot_ptr, int idj, site_mt *site_info, spec_data *dopant)
                   
                                  /* To be pointed at potpar array      */
              
                      
                  
{
char       atom1[4], atom2[4];
double     chg1, chg2;
double     p_tmp;
pot_mt     pot;
int        i,idi, n_items;
char       name[L_name],             /* Temporary species name             */
           line[L_name],             /* Store for input line from file     */
           pline[L_name];            /* Used in pot'l paramater parsing    */
int        ptype=-1;               /* Potential type index               */
int        nerrs = 0;              /* Accumulated error count            */
site_mt    spi;
spec_data  *spj;

FILE       *Fpot;

    if( (Fpot = fopen(potfile,"r")) == NULL)
        return ptype;

    n_items = sscanf(get_line(line,L_name,Fpot), "%s", name);

    if( n_items <= 0 )
       message(NULLI,NULLP,FATAL,SYSEOF,"potential type specification");

    for(i = 0; potspec[i].name; i++)             /* Is 'name' a known type? */
       if(strcmp(strlower(name), potspec[i].name) == 0)
          break;

    if(! potspec[i].name)                        /* Did the loop find 'name'? */
       message(&nerrs,line,FATAL,UNKPOT,name);   /* no                        */
    ptype = i;                                   /* yes                       */

    spj = dopant;

    while(sscanf(get_line(line,L_name,Fpot),"%s",name) > 0
                 && strcmp(strlower(name), "end") != 0)
    {
       n_items = 0;
       if(sscanf(line,"%4s %lf %4s %lf %[^#]",atom1,&chg1,atom2,&chg2,pline) <= 2)
            message(&nerrs,line,ERROR,NOPAIR);
       else
       {
                                                /* Now read in parameters   */
          (void)strcat(pline, "$");          /*   Add marker to end      */
          while(n_items < NPOTP && sscanf(pline,"%lf %[^#]", &p_tmp, pline) > 1 )
              pot.p[n_items++] = p_tmp;
       }
       for( idi = 1; idi < idj; idi++)
       {
          spi = site_info[idi];
          if( ((!strcmp(atom1,spi.name)) && (chg1 == spi.charge) &&
              (!strcmp(atom2,spj->symbol)) && (chg2 == spj->charge)) ||
                 ((!strcmp(atom1,spj->symbol)) && (chg1 == spj->charge) &&
                    (!strcmp(atom2,spi.name)) && (chg2 == spi.charge)) )
          {
                          (*pot_ptr)[idi-1] = pot;  /* Write potential to "spare" pot values */
          }
       }
    }

    fclose(Fpot);

    if(nerrs > 0)                        /* if any errors have been detected */
       message(&nerrs,NULLP,FATAL,ERRS,nerrs,(nerrs>1)?'s':' ');

    return ptype;
}
/******************************************************************************
 * prep_pot().  Convert units and add potentials for dopant.                  *
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
 * sys_spec_out().  Write a system configuration to stdout in the form of a   *
 * system specification file for MOLDY                                        *
 ******************************************************************************/
void
sys_spec_out(system_mt *system, spec_mt *species, char *molname, spec_data *dopant, site_mt *site_info, pot_mt *potpar)
{
   spec_mt      *spec;
   double       a, b, c, alpha, beta, gamma;
   mat_mp       h = system->h;
   int          i=0, imol, isite, ipos=0;
   int		specmol, id = -1;
   int          idi, idj, idij, ip;
   int          n_potpar = system->n_potpar;
   char         *specname;
   int          max_id = system->max_id;
   int          dflag = 0;

   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

/* Write header for sys_spec file */
   (void)printf("# System specification file written by RANSUB on %s\n",atime());

/* Write site data for each molecule */
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      /* Check if species matches species to be substituted */
      if( (molname != NULL) && !strcmp(strlower(spec->name), molname) )
      {
          specmol = spec->nmols - dopant->nmols; /* Subtract number of substituting species */
          id = i;
      }
      else
         specmol = spec->nmols;

      /* Check if dopant species already present in system */
      if( (dopant->name != NULL) && !strcmp(strlower(spec->name),dopant->name) )
      {
          specmol += dopant->nmols;
          dflag++;
      }

      i++;

      (void)printf("%s  %d  %s\n", spec->name, specmol,
                    spec->framework ? "framework" : "");
      for(isite=0; isite < spec->nsites; isite++)
         (void)printf("%d %9g %9g %9g %9g %9g %s\n",
                        spec->site_id[isite],
                        spec->p_f_sites[isite][0],
                        spec->p_f_sites[isite][1],
                        spec->p_f_sites[isite][2],
                        site_info[spec->site_id[isite]].mass,
                        site_info[spec->site_id[isite]].charge,
                        site_info[spec->site_id[isite]].name);
   }
   if( !dflag && id >= 0 )
   {
      (void)printf("%s  %d  %s\n", dopant->name, dopant->nmols,
              (species+id)->framework ? "framework" : "");
      (void)printf("%d %9g %9g %9g %9g %9g %s\n",
               max_id,
               (species+id)->p_f_sites[0][0],
               (species+id)->p_f_sites[0][1],
               (species+id)->p_f_sites[0][2],
      dopant->mass < 0 ? site_info[(species+id)->site_id[0]].mass:dopant->mass,
      dopant->charge == 1e6 ? site_info[(species+id)->site_id[0]].charge:dopant->charge,
      !strcmp(dopant->symbol,"") ? site_info[(species+id)->site_id[0]].name:dopant->symbol);
   }
   (void)printf("end\n");

/* Write potential parameters for pairs of site_ids */

   if( !dflag && dopant->nmols > 0 )
       max_id++;

   (void)printf("%s\n",potspec[system->ptype].name);
   for(idi = 1; idi < max_id; idi++)
   {
      for(idj = idi; idj < system->max_id; idj++)
      {
         idij = idj + idi*system->max_id;
         (void)printf("%5d %5d", idi, idj);
         for(ip = 0; ip < n_potpar; ip++)
            (void)printf(" %10g",potpar[idij].p[ip]);
         (void)putchar('\n');
      }
      if( !dflag && dopant->nmols > 0 )
      {
         if( id >=0 && !strcmp(dopant->symbol,"") )
         {
            if( idi < system->max_id)
               idij = (species+id)->site_id[0]+system->max_id*idi;
            else
               idij = (species+id)->site_id[0]+system->max_id*
                                         (species+id)->site_id[0];
         }
         potpar[idi-1] = potpar[idij];

         (void)printf("%5d %5d", idi, idj);
         for(ip = 0; ip < n_potpar; ip++)
            (void)printf(" %10g",potpar[idi-1].p[ip]);
         (void)putchar('\n');
      }
   }
   (void)printf("end\n");

/* Now we write the box dimensions */
   (void)printf("%g  %g  %g  %g  %g  %g  1  1  1\n",
          a,b,c,alpha,beta,gamma);

/* Followed by the molecules' centre of mass positions */
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      for(imol = 0; imol < spec->nmols; imol++)
      {
        specname = spec->name;
        if( molname != NULL)
           if( !strcmp(strlower(spec->name), molname))
              if( dopant->pos[ipos] == imol )
              {
                 specname = dopant->name;
                 ipos++;
              }
        (void)printf("%-*s  ", NLEN,specname);
        for( i = 0; i < 3; i++)
          (void)printf("%9g ",
             spec->c_of_m[imol][i]+0.5 - floor(spec->c_of_m[imol][i]+0.5));
        if(spec->quat != NULL)
           (void)printf("%9g %9g %9g %9g",spec->quat[imol][0],spec->quat[imol][1],
                 spec->quat[imol][2],spec->quat[imol][3]);
        (void)putchar('\n');
      }
   }
   (void)printf("end\n");

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
   int	c, cflg = 0, ans_i, sym, data_source = 0;
   char 	line[LLEN];
   extern char	*optarg;
   int		errflg = 0;
   int		intyp = 0;
   int		start, finish, inc;
   int		rflag, mflag, uflag;
   int		irec;
   int          n_elem;         /* No of records read from element data file */
   char		*filename = NULL, *dump_name = NULL;
   char		*dumplims = NULL;
   char		*molname = NULL;
   char         *elename = "elements.dat";
   char         *potname = NULL;
   char         elefile[50];
   char         potfile[50];
   FILE		*Fp;
   restrt_mt	restart_header;
   system_mt	sys;
   spec_mt	*species, *spec;
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf;
   contr_mt	control_junk;
   int          av_convert;
   int          maxmol;
   spec_data    element[NELEM];
   spec_data    dopant = {"","",-1.0, 1e6, 0,-1};

#define MAXTRY 100
   control.page_length=1000000;

   comm = argv[0];

   while( (c = getopt(argc, argv, "cr:s:m:n:u:o:w:q:z:e:y:") ) != EOF )
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
	 dopant.nmols = atoi(optarg);
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
      fputs("Usage: ransub [-r restart-file | -s sys-spec-file] ",stderr);
      fputs("[-c] [-m solvent-species] ",stderr);
      fputs("[-u solute-species] [-n no-of-substitutions] ",stderr);
      fputs("[-w mass] [-q charge] [-z symbol] [-o output-file]\n",stderr);
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
      if( !strcmp(dopant.name,"") && (dopant.nmols > 0 || molname != NULL) )
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
      {
         if( strcmp(dopant.name,"") || dopant.nmols > 0)
         {
            fputs("What is the name of the species to be replaced",stderr);
            molname = get_str("? ");
         }
         else
            mflag++;
      }

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
             message(NULLI,NULLP,WARNING,"Species \"%s\" not found\n", molname);
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

      if( dopant.nmols <= 0 )
      {
        fprintf(stderr, "How many %s species do you want to replace", molname);
	dopant.nmols = get_int("? ",0,maxmol);
      }
   }

   if( dopant.nmols < 0 )
        dopant.nmols = 0;

   if( dopant.nmols > maxmol )
      dopant.nmols = maxmol;

   (&dopant)->pos = ialloc(dopant.nmols);         

   data_source = intyp;

   if( dopant.nmols > 0 )
   {
      /* Create full path name, but don't exceed max length of string */
      strncat(strncat(elefile, ELEPATH, PATHLN-strlen(elefile)), elename, PATHLN-strlen(elefile));

      /* Read dopant species data from element data file if available */
      if( (n_elem = read_ele(element, elefile)) )
      {
         for( irec=0; irec < n_elem; irec++)
           if( !strcmp(strlower(dopant.name),strlower((element+irec)->name)) )
           {
              if( !strcmp(dopant.symbol, "") )
                 strcpy(dopant.symbol, (element+irec)->symbol);
              if( dopant.mass < 0)
                 dopant.mass = (element+irec)->mass;
              if( dopant.charge == 1e6)
                 dopant.charge = (element+irec)->charge;
           }
      }
      else
         message(NULLI,NULLP,WARNING,NOELEM,elefile);

      /* If potential parameter file exists, read in data */
      if( (potname != NULL) && (read_pot2(potfile, &potpar, sys.max_id, site_info, &dopant) < 0) )
           message(NULLI,NULLP,WARNING,NOPOTL,potfile);
   }

   switch(data_source)                  /* To read configurational data       */
   {
    case 's':                           /* Lattice_start file                 */
        lattice_start(Fp, &sys, species, qpf);
        random_pos(maxmol, dopant.nmols, (&dopant)->pos);
        sort_pos((&dopant)->pos,dopant.nmols);
        sys_spec_out(&sys, species, molname, &dopant, site_info, potpar);
      break;
    case 'r':                           /* Restart file                       */
        init_averages(sys.nspecies, restart_header.vsn,
                      control_junk.roll_interval, control_junk.roll_interval,
                      &av_convert);
        read_restart(Fp, restart_header.vsn, &sys, av_convert);
        random_pos(maxmol, dopant.nmols, (&dopant)->pos);
        sort_pos((&dopant)->pos,dopant.nmols);
        sys_spec_out(&sys, species, molname, &dopant, site_info, potpar);
      break;
    default:
      break;
    }
   return 0;
}
