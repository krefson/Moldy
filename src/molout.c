#ifndef lint
static char *RCSid = "$Header: /usr/users/moldy/CVS/moldy/src/molout.c,v 1.12.8.5 2004/04/19 04:19:31 moldydv Exp $";
#endif

#include "defs.h"
#include <stdarg.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "structs.h"
#include "ReadDCD.h"
#include "utlsup.h"
gptr	*arralloc(size_mt,int,...); 	/* Array allocator		      */

void	invert(real (*a)[3], real (*b)[3]);
void	make_sites(real (*h)[3], vec_mp c_of_m_s, 
		   quat_mp quat, vec_mp p_f_sites, 
		   real **site, int nmols, int nsites, int pbc);
void    error(char *format, ...);
void    afree(gptr *p);
char    *atime(void);
/*======================== Global vars =======================================*/
extern contr_mt		control;
/******************************************************************************
 ******************************************************************************/
void 
mat_vec_mul3(real (*m)[3], real **vec, int number)
                                /* Number of vectors to be multiplied         */
                                /* Matrix                                     */
                                /* Input vector.  Output vector same as input */
{
   int i;
   register double        a0, a1, a2;
   
   for(i = 0; i < number; i++)
   {
      a0 = vec[0][i];  a1 = vec[1][i];  a2 = vec[2][i];
      
      vec[0][i] = m[0][0]*a0 + m[0][1]*a1 + m[0][2]*a2;
      vec[1][i] = m[1][0]*a0 + m[1][1]*a1 + m[1][2]*a2;
      vec[2][i] = m[2][0]*a0 + m[2][1]*a1 + m[2][2]*a2;
   }
}
/******************************************************************************
 * Centre_mass.  Shift system centre of mass to origin (in discrete steps),   *
 ******************************************************************************/
static void
centre_mass(spec_mt *species, int nspecies, real *c_of_m)
{
   double	mass;
   spec_mt	*spec;
   int		imol;
   vec_mt	*s_c_of_m;

   mass = c_of_m[0] = c_of_m[1] = c_of_m[2] = 0.0;
   for(spec = species; spec < species + nspecies; spec++ )
   {
      s_c_of_m = spec->c_of_m;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 c_of_m[0] += spec->mass*s_c_of_m[imol][0];
	 c_of_m[1] += spec->mass*s_c_of_m[imol][1];
	 c_of_m[2] += spec->mass*s_c_of_m[imol][2];
      }
      mass += spec->nmols*spec->mass;
   }

   c_of_m[0] /= mass;
   c_of_m[1] /= mass;
   c_of_m[2] /= mass;
   c_of_m[0] = floor(c_of_m[0]+0.5);
   c_of_m[1] = floor(c_of_m[1]+0.5);
   c_of_m[2] = floor(c_of_m[2]+0.5);
}
/******************************************************************************
 * Shift.  Translate all co-ordinates.					      *
 ******************************************************************************/
static
void	shift(vec_mt (*r), int nmols, real *s)
{
   int imol;
   for(imol = 0; imol < nmols; imol++)
   {
      r[imol][0] += 0.5;
      r[imol][1] += 0.5;
      r[imol][2] += 0.5;
   }
}
/******************************************************************************
 * schakal_out().  Write a system configuration to stdout in the form of an   *
 * input data file for the graphics program SCHAKAL88.			      *
 ******************************************************************************/
static void
schakal_out(system_mt *system, mat_mp h, spec_mt *species, site_mt *site_info, 
	    char *insert, int n)
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   spec_mt	*spec;
   double	a, b, c, alpha, beta, gamma;
   mat_mt	hinv;
   int		imol, isite, is;

   invert(h,hinv);

   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

   printf("CELL %f %f %f %f %f %f\n", a, b, c, alpha, beta, gamma);
   for(spec = species; spec < species+system->nspecies; spec++)
   {

      make_sites(system->h, spec->c_of_m, spec->quat, spec->p_f_sites,
                 site, spec->nmols, spec->nsites, MOLPBC);

      mat_vec_mul3(hinv, site, spec->nsites*spec->nmols);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 puts("MOL");
	 for(is = 0; is < spec->nsites; is++)
	 {
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
	       (void)printf("ATOM %-8s %7.4f %7.4f %7.4f\n",
			    site_info[spec->site_id[is]].name,
			    site[0][isite], site[1][isite], site[2][isite]);
	    isite++;
	 }
      }
   }

   if( insert != NULL)
      (void)printf("%s\n", insert);

   (void)printf("END %d\n",n);

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*) site);
}
/******************************************************************************
 * xtl_out().  Write a system configuration to stdout in the form of an       *
 * BIOSYM XTL file.							      *
 ******************************************************************************/
static void xtl_out(system_mt *system, mat_mp h, spec_mt *species, site_mt *site_info, char *insert, int intyp)
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   double	qconv;	/* Variable for converting charge from program units */
   spec_mt	*spec;
   double	a, b, c, alpha, beta, gamma;
   mat_mt	hinv;
   int		imol, isite, is;
   int		charge;

   if( intyp == 'r' )
      qconv = CONV_Q;
   else
      qconv = 1.0;

   invert(h,hinv);

   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

   printf("TITLE %s\n",control.title);
   puts("DIMENSION 3");
   printf("CELL \n%f %f %f %f %f %f\n", a, b, c, alpha, beta, gamma);
   printf("SYMMETRY  NUMBER 1  LABEL P1\n");
   printf("SYM MAT  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0 0.0000 0.0000 0.0000\n");
   printf("ATOMS\nNAME       X          Y          Z     CHARGE   TEMP    OCCUP   SCAT\n");
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      make_sites(system->h, spec->c_of_m, spec->quat, spec->p_f_sites,
                 site, spec->nmols, spec->nsites, MOLPBC);

      mat_vec_mul3(hinv, site, spec->nsites*spec->nmols);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
	 {
            charge = (int)abs(site_info[spec->site_id[is]].charge*qconv);
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
	       (void)printf("%-4s %10.5f %10.5f %10.5f %7.4f   0.0000  1.0000   %-2s%d%c\n",
			    site_info[spec->site_id[is]].name,
			    site[0][isite], site[1][isite], site[2][isite],
			    site_info[spec->site_id[is]].charge*qconv,
                            site_info[spec->site_id[is]].name, charge, (charge >=0 ? '+':'-'));
	    isite++;
	 }
      }
   }

   if( insert != NULL)
      (void)printf("%s\n", insert);

   (void)printf("EOF\n");

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*) site);
}
/******************************************************************************
 * pdb_out().  Write a system configuration to stdout in the form of a        *
 * Brookhaven Protein Data Bank (pdb) file                                    *
 ******************************************************************************/
static void
pdb_out(system_mt *system, mat_mp h, spec_mt *species, site_mt *site_info, 
	char *insert, int intyp)
{
   double	**site = (double**)arralloc(sizeof(double),2,
                                            0,2,0,system->nsites-1);
   mat_mt	hinv;
   spec_mt	*spec;
   double	a,b,c, alpha, beta, gamma;
   double	qconv;	/* Variable for converting charge from program units */ 
   char         atom_name[3];
   double	atom_charge;
   int		imol, isite, itot=1, ispec=1;
   int		i, is;

   invert(h,hinv);
   
   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

   if( intyp == 'r' )
      qconv = CONV_Q;
   else
      qconv = 1.0;

/* Write the pdb header */
   (void)printf("HEADER     %-40s%10s%4d\n", "Moldy output", atime(), 1);
   (void)printf("TITLE      %-60s\n", control.title);
   (void)printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1        \n",
          a,b,c,alpha,beta,gamma);

   for(i = 0; i < 3; i++)
       (void)printf("SCALE%d    %10.6f%10.6f%10.6f        0.00000\n",
           i+1, hinv[i][0], hinv[i][1], hinv[i][2]);

   for(spec = species; spec < species+system->nspecies; ispec++, spec++)
   {
     make_sites(h, spec->c_of_m, spec->quat, spec->p_f_sites,
                site, spec->nmols, spec->nsites, MOLPBC);

     isite = 0;
     for(imol = 0; imol < spec->nmols; imol++)
     {
       for(is = 0; is < spec->nsites; is++)
       {
	  strncpy(atom_name, site_info[spec->site_id[is]].name, 2); atom_name[2] = 0;
	 atom_charge = site_info[spec->site_id[is]].charge*qconv;
         if(fabs(site_info[spec->site_id[is]].mass) != 0)
         {
            (void)printf("HETATM%5d %2s%-2d NON A   1    %8.3f%8.3f%8.3f\
 %5.2f %5.2f          %2s",
			 itot, atom_name, ispec, site[0][isite], site[1][isite], site[2][isite],
			 1.0,0.0, atom_name);

            if( atom_charge != 0.0 )
               (void)printf("%1.0f%c\n", fabs(atom_charge),atom_charge<0?'-':'+');
            else
               (void)printf("\n");
         }
         isite++;
         itot++;
       }
     }
   }
   (void)printf("TER   %5d      NON A   1\n",itot);
   (void)printf("END\n");
   if( insert != NULL)
      (void)printf("%s\n", insert);

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
}
/******************************************************************************
 * xyz_out().  Write a system configuration to stdout in the form of an       *
 * input data file for the graphics program XYZ (rasmol -xyz file)	      *
 ******************************************************************************/
static void
xyz_out(system_mt *system, mat_mp h, spec_mt *species, site_mt *site_info, 
	char *insert)
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   spec_mt	*spec;
   int		imol, isite, is;

/* We count the number of atoms */
   isite=0;
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
         {
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
              isite++;
         }
      }
   }
/* Now we write the Xyz header */
   (void)printf("%d\n",isite);
/* It would be nice to have here the real title */
   (void)printf("%s\n",control.title[0]?control.title:"Moldy output");
   
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      make_sites(h, spec->c_of_m, spec->quat, spec->p_f_sites,
		 site, spec->nmols, spec->nsites, MOLPBC);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
	 {
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
	       (void)printf("%-8s %7.4f %7.4f %7.4f\n",
			    site_info[spec->site_id[is]].name,
			    site[0][isite], site[1][isite], site[2][isite]);
	    isite++;
	 }
      }
   }

   if( insert != NULL)
      (void)printf("%s\n", insert);

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*) site);
}
/******************************************************************************
 * dcd_out().  Write a system configuration to stdout in the form of a        *
 * DCD data file for the graphics program VMD                                 *
 ******************************************************************************/
static void
dcd_out(system_mt *system, mat_mp h, spec_mt *species, site_mt *site_info, 
	int n, int irec, int inc)
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   float	**sitef = (float**)arralloc(sizeof(float),2,
					    0,2,0,system->nsites-1);
   spec_mt	*spec;
   int		isite, is, i, imol, isitem;

   isitem=0; 
   for(spec = species; spec < species+system->nspecies; spec++)
   {
      make_sites(h, spec->c_of_m, spec->quat, spec->p_f_sites,
		 site, spec->nmols, spec->nsites, MOLPBC);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
	 {
	    if(fabs(site_info[spec->site_id[is]].mass) != 0)
	    {
	       for(i=0; i<3; i++)
		  sitef[i][isitem] = site[i][isite];
	       isitem++;
	    }
	    isite++;
	 }
      }
   }
/* On first call write the DCD header. Always write to stdout. */
   if( n == 0 )
      write_dcdheader(stdout, control.title, isitem, irec, 0, inc, control.step);
   
   write_dcdstep(stdout, isitem, sitef[0], sitef[1], sitef[2]);

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*)site);
   afree((gptr*)sitef);
}
/******************************************************************************
 * atoms_out().  Write a system configuration to stdout in the form of an     *
 * binary atomic co-ordinates.						      *
 ******************************************************************************/
static void
atoms_out(system_mt *system, mat_mp h, spec_mt *species)
{
   double	**site = (double**)arralloc(sizeof(double),2,
					    0,2,0,system->nsites-1);
   spec_mt	*spec;
   float	fsite[3];
   mat_mt	hinv;
   int		imol, isite, is,i;

   invert(h,hinv);

   for(spec = species; spec < species+system->nspecies; spec++)
   {
      make_sites(h, spec->c_of_m, spec->quat, spec->p_f_sites,
		 site, spec->nmols, spec->nsites, MOLPBC);

      mat_vec_mul3(hinv, site, spec->nsites*spec->nmols);

      isite = 0;
      for(imol = 0; imol < spec->nmols; imol++)
      {
	 for(is = 0; is < spec->nsites; is++)
	 {
	    for(i=0; i<3; i++)
	       fsite[i]=site[i][isite];
            fwrite((char*)fsite, sizeof fsite, 1, stdout);
	    isite++;
	 }
      }
   }
   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*) site);
}
/******************************************************************************
 * cssr_out().  Write a system configuration to stdout in the form of         *
 * SERC Daresbury Lab's Cambridge Structure Search and Retrieval (cssr) file  *
 ******************************************************************************/
static void
cssr_out(system_mt *system, mat_mp h, spec_mt *species, 
	 site_mt *site_info, char *insert, int intyp)
{
   double	**site = (double**)arralloc(sizeof(double),2,
                                            0,2,0,system->nsites-1);
   mat_mt       hinv;
   spec_mt	*spec;
   double	a,b,c, alpha, beta, gamma;
   double	qconv;	/* Variable for converting charge from program units */ 
   char         atomname[5];
   int		imol, isite, itot=1, ispec=1;
   int		divd=0, is;
   
   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

   invert(h,hinv);

   if( intyp == 'r' )
      qconv = CONV_Q;
   else
      qconv = 1.0;

/* We count the number of atoms */
   isite=0;
   for(spec = species; spec < species+system->nspecies; spec++)
      for(imol = 0; imol < spec->nmols; imol++)
         for(is = 0; is < spec->nsites; is++)
            if(fabs(site_info[spec->site_id[is]].mass) != 0)
                isite++;

/* Exit if no of atoms exceeds max allowed by CSSR format */
   if( isite > 9999 )
      error("Too many atoms (%d) for CSSR format. Process aborted.", isite);

/* Write the cssr header */
   (void)printf("%38c %7.3f %7.3f %7.3f\n",' ',a,b,c);
   (void)printf("%21c %7.3f %7.3f %7.3f    SPGR =  1 P 1\n",' ',alpha,beta,gamma);
   (void)printf("%4d   0 %60s\n", isite, control.title);

   if( insert != NULL)
      (void)printf("       %53s\n", insert);
   else
      (void)printf("\n");

   for(spec = species; spec < species+system->nspecies; ispec++, spec++)
   {
     make_sites(h, spec->c_of_m, spec->quat, spec->p_f_sites,
                site, spec->nmols, spec->nsites, MOLPBC);

     mat_vec_mul3(hinv, site, spec->nsites*spec->nmols);

     isite = 0;
     for(imol = 0; imol < spec->nmols; imol++)
     {
       for(is = 0; is < spec->nsites; is++)
       {
            
         strncpy( atomname, site_info[spec->site_id[is]].name, 4);
         divd = pow(10, 4-strlen(atomname));
         if( divd > 1 )
            sprintf(atomname,"%s%d",site_info[spec->site_id[is]].name,itot%divd);

         (void)printf("%4d %-4s  %9.5f %9.5f %9.5f ",
              itot, atomname, site[0][isite], site[1][isite], site[2][isite]);
         (void)printf("   0   0   0   0   0   0   0   0 %7.3f\n",
              site_info[spec->site_id[is]].charge*qconv);
         isite++;
         itot++;
       }
     }
   }

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
}
/******************************************************************************
 * arc_out().  Write multiple system configurations to stdout in the form of  *
 * an Insight II ARC data file.                                               *
 ******************************************************************************/
static void
arc_out(system_mt *system, mat_mp h, spec_mt *species, site_mt *site_info, int intyp, int n)
{
   double       **site = (double**)arralloc(sizeof(double),2,
                                            0,2,0,system->nsites-1);
   spec_mt      *spec;
   double       a,b,c, alpha, beta, gamma;
   double       qconv;  /* Variable for converting charge from program units */
   char         atom_name[5], name_charge[7], *elem_sym;
   double       atom_charge;
   int          imol, isite, divd;
   int          i, is, ispec=1;

   a = sqrt(SQR(h[0][0]) + SQR(h[1][0]) + SQR(h[2][0]));
   b = sqrt(SQR(h[0][1]) + SQR(h[1][1]) + SQR(h[2][1]));
   c = sqrt(SQR(h[0][2]) + SQR(h[1][2]) + SQR(h[2][2]));
   alpha = 180/PI*acos((h[0][1]*h[0][2]+h[1][1]*h[1][2]+h[2][1]*h[2][2])/b/c);
   beta  = 180/PI*acos((h[0][0]*h[0][2]+h[1][0]*h[1][2]+h[2][0]*h[2][2])/a/c);
   gamma = 180/PI*acos((h[0][0]*h[0][1]+h[1][0]*h[1][1]+h[2][0]*h[2][1])/a/b);

   if( intyp == 'r' )
      qconv = CONV_Q;
   else
      qconv = 1.0;

/* First we write the ARC header */
   if( n == 0)
      (void)printf("!BIOSYM archive 3\nPBC=ON\n");

   (void)printf("Frame %d\n",n);
   (void)printf("!DATE %10s\n", atime());

/* Next we write the cell size parameters and space group (P 1 by default) */
   (void)printf("PBC%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f (P1)\n",
          a,b,c,alpha,beta,gamma);

   for(spec = species; spec < species+system->nspecies; spec++)
   {
     make_sites(h, spec->c_of_m, spec->quat, spec->p_f_sites,
           site, spec->nmols, spec->nsites, MOLPBC);

     isite = 0;
     for(imol = 0; imol < spec->nmols; imol++)
     {
         for(is = 0; is < spec->nsites; is++)
         {
            elem_sym = site_info[spec->site_id[is]].name;
            atom_charge = site_info[spec->site_id[is]].charge*qconv;
            divd = pow(10, 5-strlen(elem_sym));
            if( divd > 1)
            {
               sprintf(atom_name,"%s%d",elem_sym,(imol+1)%divd);
               sprintf(name_charge,"%s%-5.0f",elem_sym,fabs(atom_charge));
            }
            else
            {
               strncpy(atom_name, site_info[spec->site_id[is]].name, 5);
               strncpy(name_charge, site_info[spec->site_id[is]].name, 5);
            }
            if(fabs(site_info[spec->site_id[is]].mass) != 0)
               (void)printf("%-5s %14.9f %14.9f %14.9f XXX  %-2d     %-7s %-2s %6.3f\n",
                    atom_name, site[0][isite], site[1][isite], site[2][isite],
                    ispec, name_charge, elem_sym, atom_charge);
            isite++;
         }
     }
   }
   (void)printf("end\nend\n");

   if( ferror(stdout) )
      error("Error writing output - \n%s\n", strerror(errno));
   afree((gptr*) site);
}
/******************************************************************************
 * moldy_out.  Select output routine and handle file open/close.	      *
 * Translate system relative to either centre of mass of posn of framework.   *
 ******************************************************************************/
void
moldy_out(int n, int irec, int inc, system_mt *system, mat_mp h, 
	  spec_mt *species, site_mt *site_info, int outsw, int intyp, char *insert)
{
   spec_mp	spec, frame_spec  = NULL;
   vec_mt	c_of_m;
   
   for(spec = species; spec < species+system->nspecies; spec++)
      if( spec->framework )
	 frame_spec = spec;

   if( frame_spec != NULL )
      shift(system->c_of_m, system->nmols, frame_spec->c_of_m[0]);
   else
   {
      centre_mass(species, system->nspecies, c_of_m);
      shift(system->c_of_m, system->nmols, c_of_m);
   }

   switch (outsw)
   {
    case CSSR:
      cssr_out(system, h, species, site_info, insert, intyp);
      break;
    case DCD:
      dcd_out(system, h, species, site_info, n, irec, inc);
      break;
    case SHAK:
      schakal_out(system, h, species, site_info, insert, n);
      break;
    case PDB:
      pdb_out(system, h, species, site_info, insert, intyp);
      break;
    case XYZ:
      xyz_out(system, h, species, site_info, insert);
      break;
    case ARC:
      arc_out(system, h, species, site_info, intyp, n);
      break;
    case XTL:
      xtl_out(system, h, species, site_info, insert, intyp);
      break;
    case OUTBIN:
      atoms_out(system, h, species);
      break;
   }
}
