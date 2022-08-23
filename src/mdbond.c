/* MDBOND - Program for calculating bond lengths and angles
Copyright (C) 1999 Craig Fisher
For use with dump files from MOLecular DYnamics simulation code, Moldy,
by Keith Refson
 
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
 * mdbond    	Code for calculating intermolecular distances i-j and angles j/i\k    *       
 *              Determines all distances within bond cutoffs between c_of_m's         *
 *              (including mirror images) and angles within angle cutoffs for pairs   *
 *              of molecules within same bond cutoffs                                 *
 ************************************************************************************** 
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
#include "list.h"
#include "utlsup.h"

void    make_sites(real (*h)[3], vec_mp c_of_m_s,
                   quat_mp quat, vec_mp p_f_sites,
                   real **site, int nmols, int nsites, int pbc);

/*======================== Global variables ==================================*/
int ithread=0, nthreads=1;

/*
 * Structures for bond and angle data. 
 */
typedef struct
{
   char 	atom1[3];
   char 	atom2[3];
   int 		number1;
   int	 	number2;
   double	length;
} BOND;

typedef struct
{
   char 	atom1[3];
   char 	atom2[3];
   char		atom3[3];
   int 		number1;
   int	 	number2;
   int	 	number3;
   double	length1;
   double	length2;
   double	value;
} ANGLE;

/******************************************************************************
 * morethan_BOND(). Compare distances stored in BOND structure types          *
 ******************************************************************************/
NODE  *morethan_BOND(ROOT **root, BOND *data)
{
   NODE *node = NULL;
   BOND	*bd;

   if( VALID(root))
   {
       node = (*root)->head;
       bd = node->data;

       while( (node != NULL) && (bd->length < data->length) )
       {
           node = node->next;
           if( node != NULL )
              bd = node->data;
       }
   }
   return (node);
}
/******************************************************************************
 * morethan_ANGLE(). Compare angles stored in ANGLE structure types           *
 ******************************************************************************/
NODE *morethan_ANGLE(ROOT **root, ANGLE *data)
{
   NODE *node = NULL;
   ANGLE *ang;

   if( VALID(root))
   {
       node = (*root)->head;
       ang = node->data;

       while((node != NULL) && (ang->value < data->value) )
       {
           node = node->next;
           if( node != NULL )
              ang = node->data;
       }
   }
   return (node);
}
/******************************************************************************
 * angle_calc(). Calculate angle (in degrees) between two vectors             *
 ******************************************************************************/
double angle_calc(real *vec1, real *vec2)
{
double    dp, angle;                    /* Dot product and angle     */
double	  a2, b2;			/* Distances a and b */

   a2 = DOTPROD(vec1,vec1);
   b2 = DOTPROD(vec2,vec2);
   dp = DOTPROD(vec1,vec2)/sqrt(a2*b2);

   if( dp <= -1.0)
      return(180.00);
   else if( dp >= 1.0)
      return(0.0);
   else
   {
      angle = 1.0/DTOR*acos(dp);
      return(angle);
   }
}
/******************************************************************************
 * bond_calc().  Read file and determine bonds and angles within limits       *
 ******************************************************************************/
void
bond_calc(system_mt *system, spec_mt *species, site_mt *site_info, ROOT **broot,
       	  ROOT **aroot, char *spec_mask, int *blim, int *alim, int pbc, int mflag)
{
   BOND		*bond;
   ANGLE        *angle;
   spec_mt	*spec1, *spec2, *spec3;
   register double	dist1, dist2;
   register double	tmp_angle;
   vec_mt	point1, point2, point3;
   vec_mt	vec1, vec2;
   vec_mt	shift, frac;
   vec_mt	shift2, frac2;
   vec_mt	a, min, max;
   register int		i, j, k, u;
   register int		is, js, ks;
   register int		nspec1, nspec2, nspec3;
   register int		nmoli, nmolj, nmolk;
   register int		flag;
   register int	        nsites1, nsites2, nsites3=0;
   NODE		*node;
   mat_mp	h = system->h;
   double       **site1 = (double**)arralloc(sizeof(double),2,
                           0,2,0,system->nsites-1);
   double       **site2 = (double**)arralloc(sizeof(double),2,
                           0,2,0,system->nsites-1);
   double       **site3 = (double**)arralloc(sizeof(double),2,
                           0,2,0,system->nsites-1);

/* Determine no of cells to search through */
   for( u = 0; u < 3; u++)
   {
      a[u] = sqrt(SQR(h[0][u]) + SQR(h[1][u]) + SQR(h[2][u]));
      max[u] = (pbc ? ceil(blim[1]*BOND_SCALE/a[u]): 0.0);
      min[u] = -1.0 * max[u];
   }

/* Scan through selected species and determine distances and angles within limits */
   nmoli = 0;

   for(spec1 = species, nspec1 = 0; spec1 < species+system->nspecies; spec1++, nspec1++)
     if( spec_mask[nspec1] )
     {
	if( !mflag )
	   nsites1 = 1;
	else
        {
           make_sites(system->h, spec1->c_of_m, spec1->quat, spec1->p_f_sites,
                      site1, spec1->nmols, spec1->nsites, pbc?0:1);
           nsites1 = spec1->nsites;
	}
	
        for(i=0; i < spec1->nmols; i++)
        {
           nmoli++;
           nmolj = 0;

           for(spec2 = species, nspec2 = 0; spec2 <= spec1; spec2++, nspec2++)
              if( spec_mask[nspec2] )
              {
	         if( !mflag )
	            nsites2 = 1;
	         else
                 {
                    make_sites(system->h, spec2->c_of_m, spec2->quat, spec2->p_f_sites,
                               site2, spec2->nmols, spec2->nsites, pbc?0:1);
                    nsites2 = spec2->nsites;
              	 }

	         for(j=0; j < (spec2==spec1?i:spec2->nmols); j++)
                 {
                    nmolj++;
                    flag = 0;

                    for(is = 0; is < nsites1; is++)
	            {
		       if( !mflag )
                          for(u = 0; u < 3; u++)
		             point1[u] = spec1->c_of_m[i][u];
		       else
                          for(u = 0; u < 3; u++)
                             point1[u] = site1[u][i*nsites1+is];

                       for( frac[0] = min[0]; frac[0] <= max[0]; frac[0]++) 
                       for( frac[1] = min[1]; frac[1] <= max[1]; frac[1]++) 
                       for( frac[2] = min[2]; frac[2] <= max[2]; frac[2]++) 
                       {
                          mat_vec_mul(h, &frac, &shift, 1);

                          for(js = 0; js < nsites2; js++)
		          {
		             if( !mflag )
                                for(u = 0; u < 3; u++)
		                   point2[u] = spec2->c_of_m[j][u] + shift[u];
		             else
                                for(u = 0; u < 3; u++)
                                   point2[u] = site2[u][j*nsites2+js] + shift[u]; 

                             dist1 = DISTANCE(point2, point1);

                             if( (dist1 >= blim[0]*BOND_SCALE) && (dist1 <= blim[1]*BOND_SCALE) )
                             { 
                                flag = 1;
                                bond = NEW(BOND); /* Create new bond record */
		                if( !mflag && spec2->nsites > 1 )
                                   strncpy((bond->atom1),spec2->name, 2);
				else
                                   strncpy((bond->atom1),site_info[spec2->site_id[js]].name, 3);
		                if( !mflag && spec1->nsites > 1 )
                                   strncpy((bond->atom2),spec1->name, 2);
				else
                                   strncpy((bond->atom2),site_info[spec1->site_id[is]].name, 3);
                                bond->number1 = nmolj;
                                bond->number2 = nmoli;
                                bond->length = dist1;
  	                        node = morethan_BOND(broot,bond);
                                if( node == NULL )
                                {
                                   if( insert_data(broot,bond,1) < 0 )
                                      error("Error creating first node in bond list - \n%s\n",strerror(errno));
                                }
                                else
                                {
                                   if( insert_at_position(broot,node,bond,0) < 0 )
                                      error("Error inserting node in bond list - \n%s\n",strerror(errno));
                                }

                                nmolk = nmolj;

                                /* Calculate angle about spec1 molecule */
                                for(spec3 = spec2, nspec3 = nspec2; spec3 < species+system->nspecies; spec3++, nspec3++)
	                           if( spec_mask[nspec3] )
		                   {
	                              if( !mflag )
	                                 nsites3 = 1;
	                              else
                                      {
                                         make_sites(system->h, spec3->c_of_m, spec3->quat, spec3->p_f_sites,
                                                 site3, spec3->nmols, spec3->nsites, pbc?0:1);
                                         nsites3 = spec3->nsites;
              	                      }

                                      for( k=(spec3==spec2?j+1:0); k < spec3->nmols; k++)
                                      {
                                         nmolk++;

                                         if( (nmolk != nmoli) && (nmolj != nmoli) )
                                         {
                                            for( frac2[0] = min[0]; frac2[0] <= max[0]; frac2[0]++) 
                                            for( frac2[1] = min[1]; frac2[1] <= max[1]; frac2[1]++) 
                                            for( frac2[2] = min[2]; frac2[2] <= max[2]; frac2[2]++) 
                                            {
                                               mat_vec_mul(h, &frac2, &shift2, 1);
					 
                                               for(ks = 0; ks < nsites3; ks++)
		                               {
		                                  if( !mflag )
                                                     for(u = 0; u < 3; u++)
		                                        point3[u] = spec3->c_of_m[k][u] + shift2[u];
		                                  else
                                                     for(u = 0; u < 3; u++)
                                                        point3[u] = site3[u][k*nsites3+ks] + shift2[u]; 

                                                  dist2 = DISTANCE(point3, point1);

                                                  if( (dist2 >= blim[0]*BOND_SCALE ) && (dist2 <= blim[1]*BOND_SCALE) )
                                                  {
                                                     for(u = 0; u < 3; u++)
                                                     {
                                                        vec1[u] = point2[u] - point1[u];  /* Vector between atoms 2 and 1 */
                                                        vec2[u] = point3[u] - point1[u];  /* Vector between atoms 3 and 1 */
                                                     }
                                                     tmp_angle = angle_calc(vec1, vec2);
                       
                                                     if( (tmp_angle >= alim[0]) && (tmp_angle <= alim[1]) )
                                                     {
                                                        angle = NEW(ANGLE); /* Create new angle record */
		                                        if( !mflag && spec1->nsites > 1)
                                                           strncpy((angle->atom1),spec1->name, 2);
						        else
                                                           strncpy((angle->atom1),site_info[spec1->site_id[is]].name, 3);
		                                        if( !mflag && spec2->nsites > 1)
                                                           strncpy((angle->atom2),spec2->name, 2);
						        else
                                                           strncpy((angle->atom2),site_info[spec2->site_id[js]].name, 3);
		                                        if( !mflag && spec3->nsites > 1)
                                                           strncpy((angle->atom3),spec3->name, 2);
						        else
                                                           strncpy((angle->atom3),site_info[spec3->site_id[ks]].name, 3);
                                                        angle->number1 = nmoli;
                                                        angle->number2 = nmolj;
                                                        angle->number3 = nmolk;
                                                        angle->length1 = dist1;
                                                        angle->length2 = dist2;
	                                                angle->value = tmp_angle;
	                                                node = morethan_ANGLE(aroot,angle);
                                                        if( node == NULL )
	                       	                        {
                                                           if( insert_data(aroot,angle,1) < 0 )
                                                              error("Error creating first node in angle list - \n%s\n",strerror(errno));
			 	 	                }
                                                        else
                                                        {
                                                           if( insert_at_position(aroot,node,angle,0) < 0 )
                                                              error("Error inserting node in angle list - \n%s\n",strerror(errno));
                                                        }
                                                     }
					          }
                                               }
                                            }
			                 }
			              }
		                   }
                                   else
                                     nmolk += spec3->nmols;
		             }  
		          }
		       }
	            }
                    if( flag )           /* Calculate angle about spec2 molecule */ 
                    {
                       nmolk = nmoli;
                       for(js = 0; js < nsites2; js++)
		       {
		          if( !mflag )
                             for(u = 0; u < 3; u++)
		                point2[u] = spec2->c_of_m[j][u];
		          else
                             for(u = 0; u < 3; u++)
                                point2[u] = site2[u][j*spec2->nsites+js]; 

                          for(spec3 = spec1, nspec3 = 0; spec3 < species+system->nspecies; spec3++, nspec3++)
		             if( spec_mask[nspec3] )
                             {
                                for( k=(spec3==spec1?i+1:0); k < spec3->nmols; k++)
                                {
                                   nmolk++;

                                   if( (nmolk != nmolj) && (nmoli != nmolj) )
                                   {
                                      for( frac[0] = min[0]; frac[0] <= max[0]; frac[0]++) 
                                      for( frac[1] = min[1]; frac[1] <= max[1]; frac[1]++) 
                                      for( frac[2] = min[2]; frac[2] <= max[2]; frac[2]++) 
                                      {
                                         mat_vec_mul(h, &frac, &shift, 1);
 
                                         for(is = 0; is < nsites1; is++)
	                                 {
		                            if( !mflag )
                                               for(u = 0; u < 3; u++)
		                                  point1[u] = spec1->c_of_m[i][u] + shift[u];
		                            else
                                               for(u = 0; u < 3; u++)
                                                  point1[u] = site1[u][i*nsites1+is] + shift[u];

                                            dist1 = DISTANCE(point1, point2);

                                            if( (dist1 >= blim[0]*BOND_SCALE ) && (dist1 <= blim[1]*BOND_SCALE) )
                                            {
                                               for( frac2[0] = min[0]; frac2[0] <= max[0]; frac2[0]++) 
                                               for( frac2[1] = min[1]; frac2[1] <= max[1]; frac2[1]++) 
                                               for( frac2[2] = min[2]; frac2[2] <= max[2]; frac2[2]++) 
                                               {
                                                  mat_vec_mul(h, &frac2, &shift2, 1);

                                                  for(ks = 0; ks < nsites3; ks++)
		                                  {
		                                     if( !mflag )
                                                        for(u = 0; u < 3; u++)
		                                           point3[u] = spec3->c_of_m[k][u] + shift2[u];
		                                     else
                                                        for(u = 0; u < 3; u++)
                                                           point3[u] = site3[u][k*nsites3+ks] + shift2[u]; 
   
                                                     dist2 = DISTANCE(point3, point2);
   
                                                     if( (dist2 >= blim[0]*BOND_SCALE) && (dist2 <= blim[1]*BOND_SCALE) )
                                                     {
                                                        for(u = 0; u < 3; u++)
                                                        {
                                                           vec1[u] = point1[u] - point2[u];  /* Vector between atoms 1 and 2 */
                                                           vec2[u] = point3[u] - point2[u];  /* Vector between atoms 3 and 2 */
                                                        }
   
                                                        tmp_angle = angle_calc(vec1, vec2);

                                                        if( (tmp_angle >= alim[0]) && (tmp_angle <= alim[1]) )
                                                        {
                                                           angle = NEW(ANGLE);  /* Calculate new angle record */
		                                           if( !mflag && spec2->nsites > 1)
                                                              strncpy((angle->atom1),spec2->name, 2);
						           else
                                                              strncpy((angle->atom1),site_info[spec2->site_id[js]].name, 3);
		                                           if( !mflag && spec1->nsites > 1)
                                                              strncpy((angle->atom2),spec1->name, 2);
						           else
                                                              strncpy((angle->atom2),site_info[spec1->site_id[is]].name, 3);
		                                           if( !mflag && spec3->nsites > 1)
                                                              strncpy((angle->atom3),spec3->name, 2);
						           else
                                                              strncpy((angle->atom3),site_info[spec3->site_id[ks]].name, 3);
                                                           angle->number1 = nmolj;
                                                           angle->number2 = nmoli;
                                                           angle->number3 = nmolk;
                                                           angle->length1 = dist1;
                                                           angle->length2 = dist2;
                                                           angle->value = tmp_angle;
		                                           node = morethan_ANGLE(aroot,angle);
                                                           if( node == NULL )
                                                           {
                                                              if( insert_data(aroot,angle,1) < 0 )
                                                                 error("Error creating first node in angle list - \n%s\n",strerror(errno));
                                                           }
                                                           else
                                                           {
                                                              if( insert_at_position(aroot,node,angle,0) < 0 )
                                                                 error("Error inserting node in angle list - \n%s\n",strerror(errno));
                                                           }
                                                        }
					             }
                                                  }
                                               }
					    }
			                 }
                                      }
                                   }
		                }
                             }
                             else
                               nmolk += spec3->nmols;
                       }
	            }
                 }
	      }
              else
                nmolj += spec2->nmols;
        }
     }
     else
       nmoli += spec1->nmols;
}
/******************************************************************************
 * data_out().  Output bonds and angles to file with same format as Shell by  *
 *              N. Allan and G. Barrera, Bristol Univ.                        *
 ******************************************************************************/
void data_out(ROOT **broot, ROOT **aroot, int xflag)
{
   NODE         *node;
   BOND         *bd;
   ANGLE        *ang;

   puts("Bonds:");
   puts("               i           j                 R_ij");
   puts("               -           -                 ----");

   if(VALID(broot))
   {
      node = (*broot)->head;
      do
      {
          bd = node->data;
	  printf("            %4d -%4s  %4d -%4s   %12.7f\n",
              bd->number1,bd->atom1,bd->number2,bd->atom2,bd->length);
          node = node->next;
      } while(node != NULL);
   }

   puts("\nAngles:");
   puts("   i           j           k         Angle j/i\\k           R_ij          R_ik");
   puts("   -           -           -         -----------           ----          ----");

   if(VALID(aroot))
   {
      node = (*aroot)->head;
      do
      {
          ang = node->data;
	  if( !(xflag == 1 && fmod(floor(1e6*ang->value),9e7) == 0.0))
             printf("%4d -%4s  %4d -%4s  %4d -%4s   %11.6f    %12.7f  %12.7f\n",
                ang->number1,ang->atom1,ang->number2,ang->atom2,ang->number3,
                   ang->atom3,ang->value,ang->length1,ang->length2);
          node = node->next;
      } while(node != NULL);
   }
}
/******************************************************************************
 * main().   Driver program for calculating bond lengths from MOLDY files     *     
 * Acceptable inputs are sys-spec files, restart files or dump files.  	      *
 * Call: mdbond [-s sys-spec-file] [-r restart-file] [-d dump-file].	      *
 * If neither specified on command line, user is interrogated.		      *
 ******************************************************************************/
int
main(int argc, char **argv)
{
   int	c, cflg, ans_i, data_source = 0;
   char 	line[80];
   extern char	*optarg;
   int		errflg = 0;
   int		intyp = 0;
   int		start = 0, finish = 0, inc = 1;
   int		tflag = 0;
   boolean	bflag = false, aflag = false;
   int		irec;
   char         *bondlims = NULL, *anglims = NULL;
   char		*filename = NULL, *dump_base = NULL;
   char		*dump_names = NULL;
   char		*dumplims = NULL, *tempname = NULL;
   char		*dumpcommand;
   int		dump_size;
   float	*dump_buf;
   FILE		*Fp = NULL, *Hp, *Dp;
   restrt_mt	restart_header;
   system_mt	sys;
   spec_mt	*species;
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf = NULL;
   int          av_convert;
   char		*spec_list = NULL;
   char		*spec_mask = NULL;
   int		pbc = 0;    /* No periodic boundary conditions (default) */
   int		mflag = 0;  /* Bond lengths calculated between species cofms, not molecular sites */
   int		xflag = 0;  /* Include right angles (default) */
   int		dump_level = 0;

   int          blim[3], alim[3];         /* Min and max values for bonds and angles */
   ROOT         *root_bond = NULL;        /* Root of bond linked list */
   ROOT         *root_angle = NULL;       /* Root of angle linked list */
   int          arglen, ind, genflg=0;
   int		verbose = 0;

#define MAXTRY 100
   control.page_length=1000000;

   /* Set up default limits */
   sprintf(line, "%d-%d:%d", BOND_MIN, BOND_MAX, BOND_INC);
   bondlims = mystrdup(line);
   sprintf(line, "%d-%d:%d", ANGLE_MIN, ANGLE_MAX, ANGLE_INC);
   anglims = mystrdup(line);

   comm = argv[0];

   while( (c = getopt(argc, argv, "r:s:d:t:g:o:b:a:pxjv?") ) != EOF )
      switch(c)
      {
       case 'r':
       case 's':
	 if( intyp )
	    errflg++;
	 intyp = data_source = c;
	 filename = optarg;
	 break;
       case 'd':
	 dump_base = optarg;
	 break;
       case 't':
         if( tflag++ == 0)
           dumplims = mystrdup(optarg);
         else
           errflg++;
         break;
         break;
       case 'g':
	 spec_list = mystrdup(optarg);
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
            error(NOOUTF, optarg);
	 break;
       case 'b':
         (void)free(bondlims);
         bondlims = mystrdup(optarg);
         bflag = true;
         break;
       case 'a':
         (void)free(anglims);
	 anglims = mystrdup(optarg);
         aflag = true;
         break;
       case 'p': /* Apply periodic boundary conditions */
	 pbc = 1;
         break;
       case 'x': /* Don't print right angles */
	 xflag = 1;
         break;
       case 'j': /* Calculate bonds and angles between c_of_ms */
	 mflag = 1;
         break;
       case 'v':
         verbose++;
         break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fprintf(stderr,"Usage: %s -s sys-spec-file|-r restart-file ",comm);
      fputs("[-d dump-files] [-t s[-f[:n]]] [-g species] ",stderr);
      fputs("[-b bond-limits] [-a angle-limits] [-p] [-x] [-j] [-v] [-o output-file]\n",stderr);
      exit(2);
   }

   if( dump_base )
      data_source = 'd';

   if(intyp == 0)
   {
      fputs("How do you want to specify the simulated system?\n", stderr);
      fputs("Do you want to use a system specification file (1)", stderr);
      fputs(" or a restart file (2)\n", stderr);
      if( (ans_i = get_int("? ", 1, 2)) == EOF )
	 exit(2);
      intyp = ans_i-1 ? 'r': 's';

      if( (filename = get_str("File name? ")) == NULL )
	 exit(2);
   }

   switch(intyp)
   {
    case 's':
      if( (Fp = fopen(filename,"r")) == NULL)
	 error(NOSYSSPEC, filename);
      cflg = check_control(Fp);
      if( cflg )
      {
	 do
	 {
	    fscanf(Fp, "%s",line);
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
	 error(NORESTART, filename, strerror(errno)); 
      re_re_header(Fp, &restart_header, &control);
      re_re_sysdef(Fp, restart_header.vsn, &sys, &species, &site_info, &potpar);
      break;
    default:
      error("Internal error - invalid input type", "");
   }
   allocate_dynamics(&sys, species);

   spec_mask = (char*)calloc(sys.nspecies+1,sizeof(char));

/* Check species selection list */
   if( spec_list == NULL)
     {
     spec_list = malloc((int)log10(sys.nspecies) + 4);
     sprintf(spec_list,"1-%d",sys.nspecies);
     }

   if( tokenise(mystrdup(spec_list), spec_mask, sys.nspecies) == 0 )
      error(INVSPECIES, spec_list, sys.nspecies);

#ifdef DEBUG
   {
      int i;
      fputs("species mask ", stderr);
      for(i = 0; i < sys.nspecies; i++)
         if(spec_mask[i])
            fputs("1", stderr);
         else
            fputs("0", stderr);
      fputs("\n", stderr);
   }
#endif

   /* Input and check bond length limits where necessary */
   while( bondlims)
   {
      if( forstr(bondlims, &(blim[0]), &(blim[1]), &(blim[2])) )
      {
         fputs("Invalid range for bond lengths \"", stderr);
         fputs(bondlims, stderr);
         fputs("\"\n", stderr);
      }

      if( blim[0] >= blim[1] || blim[0] < 0 )
      {
         fputs("Bond length limits must satisfy max >= min and min >= 0\n", stderr);

         (void)free(bondlims);
         bondlims = NULL;
         fputs("Please specify range of bond limits in form min-max\n", stderr);
         bondlims = get_str("min-max? ");
      }
      else
         break;
   }

   if( !bflag )
      message(NULLI,NULLP,INFO,MAXBOND,
              blim[0]*BOND_SCALE,blim[1]*BOND_SCALE,blim[2]*BOND_SCALE);

#if DEBUG
fprintf(stderr, "bond range %d-%d:%d\n", blim[0], blim[1], blim[2]);
#endif

   /* Input and check angle limits where necessary */
   while (anglims)
   {
      if( forstr(anglims, &(alim[0]), &(alim[1]), &(alim[2])) )
      {
         fputs("Invalid range for angles \"", stderr);
         fputs(anglims, stderr);
         fputs("\"\n", stderr);
      }
      if( alim[0] >= alim[1] || alim[0] < 0 )
      {
         fputs("Angle limits must satisfy max >= min and min >= 0\n", stderr);

         (void)free(anglims);
         anglims = NULL;
         fputs("Please specify range of angle limits in form min-max\n", stderr);
         anglims = get_str("min-max? ");
      }
      else
         break;
   }

   if( !aflag )
      message(NULLI, NULLP, INFO, MAXANGLE, alim[0], alim[1], alim[2]);

#if DEBUG
fprintf(stderr, "angle range %d-%d:%d\n", alim[0], alim[1], alim[2]);
#endif


   switch(data_source)                  /* To read configurational data       */
   {
    case 's':                           /* Lattice_start file                 */
      lattice_start(Fp, &sys, species, qpf);
      break;
    case 'r':                           /* Restart file                       */
      init_averages(sys.nspecies, restart_header.vsn,
                    control.roll_interval, control.roll_interval,
                    &av_convert);
      control.rdf_interval = 0; /* Skip rdf data */
      read_restart(Fp, restart_header.vsn, &sys, av_convert);
      break;
    case 'd':
      /* Prepare dump file name for reading */
      if( strstr(dump_base,"%d") )
         genflg++;

      if( genflg == 0 && optind < argc ) {
         arglen = strlen(dump_base);
         for(ind=optind; ind < argc; ind++) {
	    arglen += strlen(argv[ind]) + 1;
         }
         dump_names=malloc(arglen);
         dump_names[0] = 0;
         strcat(dump_names, dump_base);
         if(optind < argc) strcat(dump_names," ");
         for(ind=optind; ind < argc; ind++) {
            strcat(dump_names,argv[ind]);
            if(ind < argc-1) strcat(dump_names," ");
         }
      }
      else
         dump_names = dump_base;

      if( (dumpcommand = malloc(256+strlen(dump_names))) == 0)
         error(COMMEM, 256+strlen(dump_names));

#if defined (HAVE_POPEN)
      sprintf(dumpcommand,"dumpext -c -1 %s%s", verbose?"-v ":"", dump_names);
      if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
      if( (Hp = popen(dumpcommand,"r")) == 0)
         error(DUMPCOMM, dumpcommand);
#else
      tempname = tmpnam((char*)0);
      sprintf(dumpcommand,"dumpext -c -1 -o %s %s%s", tempname, verbose?"-v ":"", dump_names);
      if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
      system(dumpcommand);
      if( (Hp = fopen(tempname,"rb")) == 0)
         error(FILEOPEN,tempname);
#endif
      finish = dump_info(Hp, &dump_level);

      (void)free(dumpcommand);

      if( verbose ) message(NULLI, NULLP, INFO, HEADER);
 
#if defined (HAS_POPEN)
      pclose(Hp);
#else
      fclose(Hp);
      remove(tempname);
#endif

  /* Check dump file contains necessary data */
      if( !(dump_level & 1) )
        error(NOCOMP, "C of M positions", dump_level);

  /*
   *  Ensure that the dump limits start, finish, inc are set up.
   */
      if( tflag )
        if( forstr(dumplims, &start, &finish, &inc) )
           error(INVSLICES, dumplims);

  /*
   * Allocate buffer for data
   */
      dump_size = DUMP_SIZE(~0, sys.nmols, sys.nmols_r)*sizeof(float);

      if( (dump_buf = (float*)malloc(dump_size)) == 0)
         error(BUFFMEM, dump_size);
      if( (dumpcommand = malloc(256+strlen(dump_names))) == 0)
         error(COMMEM, 256+strlen(dump_names));

#if defined (HAVE_POPEN)
      sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d %s%s",
         sys.nmols, sys.nmols_r, start, finish, inc, verbose?"-v ":"", dump_names);
                                                                                
      if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
      if( (Dp = popen(dumpcommand,"r")) == 0)
         error(DUMPCOMM, dumpcommand);
#else
      tempname = tmpnam((char*)0);
      sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d -o %s %s%s",
         sys.nmols, sys.nmols_r, start, finish, inc, tempname, verbose?"-v ":"", dump_names);
      if( verbose ) message(NULLI, NULLP, INFO, EXEC, dumpcommand);
      system(dumpcommand);
      if( (Dp = fopen(tempname,"rb")) == 0)
         error("Failed to open \"%s\"",tempname);
#endif

      for(irec = start; irec <= finish; irec+=inc)
      {
         if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
         {
            if( !strcmp(strerror(errno),"Success") )
	       error(DUMPREC, irec, dump_base, strerror(errno));
            else
               error(DUMPREC0, irec, dump_base, strerror(errno));
         }

         dump_to_moldy(dump_buf, &sys);  /* read dump data */
	 if( !mflag )
            mat_vec_mul(sys.h, sys.c_of_m, sys.c_of_m, sys.nmols);

#ifdef DEBUG
      fprintf(stderr,"Successfully read dump record %d from file \"%s\"\n",
          irec%control.maxdumps, dump_names);
#endif
         /* Perform bond/angle calculations for each slice of dump file */
         bond_calc(&sys, species, site_info, &root_bond, &root_angle, spec_mask, blim, alim, pbc, mflag); 
         printf("- Time slice %d -\n",irec);
         data_out(&root_bond, &root_angle, xflag);
         putchar('\n');
         if( root_bond != NULL && delete_list(&root_bond))
            error("Error releasing bond list data for slice %d - \n%s\n",
               irec, strerror(errno));
         if( root_angle != NULL && delete_list(&root_angle))
            error("Error releasing angle list data for slice %d - \n%s\n",
               irec, strerror(errno));
      } 
#if defined (HAVE_POPEN)
      pclose(Dp);
#else
      fclose(Dp);
      remove(tempname);
#endif
      break;
    default:
      break;
   }
   if( data_source == 's' || data_source == 'r' )
   {
/* Convert molecule positions from frac coords to Cartesian coords */
      if( !mflag )
         mat_vec_mul(sys.h, sys.c_of_m, sys.c_of_m, sys.nmols);

      bond_calc(&sys, species, site_info, &root_bond, &root_angle, spec_mask, blim, alim, pbc, mflag);
      data_out(&root_bond, &root_angle, xflag);
   }
   if( verbose ) message(NULLI, NULLP, INFO, COMPLETE);
   return 0;
}
