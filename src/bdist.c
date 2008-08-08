/* BDIST
Copyright (C) 1999, 2003 Craig Fisher */
 
/**************************************************************************************
 * bdist   code for calculating bond distributions from Shell-style bond output files *
 ************************************************************************************** 
 *  Revision Log
 *  $Log: bdist.c,v $
 *  Revision 1.10  2005/02/04 14:49:41  cf
 *  Common utility messages/errors moved to utlsup.h.
 *
 *  Revision 1.9  2004/12/07 13:00:02  cf
 *  Merged with latest utilities.
 *
 *  Revision 1.7.10.3  2004/12/06 19:03:40  cf
 *  Added stdlib.h to headers.
 *  Set uninitialized node to NULL.
 *  Removed unused variables.
 *
 *  Revision 1.7.10.2  2003/11/07 09:01:08  moldydv
 *  Fixed bug in bond length units.
 *
 *  Revision 1.7.10.1  2003/07/29 09:32:04  moldydv
 *  Able to handle multiple time slice output from mdbond.
 *  Added option -p to output coordination spheres of each molecule.
 *
 *  Revision 1.7  2002/09/19 09:26:27  kr
 *  Tidied up header declarations.
 *  Changed old includes of string,stdlib,stddef and time to <> form
 *
 *  Revision 1.6  2000/04/27 17:57:06  keith
 *  Converted to use full ANSI function prototypes
 *
 *  Revision 1.5  1999/11/01 17:16:43  keith
 *  Fixed lint complaints.
 *
 *  Revision 1.4  1999/09/24 10:51:54  keith
 *  Minor changes to Usage message.
 *
 *  Revision 1.4  1999/09/24 16:44:30  craig
 *  Minor changes to Usage message.
 *
 *  Revision 1.3  1999/09/21 15:34:10  keith
 *  Minor portability modifications.
 *
 *  Revision 1.2  1999/09/14 13:30:35  keith
 *  Fixed "return of ptr to stack var" error.
 *
 *  Revision 1.1  1999/07/22 14:02:26  keith
 *  Initial revision
 *
 *  Revision 1.0  1999/06/24 11:15:15  craig 
 *  Initial revision
 *
 */
#include "defs.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "structs.h"
#include "messages.h"
#include "list.h"
#include "utlsup.h"

#define LLENGTH 80
#define DISTRIB 0
#define COORD 1
/*
 * Structures for bond, angle and count data. Each bond and angle record has a linked list
 * of count data.
 */
typedef struct
{
   char	atom1[3];
   char atom2[3];
   int  number1;
   int  number2;
   double length;
} BOND;

typedef struct
{
   char atom1[3];
   char atom2[3];
   char atom3[3];
   int  number1;
   int  number2;
   int  number3;
   double length1;
   double length2;
   double value;
} ANGLE;

typedef struct
{
   char	atom1[3];
   char atom2[3];
   ROOT	*root_count;
} BOND_COUNTER;

typedef struct
{
   char atom1[3];
   char atom2[3];
   char atom3[3];
   ROOT *root_count;
} ANG_COUNTER;

typedef struct
{
   int  xaxis;
   int  yaxis; 
} COUNT;

typedef struct
{
   char spec[3];
   int  id;
   ROOT *root_nabor;
} POLY;

typedef struct
{
   char spec[3];
   int  count;
} NABOR;

/****************************************************************************
 * search_list_atom(). Scan list for specific atom                          *
 ****************************************************************************/
NODE *search_list_atom(ROOT **root, int atomno, int flag)
{
   NODE    *node = NULL;
   POLY    *data;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         data = node->data;
	 if( flag ) { /* Next greater atom no */
            if( data->id > atomno )
	       break;
	 } else { /* Exact match with atom no */
            if( data->id == atomno )
	       break;
	 }
         node = node->next;
      }
   }

   return node;
}
/****************************************************************************
 * search_list_spec(). Scan list for species match                          *
 ****************************************************************************/
NODE *search_list_spec(ROOT **root, char *atom)
{
   NODE    *node = NULL;
   NABOR   *data;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         data = node->data;
         if( !strcmp(atom, data->spec) )
            break;
         node = node->next;
      }
   }
   return node;
}
/****************************************************************************
 * max_bond(). return maximum distance (in angstroms) from bond list        *
 ****************************************************************************/
double max_bond(ROOT **root)
{
   NODE    *node = NULL;
   BOND    *bd;
   double  maximum = 0;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         bd = node->data;
         maximum = MAX(bd->length, maximum);
         node = node->next;
      }
   }
   return maximum;
}
/****************************************************************************
 * max_angle(). return maximum angle (in degrees) from angle list           *
 ****************************************************************************/
double max_angle(ROOT **root)
{
   NODE    *node = NULL;
   ANGLE   *ang;
   double  maximum = 0;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         ang = node->data;
         maximum = MAX(ang->value, maximum);
         node = node->next;
      }
   }
   return maximum;
}
/****************************************************************************
 * search_bond(). compare atoms in bond list with those read from file       *
 ****************************************************************************/
NODE *search_bond(ROOT **root, char *atom1, char *atom2)
{
   NODE    *node = NULL;
   BOND    *bd;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         bd = node->data;
         if( (!strcmp(atom1, bd->atom1) && !strcmp(atom2, bd->atom2)) ||
             (!strcmp(atom1, bd->atom2) && !strcmp(atom2, bd->atom1)))
             break;
         node = node->next;
      }
   }
   return node;
}
/****************************************************************************
 * search_angle(). compare atoms in angle list with those read from file    *
 ****************************************************************************/
NODE *search_angle(ROOT **root, char *atom1, char *atom2, char *atom3)
{
   NODE    *node = NULL;
   ANGLE   *ang;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         ang = node->data;
         if( !strcmp(atom1, ang->atom1) )
            if( (!strcmp(atom2, ang->atom2) && !strcmp(atom3, ang->atom3)) ||
                (!strcmp(atom2, ang->atom3) && !strcmp(atom3, ang->atom2)))
                   break;
         node = node->next;
      }
   }
   return node;
}
/****************************************************************************
 * search_count(). compare counts in linked list with those read from file  *
 ****************************************************************************/
NODE *search_count(ROOT **root, int x)
{
   NODE    *node = NULL;
   COUNT   *ct;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         ct = node->data;
         if( ct->xaxis == x )
             break;
         node = node->next;
      }
   }
   return node;
}
/****************************************************************************
 * count_poly(). Calculate no of neighbours of each molecule                *
 ****************************************************************************/
void count_poly(ROOT **bond_root, ROOT **poly_root, int *blim)
{
   NODE         *node_bond;
   NODE         *node_poly, *next_node;
   NODE         *node_nabor;
   BOND         *bd;
   POLY         *poly;
   NABOR        *nabor;

   if(VALID(bond_root))
   {
      node_bond = (*bond_root)->head;
      do
      {
         bd = node_bond->data;
         if(VALID(&poly_root) && bd->length >= blim[0]*BOND_SCALE && bd->length <= blim[1]*BOND_SCALE)
         {
            node_poly = search_list_atom(poly_root, bd->number1, 0);
            if( node_poly == NULL)
            {
               poly = NEW(POLY);
	       strncpy(poly->spec, bd->atom1, 3);
               poly->id = bd->number1;
	       next_node = search_list_atom(poly_root, poly->id, 1);
	       if( next_node != NULL)
                  insert_at_position(poly_root, next_node, poly, 0);
	       else
		  insert_data(poly_root, poly, 1);
            }
            else
	       poly = node_poly->data;
	       
            node_nabor = search_list_spec(&(poly->root_nabor), bd->atom2);
            if( node_nabor == NULL)
            {
               nabor = NEW(NABOR);
	       strncpy(nabor->spec, bd->atom2, 3);
               nabor->count = 1;
               insert_data(&(poly->root_nabor), nabor, 1);
            }
            else
            {
               nabor = node_nabor->data;
               nabor->count++;
            }
           
            node_poly = search_list_atom(poly_root, bd->number2, 0);
            if( node_poly == NULL ) 
            {
               poly = NEW(POLY);
	       strncpy(poly->spec, bd->atom2, 3);
	       poly->id = bd->number2;
	       next_node = search_list_atom(poly_root, poly->id, 1);
	       if( next_node != NULL)
                  insert_at_position(poly_root, next_node, poly, 0);
	       else
                  insert_data(poly_root, poly, 1);
	    }
	    else
	       poly = node_poly->data;

            node_nabor = search_list_spec(&(poly->root_nabor), bd->atom1);
            if( node_nabor == NULL)
            {
               nabor = NEW(NABOR);
	       strncpy(nabor->spec, bd->atom1, 3);
               nabor->count = 1;
               insert_data(&(poly->root_nabor), nabor, 1);
            }
            else
            {
               nabor = node_nabor->data;
               nabor->count++;
            }
	 }
         node_bond = node_bond->next;
      } while(node_bond != NULL);
   }
}
/****************************************************************************
 * bond_distribution(). create bond distribution list 		            *
 ****************************************************************************/
void bond_distribution(ROOT **bond_root, ROOT **bct_root, int *blim)
{
   NODE    *bd_node, *bct_node, *node_count;
   BOND    *bd;
   BOND_COUNTER	*bct;
   COUNT   *count;
   int	   u;

   /* Count no of bonds within each interval */
   if(VALID(bond_root))
   {
      bd_node = (*bond_root)->head;
      do
      {
         bd = bd_node->data;
         if(( bd->length >= blim[0]*BOND_SCALE ) && (bd->length <= blim[1]*BOND_SCALE))
         {
            bct_node = search_bond(bct_root, bd->atom1, bd->atom2);
            if( bct_node == NULL )
            {
               bct = NEW(BOND_COUNTER);
               strncpy((bct->atom1), bd->atom1, 3);
               strncpy((bct->atom2), bd->atom2, 3);
               bct->root_count = NULL;
               insert_data(bct_root, bct, 1);
            }
            else
               bct = bct_node->data; 

            u = (int)floor((bd->length/BOND_SCALE - blim[0])/blim[2]+0.5);
	    node_count = search_count(&(bct->root_count),u);
            if( node_count == NULL )
            {
               count = NEW(COUNT);
               count->xaxis = u;
               count->yaxis = 1;
               insert_data(&(bct->root_count), count, 1);
            }
	    else
            {
	       count = node_count->data;
	       count->yaxis++;
	    }
         }
         bd_node = bd_node->next;
      } while(bd_node != NULL);
   }
}
/****************************************************************************
 * angle_distribution(). create angle distribution list 		    *
 ****************************************************************************/
void angle_distribution(ROOT **angle_root, ROOT **act_root, int *alim)
{
   NODE    *ang_node, *act_node, *node_count;
   ANGLE   *ang;
   ANG_COUNTER  *act;
   COUNT   *count;
   int	   u;

   /* Count no of angles within each interval */
   if(VALID(angle_root))
   {
      ang_node = (*angle_root)->head;
      do
      {
         ang = ang_node->data;
         if(( ang->value >= alim[0]) && ( ang->value <= alim[1]))
         { 
            act_node = search_angle(act_root, ang->atom1, ang->atom2, ang->atom3);
            if( act_node == NULL )
            {
               act = NEW(ANG_COUNTER);
               strncpy((act->atom1), ang->atom1, 3);
               strncpy((act->atom2), ang->atom2, 3);
               strncpy((act->atom3), ang->atom3, 3);
               act->root_count = NULL;
               insert_data(act_root, act, 1);
            } 
            else
               act = act_node->data;

            u = (int)floor((ang->value - alim[0])/alim[2] - 0.5);
            node_count = search_count(&(act->root_count),u);
            if( node_count == NULL )
            {
               count = NEW(COUNT);
               count->xaxis = u;
               count->yaxis = 1;
               insert_data(&(act->root_count), count, 1);
            }
	    else
            {
	       count = node_count->data;
	       count->yaxis++;
	    }
         }
         ang_node = ang_node->next;
      } while(ang_node != NULL);
   }
}
/****************************************************************************
 * write_BOND(). print out bond distribution list                           *
 ****************************************************************************/
void write_BOND(ROOT **root, int *blim, boolean zero)
{
   NODE         *node;
   NODE         *node_count;
   BOND_COUNTER	*bd;
   COUNT	*ct;
   ROOT		*ct_root;

   int		u, num_intervals;
   double	xvalue;

   printf("Bond Distributions (per %4.2f Angstroms):\n",blim[2]*BOND_SCALE);

   /* Calculate total no of bond intervals */
   num_intervals = (int)ceil(1.0*(blim[1] - blim[0]) / blim[2]);

   if(VALID(root))
   {
      node = (*root)->head;
      do
      {
          bd = node->data;
          printf("\n%s %s\n", bd->atom1, bd->atom2);
          ct_root = bd->root_count;
          if(VALID(&ct_root))
          {
            node_count = ct_root->head;
            u = 0;
            do
            {
               ct = node_count->data;
               if (zero) /* Include bins with zero counts */
                 while (u<ct->xaxis)
                 {
                   xvalue = (blim[0] + u*blim[2])*BOND_SCALE;
                   printf("%5.2f     0\n", xvalue);
                   u++;
                 }
               xvalue = (blim[0] + (ct->xaxis)*blim[2])*BOND_SCALE;
               u++;
               printf("%5.2f     %d\n", xvalue, ct->yaxis);
               node_count = node_count->next;
            } while(node_count != NULL);
          }
          if (zero) /* Add any trailing bins with zero counts */
            while (u < num_intervals)
            {
              xvalue = (blim[0] + u*blim[2])*BOND_SCALE;
              printf("%5.2f     0\n", xvalue);
              u++;
            }
          node = node->next;
      } while(node != NULL);
   }
}
/****************************************************************************
 * write_ANGLE(). print out angle distribution list                         *
 ****************************************************************************/
void write_ANGLE(ROOT **root, int *alim, boolean zero)
{
   NODE         *node;
   NODE		*node_count;
   ANG_COUNTER	*ang;
   COUNT	*ct;
   ROOT		*ct_root;

   int		u, num_intervals;
   double	xvalue;

   printf("\nAngle Distributions (per %d degree%s):\n",alim[2],alim[2]==1?"":"s");

   /* Calculate total no of angle intervals */
   num_intervals = (int)ceil(1.0*(alim[1] - alim[0]) / alim[2]);

   if(VALID(root))
   {
      node = (*root)->head;
      do
      {
          ang = node->data;
          printf("\n%s %s %s\n", ang->atom1, ang->atom2, ang->atom3);
          ct_root = ang->root_count;
          if(VALID(&ct_root))
          {
            node_count = ct_root->head;
            u = 0;
            do
            {
               ct = node_count->data;
               if (zero) /* Include bins with zero counts */
                 while (u<ct->xaxis)
                 {
                   xvalue = (alim[0] + u*alim[2]);
                   printf("%5.1f     0\n", xvalue);
                   u++;
                 }
               xvalue = (alim[0] + (ct->xaxis)*alim[2]);
               printf("%5.1f     %d\n", xvalue, ct->yaxis);
               u++;
               node_count = node_count->next;
            } while(node_count != NULL);
          }
          if (zero) /* Add any trailing bins with zero counts */
            while (u < num_intervals)
            {
              xvalue = (alim[0] + u*alim[2]);
              printf("%5.1f     0\n", xvalue);
              u++;
            }
          node = node->next;
      } while(node != NULL);
      putchar('\n');
   }
}
/****************************************************************************
 * write_poly(). List no of neighbours of each molecule                     *
 ****************************************************************************/
void write_poly(ROOT **poly_root)
{
   NODE         *node_poly;
   NODE         *node_nabor;
   POLY         *poly;
   NABOR        *nabor;
   ROOT         *nb_root;

   if(VALID(poly_root))
   {
      node_poly = (*poly_root)->head;
      do
      {
         poly = node_poly->data;
         nb_root = poly->root_nabor;
         printf("%4d %2s",poly->id, poly->spec);
         if(VALID(&nb_root))
         {
            putchar(' ');
            node_nabor = nb_root->head;
            do
            {
               nabor = node_nabor->data;
               printf("    %2s: %-5d", nabor->spec, nabor->count);
               node_nabor = node_nabor->next;
            } while(node_nabor != NULL);
         }
         putchar('\n');
         node_poly = node_poly->next;
      } while(node_poly != NULL);
      putchar('\n');
   }
}
/******************************************************************************
 * main().   Driver program for accumulating and processing bond data         *
 ******************************************************************************/
int main(int argc, char **argv)
{
   extern char	*optarg;
   int		u, errflg = 0;
   char		*filename = NULL;
   char		dummy[80];
   char		*bondlims = NULL, *anglims = NULL;
   char		atom1[3], atom2[3], atom3[3];   /* Atom labels */
   double	alpha, rij, r1, r2;             /* Angle and bond length variables */
   int		blim[3], alim[3]; 
   boolean	bflag = false, aflag = false;   /* Flags for input of limits */
   int		a, b, c;
   int		outsw = 0;
   boolean      inflag = false, zero_flag = false;
   int		i, lineno=0;
   int		nslice=-1;
   char		keyword[LLENGTH], line[LLENGTH];

   BOND		*bond;
   ANGLE	*angle;
   ROOT         *root_bond = NULL;          /* Root of linked list */
   ROOT		*root_angle = NULL;         /* Root of linked list */
   ROOT         *root_bct = NULL;           /* Root of linked list */
   ROOT		*root_act = NULL;           /* Root of linked list */
   ROOT         *root_poly = NULL;          /* Root of linked list */

   FILE		*Fp;

   /* Set up default limits */
   sprintf(dummy, "%d-%d:%d", BOND_MIN, BOND_MAX, BOND_INC);
   bondlims = mystrdup(dummy);
   sprintf(dummy, "%d-%d:%d", ANGLE_MIN, ANGLE_MAX, ANGLE_INC);
   anglims = mystrdup(dummy);

   comm = argv[0];
   if( strstr(comm, "bdist") )
     outsw = DISTRIB;
   else
     outsw = COORD;

   while( (u = getopt(argc, argv, "o:a:b:pz?") ) != EOF )
      switch(u)
      {
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
            error(NOOUTF, optarg);
	 break;
       case 'a':
         (void)free(anglims);
	 anglims = mystrdup(optarg);
         aflag = true;
	 break;
       case 'b':
         (void)free(bondlims);
	 bondlims = mystrdup(optarg);
         bflag = true;
	 break;
       case 'p':
	 outsw = COORD;
	 break;
       case 'z':
	 zero_flag = true;  /* Include zero values in output */
	 break;
       default:
       case '?': 
	 errflg++;
      }

   if( errflg )
   {
      fprintf(stderr, "Usage: %s [-b bond-limits] ",comm);
      fprintf(stderr,"[-a angle-limits] %s[-o output-file] input-file\n",outsw==COORD?"":"[-p] ");
      exit(2);
   }

   if( optind <= argc)
   {
      filename = argv[optind];
      inflag = true;
   }

   if( !inflag)
     filename = get_str("Input filename? ");

   if( (Fp = fopen(filename,"r")) == NULL)
   {
      perror("Error: Couldn't open file for reading");
      exit(2);
   }

   /* Input and check bond length limits where necessary */
   while (bondlims)
   {
      if( forstr(bondlims, &(blim[0]), &(blim[1]), &(blim[2])) )
      {
         fputs("Invalid range for bond lengths \"", stderr);
         fputs(bondlims, stderr);
         fputs("\"\n", stderr);
      }

      if( blim[0] > blim[1] || blim[0] < 0 || blim[2] <= 0 )
      {
         fputs("Bond length limits must satisfy", stderr);
         fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
         (void)free(bondlims);
         bondlims = NULL;
         fputs("Please specify range of bond limits in form", stderr);
         fputs(" start-finish:increment\n", stderr);
         bondlims = get_str("s-f:n? ");
      }
      else
        break;
   }

   /* Input and check angle limits where necessary */
   while (anglims)
   {
      if( forstr(anglims, &(alim[0]), &(alim[1]), &(alim[2])) )
      {
         fputs("Invalid range for angles \"", stderr);
         fputs(anglims, stderr);
         fputs("\"\n", stderr);
      }

      if( alim[0] > alim[1] || alim[0] < 0 || alim[2] <= 0 )
      {
         fputs("Angle limits must satisfy", stderr);
         fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
         (void)free(anglims);
         anglims = NULL;
         fputs("Please specify range of angle limits in form", stderr);
         fputs(" start-finish:increment\n", stderr);
         anglims = get_str("s-f:n? ");
      }
      else
         break;
   }

   while (!feof(Fp))
   {
      /* Start reading bond list file */
      fgets(dummy,sizeof(dummy),Fp);
      lineno++;

      sscanf(dummy,"- Time slice %d -",&nslice);

#ifdef DEBUG
      fprintf(stderr,"Time slice %d\n",nslice);
#endif
      if( !strncmp(dummy,"Bonds:",6) )
      {
         for( i = 0; i < 2; i++, lineno++)
            fgets(dummy,sizeof(dummy),Fp);

         /* Read bond length data from input file */
         while ( !feof(Fp))
         { 
            sscanf(get_line(line,LLENGTH,Fp,1),"%s",keyword);
            lineno++;

            if( !strncmp(keyword,"Angles:",7) )
            {
	       lineno++;
               break;
	    }
            else
               if( sscanf(line,"%d - %s %d - %s %lf",&a,atom1,&b,atom2,&rij) < 5)
                  error("Error in line %d of \"%s\" -- should have 5 parameters", lineno, filename);

            bond = NEW(BOND);
            bond->number1 = a;
            bond->number2 = b;
            strncpy((bond->atom1), atom1, 3);
            strncpy((bond->atom2), atom2, 3);
            bond->length = rij;
            insert_data(&root_bond, bond, 1);
         }

         /* Read angle data from input file */
         for( i = 0; i < 2; i++, lineno++)
            fgets(dummy,sizeof(dummy),Fp);
      
         while ( !feof(Fp) )
         {
            if( strcmp(get_line(line,LLENGTH,Fp,0),"") )
            {
               if( sscanf(line,"%d - %s %d - %s %d - %s %lf %lf %lf\n",&a,atom1,&b,atom2,&c,atom3,&alpha,&r1,&r2) < 9)
                  error("Error in line %d of \"%s\" -- should have 9 parameters", lineno, filename);
               else
                  lineno++;

               angle = NEW(ANGLE);
               angle->number1 = a;
               angle->number2 = b;
               angle->number3 = c;
               strncpy((angle->atom1), atom1, 3);
               strncpy((angle->atom2), atom2, 3);
               strncpy((angle->atom3), atom3, 3);
               angle->length1 = r1;
               angle->length2 = r2;
               angle->value = alpha;
               insert_data(&root_angle, angle, 1);
            }
	    else
	    {
	       lineno++;
	       break;
	    }
         }
#ifdef DEBUG
	 fprintf(stderr,"Line no %d\n",lineno);
#endif
	 if( nslice > -1 )
            printf("- Time slice %d -\n\n",nslice);

/* Set maximum bond distance if not specified by user */
         if( !bflag )
            {
            blim[1] = (int)ceil(max_bond(&root_bond)/BOND_SCALE); /* Default is maximum distance in list */
            message(NULLI,NULLP,INFO,MAXBOND,
                    blim[0]*BOND_SCALE,blim[1]*BOND_SCALE,blim[2]*BOND_SCALE);
            }


/* Set maximum angle if not specified by user */
         if( !aflag )
            {
            alim[1] = (int)ceil(max_angle(&root_angle)); /* Default is maximum angle in list */
            message(NULLI, NULLP, INFO, MAXANGLE, alim[0], alim[1], alim[2]);
            }
#if DEBUG
fprintf(stderr, "bond range  %d-%d:%d\n", blim[0], blim[1], blim[2]);
fprintf(stderr, "angle range %d-%d:%d\n", alim[0], alim[1], alim[2]);
#endif

         if( outsw == DISTRIB)
         {
            /* Calculate bond and angle distributions */
            bond_distribution(&root_bond, &root_bct, blim);
            angle_distribution(&root_angle, &root_act, alim);
            /* Output data in column format */
            write_BOND(&root_bct, blim, zero_flag);
            write_ANGLE(&root_act, alim, zero_flag);
            delete_list(&root_bct);
            delete_list(&root_act);
         }
         else
         {
            /* Calculate and output no of bonds for each molecule */
            count_poly(&root_bond, &root_poly, blim);
            write_poly(&root_poly);
            delete_list(&root_poly);
         }
         delete_list(&root_bond);
         delete_list(&root_angle);
      }
   }
   fclose(Fp);
   return 0;
}
