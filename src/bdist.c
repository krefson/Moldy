/* BDIST
Copyright (C) 1999 Craig Fisher */
 
/**************************************************************************************
 * bdist   code for calculating bond distributions from Shell-style bond output files *
 ************************************************************************************** 
 *  Revision Log
 *  $Log: bdist.c,v $
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
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include "stdlib.h"
#include "string.h"
#include "list.h"

/*
 * Default limits for bond intervals and angle intervals - integers only
 */
#define BOND_MIN  0
#define BOND_MAX  20          /* Bond intervals in tenths of Angstroms */
#define BOND_INC  5

#define ANGLE_MIN  0
#define ANGLE_MAX  180        /* Angle intervals in degrees */
#define ANGLE_INC  1

/*
 * Structures for bond, angle and count data. Each bond and angle record has a linked list
 * of count data.
 */
typedef struct
{
   char	atom1[3];
   char atom2[3];
   ROOT	*root_count;
} BOND;

typedef struct
{
   char atom1[3];
   char atom2[3];
   char atom3[3];
   ROOT *root_count;
} ANGLE;

typedef struct
{
   int  xaxis;
   int  yaxis; 
} COUNT;

/*
 * strcpy with call to memory allocation
 */
static char * mystrdup(s)
char *s;
{
   char * t = NULL;
   if(s) t=malloc(strlen(s)+1);
   return t?strcpy(t,s):0;
}
/******************************************************************************
 * forstr.  Parse string str of format s-f:n  (final 2 parts optional),       *
 *          returning integer values of s,f,n.                                *
 ******************************************************************************/
int forstr(instr, start, finish, inc)
char    *instr;
int     *start, *finish, *inc;
{
   char *p, *pp, *str = mystrdup(instr);
   long strtol();
  
   if( (p = strchr(str,':')) != NULL)
   {
      *inc = strtol(p+1, &pp, 0);
      if( pp == p+1 )
         goto limerr;
      *p = 0;
   }
   if( (p = strchr(str,'-')) != NULL)
   {
      *p = 0;
      *start = strtol(str, &pp, 0);
      if( pp == str )
         goto limerr;
      *finish = strtol(p+1, &pp, 0);
      if( pp == p+1 )
         goto limerr;
   }
   else
   {
      *finish = strtol(str, &pp, 0);
      if( pp == str )
         goto limerr;
   }
   return 0;
 limerr:
   return -1;
}
/******************************************************************************
 * get_str().  Read a string from stdin, issuing a prompt.                    *
 ******************************************************************************/
char    *get_str(prompt)
char    *prompt;
{
   char         ans_str[80];
   char         *str = malloc(80);
   int          ans_flag;

   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, "%s", str) == 1)
         ans_flag++;
   }
   if( ans_flag )
      return(str);
   else
      return(NULL);
}
/****************************************************************************
 * check_bond(). compare atoms in bond list with those read from file       *
 ****************************************************************************/
int check_bond(root, atom1, atom2)
ROOT	**root;
char	*atom1;
char	*atom2;
{
   NODE    *node;
   BOND    *bd;
   int	   i=-1;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         i++;
         bd = node->data;
         if( (!strcmp(atom1, bd->atom1) && !strcmp(atom2, bd->atom2)) ||
             (!strcmp(atom1, bd->atom2) && !strcmp(atom2, bd->atom1)))
             break;
         node = node->next;
      }
      if( node == NULL)
         i = -1;
   }
   return i;
}
/****************************************************************************
 * check_angle(). compare atoms in angle list with those read from file     *
 ****************************************************************************/
int check_angle(root, atom1, atom2, atom3)
ROOT	**root;
char	*atom1;
char	*atom2;
char	*atom3;
{
   NODE    *node;
   ANGLE   *ang;
   int	   i=-1;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         i++;
         ang = node->data;
         if( !strcmp(atom1, ang->atom1) )
            if( (!strcmp(atom2, ang->atom2) && !strcmp(atom3, ang->atom3)) ||
                (!strcmp(atom2, ang->atom3) && !strcmp(atom3, ang->atom2)))
                   break;
         node = node->next;
      }
      if( node == NULL)
         i = -1;
   }
   return i;
}
/****************************************************************************
 * check_count(). compare counts in linked list with those read from file   *
 ****************************************************************************/
int check_count(root, x)
ROOT	**root;
int	x;
{
   NODE    *node;
   COUNT   *ct;
   int	   i = -1;

   if( VALID(root) )
   {
      node = (*root)->head;
      while(node != NULL)
      {
         i++;
         ct = node->data;
         if( ct->xaxis == x )
         {
             ct->yaxis = ct->yaxis + 1;
             break;
         }
         node = node->next;
      }
      if( node == NULL )
         i = -1;
   }
   return i;
}
/****************************************************************************
 * display_BOND(). print out bond distribution list                         *
 ****************************************************************************/
void display_BOND(root, bstart, binc)
ROOT	**root;
int	bstart;
int     binc;
{
   NODE         *node;
   NODE         *node_count;
   BOND		*bd;
   COUNT	*ct;
   ROOT		*ct_root;

   double	xvalue;

   printf("Bond Distributions (per %4.2f Ang.):\n",binc/10.0);

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
            do
            {
               ct = node_count->data;
               xvalue = (bstart + (ct->xaxis+0.5)*binc)/10.0;
               printf("%5.2f     %d\n", xvalue, ct->yaxis);
               node_count = node_count->next;
            } while(node_count != NULL);
          }
          node = node->next;
      } while(node != NULL);
   }
}
/****************************************************************************
 * display_ANGLE(). print out angle distribution list                       *
 ****************************************************************************/
void display_ANGLE(root, astart, ainc)
ROOT	**root;
int	astart;
int	ainc;
{
   NODE         *node;
   NODE		*node_count;
   ANGLE	*ang;
   COUNT	*ct;
   ROOT		*ct_root;

   double	xvalue;

   printf("\n\nAngle Distributions (per %d deg.):\n",ainc);

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
            do
            {
               ct = node_count->data;
               xvalue = (astart + (ct->xaxis+0.5)*ainc);
               printf("%5.1f     %d\n", xvalue, ct->yaxis);
               node_count = node_count->next;
            } while(node_count != NULL);
          }
          node = node->next;
      } while(node != NULL);
   }
}
/******************************************************************************
 * main().   Driver program for accumulating distribution data                *
 ******************************************************************************/
int main(argc, argv)
int	argc;
char	*argv[];
{
   extern char	*optarg;
   int		errflg = 0;
   char		*filename = NULL;
   char		dummy[80];
   char		*bondlims = NULL, *anglims = NULL;
   char		atom1[3], atom2[3], atom3[3];   /* Atom labels */
   double	alpha, rij, r1, r2;  /* Angle and bond length variables */
   int		astart, afin, ainc, bstart, bfin, binc; 
   int		bflag=0, aflag=0;            /* Flags for input of limits */
   char		a[7];
   int		b, c;
   int		inflag = 0;
   int		u;                  /* Index no of current bond/angle */
   int		nbond, nangle;      /* Total no of bond and angle records */
   int		i, check;

   BOND		*bond;
   ANGLE	*angle;
   COUNT	*count;
   ROOT		*root_bond = NULL;          /* Root of linked list */
   ROOT		*root_angle = NULL;         /* Root of linked list */
   NODE		*node;

   FILE		*Fp;

   while( (u = getopt(argc, argv, "i:o:a:b:") ) != EOF )
      switch(u)
      {
       case 'i':
	 filename = optarg;
         inflag++;
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    perror("failed to open file \"%s\" for output");
	 break;
       case 'a':
	 anglims = mystrdup(optarg);
	 break;
       case 'b':
	 bondlims = mystrdup(optarg);
	 break;
       default:
       case '?': 
	 errflg++;
      }

   if( errflg )
   {
      fputs("Usage: bdist [-i input-file] [-b bond-limits] ",stderr);
      fputs("[-a angle-limits] [-o output-file]\n",stderr);
      exit(2);
   }

   if( inflag == 0)
     filename = get_str("Input filename? ");

   if( (Fp = fopen(filename,"r")) == NULL)
   {
      perror("Error: Couldn't open file for reading");
      exit(2);
   }

   /* Set default values for bond limits (x10) */
   bstart = BOND_MIN;
   bfin = BOND_MAX;
   binc = BOND_INC;

   if( bondlims == NULL )
       bflag++;

   /* Input and check bond length limits where necessary */
   while (!bflag)
   {
      if( forstr(bondlims, &bstart, &bfin, &binc) )
      {
         fputs("Invalid range for bond lengths \"", stderr);
         fputs(bondlims, stderr);
         fputs("\"\n", stderr);
      }
      else
         bflag++;
      if( bstart > bfin || bstart < 0 || binc <= 0 )
      {
         fputs("Bond length limits must satisfy", stderr);
         fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
         bflag = 0;
      }
      if( !bflag)
      {
         bstart = BOND_MIN;
         bfin = BOND_MAX;
         binc = BOND_INC;
         (void)free(bondlims);
         bondlims = NULL;
         fputs("Please specify range of bond limits in form", stderr);
         fputs(" start-finish:increment\n", stderr);
         bondlims = get_str("s-f:n? ");
      }
   }

   /* Set default values for angle limits */
   astart = ANGLE_MIN;
   afin = ANGLE_MAX;
   ainc = ANGLE_INC;

   if( anglims == NULL )
       aflag++;

   /* Input and check angle limits where necessary */
   while (!aflag)
   {
      if( forstr(anglims, &astart, &afin, &ainc) )
      {
         fputs("Invalid range for angles \"", stderr);
         fputs(anglims, stderr);
         fputs("\"\n", stderr);
      }
      else
         aflag++;
      if( astart > afin || astart < 0 || ainc <= 0 )
      {
         fputs("Angle limits must satisfy", stderr);
         fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
         aflag=0;
      }
      if( !aflag)
      {
         astart = ANGLE_MIN;
         afin = ANGLE_MAX;
         ainc = ANGLE_INC;
         (void)free(anglims);
         anglims = NULL;
         fputs("Please specify range of angle limits in form", stderr);
         fputs(" start-finish:increment\n", stderr);
         anglims = get_str("s-f:n? ");
      }
   }

   /* Calculate total no of bond and angle intervals */
   nbond = floor(1.0*(bfin - bstart) / binc);
   nangle = floor(1.0*(afin - astart) / ainc);

   /* Start reading bond list file */
   fgets(dummy,sizeof(dummy),Fp);

   if( strncmp(dummy,"Bonds:",6) )
   {
       perror("Error: Incorrect file format");
       exit(0);
   }

   for( i = 0; i < 2; i++)
      fgets(dummy,sizeof(dummy),Fp);

   /* Read bond length data from input file */
   while (1)
   { 
      fscanf(Fp,"%s - %s %d - %s %lf",a,atom1,&b,atom2,&rij);

      if( !strncmp(a,"Angles:",7) )
         break;

      u = floor((10.0*rij - bstart)/binc);

      if(( u >=0 ) && (u < nbond))
      { 
         check = check_bond(&root_bond, atom1, atom2);
         if( check < 0)
         {
            bond = NEW(BOND);
            strncpy((bond->atom1), atom1, 3);
            strncpy((bond->atom2), atom2, 3);
            insert_data(&root_bond,bond,1);
         }
         else
         {
            node = select_node(&root_bond, check);
            bond = node->data; 
            if (node == NULL)
            {
               perror("Error: Failed to find bond record");
               exit(2);
            }
         }
         if( check_count(&(bond->root_count),u) < 0)
         {
            count = NEW(COUNT);
            count->xaxis = u;
            count->yaxis = 1;
            insert_data(&(bond->root_count),count,1);
         }
      }
   }

   /* Read angle data from input file */
   for( i = 0; i < 2; i++)
      fgets(dummy,sizeof(dummy),Fp);

   while ( !feof(Fp))
   { 
      fscanf(Fp,"%s - %s %d - %s %d - %s %lf %lf %lf\n",
             &a,&atom1,&b,atom2,&c,&atom3,&alpha,&r1,&r2);

      u = floor((alpha - astart)/ainc);

      if(( u >=0 ) && (u < nangle))
      { 
         check = check_angle(&root_angle, atom1, atom2, atom3);
         if( check < 0)
         {
            angle = NEW(ANGLE);
            strncpy((angle->atom1), atom1, 3);
            strncpy((angle->atom2), atom2, 3);
            strncpy((angle->atom3), atom3, 3);
            insert_data(&root_angle,angle,1);
         } 
         else
         { 
            node = select_node(&root_angle, check);
            angle = node->data;
            if (node == NULL)
            {
               perror("Error: Failed to find angle record");
               exit(2);
            }
         }
         if( check_count(&(angle->root_count),u) < 0)
         {
            count = NEW(COUNT);
            count->xaxis = u;
            count->yaxis = 1;
            insert_data(&(angle->root_count),count,1);
         }
      }
   }

   /* Write data to file in column format */
   display_BOND(&root_bond, bstart, binc);
   display_ANGLE(&root_angle, astart, ainc);

   return 0;
}
