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
 *  Revision Log
 *  $Log: mdbond.c,v $
 *  Revision 1.0  1999/06/24 11:05:12  craig 
 *  Initial revision
 *
 */
#include "defs.h"
#ifdef HAVE_STDARG_H
#include <stdarg.h>
#else
#include <varargs.h>
#endif
#include <errno.h>
#include <math.h>
#include "stdlib.h"
#include "stddef.h"
#include "string.h"
#include <stdio.h>
#include "structs.h"
#include "messages.h"
#include "list.h"
#ifdef HAVE_STDARG_H
gptr	*arralloc(size_mt,int,...); 	/* Array allocator		      */
#else
gptr	*arralloc();	        	/* Array allocator		      */
#endif

/*
 * Default limits for bond intervals and angle intervals - integers only
 */
#define BOND_MIN  2
#define BOND_MAX  20          /* Interparticle distances in tenths of Angstroms */
#define BOND_INC  5           /* BOND_INC dummy variable */

#define ANGLE_MIN  0
#define ANGLE_MAX  180        /* Angle intervals in degrees */
#define ANGLE_INC  1          /* ANGLE_INC dummy variable */

#define DOTPROD(x,y)   ((x[0]*y[0])+(x[1]*y[1])+(x[2]*y[2])) 
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

void	invert();
void	mat_vec_mul();
void	make_sites();
char	*strlower();
void	read_sysdef();
void	initialise_sysdef();
void	re_re_header();
void	re_re_sysdef();
void	allocate_dynamics();
void	lattice_start();
void	read_restart();
void	init_averages();
int	getopt();
gptr	*talloc();
FILE	*popen();
/*======================== Global vars =======================================*/
int ithread=0, nthreads=1;
/******************************************************************************
 * Dummies of 'moldy' routines so that mdbond may be linked with moldy library*
 ******************************************************************************/
void 	init_rdf()
{}
gptr *rdf_ptr()
{return 0;}
void new_lins()
{}
int lines_left()
{return 0;}
void new_page()
{}
void	new_line()
{
   (void)putchar('\n');
}
void	banner_page()
{}
void	note()
{}
void	conv_potentials()
{}
void	conv_control()
{}
/******************************************************************************
 *  message.   Deliver error message to possibly exiting.  It can be called   *
 *             BEFORE output file is opened, in which case out to stderr.     *
 ******************************************************************************/
#ifdef HAVE_STDARG_H
#   undef  va_alist
#   define      va_alist int *nerrs, ...
#   ifdef va_dcl
#      undef va_dcl
#   endif
#   define va_dcl /* */
#endif
/*VARARGS*/
void    message(va_alist)
va_dcl
{
   va_list      ap;
   char         *buff;
   int          sev;
   char         *format;
   static char  *sev_txt[] = {" *I* "," *W* "," *E* "," *F* "};
#ifdef HAVE_STDARG_H
   va_start(ap, nerrs);
#else
   int          *nerrs;

   va_start(ap);
   nerrs = va_arg(ap, int *);
#endif
   buff  = va_arg(ap, char *);
   sev   = va_arg(ap, int);
   format= va_arg(ap, char *);

   (void)fprintf(stderr, "mdbond: ");
   (void)vfprintf(stderr, format, ap);
   va_end(ap);
   fputc('\n',stderr);

   if(buff != NULL)                     /* null ptr means don't print buffer  */
   {
      (void)fprintf(stderr,"     buffer contents=\"%s\"",buff);
      fputc('\n',stderr);
   }
   if(sev >= ERROR && nerrs != NULL)
      (*nerrs)++;
   if(sev == FATAL)
      exit(3);
}

/******************************************************************************
 *  error.   Deliver error message to possibly exiting. 		      *
 ******************************************************************************/
#ifdef HAVE_STDARG_H
#undef  va_alist
#define	va_alist char *format, ...
#ifdef  va_dcl
#   undef  va_dcl
#endif
#define va_dcl /* */
#endif
/*VARARGS*/
void	error(va_alist)
va_dcl
{
   va_list	ap;
#ifdef HAVE_STDARG_H
   va_start(ap, format);
#else
   char		*format;

   va_start(ap);
   format= va_arg(ap, char *);
#endif

   (void)fprintf(stderr, "mdbond: ");
   (void)vfprintf(stderr, format, ap);
   fputc('\n',stderr);
   va_end(ap);

   exit(3);
}
/* strcpy with call to memory allocation */
static char * mystrdup(s)
char *s;
{
   char * t = NULL;
   if(s) t=malloc(strlen(s)+1);
   return t?strcpy(t,s):0;
}
/******************************************************************************
 * get_int().  Read an integer from stdin, issuing a prompt and checking      *
 * validity and range.  Loop until satisfied, returning EOF if appropriate.   *
 ******************************************************************************/
int get_int(prompt, lo, hi)
char	*prompt;
int	lo, hi;
{
   char		ans_str[80];
   int		ans_i, ans_flag;

   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, "%d", &ans_i) == 1 && ans_i >= lo && ans_i <= hi)
	 ans_flag++;
   }
   if( ans_flag )
      return(ans_i);
   else
      return(EOF);
}
/******************************************************************************
 * get_sym().  Read a character from stdin and match to supplied set	      *
 ******************************************************************************/
int get_sym(prompt, cset)
char	*prompt;
char	*cset;
{
   char		ans_c, ans_str[80];
   int		ans_flag;

   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, " %c", &ans_c) == 1 && strchr(cset, ans_c))
	 ans_flag++;
   }
   if( ans_flag )
      return(ans_c);
   else
      return(EOF);
}
/******************************************************************************
 * get_str().  Read a string from stdin, issuing a prompt.		      *
 ******************************************************************************/
char	*get_str(prompt)
char	*prompt;
{
   char		ans_str[80];
   char		*str = malloc(80);
   int		ans_flag;

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
/******************************************************************************
 * forstr.  Parse string str of format s-f:n (final 2 parts optional),        *
 *          returning integer values of s,f,n.                                *
 ******************************************************************************/
int
forstr(instr, start, finish, inc)
char	*instr;
int	*start, *finish, *inc;
{
   char	*p, *pp, *str = mystrdup(instr);
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
   if( *start > *finish || *start < 0 || *inc <= 0 )
   {
      fputs("Limits must satisfy", stderr);
      fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
      goto limerr;
   }
   return 0;
 limerr:
   return -1;
}
/******************************************************************************
 * dump_to_moldy.  Fill the 'system' arrays with the dump data in 'buf' (see  *
 * dump.c for format), expanding floats to doubles if necessary.              *
 ******************************************************************************/
#define DUMP_SIZE(level)  (( (level & 1) + (level>>1 & 1) + (level>>2 & 1) ) * \
			            (3*sys.nmols + 4*sys.nmols_r + 9)+ \
			     (level>>3 & 1) * \
			            (3*sys.nmols + 3*sys.nmols_r + 9) +\
			     (level & 1))
void
dump_to_moldy(buf, system)
float	*buf;
system_mt *system;
{
   int i;
   float	*c_of_m = buf;
   float	*quat   = buf+3*system->nmols;
   float	*h      = buf+3*system->nmols + 4*system->nmols_r;
   mat_mt hinv;

/* $dir no_recurrence */
   for(i = 0; i < system->nmols; i++)
   {
      system->c_of_m[i][0] = c_of_m[3*i];
      system->c_of_m[i][1] = c_of_m[3*i+1];
      system->c_of_m[i][2] = c_of_m[3*i+2]; 
   }
/* $dir no_recurrence */
   for(i = 0; i < system->nmols_r; i++)
   {
      system->quat[i][0] = quat[4*i];
      system->quat[i][1] = quat[4*i+1];
      system->quat[i][2] = quat[4*i+2];
      system->quat[i][3] = quat[4*i+3];
   }
/* $dir no_recurrence */
   for(i = 0; i < 3; i++)
   {
      system->h[i][0] = h[3*i];
      system->h[i][1] = h[3*i+1];
      system->h[i][2] = h[3*i+2];
   }
}
/******************************************************************************
 * morethan_BOND(). Compare distances stored in BOND structure types          *
 ******************************************************************************/
NODE  *morethan_BOND(root, data)
ROOT    **root;
BOND     *data;
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
           bd = node->data;
       }
   }
   return (node);
}
/******************************************************************************
 * morethan_ANGLE(). Compare angles stored in ANGLE structure types           *
 ******************************************************************************/
NODE *morethan_ANGLE(root, data)
ROOT    **root;
ANGLE    *data;
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
           ang = node->data;
       }
   }
   return (node);
}
/******************************************************************************
 * angle_calc(). Calculate angle (in degrees) between two vectors             *
 ******************************************************************************/
double angle_calc(vec1, vec2, a, b)
vec_mt    vec1, vec2;
double    a, b;
{
double    dp, angle;

   dp = DOTPROD(vec1,vec2)/a/b;

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
bond_calc(system, species, site_info, broot, aroot, sp_range, blim, alim)
system_mt	*system;
spec_mt		species[];
site_mt         site_info[];
ROOT		**broot, **aroot;
int		blim[3], alim[3];
int		sp_range[3];
{
   BOND		*bond;
   ANGLE        *angle;
   spec_mt	*spec1, *spec2, *spec3;
   double	dist1, dist2;
   double	tmp_angle;
   vec_mt	point1, point2, point3;
   vec_mt	vec1, vec2;
   vec_mt	shift, frac;
   vec_mt	shift2, frac2;
   vec_mt	a, min, max;
   int		i, j, k, u;
   int		nmoli, nmolj, nmolk;
   int		flag;
   NODE		*node;
   mat_mp	h = system->h;

/* Determine no of cells to search through */
   for( u = 0; u < 3; u++)
   {
      a[u] = sqrt(SQR(h[0][u]) + SQR(h[1][u]) + SQR(h[2][u]));
      max[u] = ceil(blim[1]/10.0/a[u]);
      min[u] = -1.0 * max[u];
   }

/* Scan through selected species and determine distances and angles within limits */
   nmoli = 0;

   for(spec1 = species+sp_range[0]; spec1 <= species+sp_range[1]; spec1 += sp_range[2])
     for(i=0; i < spec1->nmols; i++)
     {
        nmoli++;
        nmolj = 0;

        for(spec2 = species+sp_range[0]; spec2 <= spec1; spec2 += sp_range[2])
           for(j=0; j < (spec2==spec1?i:spec2->nmols); j++)
           {
              nmolj++;
              flag = 0;

              for(u = 0; u < 3; u++)
                 point1[u] = spec1->c_of_m[i][u];

              for( frac[0] = min[0]; frac[0] <= max[0]; frac[0]++) 
              for( frac[1] = min[1]; frac[1] <= max[1]; frac[1]++) 
              for( frac[2] = min[2]; frac[2] <= max[2]; frac[2]++) 
              {
                  mat_vec_mul(h, frac, shift, 1);

                  for( u = 0; u < 3; u++)
                     point2[u] = spec2->c_of_m[j][u] + shift[u]; 

                  dist1 = DISTANCE(point2, point1);

                  if( (dist1 >= blim[0]/10.0) && (dist1 <= blim[1]/10.0) )
                  { 
                      flag = 1;
                      bond = NEW(BOND); /* Create new bond record */
                      strncpy((bond->atom1),site_info[spec2->site_id[0]].name, 3);
                      strncpy((bond->atom2),site_info[spec1->site_id[0]].name, 3);
                      bond->number1 = nmolj;
                      bond->number2 = nmoli;
	              bond->length = dist1;
	  	      node = morethan_BOND(broot,bond);
	              if( (NUM(broot) == 0) || (node == NULL) )
                         insert_data(broot,bond,1);
                      else
                         insert_at_position(broot,node,bond,0);

                      nmolk = nmolj;

                      /* Calculate angle about spec1 molecule */
                      for(spec3 = spec2; spec3 <= species+sp_range[1]; spec3 += sp_range[2])
                         for( k=(spec3==spec2?j+1:0); k < spec3->nmols; k++)
                         {
                            nmolk++;

                            if( (nmolk != nmoli) && (nmolj != nmoli) )
                            {
                               for( frac2[0] = min[0]; frac2[0] <= max[0]; frac2[0]++) 
                               for( frac2[1] = min[1]; frac2[1] <= max[1]; frac2[1]++) 
                               for( frac2[2] = min[2]; frac2[2] <= max[2]; frac2[2]++) 
                               {
                                  mat_vec_mul(h, frac2, shift2, 1);

                                  for( u = 0; u < 3; u++)
                                     point3[u] = spec3->c_of_m[k][u] + shift2[u]; 

                                  dist2 = DISTANCE(point3, point1);

                                  if( (dist2 >= blim[0]/10.0 ) && (dist2 <= blim[1]/10.0) )
                                  {
                                     for(u = 0; u < 3; u++)
                                     {
                                        vec1[u] = point2[u] - point1[u];  /* Vector between atoms 2 and 1 */
                                        vec2[u] = point3[u] - point1[u];  /* Vector between atoms 3 and 1 */
                                     }
                                     tmp_angle = angle_calc(vec1, vec2, dist1, dist2);
                          
                                     if( (tmp_angle >= alim[0]) && (tmp_angle <= alim[1]) )
                                     {
                                        angle = NEW(ANGLE); /* Create new angle record */
                                        strncpy((angle->atom1),site_info[spec1->site_id[0]].name, 3);
                                        strncpy((angle->atom2),site_info[spec2->site_id[0]].name, 3);
                                        strncpy((angle->atom3),site_info[spec3->site_id[0]].name, 3);
                                        angle->number1 = nmoli;
                                        angle->number2 = nmolj;
                                        angle->number3 = nmolk;
                                        angle->length1 = dist1;
		                        angle->length2 = dist2;
		                        angle->value = tmp_angle;
		                        node = morethan_ANGLE(aroot,angle);
			                if( (NUM(aroot) == 0) || (node == NULL) )
                                           insert_data(aroot,angle,1);
                                        else
                                           insert_at_position(aroot,node,angle,0);
                                     }
                                  }
                               }
                            }
                         }
                  }
              }
              if( flag )           /* Calculate angle about spec2 molecule */ 
              {
                 nmolk = nmoli;
                 for(u = 0; u < 3; u++)
                     point2[u] = spec2->c_of_m[j][u];

                 for(spec3 = spec1; spec3 <= species+sp_range[1]; spec3 += sp_range[2])
                    for( k=(spec3==spec1?i+1:0); k < spec3->nmols; k++)
                    {
                        nmolk++;

                        if( (nmolk != nmolj) && (nmoli != nmolj) )
                        {
                            for( frac[0] = min[0]; frac[0] <= max[0]; frac[0]++) 
                            for( frac[1] = min[1]; frac[1] <= max[1]; frac[1]++) 
                            for( frac[2] = min[2]; frac[2] <= max[2]; frac[2]++) 
                            {
                                mat_vec_mul(h, frac, shift, 1);
 
                                for( u = 0; u < 3; u++)
                                   point1[u] = spec1->c_of_m[i][u] + shift[u]; 
 
                                dist1 = DISTANCE(point1, point2);
  
                                if( (dist1 >= blim[0]/10.0 ) && (dist1 <= blim[1]/10.0) )
                                {
                                    for( frac2[0] = min[0]; frac2[0] <= max[0]; frac2[0]++) 
                                    for( frac2[1] = min[1]; frac2[1] <= max[1]; frac2[1]++) 
                                    for( frac2[2] = min[2]; frac2[2] <= max[2]; frac2[2]++) 
                                    {
                                        mat_vec_mul(h, frac2, shift2, 1);

                                        for( u = 0; u < 3; u++)
                                           point3[u] = spec3->c_of_m[k][u] + shift2[u]; 
   
                                        dist2 = DISTANCE(point3, point2);
   
                                        if( (dist2 >= blim[0]/10.0 ) && (dist2 <= blim[1]/10.0) )
                                        {
                                            for(u = 0; u < 3; u++)
                                            {
                                                vec1[u] = point1[u] - point2[u];  /* Vector between atoms 1 and 2 */
                                                vec2[u] = point3[u] - point2[u];  /* Vector between atoms 3 and 2 */
                                            }
   
                                            tmp_angle = angle_calc(vec1, vec2, dist1, dist2);
                           
                                            if( (tmp_angle >= alim[0]) && (tmp_angle <= alim[1]) )
                                            {
                                                angle = NEW(ANGLE);  /* Calculate new angle record */
                                                strncpy((angle->atom1),site_info[spec2->site_id[0]].name, 3);
                                                strncpy((angle->atom2),site_info[spec1->site_id[0]].name, 3);
                                                strncpy((angle->atom3),site_info[spec3->site_id[0]].name, 3);
                                                angle->number1 = nmolj;
                                                angle->number2 = nmoli;
                                                angle->number3 = nmolk;
                                                angle->length1 = dist1;
                                                angle->length2 = dist2;
		                                angle->value = tmp_angle;
		                                node = morethan_ANGLE(aroot,angle);
			                        if( (NUM(aroot) == 0) || (node == NULL) )
                                                   insert_data(aroot,angle,1);
                                                else
                                                   insert_at_position(aroot,node,angle,0);
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
/******************************************************************************
 * data_out().  Output bonds and angles to file with same format as Shell by  *
 *              N. Allan and G. Barrera, Bristol Univ.                        *
 ******************************************************************************/
void data_out(broot, aroot)
ROOT    **broot;
ROOT    **aroot;
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
	  printf("             %3d -%4s   %3d -%4s   %12.7f\n",
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
          printf(" %3d -%4s   %3d -%4s   %3d -%4s   %11.6f    %12.7f  %12.7f\n",
              ang->number1,ang->atom1,ang->number2,ang->atom2,ang->number3,
                  ang->atom3,ang->value, ang->length1,ang->length2);
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
contr_mt		control;

int
main(argc, argv)
int	argc;
char	*argv[];
{
   int	c, cflg = 0, ans_i, sym = 0, data_source = 0;
   char 	line[80];
   extern char	*optarg;
   int		errflg = 0;
   int		intyp = 0;
   int		start, finish, inc;
   int		rflag = 0, sflag = 0;
   int		bflag = 0, aflag = 0;
   int		irec;
   char         *bondlims = NULL, *anglims = NULL;
   char		*filename = NULL, *dump_name = NULL;
   char		*insert = NULL, *speclims = NULL;
   char		*dumplims = NULL, *tempname;
   char		dumpcommand[256];
   int		dump_size;
   float	*dump_buf;
   FILE		*Fp, *Dp;
   restrt_mt	restart_header;
   system_mt	sys;
   spec_mt	*species;
   site_mt	*site_info;
   pot_mt	*potpar;
   quat_mt	*qpf;
   contr_mt	control_junk;
   int          av_convert;

   int		sp_range[3];              /* Range and increment for species selection */
   int          blim[3], alim[3];         /* Min and max values for bonds and angles */
   ROOT         *root_bond = NULL;        /* Root of bond linked list */
   ROOT         *root_angle = NULL;       /* Root of angle linked list */

#define MAXTRY 100
   control.page_length=1000000;

   while( (c = getopt(argc, argv, "r:s:d:t:g:o:b:a:") ) != EOF )
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
	 dump_name = optarg;
	 break;
       case 't':
	 dumplims = mystrdup(optarg);
         break;
       case 'g':
	 speclims = mystrdup(optarg);
	 break;
       case 'o':
	 if( freopen(optarg, "w", stdout) == NULL )
	    error("failed to open file \"%s\" for output", optarg);
	 break;
       case 'b':
         bondlims = mystrdup(optarg);
         break;
       case 'a':
	 anglims = mystrdup(optarg);
         break;
       default:
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fputs("Usage: mdbond [-s sys-spec-file|-r restart-file] ",stderr);
      fputs("[-d dump-file/s] [-t timeslice]] [-g species] ",stderr);
      fputs("[-b bond limits] [-a angle limits] [-o output-file]\n",stderr);
      exit(2);
   }

   if( dump_name )
      data_source = 'd';

   if(intyp == 0)
   {
      fputs("How do you want to specify the simulated system?\n", stderr);
      fputs("Do you want to use a system specification file (1)", stderr);
      fputs(" or a restart file (2)\n", stderr);
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
	 error("Couldn't open restart file \"%s\" for reading -\n%s\n", 
	       filename, strerror(errno));
      re_re_header(Fp, &restart_header, &control_junk);
      re_re_sysdef(Fp, &sys, &species, &site_info, &potpar);
      break;
    default:
      error("Internal error - invalid input type", "");
   }
   allocate_dynamics(&sys, species);

   /* Use default values for species selection limits */
   sp_range[0] = 0;
   sp_range[1] = sys.nspecies;
   sp_range[2] = 1;

   if( speclims == NULL)
        sflag ++;
    
   while (!sflag) 
   {
      if( forstr(speclims, &(sp_range[0]), &(sp_range[1]), &(sp_range[2])))
      {  
          fputs("Invalid range for molecule selection \"", stderr);
	  fputs(speclims, stderr);
	  fputs("\"\n", stderr);
      }
      else
          sflag++;
      if( sp_range[0] > sp_range[1] || sp_range[0] < 0 || sp_range[2] <= 0 )
      {
         fputs("Species limits must satisfy", stderr);
         fputs(" finish >= start, start >= 0 and increment > 0\n", stderr);
         sflag = 0;
      }
      if( sp_range[1] > sys.nspecies)
      {
         fputs("Molecule selection exceeds no. of species\n",stderr);
         sflag = 0;
      }
      if( !sflag )
      {
         sp_range[0] = 0;
         sp_range[1] = sys.nspecies;
         sp_range[2] = 1;
         (void)free(speclims);
         speclims = NULL;
         fputs("Please specify molecule selection in form", stderr);
         fputs(" start-finish:increment\n", stderr);
         speclims = get_str("s-f:n? ");
      }
   }
   /* Set default values for bond limits (x10) */
   blim[0] = BOND_MIN;
   blim[1] = BOND_MAX;
   blim[2] = BOND_INC;

   if( bondlims == NULL )
      bflag++;

   /* Input and check bond length limits where necessary */
   while( !bflag)
   {
      if( forstr(bondlims, &(blim[0]), &(blim[1]), &(blim[2])) )
      {
          fputs("Invalid range for bond lengths \"", stderr);
          fputs(bondlims, stderr);
          fputs("\"\n", stderr);
      }
      else
         bflag++;
      if( blim[0] > blim[1] || blim[0] < 0 )
      {
         fputs("Bond length limits must satisfy", stderr);
         fputs(" finish >= start and start >= 0\n", stderr);
         bflag = 0;
      }
      if( !bflag)
      {
         blim[0] = BOND_MIN;
         blim[1] = BOND_MAX;
         blim[2] = BOND_INC;
         (void)free(bondlims);
         bondlims = NULL;
         fputs("Please specify range of bond limits in form", stderr);
         fputs(" start-finish\n", stderr);
         bondlims = get_str("s-f:n? ");
      }
   }

   /* Set default values for angle limits */
   alim[0] = ANGLE_MIN;
   alim[1] = ANGLE_MAX;
   alim[2] = ANGLE_INC;

   if( anglims == NULL )
       aflag++;

   /* Input and check angle limits where necessary */
   while (!aflag)
   {
      if( forstr(anglims, &(alim[0]), &(alim[1]), &(alim[2])) )
      {
         fputs("Invalid range for angles \"", stderr);
         fputs(anglims, stderr);
         fputs("\"\n", stderr);
      }
      else
         aflag++;
      if( alim[0] > alim[1] || alim[0] < 0 <= 0 )
      {
         fputs("Angle limits must satisfy", stderr);
         fputs(" finish >= start and start >= 0\n", stderr);
         aflag=0;
      }
      if( !aflag)
      {
         alim[0] = ANGLE_MIN;
         alim[1] = ANGLE_MAX;
         alim[2] = ANGLE_INC;
         (void)free(anglims);
         anglims = NULL;
         fputs("Please specify range of angle limits in form", stderr);
         fputs(" start-finish\n", stderr);
         anglims = get_str("s-f:n? ");
      }
   }

   switch(data_source)                  /* To read configurational data       */
   {
    case 's':                           /* Lattice_start file                 */
      lattice_start(Fp, &sys, species, qpf);
      break;
    case 'r':                           /* Restart file                       */
      init_averages(sys.nspecies, restart_header.vsn,
                    control_junk.roll_interval, control_junk.roll_interval,
                    &av_convert);
      read_restart(Fp, restart_header.vsn, &sys, av_convert);
      break;
    case 'd':
      if( dump_name == 0 )
      {
         fputs("Enter canonical name of dump files (as in control)\n",stderr);
         if( (dump_name = get_str("Dump file name? ")) == NULL)
            exit(2);
      }
  /*
   *  Ensure that the dump limits start, finish, inc are set up,
   *  either on command line or by user interaction.
   */
      do
      {
         rflag = 0;
         start = finish = 0;
         inc = 1;
         if( dumplims == NULL )
         {
             fputs("Please specify range of dump records in form", stderr);
             fputs(" start-finish:increment\n", stderr);
             dumplims = get_str("s-f:n? ");
         }
         if( forstr(dumplims, &start, &finish, &inc) )
         {
            rflag++;
            fputs("Invalid range for dump records \"", stderr);
            fputs(dumplims, stderr);
            fputs("\"\n", stderr);
         }
         if( rflag )
         {
            (void)free(dumplims);
            dumplims = NULL;
         } 
      } while(rflag);
  /*
   * Allocate buffer for data
   */
      dump_size = DUMP_SIZE(~0)*sizeof(float);

      if( (dump_buf = (float*)malloc(dump_size)) == 0)
         error("malloc failed to allocate dump record buffer (%d bytes)",
            dump_size);
#if defined (HAVE_POPEN)
     sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d %s",
        sys.nmols, sys.nmols_r, start, finish, inc, dump_name);
   
      if( (Dp = popen(dumpcommand,"r")) == 0)
         error("Failed to execute \'dumpext\" command - \n%s",
            strerror(errno));
#else
      tempname = tmpnam((char*)0);
      sprintf(dumpcommand,"dumpext -R%d -Q%d -b -c 0 -t %d-%d:%d -o %s %s",
          sys.nmols,sys.nmols_r, start, finish, inc, tempname, dump_name);
      system(dumpcommand);
      if( (Dp = fopen(tempname,"rb")) == 0)
         error("Failed to open \"%s\"",tempname);
#endif

      for(irec = start; irec <= finish; irec+=inc)
      {
         if( fread(dump_buf, dump_size, 1, Dp) < 1 || ferror(Dp) )
            error("Error reading record %d in dump file - \n%s\n",
               irec, strerror(errno));

         dump_to_moldy(dump_buf, &sys);  /*read dump data */

#ifdef DEBUG
      fprintf(stderr,"Sucessfully read dump record %d from file  \"%s\"\n",
          %header.maxdumps, dump_name);
#endif
         /* Perform bond/angle calculations for each slice of dump file */
         bond_calc(sys, species, site_info, &root_bond, &root_angle, sp_range, blim, alim); 
         printf("- Time slice %d -\n",irec);
         data_out(&root_bond, &root_angle);
         putchar('\n');
         if( delete_list(&root_bond))
            error("Error releasing bond list data for slice %d - \n%s\n",
               irec, strerror(errno));
         if( delete_list(&root_angle))
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
      mat_vec_mul(sys.h, sys.c_of_m, sys.c_of_m, sys.nmols);

      bond_calc(sys, species, site_info, &root_bond, &root_angle, sp_range, blim, alim); 
      data_out(&root_bond,&root_angle);
   }
   return 0;    
}