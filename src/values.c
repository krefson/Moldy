/******************************************************************************
 *  values.  	Calculate various thermodynamic quantities and their averages *
 *  		ant print periodic output.	Contents:		      *
 * init_averages()	Allocate and set up averages database		      *
 * av_ptr()		Return pointer and size of database		      *
 * add_average()	Update an average with a new datum		      *
 * value()		Return current value of specified quantity	      *
 * roll_av()		Return rolling averagee of specified quantity	      *
 * roll_sd()		Return standard deviation of rolling average	      *
 * values()		Calculate all values and update database	      *
 * print_frame()	Print periodic output in formatted form		      *
 * output()		Call print_frame for values, rolling averages & sd's  *
 * averages()		Calculate and print averages()			      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	values.c,v $
 * Revision 1.12  92/03/24  12:41:31  keith
 * Fixed bug introduced in last revision which calculated
 * averages wrongly.
 * Added code to zero averages database if reset-averages
 * is set OR if the info in restart averages database is
 * out of date.  (ie ig begin-averages is within current run).
 * 
 * Revision 1.11  92/03/19  15:45:42  keith
 * Added support for dynamic allocation of rolling average arrays,
 * conversion of existing restart files is done on fly.
 * 
 * Revision 1.10  91/08/19  16:48:51  keith
 * Modifications for better ANSI/K&R compatibility and portability
 * --Changed sources to use "gptr" for generic pointer -- typedefed in "defs.h"
 * --Tidied up memcpy calls and used struct assignment.
 * --Moved defn of NULL to stddef.h and included that where necessary.
 * --Eliminated clashes with ANSI library names
 * --Modified defs.h to recognise CONVEX ANSI compiler
 * --Modified declaration of size_t and inclusion of sys/types.h in aux.c
 *   for GNU compiler with and without fixed includes.
 * 
 * 
 * Revision 1.9  91/03/12  15:43:29  keith
 * Tidied up typedefs size_t and include file <sys/types.h>
 * Added explicit function declarations.
 * 
 * Revision 1.8  90/05/16  18:40:53  keith
 * Renamed own freer from cfree to tfree.
 * 
 * Revision 1.7  90/05/16  14:20:40  keith
 * *** empty log message ***
 * 
 * Revision 1.6  90/05/08  17:19:09  keith
 * Fixed bug which got indexing of av_info[] wrong in test for -ve variance.
 * 
 * Revision 1.5  90/04/16  18:17:19  keith
 * Rearranged expression as workaround for CRAY CC4.1 bug.
 * 
 * Revision 1.4  89/06/22  15:45:29  keith
 * Tidied up loops over species to use one pointer as counter.
 * 
 * Revision 1.3  89/06/01  21:25:40  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.2  89/05/19  10:35:28  keith
 * Fixed bug which printed to stdout rather than control.out.
 * 
 * Revision 1.1  89/04/20  16:00:58  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: /home/eeyore/keith/md/moldy/RCS/values.c,v 1.12 92/03/24 12:41:31 keith Exp $";
#endif
/*========================== Program include files ===========================*/
#include	"defs.h"
#include	"structs.h"
#include	"messages.h"
/*========================== Library include files ===========================*/
#include	<stdio.h>
#include	<math.h>
#include	"stddef.h"
#include	"string.h"
/*========================== External function declarations ==================*/
void	mat_vec_mul();
void	q_conj_mul();			/* Quaternion multiply conjugated     */
double	det();				/* Determinant of 3x3 matrix	      */
double	vdot();				/* Vector dot product		      */
double	sum();				/* Vector sum			      */
void	zero_real();
void	energy_dyad();
double	trans_ke();
double	rot_ke();
double	precision();			/* Machine precision constant	      */
void	message();			/* Error message and exit handler     */
void	note();
char	*atime();
void	new_line();
void	new_lins();
void	new_page();
gptr            *talloc();	       /* Interface to memory allocator       */
void            tfree();	       /* Free allocated memory	      	      */
/*========================== External data references ========================*/
extern	contr_t	control;
extern	int	out_line, out_page;	/* Current line and page in output    */
/*========================== Structs local to module =========================*/
#define	MAX_ROLL_INTERVAL	100
typedef struct
{   double	value,
   		sum,
   		sum_sq,
		mean,
		sd,
   		roll[MAX_ROLL_INTERVAL],
   		roll_mean,
   		roll_sd;
} old_av_t;

typedef	union
{
   old_av_t		av;
   struct
   {	
      int	av, roll;
   }		cnt;
} old_av_u_t;

typedef struct
{   double	value,
   		sum,
   		sum_sq,
		mean,
		sd,
   		roll[1];
} av_t;

typedef	struct
{	
   int		nav, 
   		nroll, 
		iroll, 
		pad;
   double align;
} av_head_t;

typedef struct
{
   av_n		id;
   char		*name, *unit;
   int		field_width;
   char		*format;
   int		mult;
   av_t		**p;
} av_info_t;
/*========================== Global variables ================================*/
av_info_t av_info[] = { {tke_n, "Trans KE",   CONV_E_N,	11, "%11.5g",-1, NULL},
			{rke_n, "Rot KE",     CONV_E_N,	11, "%11.5g",-1, NULL},
			{pe_n,  "Pot Energy", CONV_E_N,	11, "%11.5g",NPE,NULL},
			{e_n,   "Tot Energy", CONV_E_N,	11, "%11.5g", 1, NULL},
			{tt_n,  "TTemp",      CONV_T_N,	 6, "%6.1f", -1, NULL},
			{rt_n,  "RTemp",      CONV_T_N,	 6, "%6.1f", -1, NULL},
			{t_n,   "Temp",	      CONV_T_N,	 6, "%6.1f", 1,  NULL},
			{h0_n,  "h(1,*)",     LUNIT_N,	 6, "%6.2f", 3,  NULL},
			{h1_n,  "h(2,*)",     LUNIT_N,	 6, "%6.2f", 3,  NULL},
			{h2_n,  "h(3,*)",     LUNIT_N,	 6, "%6.2f", 3,  NULL},
			{stress0_n,"Stress",  CONV_P_N,	10, "%10.3g",3,  NULL},
			{stress1_n,"Stress",  CONV_P_N,	10, "%10.3g",3,  NULL},
			{stress2_n,"Stress",  CONV_P_N,	10, "%10.3g",3,  NULL},
			{press_n,"Pressure",  CONV_P_N,	10, "%10.3g",1,  NULL},
			{vir_n, "Virial",     CONV_E_N,	11, "%11.5g",1,  NULL},
			{msqf_n,"<F**2>",     CONV_F_N,	10, "%10.5g",-3, NULL},
			{msqt_n,"<N**2>",     CONV_N_N,	10, "%10.5g",-3, NULL},
			{dip_n, "Dip Mom",    CONV_D_N,	 8, "%8.2g", 3,  NULL}};
static  int	av_size;		/* Size of averages database          */
static  int	av_t_size;		/* Size of entry inaverages database  */
static  av_head_t  *av_head;
static	av_t	*av;			/* Dynamic array of averages structs  */
static	int	navs = 0;		/* Size of array av                   */
static	int	max_row = 0;		/* Largest number of components       */
static  int	max_col = (int)press_n;	/* Number to print across page        */
static  int	av_convert = 0;
static  int     av_tmp_size;
static  gptr    *av_tmp;
/*========================== Macros ==========================================*/
#define NAVT			(int)end
#define INC(av_p)    (av_p = (av_t*)((char*)av_p + av_t_size))
/*============================================================================*/
/******************************************************************************
 *  init_averages  Allocate space for and initialise the averages database.   *
 *  The first word is the 'count' of numbers summed to date, (Pointed at by   *
 *  global av_cnt) which is followed by an array of av_t structs of dimension *
 *  navs. 								      *
 *  Struct array av_info is set up with .av_p pointing to the appropriate     *
 *  entry in av.  The multiplicity (number of components) is also set by the  *
 *  following rule - a positive entry is the true multiplicity and a negative *
 *  one is multiplied by the number of species for the true multiplicity.     *
 ******************************************************************************/
void	init_averages(nspecies, vsn, roll_interval, old_roll_interval)
int	nspecies;
char	*vsn;
int	roll_interval, old_roll_interval;
{
   av_t		*av_p;
   int		i, imult;
   int		major, minor, cmajor=1, cminor=13;

   for(i = 0; i < (int)end; i++)	/* cycle over enum types av_n         */
   {
      if(av_info[i].mult < 0)		/* Set true multiplicity of each type */
         av_info[i].mult = -nspecies * av_info[i].mult;
      navs += av_info[i].mult;		/* Count total			      */
      if(i < max_col && av_info[i].mult > max_row)
         max_row = av_info[i].mult;
   }

   /* Determine size of database, Allocate space and set pointers             */
   av_t_size = sizeof(av_t)+(roll_interval-1)*sizeof(double);
   av_size = sizeof(av_head_t) + navs*av_t_size;
   av_head  = (av_head_t*)aalloc( av_size, char);
   av_head->nav = av_head->nroll = av_head->iroll = 0;
   av       = (av_t*)(av_head+1);
   zero_double(av, navs*av_t_size/sizeof(double));

   av_p = av;
   for(i = 0; i < (int)end; i++)	/* Set up pointers to area of array   */
   {					/* reserved for each type, size=mult. */
      av_info[i].p = aalloc(av_info[i].mult, av_t*);
      for(imult = 0; imult < av_info[i].mult; imult++)
      {
	 av_info[i].p[imult] = av_p;
	 INC(av_p);
      }
   } 
   /*
    * Do we have to do any conversion on averages read from restart file?
    * We just allocate buffers and set flags here.
    */
   if( vsn )
   {
      /*
       * First check whether restart was written by 1.13 or earlier.
       */
      if( sscanf(vsn, "%d.%d", &major, &minor) < 2 )
	 message(NULLI, NULLP, FATAL, INRVSN, vsn);
      if( major < cmajor || (major==cmajor && minor <= cminor ) )
      {
	 av_tmp_size = (navs+1)*sizeof(old_av_u_t);
	 av_tmp = aalloc(av_tmp_size, char);
	 av_convert = 1;
      }
      /*
       * Has size of rolling average store changed?
       */
      else if (roll_interval != old_roll_interval )
      {
	 av_tmp_size =  sizeof(av_head_t) 
	             + navs*(sizeof(av_t)+(old_roll_interval-1)*sizeof(double));
	 av_tmp = aalloc(av_tmp_size, char);
	 av_convert = 2;	 
      }
   }
}
/******************************************************************************
 * convert_averages.  Update averages database if roll_interval changed or    *
 * if restart file written using old "static" scheme.			      *
 * Also a convenient place to implement reset_averages.			      *
 ******************************************************************************/
void	convert_averages(roll_interval, old_roll_interval)
int	roll_interval, old_roll_interval;
{
   int iav, old_nroll, old_iroll, rbl, prev_av_t_size;
   old_av_u_t *old_av_p=(old_av_u_t *)av_tmp;
   av_head_t	*prev_av_head = (av_head_t *)av_tmp;
   av_t	        *av_p, *prev_av_p;

   switch(av_convert)
   {
    case 0:					/* No conversion needed       */
      break;
    case 1:					/* Convert from static scheme */
      old_nroll = old_av_p[0].cnt.roll;
      old_iroll = control.istep % old_roll_interval;
      av_head->nroll = MIN(old_nroll,roll_interval);
      av_head->iroll = av_head->nroll % roll_interval;
      av_head->nav   = old_av_p[0].cnt.av;
      rbl = MIN(old_iroll, av_head->nroll);
      old_av_p++;
      av_p = av_info[0].p[0];
      for(iav = 0; iav < navs; iav++)
      {
	 av_p->value  = old_av_p->av.value;
	 av_p->sum    = old_av_p->av.sum;
	 av_p->sum_sq = old_av_p->av.sum_sq;
	 av_p->mean   = old_av_p->av.mean;
	 av_p->sd     = old_av_p->av.sd;
	 (void)memcpy((gptr*)(av_p->roll + av_head->nroll - rbl),
		      (gptr*)(old_av_p->av.roll + old_iroll - rbl),
		      rbl*sizeof(double));
	 (void)memcpy((gptr*)av_p->roll, 
		      (gptr*)(old_av_p->av.roll+old_nroll-av_head->nroll+rbl),
		      (av_head->nroll-rbl)*sizeof(double));
	 INC(av_p);
	 old_av_p++;
      }
      break;
    case 2:					/* Change roll_interval       */
      prev_av_t_size = sizeof(av_t)+(old_roll_interval-1)*sizeof(double);
      old_nroll = prev_av_head->nroll;
      old_iroll = prev_av_head->iroll;
      av_head->nroll = MIN(old_nroll,roll_interval);
      av_head->iroll = av_head->nroll % roll_interval;
      av_head->nav   = prev_av_head->nav;
      rbl = MIN(old_iroll, av_head->nroll);
      prev_av_p = (av_t *)(prev_av_head+1);
      av_p = av_info[0].p[0];
      for(iav = 0; iav < navs; iav++)
      {
	 /*
	  * Can do a struct copy -- will only pick up 1st roll entry
	  */
	 *av_p = *prev_av_p;
	 (void)memcpy((gptr*)(av_p->roll + av_head->nroll - rbl),
		      (gptr*)(prev_av_p->roll + old_iroll - rbl),
		      rbl*sizeof(double));
	 (void)memcpy((gptr*)av_p->roll, 
		      (gptr*)(prev_av_p->roll+old_nroll-av_head->nroll+rbl),
		      (av_head->nroll-rbl)*sizeof(double));

	 INC(av_p);
	 prev_av_p = (av_t*)((char*)prev_av_p + prev_av_t_size);
      }
      break;
   }
   av_convert = 0;
   /*
    *  Reset averages and counters to zero if a) requested
    *  or b) we have not yet reached begin_average. (The latter 
    *  avoids the situation of non-contiguous averages).
    *  Yes I know we just set them up, but this way is clearest.
    */
   if( control.reset_averages || control.istep+1 <= control.begin_average )
   {
      av_p = av_info[0].p[0];
      for(iav = 0; iav < navs; iav++)
      {
	 av_p->sum = av_p->sum_sq = av_p->mean = av_p->sd = 0.0;
	 INC(av_p);
      }
      av_head->nav = 0;
   }
}
/******************************************************************************
 * av_ptr   Return a pointer to averages database and its size (for restart)  *
 ******************************************************************************/
gptr	*av_ptr(size)
size_t	*size;
{
   switch(av_convert)
   {
    case 0:
      *size = av_size;
      if(av_head != NULL)
	 return((gptr*)av_head);
      break;
    case 1:
    case 2:
      *size = av_tmp_size;
      if(av_tmp != NULL)
	 return(av_tmp);
   }
   message(NULLI, NULLP, FATAL, AVNOC, "av_ptr");
   return(NULL);					/* To satisfy lint    */
}
/******************************************************************************
 * add_average  update the averages database with new datum                   *
 ******************************************************************************/
static void	add_average(datum, type, offset)
double	datum;				/* Datum to store and accumulate sums */
av_n	type;				/* What kind (ie where to store)      */
int	offset;				/* Sub-type or which component        */
{
   av_t		*av_p;
   if(offset < 0 || offset > av_info[(int)type].mult - 1)
      message(NULLI, NULLP, FATAL, AVBNDS, offset, av_info[(int)type].name);

   av_p = av_info[(int)type].p[offset];
 
   av_p->value = datum;
   if(control.istep >= control.begin_average)
   {
      av_p->sum += datum;
      av_p->sum_sq += datum * datum;
   }
   av_p->roll[av_head->iroll] = datum;
}
/******************************************************************************
 * values   Calculate the values of the thermodynamic quantities, maintain and*
 * if necessary print them, their averages and rolling averages.              *
 ******************************************************************************/
void	values(system, species, meansq_f_t, pe, dipole, stress_vir)
system_p	system;		/* record of system info		      */
spec_t	species[];	/* Records of info for each species	      */
vec_t		meansq_f_t[][2];/* mean square forces and torques             */
double		pe[];		/* potential energy real/reciprocal space     */
vec_t		dipole;		/* dipole moment of whole system              */
mat_t		stress_vir;	/* 'Potential' part of stress, or virial      */
{
   spec_p	spec;
   int		ispec, ipe;
   double	e, tot_ke = 0.0, tot_pe = 0.0;
   int		i, j, k;
   mat_t	ke_dyad,
                stress;
   double	vol = det(system->h);

   for(ipe = 0; ipe < NPE; ipe++)
   {
      add_average(CONV_E * pe[ipe], pe_n, ipe);
      tot_pe += pe[ipe];
   }

   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
      ispec = spec-species;
      e = trans_ke(system->h, spec->vel, spec->mass, spec->nmols);
      add_average(CONV_E * e, tke_n, ispec);	/* c of m  kinetic energy     */
      tot_ke += e;

      add_average(e/(1.5*spec->nmols*kB), tt_n, ispec);
						/* c of mass temp.	      */
      if(spec->rdof > 0)			/* Only if polyatomic species */
      {
         e = rot_ke(spec->quat, spec->qdot, spec->inertia, spec->nmols);
         add_average(CONV_E * e, rke_n, ispec);	/* Rotational kinetic energy  */
         tot_ke += e;
         add_average(e/(0.5*kB*spec->rdof*spec->nmols), rt_n, ispec);
      }						/* and temperature            */
   }   						/* Overall temperature        */
   add_average(CONV_E*(tot_ke+tot_pe), e_n, 0);	/* Total energy               */
   add_average(tot_ke/(0.5*kB*system->d_of_f), t_n, 0);
   for(i = 0; i < 3; i++)			/* Non-zero (upper triangle)  */
   {
      add_average(system->h[i][0], h0_n, i);
      add_average(system->h[i][1], h1_n, i);
      add_average(system->h[i][2], h2_n, i);
   }

   zero_real(ke_dyad[0],9);
   for (spec = species; spec < &species[system->nspecies]; spec++)
      energy_dyad(ke_dyad, system->h, spec->vel, spec->mass, spec->nmols);

   k = 0;
   for(i = 0; i < 3; i++)
   {
      for(j = i; j < 3; j++)
      {
         stress[j][i] = stress[i][j] = CONV_P * 0.5/vol * 
          (ke_dyad[i][j] + ke_dyad[j][i] + stress_vir[i][j] + stress_vir[j][i]);
      }
      add_average(stress[i][0], stress0_n, i);
      add_average(stress[i][1], stress1_n, i);
      add_average(stress[i][2], stress2_n, i);
   }
   add_average((stress[0][0] + stress[1][1] + stress[2][2])/3.0, press_n, 0);
   add_average((stress_vir[0][0] + stress_vir[1][1] + stress_vir[2][2])
               *CONV_V/3.0, vir_n, 0);
   
   k = 0;
   for(ispec = 0; ispec < system->nspecies; ispec++)
      for(i = 0; i < 3; i++)
      {
         add_average(CONV_F*CONV_F*meansq_f_t[ispec][0][i], msqf_n, k);
         add_average(CONV_N*CONV_N*meansq_f_t[ispec][1][i], msqt_n, k++);
      }

   for(i = 0; i < 3; i++)
      add_average(CONV_D * dipole[i], dip_n, i);
   /*
    * Update counters.
    */
   if(control.istep >= control.begin_average)
      (av_head->nav)++;
   if(av_head->nroll < control.roll_interval)
      (av_head->nroll)++;
   av_head->iroll = (av_head->iroll+1) % control.roll_interval;
}
/******************************************************************************
 *  value,  roll_av,  roll_sd.   Functions returning the value, rolling       *
 *  average, or s.d for rolling average indicated by pointer p.               *
 ******************************************************************************/
double value(type, comp)
av_n	type;
int	comp;
{
   return(av_info[(int)type].p[comp]->value);
}

double	roll_av(type, comp)
av_n	type;
int	comp;
{
   int	i;
   double	mean = 0.0;

   for(i = 0; i < av_head->nroll; i++)
      mean += av_info[(int)type].p[comp]->roll[i];
   return(mean/ av_head->nroll);
}

double	roll_sd(type, comp)
av_n	type;
int	comp;
{
   int	i;
   double	*roll, ssq = 0.0, mean = roll_av(type, comp), var,
		bottom = -32.0*sqrt((double)control.roll_interval)*precision();

   roll = av_info[(int)type].p[comp]->roll;
   for(i = 0; i < av_head->nroll; i++)
      ssq += roll[i] * roll[i];

   var = ssq/ av_head->nroll - mean*mean;
   if( var * av_head->nroll < ssq * bottom)
      message(NULLI, NULLP, WARNING, NEGVAR, "roll_sd", var,
	      av_info[(int)type].name);

   return(var > 0.0 ? sqrt(var): 0.0);
}
/******************************************************************************
 *  Print frame.   Print the values from the structs pointed at by            *
 *  av_info[i].p in a reasonable format.  Function parameter allows various   *
 *  info from struct to be printed in the same format.  av_info contains      *
 *  the field width and format to use for each data type.                     *
 ******************************************************************************/
void print_frame(header_sym, header_text, f)
char	header_sym;
char	*header_text;
double	(*f)();
{
   int	row, col, icol;
   static int	out_width = 1;
   static boolean	initial = true;
   if(initial)
   {
      initial = false;
      out_line = control.page_length;     	/* Force header on 1st call   */
      for(icol = 0; icol < max_col; icol++)	/* Count total width	      */
         out_width  += av_info[icol].field_width + 1;
      /* if(out_width > control.page_width)*/
   }
   if(out_line + max_row + 1 > control.page_length)	/* If near end of page*/
   {
      new_page();
      for(icol = 0; icol < max_col; icol++)             /* Print column titles*/
         (void)printf(" %*s", av_info[icol].field_width, av_info[icol].name);
      new_line();
   }
   col = 0;					/* Print row of 'header_sym'  */
   while(col++ < 8)				/* with 'header_text' in the  */
      (void)putchar(header_sym);		/* middle.		      */
   col += printf(" %s ", header_text);
   while(col++ < out_width)
      (void)putchar(header_sym);
   new_line();
   
   for(row = 0; row < max_row; row++)		/* Print 'max_col' fields     */
   {						/* across page, up to max_row */
      for(icol = 0; icol < max_col; icol++)	/* down.  Print value returned*/
      {						/* by (*f) in field or fill   */
         (void)putchar(' ');			/* with spaces                */
         if(row < av_info[icol].mult)
            (void)printf( av_info[icol].format, (*f)((av_n)icol, row));
         else
            for(col = 0; col < av_info[icol].field_width; col++)
               (void)putchar(' ');
      }
      new_line();
   }
}
/******************************************************************************
 *  output    Main output routine to be called periodically.                  *
 *  Calls print_frame which does most of the real work.                       *
 ******************************************************************************/
void	output()
{
   char	s[64];

   (void)sprintf(s, "Timestep %d      Current values", control.istep);
   print_frame('=', s, value);

   (void)sprintf(s,"Rolling averages over last %d timesteps", av_head->nroll);
   print_frame('-', s, roll_av);
   print_frame('-', "Standard deviations", roll_sd);
   (void)fflush(stdout);
}
/******************************************************************************
 *  averages   calculate and print averages, reset counter.		      *
 ******************************************************************************/
void	averages()
{
   int	i, iav, col;
   double	variance,
		bottom = -32.0*sqrt((double)av_head->nav)*precision();
   av_t	*av_p;

   if(av_head == NULL)
      message(NULLI, NULLP, FATAL, AVNOC, "averages");

   if(av_head->nav == 0)
   {
      note("no sums accumulated for averages");
      return;
   }   
     
   if(out_line + NAVT +4 > control.page_length)	/* If near end of page*/
      new_page();
   else
      new_lins(2);
   (void)printf( "   Averages over last %d timesteps",av_head->nav);
   new_line();

   av_p = av;
   for(iav = 0; iav < NAVT; iav++)
   {
      for(i = 0; i < av_info[iav].mult; i++)
      {
	 av_p = av_info[iav].p[i];
	 av_p->mean = av_p->sum / av_head->nav;
	 variance = av_p->sum_sq/ av_head->nav - av_p->mean * av_p->mean;
	 if(variance * av_head->nav  < av_p->sum_sq * bottom)
	    message(NULLI, NULLP, WARNING, NEGVAR, "averages", variance,
		    av_info[(int)iav].name);
	 av_p->sd   = variance > 0.0 ? sqrt(variance) : 0.0;
	 av_p->sum = 0.0;  av_p->sum_sq = 0.0;
      }
   }
   av_head->nav = 0;

   for(iav = 0; iav < NAVT; iav++)
   {
      col = 0;
      col += printf( "   %-16s= ",av_info[iav].name);
      for(i = 0; i < av_info[iav].mult; i++)
      {
         if(col + 2 * av_info[iav].field_width + 6 > control.page_width)
         {
            new_line();
            col = printf("                     ");
         }
         col += printf( av_info[iav].format, av_info[iav].p[i]->mean);
         col += printf(" +/- ");
         col += printf(av_info[iav].format, av_info[iav].p[i]->sd);
         if(i < av_info[iav].mult - 1)
            col += printf(",");
      }
      if(col + 10 > control.page_width) new_line();
      col += printf("  %s", av_info[iav].unit);
      new_line();
   }
   (void)fflush(stdout);
}
