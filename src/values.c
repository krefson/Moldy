/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1988, 1992, 1993 Keith Refson
 
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
what you give them.   Help stamp out software-hoarding!  */
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
 */
/*========================== Program include files ===========================*/
#include	"defs.h"
#include	"structs.h"
#include	"messages.h"
/*========================== Library include files ===========================*/
#include	<math.h>
#include	<string.h>
#include	<stdio.h>
/*========================== External function declarations ==================*/
void	mat_vec_mul(real (*m)[3], vec_mp in_vec, vec_mp out_vec, int number);
void	q_conj_mul(quat_mp p, quat_mp q, quat_mp r, int n);
					/* Quaternion multiply conjugated     */
double	det(real (*a)[3]);		/* Determinant of 3x3 matrix	      */
double	vdot(int n, real *x, int ix, real *y, int iy); /* Vector dot product  */
double	sum(register int n, register double *x, register int ix); /* Vect. sum*/
void	zero_real(real *r, int n);
void	zero_double(double *r, int n);
void	zero_dbls(double *r, size_mt n);
void	energy_dyad(mat_mt ke_dyad, mat_mt h, real s, 
		    vec_mp vels, double mass, int nmols);
double  ke_cell(mat_mt hmom, real w);
double	trans_ke(real (*h)[3], vec_mt (*vel_s), real s, double mass, int nmols);
double	rot_ke(quat_mt (*omega_p), real, real *inertia, int nmols);
double	precision(void);		/* Machine precision constant	      */
char	*atime(void);
void	new_line(void);
void	new_lins(int n);
void	new_page(void);
gptr    *talloc(int n, size_mt size, int line, char *file);
				       /* Interface to memory allocator       */
void    tfree(gptr *p);		       /* Free allocated memory	      	      */
int	lines_left(void);		/* on current page of output	      */
void		note(char *, ...);	/* Write a message to the output file */
void		message(int *, ...);	/* Write a warning or error message   */
/*========================== External data references ========================*/
extern	contr_mt	control;            /* Main simulation control parms. */
/*========================== Structs local to module =========================*/

typedef struct
{   double	value,
   		sum,
   		sum_sq,
		mean,
		sd,
   		roll[1];
} av_mt;

typedef struct
{
   av_n		id;
   char		name[12], *unit;
   int		field_width;
   char		*format;
   int		mult;
   av_mt		**p;
} av_info_t;
/*========================== Global variables ================================*/
static
av_info_t av_info[] = { {tke_n, "Trans KE",   CONV_E_N,	11, "%11.5g",-1, 0},
			{rke_n, "Rot KE",     CONV_E_N,	11, "%11.5g",-1, 0},
			{pe_n,  "Pot Energy", CONV_E_N,	11, "%11.5g",NPE,0},
			{e_n,   "Tot Energy", CONV_E_N,	11, "%11.5g", 2, 0},
			{tt_n,  "TTemp",      CONV_T_N,	 6, "%6.1f", -1, 0},
			{rt_n,  "RTemp",      CONV_T_N,	 6, "%6.1f", -1, 0},
			{t_n,   "Temp",	      CONV_T_N,	 6, "%6.1f", 1,  0},
			{h0_n,  "h(1,*)",     LUNIT_N,	 6, "%6.2f", 3,  0},
			{h1_n,  "h(2,*)",     LUNIT_N,	 6, "%6.2f", 3,  0},
			{h2_n,  "h(3,*)",     LUNIT_N,	 6, "%6.2f", 3,  0},
			{stress0_n,"Stress",  CONV_P_N,	10, "%10.3g",3,  0},
			{stress1_n,"Stress",  CONV_P_N,	10, "%10.3g",3,  0},
			{stress2_n,"Stress",  CONV_P_N,	10, "%10.3g",3,  0},
			{press_n,"Pressure",  CONV_P_N,	10, "%10.3g",1,  0},
			{vir_n, "Virial",     CONV_E_N,	11, "%11.5g",1,  0},
			{msqf_n,"<F**2>",     CONV_F_N,	10, "%10.5g",-3, 0},
			{msqt_n,"<N**2>",     CONV_N_N,	10, "%10.5g",-3, 0},
			{dip_n, "Dip Mom",    CONV_D_N,	 8, "%8.2g", 3,  0}};
static  size_mt	av_size;		/* Size of averages database          */
static  size_mt	av_mt_size;		/* Size of entry inaverages database  */
static  av_head_mt  *av_head;
static	av_mt	*av;			/* Dynamic array of averages structs  */
static	int	navs = 0;		/* Size of array av                   */
static  int	e_n_fill = 0;
static	int	max_row = 0;		/* Largest number of components       */
static  int	max_col = (int)press_n;	/* Number to print across page        */
static  size_mt av_tmp_size;
static  gptr    *av_tmp;
/*========================== Macros ==========================================*/
#define NAVT			(int)end
#define INC(av_mp)    (av_mp = (av_mt*)((double*)av_mp + av_mt_size/sizeof(double)))
#define DEC(av_mp)    (av_mp = (av_mt*)((double*)av_mp - av_mt_size/sizeof(double)))
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
void	init_averages(int nspecies, char *vsn, long int roll_interval, 
		      long int old_roll_interval, int *av_convert)
{
   av_mt		*av_mp;
   int		i, imult;
   int		major, minor;

   for(i = 0; i < NAVT; i++)	/* cycle over enum types av_n         */
   {
      if(av_info[i].mult < 0)		/* Set true multiplicity of each type */
         av_info[i].mult = -nspecies * av_info[i].mult;
      if( i == e_n )
	 e_n_fill = navs;
      navs += av_info[i].mult;		/* Count total			      */
      if(i < max_col && av_info[i].mult > max_row)
         max_row = av_info[i].mult;
   }
   if( control.const_pressure > 0 )
   {
      if( control.const_temp > 0 )
	 strcpy(av_info[3].name,"Energy E,G");
      else
	 strcpy(av_info[3].name,"Energy E,H");
   }
   else
   {
      if( control.const_temp > 0 )
	 strcpy(av_info[3].name,"Energy E,F");
      else
	 strcpy(av_info[3].name,"Energy E");
   }

   /*
    * Do we have to do any conversion on averages read from restart file?
    * We just allocate buffers and set flags here.
    */
   *av_convert = 0;
   if( vsn )
   {
      /*
       * First check whether restart was written by 2.19 or earlier.
       */
      if( sscanf(vsn, "%d.%d", &major, &minor) < 2 )
	 message(NULLI, NULLP, FATAL, INRVSN, vsn);
      if( major < 2 || (major==2 && minor <= 19 ) )
      {
	 av_tmp_size =  sizeof(av_head_mt) 
	             + (navs-1)*(sizeof(av_mt)+(old_roll_interval-1)*sizeof(double));
	 av_tmp = balloc(1,av_tmp_size);
	 *av_convert = 1;
      }
      /*
       * Has size of rolling average store changed?
       */
      if (roll_interval != old_roll_interval )
      {
	 if( *av_convert == 1 )
	 {
	    message(NULLI, NULLP, WARNING, FIXROL);
	    roll_interval = control.roll_interval = old_roll_interval;
	 }
	 else
	 {
	    av_tmp_size =  sizeof(av_head_mt) 
	       + navs*(sizeof(av_mt)+(old_roll_interval-1)*sizeof(double));
	    av_tmp = balloc(1, av_tmp_size);
	    *av_convert = 2;
	 }
      }
   }
   /* Determine size of database, Allocate space and set pointers             */
   av_mt_size = sizeof(av_mt)+(roll_interval-1)*sizeof(double);
   av_size = sizeof(av_head_mt) + navs*av_mt_size;
   av_head  = (av_head_mt*)balloc(1, av_size);
   av_head->nav = av_head->nroll = av_head->iroll = 0;
   av       = (av_mt*)(av_head+1);
   zero_dbls((double*)av, navs*av_mt_size/sizeof(double));

   av_mp = av;
   for(i = 0; i < (int)end; i++)	/* Set up pointers to area of array   */
   {					/* reserved for each type, size=mult. */
      av_info[i].p = aalloc(av_info[i].mult, av_mt*);
      for(imult = 0; imult < av_info[i].mult; imult++)
      {
	 av_info[i].p[imult] = av_mp;
	 INC(av_mp);
      }
   } 
}
/******************************************************************************
 * convert_averages.  Update averages database if roll_interval changed or    *
 * if restart file written using old "static" scheme.			      *
 * Also a convenient place to implement reset_averages.			      *
 ******************************************************************************/
void	convert_averages(long roll_interval, long old_roll_interval, 
			 int av_convert)
{
   int iav, old_nroll, old_iroll, rbl;
   size_mt prev_av_mt_size;
   av_head_mt	*prev_av_head = (av_head_mt *)av_tmp;
   av_mt	        *av_mp, *prev_av_mp;

   if( av_convert == 1 )
   {
      av_mp = av_info[0].p[0];
      prev_av_mp = (av_mt *)(prev_av_head+1);

      for(iav = 0; iav < e_n_fill; iav++)
      {
	 memcp(av_mp, prev_av_mp, av_mt_size);
	 INC(av_mp);
	 INC(prev_av_mp);
      }
      memcp(av_mp, prev_av_mp, av_mt_size);
      INC(av_mp);
      for(iav++; iav < navs; iav++)
      {
	 memcp(av_mp, prev_av_mp, av_mt_size);
	 INC(av_mp);
	 INC(prev_av_mp);
      }
      control.reset_averages++;
   }
   
   if( av_convert == 2 )
   {
      prev_av_mt_size = sizeof(av_mt)+(old_roll_interval-1)*sizeof(double);
      old_nroll = prev_av_head->nroll;
      old_iroll = prev_av_head->iroll;
      av_head->nroll = MIN(old_nroll,roll_interval);
      av_head->iroll = av_head->nroll % roll_interval;
      av_head->nav   = prev_av_head->nav;
      rbl = MIN(old_iroll, av_head->nroll);
      prev_av_mp = (av_mt *)(prev_av_head+1);
      av_mp = av_info[0].p[0];
      for(iav = 0; iav < navs; iav++)
      {
	 /*
	  * Can do a struct copy -- will only pick up 1st roll entry
	  */
	 *av_mp = *prev_av_mp;
	 memcp(av_mp->roll+av_head->nroll-rbl, prev_av_mp->roll+old_iroll-rbl,
	       rbl*sizeof(double));
	 memcp(av_mp->roll, prev_av_mp->roll+old_nroll-av_head->nroll+rbl,
	       (av_head->nroll-rbl)*sizeof(double));

	 INC(av_mp);
	 prev_av_mp = (av_mt*)((double*)prev_av_mp + prev_av_mt_size/sizeof(double));
      }
   } 
   /*
    *  Reset averages and counters to zero if a) requested
    *  or b) we have not yet reached begin_average. (The latter 
    *  avoids the situation of non-contiguous averages).
    *  Yes I know we just set them up, but this way is clearest.
    */
   if( control.reset_averages || control.istep+1 <= control.begin_average )
   {
      av_mp = av_info[0].p[0];
      for(iav = 0; iav < navs; iav++)
      {
	 av_mp->sum = av_mp->sum_sq = av_mp->mean = av_mp->sd = 0.0;
	 INC(av_mp);
      }
      av_head->nav = 0;
   }
}
/******************************************************************************
 * av_ptr   Return a pointer to averages database and its size (for restart)  *
 ******************************************************************************/
gptr	*av_ptr(size_mt *size, int av_convert)
{
   if(av_convert)
   {
      *size = av_tmp_size;
      if(av_tmp != 0)
	 return(av_tmp);
   }
   else
   {
      *size = av_size;
      if(av_head != 0)
	 return((gptr*)av_head);
   }
   message(NULLI, NULLP, FATAL, AVNOC, "av_ptr");
   return 0;						/* To satisfy lint    */
}
/******************************************************************************
 * add_average  update the averages database with new datum                   *
 ******************************************************************************/
static void	add_average(double datum,/* Datum to store and accumulate sums*/ 
			    av_n type, 	 /* What kind (ie where to store)     */
			    int offset)	 /* Sub-type or which component       */
{
   av_mt		*av_mp;
   if(offset < 0 || offset > av_info[(int)type].mult - 1)
      message(NULLI, NULLP, FATAL, AVBNDS, offset, av_info[(int)type].name);

   av_mp = av_info[(int)type].p[offset];
 
   av_mp->value = datum;
   if(control.istep >= control.begin_average)
   {
      av_mp->sum += datum;
      av_mp->sum_sq += datum * datum;
   }
   av_mp->roll[av_head->iroll] = datum;
}
/******************************************************************************
 * values   Calculate the values of the thermodynamic quantities, maintain and*
 * if necessary print them, their averages and rolling averages.              *
 ******************************************************************************/
void	values(system_mt *system,        /* record of system info             */
	       spec_mt   *species, 	 /* Records of info for each species  */
	       vec_mt   (*meansq_f_t)[2],/* mean square forces and torques    */ 
	       double    *pe, 		 /* potential energy real,reciprocal  */
	       real      *dipole, 	 /* dipole moment of whole system     */
	       real     (*stress_vir)[3])/* 'Potential' part of stress virial */
{
   spec_mp	spec;
   int		ispec, ipe;
   double	e, tot_ke = 0.0, tot_pe = 0.0;
   double	ske, gktls, hke, hpe;
   int		i, j, k;
   int		tdof = 0;
   mat_mt	ke_dyad,
                stress;
   double	vol = det(system->h);
   static	double *tkep1, *tkem1, *tkem3, *tkem5;
   static	double *rkep1, *rkem1, *rkem3, *rkem5;
   static	double skep1, skem1,skem3, skem5;
   static	double hkep1, hkem1,hkem3, hkem5;
   static	boolean firstcall = true;
   static	double tsold;

   if( firstcall )
   {
      firstcall = false;
      tkep1 = (double*)aalloc(system->nspecies, double);
      tkem1 = (double*)aalloc(system->nspecies, double);
      tkem3 = (double*)aalloc(system->nspecies, double);
      tkem5 = (double*)aalloc(system->nspecies, double);
      rkep1 = (double*)aalloc(system->nspecies, double);
      rkem1 = (double*)aalloc(system->nspecies, double);
      rkem3 = (double*)aalloc(system->nspecies, double);
      rkem5 = (double*)aalloc(system->nspecies, double);
      for(ispec=0; ispec < system->nspecies; ispec++)
      {
	 tkep1[ispec] = rkep1[ispec] = 
	 tkem1[ispec] = rkem1[ispec] = 
	 tkem3[ispec] = rkem3[ispec] = 
	 tkem5[ispec] = rkem5[ispec] = -1.0;
      }
      skep1 = skem1 = skem3 = skem5 = -1.0;
      hkep1 = hkem1 = hkem3 = hkem5 = -1.0;
      tsold = system->ts;
   }

    for (spec = species; spec < &species[system->nspecies]; spec++)
      if( ! spec->framework )
        tdof += 3*spec->nmols;
   /*
    * Mazur's leapfrog KE predictor. [J. Comp. Phys (1997) 136, 354-365) Eq. 30]
    */
#ifndef KE1
#define KEINT(a,b,c,d) (1.0/48.0*(15*a+45*b-15*c+3*d))
#else
#define KEINT(a,b,c,d) (a)
#endif

   for(ipe = 0; ipe < NPE; ipe++)
   {
      add_average(CONV_E * pe[ipe], pe_n, ipe);
      tot_pe += pe[ipe];
   }

   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
      ispec = spec-species;
      tkem5[ispec] = tkem3[ispec]; 
      tkem3[ispec] = tkem1[ispec]; 
      tkem1[ispec] = tkep1[ispec];
      tkep1[ispec] = trans_ke(system->h, spec->mom, 0.5*(tsold+system->ts),  spec->mass, spec->nmols);
      if(tkem1[ispec] < 0.0)
	 tkem1[ispec] = tkem3[ispec] = tkem5[ispec] = tkep1[ispec];
      e = KEINT(tkep1[ispec], tkem1[ispec], tkem3[ispec], tkem5[ispec]);
      add_average(CONV_E * e, tke_n, ispec);	/* c of m  kinetic energy     */
      tot_ke += e;

      add_average(e/(1.5*spec->nmols*kB)*tdof/(tdof-3), tt_n, ispec);
						/* c of mass temp.	      */
      if(spec->rdof > 0)			/* Only if polyatomic species */
      {
	 rkem5[ispec] = rkem3[ispec]; 
	 rkem3[ispec] = rkem1[ispec]; 
	 rkem1[ispec] = rkep1[ispec];
         rkep1[ispec] = rot_ke(spec->amom, 0.5*(tsold+system->ts), spec->inertia, spec->nmols);
	 if(rkem1[ispec] < 0.0)
	    rkem1[ispec] = rkem3[ispec] = rkem5[ispec] = rkep1[ispec];
	 e = KEINT(rkep1[ispec], rkem1[ispec], rkem3[ispec], rkem5[ispec]);
         add_average(CONV_E * e, rke_n, ispec);	/* Rotational kinetic energy  */
         tot_ke += e;
         add_average(e/(0.5*kB*spec->rdof*spec->nmols), rt_n, ispec);
      }						/* and temperature            */
   }   						/* Overall temperature        */

   tsold = system->ts;
   skem5 = skem3;
   skem3 = skem1;
   skem1 = skep1;
   skep1 = SQR(system->tsmom)/(2.0*control.ttmass);
   if( skem1 < 0.0 )
      skem1 = skem3 = skem5 = skep1;
   ske = KEINT(skep1, skem1, skem3, skem5);

   gktls = system->d_of_f*kB*control.temp*log(system->ts);

   hkem5 = hkem3;
   hkem3 = hkem1;
   hkem1 = hkep1;
   hkep1 = ke_cell(system->hmom, control.pmass);
   if( hkem1 < 0.0 )
      hkem1 = hkem3 = hkem5 = hkep1;
   hke = KEINT(hkep1, hkem1, hkem3, hkem5);

   hpe = control.pressure*vol;

   add_average(CONV_E*(tot_ke+tot_pe+ske+gktls + hke + hpe), e_n, 0); 
   add_average(CONV_E*(tot_ke+tot_pe), e_n, 1); 

   add_average(tot_ke/(0.5*kB*system->d_of_f), t_n, 0);
   for(i = 0; i < 3; i++)			/* Non-zero (upper triangle)  */
   {
      add_average(system->h[i][0], h0_n, i);
      add_average(system->h[i][1], h1_n, i);
      add_average(system->h[i][2], h2_n, i);
   }

   zero_real(ke_dyad[0],9);
   for (spec = species; spec < &species[system->nspecies]; spec++)
      energy_dyad(ke_dyad, system->h, system->ts, spec->mom, spec->mass, spec->nmols);

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
         add_average(CONV_N*CONV_N*meansq_f_t[ispec][1][i], msqt_n, k);
	 k++;
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
double value(av_n type, int comp)
{
   return(av_info[(int)type].p[comp]->value);
}

double	roll_av(av_n type, int comp)
{
   int	i;
   double	mean = 0.0;

   for(i = 0; i < av_head->nroll; i++)
      mean += av_info[(int)type].p[comp]->roll[i];
   return(mean/ av_head->nroll);
}

static
double	roll_sd(av_n type, int comp)
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
static
void print_frame(int header_sym, char *header_text, double (*f)(av_n, int))
{
   int	row, col, icol;
   static int	out_width = 1;
   static boolean	initial = true;
   if(initial)
   {
      for(icol = 0; icol < max_col; icol++)	/* Count total width	      */
         out_width  += av_info[icol].field_width + 1;
   }
   if( initial || lines_left() < max_row + 1 )	/* If near end of page*/
   {
      new_page();
      for(icol = 0; icol < max_col; icol++)             /* Print column titles*/
         (void)printf(" %*s", av_info[icol].field_width, av_info[icol].name);
      new_line();
      initial = false;
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
void	output(void)
{
   char	s[64];

   (void)sprintf(s, "Timestep %ld      Current values", control.istep);
   print_frame('=', s, value);

   (void)sprintf(s,"Rolling averages over last %d timesteps", av_head->nroll);
   print_frame('-', s, roll_av);
   print_frame('-', "Standard deviations", roll_sd);
   (void)fflush(stdout);
}
/******************************************************************************
 *  averages   calculate and print averages, reset counter.		      *
 ******************************************************************************/
void	averages(void)
{
   int	i, iav, col;
   double	variance,
		bottom = -32.0*sqrt((double)av_head->nav)*precision();
   av_mt	*av_mp;

   if(av_head == 0)
      message(NULLI, NULLP, FATAL, AVNOC, "averages");

   if(av_head->nav == 0)
   {
      note("no sums accumulated for averages");
      return;
   }   
     
   if(lines_left() < NAVT + 4 )	/* If near end of page*/
      new_page();
   else
      new_lins(2);
   (void)printf( "   Averages over last %d timesteps",av_head->nav);
   new_line();

   for(iav = 0; iav < NAVT; iav++)
   {
      for(i = 0; i < av_info[iav].mult; i++)
      {
	 av_mp = av_info[iav].p[i];
	 av_mp->mean = av_mp->sum / av_head->nav;
	 variance = av_mp->sum_sq/ av_head->nav - av_mp->mean * av_mp->mean;
	 if(variance * av_head->nav  < av_mp->sum_sq * bottom)
	    message(NULLI, NULLP, WARNING, NEGVAR, "averages", variance,
		    av_info[(int)iav].name);
	 av_mp->sd   = variance > 0.0 ? sqrt(variance) : 0.0;
	 av_mp->sum = 0.0;  av_mp->sum_sq = 0.0;
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
