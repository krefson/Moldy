/******************************************************************************
 * startup	Functions to perform operations to set up and initialise the  *
 *		simulation and control the reading of control parameters,     *
 *		system specification and restart files.  Contents:	      *
 * start_up()		Master start-up control function		      *
 * allocate_dynamics()	Set up dynamic variable arrays by memory allocation   *
 * initialise_sysdef()	Calculate 'whole-system' quantities and more	      *
 * check_sysdef()	Check new system specification is consistent with old *
 * interp()		Scale derivatives when changing timestep	      *
 * skew_start()		Set up skew-cyclic initial configuration	      *
 * random_start()	Set up random initial configuration (not used)	      *
 *   			(N.B) lattice_start in "input.c"		      *
 * thermalise()		Set up Maxwell-Boltzmann velocity distribution	      *
 * select_spec()	Randomly choose a species			      *
 * random_quat()	Generate quaternion  for uniform random rotation      *
 * gauss_rand()		Return random sample from univariant gaussian         *
 ******************************************************************************
 *      Revision Log
 *       $Log:	startup.c,v $
 * Revision 1.4  89/06/23  15:35:10  keith
 * print-control option deleted.
 * 
 * Revision 1.3  89/06/22  15:45:19  keith
 * Tidied up loops over species to use one pointer as counter.
 * 
 * Revision 1.2  89/06/14  17:57:48  keith
 * Control.out eliminated, use printf and freopen instead to direct output.
 * 
 * Revision 1.1  89/04/27  15:16:41  keith
 * Initial revision
 * 
 * 
 */
#ifndef lint
static char *RCSid = "$Header: startup.c,v 1.5 89/07/05 18:41:16 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#include	<stdio.h>
#include	<math.h>
#include 	"string.h"
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== Library declarations ============================*/
void		 cfree();
/*========================== External function declarations ==================*/
void		message();
void		read_control();
void		read_sysdef();
void		lattice_start();
void		conv_control();
void		conv_potentials();
void		init_averages();
void		init_rdf();
void		re_re_header();
void		re_re_sysdef();
void		read_restart();
void		banner_page();
void		zero_real();
void		jacobi();
void		transpose();
void		mat_vec_mul();
void		q_mul();
void		rot_to_q();
void		print_sysdef();
char		*atime();
double		mdrand();
void		smdrand();
/*========================== External data references ========================*/
extern	unit_t		input_unit;
extern	contr_t		control;
extern	restrt_t	restart_header;
/*========================== Global variables ================================*/
static	unit_t	prog_unit = {MUNIT, LUNIT, TUNIT, QUNIT};
#ifdef	DEBUG
static	char	afmt[] = "    %8s = %8X %8s = %8X %8s = %8X %8s = %8X\
 %8s = %8X %8s = %8X\n";
#endif
/*============================================================================*/
/******************************************************************************
 * gauss_rand. Return a random variable from a gaussian distribution with uni *
 * variance.  After Press, Plannery, Teulkolsk & Vetterling p202.             *
 ******************************************************************************/
double	gauss_rand()
{
   static int		set = 1;
   static double	gset;
   double		v1, v2, r, fac;
   if(set)
   {
      do
      {
         v1 = 2.0*mdrand()-1.0;  v2 = 2.0*mdrand() - 1.0;
         r = v1*v1 + v2*v2;
      } while(r > 1.0);
      fac = sqrt(-2.0 * log(r) / r);
      gset = v1 * fac;
      set = 0;
      return(v2*fac);
   }
   else
   {
      set = 1;
      return(gset);
   }
}
/******************************************************************************
 * random_quat.   Generate quaternions to represent uniform random rotations. *
 * If q=(cos(psi/2), l sin(psi/2)) , l is a unit vector, uniformity requires  *
 * that q(0) be distributed linearly ie p(q(0)) = q(0). L is a random unit    *
 * vector, generated from spherical co-ordinates theta and phi with           *
 * distribution p(theta) = sin(theta) 					      *
 ******************************************************************************/
void	random_quat(q, n)
quat_p	q;				/* First quaternion		(out) */
int	n;				/* Number to be generated.       (in) */
{
   double	phi, cos_theta, sin_theta, st2;
   while(n-- > 0)
   {
      phi = 2.0*PI*mdrand();		/* Phi is uniform on [0, 2pi)	      */
      cos_theta = 1.0 - 2.0*mdrand();	/* 0 <= theta < pi, p(theta)=sin()    */
      sin_theta = sqrt(1.0 - SQR(cos_theta));
      (*q)[0] = sqrt(mdrand());
      st2 = sqrt(1.0 - SQR((*q)[0]));
      (*q)[1] = st2*sin_theta*sin(phi);
      (*q)[2] = st2*sin_theta*cos(phi);
      (*q)[3] = st2*cos_theta;
      q++;
   }
}
/******************************************************************************
 *  select_spec  Choose a species at random, weighted by proportion of whole  *
 ******************************************************************************/
spec_p	select_spec(system, species)
system_p	system;
spec_t		species[];
{
   int		sel = mdrand() * system->nmols;
   spec_p	spec;

   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
      if(sel < spec->nmols)
         return(spec);
      sel -= spec->nmols;
   }
   return((spec_p)-1);
}
/******************************************************************************
 *  skew_start.  Make a starting configuration with all molecules arranged at *
 *  regular intervals on a line passing at an angle through the MD cell.  This*
 *  permits a reasonable spacing between molecules without restricting the    *
 *  number allowed as in a lattice start.                                     *
 ******************************************************************************/
void	skew_start(system, species)
system_p	system;
spec_t	species[];
{
   int		ispec, imol;		/* Counters for species, molecules etc*/
   spec_p	spec;
   double	mass = 0.0;		/* Whole system mass		      */
   double	n_third = exp(log((double)system->nmols)/3.0);
   int		nz = 1, ny = (int)(n_third+0.5), nx = (int)(SQR(n_third)+0.5);
   int		*nmols = ialloc(system->nspecies), nm;
   double	delta_x = (double)nx / system->nmols,
                delta_y = (double)ny / system->nmols,
                delta_z = (double)nz / system->nmols;

   for (spec = species; spec < &species[system->nspecies]; spec++)
      mass += spec->mass * spec->nmols;
      
   system->h[0][0] = system->h[1][1] = system->h[2][2] 
                   = exp((log(mass) - log(control.density))/3.0);
 					/* L = cube root of mass/density      */
   for(imol = 0; imol < system->nmols; imol++)
   {
      do
      {
         spec = select_spec(system, species);		/* Choose species     */
	 ispec = spec-species;
      }
      while(nmols[ispec] >= spec->nmols);	/* Repeat if all set  */

      nm = nmols[ispec];
      spec->c_of_m[nm][0] = imol*delta_x;
      spec->c_of_m[nm][1] = imol*delta_y;
      spec->c_of_m[nm][2] = imol*delta_z;
      nmols[ispec]++;
   }

   random_quat(system->quat, system->nmols_r);
   cfree((char*)nmols);
}
/******************************************************************************
 * random_start.  This function generates a completely random starting        *
 * configuration.  Molecules are placed at random locations and orientations  *
 * in md box, whose size is chosed to give the required density.	      *
 ******************************************************************************/
#ifdef NOT_FOR_NOW
void	random_start(system, species)
system_p	system;
spec_t	species[];
{
   int		imol, i;		/* Counters for species, molecules etc*/
   double	mass = 0.0;		/* Whole system mass		      */

   for (spec = species; spec < &species[system->nspecies]; spec++)
      mass += spec->mass * spec->nmols;
      
   system->h[0][0] = system->h[1][1] = system->h[2][2] 
                   = exp((log(mass) - log(control.density))/3.0);
 					/* L = cube root of mass/density      */
   for(imol = 0; imol < system->nmols; imol++)
      for(i = 0; i < 3; i++)		/* Centre of mass co-ords -1 < x < 1  */
         system->c_of_m[imol][i] = mdrand() - 1.0;
   random_quat(system->quat, system->nmols_r);
}
#endif
/******************************************************************************
 *  thermalise  set velocities and quaternion derivatives to values sampled   *
 *  at random from the Boltzmann distribution for the required temperature.   *
 ******************************************************************************/
void	thermalise(system, species)
system_p	system;
spec_t	species[];
{
   int		imol, i;		/* Counters for species, molecules etc*/
   spec_p	spec;			/* Pointer to species[ispec]	      */
   double	omega_sq;		/* |omega|squared / 4		      */
   double	root_ktm, root_kti[3];	/* Gaussian widths of MB distribution */
   double	total_mass = 0;
   vec_t	momentum;	      	/* Whole system momentum	      */

   zero_real(momentum, 3);
   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
      root_ktm = sqrt(kB * control.temp / spec->mass) / system->h[0][0];
      total_mass += spec->mass*spec->nmols;
      for(imol = 0; imol < spec->nmols; imol++)
         for(i = 0; i < 3; i++)		/* Centre of mass co-ords -1 < x < 1  */
	 {
            spec->vel[imol][i]    = root_ktm * gauss_rand();
	    momentum[i] += spec->mass*spec->vel[imol][i];
	 }

      if(spec->rdof > 0)
      {
         for(i = 0; i < 3; i++)
	    if(spec->inertia[i] != 0.0)
               root_kti[i] = 0.5*sqrt(kB * control.temp / spec->inertia[i]);
	    else
	       root_kti[i] = 0.0;

         for(imol = 0; imol < spec->nmols; imol++)
         {
            for(i = 0; i < 3; i++)	/* Centre of mass co-ords -1 < x < 1  */
               spec->qdot[imol][i+1] = root_kti[i] * gauss_rand();
            omega_sq = -SUMSQ2(spec->qdot[imol]);
            for(i = 0; i < 4; i++)
               spec->qddotvo[imol][i] = spec->qddoto[imol][i] =
               spec->qddot[imol][i] = omega_sq*spec->quat[imol][i];
	 }
      }
   }
   q_mul(system->quat, system->qdot, system->qdot, system->nmols_r);
   for(i = 0; i < 3; i++)		/* Subtract any whole system velocity*/
      for(imol = 0; imol < system->nmols; imol++)
         system->vel[imol][i] -= momentum[i] / total_mass;
}
/******************************************************************************
 *  initialise_sysdef    Computes various quanities and completes the set-up  *
 *  of the system once it has been read in.  Quantities computed are :-       *
 *  Total numbers of molecules, sites, degrees of freedom in whole system,    *
 *  mass and inertia tensor for each species.  The molecular site co-         *
 *  ordinates are re-expressed wrt the molecular centre of mass, the inertia  *
 *  tensor is diagonalised and the site co-ordinates rotated to the principal *
 *  frame.								      *
 ******************************************************************************/
void	initialise_sysdef(system, species, site_info, qpf)
system_p	system;
spec_t		species[];
site_t		site_info[];
quat_t		qpf;			/* Quaternion rotation to princ.frame*/
{
   vec_t	c_of_m;			/* Co-ordinates of centre of mass    */
   vec_t	dipole;			/* Molecular dipole moment	     */
   mat_t	inertia,		/* Inertia tensor		     */
   		v;			/* Transformation matrix to prin. fr.*/
   spec_p	spec;			/* Used for looping over species     */
   double	mass;			/* Temporary for site mass	     */
   int		nz;			/* Count of zero moments of inertia  */
   int		i, j, isite, id; 	/* Various loop counters	     */
   int		nrot;			/* Number of jacobi rotations        */
   boolean	flag;			/* Used to test for charges	     */

   system->nsites  = 0;  system->nmols  = 0;
   system->nmols_r = 0;  system->d_of_f = 0;

   for (spec = species; spec < &species[system->nspecies]; spec++)
   {					/* Loop over molecular species       */
      system->nmols  += spec->nmols;
      system->nsites += spec->nmols * spec->nsites;
      system->d_of_f += spec->nmols * 3;

      zero_real(c_of_m,3);		/* Initialise C_of_M for this species*/
      zero_real(dipole,3);		/* And dipole moment		     */
      for(isite=0; isite < spec->nsites; isite++) /* Calculate (sum m*r) and */
      {					/* molecular mass.		     */
         for(i=0; i<3; i++)
            c_of_m[i] += spec->p_f_sites[isite][i] 
                         * site_info[spec->site_id[isite]].mass;
         spec->mass += site_info[spec->site_id[isite]].mass;
      }

      if(spec->mass < 1.0)		/* Lighter than 1 amu ?              */
         message(NULLI,NULLP,FATAL,ZMASS,spec->name,spec->mass);

      for(i=0; i < 3; i++)		/* Finish calculation of c. of mass. */
         c_of_m[i] /= spec->mass;

      for(isite=0; isite < spec->nsites; isite++)
         for(i=0; i < 3; i++)		/* Subtract c_of_m from co-ordinates */
            spec->p_f_sites[isite][i] -= c_of_m[i];

      if(spec->nsites > 1)		/* If molecule is polyatomic         */
      {
         zero_real(inertia[0],9);	/* Initialise inertia tensor	     */
         for(isite=0; isite < spec->nsites; isite++)
         {
            mass = site_info[spec->site_id[isite]].mass;
            for(i=0; i < 3; i++)	/* Calculate inertia tensor	     */
            {
               inertia[i][i] += mass * SUMSQ(spec->p_f_sites[isite]);
               for(j=0; j < 3; j++)
         	  inertia[i][j] -= mass * spec->p_f_sites[isite][i] 
         	                        * spec->p_f_sites[isite][j];
            }
         }
#ifdef	DEBUG
         printf(" *D* Molecule type %d, mass = %g, C of M = (%g,%g,%g)\n",
                spec-species, spec->mass, c_of_m[0], c_of_m[1], c_of_m[2]);
         print_mat(inertia, " *D* Inertia Tensor");
#endif
         jacobi(inertia, 3, spec->inertia, v, &nrot);
         transpose(v,v);			/* Rotation matrix to pr. fr.*/
	 rot_to_q(v, qpf);			/* make equivalent quaternion*/
#ifdef	DEBUG
         print_mat(v," *D* Rotation Mat.");
#endif
         mat_vec_mul(v,spec->p_f_sites, spec->p_f_sites, spec->nsites);

         nz = (spec->inertia[0]==0) + (spec->inertia[1]==0)  /* Count zero   */
                                    + (spec->inertia[2]==0); /* moments.     */
         spec->rdof = 3-nz;			/* Rotational deg. of freedom*/
         if(spec->rdof > 0)			/* Count molecules with      */
            system->nmols_r += spec->nmols;     /* rotational freedom.       */
         system->d_of_f += spec->rdof * spec->nmols;/* Count total d of f    */

	 for(i = 0; i < 3; i++)			/* Calculate molecular dipole*/
	    for(isite = 0; isite < spec->nsites; isite++)
	       dipole[i] += spec->p_f_sites[isite][i] * 
	       		    site_info[spec->site_id[isite]].charge;
	 spec->dipole = sqrt(SUMSQ(dipole));
      }
   }
#ifdef	DEBUG
   printf(" *D* Totals: nsites = %d, nmols = %d, nmols_r = %d, dof = %d\n",
          system->nsites, system->nmols, system->nmols_r, system->d_of_f);
#endif
   
   flag = false;			/* Test to see if any charges present */
   for(id = 1; id < system->max_id; id++)
      flag = site_info[id].charge != 0.0;
   if(!flag) control.alpha = -1.0;	/* Don't call Ewald sum if not	      */
}
/******************************************************************************
 *  allocate_dynamics	 Allocate memory for the dynamic MD variables         *
 ******************************************************************************/
static void	allocate_dynamics(system, species)
system_p	system;
spec_t	species[];
{
   spec_p	spec;			/* Alias for species[ispec]           */
   int		nmol_cum = 0,		/* Cumulative number of molecules     */
   		nmolr_cum = 0;		/* As above excluding point atoms     */

   system->c_of_m = ralloc(system->nmols);
   system->vel    = ralloc(system->nmols);
   system->velp   = ralloc(system->nmols);
   system->acc    = ralloc(system->nmols);
   system->acco   = ralloc(system->nmols);
   system->accvo  = ralloc(system->nmols);
#ifdef	DEBUG
   printf(" *D* System Dynamic variables (all %d x 3 reals)\n",system->nmols);
   printf(afmt,"c_of_m",system->c_of_m,"vel",system->vel,"velp",system->velp,
          "acc",system->acc,"acco",system->acco,"accvo",system->accvo);
#endif
   if(system->nmols_r > 0)
   {
      system->quat   = qalloc(system->nmols_r);
      system->qdot   = qalloc(system->nmols_r);
      system->qdotp  = qalloc(system->nmols_r);
      system->qddot  = qalloc(system->nmols_r);
      system->qddoto = qalloc(system->nmols_r);
      system->qddotvo= qalloc(system->nmols_r);
#ifdef	DEBUG
   printf(" *D* System Dynamic variables (all %d x 4 reals)\n",system->nmols);
   printf(afmt,"quat",system->quat, "qdot",system->qdot, "qdotp",system->qdotp,
       "qddot",system->qddot,"qddoto",system->qddoto,"qddotvo",system->qddotvo);
#endif
   }
   system->h       = ralloc(3);
   system->hdot    = ralloc(3);
   system->hdotp   = ralloc(3);
   system->hddot   = ralloc(3);
   system->hddoto  = ralloc(3);
   system->hddotvo = ralloc(3);
#ifdef	DEBUG
   printf(" *D* System Dynamic variables (all 9 reals)\n");
   printf(afmt,"h",system->h, "hdot",system->hdot, "hdotp",system->hdotp,
       "hddot",system->hddot,"hddoto",system->hddoto,"hddotvo",system->hddotvo);
#endif

   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
      spec->c_of_m = system->c_of_m + nmol_cum;
      spec->vel    = system->vel    + nmol_cum;
      spec->velp   = system->velp   + nmol_cum;
      spec->acc    = system->acc    + nmol_cum;
      spec->acco   = system->acco   + nmol_cum;
      spec->accvo  = system->accvo  + nmol_cum;
#ifdef	DEBUG
      printf(" *D* Species %d Dynamic variables (all %d x 3 reals)\n",
	     spec-species, spec->nmols);
      printf(afmt,"c_of_m",spec->c_of_m,"vel",spec->vel,"velp",spec->velp,
          "acc",spec->acc,"acco",spec->acco,"accvo",spec->accvo);
#endif
      if(spec->rdof > 0)
      {
         spec->quat    = system->quat    + nmolr_cum;
         spec->qdot    = system->qdot    + nmolr_cum;
         spec->qdotp   = system->qdotp   + nmolr_cum;
         spec->qddot   = system->qddot   + nmolr_cum;
         spec->qddoto  = system->qddoto  + nmolr_cum;
         spec->qddotvo = system->qddotvo + nmolr_cum;
         nmolr_cum += spec->nmols;
#ifdef	DEBUG
         printf(" *D* Species %d Dynamic variables (all %d x 4 reals)\n",
		spec-species, spec->nmols);
         printf(afmt,"quat",spec->quat, "qdot",spec->qdot, "qdotp",spec->qdotp,
             "qddot",spec->qddot,"qddoto",spec->qddoto,"qddotvo",spec->qddotvo);
#endif
      }
      nmol_cum += spec->nmols;
   }
}
/******************************************************************************
 *  Interpolate_derivatives & interp.    Interp is a quadratic interpolation  *
 *  routine for scaling derivatives for a new timestep.                       *
 *  Interpolate_derivatives calls it for all the dynamic variables.	      *
 ******************************************************************************/
void	interp(ratio, x, xo, xvo, n)
double	ratio;				/* Between old and new timesteps      */
real	*x, *xo, *xvo;			/* Pointers to 1st dynamic variable   */
int	n;				/* Size of arrays x, xo, xvo          */
{
   double	c1 = 1.0 - 1.5*ratio + 0.5*SQR(ratio),
   		c2 = 2.0*ratio - SQR(ratio),
   		c3 = -0.5*ratio + 0.5*SQR(ratio),
   		c4 = 1.0 - 3.0*ratio + 2.0*SQR(ratio),
   		c5 = 4.0*ratio - 4.0*SQR(ratio),
   		c6 = -1.0*ratio + 2.0*SQR(ratio);
   double	tmp;
   
   while(n-- > 0)
   {
      tmp	= *xo;
      *xo	= c1 * *x + c2 * *xo + c3 * *xvo;
      *xvo	= c4 * *x + c5 * tmp + c6 * *xvo;
      x++; xo++; xvo++;
   }
}
void	interpolate_derivatives(sys, step, step1)
system_p	sys;
double		step, step1;
{
   double	ratio = step1/step;
   message(NULLI, NULLP, INFO, NEWTS, step, step1);
   interp(ratio, sys->acc[0],   sys->acco[0],  sys->accvo[0],   3*sys->nmols);
   interp(ratio, sys->qddot[0], sys->qddoto[0],sys->qddotvo[0], 4*sys->nmols_r);
   interp(ratio, sys->hddot[0], sys->hddoto[0],sys->hddotvo[0], 9);
}
/******************************************************************************
 *  check_sysdef.   		Read in system specification from the restart *
 *  file and check its consistency with the previously defined spec.  This is *
 *  for use when replacing the system spec in the middle of a run.  The system*
 *  spec from restart is discarded.  Only the minimum checks are applied for  *
 *  the dynamic variables stored in the restart file to make sense with the   *
 *  new system spec.  Checked are a) Number of species, b) Number of molecules*
 *  of each species and c) that the molecules have the same rotational degrees*
 *  of freedom.  The number of sites can change subject to c).                *
 ******************************************************************************/
void	check_sysdef(restart, system, species)
FILE		*restart;		/* Restart file pointer		      */
system_p	system;			/* NEW 'system' struct		      */
spec_t	species[];		/* NEW 'species' struct array	      */
{
   system_t	sys_tmp;		/* Local temporaries of system,       */
   spec_p	spec_tmp, spec;		/* species, site_info and potpar      */
   site_p	site_tmp;		/* used when overwriting restart      */
   pot_p	pot_tmp;		/* sysdef.			      */

   re_re_sysdef(restart, &sys_tmp, &spec_tmp, &site_tmp, &pot_tmp);

   if(system->nspecies != sys_tmp.nspecies)
      message(NULLI, NULLP, FATAL, NSPCON, system->nspecies, sys_tmp.nspecies);
   for (spec = species; spec < &species[system->nspecies]; spec++)
   {
      if(spec->nmols != spec_tmp->nmols)
         message(NULLI, NULLP, FATAL, NMLCON, spec_tmp->name,
                 spec->nmols, spec_tmp->nmols);
      if(spec->rdof != spec_tmp->rdof)
         message(NULLI, NULLP, FATAL, NDFCON, spec_tmp->name,
                 spec->rdof, spec_tmp->rdof);
      cfree((char*)spec_tmp->p_f_sites);
      cfree((char*)spec_tmp->site_id);
   }
   cfree((char*)(spec_tmp-system->nspecies));
   cfree((char*)site_tmp);
   cfree((char*)pot_tmp);
}	
/******************************************************************************
 *  startup	This function sets up everything that is needed to start a    *
 *  run.  It controls the reading in of the control, system specification and *
 *  restart files, conversions to program units, calculation of molecular     *
 *  mass & moments of inertia and evaluation of 'whole system' quantities eg  *
 *  number of molecules.  						      *
 ******************************************************************************/
void	start_up(contr_name, system, species, site_info, potpar)
char		*contr_name;		/* Name of control file "" for stdin  */
system_p	system;			/* Pointer to system struct	      */
spec_p		*species;		/* Pointer to species array	      */
site_p		*site_info;		/* Pointer to site_info array	      */
pot_p		*potpar;		/* Pointer to pot'l parameter array   */
{
   FILE		*contr_file,		/* File pointer for control read      */
   		*sysdef,		/* File pointer for sysdef file read  */
   		*backup = NULL,		/* File pointer for backup file read  */
   		*restart = NULL;	/* File pointer for restart file read */
   double	old_step;		/* Timestep read from restart file    */
   int		old_dump_interval;	/* To check if altered on restart     */
   int		old_max_dumps;		/* To check if altered on restart     */
   boolean	flag;			/* Used to test 'fseek'		      */
   long		pos;			/* Where control info starts on input */
   restrt_t	backup_header;		/* To read backup file header into    */
   contr_t	backup_control;		/* Control struct from backup file    */
   quat_t	qpf;			/* Quat of rotation to principal frame*/

   if(contr_name[0] == '\0')		/* Null name - read control from      */
      contr_file = stdin;		/* standard input.		      */
   else					/* Open named file for reading control*/
   {
      contr_file = fopen(contr_name,"r");
      if(contr_file == NULL)
         message(NULLI, NULLP, FATAL, OCFAIL, contr_name);
   }
   pos = ftell(contr_file);		/* Current file pos needed for cray   */
   read_control(contr_file);		/* Do keyword read of control params  */
   conv_control(&prog_unit,true);	/* Convert to program units           */

   if(control.restart_file[0] != '\0')	/* Open restart file, get backup name */
   {
      if((restart = fopen(control.restart_file,"rb")) == NULL)
         message(NULLI, NULLP, FATAL, ORFAIL, control.restart_file);
      re_re_header(restart, &restart_header, &control);
   }
   if( control.backup_file[0] != '\0' &&
      ( backup = fopen(control.backup_file,"rb")) != NULL )  /* Backup exists */
   {
      re_re_header(backup, &backup_header, &backup_control);
      if(restart && (backup_header.timestamp < restart_header.timestamp))
         backup = NULL;			/* Backup older than restart-don't use*/
   }

   if( backup )
   /* We are restarting from backup file. File is already open and header and *
    * control struct have been read into backup_header and backup_control resp*
    * 1) Set number of steps to number yet to do.			      *
    * 2) Copy backup header and control into the correct structs.	      *
    * 3) Read rest of backup file as for ordinary restart, omitting reread of *
    *    control and scaling since nothing can have changed.		      */
   {
      (void)memcpy((char*)&control, (char*)&backup_control, sizeof(contr_t));
      (void)memcpy((char*)&restart_header,
		   (char*)&backup_header, sizeof(restrt_t));
      re_re_sysdef(backup, system, species, site_info, potpar);
      allocate_dynamics(system, *species);/* Memory for dynamic variables     */
      init_averages(system->nspecies);	/* Create data structure for averages */
      if(control.rdf_interval > 0)
         init_rdf(system);		/* Prepare to calculate rdf	      */
      read_restart(backup, system);	/* Saved dynamic vars and averages    */
      (void)fclose(backup);
      message(NULLI, NULLP, INFO, BACKUP, control.backup_file);
   }
   else if( !restart )			
   /* Initiate a new run from scratch, ie not restart			      *
    * 1) Read the system specification.  This is either in a file of its own  *
    *    or follows the control info.                                         *
    * 2) Convert potential parameters, site masses and charges to prog. units.*
    * 3) Calculate molecular masses, moments of inertia and nmols etc.        *
    * 4) Allocate memory for the MD dynamic variables.			      *
    * 5) Call the dynamic var initialisation routine.			      *
    * 6) Set up averages data structures.				      */
   {
      if( control.sysdef[0] == '\0' || strcmp(control.sysdef,contr_name) == 0 )
         sysdef = contr_file;		/* Sys def'n is tacked onto control   */
      else				/* Sys def'n is in separate file      */
      {					/* Open system specification file     */
         sysdef = fopen(control.sysdef,"r");
         if(sysdef == NULL)
            message(NULLI, NULLP, FATAL, ODFAIL, control.sysdef);
      }
					/* Read system specification file     */
      read_sysdef(sysdef, system, species, site_info, potpar);

      conv_potentials(&input_unit, &prog_unit, *potpar, system->n_potpar,
                         system->ptype, *site_info, system->max_id);
      initialise_sysdef(system, *species, *site_info, qpf);
      allocate_dynamics(system, *species);	/* Allocate dynamic arrays    */

      smdrand(control.seed);			/* Seed random number generato*/
      if( control.lattice_start )		/* Choose startup method      */
	 lattice_start(sysdef, system, *species, qpf); /* Lattice - from file */
      else	
         skew_start(system, *species);		/* Start from skew lattice    */
      thermalise(system, *species);		/* Initialise velocities      */

      init_averages(system->nspecies);		/* Prepare to collect averages*/
      if( control.rdf_interval > 0 )
         init_rdf(system);			/* Prepare to calculate rdf   */

      (void)strcpy(restart_header.init_date, atime());
      (void)strcpy(restart_header.title,control.title);

      if(sysdef != contr_file)
         (void)fclose(sysdef);			/* Finished with sys spec file*/
   }
   else
   /* Continue from a restart file.  The restart file contains a header, the  *
    * saved 'control' struct, system specification, dynamic variables and     *
    * averages.  The old 'control' values act as defaults for this run. To do *
    * this they are read into 'control', and the control file is rewound and  *
    * reread.  To keep the units consistent, the old control is converted to  *
    * the NEW input units before rereading, and then back to program units.   *
    * Also the scale_end and begin_average parameters are incremented by the  *
    * initial timestep only if they are explicitly specified.                 *
    * The next thing to be read is the system spec.  Normally this is just    *
    * taken from the restart file, but a facility is provided to allow its    *
    * replacement with a setup from a new system spec 'source' file.  Its     *
    * consistency with the old sys-spec is only checked in so far as the size *
    * and shape of the dynamic variable arrays are identical. Beware.         *
    * Memory is allocated for the dynamic variables and averages which are    *
    * read from the input file.  Finally if the timestep has been altered the *
    * accelerations at previous timesteps are adjusted to the new step.       */
   {
      old_step = control.step;				/* Needed for scaling */
      old_dump_interval	     = control.dump_interval;	/* Check if these par */
      old_max_dumps	     = control.maxdumps;	/* -amaters altered.  */
      conv_control(&prog_unit, false);
      control.scale_end     -= control.istep;		/* These parameters   */
      control.begin_average -= control.istep;		/* are respecified    */
      control.begin_rdf     -= control.istep;		/* RELATIVE to current*/
      control.begin_dump    -= control.istep;		/* time step.	      */
      control.nsteps	    -= control.istep;

      flag = fseek(contr_file,pos,0);			/* Rewind control file*/
      if( flag )
         message(NULLI,NULLP,FATAL,SEFAIL,contr_name[0]?contr_name:"stdin");
      read_control(contr_file);				/* Reread control file*/
      conv_control(&prog_unit,true);			/* Back to prog units */

      control.scale_end     += control.istep;
      control.begin_average += control.istep;
      control.begin_rdf     += control.istep;
      control.begin_dump    += control.istep;
      control.nsteps	    += control.istep;
      if( control.maxdumps != old_max_dumps ||		/* Need to restart    */
	  control.dump_interval != old_dump_interval )	/* dump seq if changed*/
      {
         control.begin_dump = control.istep + 1;	/* Set new beginning  */
      }

      if( !control.new_sysdef )		/* Usual case, get sysdef from restart*/
         re_re_sysdef(restart, system, species, site_info, potpar);
      else				/* Get sysdef from new sys-spec file  */
      {
         if( control.sysdef[0]=='\0' || strcmp(control.sysdef,contr_name) == 0 )
            sysdef = contr_file;	/* Sys def'n is tacked onto control   */
         else				/* Sys def'n is in separate file      */
         {				/* Open system specification file     */
            sysdef = fopen(control.sysdef,"r");
            if( sysdef == NULL )
               message(NULLI, NULLP, FATAL, OSFAIL, control.sysdef);
         }
	 /* Read in and set up new system spec (just as in case of new run)   */
         read_sysdef(sysdef, system, species, site_info, potpar);
         conv_potentials(&input_unit, &prog_unit, *potpar, system->n_potpar,
                            system->ptype, *site_info, system->max_id);
#ifdef	DEBUG
         printf(" *D* Read and converted new system specification\n");
#endif
         initialise_sysdef(system, *species, *site_info, qpf);
         check_sysdef(restart, system, *species);/* Consistent with saved one?*/
         if(sysdef != contr_file)
            (void)fclose(sysdef);
         control.reset_averages = 0;	/* Averages invalid if sysdef changed */
      }
      allocate_dynamics(system, *species);/* Memory for dynamic variables     */
      init_averages(system->nspecies);	/* Create data structure for averages */
      if(control.rdf_interval > 0)
         init_rdf(system);		/* Prepare to calculate rdf	      */
      read_restart(restart, system);	/* Saved dynamic vars and averages    */
      if(control.step != old_step)
         interpolate_derivatives(system, old_step, control.step);
      (void)fclose(restart);
      message(NULLI, NULLP, INFO, RESUCC, control.restart_file);
   }               
   (void)fclose(contr_file);

   if( control.out_file[0] != NULL )	/* Open output file (or use stdout)   */
   {
      (void)fflush(stdout);		/* Purge buffer before opening file   */
      if( freopen(control.out_file, "a", stdout) == NULL )
         message(NULLI, NULLP, FATAL, OOFAIL, control.out_file);
   }
   banner_page(system, *species);
}
