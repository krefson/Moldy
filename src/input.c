/******************************************************************************
 * Input	Functions for reading and verifying the input files (except   *
 *		the restart file). Contents:				      *
 * Strlower();		Convert a string to lowercase (Internal use only)     *
 * Get_line();		Read next input line.         (Internal use only)     *
 * Read_sysdef()       	Read the system specification file     	       	      *
 * Print_sysdef()       Inverse of read_sysdef()       	       	       	      *
 * Lattice_start()	Read initial crystal structure and set it up	      *
 * Read_control()       Read control file				      *
 ******************************************************************************
 *      Revision Log
 *       $Log:	input.c,v $
 * Revision 1.1  89/04/20  16:00:42  keith
 * Initial revision
 * 
 */
#ifndef lint
static char *RCSid = "$Header: input.c,v 1.1 89/04/20 16:00:42 keith Exp $";
#endif
/*========================== Library include files ===========================*/
#include	<ctype.h>
#include	<math.h>
#include 	"string.h"
/*========================== Program include files ===========================*/
#include	"structs.h"
#include	"messages.h"
/*========================== External function declarations ==================*/
void		q_mul_1();
void		message();
/*========================== External data references ========================*/
extern	contr_t	control;		/* Main simulation control record     */
extern	unit_t	input_unit;		/* Unit specification (see Convert.c) */
/*========================== Structs local to module =========================*/
typedef	struct				/* Struct template for keyword	      */
{					/* in read_control.		      */
   char	*key,
	*format,
	*ptr;
}	match_t;
/*========================== Potential type specification ====================*/
static char	*types[] = {"lennard-jones","buckingham","mcy"};
static int	npotp[]  = {2,		    3,	 	 4};
/*========================== Macros ==========================================*/
#define	NPOTT	(sizeof npotp / sizeof(int))
#define		LLEN		132
		/* Flags to indicate status of potpar and site_info records   */
#define		S_USED		0x01
#define		S_MASS		0x02
#define		S_CHARGE	0x04
#define		S_NAME		0x08
		/* Copy a struct of type pot_t			      */
#define	potcpy(s1, s2)	(void)memcpy((char*)(s1),(char*)(s2), sizeof(pot_t))
/*========================== Control file keyword template ===================*/
					/* format SFORM is defined as %NAMLENs*/
					/* in structs.h, to avoid overflow */
static	match_t	match[] = {
		{"title",	    	SFORM,	 control.title},
                {"nsteps",	    	"%d",	(char*)&control.nsteps},
                {"step",	    	"%lf",	(char*)&control.step},
                {"print-sys-spec",  	"%d",	(char*)&control.print_sysdef},
                {"new-sys-spec",    	"%d",	(char*)&control.new_sysdef},
		{"lattice-start",	"%d",	(char*)&control.lattice_start},
                {"sys-spec-file",   	SFORM,	 control.sysdef},
                {"restart-file",    	SFORM,	 control.restart_file},
                {"save-file",	   	SFORM,	 control.save_file},
                {"dump-file",	    	SFORM,	 control.dump_file},
                {"backup-file",	    	SFORM,	 control.backup_file},
		{"temp-file",		SFORM,	 control.temp_file},
                {"out-file",	    	SFORM,	 control.out_file},
		{"nbins",		"%d",	(char*)&control.nbins},
		{"seed",		"%d",	(char*)&control.seed},
                {"page-width",	    	"%d",	(char*)&control.page_width},
                {"page-length",	   	"%d",	(char*)&control.page_length},
                {"scale-interval",    	"%d",	(char*)&control.scale_interval},
		{"scale-separately",	"%d", (char*)&control.scale_separately},
                {"const-pressure",  	"%d",	(char*)&control.const_pressure},
                {"reset-averages",  	"%d",	(char*)&control.reset_averages},
                {"scale-end",	   	"%d",	(char*)&control.scale_end},
                {"begin-average",   	"%d",	(char*)&control.begin_average},
                {"average-interval",	"%d", (char*)&control.average_interval},
		{"begin-dump",		"%d",   (char*)&control.begin_dump},
		{"dump-interval",  	"%d",	(char*)&control.dump_interval},
		{"dump-level",		"%d",	(char*)&control.dump_level},
		{"ndumps",		"%d",	(char*)&control.maxdumps},
		{"backup-interval",	"%d",  (char*)&control.backup_interval},
                {"roll-interval",   	"%d",	(char*)&control.roll_interval},
                {"print-interval",  	"%d",	(char*)&control.print_interval},
		{"begin-rdf",		"%d",	(char*)&control.begin_rdf},
		{"rdf-interval",	"%d",	(char*)&control.rdf_interval},
		{"rdf-out",		"%d",	(char*)&control.rdf_out},
                {"temperature",	    	"%lf",	(char*)&control.temp},
                {"pressure",	    	"%lf",	(char*)&control.pressure},
                {"w",		    	"%lf",	(char*)&control.pmass},
                {"cutoff",	    	"%lf",	(char*)&control.cutoff},
		{"subcell",		"%lf",	(char*)&control.subcell},
                {"density",	    	"%lf",	(char*)&control.density},
                {"alpha",	    	"%lf",	(char*)&control.alpha},
		{"k-cutoff",		"%lf",	(char*)&control.k_cutoff},
		{"rdf-limit",		"%lf",	(char*)&control.limit},
		{"cpu-limit",		"%lf",	(char*)&control.cpu_limit},
                {"mass-unit",	    	"%lf",	(char*)&input_unit.m},
                {"length-unit",	    	"%lf",	(char*)&input_unit.l},
                {"time-unit",		"%lf",	(char*)&input_unit.t},
                {"charge-unit",		"%lf",	(char*)&input_unit.q}
                          };
#define		NMATCH	(sizeof match / sizeof(match_t))
/*=============================================================================
 |   Start of functions							      |
 =============================================================================*/
/******************************************************************************
 *  get_line  read an input line skipping blank and comment lines	      *
 ******************************************************************************/
char	*get_line(line, len, file)
char	*line;
int	len;
FILE	*file;
{
   char	*s, *t;
   do
   {
      s = fgets(line, len, file);		/* Read one line of input     */
      if(s == NULL) break;			/* exit if end of file        */
      t = s + strlen(s) - 1;
      while(t >= s && (*t == ' ' || *t == '\t' || *t == '\n'))
         *t-- = '\0';				/* Strip trailing white space */
   }
   while(*s == '\0' || *s == '#');		/* Repeat if blank or comment */
   if(s == NULL)
      *line = '\0';				/* Return null at eof         */
   return(line);
}
/******************************************************************************
 * strlower   convert a string to lowercase anr return a pointer to it        *
 ******************************************************************************/
char	*strlower(s)
char	*s;
{
   char	*t;
   for(t = s; *t != '\0'; t++)
      *t = isupper(*t) ? tolower(*t) : *t;
   return(s);
}

/******************************************************************************
 *  read_sysdef    Read the system specification file which must be open and  *
 *  pointed to by parameter 'file'.  Set up the structures system and species *
 *  and arrays site_info and potpar (allocating space) and read in and check  *
 *  the supplied values.  The reading is done in two passes.  Pass 1 simply   *
 *  counts the number of species, number of sites on each species and the     *
 *  largest site identifier index in order to allocate the dynamic arrays.    *
 *  Pass 2 does the actual reading and checking.                              *
 ******************************************************************************/
void	read_sysdef(file, system, spec_pp, site_info, pot_ptr)
FILE		*file;			/* File pointer to read info from     */
system_p	system;			/* Pointer to system array (in main)  */
spec_p		*spec_pp;		/* Pointer to be set to species array */
site_p		*site_info;		/* To be pointed at site_info array   */
pot_p		*pot_ptr;		/* To be pointed at potpar array      */
{
   int		nspecies = 0,		/* Number of distinct species         */
   		max_id = 0,		/* Largest site identifier index      */
   		id, idi, idj,		/* Temp. site identifier index        */
   		ispec, isite,		/* species and site counters	      */
		sflag,			/* Temporary flag		      */
		i,			/* Counter			      */
		n_potpar,		/* Number of parameters for this pot'l*/
   		n_items;		/* How many items scanf found in input*/
   struct list_t {int n; struct list_t *p;};/* Template for linked list nsites*/
   struct list_t nsites_base,		/* Head of list (contains no datum)   */
   		*nsites,		/* List entry for current species     */
                *last = &nsites_base;	/* List entry for previous species    */
   int		nerrs = 0;		/* Accumulated error count	      */
   int		flag;			/* Used to test 'fseek' result        */
   long		start_pos = ftell(file);/* Rewind marker for second pass      */
   char		name[LLEN],		/* Species name temporary             */
   		line[LLEN];		/* Store for input line from file     */
   double	mass, charge, p_tmp;	/* Local temporaries		      */
   double	p_f_sites[3];		/* Local temporary		      */
   pot_p	pp1;			/* Used for acces to potpar ij and ji */
   spec_p	species, spec;		/* Local pointer to species array     */
   site_p	s_ptr;			/* Local pointer to site info array   */
   static pot_t	pot = {S_USED};	/* Local storage for potentials       */

   message(&nerrs,NULLP,INFO,SYSRD);
   /* First pass - read system definition and count nspecies, nsites, max_id  */
   (void)get_line(line,LLEN,file);		/* Read first line.	      */
   while(sscanf(line, "%s", name) > 0 && strcmp(strlower(name), "end") != 0)
   {						/* Loop, parsing 'line' for   */
      nspecies++;				/* name of new species.       */
      nsites = aalloc(1, struct list_t); 	/* Make new list element      */
      last->p = nsites;				/* Link it in		      */
      last = nsites;				/* Backwards pointer for link */
      while(sscanf(get_line(line,LLEN,file), "%d", &id) > 0)
      {						/* Loop, reading and parsing  */
         nsites->n++;				/* for integer ie new site id.*/
         max_id = MAX(max_id, id);		/* Count nsites, greatest id. */
      }						/* Leave 'line' if parse fails*/
   }
   if(nspecies == 0)				/* Empty file??		      */
      message(&nerrs,NULLP,FATAL,NOSPEC);

   /* Allocate arrays of species and site info records */
   max_id++;
   system->max_id = max_id;
   *spec_pp    = aalloc(nspecies, spec_t );
   *site_info  = aalloc(max_id, site_t );
   *pot_ptr    = aalloc(max_id * max_id, pot_t );
   species = *spec_pp;   			/* Local pointer for neatness.*/
   system->nspecies = nspecies;

   flag = fseek(file, start_pos, 0);		/* Prepare to reread input.   */
   if(flag)
      message(NULLI, NULLP, FATAL, SEFAIL, "control file");
   nsites = &nsites_base;
   /* Pass 2.  read system definition and set up species and site_info arrays */
   for(ispec = 0, spec = species; ispec < nspecies; ispec++, spec++)
   {						/* Loop over all species.     */
      n_items = sscanf(get_line(line,LLEN,file),"%s %d", name, &spec->nmols);
      name[sizeof spec->name-1] = '\0';		/* Truncate before copying    */
      (void)strcpy(spec->name, name);		/* to avoid overflow.         */
      nsites = nsites->p;			/* Find next element of list  */
      spec->nsites = nsites->n;			/* which contains nsites.     */
      if(n_items < 2)				/* Number of mols not supplied*/
         message(&nerrs,line, ERROR, NONUM, name);
      if(spec->nmols <= 0)			/* Can't have <=0 molecules   */
         message(&nerrs,line,ERROR, NOMOLS, spec->nmols, name);
      if(spec->nsites <=0)			/* or ghost molecules!        */
         message(&nerrs,NULLP,ERROR,NOSITE,spec->nsites,name);
      spec->p_f_sites = ralloc(spec->nsites);	/* Allocate space and set     */
      spec->site_id   = ialloc(spec->nsites+1);	/* pointers for each species. */

      for(isite = 0; isite < spec->nsites; isite++)
      {						/* Loop over sites on molecule*/
        n_items =sscanf(get_line(line,LLEN,file), "%d %lf %lf %lf %lf %lf %s",
                        &id,			/* Get and parse line of input*/
                        p_f_sites,
                        p_f_sites + 1,
                        p_f_sites + 2,
                        &mass,  &charge,  name);
        if(id <= 0)				/* Test for valid site index. */
        {
           message(&nerrs,line, ERROR, INVSID, id);
           id = 0;
        }
        spec->site_id[isite] = id;		/* Put id into rightful place.*/
        name[sizeof s_ptr->name - 1] = '\0';	/* Truncate site name.        */
        s_ptr = *site_info + id;		/* Reference (*site_info)[id].*/
        s_ptr->flag |= S_USED;
        switch (n_items)			/* Handle input items in      */
        {					/* reverse order.             */
           case 7:				/* Site name supplied.	      */
              if(s_ptr->flag & S_NAME && strcmp(name, s_ptr->name) != 0)
                 message(&nerrs,line,ERROR,NCONF,id,s_ptr->name);
              else
                 (void)strcpy(s_ptr->name, name);
              s_ptr->flag |= S_NAME;
           case 6:				/* Site charge supplied.      */
              if(s_ptr->flag & S_CHARGE && charge != s_ptr->charge)
                 message(&nerrs,line,ERROR,CCONF,id, s_ptr->charge);
              else
                 s_ptr->charge = charge;
              s_ptr->flag |= S_CHARGE;
           case 5:				/* Site mass supplied.        */
              if(s_ptr->flag & S_MASS && mass != s_ptr->mass)
                 message(&nerrs,line,ERROR,MCONF,id, s_ptr->mass);
              else if(mass < 0.0)
                 message(&nerrs,NULLP,ERROR,INVMAS,id,mass);
              else
                 s_ptr->mass = mass;
              s_ptr->flag |= S_MASS;
           case 4:				/* All site co-ordinates      */
	      for( i = 0; i < 3; i++ )
	         spec->p_f_sites[isite][i] = p_f_sites[i];
              break;				/* Normal exit from switch    */
           case 3:                      	/* One or more site 	      */
           case 2:                      	/* co-ordinates are missing.  */
           case 1:
              message(&nerrs,line,ERROR,MISSCO,n_items-1,id);
              break;
           default:				/* Should never occur.	      */
              message(&nerrs,NULLP, FATAL, BROKEN, n_items);
        }
      }
   }
   /* Check that all sites have been fully defined, and for gaps in ordering. */
   for(id = 1; id < max_id; id++)
   {
      sflag = (*site_info)[id].flag;
      if(sflag & S_USED)
      {
         if(! (sflag & S_MASS))
            message(&nerrs,NULLP,ERROR,NOMASS,id);
         if(! (sflag & S_CHARGE))
            message(&nerrs,NULLP,ERROR,NOCGRG,id);
         if(! (sflag & S_NAME))
            message(&nerrs,NULLP,ERROR,NONAME,id);
      }
      else
         message(&nerrs,NULLP, WARNING, NOTUSD,id);
   }
   (void)get_line(line,LLEN,file);		/* read "end" -left by pass 1 */

   /* Next line is keyword indicating type of potentials to be used	      */
   n_items = sscanf(get_line(line,LLEN,file), "%s", name);
   for(i = 0; i < NPOTT; i++)			/* Is 'name' a known type?    */
      if(strcmp(strlower(name), types[i]) == 0)
         break;
   if(i == NPOTT)				/* Did the loop find 'name'?  */
      message(&nerrs,line,FATAL,UNKPOT,name);	/* no			      */
   system->ptype = i;				/* yes		              */
   n_potpar = system->n_potpar = npotp[i];
   						/* Now read in parameters     */
   while(sscanf(get_line(line,LLEN,file),"%s",name) > 0
                    && strcmp(strlower(name), "end") != 0)
   {
      n_items = 0;
      if(sscanf(line,"%d %d %[^#]",&idi,&idj,line) <= 2) /* Not enough values */
         message(&nerrs,line,ERROR,NOPAIR);
      else
      {						/* Parse potential parameters */
	 (void)strcat(line, "$");		/* Add marker to end	      */
         while( n_items < NPOTP && sscanf(line,"%lf %[^#]", &p_tmp, line) > 1 )
	    pot.p[n_items++] = p_tmp;
      }
      if (n_items < n_potpar)
         message(&nerrs,line,ERROR,NOPOTP,n_potpar);
      else				/* Test id pair and if OK store values*/
      {
         if(idi < 1 || idi >= max_id)
            message(&nerrs,line,ERROR,IDOUTR,idi);
         if(idj < 1 || idj >= max_id)
            message(&nerrs,line,ERROR,IDOUTR,idj);
         if(!(  (*site_info)[idi].flag & S_USED
              &&(*site_info)[idj].flag & S_USED))
            message(&nerrs,line,WARNING,EXTPOT);
         pp1 = (*pot_ptr) + (idi + idj * system->max_id);
         if(pp1->flag & S_USED)		/* pot'l for this id pair already set?*/
            message(&nerrs,line,ERROR,DUPPOT);
         else				/* Put values into pp1  and  pp2      */
         {
            potcpy(*pot_ptr + idi + idj * system->max_id, &pot);
            potcpy(*pot_ptr + idj + idi * system->max_id, &pot);
         }
      }
   }
   /* Check whether potentials have been defined for all 'used' site id's     */
   for(idi = 0; idi < max_id; idi++)
      for(idj = idi; idj < max_id; idj++)
         if( (   (*site_info)[idi].flag & S_USED    /* True if sites idi, idj */
              && (*site_info)[idj].flag & S_USED)   /* both used but pot'l not*/
            && !((*pot_ptr)[idi + max_id * idj].flag & S_USED))
            message(&nerrs,NULLP,WARNING,NOPOT,idi,idj);

   if(nerrs > 0)			/* if any errors have been detected   */
      message(&nerrs,NULLP,FATAL,ERRS,nerrs);
   else
      message(&nerrs,NULLP,INFO,SUCCES);
}
/******************************************************************************
 *  print sysdef   Print out the definition of the system, in the format that *
 *  read_sysdef can interpret.						      *
 ******************************************************************************/
void	print_sysdef(system, species, site_info, potpar, file)
system_p	system;			/* Pointer to system array (in main)  */
spec_t		species[];		/* Pointer to species array 	      */
site_p		site_info;		/* pointer to site_info array	      */
pot_t		potpar[];		/* Potential parameter array	      */
FILE		*file;			/* File pointer to write info to      */
{
   spec_p	spec;
   int	ispec, isite, idi, idj, idij, ip;
   int	n_potpar = npotp[system->ptype];
   for(ispec = 0, spec = species; ispec < system->nspecies; ispec++, spec++)
   {
      (void)fprintf(file," %-16s  %d\n", spec->name, spec->nmols);
      for(isite=0; isite < spec->nsites; isite++)
         (void)fprintf(file," %6d %9g %9g %9g %9g %9g %s\n",
         		spec->site_id[isite],
	         	spec->p_f_sites[isite][0],
	         	spec->p_f_sites[isite][1],
	         	spec->p_f_sites[isite][2],
	         	site_info[spec->site_id[isite]].mass,
	         	site_info[spec->site_id[isite]].charge,
	         	site_info[spec->site_id[isite]].name);
   }
   (void)fprintf(file," end\n");
   (void)fprintf(file," %s potential parameters\n",types[system->ptype]);
   for(idi = 1; idi < system->max_id; idi++)
      for(idj = idi; idj < system->max_id; idj++)
      {
         idij = idj + idi * system->max_id;
         if(potpar[idij].flag & S_USED)
         {
            (void)fprintf(file," %6d %6d", idi, idj);
            for(ip = 0; ip < n_potpar; ip++)
               (void)fprintf(file,"%9g",potpar[idij].p[ip]);
            (void)fprintf(file,"\n");
         }
      }
   (void)fprintf(file," end\n");
}
/******************************************************************************
 * lattice_start   Initialse the simulation co-ordinates on a lattice. Read   *
 * from the end of the system specification file.  The format is one line     *
 * specifying the unit cell (6 x floating point + 3 x int # cells)	      *
 *    a  b  c  alpha  beta  gamma  nx ny nz                                   *
 * followed by one line for each molecule in the unit cell:                   *
 *    species  x  y  z  q0  q1  q2  q3                                        *
 * 'species'  is the name,  x, y, z are FRACTIONAL co-ords and 4 quaternions. *
 ******************************************************************************/
void	lattice_start(file, system, species, qpf)
FILE	*file; 				/* File to read info from	      */
system_p system;			/* System info struct		      */
spec_p	species;			/* Array of species info structs      */
quat_t	qpf;				/* Princ frame rotation quaternion    */
{
   typedef struct init_s {int species;  struct init_s *next;
                  double r[3], q[4];} init_t; 	/* For linked list of coords  */
   init_t	*cur, *init = NULL;		/* Current and header of list */
   double	a, b, c, alpha, beta, gamma;	/* Unit cell lengths, angles  */
   int		ix, iy, iz, nx, ny, nz;		/* Number of unit cells in MDC*/
   char		line[LLEN], name[LLEN];
   int		n_items, nerrs = 0, ispec, imol, i;
   int		*nmols = ialloc(system->nspecies);
   real		ca, cb, cg, sg;
   quat_t	q;

   n_items = sscanf(get_line(line,LLEN,file),"%lf%lf%lf%lf%lf%lf%d%d%d",
		    &a, &b, &c, &alpha, &beta, &gamma, &nx, &ny, &nz);
   if(n_items < 9)
      message(&nerrs, line, ERROR, NOCELL);
   if( ! (a > 0 && b > 0 && c > 0 && nx > 0 && ny > 0 && nz > 0 &&
	  alpha > 0 && alpha < 180.0 && beta > 0 && beta < 180.0 &&
	  gamma > 0 && gamma < 180.0))
      message(&nerrs, line,  ERROR, INVCEL);

   ca = cos(alpha*DTOR); cb = cos(beta*DTOR); cg = cos(gamma*DTOR);
   sg = sin(gamma*DTOR);

   system->h[0][0] = nx*a;			/* Set up MD cell matrix      */
   system->h[0][1] = ny*b * cg;			/* from lengths and angles.   */
   system->h[1][1] = ny*b * sg;
   system->h[0][2] = nz*c * cb;
   system->h[1][2] = nz*c / sg * sqrt(ca*ca + cb*cb * cg*cg - 2*ca*cb*cg);
   system->h[2][2] = nz*c / sg * sqrt(1 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg);

   while(sscanf(get_line(line,LLEN,file), "%s", name) > 0 &&
	 strcmp(strlower(name), "end") != 0)	/* Cycle over lines in file   */
   {
      cur =aalloc(1, init_t );			/* Get struct for new molecule*/
      cur->next = init;  init = cur;		/* Link it into list	      */
      n_items = sscanf(line, "%s%lf%lf%lf%lf%lf%lf%lf",
		       name, &cur->r[0], &cur->r[1], &cur->r[2],
		       &cur->q[0], &cur->q[1], &cur->q[2], &cur->q[3]);
      if(n_items > 1)				/* Have name of molecule      */
      {
	 for(ispec = 0; ispec < system->nspecies; ispec++)
	    if(strcmp(strlower(name),strlower(species[ispec].name)) == 0)
	       break;
	 if(ispec >= system->nspecies)		/* Didn't find it	      */
	    message(&nerrs,NULLP,ERROR,UNKSPE,name);
	 else					/* Found it - check values    */
	 {
	    cur->species = ispec;
	    if(n_items < 4)
	       message(&nerrs, line, ERROR, FEWCOO, name);
	    if(n_items < 8 && species[cur->species].rdof != 0)
	       message(&nerrs, line, ERROR, FEWQUA, name);
	    if(cur->r[0] < 0 || cur->r[0] >= 1 ||
	       cur->r[1] < 0 || cur->r[1] >= 1 ||
	       cur->r[2] < 0 || cur->r[2] >= 1)
	       message(&nerrs,NULLP,ERROR,FRACCO,cur->r[0],cur->r[1],cur->r[2]);
	    if(fabs(1.0 - SQR(cur->q[0]) - SQR(cur->q[1])
		        - SQR(cur->q[2]) - SQR(cur->q[3])) > 1e-4)
	       message(&nerrs,NULLP,ERROR,QNORM,
		       cur->q[0],cur->q[1],cur->q[2],cur->q[3]);
	    nmols[cur->species] += nx*ny*nz;
	 }
      }
   }
   for(ispec = 0; ispec < system->nspecies; ispec++)
   {
      if(nmols[ispec] != species[ispec].nmols)
         message(&nerrs,NULLP,ERROR,NIMOLS,species[ispec].name,
		 nmols[ispec],species[ispec].nmols);
      nmols[ispec] = 0;
   }

   if(nerrs > 0)				/* Is file all correct?       */
      message(NULLI, NULLP, FATAL, INITER, nerrs);

   for(cur = init; cur != NULL; cur = cur->next)
   {
      ispec = cur->species;
      for(ix = 0; ix < nx; ix++)
         for(iy = 0; iy < ny; iy++)
            for(iz = 0; iz < nz; iz++)
	    {
	       imol = nmols[ispec];
	       species[ispec].c_of_m[imol][0] = (cur->r[0]+ix)/nx - 0.5;
	       species[ispec].c_of_m[imol][1] = (cur->r[1]+iy)/ny - 0.5;
	       species[ispec].c_of_m[imol][2] = (cur->r[2]+iz)/nz - 0.5;
	       for( i = 0; i < 4; i++ )
	          q[i] = cur->q[i];		/* Convert type to 'real'     */
	       q_mul_1(q, qpf, species[ispec].quat[imol]);
	       nmols[ispec]++;
	    }
   }
}

/******************************************************************************
 *  read_control.   Read the control parameters from the standard input.      *
 *  Input lines are of the form " key = value ", one per line.  The struct    *
 *  'match' is searched for a matching key, and if found converts the value   *
 *  according to the format string and stores it at the value of the pointer  *
 *  in 'match'.	"name=" with no value means assign a null string. 	      *
 ******************************************************************************/
void	read_control(file)
FILE	*file;
{
   char		line[LLEN],
   		name[LLEN],
   		value[LLEN];
   int		i, n_items;
   int		nerrs = 0;

   while( *get_line(line,LLEN,file) != '\0' )
   {
      n_items = sscanf(line, " %[^= ] = %127[^#]",name, value);
      if(!strcmp(strlower(name),"end"))
         break;
      if( !strcmp(name,"?") )
      {
	 for( i = 0; i < NMATCH; i++ )
	 {
	    (void)printf(" %s =",match[i].key);
	    (void)printf(match[i].format,match[i].ptr);
	    (void)printf("\n");
	 }
	 continue;
      }
      if( n_items < 1 )
         message(&nerrs,line,ERROR,NOVAL,name);
      else
      {
         for( i = 0; i < NMATCH; i++ )		/* Search table for key	      */
            if( !strcmp(strlower(name), match[i].key) )
               break;				/* Found it          	      */
	 if( i == NMATCH )			/* Reached end without success*/
            message(&nerrs,line,ERROR,NOTFND,name);
         else					/* Found it, so convert value */
         {
            if( n_items == 1 && strcmp(match[i].format,SFORM) == 0 )
                *match[i].ptr = '\0';		/* name=<empty> - assign null */
	    else
	    {
               n_items = sscanf(value, match[i].format, match[i].ptr);
               if( n_items < 1 )
                  message(&nerrs,NULLP,ERROR,BADVAL,value,name);
	    }
         }
      }
   }
   if( nerrs > 0 )
      message(&nerrs,NULLP,FATAL,ERRCON,nerrs,(nerrs>1)?'s':' ');
   else
      message(&nerrs,NULLP,INFO,SUCCON);
}
