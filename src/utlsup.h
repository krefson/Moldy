/*
 * utlsup definitions
 */
#define MAX_SPECIES 50                     /* Maximum no of different species          */
#define DOTPROD(x,y)   ((x[0]*y[0])+(x[1]*y[1])+(x[2]*y[2])) /* Scalar product between x and y */
#define xtoupper(c) (islower(c) ? toupper(c) : c)
/*
 * file formats
 */
#define SHAK   0
#define XYZ 1
#define OUTBIN 2
#define DCD 3
#define PDB 4
#define CSSR 5
#define ARC 6
#define XTL 7
/*
 * utlsup functions
 */
extern void init_rdf (system_mt *sys); 
extern gptr *rdf_ptr (int *size); 
extern void new_lins (int ); 
extern int lines_left (void); 
extern void new_page (void); 
extern void new_line (void); 
extern void banner_page (system_mt *sys, spec_mt *spec, restrt_mt *rh); 
extern void note (char *format, ...); 
extern void message (int *nerrs, ...); 
extern void error (char *format, ...); 
extern char *mystrdup (char *s); 
extern int get_int (char *prompt, int lo, int hi); 
extern real get_real (char *prompt, real lo, real hi); 
extern int get_sym (char *prompt, char *cset); 
extern char *get_str (char *prompt); 
extern int forstr (char *instr, int *start, int *finish, int *inc); 
extern void dump_to_moldy (float *buf, system_mt *system); 
extern void invert (real (*a)[3], real (*b)[3]);
extern void mat_vec_mul (real (*m)[3], vec_mp in_vec, vec_mp out_vec, int number);
extern void traj_con (system_mt *system, vec_mt (*prev_cofm), int n); 
extern void traj_con2 (spec_mt *species, vec_mt (*prev_cofm), vec_mt (*traj_cofm), int nspecies); 
extern int range_in (system_mt *system, real (*range)[3]); 
extern char *comm;
extern FILE  *open_dump(char *fname, char *mode);
extern int close_dump(FILE *dumpf);
extern int rewind_dump(FILE *dumpf, int xdr);
extern int read_dump_header(char *fname, FILE *dumpf, dump_mt *hdr_p, boolean *xdr_write,
		     size_mt sysinfo_size, dump_sysinfo_mt *dump_sysinfo);
int     tokenise(char *fields, char *mask, int len);
/*
 * Moldy functions
 */
gptr	*arralloc(size_mt,int,...); 	/* Array allocator */
gptr	*talloc(int n, size_mt size, int line, char *file);
void    tfree(gptr *p);
char    *atime(void);

void    moldy_out(int n, int irec, int inc, system_mt *system, 
		  mat_mp h, spec_mt *species, site_mt *site_info, 
		  int outsw, int intyp, char *insert);
char	*strlower(char *s);
void	read_sysdef(FILE *file, system_mp system, spec_mp *spec_pp, 
		    site_mp *site_info, pot_mp *pot_ptr);
void	initialise_sysdef(system_mp system, spec_mt *species, 
			  site_mt *site_info, quat_mt (*qpf));
void	re_re_header(FILE *restart, restrt_mt *header, contr_mt *contr);
void	re_re_sysdef(FILE *restart, char *vsn, system_mp system, 
		     spec_mp *spec_ptr, site_mp *site_info, pot_mp *pot_ptr);
void	allocate_dynamics(system_mp system, spec_mt *species);
void	lattice_start(FILE *file, system_mp system, spec_mp species, 
		      quat_mt (*qpf));
void	read_restart(FILE *restart, char *vsn, system_mp system, 
		     int av_convert);
void	init_averages(int nspecies, char *vsn, long int roll_interval, 
		      long int old_roll_interval, int *av_convert);
int	getopt(int, char *const *, const char *);
gptr	*talloc(int n, size_mt size, int line, char *file);
void	zero_real(real *r, int n);
void    conv_potentials(const unit_mt *unit_from, const 
			unit_mt *unit_to, pot_mt *potpar, int npotpar,
			int ptype, site_mt *site_info, int max_id);
extern char    *get_line(char *line, int len, FILE *file, int skip);
int	get_tokens(char *buf, char **linev, char *sep);
