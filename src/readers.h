/* reader functions */
char	*read_ftype(char *filename);
int     read_ele(spec_data *element, char *filename);
int     read_pdb(char *, mat_mp, char (*)[NLEN], vec_mp, double *, char *, char *);
int     read_cssr(char *, mat_mp, char (*)[NLEN], vec_mp, double *, char *, char *);
int     read_shak(char *, mat_mp, char (*)[NLEN], vec_mp, double *, char *, double *);
int     read_xtl(char *, mat_mp, char (*)[NLEN], vec_mp, double *, char *, char *);
int     read_xyz(char *, mat_mp, char (*)[NLEN], vec_mp, char *);
int     read_shelx(char *, mat_mp, char (*)[NLEN], vec_mp, char *, char *);
int     read_pot(char *, pot_mp *, site_mt *, int, int);
int     is_symbol(const char *);

/* sginfo functions */
int     sgexpand(int , int , vec_mt *, char (*)[NLEN], double *, char *);

/* Moldy functions */
double  det(mat_mt a);
