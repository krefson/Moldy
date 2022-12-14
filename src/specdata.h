/* Variable sizes can be adjusted by user where necessary */
#define NELEM 150                       /* Max no of records in element data file   */
#define NLEN 16                         /* Maximum length of species' name          */
#define LLEN 132                        /* Max length of single line in file */
#define MAX_ATOMS 19999                 /* Maximum no of atoms to be read in */
#define MAX_SPECIES 50                  /* Maximum no of different species */

/* Items below here should not need to be altered */
#define NOELEM  "Couldn't open element data file \"%s\" for reading"
#define NOPOTL  "Couldn't open potential parameter file \"%s\" for reading"
#define NOSUB   "Species \"%s\" cannot be found for substitution"

#define TITLE_SIZE  80

typedef struct              /* Element/species data */
{
   char         name[NLEN], /* Name of species */
                symbol[4];  /* Symbol of species */
   double       mass,       /* Mass of species */
                charge;     /* Charge of species */
} spec_data;

typedef struct              /* Element/species data */
{
   char         *name, /* Name of species */
                *symbol;  /* Symbol of species */
   double       mass,       /* Mass of species */
                charge;     /* Charge of species */
} ele_data;
