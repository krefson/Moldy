/* Variable sizes can be adjusted by user where necessary */
#define NELEM 150                       /* Max no of records in element data file   */
#define NLEN 32                         /* Maximum length of species' name          */
#define LLEN 132                        /* Max length of single line in file        */
#define MAX_ATOMS 9999			/* Maximum no of atoms to be read in        */
#define PATHLN 50                       /* Max length of filenames (including path) */

/* Paths should end in a slash so that filenames can be concatenated */
#define ELEPATH  "./"
#define POTPATH  "./"

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
