#define NELEM 150                       /* Max no of records in element data file   */
#define NSPEC 50                        /* Maximum no of different species          */
#define NLEN 32                         /* Maximum length of species' name          */
#define LLEN 132                        /* Max length of single line in file        */
#define MAX_ATOMS 9999			/* Maximum no of atoms to be read in        */
#define PATHLN 50                       /* Max length of filenames (including path) */

#define ELEPATH  ""
#define POTPATH  ""

#define NOELEM  "Couldn't open element data file \"%s\" for reading"
#define NOPOTL  "Couldn't open potential parameter file \"%s\" for reading"
#define NOSUB   "Species \"%s\" cannot be found for substitution"

#define  PDB              0
#define  CSSR             1
#define  SHAK             2

#define TITLE_SIZE  80

typedef struct              /* Element/species data */
{
   char         name[NLEN], /* Name of species */
                symbol[4];  /* Symbol of species */
   double       mass,       /* Mass of species */
                charge;     /* Charge of species */
} spec_data;

