#include "defs.h"
#include	"stddef.h"
#include 	"string.h"
#include 	<stdio.h>
#include	"structs.h"

#undef MIN
#define MIN(x,y) ( (x) > (y) ? (y) : (x))

char *malloc(), *calloc();

/*========================== Macros ==========================================*/
#define DUMP_SIZE(level, n, n_r)  \
   			 (( (level & 1) + (level>>1 & 1) + (level>>2 & 1) ) * \
                                    (3*n + 4*n_r + 9)+ \
                             (level>>3 & 1) * \
                                    (3*n + 3*n_r + 9) +\
                             (level & 1))
/*============================================================================*/
typedef struct list_t
{
   struct list_t	*next;
   int			i;
   char *p;
   int num;
} list_t;

typedef struct cpt_t
{
   int ncpt, offset, size, mols;
} cpt_t;


void
insert(entry, head)
list_t	*entry, *head;
{
   while( head->next != NULL && entry->i > head->next->i)
      head = head->next;

   entry->next = head->next;
   head->next  = entry;
}

void
print_list(head)
list_t	*head;
{
   if(head == NULL)
      return;
   fprintf(stderr,"%-8d%s\n",head->i, head->p);
   print_list(head->next);
}

/******************************************************************************
 * forstr.  Parse string str of format s-f:n  (final 2 parts optional),       *
 *          returning integer values of s,f,n. f defaults to s and n to 1     *
 ******************************************************************************/
int
forstr(str, start, finish, inc)
char	*str;
int	*start, *finish, *inc;
{
   char	*p, *pp;
   long strtol();
   
   if( (p = strchr(str,':')) != NULL)
   {
      *inc = strtol(p+1, &pp, 0);
      if( pp == p+1 )
	 goto limerr;
      *p = 0;
   }
   else
      *inc = 1;
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
      *start = *finish = strtol(str, &pp, 0);
      if( pp == str )
	 goto limerr;
   }
   return 0;
 limerr:
   return -1;
}
/******************************************************************************
 * Seq gen  Generate a sequence of file names from a "printf" style spec and  *
 * range of filenames.							      *
 ******************************************************************************/
char **seqgen(fmt, range)
char	*fmt, *range;
{
   int	start, finish, inc;
   int  i, iseq = 0;
   char **seq;
   char	buf[256];

   if( forstr(range, &start, &finish, &inc) != 0 )
      return NULL;

   if( (seq=(char**)calloc(1+(finish-start+1)/inc,sizeof(char*))) == NULL)
      return NULL;

   for(i = start; i <= finish; i++)
   {
      seq[iseq++] = strdup((sprintf(buf,fmt,i),buf));
   }
   seq[iseq] = NULL;
   return seq;
}
/******************************************************************************
 * Put.  Write data in text or binary form.				      *
 ******************************************************************************/
void
put(buf, n, bflg)
float *buf;
int   n;
int   bflg;
{
#ifdef DEBUG2
   fprintf(stderr,"Put: %d at %8x (%s)\n", n, buf, bflg? "binary":"text");
#endif
   if( bflg )
      fwrite((char*)buf, sizeof(float), n, stdout);
   else
      while(n-- > 0)
	 fprintf(stdout,"%7g ",*buf++);

}
/******************************************************************************
 * Extract.  Process one dump file, extracting and outputting data.	      *
 ******************************************************************************/
void
extract(dump_name, cpt_mask, molecules, cpt, ncpt, tslice, num, inc, bflg, nmols)
char	*dump_name;
int	cpt_mask;
list_t	*molecules;
cpt_t	cpt[];
int	ncpt, tslice, num, inc;
int	bflg, nmols;
{
   FILE		*dump_file;
   dump_t	header;
   float	*buf = (float*)calloc(4*nmols,sizeof(float));/* nmols > nmols_r */
   long		dump_base;
   int		icpt;
   list_t	*mol;
   
   if( (dump_file = fopen(dump_name, "rb")) == NULL)
   {
      fprintf(stderr, "Failed to open dump file \"%s\"\n", dump_name);
      exit(2);
   }
#ifdef DEBUG
   fprintf(stderr,"Working on file \"%s\" (%d-%d)\n", dump_name, tslice, num);
#endif
   fread((char*)&header, sizeof(dump_t), 1, dump_file);
   if( ferror(dump_file) || feof(dump_file) )
   {
      fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
      exit(2);
   }

   dump_base = sizeof(dump_t)+tslice*header.dump_size*sizeof(float);
   while(tslice < num && tslice < header.dump_interval * header.ndumps)
   {
#ifdef DEBUG
      fprintf(stderr,"Timeslice %d of %d\n", tslice, num);
#endif
      for( icpt = 0; icpt < ncpt; icpt++)
      {
	 if( cpt_mask & (1 << icpt) )
	 {
	    fseek(dump_file, dump_base+cpt[icpt].offset*sizeof(float), 0);
	    fread((char*)buf, sizeof(float), cpt[icpt].size, dump_file);
	    if( cpt[icpt].mols )
	       for(mol = molecules; mol != 0; mol = mol->next)
		  put(buf+mol->i*cpt[icpt].ncpt, mol->num*cpt[icpt].ncpt, bflg);
	    else
	       put(buf, cpt[icpt].size, bflg);
	 }
      }
      if( ! bflg )
	 putchar('\n');
      tslice += inc;
      dump_base += inc*sizeof(float) * header.dump_size;
   }
   (void)free((char*)buf);
}

main(argc, argv)
int	argc;
char	*argv[];
{
   int	c;
   extern int	optind;
   extern char	*optarg;
   int		errflg = 0, genflg = 0, tsflg = 0, bflg = 0;
   int          nmols, nmols_r;
   int		cpt_mask;
   char		*dump_name;
   char		*tsrange, *filerange;
   char		**filelist;
   FILE		*dump_file;
   int		nfiles = 0;
   int		start,finish,inc;
   int		tslice, numslice, maxslice;
   int		find, offset, icpt;
   
   static cpt_t cpt[] = {{3, 0, 3, 1},
			  {4, 0, 4, 1},
			  {9, 0, 1, 0},
			  {1, 0, 1, 0},
			  {3, 0, 3, 1},
			  {4, 0, 4, 1},
			  {9, 0, 1, 0},
			  {3, 0, 3, 1},
			  {4, 0, 4, 1},
			  {9, 0, 1, 0},
			  {3, 0, 3, 1},
			  {3, 0, 3, 1},
			  {9, 0, 1, 0} };
#define NCPT (sizeof(cpt)/sizeof(cpt_t))

   static int level_mask[16] = {  0x0000,0x000f,0x0070,0x007f,
				  0x0380,0x038f,0x03f0,0x03ff,
				  0x1c00,0x1c0f,0x1c70,0x1c7f,
				  0x1f80,0x1f8f,0x1ff0,0x1fff};

   dump_t	proto_header, header;

   list_t	f_head;
   list_t	mol_head;
   list_t	*cur;

   mol_head.next = NULL;

   while( (c = getopt(argc, argv, "c:bR:Q:n:t:m:") ) != EOF )
      switch(c)
      {
       case 'c':
	 cpt_mask = strtol(optarg,(char*)NULL,0);
	 break;
       case 'b':
	 bflg++;
	 break;
       case 'R':
	 nmols = strtol(optarg,(char*)NULL,0);
	 break;
       case 'Q':
	 nmols_r = strtol(optarg,(char*)NULL,0);
	 break;
       case 'n':
	 filerange = optarg;
	 genflg++;
	 break;
       case 't':
	 if( tsflg++ == 0)
	    tsrange = optarg;
	 else
	    errflg++;
	 break;
       case 'm':
	 if( forstr(optarg, &start, &finish, &inc) )
	    errflg++;
	 else
	 {
	    cur = (list_t *)calloc(1, sizeof(list_t));
	    cur->i = start;  cur->num=finish-start+1;
	    insert(cur, &mol_head);
	 }	 
	 break;
       case '?':
	 errflg++;
      }


   if(optind >= argc)
   {
      fprintf(stderr,"%s: no dump files given\n",argv[0]);
      errflg++;
   }
   if( errflg )
   {
      fprintf(stderr,
	   "Usage: dumpextract [-Rn] [-Qn] [-b] [-n range] [-c cpts]\
 [-t timeslices] [-m molecules] dumpfiles\n");
      exit(2);
   }

   /*
    *  Generate list of dump files if required
    */
   if(genflg)
   {
      if( (filelist = seqgen(argv[optind], filerange)) == NULL)
      {
	 fprintf(stderr,"%s: invalid dump range \"%s\"\n", argv[0]);
	 errflg++;
      }
   }
   else
      filelist = argv+optind;

   /*
    *  Check all dump files for correctness and build ordered list
    */
   for(find = 0; filelist[find] != NULL; find++)
   {
      dump_name = filelist[find];
      if( (dump_file = fopen(dump_name, "rb")) == NULL)
      {
	 fprintf(stderr, "Failed to open dump file \"%s\"\n", dump_name);
	 exit(2);
      }
      fread((char*)&header, sizeof(dump_t), 1, dump_file);
      if( ferror(dump_file) || feof(dump_file) )
      {
	 fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
	   exit(2);
      }
      
      if( nfiles++ == 0 )
	 proto_header = header;
      else if( strncmp(header.title, proto_header.title, L_name) ||
	    header.dump_interval != proto_header.dump_interval ||
	    header.dump_level != proto_header.dump_level       ||
	    header.dump_size != proto_header.dump_size         ||
	    header.dump_init != proto_header.dump_init )
      {
         fprintf(stderr,"Dump headers don't match: file\"%s\"\n", dump_name);
	 exit(2);
      };

      (void)fclose(dump_file);
      cur = (list_t *)calloc(1, sizeof(list_t));
      cur->p = dump_name;
      cur->i = header.istep/header.dump_interval;
      cur->num = header.ndumps;
      insert(cur, &f_head);
#ifdef DEBUG
      fprintf(stderr,"File \"%s\" \nslice %5d length %5d\n",
	              dump_name, cur->i, cur->num);
#endif
   }

   if( (cpt_mask & level_mask[proto_header.dump_level]) != cpt_mask )
   {
      fprintf(stderr,"Sorry the components requested by mask %d",cpt_mask);
      fprintf(stderr," are not contained in a dump of level %d\n",
	      header.dump_level);
      exit(2);
   }

   tslice = start = f_head.next->i;
   for(cur = f_head.next; cur; cur = cur->next)
   {
      if( cur->i != tslice )
      {
	 fprintf(stderr,"Dump file \"%s\" out of sequence at slice %d\n",
		        cur->p, tslice);
	 exit(2);
      }
      cur->i -= start;
      tslice += cur->num;
#ifdef DEBUG
      fprintf(stderr,"File \"%s\" \nslice %5d length %5d\n",
	              cur->p, cur->i, cur->num);
#endif
   }
   maxslice = tslice - start;

   if(DUMP_SIZE(header.dump_level, nmols,nmols_r) != proto_header.dump_size)
   {
      fprintf(stderr, "Number of molecules (%d/%d) ",nmols,nmols_r);
      fprintf(stderr, "incompatible with dump size (%d) and level(%d)\n",
	      proto_header.dump_size, proto_header.dump_level);
      exit(2);
   }
   /*
    *  Set up timestep mask
    */
   if( ! tsflg )
   {
      numslice = maxslice;
      tslice = 0;
      inc = 1;
   }
   else if( forstr(tsrange, &tslice, &numslice, &inc) )
   {
      fprintf(stderr,"Incorrect time slice selection \"%s\"\n", tsrange);
      exit(2);
   }
   else
   {
      numslice += 1;
      if(tslice < 0 || numslice > maxslice)
      {
	 fprintf(stderr, "Error in dump sequence - step %d not found\n", 
		 numslice);
	 exit(2);
      }
   }

   /*
    *  Molecule mask
    */
   if(mol_head.next == 0)
   {
      cur = (list_t*)calloc(1,sizeof(list_t));
      cur->i = 0;
      cur->num = nmols;
      mol_head.next = cur;
   }
   /*
    *  Set up component size and offsets
    */
#ifdef DEBUG
   fprintf(stderr,"Size Offset\n");
#endif
   offset=0;
   cpt[0].size = cpt[4].size = cpt[7].size = cpt[10].size = nmols;
   cpt[1].size = cpt[5].size = cpt[8].size = cpt[11].size = nmols_r;
   for(icpt = 0; icpt < NCPT; icpt++)
   {
      cpt[icpt].size *= cpt[icpt].ncpt;
      cpt[icpt].offset = offset;
      if( (1<<icpt) & level_mask[proto_header.dump_level] )
	 offset += cpt[icpt].size;
#ifdef DEBUG
      fprintf(stderr,"%6d %6d\n",cpt[icpt].size,cpt[icpt].offset);
#endif
   }
   /*
    *  Do the extraction
    */
   for( cur = f_head.next; cur; cur = cur->next)
   {
      if( cur->i <= tslice && tslice < MIN(cur->i + cur->num, numslice) )
      {
	 extract(cur->p, cpt_mask, mol_head.next, cpt, NCPT, tslice-cur->i,
		 MIN(cur->num,numslice-cur->i), inc, bflg, nmols);
	 tslice += (cur->i + cur->num - tslice - 1) / inc * inc + inc;
      }
   }

   return 0;
}
