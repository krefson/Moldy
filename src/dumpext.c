#include 	"string.h"
#include 	<stdio.h>
#include	"structs.h"

char	*calloc();

typedef struct list_t
{
   struct list_t	*next;
   int			i;
   char			*p;
} list_t;

static int	nmols, nmols_r;

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
   char outbuf[32];
   char *gcvt();
   if( bflg )
      fwrite((char*)buf, sizeof(float), n, stdout);
   else
      while(n-- > 0)
      {
	 fputs(gcvt(*buf++, 7, outbuf), stdout);
	 putchar(' ');
      }

}
/******************************************************************************
 * Extract.  Process one dump file, extracting and outputting data.	      *
 ******************************************************************************/
void
extract(dump_name, pflg, rlist, qlist, bflg)
char	*dump_name;
int	pflg;
list_t	*rlist, *qlist;
int	bflg;
{
   FILE		*dump_file;
   dump_t	header;
   list_t	*r_ptr, *q_ptr;
   float	r_buf[3], q_buf[3];
   long		idump, dump_base;
   long		q_offset = 3 * sizeof(float) * nmols,
   		h_offset = q_offset + 4 * sizeof(float) * nmols_r,
   		p_offset = h_offset + 9 * sizeof(float);
   
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

   dump_base = sizeof(dump_t);
   for( idump = 0; idump < header.ndumps; idump++)
   {
      if(pflg)
      {
	 fseek(dump_file, dump_base + p_offset, 0);
	 fread((char*)r_buf, sizeof(float), 1, dump_file);
	 put(r_buf, 1, bflg);
      }
      for(r_ptr = rlist; r_ptr; r_ptr = r_ptr->next)
      {
	 fseek(dump_file, dump_base + 3*sizeof(float)*r_ptr->i, 0);
	 fread((char*)r_buf, sizeof(float), 3, dump_file);
	 put(r_buf, 3, bflg);
      }
      for(q_ptr = qlist; q_ptr; q_ptr = q_ptr->next)
      {
	 fseek(dump_file, dump_base + q_offset + 4*sizeof(float)*q_ptr->i, 0);
	 fread((char*)q_buf, sizeof(float), 4, dump_file);
	 put(q_buf, 4, bflg);
      }
      if( ! bflg )
	 putchar('\n');
      dump_base += sizeof(float) * header.dump_size;
   }
}

main(argc, argv)
int	argc;
char	*argv[];
{
   int	c;
   extern int	optind;
   extern char	*optarg;
   int		errflg = 0, genflg = 0;
   char		*dump_name;
   char		*range;
   char		**filelist;
   FILE		*dump_file;
   int		nfiles = 0, pflg = 0, bflg = 0;
   int		start,finish,inc, i;
   int		find;

   dump_t	proto_header, header;

   list_t	r_head, q_head, f_head;
   list_t	*cur;

   r_head.next = q_head.next = f_head.next = NULL;

   while( (c = getopt(argc, argv, "pR:r:Q:q:n:b") ) != EOF )
      switch(c)
      {
       case 'p':
	 pflg++;
	 break;
       case 'b':
	 bflg++;
	 break;
       case 'R':
	 nmols = atoi(optarg);
	 break;
       case 'Q':
	 nmols_r = atoi(optarg);
	 break;
       case 'n':
	 range = optarg;
	 genflg++;
	 break;
       case 'r':
	 if( forstr(optarg, &start, &finish, &inc) )
	    errflg++;
	 else
	    for(i = start; i <= finish; i+= inc)
	    {
	       cur = (list_t *)calloc(1, sizeof(list_t));
	       cur->i = i;
	       if( cur->i < 0 || cur->i >= nmols)
	       {
		  fprintf(stderr,"Molecule index %d out of range\n", cur->i);
		  exit(2);
	       }
	       insert(cur, &r_head);
	    }
	 break;
       case 'q':
	 if( forstr(optarg, &start, &finish, &inc) )
	    errflg++;
	 else
	    for(i = start; i <= finish; i+= inc)
	    {
	       cur = (list_t *)calloc(1, sizeof(list_t));
	       cur->i = i;
	       if( cur->i < 0 || cur->i >= nmols_r)
	       {
		  fprintf(stderr, "Quaternion index %d out of range\n", cur->i);
		  exit(2);
	       }
	       insert(cur, &q_head);
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
   if(genflg)
   {
      if( (filelist = seqgen(argv[optind], range)) == NULL)
      {
	 fprintf(stderr,"%s: invalid dump range \"%s\"\n", argv[0]);
	 errflg++;
      }
   }
   else
      filelist = argv+optind;

   if( errflg )
   {
      fprintf(stderr,
	   "Usage: dumpextract [-Rn] [-ri] [-Qn] [-qi] [-p] [-n range] dumpfiles\n");
      exit(2);
   }

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
      else
      {
	 if( strncmp(header.title, proto_header.title, L_name) ||
	    header.dump_interval != proto_header.dump_interval ||
	    header.dump_level != proto_header.dump_level       ||
	    header.dump_size != proto_header.dump_size         ||
	    header.dump_init != proto_header.dump_init )
	 {
	    fprintf(stderr,"Dump headers don't match: file\"%s\"\n", dump_name);
	    exit(2);
	 }
      }
      (void)fclose(dump_file);
      cur = (list_t *)calloc(1, sizeof(list_t));
      cur->p = dump_name;
      cur->i = header.istep;
      insert(cur, &f_head);
   }
   
   for( cur = f_head.next; cur; cur = cur->next)
   {
      extract(cur->p, pflg, r_head.next, q_head.next, bflg);
   }
   return 0;
}
