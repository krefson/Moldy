#include 	<stdio.h>
#include	"structs.h"

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

void
extract(dump_name, pflg, rlist, qlist)
char	*dump_name;
int	pflg;
list_t	*rlist, *qlist;
{
   FILE		*dump_file;
   dump_t	header;
   list_t	*r_ptr, *q_ptr;
   float	r_buf[3], q_buf[3];
   int		idump, dump_base;
   int		q_offset = 3 * sizeof(float) * nmols,
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
	 fread(&r_buf, sizeof(float), 1, dump_file);
	 printf("%12.7g ", r_buf[0]);
      }
      for(r_ptr = rlist; r_ptr; r_ptr = r_ptr->next)
      {
	 fseek(dump_file, dump_base + 3*sizeof(float)*r_ptr->i, 0);
	 fread(&r_buf, sizeof(float), 3, dump_file);
	 printf("%12.7g %12.7g %12.7g", r_buf[0], r_buf[1], r_buf[2]);
      }
      for(q_ptr = qlist; q_ptr; q_ptr = q_ptr->next)
      {
	 fseek(dump_file, dump_base + q_offset + 4*sizeof(float)*q_ptr->i, 0);
	 fread(&q_buf, sizeof(float), 4, dump_file);
	 printf("%12.7g %12.7g %12.7g %12.7g", q_buf[0], q_buf[1],
		                               q_buf[2], q_buf[3]);
      }
      putchar('\n');
      dump_base += sizeof(float) * header.dump_size;
   }
}

main(argc, argv)
int	argc;
char	*argv[];
{
   int	c, val;
   extern int	optind;
   extern char	*optarg;
   int		errflg = 0;
   char		*dump_name;
   FILE		*dump_file;
   int		nfiles = 0, pflg = 0;

   dump_t	proto_header, header;

   list_t	r_head, q_head, f_head;
   list_t	*cur;

   r_head.next = q_head.next = f_head.next = NULL;

   while( (c = getopt(argc, argv, "pR:r:Q:q:") ) != EOF )
      switch(c)
      {
       case 'p':
	 pflg++;
	 break;
       case 'R':
	 nmols = atoi(optarg);
	 break;
       case 'Q':
	 nmols_r = atoi(optarg);
	 break;
       case 'r':
	 cur = (list_t *)calloc(1, sizeof(list_t));
	 cur->i = atoi(optarg);
	 if( cur->i < 0 || cur->i >= nmols)
	 {
	    fprintf(stderr, "Molecule index %d out of range\n", cur->i);
	    exit(2);
	 }
	 insert(cur, &r_head);
	 break;
       case 'q':
	 cur = (list_t *)calloc(1, sizeof(list_t));
	 cur->i = atoi(optarg);
	 if( cur->i < 0 || cur->i >= nmols_r)
	 {
	    fprintf(stderr, "Quaternion index %d out of range\n", cur->i);
	    exit(2);
	 }
	 insert(cur, &q_head);
	 break;
       case '?':
	 errflg++;
      }

   if( errflg )
   {
      fprintf(stderr,
	      "Usage: dumpextract [-Rn] [-ri] [-Qn] [-qi] [-p] dumpfiles\n");
      exit(2);
   }

   while( optind < argc )
   {
      dump_name = argv[optind];
      optind++;
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
      extract(cur->p, pflg, r_head.next, q_head.next);
   }
}
