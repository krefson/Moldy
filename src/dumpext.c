/* MOLecular DYnamics simulation code, Moldy.
Copyright (C) 1988, 1992, 1993 Keith Refson
 
This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 
In other words, you are welcome to use, share and improve this program.
You are forbidden to forbid anyone else to use, share and improve
what you give them.   Help stamp out software-hoarding!  */
#include "defs.h"
#include	"stdlib.h"
#include	"stddef.h"
#include 	"string.h"
#include 	<stdio.h>
#include	"structs.h"
#ifdef USE_XDR
#ifdef sun
#   define free xxfree
#   include	"xdr.h"
#   undef free
#else
#   include	"xdr.h"
#endif
#endif

/******************************************************************************
 * strstr replacement for pre-ANSI machines which don't have it.              *
 ******************************************************************************/
#ifndef ANSI_LIBS
char *strstr(cs, ct)
char *cs, *ct;
{
   char *end = cs+strlen(cs)-strlen(ct);
   for(; cs < end; cs++)
      if( !strcmp(cs,ct) )
	 return cs;
   return 0;      
}
#endif

int av_convert;
#undef MIN
#define MIN(x,y) ( (x) > (y) ? (y) : (x))

/*========================== Macros ==========================================*/
#define DUMP_SIZE(level, n, n_r)  \
   			 (( (level & 1) + (level>>1 & 1) + (level>>2 & 1) ) * \
                                    (3*n + 4*n_r + 9)+ \
                             (level>>3 & 1) * \
                                    (3*n + 3*n_r + 9) +\
                             (level & 1))
/*============================================================================*/
typedef struct list_mt
{
   struct list_mt	*next;
   int			i;
   char *p;
   int num;
} list_mt;

typedef struct cpt_mt
{
   int ncpt, offset, size, mols;
   char name[32];
} cpt_mt;

/******************************************************************************
 * get_int().  Read an integer from stdin, issuing a prompt and checking      *
 * validity and range.  Loop until satisfied, returning EOF if appropriate.   *
 ******************************************************************************/
int get_int(prompt, lo, hi)
char    *prompt;
int     lo, hi;
{
   char         ans_str[80];
   int          ans_i, ans_flag;

   ans_flag = 0;
   while( ! feof(stdin) && ! ans_flag )
   {
      fputs(prompt, stderr);
      fflush(stderr);
      fgets(ans_str, sizeof ans_str, stdin);
      if( sscanf(ans_str, "%d", &ans_i) == 1 && ans_i >= lo && ans_i <= hi)
         ans_flag++;
   }
   if( feof(stdin) )
   {
      fprintf(stderr,"\nExit requested\n");
      exit(3);
   }
   if( ans_flag )
      return(ans_i);
   else
      return(EOF);
}
/******************************************************************************
 * List manipulation procedures						      *
 ******************************************************************************/
void
insert(entry, head)
list_mt	*entry, *head;
{
   while( head->next != NULL && entry->i > head->next->i)
      head = head->next;

   entry->next = head->next;
   head->next  = entry;
}

void
print_list(head)
list_mt	*head;
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
extract(dump_name, cpt_mask, molecules, cpt, ncpt, tslice, num, inc, 
	bflg, nmols, xdr)
char	*dump_name;
int	cpt_mask;
list_mt	*molecules;
cpt_mt	cpt[];
int	ncpt, tslice, num, inc;
int	bflg, nmols, xdr;
{
   FILE		*dump_file;
   dump_mt	header;
   float	*buf = (float*)calloc(4*nmols,sizeof(float));/* nmols > nmols_r */
   long		dump_base;
   int		icpt;
   list_mt	*mol;
   int		errflg = 0;
#ifdef USE_XDR
   XDR          xdrs;
#endif

   
   if( (dump_file = fopen(dump_name, "rb")) == NULL)
   {
      fprintf(stderr, "Failed to open dump file \"%s\"\n", dump_name);
      exit(2);
   }
#ifdef DEBUG
   fprintf(stderr,"Working on file \"%s\" (%d-%d)\n", dump_name, tslice, num);
#endif
#ifdef USE_XDR
   /*
    * Attempt to read dump header in XDR format
    */
   if( xdr )
   {
      xdrstdio_create(&xdrs, dump_file, XDR_DECODE);
      errflg = ! xdr_dump(&xdrs, &header);
   }
   else
#endif
   {
      if( fread((char*)&header, sizeof(dump_mt), 1, dump_file) == 0 )
	 errflg = false;
   }
         
   if( errflg || ferror(dump_file) || feof(dump_file) )
   {
      fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
      exit(2);
   }

#ifdef USE_XDR
   if( xdr )
      dump_base = XDR_DUMP_SIZE+tslice*header.dump_size*XDR_FLOAT_SIZE;
   else
#endif
      dump_base = sizeof(dump_mt)+tslice*header.dump_size*sizeof(float);
   while(tslice < num && tslice < header.ndumps)
   {
#ifdef DEBUG
      fprintf(stderr,"Timeslice %d of %d\n", tslice, num);
#endif
      for( icpt = 0; icpt < ncpt; icpt++)
      {
	 if( cpt_mask & (1 << icpt) )
	 {
#ifdef USE_XDR
	    if( xdr )
	    {
	       xdr_setpos(&xdrs, dump_base+cpt[icpt].offset*XDR_FLOAT_SIZE);
	       xdr_vector(&xdrs, (char*)buf, (u_int)cpt[icpt].size, XDR_FLOAT_SIZE, 
			  xdr_float);
	    }
	    else
#endif
	    {
	       fseek(dump_file, dump_base+cpt[icpt].offset*sizeof(float), 0);
	       fread((char*)buf, sizeof(float), cpt[icpt].size, dump_file);
	    }
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
#ifdef USE_XDR
      if( xdr )
	 dump_base += inc*XDR_FLOAT_SIZE * header.dump_size;
      else
#endif
	 dump_base += inc*sizeof(float) * header.dump_size;
   }
#ifdef USE_XDR
   if( xdr )
      xdr_destroy(&xdrs);
#endif
   (void)fclose(dump_file);
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
   int          nmols=0, nmols_r=0;
   int		xcpt=0;
   char		*dump_name=0, *out_name=0;
   char		*tsrange, *filerange;
   char		**filelist;
   FILE		*dump_file;
   int		nfiles = 0;
   int		start,finish,inc;
   int		tslice, numslice, maxslice;
   int		find, offset, icpt;
   int		xdr = 0;
#ifdef USE_XDR
   XDR          xdrs;
#endif
   
   static cpt_mt cpt[] = {{3, 0, 3, 1, "C of M positions"},
			 {4, 0, 4, 1, "quaternions"},
			 {9, 0, 1, 0, "unit cell matrix"},
			 {1, 0, 1, 0, "potential energy"},
			 {3, 0, 3, 1, "C of M velocities"},
			 {4, 0, 4, 1, "quaternion derivatives"},
			 {9, 0, 1, 0, "unit cell velocities"},
			 {3, 0, 3, 1, "C of M accelerations"},
			 {4, 0, 4, 1, "quaternion accelerations"},
			 {9, 0, 1, 0, "unit cell accelerations"},
			 {3, 0, 3, 1, "C of M forces"},
			 {3, 0, 3, 1, "torques"},
			 {9, 0, 1, 0, "stress tensor"} };
#define NCPT (sizeof(cpt)/sizeof(cpt_mt))

   static int level_mask[16] = {  0x0000,0x000f,0x0070,0x007f,
				  0x0380,0x038f,0x03f0,0x03ff,
				  0x1c00,0x1c0f,0x1c70,0x1c7f,
				  0x1f80,0x1f8f,0x1ff0,0x1fff};

   dump_mt	proto_header, header;

   list_mt	f_head;
   list_mt	mol_head;
   list_mt	*cur;

   mol_head.next = NULL;
   f_head.next = NULL;

   while( (c = getopt(argc, argv, "c:bR:Q:n:t:m:o:") ) != EOF )
      switch(c)
      {
       case 'c':
	 xcpt = strtol(optarg,(char**)0,0);
	 break;
       case 'b':
	 bflg++;
	 break;
       case 'R':
	 nmols = strtol(optarg,(char**)0,0);
	 break;
       case 'Q':
	 nmols_r = strtol(optarg,(char**)0,0);
	 break;
       case 'n':
	 filerange = optarg;
	 genflg++;
	 break;
       case 'o':
	 out_name = optarg;
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
	    cur = (list_mt *)calloc(1, sizeof(list_mt));
	    cur->i = start;  cur->num=finish-start+1;
	    insert(cur, &mol_head);
	 }	 
	 break;
       case '?': 
       case 'h':
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
	   "Usage: dumpextract [-Rn] [-Qn] [-b] [-n dumpfile-range] [-c cpt]\
 [-t timeslices] [-m molecules] [-o output-file] dumpfiles\n");
      exit(2);
   }
   /*
    *  Interactive input of parameters not supplied as argument
    */
   if( ! nmols )
      nmols = get_int("Number of molecules? ",1,1000000);
   if( ! nmols_r )
      nmols_r = get_int("Number of polyatomic molecules? ",0,1000000);
   if( ! xcpt )
   {
      fprintf(stderr,"Which quantity do you require?\n");
      for(icpt = 0; icpt < NCPT; icpt++)
	 fprintf(stderr,"\t%-32s %d\n",cpt[icpt].name,icpt+1);
      xcpt=get_int("Quantity index (1-13)? ",1,NCPT);
   }

   /*
    *  Molecule mask
    */
   if(mol_head.next == 0)
   {
      cur = (list_mt*)calloc(1,sizeof(list_mt));
      cur->i = 0;
      cur->num = nmols;
      mol_head.next = cur;
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
#ifdef USE_XDR
      /*
       * Attempt to read dump header in XDR format
       */
      xdrstdio_create(&xdrs, dump_file, XDR_DECODE);
      if( xdr_dump(&xdrs, &header) )
      {
	 header.vsn[sizeof header.vsn - 1] = '\0';
	 if( strstr(header.vsn,"(XDR)") )
	 {
	    errflg = 0;
	    xdr = 1;
	 }
      }
      else
	 errflg = 1;
#endif
      /*
       * If we failed, try to read header as native struct image.
       */
      if( ! xdr )
      {
	 if( fseek(dump_file, 0L, 0) == 0 &&
             fread((char*)&header, sizeof(dump_mt), 1, dump_file) == 1) 
	    errflg = 0;
      }
      if( errflg || ferror(dump_file) || feof(dump_file) )
      {
	 fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
	 exit(2);
      }
      
      if( nfiles++ == 0 )
	 proto_header = header;
      else if( strncmp(header.title, proto_header.title, L_name) ||
	      strncmp(header.vsn, proto_header.vsn, sizeof header.vsn) ||
	    header.dump_interval != proto_header.dump_interval ||
	    header.dump_level != proto_header.dump_level       ||
	    header.dump_size != proto_header.dump_size         ||
	    header.dump_init != proto_header.dump_init )
      {
         fprintf(stderr,"Dump headers don't match: file\"%s\"\n", dump_name);
	 exit(2);
      };

#ifdef USE_XDR
      if( xdr )
	 xdr_destroy(&xdrs);
#endif
      (void)fclose(dump_file);
      cur = (list_mt *)calloc(1, sizeof(list_mt));
      cur->p = dump_name;
      cur->i = header.istep/header.dump_interval;
      cur->num = header.ndumps;
      insert(cur, &f_head);
#ifdef DEBUG
      fprintf(stderr,"File \"%s\" \nslice %5d length %5d\n",
	              dump_name, cur->i, cur->num);
#endif
   }

   if( ! (1 << xcpt-1 & level_mask[proto_header.dump_level]) )
   {
      fprintf(stderr,"Sorry the component requested (%s)",cpt[xcpt-1].name);
      fprintf(stderr," is not contained in a dump of level %d\n",
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
    *  Safety.  Have we correctly calculated the dump record length?
    */
   if( offset != proto_header.dump_size )
   {
      fprintf(stderr,"Internal error in dump record length (%d,%d)\n",
	      offset,proto_header.dump_size);
      exit(3);
   }
   /*
    * Open output file  if requested
    */
   if( out_name )
      if( ! freopen(out_name, bflg?"wb":"w", stdout) )
      {
	 fprintf(stderr,"Failed to open file \"%s\" for output - ",out_name);
	 perror("");
	 exit(4);
      }
   /*
    *  Do the extraction
    */
   for( cur = f_head.next; cur; cur = cur->next)
   {
      if( cur->i <= tslice && tslice < MIN(cur->i + cur->num, numslice) )
      {
	 extract(cur->p, 1<<xcpt-1, mol_head.next, cpt, NCPT, tslice-cur->i,
		 MIN(cur->num,numslice-cur->i), inc, bflg, nmols, xdr);
	 tslice += (cur->i + cur->num - tslice - 1) / inc * inc + inc;
      }
   }

   return 0;
}
