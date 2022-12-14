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
#include	<stdlib.h>
#include	<stddef.h>
#include 	<string.h>
#include        <time.h>
#include 	<stdio.h>
#include	"structs.h"
#ifdef USE_XDR
#   include	"xdr.h"
#endif
#include        "utlsup.h"
#include	"messages.h"

#ifdef USE_XDR
   XDR          xdrs;
#endif

static int verbose;
/******************************************************************************
 * strstr replacement for pre-ANSI machines which don't have it.              *
 ******************************************************************************/
#ifndef STDC_HEADERS
char *strstr(cs, ct)
char *cs, *ct;
{
   char *end = cs+strlen(cs)-strlen(ct);
   for(; cs <= end; cs++)
      if( !strcmp(cs,ct) )
	 return cs;
   return 0;      
}
#endif

int av_convert;
#undef MIN
#define MIN(x,y) ( (x) > (y) ? (y) : (x))

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
 * List manipulation procedures						      *
 ******************************************************************************/
void
insert(list_mt *entry, list_mt *head)
{
   while( head->next != NULL && entry->i > head->next->i)
      head = head->next;

   entry->next = head->next;
   head->next  = entry;
}

void
print_list(list_mt *head)
{
   if(head == NULL)
      return;
   fprintf(stderr,"%-8d%s\n",head->i, head->p);
   print_list(head->next);
}

/******************************************************************************
 * Put.  Write data in text or binary form.				      *
 ******************************************************************************/
void
put(float *buf, int n, int bflg)
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

static
int read_dump_header(char *fname, FILE *dumpf, dump_mt *hdr_p, boolean *xdr_write,
		     size_mt sysinfo_size, dump_sysinfo_mt *dump_sysinfo)
{
   int      errflg = true;	/* Provisionally !!   */
   char     vbuf[sizeof hdr_p->vsn + 1];
   int	    vmajor,vminor;

   *xdr_write = false;
#ifdef USE_XDR
   /*
       * Attempt to read dump header in XDR format
       */
   if( xdr_dump(&xdrs, hdr_p) )
   {
      strncpy(vbuf,hdr_p->vsn,sizeof hdr_p->vsn);
      vbuf[sizeof hdr_p->vsn] = '\0';
      if( strstr(vbuf,"(XDR)") )
      {
	 errflg = false;
	 *xdr_write = true;
      }
   }
#endif
   /*
    * If we failed, try to read header as native struct image.
    */
   if( ! *xdr_write )
   {
      if( fseek(dumpf, 0L, 0) ) 
	 message(NULLI, NULLP, WARNING, SEFAIL, fname, strerror(errno));
      else if( fread((gptr*)&*hdr_p, sizeof(dump_mt), 1, dumpf) == 0 )
	 message(NULLI, NULLP, WARNING, DRERR, fname, strerror(errno));
      else
	 errflg = false;
   }
   if( ! errflg )
   {
      /*
       * Parse header version
       */
      errflg = true;
      if( sscanf(hdr_p->vsn, "%d.%d", &vmajor, &vminor) < 2 )
	 message(NULLI, NULLP, WARNING, INDVSN, hdr_p->vsn);
      if( vmajor < 2 || vminor <= 17)
	 message(NULLI, NULLP, WARNING, OLDVSN, hdr_p->vsn);
      else
	 errflg = false;
   }
   if( errflg ) return errflg;

   
   if( dump_sysinfo == 0)
      return errflg;
   else if ( sysinfo_size == sizeof(dump_sysinfo_mt) )
   {
      /*
       * Now check for sysinfo and read fixed part of it.  This is needed to
       * determine species count to read the whole thing.
       */
#ifdef USE_XDR
      if( *xdr_write ) {
	 if( ! xdr_dump_sysinfo_hdr(&xdrs, dump_sysinfo) )
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
	 errflg = false;
      } else
#endif
      {
	 if( fread((gptr*)dump_sysinfo,sizeof(dump_sysinfo_mt), 1, dumpf) == 0)
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
	 errflg = false;
      }
   }
   else
   {
      /*
       * Now check for sysinfo and read it all.  N.B.  Buffer must be 
       * allocated to full expected size by prior call to read_dump_header.
       */
#ifdef USE_XDR
      if( *xdr_write ) {
	 if( ! xdr_dump_sysinfo(&xdrs, dump_sysinfo, vmajor, vminor) )
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
      if (sizeof(dump_sysinfo_mt) 
	  + sizeof(mol_mt) * (dump_sysinfo->nspecies-1) > sysinfo_size)
      {
	 /*
	  * We have already overrun the end of the "dump_sysinfo" buffer.
	  * Perhaps we can exit gracefully before crashing?
          */
	 message(NULLI, NULLP, FATAL, RDHERR,  sizeof(dump_sysinfo_mt) 
		 + sizeof(mol_mt) * (dump_sysinfo->nspecies-1), sysinfo_size);
      }
	 errflg = false;
      } else
#endif
      {
	 if( fread((gptr*)dump_sysinfo, sysinfo_size, 1, dumpf) == 0)
	    message(NULLI, NULLP, FATAL, DRERR, fname, strerror(errno));
	 errflg = false;
      }
   }
   return errflg;
}

static
void read_dump_record(gptr *dump_buf, FILE *dumpf, size_mt dump_size, 
		       char *cur_file, boolean xdr_write)
{
#ifdef USE_XDR
   if( xdr_write )
   {
      if( ! xdr_vector(&xdrs, dump_buf, dump_size, sizeof(float), 
		     (xdrproc_t)xdr_float) )
	 message(NULLI, NULLP, FATAL, DRERR, cur_file, strerror(errno));
   }
   else
#endif
   {
      if( fread(dump_buf, sizeof(float), dump_size, dumpf) ==0 )
	 message(NULLI, NULLP, FATAL, DRERR, cur_file, strerror(errno));
   }
}

/******************************************************************************
 * Extract.  Process one dump file, extracting and outputting data.	      *
 ******************************************************************************/
void
extract(char *dump_name, int cpt_mask, list_mt *molecules, cpt_mt *cpt, int ncpt, 
	int tslice, int num, int inc, int bflg, int nmols, int xdr, size_mt sysinfo_size)
{
   FILE		*dump_file;
   dump_mt	header;
   float	*buf = (float*)calloc(4*nmols,sizeof(float));/* nmols > nmols_r */
   long		dump_base;
   int		icpt, start, nitems;
   list_mt	*mol;
   dump_sysinfo_mt *dump_sysinfo =  (dump_sysinfo_mt*)malloc(sysinfo_size);
   
   if( (dump_file = open_dump(dump_name, "rb")) == NULL)
   {
      fprintf(stderr, "Failed to open dump file \"%s\"\n", dump_name);
      exit(2);
   }
#ifdef DEBUG
   fprintf(stderr,"Working on file \"%s\" (%d-%d)\n", dump_name, tslice, num);
#endif
         
   if( read_dump_header(dump_name, dump_file, &header, &xdr, 
			sysinfo_size, dump_sysinfo) 
       || ferror(dump_file) || feof(dump_file) )
   {
      fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
      exit(2);
   }

#ifdef USE_XDR
   if( xdr )
      dump_base = XDR_DUMP_SIZE + sysinfo_size
	                        + tslice*header.dump_size*XDR_FLOAT_SIZE;
   else
#endif
      dump_base = sizeof(dump_mt) + sysinfo_size 
	                        + tslice*header.dump_size*sizeof(float);
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
			  (xdrproc_t)xdr_float);
	    }
	    else
#endif
	    {
	       fseek(dump_file, dump_base+cpt[icpt].offset*sizeof(float), 0);
	       fread((char*)buf, sizeof(float), cpt[icpt].size, dump_file);
	    }
	    if( cpt[icpt].mols )
	       for(mol = molecules; mol != 0; mol = mol->next)
	       {
		  start  = mol->i   * cpt[icpt].ncpt;
		  nitems = mol->num * cpt[icpt].ncpt;
		  if( start+nitems > cpt[icpt].size )
		     nitems = cpt[icpt].size - start;
		  put(buf+start, nitems, bflg);
	       }
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
   (void)close_dump(dump_file);
   (void)free((char*)buf);
   (void)free((char*)dump_sysinfo);
}

void	print_header(dump_mt *header, dump_sysinfo_mt *sysinfo)
{
   int ispec;
   printf("Title\t\t\t\t= \"%s\"\n",header->title);
   printf("RCS revision\t\t\t= %.*s\n", (int)strlen(header->vsn), header->vsn);
   printf("Istep\t\t\t\t= %ld\n",header->istep);
   printf("Dump interval\t\t\t= %ld\n", header->dump_interval);
   printf("Dump level\t\t\t= %d\n", header->dump_level);
   printf("Max dumps\t\t\t= %d\n", header->maxdumps);
   printf("Dump size\t\t\t= %d\n", header->dump_size);
   printf("Number of dumps\t\t\t= %d\n", header->ndumps);
   printf("Timestamp\t\t\t= %s", ctime((time_t*)&header->timestamp));
   printf("Restart timestamp\t\t= %s", ctime((time_t*)&header->restart_timestamp));
   printf("Dump start\t\t\t= %s", ctime((time_t*)&header->dump_init));
   printf("Time between dumps\t\t= %f %s\n", sysinfo->deltat, TUNIT_N);
   printf("Number of particles\t\t= %d\n", sysinfo->nmols);
   printf("Number of polyatomics\t\t= %d\n", sysinfo->nmols_r);
   printf("Number of species\t\t= %d\n", sysinfo->nspecies);
   for(ispec = 0; ispec < sysinfo->nspecies; ispec++)
   {
      printf("Species %d name\t\t\t= %s\n", ispec+1, sysinfo->mol[ispec].name);
      if( sysinfo->mol[ispec].charge != 0.0)
        printf("  Number of ions\t\t= %d\n",sysinfo->mol[ispec].nmols);
      else
        printf("  Number of molecules\t\t= %d\n",sysinfo->mol[ispec].nmols);
      if( sysinfo->mol[ispec].framework )
      {
         if( sysinfo->mol[ispec].charge != 0.0)
	   printf("  Ion is a framework\n");
         else
	   printf("  Molecule is a framework\n");
      }
      else
	 printf("  Rotational deg. of freedom\t= %d\n", sysinfo->mol[ispec].rdof);
      printf("  Mass\t\t\t\t= %f amu\n", sysinfo->mol[ispec].mass);
      if (sysinfo->mol[ispec].rdof > 0 )
	 printf("  Moments of inertia\t\t= %f  %f  %f %s\n", sysinfo->mol[ispec].inertia[0],
		sysinfo->mol[ispec].inertia[1],sysinfo->mol[ispec].inertia[2], IUNIT_N);
      printf("  Charge\t\t\t= %f %s\n", CONV_Q*sysinfo->mol[ispec].charge, CONV_Q_N); 
      printf("  Dipole moment\t\t\t= %f %s\n", CONV_D*sysinfo->mol[ispec].dipole, CONV_D_N);
   }
}

int
main(int argc, char **argv)
{
   int	c;
   extern int	optind;
   extern char	*optarg;
   int		errflg = 0, genflg = 0, tsflg = 0, bflg = 0;
   int          nmols=-1, nmols_r=-1;
   int		xcpt= -2;
   char		*dump_name=0, *dump_base=0, *out_name=0;
   char		cur_dump[256];
   char		*tsrange;
   FILE		*dump_file;
   int		nfiles = 0;
   int		start,finish,inc;
   int		tslice, numslice, maxslice;
   int		offset, icpt;
   int		idump0;
   int		xdr = 0;

   dump_sysinfo_mt *dump_sysinfo;
   size_mt	sysinfo_size;
   
   static cpt_mt cpt[] = {{3, 0, 3, 1, "C of M positions"},
			  {4, 0, 4, 1, "quaternions"},
			  {9, 0, 1, 0, "unit cell matrix"},
			  {1, 0, 1, 0, "Thermostat variable"},
			  {1, 0, 1, 0, "potential energy"},
			  {3, 0, 3, 1, "C of M velocities"},
			  {3, 0, 3, 1, "angular velocities"},
			  {9, 0, 1, 0, "unit cell velocities"},
			  {1, 0, 1, 0, "Thermostat momentum"},
			  {3, 0, 3, 1, "C of M forces"},
			  {3, 0, 3, 1, "torques"},
			  {9, 0, 1, 0, "stress tensor"} };
#define NCPT (int)(sizeof(cpt)/sizeof(cpt_mt))

   static int level_mask[16] = {  0x0000,0x001f,0x01e0,0x01ff,
				  0x0000,0x001f,0x01e0,0x01ff,
				  0x0e00,0x0e1f,0x0fe0,0x0fff,
				  0x0e00,0x0e1f,0x0fe0,0x0fff};

   dump_mt	proto_header, header;

   list_mt	f_head;
   list_mt	mol_head;
   list_mt	*cur;

   mol_head.next = NULL;
   f_head.next = NULL;

   verbose = 0;

   while( (c = getopt(argc, argv, "c:br:R:q:Q:t:m:o:v") ) != EOF )
      switch(c)
      {
       case 'c':
	 xcpt = strtol(optarg,(char**)0,0);
	 break;
       case 'b':
	 bflg++;
	 break;
       case 'r':
       case 'R':
	 nmols = strtol(optarg,(char**)0,0);
	 break;
       case 'q':
       case 'Q':
	 nmols_r = strtol(optarg,(char**)0,0);
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
	    if( inc != 1 )
	    {
	       fprintf(stderr,"\":%d\" modifier not allowed for molecule ranges\n",inc);
	       errflg++;
	    }
	    cur = (list_mt *)calloc(1, sizeof(list_mt));
	    cur->i = start;  cur->num=finish-start+1;
	    insert(cur, &mol_head);
	 }	 
	 break;
       case 'v':
         verbose++;
	 break;
       case '?': 
       case 'h':
	 errflg++;
	 break;
      }


   if(optind >= argc)
   {
      fprintf(stderr,"%s: no dump files given\n",argv[0]);
      errflg++;
   }
   if( errflg )
   {
      fprintf(stderr,
	   "Usage: dumpext [-Rn] [-Qn] [-b] [-c cpt]\
 [-t timeslices] [-m molecules] [-v] [-o output-file] dumpfiles\n");
      exit(2);
   }
   /*
    *  Interactive input of parameters not supplied as argument
    */
#if 0
   if( nmols <= 0)
      nmols = get_int("Number of molecules? ",1,1000000);
   if( nmols_r < 0)
      nmols_r = get_int("Number of polyatomic molecules? ",0,1000000);
#endif
   if( xcpt < -1 )
   {
      fprintf(stderr,"Which quantity do you require?\n");
      fprintf(stderr,"\t%-32s %d\n","Header information",-1);
      fprintf(stderr,"\t%-32s %d\n","All data components",0);
      for(icpt = 0; icpt < NCPT; icpt++)
	 fprintf(stderr,"\t%-32s %d\n",cpt[icpt].name,icpt+1);
      xcpt=get_int("Quantity index (0-12)? ",-1,NCPT);
   }
   /*
    *  Generate list of dump files if required
    */
   if (strstr(argv[optind],"%d") )
   {
      dump_base = argv[optind];
      genflg++;
#define MAXTRY 500
      idump0 = -1;
      if ( verbose )  fprintf(stderr,"Searching for dump file matching \"%s\"\n",dump_base);
      do                      /* Search for a dump file matching pattern */
        sprintf(cur_dump, dump_base, ++idump0);
      while( (dump_file = fopen(cur_dump, "rb")) == NULL && idump0 < MAXTRY);
      if( dump_file == NULL )        /* If we didn't find one . .               */
      {
        fprintf(stderr,"No dump files found matching \"%s\".\n",dump_base);
        exit(2);
      } else if (verbose) {
        fprintf(stderr,"... found \"%s\"\n", cur_dump);
      }
      (void)fclose(dump_file);
   }
   else
      idump0 = optind;
   
   /*
    *  Check all dump files for correctness and build ordered list
    */
   while(1)
   {
      if( genflg )
      {
	 sprintf(cur_dump, dump_base, idump0++);
	 dump_name = cur_dump;
      }
      else
      {
	 dump_name = argv[idump0++];
	 if( dump_name == 0 )
	    break;
      }

      if( verbose )
         fprintf(stderr, "Checking dump file \"%s\"\n",dump_name);
      if( (dump_file = open_dump(dump_name, "rb")) == NULL)
      {
	 if( genflg )
         {
           if( verbose )
              fputs("End of dump sequence reached\n",stderr);
	    break;		/* Exit loop if at end of sequence */
         }
	 fprintf(stderr, "Failed to open dump file \"%s\"\n", dump_name);
	 exit(2);
      }

      /*
       * Read dump header. On first call we need to read only the first
       * part of sysinfo to determine nspecies and consequently the size
       * of the buffer needed to hold all of it.
       */
      sysinfo_size = sizeof(dump_sysinfo_mt);
      dump_sysinfo = (dump_sysinfo_mt*)malloc(sysinfo_size);
      if( read_dump_header(dump_name, dump_file, &header, &xdr, 
			   sysinfo_size, dump_sysinfo) 
	                       || ferror(dump_file) || feof(dump_file) )
      {
	 fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
	 exit(2);
      }
      sysinfo_size = sizeof(dump_sysinfo_mt) 
	                   + sizeof(mol_mt) * (dump_sysinfo->nspecies-1);
      
     /* Check that dump headers match */
      if( nfiles++ == 0 )
      {
	 proto_header = header;
      }
      else if( strncmp(header.title, proto_header.title, L_name) ||
	      strncmp(header.vsn, proto_header.vsn, sizeof header.vsn) ||
	    header.dump_interval != proto_header.dump_interval ||
	    header.dump_level != proto_header.dump_level       ||
	    header.dump_size != proto_header.dump_size         ||
	    header.dump_init != proto_header.dump_init )
      {
         fprintf(stderr,"Dump headers don't match: file\"%s\"\n", dump_name);
	 exit(2);
      } else if (verbose) {
        fprintf(stderr,"Dump file \"%s\" contains slices %ld to %ld\n",dump_name,
              header.istep/header.dump_interval,header.istep/header.dump_interval+header.ndumps-1);
      }
     /*
      * Allocate space for and read dump sysinfo.
      */
      (void)free(dump_sysinfo);
      dump_sysinfo = (dump_sysinfo_mt*)malloc(sysinfo_size);
     /*
      * Rewind and reread header, this time including sysinfo.
      */
      (void)rewind_dump(dump_file, xdr);
      if( read_dump_header(dump_name, dump_file, &header, &xdr, 
        sysinfo_size, dump_sysinfo) 
          || ferror(dump_file) || feof(dump_file) )
      {
	fprintf(stderr, "Failed to read dump header \"%s\"\n", dump_name);
	exit(2);
      }
      nmols   = dump_sysinfo->nmols;
      nmols_r = dump_sysinfo->nmols_r;

      (void)close_dump(dump_file);

      cur = (list_mt *)calloc(1, sizeof(list_mt));
      cur->p = mystrdup(dump_name);
      cur->i = header.istep/header.dump_interval;
      cur->num = header.ndumps;
      insert(cur, &f_head);
      if( verbose )
         fprintf(stderr, "Dump file \"%s\" added to list to be processed\n",dump_name);
#ifdef DEBUG
      fprintf(stderr,"File \"%s\" \nslice %5d length %5d\n",
	              dump_name, cur->i, cur->num);
#endif
     if( xcpt == -1 )
     {
       if( bflg ) /* Binary output */
       {
	 fwrite(&header, sizeof(header), 1, stdout);
	 fwrite(dump_sysinfo, sizeof(*dump_sysinfo), 1, stdout);
       }
       else
         {  /* Text output */
         if( nfiles > 1 )
           printf("\n");
         printf("File name\t\t\t= %s\n", dump_name);
 	 print_header(&header, dump_sysinfo);
         }
     }
   }

   if( xcpt == -1 )
     exit(0);

   /*
    * Check molecule selections are in range.
    */
   for(cur=mol_head.next; cur; cur = cur->next)
   {
      if( cur->i < 0 || cur->i+cur->num > nmols)
      {
	 fprintf(stderr, "Error in molecule selection: \"%d-%d\" out of range.\n",
		cur->i, cur->i+cur->num-1);
	 exit(2);
      }
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

   if( xcpt > 0 && ! (1 << (xcpt-1) & level_mask[proto_header.dump_level]) )
   {
      fprintf(stderr,"The component requested (%s)",cpt[xcpt-1].name);
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
		 numslice-1);
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
   cpt[0].size = cpt[5].size = cpt[9].size = nmols;
   cpt[1].size = cpt[6].size = cpt[10].size = nmols_r;
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
    * Open output file if requested
    */
   if( out_name )
      if( ! freopen(out_name, bflg?"wb":"w", stdout) )
      {
	 fprintf(stderr,"Failed to open file \"%s\" for output - ",out_name);
	 exit(4);
      }

   /*
    *  Do the extraction
    */
   for( cur = f_head.next; cur; cur = cur->next)
   {
      if( cur->i <= tslice && tslice < MIN(cur->i + cur->num, numslice) )
      {
	 extract(cur->p, xcpt?1<<(xcpt-1):~0, mol_head.next, cpt, NCPT, 
		 tslice-cur->i,
		 MIN(cur->num,numslice-cur->i), inc, bflg, nmols, xdr,
		 sysinfo_size);
	 tslice += (cur->i + cur->num - tslice - 1) / inc * inc + inc;
      }
   }

   if( verbose )
      fprintf(stderr, "Data successfully extracted from \"%s\"\n",argv[optind]);

   return 0;
}
