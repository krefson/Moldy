/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: ReadDCD.C,v $
 *      $Author: dalke $        $Locker:  $                $State: Exp $
 *      $Revision: 1.6 $      $Date: 1997/03/13 17:38:56 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * C routines to read and write binary DCD files (which use the goofy
 * FORTRAN UNFORMATTED format).  These routines are courtesy of
 * Mark Nelson (autographs available upon request and $1 tip).
 *
 ***************************************************************************/

#include <string.h>
#include "ReadDCD.h"

/************************************************************************/
/*									*/
/*				MACRO CHECK_FREAD			*/
/*									*/
/*	CHECK_FREAD is used to check the status of a file read.  It     */
/*   is passed the return code from read and a string to print out if   */
/*   an error is detected.  If an error is found, an error message is   */
/*   printed out and the program terminates.  This was made into a      */
/*   macro because it had to be done over and over and over . . .	*/
/*									*/
/************************************************************************/

#define CHECK_FREAD(X, msg)  if (X==-1) \
			     { \
				return(DCD_BADREAD); \
			     }

/************************************************************************/
/*									*/
/*			MACRO CHECK_FEOF				*/
/*									*/
/*	CHECK_FEOF checks for an abnormal EOF on a file read.  It is    */
/*   passed the return status from read and a message to print on error.*/
/*   if an EOF is detected, the error message is printed and the program*/
/*   terminates.  This is used for reads that should find something in  */
/*   the file.								*/
/*									*/
/************************************************************************/

#define CHECK_FEOF(X, msg)  if (X==0) \
			     { \
				return(DCD_BADEOF); \
			     }


/*			GLOBAL VARIABLES				*/
extern int errno;

void pad(char *s, int len)
{
	int curlen;
	int i;

	curlen=strlen(s);

	if (curlen>len)
	{
		s[len]=0;
		return;
	}

	for (i=curlen; i<len; i++)
	{
		s[i]=' ';
	}

	s[i]=0;
}

/************************************************************************/
/*									*/
/*			FUNCTION open_dcd_read				*/
/*									*/
/*   INPUTS:								*/
/*	filename - filename to read in					*/
/*									*/
/*   OUTPUTS:								*/
/*	an open filedesriptor (integer) is returned if the call is      */
/*   successful, otherwise a negative number is returned		*/
/*									*/
/*	open_dcd_read opens a dcd file for input.  It first does a check*/
/*   to see if the file really exists.  If the file does not exist,     */
/*   a DCD_DNE is returned, if the file exists but can' be opened for   */
/*   some reason, a DCD_OPENFAILED is returned.  If the open is         */
/*   successful, the file descriptor is returned.			*/
/*									*/
/************************************************************************/

int open_dcd_read(const char *filename)

{
	struct stat stbuf;	/*  Stat structure to check file	*/
	int dcdfd;		/*  file descriptor for dcd file	*/

	/*  Do a stat just to see if the file really exists	*/
	if (stat(filename, &stbuf) != 0)
	{
		if (errno == ENOENT)
		{
			return(DCD_DNE);
		}
	}

	/*  Try and open the file				*/
	dcdfd=open(filename, O_RDONLY);

	if (dcdfd == -1)
	{
		return(DCD_OPENFAILED);
	}

	return(dcdfd);
}

/****************************************************************/
/*								*/
/*			FUNCTION read_dcdheader			*/
/*								*/
/*   INPUTS:							*/
/*	fd - file descriptor for the dcd file			*/
/*	N - Number of atoms					*/
/*	NSET - Number of sets of coordinates			*/
/*	ISTART - Starting timestep of DCD file			*/
/*	NSAVC - Timesteps between DCD saves			*/
/*	DELTA - length of a timestep				*/
/*								*/
/*   OUTPUTS:							*/
/*	N, NSET, ISTART, NSAVC, and DELTA are returned populated*/
/*   a 0 is returned if the call is successful, or a negative   */
/*   number if errors are detected				*/
/*								*/
/*	read_header reads in the header from a DCD file and     */
/*   returns the timestep size and the number of atoms in the   */
/*   system.  A 0 is returned if the header is successfully     */
/*   read, or a negative error code is returned otherwise.      */
/*								*/
/****************************************************************/

int read_dcdheader(int fd, int *N, int *NSET, int *ISTART, 
		   int *NSAVC, double *DELTA, int *NAMNF, 
		   int **FREEINDEXES)

{
	int input_integer;	/*  Integer buffer space	*/
	char bigbuf[256];	/*  A large string buffer	*/
	int ret_val;		/*  Return value from read	*/
	int i;			/*  Loop counter		*/
	char HDR[5];		/*  Header = "CORD"		*/
	int I;
	int NTITLE;

	/*  First thing in the file should be an 84		*/
	ret_val = read(fd, &input_integer, sizeof(int));

	CHECK_FREAD(ret_val, "reading first int from dcd file");
	CHECK_FEOF(ret_val, "reading first int from dcd file");

	if (input_integer != 84)
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in the string "CORD"			*/
	ret_val = read(fd, HDR, 4);

	CHECK_FREAD(ret_val, "reading CORD from dcd file");
	CHECK_FEOF(ret_val, "reading CORD from dcd file");

	HDR[ret_val]=0;

	/*
// temp. fix until the read solution
//	if (strcmp(HDR, "CORD") != 0)
//	{
//		return(DCD_BADFORMAT);
//	}*/

	/*  Read in the number of Sets of coordinates, NSET  */
	ret_val = read(fd, NSET, sizeof(int));

	CHECK_FREAD(ret_val, "reading NSET");
	CHECK_FEOF(ret_val, "reading NSET");

	/*  Read in ISTART, the starting timestep	     */
	ret_val = read(fd, ISTART, sizeof(int));

	CHECK_FREAD(ret_val, "reading ISTART");
	CHECK_FEOF(ret_val, "reading ISTART");

	/*  Read in NSAVC, the number of timesteps between   */
	/*  dcd saves					     */
	ret_val = read(fd, NSAVC, sizeof(int));

	CHECK_FREAD(ret_val, "reading NSAVC");
	CHECK_FEOF(ret_val, "reading NSAVC");


	/*  Skip blank integers				     */
	for (i=0; i<5; i++)
	{
		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading I");
		CHECK_FEOF(ret_val, "reading I");
	}

	/*  Read NAMNF, the number of free atoms	     */
	ret_val = read(fd, NAMNF, sizeof(int));

	CHECK_FREAD(ret_val, "reading NAMNF");
	CHECK_FEOF(ret_val, "reading NAMNF");

	/*  Read in the timestep, DELTA				*/
	ret_val = read(fd, DELTA, sizeof(double));

	CHECK_FREAD(ret_val, "reading DELTA");
	CHECK_FEOF(ret_val, "reading DELTA");

	/*  Skip blank integers					*/
	for (i=0; i<9; i++)
	{
		ret_val = read(fd, &I, sizeof(int));

		CHECK_FREAD(ret_val, "reading I");
		CHECK_FEOF(ret_val, "reading I");
	}

	/*  Get the end size of the first block			*/
	ret_val = read(fd, &input_integer, sizeof(int));

	CHECK_FREAD(ret_val, "reading second 84 from dcd file");
	CHECK_FEOF(ret_val, "reading second 84 from dcd file");

	if (input_integer != 84)
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in the size of the next block			*/
	ret_val = read(fd, &input_integer, sizeof(int));

	CHECK_FREAD(ret_val, "reading size of title block");
	CHECK_FEOF(ret_val, "reading size of title block");

	if ( ((input_integer-4)%80) == 0)
	{
		/*  Read NTITLE, the number of 80 characeter    */
		/*  title strings there are			*/
		ret_val = read(fd, &NTITLE, sizeof(int));

		CHECK_FREAD(ret_val, "reading NTITLE");
		CHECK_FEOF(ret_val, "reading NTITLE");

		for (i=0; i<NTITLE; i++)
		{
			ret_val = read(fd, bigbuf, 80);

			CHECK_FREAD(ret_val, "reading TITLE");
			CHECK_FEOF(ret_val, "reading TITLE");
		}

		/*  Get the ending size for this block		*/
		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading size of title block");
		CHECK_FEOF(ret_val, "reading size of title block");
	}
	else
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in an 4				*/
	ret_val = read(fd, &input_integer, sizeof(int));

	CHECK_FREAD(ret_val, "reading an 4");
	CHECK_FEOF(ret_val, "reading an 4");

	if (input_integer != 4)
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in the number of atoms			*/
	ret_val = read(fd, N, sizeof(int));

	CHECK_FREAD(ret_val, "reading number of atoms");
	CHECK_FEOF(ret_val, "reading number of atoms");

	/*  Read in an 4				*/
	ret_val = read(fd, &input_integer, sizeof(int));

	CHECK_FREAD(ret_val, "reading an 4");
	CHECK_FEOF(ret_val, "reading an 4");

	if (input_integer != 4)
	{
		return(DCD_BADFORMAT);
	}

	if (*NAMNF != 0)
	{
		(*FREEINDEXES) = (int *) calloc(((*N)-(*NAMNF)), sizeof(int));

		if (*FREEINDEXES == NULL)
			return(DCD_BADMALLOC);
	
		/*  Read in an size				*/
		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading size of index array");
		CHECK_FEOF(ret_val, "reading size of index array");

		if (input_integer != ((*N)-(*NAMNF))*4)
		{
			return(DCD_BADFORMAT);
		}
		
		ret_val = read(fd, (*FREEINDEXES), ((*N)-(*NAMNF))*sizeof(int));

		CHECK_FREAD(ret_val, "reading size of index array");
		CHECK_FEOF(ret_val, "reading size of index array");

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading size of index array");
		CHECK_FEOF(ret_val, "reading size of index array");

		if (input_integer != ((*N)-(*NAMNF))*4)
		{
			return(DCD_BADFORMAT);
		}
	}

	return(0);
}

/************************************************************************/
/*									*/
/*			FUNCTION read_dcdstep				*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor to use					*/
/*	N - Number of atoms						*/
/*	X - array of X coordinates					*/
/*	Y - array of Y coordinates					*/
/*	Z - array of Z coordinates					*/
/*	num_fixed - Number of fixed atoms				*/
/*	first - flag 1->first time called				*/
/*	indexes - indexes for free atoms				*/
/*									*/
/*   OUTPUTS:								*/
/*	read step populates the arrays X, Y, Z and returns a 0 if the   */
/*   read is successful.  If the read fails, a negative error code is   */
/*   returned.								*/
/*									*/
/*	read_step reads in the coordinates from one step.  It places    */
/*   these coordinates into the arrays X, Y, and Z.			*/
/*									*/
/************************************************************************/

int read_dcdstep(int fd, int N, float *X, float *Y, float *Z, int num_fixed,
		 int first, int *indexes)

{
	int ret_val;		/*  Return value from read		*/
	int input_integer;	/*  Integer buffer space		*/
	int i;			/*  Loop counter			*/
	static float *tmpX;

	if (first && num_fixed)
	{
		tmpX = (float *) calloc(N-num_fixed, sizeof(float));

		if (tmpX==NULL)
		{
			return(DCD_BADMALLOC);
		}
	}

	/*  Get the first size from the file				*/
	ret_val = read(fd, &input_integer, sizeof(int));

	CHECK_FREAD(ret_val, "reading number of atoms at begining of step");

	/*  See if we've reached the end of the file			*/
	if (ret_val == 0)
	{
		free(tmpX);

		return(-1);
	}

	if ( (num_fixed==0) || first)
	{
		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, X, 4*N);
	
		CHECK_FREAD(ret_val, "reading X array");
		CHECK_FEOF(ret_val, "reading X array");

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading number of atoms after X array");
		CHECK_FEOF(ret_val, "reading number of atoms after X array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading number of atoms after X array");
		CHECK_FEOF(ret_val, "reading number of atoms after X array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, Y, 4*N);

		CHECK_FREAD(ret_val, "reading Y array");
		CHECK_FEOF(ret_val, "reading Y array");

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading number of atoms after Y array");
		CHECK_FEOF(ret_val, "reading number of atoms after Y array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading number of atoms after Y array");
		CHECK_FEOF(ret_val, "reading number of atoms after Y array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, Z, 4*N);

		CHECK_FREAD(ret_val, "reading Z array");
		CHECK_FEOF(ret_val, "reading Z array");

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading number of atoms after Z array");
		CHECK_FEOF(ret_val, "reading number of atoms after Z array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}
	}
	else
	{
		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, tmpX, 4*(N-num_fixed));
	
		CHECK_FREAD(ret_val, "reading Xtmp array");
		CHECK_FEOF(ret_val, "reading Xtmp array");

		for (i=0; i<N-num_fixed; i++)
		{
			X[indexes[i]-1]=tmpX[i];
		}

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading number of atoms after X array");
		CHECK_FEOF(ret_val, "reading number of atoms after X array");

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, &input_integer, sizeof(int));

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, tmpX, 4*(N-num_fixed));
	
		CHECK_FREAD(ret_val, "reading Xtmp array");
		CHECK_FEOF(ret_val, "reading Xtmp array");

		for (i=0; i<N-num_fixed; i++)
		{
			Y[indexes[i]-1]=tmpX[i];
		}

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading number of atoms after Y array");
		CHECK_FEOF(ret_val, "reading number of atoms after Y array");

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, &input_integer, sizeof(int));

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = read(fd, tmpX, 4*(N-num_fixed));
	
		CHECK_FREAD(ret_val, "reading Xtmp array");
		CHECK_FEOF(ret_val, "reading Xtmp array");

		for (i=0; i<N-num_fixed; i++)
		{
			Z[indexes[i]-1]=tmpX[i];
		}

		ret_val = read(fd, &input_integer, sizeof(int));

		CHECK_FREAD(ret_val, "reading number of atoms after Z array");
		CHECK_FEOF(ret_val, "reading number of atoms after Z array");

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}
	}

	return(0);
}


/*********************************************************************/
/*								     */
/*			FUNCTION open_dcd_write			     */
/*								     */
/*   INPUTS:							     */
/*	dcdfile - Name of the dcd file				     */
/*								     */
/*   OUTPUTS:							     */
/*	returns an open file descriptor for writing		     */
/*								     */
/*	This function will open a dcd file for writing.  It takes    */
/*   the filename to open as its only argument.	 It will return a    */
/*   valid file descriptor if successful or DCD_OPENFAILED if the    */
/*   open fails for some reason.  If the file specifed already       */
/*   exists, DCD_FILEEXISTS is returned.			     */
/*								     */
/*********************************************************************/

int open_dcd_write(char *dcdname)

{
	struct stat sbuf;
	int dcdfd;

	if (stat(dcdname, &sbuf) == 0) 
	{
		return(DCD_FILEEXISTS);
	} 
	else
	{
		if ( (dcdfd = open(dcdname, O_WRONLY|O_CREAT|O_EXCL,
					S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0)
		{
			return(DCD_OPENFAILED);
		}
	}

	return(dcdfd);
}

/************************************************************************/
/*									*/
/*				FUNCTION write_dcdstep			*/
/*									*/
/*   INPUTS:								*/
/*	fd - file descriptor for the DCD file to write to		*/
/*	N - Number of atoms						*/
/*	X - X coordinates						*/
/*	Y - Y coordinates						*/
/*	Z - Z coordinates						*/
/*									*/
/*   OUTPUTS:								*/
/*	none								*/
/*									*/
/*	write_dcdstep writes the coordinates out for a given timestep   */
/*   to the specified DCD file.						*/
/*									*/
/************************************************************************/

int write_dcdstep(int fd, int N, float *X, float *Y, float *Z)

{
	int out_integer;

	out_integer = N*4;
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) X, out_integer);
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) Y, out_integer);
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) Z, out_integer);
	write(fd, (char *) &out_integer, sizeof(int));
	
	return 1;
}

/*****************************************************************************/
/*									     */
/*				FUNCTION write_dcdheader		     */
/*									     */
/*   INPUTS:								     */
/*	fd - file descriptor for the dcd file				     */
/*	N - Number of atoms						     */
/*	NSET - Number of sets of coordinates				     */
/*	ISTART - Starting timestep of DCD file				     */
/*	NSAVC - Timesteps between DCD saves				     */
/*	DELTA - length of a timestep					     */
/*									     */
/*   OUTPUTS:								     */
/*	none								     */
/*									     */
/*	This function prints the "header" information to the DCD file.  Since*/
/*   this is duplicating an unformatted binary output from FORTRAN, its ugly.*/
/*   So if you're squimish, don't look.					     */
/*									     */
/*****************************************************************************/

int write_dcdheader(int fd, char *filename, int N, int NSET, int ISTART, 
		   int NSAVC, double DELTA)
{
	int	out_integer;
	char	title_string[200];
	int	user_id;
	struct  passwd *pwbuf;
	time_t 	cur_time;
	struct  tm *tmbuf;
	char    time_str[11];

	out_integer = 84;
	write(fd, (char *) & out_integer, sizeof(int));
	strcpy(title_string, "CORD");
	write(fd, title_string, 4);
	write(fd, (char *) &NSET, sizeof(int));
	write(fd, (char *) &ISTART, sizeof(int));
	write(fd, (char *) &NSAVC, sizeof(int));
	out_integer=0;
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &DELTA, sizeof(double));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	write(fd, (char *) &out_integer, sizeof(int));
	out_integer = 84;
	write(fd, (char *) & out_integer, sizeof(int));

	out_integer = 164;
	write(fd, (char *) & out_integer, sizeof(int));
	out_integer = 2;
	write(fd, (char *) & out_integer, sizeof(int));

	sprintf(title_string, "REMARKS FILENAME=%s CREATED BY VMD", filename);
	pad(title_string, 80);
	write(fd, title_string, 80);

	user_id=(int)getuid();
	pwbuf=getpwuid(user_id);
	cur_time=time(NULL);
	tmbuf=localtime(&cur_time);
	strftime(time_str, 10, "%m/%d/%y", tmbuf);

	sprintf(title_string, "REMARKS DATE: %s CREATED BY USER: %s",
	   time_str, pwbuf->pw_name);
	pad(title_string, 80);
	write(fd, title_string, 80);
	out_integer = 164;
	write(fd, (char *) & out_integer, sizeof(int));
	out_integer = 4;
	write(fd, (char *) & out_integer, sizeof(int));
	out_integer = N;
	write(fd, (char *) & out_integer, sizeof(int));
	out_integer = 4;
	write(fd, (char *) & out_integer, sizeof(int));
	
	return 1;
}

/****************************************************************/
/*								*/
/*			FUNCTION close_dcd_read			*/
/*								*/
/*   INPUTS:							*/
/*	fd - file descriptor to close				*/
/*	num_fixed - the number of fixed atoms			*/
/*	indexes - Array of indexes to be deallocated		*/
/*								*/
/*   OUTPUTS:							*/
/*	the file pointed to by fd is closed and the memory      */
/*   pointed to by indexes is freed if any was allocated	*/
/*								*/
/*	close_dcd_read closes a dcd file that was opened for    */
/*   reading.  It also deallocated memory used for the indexes  */
/*   used for the free atom list, if there were fixed atoms.    */
/*								*/
/****************************************************************/

void close_dcd_read(int fd, int num_fixed, int *indexes)

{
	close(fd);

	if (num_fixed)
	{
		free(indexes);
	}
}

/****************************************************************/
/*								*/
/*			FUNCTION close_dcd_write		*/
/*								*/
/*   INPUTS:							*/
/*	fd - file descriptor to close				*/
/*								*/
/*   OUTPUTS:							*/
/*	the file pointed to by fd				*/
/*								*/
/*	close_dcd_write close a dcd file that was opened for    */
/*   writing							*/
/*								*/
/****************************************************************/

void close_dcd_write(int fd)

{
	close(fd);
}

