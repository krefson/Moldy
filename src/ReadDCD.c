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
 *      $RCSfile: ReadDCD.c,v $
 *      $Author: keith $        $Locker:  $                $State: Exp $
 *      $Revision: 1.7 $      $Date: 1997/11/27 15:58:19 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * C routines to read and write binary DCD files (which use the goofy
 * FORTRAN UNFORMATTED format).  These routines are courtesy of
 * Mark Nelson (autographs available upon request and $1 tip).
 *
 * Converted to ANSI C for use with Moldy - Keith Refson, Nov 1997         
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

#define CHECK_FREAD(X, msg)  if (X<=0 && ferror(fp)) \
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

#define CHECK_FEOF(X, msg)  if (X<=0 && feof(fp)) \
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

int read_dcdheader(FILE *fp, int *N, int *NSET, int *ISTART, int *NSAVC, 
		   double *DELTA, int *NAMNF, int **FREEINDEXES)
{
	int input_integer;	/*  Integer buffer space	*/
	char bigbuf[256];	/*  A large string buffer	*/
	int ret_val;		/*  Return value from read	*/
	int i;			/*  Loop counter		*/
	char HDR[5];		/*  Header = "CORD"		*/
	int I;
	int NTITLE;

	/*  First thing in the file should be an 84		*/
	ret_val = fread( &input_integer, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading first int from dcd file");
	CHECK_FEOF(ret_val, "reading first int from dcd file");

	if (input_integer != 84)
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in the string "CORD"			*/
	ret_val = fread(HDR, 1, 4, fp);

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
	ret_val = fread(NSET, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading NSET");
	CHECK_FEOF(ret_val, "reading NSET");

	/*  Read in ISTART, the starting timestep	     */
	ret_val = fread(ISTART, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading ISTART");
	CHECK_FEOF(ret_val, "reading ISTART");

	/*  Read in NSAVC, the number of timesteps between   */
	/*  dcd saves					     */
	ret_val = fread(NSAVC, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading NSAVC");
	CHECK_FEOF(ret_val, "reading NSAVC");


	/*  Skip blank integers				     */
	for (i=0; i<5; i++)
	{
		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading I");
		CHECK_FEOF(ret_val, "reading I");
	}

	/*  Read NAMNF, the number of free atoms	     */
	ret_val = fread(NAMNF, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading NAMNF");
	CHECK_FEOF(ret_val, "reading NAMNF");

	/*  Read in the timestep, DELTA				*/
	ret_val = fread(DELTA, sizeof(double), 1, fp);

	CHECK_FREAD(ret_val, "reading DELTA");
	CHECK_FEOF(ret_val, "reading DELTA");

	/*  Skip blank integers					*/
	for (i=0; i<9; i++)
	{
		ret_val = fread(&I, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading I");
		CHECK_FEOF(ret_val, "reading I");
	}

	/*  Get the end size of the first block			*/
	ret_val = fread(&input_integer, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading second 84 from dcd file");
	CHECK_FEOF(ret_val, "reading second 84 from dcd file");

	if (input_integer != 84)
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in the size of the next block			*/
	ret_val = fread(&input_integer, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading size of title block");
	CHECK_FEOF(ret_val, "reading size of title block");

	if ( ((input_integer-4)%80) == 0)
	{
		/*  Read NTITLE, the number of 80 characeter    */
		/*  title strings there are			*/
		ret_val = fread(&NTITLE, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading NTITLE");
		CHECK_FEOF(ret_val, "reading NTITLE");

		for (i=0; i<NTITLE; i++)
		{
			ret_val = fread(bigbuf, 80, 1, fp);

			CHECK_FREAD(ret_val, "reading TITLE");
			CHECK_FEOF(ret_val, "reading TITLE");
		}

		/*  Get the ending size for this block		*/
		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading size of title block");
		CHECK_FEOF(ret_val, "reading size of title block");
	}
	else
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in an 4				*/
	ret_val = fread(&input_integer, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading an 4");
	CHECK_FEOF(ret_val, "reading an 4");

	if (input_integer != 4)
	{
		return(DCD_BADFORMAT);
	}

	/*  Read in the number of atoms			*/
	ret_val = fread(N, sizeof(int), 1, fp);

	CHECK_FREAD(ret_val, "reading number of atoms");
	CHECK_FEOF(ret_val, "reading number of atoms");

	/*  Read in an 4				*/
	ret_val = fread(&input_integer, sizeof(int), 1, fp);

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
		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading size of index array");
		CHECK_FEOF(ret_val, "reading size of index array");

		if (input_integer != ((*N)-(*NAMNF))*sizeof(int))
		{
			return(DCD_BADFORMAT);
		}
		
		ret_val = fread((*FREEINDEXES), sizeof(int), (*N)-(*NAMNF), fp);

		CHECK_FREAD(ret_val, "reading size of index array");
		CHECK_FEOF(ret_val, "reading size of index array");

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading size of index array");
		CHECK_FEOF(ret_val, "reading size of index array");

		if (input_integer != ((*N)-(*NAMNF))*sizeof(int))
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

int read_dcdstep(FILE *fp, int N, float *X, float *Y, float *Z, 
		 int num_fixed, int first, int *indexes)
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
	ret_val = fread(&input_integer, sizeof(int), 1, fp);

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

		ret_val = fread(X, N, 4, fp);
	
		CHECK_FREAD(ret_val, "reading X array");
		CHECK_FEOF(ret_val, "reading X array");

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading number of atoms after X array");
		CHECK_FEOF(ret_val, "reading number of atoms after X array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading number of atoms after X array");
		CHECK_FEOF(ret_val, "reading number of atoms after X array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(Y, N, 4, fp);

		CHECK_FREAD(ret_val, "reading Y array");
		CHECK_FEOF(ret_val, "reading Y array");

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading number of atoms after Y array");
		CHECK_FEOF(ret_val, "reading number of atoms after Y array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading number of atoms after Y array");
		CHECK_FEOF(ret_val, "reading number of atoms after Y array");

		if (input_integer != 4*N)
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(Z, N, 4, fp);

		CHECK_FREAD(ret_val, "reading Z array");
		CHECK_FEOF(ret_val, "reading Z array");

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

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

		ret_val = fread(tmpX, N-num_fixed, 4, fp);
	
		CHECK_FREAD(ret_val, "reading Xtmp array");
		CHECK_FEOF(ret_val, "reading Xtmp array");

		for (i=0; i<N-num_fixed; i++)
		{
			X[indexes[i]-1]=tmpX[i];
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading number of atoms after X array");
		CHECK_FEOF(ret_val, "reading number of atoms after X array");

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(tmpX, N-num_fixed, 4, fp);
	
		CHECK_FREAD(ret_val, "reading Xtmp array");
		CHECK_FEOF(ret_val, "reading Xtmp array");

		for (i=0; i<N-num_fixed; i++)
		{
			Y[indexes[i]-1]=tmpX[i];
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading number of atoms after Y array");
		CHECK_FEOF(ret_val, "reading number of atoms after Y array");

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}

		ret_val = fread(tmpX, N-num_fixed, 4, fp);
	
		CHECK_FREAD(ret_val, "reading Xtmp array");
		CHECK_FEOF(ret_val, "reading Xtmp array");

		for (i=0; i<N-num_fixed; i++)
		{
			Z[indexes[i]-1]=tmpX[i];
		}

		ret_val = fread(&input_integer, sizeof(int), 1, fp);

		CHECK_FREAD(ret_val, "reading number of atoms after Z array");
		CHECK_FEOF(ret_val, "reading number of atoms after Z array");

		if (input_integer != 4*(N-num_fixed))
		{
			return(DCD_BADFORMAT);
		}
	}

	return(0);
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

int write_dcdstep(FILE *fp, int N, float *X, float *Y, float *Z)
{
	int out_integer;

	out_integer = N*4;
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) X, out_integer, 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) Y, out_integer, 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) Z, out_integer, 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	
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

int write_dcdheader(FILE *fp, char *filename, int N, int NSET, int ISTART, int NSAVC, double DELTA)
{
	int	out_integer;
	char	title_string[200];
	time_t 	cur_time;
	struct  tm *tmbuf;
	char    time_str[11];

	out_integer = 84;
	fwrite((char *) & out_integer, sizeof(int), 1, fp);
	strcpy(title_string, "CORD");
	fwrite(title_string, 4, 1, fp);
	fwrite((char *) &NSET, sizeof(int), 1, fp);
	fwrite((char *) &ISTART, sizeof(int), 1, fp);
	fwrite((char *) &NSAVC, sizeof(int), 1, fp);
	out_integer=0;
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &DELTA, sizeof(double), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	fwrite((char *) &out_integer, sizeof(int), 1, fp);
	out_integer = 84;
	fwrite((char *) & out_integer, sizeof(int), 1, fp);

	out_integer = 164;
	fwrite((char *) & out_integer, sizeof(int), 1, fp);
	out_integer = 2;
	fwrite((char *) & out_integer, sizeof(int), 1, fp);

	sprintf(title_string, "REMARKS FILENAME=%s CREATED BY VMD", filename);
	pad(title_string, 80);
	fwrite(title_string, 80, 1, fp);

	cur_time=time(NULL);
	tmbuf=localtime(&cur_time);
	strftime(time_str, 10, "%m/%d/%y", tmbuf);

	sprintf(title_string, "REMARKS DATE: %s CREATED BY USER: %s",
	   time_str, "mdshak");
	pad(title_string, 80);
	fwrite(title_string, 80, 1, fp);
	out_integer = 164;
	fwrite((char *) & out_integer, sizeof(int), 1, fp);
	out_integer = 4;
	fwrite((char *) & out_integer, sizeof(int), 1, fp);
	out_integer = N;
	fwrite((char *) & out_integer, sizeof(int), 1, fp);
	out_integer = 4;
	fwrite((char *) & out_integer, sizeof(int), 1, fp);
	
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

void close_dcd_read(FILE *fp, int num_fixed, int *indexes)
{
	fclose(fp);

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

void close_dcd_write(FILE *fp)
{
	fclose(fp);
}

