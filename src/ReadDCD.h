/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 * C routines to read and write binary DCD files (which use the goofy
 * FORTRAN UNFORMATTED format).  These routines are courtesy of
 * Mark Nelson (autographs available upon request and $1 tip).
 *
 ***************************************************************************/
#ifndef READ_DCD_H
#define READ_DCD_H

#include <stdio.h>
#include "stdlib.h"
#include <errno.h>
#include "time.h"

/*  DEFINE ERROR CODES THAT MAY BE RETURNED BY DCD ROUTINES		*/
#define DCD_DNE		-2	/*  DCD file does not exist		*/
#define DCD_OPENFAILED	-3	/*  Open of DCD file failed		*/
#define DCD_BADREAD 	-4	/*  read call on DCD file failed	*/
#define DCD_BADEOF	-5	/*  premature EOF found in DCD file	*/
#define DCD_BADFORMAT	-6	/*  format of DCD file is wrong		*/
#define DCD_FILEEXISTS  -7	/*  output file already exists		*/
#define DCD_BADMALLOC   -8	/*  malloc failed			*/


/*			FUNCTION ALLUSIONS				*/
int read_dcdheader(FILE *fp, int *N, int *NSET, int *ISTART, int *NSAVC, double *DELTA, int *NAMNF, int **FREEINDEXES);	
				/*  Read the DCD header			*/
int read_dcdstep(FILE *fp, int N, float *X, float *Y, float *Z, int num_fixed, int first, int *indexes);	
				/*  Read a timestep's values		*/
int write_dcdstep(FILE *fp, int N, float *X, float *Y, float *Z);
				/*  Write out a timesteps values	*/
int write_dcdheader(FILE *fp, char *filename, int N, int NSET, int ISTART, int NSAVC, double DELTA);	
				/*  Write a dcd header			*/
void close_dcd_read(FILE *fp, int num_fixed, int *indexes);
				/*  Close a dcd file open for reading   */
void close_dcd_write(FILE *fp);	/*  Close a dcd file open for writing   */

#endif /* ! DCDLIB_H */

