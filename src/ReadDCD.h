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
 *      $RCSfile: ReadDCD.h,v $
 *      $Author: billh $        $Locker:  $                $State: Exp $
 *      $Revision: 1.3 $      $Date: 1995/05/11 23:42:25 $
 *
 ***************************************************************************
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
#include <malloc.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <pwd.h>
#include <time.h>

/*  DEFINE ERROR CODES THAT MAY BE RETURNED BY DCD ROUTINES		*/
#define DCD_DNE		-2	/*  DCD file does not exist		*/
#define DCD_OPENFAILED	-3	/*  Open of DCD file failed		*/
#define DCD_BADREAD 	-4	/*  read call on DCD file failed	*/
#define DCD_BADEOF	-5	/*  premature EOF found in DCD file	*/
#define DCD_BADFORMAT	-6	/*  format of DCD file is wrong		*/
#define DCD_FILEEXISTS  -7	/*  output file already exists		*/
#define DCD_BADMALLOC   -8	/*  malloc failed			*/


/*			FUNCTION ALLUSIONS				*/
int open_dcd_read(const char *);      /*  Open a DCD file for reading 	*/
int read_dcdheader(int, int*, int*, int*, int*, double*, int*, int**);	
				/*  Read the DCD header			*/
int read_dcdstep(int, int, float*, float*, float*, int, int, int*);	
				/*  Read a timestep's values		*/
int open_dcd_write(char *);     /*  Open a DCD file for writing		*/
int write_dcdstep(int, int, float *, float *, float *);
				/*  Write out a timesteps values	*/
int write_dcdheader(int, char*, int, int, int, int, double);	
				/*  Write a dcd header			*/
void close_dcd_read(int, int, int *);
				/*  Close a dcd file open for reading   */
void close_dcd_write(int);	/*  Close a dcd file open for writing   */

#endif /* ! DCDLIB_H */

