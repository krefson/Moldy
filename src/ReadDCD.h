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
 *      $Author: keith $        $Locker:  $                $State: Exp $
 *      $Revision: 1.4 $      $Date: 1997/11/27 15:58:19 $
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
int read_dcdheader();	
				/*  Read the DCD header			*/
int read_dcdstep();	
				/*  Read a timestep's values		*/
int write_dcdstep();
				/*  Write out a timesteps values	*/
int write_dcdheader();	
				/*  Write a dcd header			*/
void close_dcd_read();
				/*  Close a dcd file open for reading   */
void close_dcd_write();	/*  Close a dcd file open for writing   */

#endif /* ! DCDLIB_H */

