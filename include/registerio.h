/**
 * @file  registerio.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


/**********************************************************
Routines for handling register.dat files
***********************************************************/


#ifndef REGISTERIO_H_INC
#define REGISTERIO_H_INC

#include <stdio.h>
#include "matrix.h"

int regio_read_register(char *regfile, char **subject, float *inplaneres,
                        float *betplaneres, float *intensity,  MATRIX **R,
                        int *float2int);

int regio_print_register(FILE *fp, char *subject, float inplaneres,
                         float betplaneres, float intensity, MATRIX *R,
                         int float2int);

int regio_write_register(char *regfile, char *subject, float inplaneres,
                         float betplaneres, float intensity, MATRIX *R,
                         int float2int);

int regio_read_mincxfm(char *xfmfile, MATRIX **R, char **fileinfo);
int regio_write_mincxfm(char *xfmfile, MATRIX *R, char *fileinfo);


int regio_read_xfm4(char *xfmfile, MATRIX **R);
int regio_read_xfm(char *xfmfile, MATRIX **R);


#endif /*BF_H_INC*/
