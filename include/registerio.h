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
 *    $Date: 2007/12/21 04:30:08 $
 *    $Revision: 1.6 $
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

#include "mri.h"
#include "mrisurf.h"
int
regio_write_surfacexform_to_register_dat(MATRIX *B, char *fname, 
                                         MRI_SURFACE *mris, MRI *mri, 
                                         char *subject, int float2int);
MATRIX *
regio_read_surfacexform_from_register_dat(char *fname, MRI_SURFACE *mris, 
                                          MRI *mri, char **subject);


#endif /*BF_H_INC*/
