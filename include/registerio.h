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
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.9 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


/**********************************************************
Routines for handling register.dat files
***********************************************************/


#ifndef REGISTERIO_H_INC
#define REGISTERIO_H_INC

#include <stdio.h>
#include "matrix.h"
#include "mri.h"
#include "mrisurf.h"

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
int regio_write_surfacexform_to_register_dat(MATRIX *B, char *fname, 
                                             MRI_SURFACE *mris, MRI *mri, 
                                             char *subject, int float2int);
MATRIX *regio_read_surfacexform_from_register_dat(char *fname, 
                                                  MRI_SURFACE *mris, 
                                                  MRI *mri, char **subject);

int
regio_write_surfacexform_to_register_dat(MATRIX *B, char *fname, 
                                         MRI_SURFACE *mris, MRI *mri, 
                                         char *subject, int float2int);
MATRIX *
regio_read_surfacexform_from_register_dat(char *fname, MRI_SURFACE *mris, 
                                          MRI *mri, char **subject);

MATRIX *regio_read_registermat(char *regfile);

#endif /*BF_H_INC*/
