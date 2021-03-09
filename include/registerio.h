/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

int regio_read_register(const char *regfile, char **subject, float *inplaneres,
                        float *betplaneres, float *intensity,  MATRIX **R,
                        int *float2int);

int regio_print_register(FILE *fp, const char *subject, float inplaneres,
                         float betplaneres, float intensity, const MATRIX *R,
                         int float2int);

int regio_write_register(const char *regfile, const char *subject, float inplaneres,
                         float betplaneres, float intensity, const MATRIX *R,
                         int float2int);

int regio_read_mincxfm(const char *xfmfile, MATRIX **R, char **fileinfo);
int regio_write_mincxfm(const char *xfmfile, const MATRIX *R, const char *fileinfo);


int regio_read_xfm4(const char *xfmfile, MATRIX **R);
int regio_read_xfm(const char *xfmfile, MATRIX **R);
int regio_write_surfacexform_to_register_dat(const MATRIX *B, const char *fname, 
                                             const MRI_SURFACE *mris, const MRI *mri, 
                                             const char *subject, int float2int);
MATRIX *regio_read_surfacexform_from_register_dat(const char *fname, 
                                                  const MRI_SURFACE *mris, 
                                                  const MRI *mri, char **subject);

MATRIX *regio_read_registermat(const char *regfile);
char *regio_read_subject(const char *regfile);

#endif /*BF_H_INC*/
