/**
 * @file  unwarpGradientNonlinearity.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.4 $
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



#ifndef unwarpGradientNonlinearity_H
#define unwarpGradientNonlinearity_H

#include"mri.h"

MRI *unwarpGradientNonlinearity(MRI *mri,
                                char *unwarp_gradientType,
                                char *unwarp_partialUnwarp,
                                char *unwarp_jacobianCorrection,
                                char *unwarp_interpType,
                                int unwarp_sincInterpHW);

int uGN_loadGradientData(char *unwarp_gradientType,
                         MATRIX *M_XYZ_2_beadIJK,
                         float **p_bead_dX,
                         float **p_bead_dY,
                         float **p_bead_dZ,
                         int *p_maxBeadI,
                         int *p_maxBeadJ,
                         int *p_maxBeadK);

int uGN_linInterp(float *bead_dX, float *bead_dY, float *bead_dZ,
                  float beadI, float beadJ, float beadK,
                  int maxBeadI, int maxBeadJ, int maxBeadK,
                  float *p_voxel_dX, float *p_voxel_dY, float *p_voxel_dZ);



#ifdef unwarpGradientNonlinearity_SRC
int unwarp_flag = 0;
char unwarp_gradientType[STRLEN];
char unwarp_partialUnwarp[STRLEN];
char unwarp_jacobianCorrection[STRLEN];
char unwarp_interpType[STRLEN];
int unwarp_sincInterpHW = 0;

#else
extern int unwarp_flag;
extern char unwarp_gradientType[STRLEN];
extern char unwarp_partialUnwarp[STRLEN];
extern char unwarp_jacobianCorrection[STRLEN];
extern char unwarp_interpType[STRLEN];
extern int unwarp_sincInterpHW;

#endif /* #ifdef unwarpGradientNonlinearity_SRC */

#endif /* #ifndef unwarpGradientNonlinearity_H */
