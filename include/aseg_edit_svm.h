/**
 * @file aseg_edit_svm.h
 * @brief prototypes and structs for using SVMs to reclassify aseg voxels
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2008/02/07 00:55:32 $
 *    $Revision: 1.1 $
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

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#ifndef ASEG_EDIT_SVM_H
#define ASEG_EDIT_SVM_H

#include "mri.h"
#include "voxlist.h"

float **build_svm_training_data(VOXEL_LIST *vl_total, int target_label, 
                                float *svm_classes, MRI **mri_smooth, 
                                MRI **mri_grad, MRI **mri_laplacian, MRI *mri_dtrans, int wsize, int nscales);
float *build_svm_inputs_at_voxel(VOXEL_LIST *vl, int i, int target_label, 
                                 MRI **mri_smooth, MRI **mri_grad, MRI **mri_laplacian, MRI *mri_dtrans, 
                                 MRI *mri_dtrans_grad,
                                 int wsize, int nscales, float *svm_inputs) ;

#define NINPUTS(wsize, nscales)    (4+(wsize*wsize*wsize)*nscales*(1+1)) // 3+1 = 3 grad components plus image intensity
#endif
