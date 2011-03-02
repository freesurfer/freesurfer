/**
 * @file aseg_edit_svm.c
 * @brief utilities for using SVMs to reclassify aseg voxels
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:42 $
 *    $Revision: 1.2 $
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

extern const char* Progname;

/*-----------------------------------------------------
  INCLUDE FILES
  -------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "error.h"
#include "mri.h"
#include "voxlist.h"
#include "diag.h"
#include "aseg_edit_svm.h"

#define NINPUTS(wsize, nscales)    (4+(wsize*wsize*wsize)*nscales*(1+1)) // 3+1 = 3 grad components plus image intensity
float **
build_svm_training_data(VOXEL_LIST *vl_total, int target_label, float *svm_classes, 
                        MRI **mri_smooth, MRI **mri_grad, MRI **mri_laplacian, MRI *mri_dtrans,
                        int wsize, int nscales)
{
  float **svm_inputs ;
  int   i, ninputs ;
  MRI   *mri_dtrans_grad ;
 
  svm_inputs = (float **)calloc(vl_total->nvox, sizeof(float *)) ;
  if (svm_inputs == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d SVM input vector", Progname, vl_total->nvox) ;

  ninputs = NINPUTS(wsize, nscales) ;
  mri_dtrans_grad = MRIsobel(mri_dtrans, NULL, NULL) ;
  for (i = 0 ; i < vl_total->nvox ; i++)
  {
    if (svm_classes)
      svm_classes[i] = ((int)vl_total->vdst[i] == target_label) ? 1 : -1 ;

    svm_inputs[i] = build_svm_inputs_at_voxel(vl_total, i, target_label,
                                              mri_smooth, mri_grad, mri_laplacian,
                                              mri_dtrans, mri_dtrans_grad, wsize, nscales, NULL) ;
    if (svm_inputs[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %dth SVM input vector size %d", 
                Progname, vl_total->nvox, ninputs) ;
  }
  MRIfree(&mri_dtrans_grad) ;
  return(svm_inputs) ;
}
float *
build_svm_inputs_at_voxel(VOXEL_LIST *vl, int i, int target_label, 
    MRI **mri_smooth, MRI **mri_grad, MRI **mri_laplacian, MRI *mri_dtrans, MRI *mri_dtrans_grad, int wsize, int nscales, float *svm_inputs)
{
  int   s, xk, yk, zk, xi, yi, zi, x, y, z, ninputs, whalf, input ;

  whalf = (wsize-1)/2 ;
  ninputs = NINPUTS(wsize, nscales) ; 
  if (svm_inputs == NULL)
  {
    svm_inputs = (float *)calloc(ninputs, sizeof(float)) ;
    if (svm_inputs == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d SVM input vector", Progname, ninputs) ;
  }

  x = vl->xi[i] ; y = vl->yi[i] ; z = vl->zi[i] ; 

  input = 0 ;
  svm_inputs[input++] = MRIgetVoxVal(mri_dtrans, x, y, z, 0) ;
  if (!finite(svm_inputs[input-1]))
    DiagBreak() ;
  svm_inputs[input++] = MRIgetVoxVal(mri_dtrans_grad, x, y, z, 0) ;
  if (!finite(svm_inputs[input-1]))
    DiagBreak() ;
  svm_inputs[input++] = MRIgetVoxVal(mri_dtrans_grad, x, y, z, 1) ;
  if (!finite(svm_inputs[input-1]))
    DiagBreak() ;
  svm_inputs[input++] = MRIgetVoxVal(mri_dtrans_grad, x, y, z, 2) ;
  if (!finite(svm_inputs[input-1]))
    DiagBreak() ;
  for (xk = -whalf ; xk <= whalf ; xk++)
    for (yk = -whalf ; yk <= whalf ; yk++)
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        for (s = 0 ; s < nscales ; s++)
        {
          xi = mri_smooth[s]->xi[x+xk] ; yi = mri_smooth[s]->yi[y+yk] ; 
          zi = mri_smooth[s]->zi[z+zk] ;
          svm_inputs[input++] = MRIgetVoxVal(mri_smooth[s], xi, yi, zi, 0) ;
          if (!finite(svm_inputs[input-1]))
            DiagBreak() ;
#if 0
          svm_inputs[input++] = MRIgetVoxVal(mri_grad[s], xi, yi, zi, 0) ;
          svm_inputs[input++] = MRIgetVoxVal(mri_grad[s], xi, yi, zi, 1) ;
          svm_inputs[input++] = MRIgetVoxVal(mri_grad[s], xi, yi, zi, 2) ;
#endif
          svm_inputs[input++] = MRIgetVoxVal(mri_laplacian[s], xi, yi, zi, 0) ;
          if (!finite(svm_inputs[input-1]))
            DiagBreak() ;
        }
      }

  for (input = 0 ; input < ninputs ; input++)
    if (!finite(svm_inputs[input]))
      DiagBreak() ;
  return(svm_inputs) ;
}
