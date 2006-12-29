/**
 * @file  gcarray.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.4 $
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


#ifndef GCARRAY_H
#define GCARRAY_H

#include "gclass.h"
#include "mri.h"

/*
  the classifiers are spaced so that there is scale/2 padding at each
  border, then a classifier center every scale pixels.
  */
typedef struct
{
  int       scale ;     /* size reduction from original full-res. image */
  int       width ;
  int       height ;
  int       depth ;
  int       swidth ;    /* dimensions of source MRI */
  int       sheight ;
  int       sdepth ;
  int       nvars ;     /* # of inputs */
  GCLASSIFY ****gcs ;
  Transform *transform ;
  Transform *inverse_transform ;
  Real      xstart ;     /* coord. of (0,0,0) in Talairach space */
  Real      ystart ;
  Real      zstart ;
}
GCARRAY ;


GCARRAY    *GCarrayAlloc(MRI *mri_template, int scale, int nvars);
GCARRAY    *GCarrayTrainAll(GCARRAY *gcarray, char *training_file_name,
                            int scale, int nvars) ;
int     GCarraySetTransform(GCARRAY *gcarray, Transform *transform,
                            Transform *inverse_transform) ;
int     GCarrayFree(GCARRAY **pgcarray) ;
int     GCarrayTrain(GCARRAY *gcarray, MRI *mri_src,MRI *mri_norm,
                     MRI *mri_target);
int     GCarrayUpdate(GCARRAY *gcarray, MRI *mri_src,MRI *mri_norm,
                      MRI *mri_target);
int     GCarrayFinish(GCARRAY *gcarray) ;
MRI     *GCarrayClassify(GCARRAY *gcarray, MRI *mri_src, MRI *mri_dst,
                         float conf, MRI *mri_probs, MRI *mri_classes) ;
int     GCarrayToVoxel(GCARRAY *gcarray, int xc, int yc, int zc,
                       int *pxv, int *pyv, int *pzv) ;
int     GCarrayVoxelToClass(GCARRAY *gcarray, int xv, int yv, int zv,
                            int *pxc, int *pyc, int *pzc) ;
GCARRAY    *GCarrayRead(char *fname) ;
int     GCarrayWrite(GCARRAY *gcarray, char *fname) ;
MRI     *GCarrayThreshold(GCARRAY *gcarray, MRI *mri_probs, MRI *mri_classes,
                          MRI *mri_dst, float threshold) ;

GCLASSIFY *MRIgetClassifier(GCARRAY *gcarray, MRI *mri, int xv, int yv,int zv);
int      GCarrayUpdateMeans(GCARRAY *gcarray, MRI *mris[], MRI *mri_target,
                            int nimages) ;
int      GCarrayComputeMeans(GCARRAY *gcarray) ;
int      GCarrayComputeCovariances(GCARRAY *gcarray) ;
int      GCarrayUpdateCovariances(GCARRAY *gcarray, MRI *mris[],
                                  MRI *mri_target, int nimages) ;

#endif
