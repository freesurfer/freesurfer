#ifndef MRI_CLASS_H
#define MRI_CLASS_H

#include "gclass.h"
#include "mri.h"
#include "volume_io.h"

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
  int       nvars ;
  GCLASSIFY ****gcs ;
  Transform *transform ;
  Transform *inverse_transform ;
  Real      xstart ;     /* coord. of (0,0,0) in Talairach space */
  Real      ystart ;
  Real      zstart ;
} MRI_CLASSIFIER, MRIC ;


MRIC    *MRIclassAlloc(int width, int height, int depth, int scale, int nvars);
int     MRIclassTrainAll(MRIC *mric, char *training_file_name) ;
int     MRIclassSetTransform(MRIC *mric, Transform *transform, 
                             Transform *inverse_transform) ;
int     MRIclassFree(MRIC **pmric) ;
int     MRIclassTrain(MRIC *mric, MRI *mri_src,MRI *mri_norm,MRI *mri_target);
int     MRIclassUpdate(MRIC *mric, MRI *mri_src,MRI *mri_norm,MRI *mri_target);
int     MRIclassFinish(MRIC *mric) ;
MRI     *MRIclassify(MRIC *mric, MRI *mri_src, MRI *mri_dst, MRI *mri_norm,
                     float conf, MRI *mri_probs, MRI *mri_classes) ;
int     MRIclassToVoxel(MRIC *mric, int xc, int yc, int zc,
                        int *pxv, int *pyv, int *pzv) ;
int     MRIvoxelToClass(MRIC *mric, int xv, int yv, int zv,
                        int *pxc, int *pyc, int *pzc) ;
MRIC    *MRIclassRead(char *fname) ;
int     MRIclassWrite(MRIC *mric, char *fname) ;
MRI     *MRIclassThreshold(MRIC *mric, MRI *mri_probs, MRI *mri_classes,
                           MRI *mri_dst, float threshold) ;

GCLASSIFY *MRIgetClassifier(MRIC *mric, MRI *mri, int xv, int yv, int zv) ;
int      MRIclassUpdateMeans(MRIC *mric, MRI *mris[], MRI *mri_target, 
                             int nimages) ;
int      MRIclassComputeMeans(MRIC *mric) ;
int      MRIclassComputeCovariances(MRIC *mric) ;
int      MRIclassUpdateCovariances(MRIC *mric, MRI *mris[], MRI *mri_target, 
                             int nimages) ;

#endif
