#ifndef MRI_CLASS_H
#define MRI_CLASS_H

#include "gclass.h"

typedef struct
{
  int       scale ;     /* size reduction from original full-res. image */
  int       width ;
  int       height ;
  int       depth ;
  int       nvars ;
  GCLASSIFY ****gcs ;
} MRI_CLASSIFIER, MRIC ;


MRIC    *MRIclassAlloc(int width, int height, int depth, int scale, int nvars);
int     MRIclassFree(MRIC **pmric) ;
int     MRIclassTrain(MRIC *mric, MRI *mri_src,MRI *mri_norm,MRI *mri_target);
MRI     *MRIclassify(MRIC *mric, MRI *mri_src, MRI *mri_dst, MRI *mri_norm,
                     float conf) ;

#endif
