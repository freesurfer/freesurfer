#ifndef MRICLASS_H
#define MRICLASS_H

#include "classify.h"
#include "backprop.h"
#include "gclass.h"
#include "artmap.h"
#include "box.h"
#include "mri.h"

#define GAUSSIAN_NCLASSES    4
#define BACKGROUND           0
#define GREY_MATTER          1
#define WHITE_MATTER         2
#define BRIGHT_MATTER        3
#define LO_LIM               70
#define HI_LIM               150

#define MAX_INPUTS           4


/* bitfield of feature types */
#define FEATURE_INTENSITY      0x00001
#define FEATURE_ZSCORE3        0x00002
#define FEATURE_ZSCORE5        0x00004
#define FEATURE_DIRECTION      0x00010
#define FEATURE_MEAN3          0x00020
#define FEATURE_MEAN5          0x00040
#define FEATURE_CPOLV_MEAN3    0x00080
#define FEATURE_CPOLV_MEAN5    0x00100
#define FEATURE_CPOLV_MEDIAN3  0x00200
#define FEATURE_CPOLV_MEDIAN5  0x00400
#define FEATURE_CPOLV          (FEATURE_CPOLV_MEAN3 | FEATURE_CPOLV_MEAN5 | \
                                FEATURE_CPOLV_MEDIAN3 | FEATURE_CPOLV_MEDIAN5)
#define MAX_FEATURE            FEATURE_CPOLV_MEDIAN5

typedef struct
{
  int         features ;    /* bit field of above features */
  int         type ;        /* what kind of classifier are we using */
  int         ninputs ;     /* # of inputs to the classifier */
  void        *parms ;      /* classifier specific parameters */
  GCLASSIFY   *gc ;         /* if type == CLASSIFIER_GAUSSIAN */
  BACKPROP    *bp ;         /* if type == CLASSIFIER_BACKPROP */
  ARTMAP      *artmap ;     /* if type == CLASSIFIER_ARTMAP */
  MRI         *mri_priors ; /* prior probabilities */
  char        prior_fname[100] ;
} MRI_CLASSIFIER, MRIC ;

MRIC   *MRICalloc(int type, int features, void *parms) ;
int    MRICfree(MRIC **pmri) ;
int    MRICtrain(MRIC *mric, char *file_name, char *prior_fname) ;
/*int    MRICclassify(MRIC *mric, MRI *mri_src, float *pprob) ;*/
MRIC   *MRICread(char *fname) ;
int    MRICwrite(MRIC *mric, char *fname) ;
MRI    *MRICclassify(MRIC *mric, MRI *mri_src, MRI *mri_dst, float conf, 
                     MRI *mri_probs, MRI *mri_classes) ;
int    MRICupdateMeans(MRIC *mric, MRI *mri_src, MRI *mri_target, BOX *bbox) ;
int    MRICcomputeMeans(MRIC *mric) ;  
int    MRICupdateCovariances(MRIC *mric, MRI *mri_src, MRI *mri_target, 
                             BOX *bbox) ;
int    MRICcomputeCovariances(MRIC *mric) ;
int    MRICcomputeInputs(MRI *mri, int x,int y,int z,float *obs, int features);
MRI    *MRICbuildTargetImage(MRI *mri_src, MRI *mri_target, MRI *mri_wm,
                             int lo_lim, int hi_lim) ;
MRI    *MRICupdatePriors(MRI *mri_target, MRI *mri_priors, int scale) ;
int    MRInormalizePriors(MRI *mri_priors) ;

#endif
