#ifndef GCA_H
#define GCA_H

#include "mri.h"
#include "transform.h"

/* the volume that the classifiers are distributed within */
#define DEFAULT_VOLUME_SIZE   256

typedef struct
{
  float   spacing ;
  int     use_gradient ;
} GCA_PARMS ;

/*
  the classifiers are spaced so that there is scale/2 padding at each
  border, then a classifier center every scale pixels.
  */



typedef struct
{
  float   mean ;
  float   var ;
  float   prior ;
} GC, GAUSSIAN_CLASSIFIER_1D ;

typedef struct
{
  int  nlabels ;
  int  max_labels ;   /* amount allocated */
  char *labels ;
  GC   *gcs ;
  int  total_training ;  /* total # of times this node was was accessed */
} GCA_NODE ;

typedef struct
{
  float     spacing ;    /* inter-node spacing */
  int       width ;
  int       height ;
  int       depth ;
  GCA_NODE  ***nodes ;
  int       ninputs ;
} GAUSSIAN_CLASSIFIER_ARRAY, GCA ;


int  GCAvoxelToNode(GCA *gca, int xv, int yv, int zv, int *pxn,int *pyn,int *pzn);
GCA  *GCAalloc(int ninputs, float spacing, int length) ;
int  GCAfree(GCA **pgca) ;
int  GCANfree(GCA_NODE *gcan) ;
int  GCAtrain(GCA *gca, MRI *mri_inputs, MRI *mri_labels, LTA *lta) ;
int  GCAwrite(GCA *gca, char *fname) ;
GCA  *GCAread(char *fname) ;
int  GCAcompleteTraining(GCA *gca) ;
MRI  *GCAlabel(MRI *mri_src, GCA *gca, MRI *mri_dst, LTA *lta) ;

#endif
