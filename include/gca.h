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
  int    xn ;         /* node coordinates */
  int    yn ;
  int    zn ;
  int    label ;
  int    n ;          /* index in gcan structure */
  float  prior ;
  float  var ;
  float  mean ;
} GCA_SAMPLE, GCAS ;

typedef struct
{
  float   mean ;
  float   var ;
  float   prior ;
} GC1D, GAUSSIAN_CLASSIFIER_1D ;

typedef struct
{
  int  nlabels ;
  int  max_labels ;   /* amount allocated */
  char *labels ;
  GC1D *gcs ;
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


int  GCAvoxelToNode(GCA *gca, MRI *mri,
                    int xv, int yv, int zv, int *pxn,int *pyn,int *pzn);
GCA  *GCAalloc(int ninputs, float spacing, int width, int height, int depth) ;
int  GCAfree(GCA **pgca) ;
int  GCANfree(GCA_NODE *gcan) ;
int  GCAtrain(GCA *gca, MRI *mri_inputs, MRI *mri_labels, LTA *lta) ;
int  GCAwrite(GCA *gca, char *fname) ;
GCA  *GCAread(char *fname) ;
int  GCAcompleteTraining(GCA *gca) ;
MRI  *GCAlabel(MRI *mri_src, GCA *gca, MRI *mri_dst, LTA *lta) ;
MRI  *GCAclassify(MRI *mri_src,GCA *gca,MRI *mri_dst,LTA *lta,int max_labels);
GCA  *GCAreduce(GCA *gca_src) ;
int  GCAnodeToVoxel(GCA *gca, MRI *mri, int xn, int yn, int zn, int *pxv, 
                    int *pyv, int *pzv) ;
float GCAcomputeLogImageProbability(GCA *gca, MRI *mri_inputs, MRI *mri_labels,
                                    LTA *lta) ;
float  GCAcomputeLogSampleProbability(GCA *gca, GCA_SAMPLE *gcas, 
                                      MRI *mri_inputs,
                                      MATRIX *m_L,int nsamples);
int    GCAsourceVoxelToNodePoint(GCA *gca, MRI *mri, LTA *lta,
                                 Real xv, Real yv, Real zv, 
                                 Real *pxn, Real *pyn, Real *pzn) ;
#if 0
int    GCAsampleStats(GCA *gca, MRI *mri, LTA *lta, int class, 
                      Real x, Real y, Real z, 
                      Real *pmean, Real *pvar, Real *pprior) ;
#endif



GCA_SAMPLE *GCAfindContrastSamples(GCA *gca, int *pnsamples, int min_spacing,
                                 float min_prior) ;
GCA_SAMPLE *GCAfindStableSamples(GCA *gca, int *pnsamples, int min_spacing,
                                 float min_prior) ;
GCA_SAMPLE *GCAfindStableSamplesByLabel(GCA *gca, int nsamples, 
                                        float min_prior) ;
int       GCAtransformSamples(GCA *gca_src, GCA *gca_dst, GCA_SAMPLE *gcas, 
                              int nsamples) ;
int        GCAwriteSamples(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, int nsamples,
                           char *fname) ;
int        GCAtransformAndWriteSamples(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, 
                                       int nsamples,char *fname,LTA *lta) ;
MRI        *GCAmri(GCA *gca, MRI *mri) ;

#define MIN_PRIOR  0.5

#endif
