#ifndef GCA_H
#define GCA_H

#include "mri.h"
#include "transform.h"

/* the volume that the classifiers are distributed within */
#define DEFAULT_VOLUME_SIZE   256
#define MAX_GCA_LABELS        1000

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
  float  log_p ;      /* current log probability of this sample */
  int    x ;          /* image coordinates */
  int    y ;
  int    z ;   
} GCA_SAMPLE, GCAS ;

#define GIBBS_NEIGHBORHOOD   6
#define GIBBS_NEIGHBORS      GIBBS_NEIGHBORHOOD
#define MAX_NBR_LABELS       9
typedef struct
{
  float   mean ;
  float   var ;
  float   prior ;
  float   *label_priors[GIBBS_NEIGHBORHOOD] ;
  char    *labels[GIBBS_NEIGHBORHOOD] ;
  short   nlabels[ GIBBS_NEIGHBORHOOD];
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
  double   T1_mean ;
  double   PD_mean ;
  double   T2_mean ;
  double   T1_var ;
  double   PD_var ;
  double   T2_var ;
  int      label ;
  int      total_training ;
} GCA_TISSUE_PARMS ;

typedef struct
{
  float     spacing ;    /* inter-node spacing */
  int       width ;
  int       height ;
  int       depth ;
  GCA_NODE  ***nodes ;
  int       ninputs ;
  GCA_TISSUE_PARMS tissue_parms[MAX_GCA_LABELS] ;
} GAUSSIAN_CLASSIFIER_ARRAY, GCA ;


int  GCAunifyVariance(GCA *gca) ;
int  GCAvoxelToNode(GCA *gca, MRI *mri,
                    int xv, int yv, int zv, int *pxn,int *pyn,int *pzn);
GCA  *GCAalloc(int ninputs, float spacing, int width, int height, int depth) ;
int  GCAfree(GCA **pgca) ;
int  GCANfree(GCA_NODE *gcan) ;
int  GCAtrain(GCA *gca, MRI *mri_inputs, MRI *mri_labels, LTA *lta, 
              GCA *gca_prune) ;
int  GCAwrite(GCA *gca, char *fname) ;
GCA  *GCAread(char *fname) ;
int  GCAcompleteTraining(GCA *gca) ;
MRI  *GCAlabel(MRI *mri_src, GCA *gca, MRI *mri_dst, LTA *lta) ;
MRI  *GCAclassify(MRI *mri_src,GCA *gca,MRI *mri_dst,LTA *lta,int max_labels);
MRI  *GCAreclassifyUsingGibbsPriors(MRI *mri_inputs, GCA *gca, MRI *mri_dst,
                                    LTA *lta, int max_iter, MRI *mri_fixed) ;
GCA  *GCAreduce(GCA *gca_src) ;
int  GCAnodeToVoxel(GCA *gca, MRI *mri, int xn, int yn, int zn, int *pxv, 
                    int *pyv, int *pzv) ;
float GCAcomputeLogImageProbability(GCA *gca, MRI *mri_inputs, MRI *mri_labels,
                                    LTA *lta) ;
float  GCAcomputeLogSampleProbability(GCA *gca, GCA_SAMPLE *gcas, 
                                      MRI *mri_inputs,
                                      MATRIX *m_L,int nsamples);
int   GCArankSamples(GCA *gca, GCA_SAMPLE *gcas, int nsamples, 
                     int *ordered_indices) ;
MRI  *GCAanneal(MRI *mri_inputs, GCA *gca, MRI *mri_dst,LTA *lta, 
                int max_iter);
int    GCAsourceVoxelToNode(GCA *gca, MRI *mri, LTA *lta,
                                 int xv, int yv, int zv, 
                                 int *pxn, int *pyn, int *pzn) ;
int  GCAnodeToSourceVoxel(GCA *gca, MRI *mri, MATRIX *m_L, 
                          int xv, int yv, int zv,
                          int *pxn, int *pyn, int *pzn) ;
#if 0
int    GCAsampleStats(GCA *gca, MRI *mri, LTA *lta, int class, 
                      Real x, Real y, Real z, 
                      Real *pmean, Real *pvar, Real *pprior) ;
#endif



MRI  *GCAannealUnlikelyVoxels(MRI *mri_inputs, GCA *gca, MRI *mri_dst,
                              LTA *lta, int max_iter, MRI *mri_fixed) ;
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
MRI        *GCAlabelMri(GCA *gca, MRI *mri, int label, LTA *lta) ;
MRI        *GCAbuildMostLikelyVolume(GCA *gca, MRI *mri) ;
MRI  *GCAlabelProbabilities(MRI *mri_inputs, GCA *gca, MRI *mri_dst, LTA *lta);
MRI  *GCAcomputeProbabilities(MRI *mri_inputs, GCA *gca, MRI *mri_labels, 
                              MRI *mri_dst, LTA *lta);

MRI   *GCAconstrainLabelTopology(GCA *gca, MRI *mri_inputs, MRI *mri_src, 
                                 MRI *mri_dst, LTA *lta) ;
MRI   *GCAexpandLabelIntoWM(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                            MRI *mri_dst, LTA *lta,MRI *mri_fixed,
                            int target_label) ;
MRI   *GCAexpandVentricle(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                          MRI *mri_dst, LTA *lta, int target_label) ;
MRI   *GCAnormalizeSamples(MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas, 
                           int nsamples, LTA *lta) ;
float GCAlabelProbability(MRI *mri_src, GCA *gca, LTA *lta,
                          int x, int y, int z, int label) ;
MRI   *GCAmaxLikelihoodBorders(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                               MRI *mri_dst, LTA *lta, int max_iter,
                               float min_ratio) ;
int   GCAaccumulateTissueStatistics(GCA *gca, MRI *mri_T1, MRI *mri_PD, 
                                    MRI *mri_parc, LTA *lta) ;
int   GCAnormalizeTissueStatistics(GCA *gca) ;
char *cma_label_to_name(int label) ;
int  GCArenormalizeLabels(MRI *mri_in, MRI *mri_labeled, GCA *gca, LTA *lta) ;
MRI   *GCArelabel_cortical_gray_and_white(GCA *gca, MRI *mri_inputs, 
                                          MRI *mri_src, MRI *mri_dst,LTA *lta);

int GCAdump(GCA *gca, MRI *mri, int x, int y, int z, LTA *lta, FILE *fp, 
            int verbose) ;
int GCArenormalizeIntensities(GCA *gca, int *labels, float *intensities, 
                              int num) ;
double  GCAlabelMean(GCA *gca, int label) ;
int     GCAregularizeConditionalDensities(GCA *gca, float smooth) ;


#define MIN_PRIOR  0.5

#endif
