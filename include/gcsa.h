#ifndef GCSA_H
#define GCSA_H

#include "mrisurf.h"
#include "mrishash.h"

#define GCSA_MAGIC 0xababcdcd


#define GIBBS_SURFACE_NEIGHBORHOOD   4
#define GIBBS_SURFACE_NEIGHBORS      GIBBS_SURFACE_NEIGHBORHOOD

typedef struct
{
  VECTOR  *v_means ;
  MATRIX  *m_cov ;
  float   prior ;
  float   *label_priors[GIBBS_SURFACE_NEIGHBORHOOD] ;
  int    *labels[GIBBS_SURFACE_NEIGHBORHOOD] ;
  short   nlabels[ GIBBS_SURFACE_NEIGHBORHOOD];
} GCS, SURFACE_GAUSSIAN_CLASSIFIER ;

typedef struct
{
  int     nlabels ;
  int     max_labels ;   /* amount allocated */
  int     *labels ;
  GCS     *gcs ;
  int     total_training ;  /* total # of times this node was was accessed */
} GCSA_NODE ;

typedef struct
{
  MRI_SURFACE      *mris ;
  GCSA_NODE        *nodes ;
  int              ninputs ;
  int              icno ;
  MRIS_HASH_TABLE  *mht ;
} GAUSSIAN_CLASSIFIER_SURFACE_ARRAY, GCSA ;


#if 1

GCSA  *GCSAalloc(int ninputs, int nvertices) ;
int   GCSAfree(GCSA **pgcas) ;
int   GCSAtrainMeans(GCSA *gcsa, MRI_SURFACE *mris) ;
int   GCSAtrainCovariances(GCSA *gcsa, MRI_SURFACE *mris);
int   GCSAnormalizeMeans(GCSA *gcsa) ;
int   GCSAnormalizeCovariances(GCSA *gcsa) ;
int   GCSAwrite(GCSA *gcsa, char *fname) ;
GCSA  *GCSAread(char *fname) ;
int   GCSAlabel(GCSA *gcsa, MRI_SURFACE *mris) ;

#else
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

#define MIN_PRIOR  0.5

#endif
#endif
