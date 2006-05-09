#ifndef GCA_H
#define GCA_H

#include "mri.h"
#include "transform.h"
#define MIN_PRIOR  0.5
#define MAX_GCA_INPUTS 100
/* GCA types *************/
#define GCA_NORMAL     0  // standard way to create GCA
#define GCA_FLASH      1  // used flash data to create GCA
#define GCA_PARAM      2  // used T1 and PD data to create GCA
#define GCA_UNKNOWN   99  // what ???

#define FILE_TAG        0xab2c
#define TAG_PARAMETERS  0x0001
#define TAG_GCA_TYPE    0x0002
#define TAG_GCA_DIRCOS  0x0003

/* the volume that the classifiers are distributed within */
#define DEFAULT_VOLUME_SIZE   256
#define MAX_GCA_LABELS        100

/* for flags field */
#define GCA_NO_FLAGS          0x0000
#define GCA_NO_MRF            0x0001
#define GCA_XGRAD             0x0002
#define GCA_YGRAD             0x0004
#define GCA_ZGRAD             0x0008
#define GCA_GRAD              (GCA_XGRAD | GCA_YGRAD | GCA_ZGRAD)
#define GCA_NO_GCS            0x0010

typedef struct
{
  float   prior_spacing ;
  float   node_spacing ;
  int     use_gradient ;
} GCA_PARMS ;

/*
  the classifiers are spaced so that there is scale/2 padding at each
  border, then a classifier center every scale pixels.
  */

typedef struct
{
  int    xp ;         /* prior coordinates */
  int    yp ;
  int    zp ;
  int    label ;
  float  prior ;
  float  *covars ;
  float  *means ;
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
  float   *means ;
  float   *covars ;
  float   **label_priors ;
  unsigned char    **labels ;
  short   *nlabels;
  short   n_just_priors ;
  int     ntraining ;
  char    regularized ;
} GC1D, GAUSSIAN_CLASSIFIER_1D ;

typedef struct
{
  short nlabels ;
  short max_labels ;
  unsigned char  *labels ;
  float *priors ;
  int   total_training ;
} GCA_PRIOR ;

typedef struct
{
  int  nlabels ;
  int  max_labels ;   /* amount allocated */
  unsigned char *labels ;
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
  float     node_spacing ;    /* inter-node spacing */
  float     prior_spacing ;    /* inter-prior spacing */
  int       node_width ;
  int       node_height ;
  int       node_depth ;
  GCA_NODE  ***nodes ;
  int       prior_width ;
  int       prior_height ;
  int       prior_depth ;
  GCA_PRIOR ***priors ;
  int       ninputs ;
  GCA_TISSUE_PARMS tissue_parms[MAX_GCA_LABELS] ;
  int    flags ;
  int       type ;
  double    TRs[MAX_GCA_INPUTS] ;  /* TR of training data (in msec) */
  double    FAs[MAX_GCA_INPUTS] ;  /* flip angles of training data (in radians) */
  double    TEs[MAX_GCA_INPUTS] ;  /* TE  of training data (in msec) */
  // direction cosine info (common to both prior and node)
  float         x_r, x_a, x_s; 
  float         y_r, y_a, y_s; 
  float         z_r, z_a, z_s; 
  float         c_r, c_a, c_s; 
  int           width;            
  int           height;
  int           depth;
  float         xsize;
  float         ysize;
  float         zsize;
  MRI          *mri_node__;       // these three MRI are helper to get
  MRI          *mri_prior__;      // appropriate transforms
  MRI          *mri_tal__;
  MATRIX       *node_i_to_r__;
  MATRIX       *node_r_to_i__;
  MATRIX       *prior_i_to_r__;
  MATRIX       *prior_r_to_i__;
  MATRIX       *tal_i_to_r__;
  MATRIX       *tal_r_to_i__;
  MATRIX       *tmp__;
} GAUSSIAN_CLASSIFIER_ARRAY, GCA ;


int  GCAsetFlashParameters(GCA *gca, double *TRs, double *FAs, double *TEs) ;
int  GCAunifyVariance(GCA *gca) ;
int  GCAvoxelToPrior(GCA *gca, MRI *mri,
                    int xv, int yv, int zv, int *pxp,int *pyp,int *pzp);
int  GCAvoxelToNode(GCA *gca, MRI *mri,
                    int xv, int yv, int zv, int *pxn,int *pyn,int *pzn);
GCA  *GCAalloc(int ninputs, float prior_spacing, float node_spacing, int width, int height, int depth, int flags) ;
int  GCAfree(GCA **pgca) ;
int  GCANfree(GCA_NODE *gcan, int ninputs) ;
int  GCAtrain(GCA *gca, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform, 
              GCA *gca_prune, int noint) ;
int  GCAtrainCovariances(GCA *gca, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform) ;
int  GCAwrite(GCA *gca, char *fname) ;
GCA  *GCAread(char *fname) ;
int  GCAcompleteMeanTraining(GCA *gca) ;
int  GCAcompleteCovarianceTraining(GCA *gca) ;
MRI  *GCAlabel(MRI *mri_src, GCA *gca, MRI *mri_dst, TRANSFORM *transform) ;
MRI  *GCAclassify(MRI *mri_src,GCA *gca,MRI *mri_dst,TRANSFORM *transform,int max_labels);
MRI  *GCAreclassifyUsingGibbsPriors(MRI *mri_inputs, GCA *gca, MRI *mri_dst,
                                    TRANSFORM *transform, int max_iter, MRI *mri_fixed,
                                    int restart, void (*update_func)(MRI *));
GCA  *GCAreduce(GCA *gca_src) ;
int  GCAnodeToVoxel(GCA *gca, MRI *mri, int xn, int yn, int zn, int *pxv, 
                    int *pyv, int *pzv) ;
int  GCApriorToVoxel(GCA *gca, MRI *mri, int xn, int yn, int zn, int *pxv, 
                    int *pyv, int *pzv) ;
float GCAcomputeLogImageProbability(GCA *gca, MRI *mri_inputs, MRI *mri_labels,
                                    TRANSFORM *transform) ;
float  GCAcomputeLogSampleProbability(GCA *gca, GCA_SAMPLE *gcas, 
                                      MRI *mri_inputs,
                                      TRANSFORM *transform,int nsamples);
float  GCAcomputeLogSampleProbabilityUsingCoords(GCA *gca, GCA_SAMPLE *gcas, 
                                      MRI *mri_inputs,
                                      TRANSFORM *transform,int nsamples);
float  GCAnormalizedLogSampleProbability(GCA *gca, GCA_SAMPLE *gcas, 
                                      MRI *mri_inputs,
                                      TRANSFORM *transform,int nsamples);
int   GCAremoveOutlyingSamples(GCA *gca, GCA_SAMPLE *gcas, 
                                      MRI *mri_inputs,
                                      TRANSFORM *transform,int nsamples, float nsigma);
int   GCArankSamples(GCA *gca, GCA_SAMPLE *gcas, int nsamples, 
                     int *ordered_indices) ;
MRI  *GCAanneal(MRI *mri_inputs, GCA *gca, MRI *mri_dst,TRANSFORM *transform, 
                int max_iter);
int    GCAsourceVoxelToNode(GCA *gca, MRI *mri, TRANSFORM *transform,
                                 int xv, int yv, int zv, 
                                 int *pxn, int *pyn, int *pzn) ;
int    GCAsourceVoxelToPrior(GCA *gca, MRI *mri, TRANSFORM *transform,
                                 int xv, int yv, int zv, 
                                 int *pxp, int *pyp, int *pzp) ;
int  GCAnodeToSourceVoxel(GCA *gca, MRI *mri, TRANSFORM *transform, 
                          int xv, int yv, int zv,
                          int *pxn, int *pyn, int *pzn) ;
int  GCAnodeToSourceVoxelFloat(GCA *gca, MRI *mri, TRANSFORM *transform, 
                          int xv, int yv, int zv,
                          float *pxn, float *pyn, float *pzn) ;
#if 0
int    GCAsampleStats(GCA *gca, MRI *mri, TRANSFORM *transform, int class, 
                      Real x, Real y, Real z, 
                      Real *pmean, Real *pvar, Real *pprior) ;
#endif



MRI  *GCAannealUnlikelyVoxels(MRI *mri_inputs, GCA *gca, MRI *mri_dst,
                              TRANSFORM *transform, int max_iter, MRI *mri_fixed) ;
GCA_SAMPLE *GCAfindContrastSamples(GCA *gca, int *pnsamples, int min_spacing,
                                 float min_prior) ;
GCA_SAMPLE *GCAfindAllSamples(GCA *gca, int *pnsamples, int *exclude_list,
															int unknown_nbr_spacing) ;
GCA_SAMPLE *GCAfindStableSamples(GCA *gca, int *pnsamples, int min_spacing,
                                 float min_prior, int *exclude_list, int unknown_nbr_spacing) ;
GCA_SAMPLE *GCAfindStableSamplesByLabel(GCA *gca, int nsamples, 
                                        float min_prior) ;
int       GCAtransformSamples(GCA *gca_src, GCA *gca_dst, GCA_SAMPLE *gcas, 
                              int nsamples) ;
int        GCAwriteSamples(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, int nsamples,
                           char *fname) ;
int        GCAtransformAndWriteSamples(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, 
                                       int nsamples,char *fname,TRANSFORM *transform) ;
int        GCAcomputeSampleCoords(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, 
                                  int nsamples,TRANSFORM *transform) ;
MRI        *GCAmri(GCA *gca, MRI *mri) ;
MRI        *GCAlabelMri(GCA *gca, MRI *mri, int label, TRANSFORM *transform) ;
MRI        *GCAbuildMostLikelyVolume(GCA *gca, MRI *mri) ;
MRI        *GCAbuildMostLikelyVolumeForStructure(GCA *gca, MRI *mri_seg, int label, int border, TRANSFORM *transform,
																								 MRI *mri_labels) ;
MRI        *GCAbuildMostLikelyVolumeFrame(GCA *gca, MRI *mri, int frame) ;
MRI *GCAbuildMostLikelyLabelVolume(GCA *gca);
MRI  *GCAlabelProbabilities(MRI *mri_inputs, GCA *gca, MRI *mri_dst, TRANSFORM *transform);
MRI  *GCAcomputeProbabilities(MRI *mri_inputs, GCA *gca, MRI *mri_labels, 
                              MRI *mri_dst, TRANSFORM *transform);

int   GCAcomputeMAPlabelAtLocation(GCA *gca, int xp, int yp, int zp, float *vals,
                             int *pmax_n, float *plog_p) ;
int   GCAcomputeMLElabelAtLocation(GCA *gca, int x, int y, int z, float *vals, int *pmax_n,float *plog_p);
MRI   *GCAconstrainLabelTopology(GCA *gca, MRI *mri_inputs, MRI *mri_src, 
                                 MRI *mri_dst, TRANSFORM *transform) ;
MRI   *GCAexpandLabelIntoWM(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                            MRI *mri_dst, TRANSFORM *transform,MRI *mri_fixed,
                            int target_label) ;
MRI   *GCAexpandVentricle(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                          MRI *mri_dst, TRANSFORM *transform, int target_label) ;
MRI   *GCAexpandCortex(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                       MRI *mri_dst, TRANSFORM *transform) ;
MRI   *GCAnormalizeSamples(MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas, 
                           int nsamples, TRANSFORM *transform, char *ctl_point_fname) ;
MRI *
GCAnormalizeSamplesAllChannels(MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas, int nsamples, 
			       TRANSFORM *transform, char *ctl_point_fname);
float GCAlabelProbability(MRI *mri_src, GCA *gca, TRANSFORM *transform,
                          int x, int y, int z, int label) ;
MRI   *GCAmaxLikelihoodBorders(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                               MRI *mri_dst, TRANSFORM *transform, int max_iter,
                               float min_ratio) ;
int   GCAaccumulateTissueStatistics(GCA *gca, MRI *mri_T1, MRI *mri_PD, 
                                    MRI *mri_parc, TRANSFORM *transform) ;
int   GCAhistogramTissueStatistics(GCA *gca, MRI *mri_T1, MRI *mri_PD, 
                                   MRI *mri_parc, TRANSFORM *transform, char *fname) ;
int   GCAnormalizeTissueStatistics(GCA *gca) ;
int  GCArenormalize(MRI *mri_in, MRI *mri_labeled, GCA *gca, TRANSFORM *transform) ;
int  GCAmapRenormalize(GCA *gca, MRI *mri, TRANSFORM *transform) ;
int  GCAmapRenormalizeWithAlignment(GCA *gca, MRI *mri, TRANSFORM *transform, FILE *logfp, char *base_name, LTA **plta) ;
int  GCArenormalizeAdaptive(MRI *mri_in, MRI *mri_labeled, GCA *gca, TRANSFORM *transform,
                            int wsize, float pthresh) ;
int  GCArenormalizeLabels(MRI *mri_in, MRI *mri_labeled, GCA *gca, TRANSFORM *transform) ;
MRI   *GCArelabel_cortical_gray_and_white(GCA *gca, MRI *mri_inputs, 
                                          MRI *mri_src, MRI *mri_dst,TRANSFORM *transform);

int GCAdump(GCA *gca, MRI *mri, int x, int y, int z, TRANSFORM *transform, FILE *fp, 
            int verbose) ;
int GCArenormalizeIntensities(GCA *gca, int *labels, float *intensities, 
                              int num) ;
int GCArenormalizeIntensitiesAbsolute(GCA *gca, int *labels, 
                                      float *intensities, int num) ;
int GCArenormalizeToExample(GCA *gca, MRI *mri_seg, MRI *mri_T1) ;

int     GCAlabelMode(GCA *gca, int label, float *modes) ;
int     GCAlabelMean(GCA *gca, int label, float *means) ;
int     GCAlabelMeanFromImage(GCA *gca, TRANSFORM *transform, MRI *mri, int label, float *means) ;
MATRIX  *GCAlabelCovariance(GCA *gca, int label, MATRIX *m_total) ;
int     GCAregularizeConditionalDensities(GCA *gca, float smooth) ;
int     GCAmeanFilterConditionalDensities(GCA *gca, float navgs) ;
int     GCArenormalizeToFlash(GCA *gca, char *tissue_parms_fname, MRI *mri) ;
int     GCAhistoScaleImageIntensities(GCA *gca, MRI *mri) ;
int     GCAhisto(GCA *gca, int nbins, int **pcounts) ;
int     GCAcomputeVoxelLikelihoods(GCA *gca, MRI *mri_in, int x, int y, int z, 
                                   TRANSFORM *transform, int *labels, double *likelihoods);
GCA_PRIOR *getGCAP(GCA *gca, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv) ;
float getPrior(GCA_PRIOR *gcap, int label) ;
int   GCApriorToNode(GCA *gca, int xp, int yp, int zp, int *pxn, int *pyn, int *pzn) ;
int   GCAfreeGibbs(GCA *gca) ;
GC1D *GCAfindPriorGC(GCA *gca, int xp, int yp, int zp,int label) ;
int  GCApriorToSourceVoxel(GCA *gca, MRI *mri, TRANSFORM *transform, int xp, int yp, int zp, 
                           int *pxv, int *pyv, int *pzv) ;
int  GCApriorToSourceVoxelFloat(GCA *gca, MRI *mri, TRANSFORM *transform, 
                                int xp, int yp, int zp, 
                                float *pxv, float *pyv, float *pzv) ;
int GCArenormalizeFromAtlas(GCA *gca, GCA *gca_template) ;
GC1D *GCAfindGC(GCA *gca, int x, int y, int z,int label) ;
GC1D *GCAfindSourceGC(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z, int label) ;
int GCAlabelExists(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z, int label) ;

VECTOR *load_mean_vector(GC1D *gc, VECTOR *v_means, int ninputs) ;
MATRIX *load_covariance_matrix(GC1D *gc, MATRIX *m_cov, int ninputs) ;
MATRIX *load_inverse_covariance_matrix(GC1D *gc, MATRIX *m_cov, int ninputs) ;
double covariance_determinant(GC1D *gc, int ninputs) ;
void load_vals(MRI *mri_inputs, float x, float y, float z, float *vals, int ninputs) ;
int    GCAisPossible(GCA *gca, MRI *mri, int label, TRANSFORM *transform, int x, int y, int z, int use_mrf) ;
double GCAcomputePosteriorDensity(GCA_PRIOR *gcap, GCA_NODE *gcan, int n, float *vals, int ninputs) ;
double GCAcomputeConditionalDensity(GC1D *gc, float *vals, int ninputs, int label) ;
double GCAmahDistIdentityCovariance(GC1D *gc, float *vals, int ninputs) ;
double GCAmahDist(GC1D *gc, float *vals, int ninputs) ;
int    GCAfreeSamples(GCA_SAMPLE **pgcas, int nsamples) ;
double GCAsampleMahDist(GCA_SAMPLE *gcas, float *vals, int ninputs) ;
int GCAnormalizePD(GCA *gca, MRI *mri_inputs, TRANSFORM *transform) ;
GCA *GCAcreateFlashGCAfromParameterGCA(GCA *gca_parms, double *TR, double *fa, double *TE, int nflash, double lambda);
GCA *GCAcreateWeightedFlashGCAfromParameterGCA(GCA *gca_parms, double *TR, double *fa, double *TE, int nflash, double *wts, double lambda);
GCA *GCAcreateFlashGCAfromFlashGCA(GCA *gca_parms, double *TR, double *fa, double *TE, int nflash);
int GCAfixSingularCovarianceMatrices(GCA *gca) ;
int GCAregularizeCovariance(GCA *gca, float regularize) ;
int GCAnormalizeMeans(GCA *gca, float target) ;
double GCAcomputeConditionalLogDensity(GC1D *gc, float *vals, int ninputs, int label) ;
double GCAcomputeNormalizedConditionalDensity(GCA *gca, int xp, int yp, int zp, float *vals, int label);
MRI    *GCArelabelNonbrain(GCA *gca, MRI *mri_inputs, MRI *mri_src, MRI *mri_dst, TRANSFORM *transform) ;
int    GCAreplaceLabels(GCA *gca, int in_label, int out_label) ;
///////////////////////////////////////////////////////////////////

// setting up global node and prior volume parameters
void GCAsetup(GCA *gca);
void GCAreinit(MRI *mri, GCA *gca); // reinit gca with mri values
void GCAcleanup(GCA *gca);
void GCAcopyDCToMRI(GCA *gca, MRI *mri); // copy direction cosine info to MRI
void GCAsetVolGeom(GCA *gca, VOL_GEOM *vg); 
int GCAregularizeCovarianceMatrices(GCA *gca, double lambda) ;
int GCAreplaceRightWithLeft(GCA *gca) ;
int GCAcomputeLabelStats(GCA *gca, int target_label, float *pvar, float *means);
GCA_NODE *GCAbuildRegionalGCAN(GCA *gca, int x, int y, int z, int wsize) ;
int set_mean_vector(GC1D *gc, VECTOR *v_means, int ninputs) ;
int set_covariance_matrix(GC1D *gc, MATRIX *m_cov, int ninputs) ;
int GCAmaxLikelihoodLabel(GCA_NODE *gcan, float *vals, int ninputs, float *plikelihood) ;
int GCAfreeRegionalGCAN(GCA_NODE **pgcan) ;
GCA *GCAcompactify(GCA *gca);
MRI *GCAreplaceImpossibleLabels(MRI *mri_inputs, GCA *gca, MRI *mri_in_labels, MRI *mri_out_labels, TRANSFORM *transform) ;
GC1D *alloc_gcs(int nlabels, int flags, int ninputs) ;
int free_gcs(GC1D *gcs, int nlabels, int ninputs) ;
int GCAmapRenormalizeByClass(GCA *gca, MRI *mri, TRANSFORM *transform) ;
extern int Ggca_x, Ggca_y, Ggca_z, Ggca_label, Ggca_nbr_label ;
extern char *G_write_probs ;
MRI *GCAmarkImpossible(GCA *gca, MRI *mri_labeled, MRI *mri_dst, TRANSFORM *transform) ;
int GCAclassMode(GCA *gca, int class, float *modes) ;
int GCAcomputeLabelMeansAndCovariances(GCA *gca, int target_label, MATRIX **p_mcov, VECTOR **p_vmeans) ;
double GCAgibbsImpossibleConfiguration(GCA *gca, 
																			 MRI *mri_labels, 
																			 int x, int y, int z, 
																			 TRANSFORM *transform) ;
MRI *GCAlabelWMandWMSAs(GCA *gca, MRI *mri_inputs, MRI *mri_src_labels, MRI *mri_dst_labels, TRANSFORM *transform);

#define GCA_DEFAULT_NOISE_PARAMETER  1

#endif
