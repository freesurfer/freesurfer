#define USE_ROUND 1
#define BEVIN_GCACOMPUTELOGSAMPLEPROBABILITY_REPRODUCIBLE

/**
 * @brief Core routines implementing the Gaussian Classifier Atlas mechanism
 *
 * Routines implementing the mechanism of encoding voxel label information
 * based on probabilistic information estimated from a training set.
 *
 * Reference:
 * "Whole Brain Segmentation: Automated Labeling of Neuroanatomical
 * Structures in the Human Brain", Fischl et al.
 * (2002) Neuron, 33:341-355.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2011-2015 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "faster_variants.h"
#include "romp_support.h"

#include "cma.h"
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "flash.h"
#include "gca.h"
#include "intensity_eig.h"
#include "macros.h"
#include "mri.h"
#include "mri2.h"
#include "mrimorph.h"
#include "mrisegment.h"
#include "numerics.h"
#include "proto.h"
#include "tags.h"
#include "talairachex.h"
#include "transform.h"
#include "utils.h"

#include "affine.h"

#include "timer.h"

#include "znzlib.h"

#if WITH_DMALLOC
#include <dmalloc.h>
#endif

#ifdef HAVE_OPENMP
//#undef HAVE_OPENMP
#endif

extern const char *Progname;

int Ggca_label = -1;
int Ggca_nbr_label = -1;
int Ggca_x = -1;
int Ggca_y = -1;
int Ggca_z = -1;
int Gxp = -1;
int Gyp = -1;
int Gzp = -1;
int Gxn = -1;  // 32;
int Gyn = -1;  // 21;
int Gzn = -1;  // 32;
char *G_write_probs = NULL;

/* this is the hack section */
static double PRIOR_FACTOR = 0.1;
#define LABEL_UNDETERMINED 255

static int total_pruned = 0;

#define MIN_DET 1e-7
/* less than this, and the covariance matrix is poorly conditioned */

#define MAX_LABELS_PER_GCAN 50
#define MAX_DIFFERENT_LABELS 5000
#define MIN_VAR (5 * 5) /* should make this configurable */
#define BIG_AND_NEGATIVE -10000000.0
#define VERY_UNLIKELY 1e-10
#define UNKNOWN_DIST 4 /* within 4 mm of some known label */
#define GCA_OLD_VERSION 2.0
#define GCA_UCHAR_VERSION 4.0  // labels were uchars in file
#define GCA_INT_VERSION 5.0    // labels are ints in file
#define DEFAULT_MAX_LABELS_PER_GCAN 4

static float get_node_prior(GCA *gca, int label, int xn, int yn, int zn);
// static int gcapBrainIsPossible(GCA_PRIOR *gcap) ;

static double gcaGibbsLogPosterior(GCA *gca,
                                   MRI *mri_labels,
                                   float *vals,
                                   int label,
                                   int x,
                                   int y,
                                   int z,
                                   GCA_PRIOR *gcap,
                                   GCA_NODE *gcan,
                                   TRANSFORM *transform);
#define INTERP_PRIOR 0
#if INTERP_PRIOR
static float gcaComputePrior(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z, int label);
#endif

static int gcaScale(GCA *gca, int *labels, int *contra_labels, float *scales, int nlabels, int dir);
static int gcaMaxPriorLabel(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z);
static GCA_NODE *gcaBuildNbhdGCAN(GCA *gca, int x, int y, int z, float sigma);
static GCA_PRIOR *gcaBuildNbhdGCAP(GCA *gca, int x, int y, int z, float sigma, GCA *gca_smooth);

static HISTOGRAM *gcaGetLabelHistogram(GCA *gca, int label, int frame, int border);
int GCAmaxLabel(GCA *gca);
static int gcaRelabelSegment(GCA *gca, TRANSFORM *transform, MRI *mri_inputs, MRI *mri_dst, MRI_SEGMENT *mseg);

static int gcapGetMaxPriorLabel(GCA_PRIOR *gcap, double *p_prior);
double compute_partial_volume_log_posterior(GCA *gca, GCA_NODE *gcan, GCA_PRIOR *gcap, float *vals, int l1, int l2);
static double gcaComputeSampleConditionalLogDensity(GCA_SAMPLE *gcas, float *vals, int ninputs, int label);
static VECTOR *load_sample_mean_vector(GCA_SAMPLE *gcas, VECTOR *v_means, int ninputs);
static MATRIX *load_sample_covariance_matrix(GCA_SAMPLE *gcas, MATRIX *m_cov, int ninputs);
static double sample_covariance_determinant(GCA_SAMPLE *gcas, int ninputs);
static GC1D *gcanGetGC(GCA_NODE *gcan, int label);
static GC1D *findGCInWindow(GCA *gca, int x, int y, int z, int label, int wsize);


static int gcaCheck(GCA *gca);
double gcaVoxelLogPosterior(GCA *gca, MRI *mri_labels, MRI *mri_inputs, int x, int y, int z, TRANSFORM *transform);
static double gcaGibbsImpossibleConfiguration(GCA *gca, MRI *mri_labels, int x, int y, int z, TRANSFORM *transform);
static GCA_SAMPLE *gcaExtractLabelAsSamples(
    GCA *gca, MRI *mri_labeled, TRANSFORM *transform, int *pnsamples, int label);
static GCA_SAMPLE *gcaExtractThresholdedRegionLabelAsSamples(GCA *gca,
                                                             MRI *mri_labeled,
                                                             TRANSFORM *transform,
                                                             int *pnsamples,
                                                             int label,
                                                             int xn,
                                                             int yn,
                                                             int zn,
                                                             int wsize,
                                                             float pthresh);
static double gcaComputeSampleConditionalDensity(GCA_SAMPLE *gcas, float *vals, int ninputs, int label);
static double gcaComputeLogDensity(GC1D *gc, float *vals, int ninputs, float prior, int label);
static double gcaComputeSampleLogDensity(GCA_SAMPLE *gcas, float *vals, int ninputs);
static int GCAupdateNode(GCA *gca, MRI *mri, int xn, int yn, int zn, float *vals, int label, GCA *gca_prune, int noint);
static int GCAupdateNodeCovariance(GCA *gca, MRI *mri, int xn, int yn, int zn, float *vals, int label);
static int GCAupdatePrior(GCA *gca, MRI *mri, int xn, int yn, int zn, int label);
static int GCAupdateNodeGibbsPriors(GCA *gca, MRI *mri, int xn, int yn, int zn, int x, int y, int z, int label);
static int different_nbr_max_labels(GCA *gca, int x, int y, int z, int wsize, int label);
static int gcaRegionStats(GCA *gca,
                          int x0,
                          int y0,
                          int z0,
                          int wsize,
                          float *priors,
                          float *vars,
                          float means[MAX_DIFFERENT_LABELS][MAX_GCA_INPUTS]);
static int gcaFindBestSample(GCA *gca, int x, int y, int z, int best_label, int wsize, GCA_SAMPLE *gcas);

static int mriFillRegion(MRI *mri, int x, int y, int z, int frame, float fill_val, int whalf);
static int gcaFindMaxPriors(GCA *gca, float *max_priors);
static int gcaFindIntensityBounds(GCA *gca, float *pmin, float *pmax);
static int dump_gcan(GCA *gca, GCA_NODE *gcan, FILE *fp, int verbose, GCA_PRIOR *gcap);
static GCA_NODE *findSourceGCAN(GCA *gca, MRI *mri_src, TRANSFORM *transform, int x, int y, int z);
static int borderVoxel(MRI *mri, int x, int y, int z);
static int GCAmaxLikelihoodBorderLabel(
    GCA *gca, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform, int x, int y, int z, float min_ratio);

float getPrior(GCA_PRIOR *gcap, int label);
GCA_PRIOR *getGCAP(GCA *gca, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv);
GCA_PRIOR *getGCAPfloat(GCA *gca, MRI *mri, TRANSFORM *transform, float xv, float yv, float zv);
static int gcaNodeToPrior(const GCA *gca, int xn, int yn, int zn, int *pxp, int *pyp, int *pzp);
static HISTOGRAM *gcaHistogramSamples(
    GCA *gca, GCA_SAMPLE *gcas, MRI *mri, TRANSFORM *transform, int nsamples, HISTOGRAM *histo, int frame);
int GCApriorToNode(const GCA *gca, int xp, int yp, int zp, int *pxn, int *pyn, int *pzn);

/* arrays for indexing 6-connected neighbors */
static int xnbr_offset[] = {1, -1, 0, 0, 0, 0};
static int ynbr_offset[] = {0, 0, 1, -1, 0, 0};
static int znbr_offset[] = {0, 0, 0, 0, 1, -1};
static int boundsCheck(int *pix, int *piy, int *piz, MRI *mri);

static void set_equilavent_classes(int *equivalent_classes);

#ifdef FASTER_MRI_EM_REGISTER
static void load_vals_xyzInt(const MRI *mri_inputs, int x, int y, int z, float *vals, int ninputs);
static double gcaComputeSampleLogDensity_1_input(GCA_SAMPLE *gcas, float val);
static double gcaComputeSampleConditionalLogDensity_1_input(GCA_SAMPLE *gcas, float val, int label);
static double GCAsampleMahDist_1_input(GCA_SAMPLE *gcas, float val);
#endif


double GCAimageLikelihoodAtNode(GCA *gca, GCA_NODE *gcan, float *vals, int label)
{
  GC1D *gc;
  float likelihood;

  gc = gcanGetGC(gcan, label);
  if (gc == NULL) return (0.0f);
  likelihood = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, label);
  return (likelihood);
}

static int initialize_ventricle_alignment(MRI *mri_seg, MRI *mri, MATRIX *m_L, char *base_name, int border, int label);

void GCAsetVolGeom(GCA *gca, VOL_GEOM *vg)
{
  vg->width = gca->width;
  vg->height = gca->height;
  vg->depth = gca->depth;
  vg->xsize = gca->xsize;
  vg->ysize = gca->ysize;
  vg->zsize = gca->zsize;
  vg->x_r = gca->x_r;
  vg->y_r = gca->y_r;
  vg->z_r = gca->z_r;
  vg->x_a = gca->x_a;
  vg->y_a = gca->y_a;
  vg->z_a = gca->z_a;
  vg->x_s = gca->x_s;
  vg->y_s = gca->y_s;
  vg->z_s = gca->z_s;
  vg->c_r = gca->c_r;
  vg->c_a = gca->c_a;
  vg->c_s = gca->c_s;
  vg->valid = 1;
}

void GCAcopyDCToMRI(GCA *gca, MRI *mri)
{
  mri->x_r = gca->x_r;
  mri->y_r = gca->y_r;
  mri->z_r = gca->z_r;
  mri->x_a = gca->x_a;
  mri->y_a = gca->y_a;
  mri->z_a = gca->z_a;
  mri->x_s = gca->x_s;
  mri->y_s = gca->y_s;
  mri->z_s = gca->z_s;
  mri->c_r = gca->c_r;
  mri->c_a = gca->c_a;
  mri->c_s = gca->c_s;
  mri->ras_good_flag = 1;
  MRIreInitCache(mri);

  // mri->i_to_r__ = extract_i_to_r(mri);
  MATRIX *tmp;
  tmp = extract_i_to_r(mri);
  AffineMatrixAlloc(&(mri->i_to_r__));
  SetAffineMatrix(mri->i_to_r__, tmp);

  mri->r_to_i__ = extract_r_to_i(mri);
  MatrixFree(&tmp);
}

void GCAcopyDCToGCA(GCA *gca, GCA *gca_dst)
{
  gca_dst->x_r = gca->x_r;
  gca_dst->y_r = gca->y_r;
  gca_dst->z_r = gca->z_r;
  gca_dst->x_a = gca->x_a;
  gca_dst->y_a = gca->y_a;
  gca_dst->z_a = gca->z_a;
  gca_dst->x_s = gca->x_s;
  gca_dst->y_s = gca->y_s;
  gca_dst->z_s = gca->z_s;
  gca_dst->c_r = gca->c_r;
  gca_dst->c_a = gca->c_a;
  gca_dst->c_s = gca->c_s;
  //
  gca_dst->node_i_to_r__ = MatrixCopy(gca->node_i_to_r__, gca_dst->node_i_to_r__);
  gca_dst->node_r_to_i__ = MatrixCopy(gca->node_r_to_i__, gca_dst->node_r_to_i__);
  gca_dst->prior_i_to_r__ = MatrixCopy(gca->prior_i_to_r__, gca_dst->prior_i_to_r__);
  gca_dst->prior_r_to_i__ = AffineMatrixCopy(gca->prior_r_to_i__, gca_dst->prior_r_to_i__);
  gca_dst->tal_i_to_r__ = MatrixCopy(gca->tal_i_to_r__, gca_dst->tal_i_to_r__);
  gca_dst->tal_r_to_i__ = MatrixCopy(gca->tal_r_to_i__, gca_dst->tal_r_to_i__);
  //
  gca_dst->mri_prior__ = MRIcopy(gca->mri_prior__, gca_dst->mri_prior__);
  gca_dst->mri_node__ = MRIcopy(gca->mri_node__, gca_dst->mri_node__);
  gca_dst->mri_tal__ = MRIcopy(gca->mri_tal__, gca_dst->mri_tal__);
}

void GCAcleanup(GCA *gca)
{
  if (gca->mri_node__) {
    MRIfree(&gca->mri_node__);
    gca->mri_node__ = 0;
  }
  if (gca->mri_prior__) {
    MRIfree(&gca->mri_prior__);
    gca->mri_prior__ = 0;
  }
  if (gca->mri_tal__) {
    MRIfree(&gca->mri_tal__);
    gca->mri_tal__ = 0;
  }
  if (gca->node_i_to_r__) {
    MatrixFree(&(gca->node_i_to_r__));
    gca->node_i_to_r__ = 0;
  }
  if (gca->node_r_to_i__) {
    MatrixFree(&(gca->node_r_to_i__));
    gca->node_r_to_i__ = 0;
  }
  if (gca->prior_i_to_r__) {
    MatrixFree(&(gca->prior_i_to_r__));
    gca->prior_i_to_r__ = 0;
  }
  if (gca->prior_r_to_i__) {
    AffineMatrixFree(&(gca->prior_r_to_i__));
    gca->prior_r_to_i__ = 0;
  }
  if (gca->tal_i_to_r__) {
    MatrixFree(&(gca->tal_i_to_r__));
    gca->tal_i_to_r__ = 0;
  }
  if (gca->tal_r_to_i__) {
    MatrixFree(&(gca->tal_r_to_i__));
    gca->tal_r_to_i__ = 0;
  }
  if (gca->tmp__) {
    MatrixFree(&gca->tmp__);
    gca->tmp__ = 0;
  }
}

// set up mri's used in GCA
void GCAsetup(GCA *gca)
{
  // set up node part /////////////////////
  if (gca->mri_node__) {
    MRIfree(&gca->mri_node__);
    gca->mri_node__ = 0;
  }
  gca->mri_node__ = MRIallocHeader(gca->node_width, gca->node_height, gca->node_depth, MRI_UCHAR, 1);
  /* Copy the voxel resolutions.  Set the defaults */
  gca->mri_node__->xsize = gca->xsize * gca->node_spacing;
  gca->mri_node__->ysize = gca->ysize * gca->node_spacing;
  gca->mri_node__->zsize = gca->zsize * gca->node_spacing;

  // this will recalculate i_to_r__ etc and
  // thus xsize etc must be set correctly
  GCAcopyDCToMRI(gca, gca->mri_node__);

  // fprintf(stdout, "node voxelToRAS\n");
  // MATRIX *mnode = extract_i_to_r(gca->mri_node__);
  // MatrixPrint(stdout, mnode);
  // MatrixFree(&mnode);

  // setup prior part //////////////////////////////////////
  if (gca->mri_prior__) {
    MRIfree(&gca->mri_prior__);
    gca->mri_prior__ = 0;
  }
  gca->mri_prior__ = MRIallocHeader(gca->prior_width, gca->prior_height, gca->prior_depth, MRI_UCHAR, 1);

  /* Copy the voxel resolutions.  Set the defaults */
  gca->mri_prior__->xsize = gca->xsize * gca->prior_spacing;
  gca->mri_prior__->ysize = gca->ysize * gca->prior_spacing;
  gca->mri_prior__->zsize = gca->zsize * gca->prior_spacing;

  GCAcopyDCToMRI(gca, gca->mri_prior__);

  // fprintf(stdout, "prior voxelToRAS\n");
  // MATRIX *mprior = extract_i_to_r(gca->mri_prior__);
  // MatrixPrint(stdout, mprior);
  // MatrixFree(&mprior);

  // set up the default talairach volume ////////////////
  if (gca->mri_tal__) {
    MRIfree(&gca->mri_tal__);
    gca->mri_tal__ = 0;
  }
  gca->mri_tal__ = MRIallocHeader(gca->width, gca->height, gca->depth, MRI_UCHAR, 1);

  /* Copy the voxel resolutions.  Set the defaults */
  gca->mri_tal__->xsize = gca->xsize;
  gca->mri_tal__->ysize = gca->ysize;
  gca->mri_tal__->zsize = gca->zsize;

  GCAcopyDCToMRI(gca, gca->mri_tal__);

  // fprintf(stdout, "tal voxelToRAS\n");
  // mtal = extract_i_to_r(gca->mri_tal__);
  // MatrixPrint(stdout, mtal);
  // MatrixFree(&mtal);
  if (gca->node_i_to_r__) {
    MatrixFree(&(gca->node_i_to_r__));
    gca->node_i_to_r__ = 0;
  }
  if (gca->node_r_to_i__) {
    MatrixFree(&(gca->node_r_to_i__));
    gca->node_r_to_i__ = 0;
  }
  if (gca->prior_i_to_r__) {
    MatrixFree(&(gca->prior_i_to_r__));
    gca->prior_i_to_r__ = 0;
  }
  if (gca->prior_r_to_i__) {
    AffineMatrixFree(&(gca->prior_r_to_i__));
    gca->prior_r_to_i__ = 0;
  }
  if (gca->tal_i_to_r__) {
    MatrixFree(&(gca->tal_i_to_r__));
    gca->tal_i_to_r__ = 0;
  }
  if (gca->tal_r_to_i__) {
    MatrixFree(&(gca->tal_r_to_i__));
    gca->tal_r_to_i__ = 0;
  }
  if (gca->tmp__) {
    MatrixFree(&(gca->tmp__));
    gca->tmp__ = 0;
  }
  gca->node_i_to_r__ = extract_i_to_r(gca->mri_node__);
  gca->node_r_to_i__ = extract_r_to_i(gca->mri_node__);
  gca->prior_i_to_r__ = extract_i_to_r(gca->mri_prior__);

  MATRIX *tmp = extract_r_to_i(gca->mri_prior__);
  AffineMatrixAlloc(&(gca->prior_r_to_i__));
  SetAffineMatrix(gca->prior_r_to_i__, tmp);
  MatrixFree(&tmp);

  gca->tal_i_to_r__ = extract_i_to_r(gca->mri_tal__);
  gca->tal_r_to_i__ = extract_r_to_i(gca->mri_tal__);
  gca->tmp__ = MatrixAlloc(4, 4, MATRIX_REAL);
}

// using the values of mri, modify gca
// now we can keep track of voxelToRAS info even for non-standard type
// This has the potential of changing direction cosines
void GCAreinit(MRI *mri, GCA *gca)
{
  // Keep the direction cosine the same to avoid used by
  // different mri
  // modify direction cosines etc.
  gca->x_r = mri->x_r;
  gca->y_r = mri->y_r;
  gca->z_r = mri->z_r;
  gca->x_a = mri->x_a;
  gca->y_a = mri->y_a;
  gca->z_a = mri->z_a;
  gca->x_s = mri->x_s;
  gca->y_s = mri->y_s;
  gca->z_s = mri->z_s;
  gca->c_r = mri->c_r;
  gca->c_a = mri->c_a;
  gca->c_s = mri->c_s;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("gca reinit c_(ras) = (%.2f, %.2f, %.2f)\n", gca->c_r, gca->c_a, gca->c_s);

  // modify width height depth
  if (gca->width != mri->width) fprintf(stdout, "gca width modified from %d to %d\n", gca->width, mri->width);
  if (gca->height != mri->height) fprintf(stdout, "gca height modified from %d to %d\n", gca->height, mri->height);
  if (gca->depth != mri->depth) fprintf(stdout, "gca depth modified from %d to %d\n", gca->depth, mri->depth);
  gca->width = mri->width;
  gca->height = mri->height;
  gca->depth = mri->depth;

  if (gca->xsize != mri->xsize) fprintf(stdout, "gca xsize modified from %.3f to %.3f\n", gca->xsize, mri->xsize);
  if (gca->ysize != mri->ysize) fprintf(stdout, "gca ysize modified from %.3f to %.3f\n", gca->ysize, mri->ysize);
  if (gca->zsize != mri->zsize) fprintf(stdout, "gca zsize modified from %.3f to %.3f\n", gca->zsize, mri->zsize);
  gca->xsize = mri->xsize;
  gca->ysize = mri->ysize;
  gca->zsize = mri->zsize;
  //
  GCAsetup(gca);
  fflush(stdout);
}

GCA_PRIOR *getGCAP(GCA *gca, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv)
{
  int xp, yp, zp;
  GCA_PRIOR *gcap = NULL;

  if (!GCAsourceVoxelToPrior(gca, mri, transform, xv, yv, zv, &xp, &yp, &zp)) {
    gcap = &gca->priors[xp][yp][zp];
  }

  return (gcap);
}

GCA_PRIOR *getGCAPfloat(GCA *gca, MRI *mri, TRANSFORM *transform, float xv, float yv, float zv)
{
  int xp, yp, zp;
  GCA_PRIOR *gcap = NULL;

  if (!GCAsourceFloatVoxelToPrior(gca, mri, transform, xv, yv, zv, &xp, &yp, &zp)) {
    gcap = &gca->priors[xp][yp][zp];
  }

  return (gcap);
}

GCA_NODE *getGCAN(GCA *gca, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv)
{
  int xn, yn, zn;
  GCA_NODE *gcan = NULL;

  if (!GCAsourceVoxelToNode(gca, mri, transform, xv, yv, zv, &xn, &yn, &zn)) {
    gcan = &gca->nodes[xn][yn][zn];
  }

  return (gcan);
}

float getPrior(GCA_PRIOR *gcap, int label)
{
  int n;

  if (gcap == NULL) {
    return (VERY_UNLIKELY);
  }

  // find the label
  for (n = 0; n < gcap->nlabels; n++)
    if (gcap->labels[n] == label) {
      break;
    }
  // cannot find it
  if (n >= gcap->nlabels) {
    if (gcap->total_training > 0) {
      return (0.1f / (float)gcap->total_training); /* make it unlikely */
    }
    else {
      return (VERY_UNLIKELY);
    }
  }
  // return found one
  if (DZERO(gcap->priors[n])) {
    return (VERY_UNLIKELY);
  }
  return (gcap->priors[n]);
}

int GCAisPossible(GCA *gca, MRI *mri, int label, TRANSFORM *transform, int x, int y, int z, int use_mrf)
{
  int n, i, found, nbr_label, xnbr, ynbr, znbr, xn, yn, zn;
  GCA_PRIOR *gcap = NULL;
  GC1D *gc;

  gcap = getGCAP(gca, mri, transform, x, y, z);
  if (gcap == NULL) {
    return (0);
  }
  // if label found
  for (found = n = 0; n < gcap->nlabels; n++)
    if (gcap->labels[n] == label) {
      found = 1;
      break;
    }
  if (found == 0) {
    return (0);
  }
  if (use_mrf == 0) {
    return (1);
  }

  if (GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn) != NO_ERROR) {
    return (1);
  }
  gc = GCAfindGC(gca, xn, yn, zn, label);
  if (gc == NULL) {
    return (1);
  }
  for (i = 0; i < GIBBS_NEIGHBORS; i++) {
    found = 0;
    xnbr = mri->xi[x + xnbr_offset[i]];  // xnbr_offset 1, -1, 0,  0, 0,  0
    ynbr = mri->yi[y + ynbr_offset[i]];  // ynbr_offset 0,  0, 1, -1, 0,  0
    znbr = mri->zi[z + znbr_offset[i]];  // znbr_offset 0,  0, 0,  0, 1, -1
    nbr_label = nint(MRIgetVoxVal(mri, xnbr, ynbr, znbr, 0));
    for (n = 0; n < gc->nlabels[i]; n++)
      if (gc->labels[i][n] == nbr_label) {
        found = 1;
        break;
      }
    if (found == 0)  // this pair never occurred here
    {
      return (0);
    }
  }
  return (1);
}


static float get_node_prior(GCA *gca, int label, int xn, int yn, int zn)
{
  int xp, yp, zp;
  GCA_PRIOR *gcap;

  if (gcaNodeToPrior(gca, xn, yn, zn, &xp, &yp, &zp) == NO_ERROR) {
    gcap = &gca->priors[xp][yp][zp];
    if (gcap == NULL) {
      return (VERY_UNLIKELY);
    }
    return (getPrior(gcap, label));
  }
  else {
    return (VERY_UNLIKELY);
  }
}

// use the same function for bounds checking
// if (ix, iy, iz) is within mri volume, returns NO_ERROR
//                                       otherwise ERROR_BADPARAM
static int boundsCheck(int *pix, int *piy, int *piz, MRI *mri)
{
  int ix = *pix;
  int iy = *piy;
  int iz = *piz;
  int errCode = NO_ERROR;  // initialize

  if (ix < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (iy < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (iz < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (ix > mri->width - 1) {
    errCode = ERROR_BADPARM;
  }
  else if (iy > mri->height - 1) {
    errCode = ERROR_BADPARM;
  }
  else if (iz > mri->depth - 1) {
    errCode = ERROR_BADPARM;
  }


  return errCode;
}

static int gcaNodeToPrior(const GCA *gca, int xn, int yn, int zn, int *pxp, int *pyp, int *pzp)
{
  // initialize errCode
  int errCode = NO_ERROR;

// see GCApriorToNode comment
// node_spacing > prior_spacing
// integer operation is floor
  *pxp = (int)((double)xn * gca->node_spacing / gca->prior_spacing);
  *pyp = (int)((double)yn * gca->node_spacing / gca->prior_spacing);
  *pzp = (int)((double)zn * gca->node_spacing / gca->prior_spacing);
  // no error should occur
  return errCode;
}

int GCApriorToNode(const GCA *gca, int xp, int yp, int zp, int *pxn, int *pyn, int *pzn)
{
  int errCode = NO_ERROR;
/////////////////////////////////////////////////////////////////////
// Note that prior and node share the common direction cosines/translation
//
//      prior has voxelToRAS = [ R | T ][prior_size | 0] = X [prior_size | 0]
//                             [ 0 | 1 ][    0      | 1]     [    0      | 1]
//      node  has voxelToRAS = [ R | T ][ node_size | 0] = X [ node_size | 0]
//                             [ 0 | 1 ][    0      | 1]     [    0      | 1]
//
//        prior   --->  RAS
//          |            |
//          |            | identity
//          V            V
//         node    ---> RAS
//                                  -1   -1
//    priorToNode = [ node_size | 0]  * X   * X * [ prior_size | 0]
//                  [     0     | 1]              [     0      | 1]
//
//                = [         -1             |     ]
//                  [ node_size * prior_size |  0  ]
//                  [           0            |  1  ]
//
// Note that node_spacing > prior_spacing and thus the following will do
// .5 becomes 0.
// but this is OK, since the correspondence is
//               |  0  |   1  |   prior
//               |     0      |   node
// i.e. 1 becomes 0

// make sure this agrees with sourceVoxelTo[Node,Prior]! (BRF)
  *pxn = (int)((double)xp * gca->prior_spacing / gca->node_spacing);
  *pyn = (int)((double)yp * gca->prior_spacing / gca->node_spacing);
  *pzn = (int)((double)zp * gca->prior_spacing / gca->node_spacing);
  // no error should occur
  return errCode;
}

static int dump_gcan(GCA *gca, GCA_NODE *gcan, FILE *fp, int verbose, GCA_PRIOR *gcap)
{
  int n, i, j, n1;
  GC1D *gc;
  float prior;
  VECTOR *v_means = NULL;
  MATRIX *m_cov = NULL;

  if (gcap == NULL) {
    return -1;
  }
  // print out gcan labels
  fprintf(fp,
          "node labels at this point -----"
          "---------------------------------\n");
  for (n = 0; n < gcan->nlabels; n++) {
    // find the same label in gcap
    for (n1 = 0; n1 < gcap->nlabels; n1++)
      if (gcap->labels[n1] == gcan->labels[n]) {
        break;
      }
    // the following won't happen
    if (n1 >= gcap->nlabels) {
      continue;
    }
    prior = getPrior(gcap, gcan->labels[n]);
    if (FZERO(prior)) {
      continue;
    }
    // get gc for this label
    gc = &gcan->gcs[n];
    v_means = load_mean_vector(gc, v_means, gca->ninputs);
    m_cov = load_covariance_matrix(gc, m_cov, gca->ninputs);
    fprintf(fp,
            "%d: label %s (%d), prior %2.3f, ntr %d, means ",
            n,
            cma_label_to_name(gcan->labels[n]),
            gcan->labels[n],
            prior,
            gc->ntraining);
    for (i = 0; i < gca->ninputs; i++) {
      fprintf(fp, "%2.1f ", VECTOR_ELT(v_means, i + 1));
    }
    fprintf(fp, ", cov:\n");
    MatrixPrint(fp, m_cov);

    if (verbose) {
      for (i = 0; i < gca->ninputs; i++) {
        for (j = 0; j < gca->ninputs; j++) {
          fprintf(fp, "%2.1f\t", *MATRIX_RELT(m_cov, i + 1, j + 1));
        }
        fprintf(fp, "\n");
      }
      // 6 neighbors
      for (i = 0; i < GIBBS_NEIGHBORS; i++) {
        fprintf(
            fp, "\tnbr %d (%d,%d,%d): %d labels\n", i, xnbr_offset[i], ynbr_offset[i], znbr_offset[i], gc->nlabels[i]);
        fflush(fp);
        for (j = 0; j < gc->nlabels[i]; j++) {
          fprintf(fp, "\t\tlabel %s, prior %2.3f\n", cma_label_to_name(gc->labels[i][j]), gc->label_priors[i][j]);
          fflush(fp);
        }
      }
    }
  }
  fprintf(fp,
          "--------------------------------"
          "--------------------------------\n");
  MatrixFree(&m_cov);
  VectorFree(&v_means);
  return (NO_ERROR);
}

////////////////////////////////////////////////////////////////
// transform from template -> node
////////////////////////////////////////////////////////////////
int GCAvoxelToNodeReal(GCA *gca, MRI *mri, double xv, double yv, double zv, double *pxn, double *pyn, double *pzn)
{
//               i_to_r
//        mri    ----->    RAS
//         |                |
//         | we need        | identity
//         V                V
//        node   <-----    RAS
//               r_to_i
  AffineMatrix voxelToNode, nodeFromRAS;
  SetAffineMatrix(&nodeFromRAS, gca->node_r_to_i__);
  AffineMM(&voxelToNode, &nodeFromRAS, mri->i_to_r__);

  AffineVector vv, nv;
  SetAffineVector(&vv, xv, yv, zv);
  AffineMV(&nv, &voxelToNode, &vv);

  float nxf, nyf, nzf;
  GetAffineVector(&nv, &nxf, &nyf, &nzf);
  *pxn = nxf;
  *pyn = nyf;
  *pzn = nzf;


  return NO_ERROR;
}

int GCAvoxelToNode(GCA *gca, MRI *mri, int xv, int yv, int zv, int *pxn, int *pyn, int *pzn)
{
  int ixn, iyn, izn, xp, yp, zp;
  int errCode = NO_ERROR;

  GCAvoxelToPrior(gca, mri, xv, yv, zv, &xp, &yp, &zp);
  GCApriorToNode(gca, xp, yp, zp, &ixn, &iyn, &izn);

  // xn, yn, zn are double.  we use voxel center as integer
  // if outofbounds, tell it
  errCode = boundsCheck(&ixn, &iyn, &izn, gca->mri_node__);
  //
  *pxn = ixn;
  *pyn = iyn;
  *pzn = izn;

  return errCode;
}

// use rounding to compute the index, not truncation
int GCAvoxelToNodeNearest(GCA *gca, MRI *mri, int xv, int yv, int zv, int *pxn, int *pyn, int *pzn)
{
  double xn, yn, zn;
  int ixn, iyn, izn;
  int errCode = NO_ERROR;

  GCAvoxelToNodeReal(gca, mri, xv, yv, zv, &xn, &yn, &zn);

  // xn, yn, zn are double.  we use voxel center as integer
  ixn = (int)nint(xn);
  iyn = (int)nint(yn);
  izn = (int)nint(zn);
  // if outofbounds, tell it
  errCode = boundsCheck(&ixn, &iyn, &izn, gca->mri_node__);
  //
  *pxn = ixn;
  *pyn = iyn;
  *pzn = izn;

  return errCode;
}

////////////////////////////////////////////////////////////////////
// transform from template -> prior
////////////////////////////////////////////////////////////////////
int GCAvoxelToPriorReal(const GCA *gca,
                        const MRI *mri,
                        const double xv,
                        const double yv,
                        const double zv,
                        double *pxp,
                        double *pyp,
                        double *pzp)
{
// printf("GCAvoxelToPriorReal1\n");
  AffineMatrix voxelToPrior;

  AffineMM(&voxelToPrior, gca->prior_r_to_i__, mri->i_to_r__);

  AffineVector vv, pv;

  SetAffineVector(&vv, xv, yv, zv);
  AffineMV(&pv, &voxelToPrior, &vv);

  float pxf, pyf, pzf;
  GetAffineVector(&pv, &pxf, &pyf, &pzf);
  *pxp = pxf;
  *pyp = pyf;
  *pzp = pzf;


  return NO_ERROR;
}

int GCAvoxelToPrior(
    const GCA *gca, const MRI *mri, const int xv, const int yv, const int zv, int *pxp, int *pyp, int *pzp)
{
  double xp, yp, zp;
  int ixp, iyp, izp;
  int errCode = NO_ERROR;

  /*!
    @bugs gca would be const if its tmp__ variable wasn't used as workspace in GCAvoxelToPriorReal
  */

  GCAvoxelToPriorReal(gca, mri, xv, yv, zv, &xp, &yp, &zp);

#if USE_ROUND
  ixp = (int)nint(xp);
  iyp = (int)nint(yp);
  izp = (int)nint(zp);
  if (ixp >= gca->prior_width) ixp = gca->prior_width - 1;
  if (iyp >= gca->prior_height) iyp = gca->prior_height - 1;
  if (izp >= gca->prior_depth) izp = gca->prior_depth - 1;
#else
  ixp = (int)floor(xp);
  iyp = (int)floor(yp);
  izp = (int)floor(zp);
#endif
  // bound check
  // if outofbounds, tell it
  errCode = boundsCheck(&ixp, &iyp, &izp, gca->mri_prior__);
  //
  *pxp = ixp;
  *pyp = iyp;
  *pzp = izp;

  return errCode;
}

// use rounding to compute the index, not truncation
int GCAvoxelToPriorNearest(GCA *gca, MRI *mri, int xv, int yv, int zv, int *pxp, int *pyp, int *pzp)
{
  double xp, yp, zp;
  int ixp, iyp, izp;
  int errCode = NO_ERROR;

  GCAvoxelToPriorReal(gca, mri, xv, yv, zv, &xp, &yp, &zp);

  ixp = (int)nint(xp);
  iyp = (int)nint(yp);
  izp = (int)nint(zp);
  // bound check
  // if outofbounds, tell it
  errCode = boundsCheck(&ixp, &iyp, &izp, gca->mri_prior__);
  //
  *pxp = ixp;
  *pyp = iyp;
  *pzp = izp;

  return errCode;
}

/////////////////////////////////////////////////////////////////////
// transform node->template
/////////////////////////////////////////////////////////////////////
int GCAnodeToVoxelReal(GCA *gca, MRI *mri, double xn, double yn, double zn, double *pxv, double *pyv, double *pzv)
{
  //               r_to_i
  //        mri    <-----    RAS
  //         ^                ^
  //         | we need        | identity
  //         |                |
  //        node   ----->    RAS
  //               i_to_r
  MATRIX *rasFromNode = gca->node_i_to_r__;
  //  extract_i_to_r(gca->mri_node__);
  MATRIX *voxelFromRAS = mri->r_to_i__;  //  extract_r_to_i(mri);
  MATRIX *nodeToVoxel = gca->tmp__;
  MatrixMultiply(voxelFromRAS, rasFromNode, gca->tmp__);

  TransformWithMatrix(nodeToVoxel, xn, yn, zn, pxv, pyv, pzv);

  // MatrixFree(&rasFromNode);
  // MatrixFree(&voxelFromRAS);
  // MatrixFree(&nodeToVoxel);

  return NO_ERROR;
}

int GCAnodeToVoxel(GCA *gca, MRI *mri, int xn, int yn, int zn, int *pxv, int *pyv, int *pzv)
{
  double xv, yv, zv;
  int ixv, iyv, izv;
  int errCode = NO_ERROR;

  GCAnodeToVoxelReal(gca, mri, xn, yn, zn, &xv, &yv, &zv);

  // addition makes overall error smaller
  //             0 picks  1 out of 0, 1, 2, 3 possible choices
  // without it, 0 pickes 0 out of 0, 1, 2, 3 possible choices
  // if you used float, then 0 picks 1.5.
  ixv = (int)floor(xv + gca->node_spacing / 2.);
  iyv = (int)floor(yv + gca->node_spacing / 2.);
  izv = (int)floor(zv + gca->node_spacing / 2.);

  // bound check
  // if outofbounds, tell it
  errCode = boundsCheck(&ixv, &iyv, &izv, mri);
  //

  *pxv = ixv;
  *pyv = iyv;
  *pzv = izv;
  return errCode;
}

//////////////////////////////////////////////////////////////////////
// transform from prior-> template
//////////////////////////////////////////////////////////////////////
int GCApriorToVoxelReal(GCA *gca, MRI *mri, double xp, double yp, double zp, double *pxv, double *pyv, double *pzv)
{
  int tid;
  MATRIX *rasFromPrior = gca->prior_i_to_r__;
  // extract_i_to_r(gca->mri_prior__);
  MATRIX *voxelFromRAS = mri->r_to_i__;  // extract_r_to_i(mri);
  static MATRIX *priorToVoxel[_MAX_FS_THREADS];
#ifdef HAVE_OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif
  priorToVoxel[tid] = MatrixMultiply(voxelFromRAS, rasFromPrior, priorToVoxel[tid]);

  // TransformWithMatrix(priorToVoxel, xp, yp, zp, pxv, pyv, pzv);
  TransformWithMatrix(priorToVoxel[tid], xp, yp, zp, pxv, pyv, pzv);

  // MatrixFree(&rasFromPrior);
  // MatrixFree(&voxelFromRAS);
  // MatrixFree(&priorToVoxel);

  return NO_ERROR;
}

int GCApriorToVoxel(GCA *gca, MRI *mri, int xp, int yp, int zp, int *pxv, int *pyv, int *pzv)
{
  double xv, yv, zv;
  int ixv, iyv, izv;
  int errCode = NO_ERROR;

  GCApriorToVoxelReal(gca, mri, xp, yp, zp, &xv, &yv, &zv);
// addition makes overall error smaller
// without it, 0 pickes 0 out of 0, 1 possible choices
  ixv = (int)floor(xv + gca->prior_spacing / 2.);
  iyv = (int)floor(yv + gca->prior_spacing / 2.);
  izv = (int)floor(zv + gca->prior_spacing / 2.);
  // bound check
  // if outofbounds, tell it
  errCode = boundsCheck(&ixv, &iyv, &izv, mri);
  //
  *pxv = ixv;
  *pyv = iyv;
  *pzv = izv;

  errCode = NO_ERROR;

  return errCode;
}

///////////////////////////////////////////////////////////////////////
// transform from source -> template space -> prior
//////////////////////////////////////////////////////////////////////
int GCAsourceVoxelToPrior(
    const GCA *gca, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv, int *pxp, int *pyp, int *pzp)
{
  float xt = 0, yt = 0, zt = 0;
  double xrt, yrt, zrt;
  int retval = NO_ERROR;
  LTA *lta;

  if (transform->type != MORPH_3D_TYPE) {
    if (transform->type == LINEAR_VOX_TO_VOX) {
      lta = (LTA *)transform->xform;
      // transform point to talairach volume point
      TransformWithMatrix(lta->xforms[0].m_L, xv, yv, zv, &xrt, &yrt, &zrt);
      /*!
      @BUGS
      This is strange - a downcast from double to float,
      when the next step for these variables is a
      call to GCAvoxelToPrior, which will immediately
      convert them to int
      */
      // we should change everything to use double, but at the moment they don't, which is why the cast is here (BRF)
      xt = xrt;
      yt = yrt;
      zt = zrt;
      // TransformSample(transform, xv, yv, zv, &xt, &yt, &zt) ;
    }
    else
      ErrorExit(ERROR_BADPARM, "GCAsourceVoxelToPrior: needs vox-to-vox transform");
  }
  else  // morph 3d type can go directly from source to template
  {
    TransformSample(transform, xv, yv, zv, &xt, &yt, &zt);
  }
// get the position in gca from talairach volume
/*!
  @BUGS
  Combine with note above - xt, yt and zt will be converted to
  integers on this call
*/
#if USE_ROUND
  {
    double xpf, ypf, zpf;
    GCAvoxelToPriorReal(gca, gca->mri_tal__, xt, yt, zt, &xpf, &ypf, &zpf);
    *pxp = nint(xpf);
    *pyp = nint(ypf);
    *pzp = nint(zpf);
  }
#else
  GCAvoxelToPrior(gca, gca->mri_tal__, xt, yt, zt, pxp, pyp, pzp);
#endif

  if (*pxp < 0 || *pyp < 0 || *pzp < 0 || *pxp >= gca->prior_width || *pyp >= gca->prior_height ||
      *pzp >= gca->prior_depth) {
    retval = (ERROR_BADPARM);
  }
  if (*pxp < 0) {
    *pxp = 0;
  }
  if (*pyp < 0) {
    *pyp = 0;
  }
  if (*pzp < 0) {
    *pzp = 0;
  }
  if (*pxp >= gca->prior_width) {
    *pxp = gca->prior_width - 1;
  }
  if (*pyp >= gca->prior_height) {
    *pyp = gca->prior_height - 1;
  }
  if (*pzp >= gca->prior_depth) {
    *pzp = gca->prior_depth - 1;
  }
  return (retval);
}

int GCAsourceFloatVoxelToPrior(
    GCA *gca, MRI *mri, TRANSFORM *transform, float xv, float yv, float zv, int *pxp, int *pyp, int *pzp)
{
  float xt = 0, yt = 0, zt = 0;
  double xrt, yrt, zrt, xrp, yrp, zrp;

  LTA *lta;
  if (transform->type != MORPH_3D_TYPE) {
    if (transform->type == LINEAR_VOX_TO_VOX) {
      lta = (LTA *)transform->xform;
      // transform point to talairach volume point
      TransformWithMatrix(lta->xforms[0].m_L, xv, yv, zv, &xrt, &yrt, &zrt);
      xt = xrt;
      yt = yrt;
      zt = zrt;
      // TransformSample(transform, xv, yv, zv, &xt, &yt, &zt) ;
    }
    else
      ErrorExit(ERROR_BADPARM, "GCAsourceVoxelToPrior: needs vox-to-vox transform");
  }
  else  // morph 3d type can go directly from source to template
  {
    TransformSample(transform, xv, yv, zv, &xt, &yt, &zt);
  }
  // get the position in gca from talairach volume
  GCAvoxelToPriorReal(gca, gca->mri_tal__, xt, yt, zt, &xrp, &yrp, &zrp);
  *pxp = nint(xrp);
  *pyp = nint(yrp);
  *pzp = nint(zrp);
  if (*pxp < 0 || *pyp < 0 || *pzp < 0 || *pxp >= gca->prior_width || *pyp >= gca->prior_height ||
      *pzp >= gca->prior_depth) {
    return (ERROR_BADPARM);
  }
  return (NO_ERROR);
}

int GCAsourceFloatVoxelToPriorReal(
    GCA *gca, MRI *mri, TRANSFORM *transform, float xv, float yv, float zv, double *pxp, double *pyp, double *pzp)
{
  float xt = 0, yt = 0, zt = 0;
  double xrt, yrt, zrt, xrp, yrp, zrp;

  LTA *lta;
  if (transform->type != MORPH_3D_TYPE) {
    if (transform->type == LINEAR_VOX_TO_VOX) {
      lta = (LTA *)transform->xform;
      // transform point to talairach volume point
      TransformWithMatrix(lta->xforms[0].m_L, xv, yv, zv, &xrt, &yrt, &zrt);
      xt = xrt;
      yt = yrt;
      zt = zrt;
      // TransformSample(transform, xv, yv, zv, &xt, &yt, &zt) ;
    }
    else
      ErrorExit(ERROR_BADPARM, "GCAsourceVoxelToPrior: needs vox-to-vox transform");
  }
  else  // morph 3d type can go directly from source to template
  {
    TransformSampleReal2(transform, xv, yv, zv, &xt, &yt, &zt);
  }
  // get the position in gca from talairach volume
  GCAvoxelToPriorReal(gca, gca->mri_tal__, xt, yt, zt, &xrp, &yrp, &zrp);
  *pxp = xrp;
  *pyp = yrp;
  *pzp = zrp;
  if (*pxp < 0 || *pyp < 0 || *pzp < 0 || *pxp >= gca->prior_width || *pyp >= gca->prior_height ||
      *pzp >= gca->prior_depth)
    return (ERROR_BADPARM);
  return (NO_ERROR);
}

/////////////////////////////////////////////////////////////////////
// transform from source -> template space -> prior
/////////////////////////////////////////////////////////////////////
int GCAsourceVoxelToNode(
    const GCA *gca, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv, int *pxn, int *pyn, int *pzn)
{
  int xp, yp, zp;

  GCAsourceVoxelToPrior(gca, mri, transform, xv, yv, zv, &xp, &yp, &zp);
  GCApriorToNode(gca, xp, yp, zp, pxn, pyn, pzn);
  return (NO_ERROR);
}

//////////////////////////////////////////////////////////////////////
// get the matrix that transforms prior coordinates into source voxel coords
////////////////////////////////////////////////////////////////////
MATRIX *GCAgetPriorToSourceVoxelMatrix(GCA *gca, const MRI *mri, TRANSFORM *transform)
{
  MATRIX *rasFromPrior = gca->prior_i_to_r__;
  // extract_i_to_r(gca->mri_prior__);
  MATRIX *voxelFromRAS = gca->mri_tal__->r_to_i__;  // extract_r_to_i(mri);
  MATRIX *m_prior2voxel, *m_prior2source_voxel, *m_tmp;
  LTA *lta = (LTA *)transform->xform;

  // go to the template voxel position
  //  GCApriorToVoxelReal(gca, gca->mri_tal__, xp, yp, zp, &xt, &yt, &zt);
  m_prior2voxel = MatrixMultiply(voxelFromRAS, rasFromPrior, NULL);
  if (transform->type != LINEAR_VOX_TO_VOX)  // from src to talairach volume
    ErrorExit(ERROR_BADPARM, "GCAgetPriorToSourceVoxelMatrix: needs vox-to-vox transform");
  m_tmp = MatrixMultiply(lta->inv_xforms[0].m_L, m_prior2voxel, NULL);
  m_prior2source_voxel = m_tmp;
  MatrixFree(&m_prior2voxel);
  return (m_prior2source_voxel);
}
//////////////////////////////////////////////////////////////////////
// transform from node->template->source
////////////////////////////////////////////////////////////////////
int GCApriorToSourceVoxelFloat(GCA *gca,
                               const MRI *mri,
                               TRANSFORM *transform,
                               const int xp,
                               const int yp,
                               const int zp,
                               float *pxv,
                               float *pyv,
                               float *pzv)
{
  int width, height, depth;
  double xt = 0, yt = 0, zt = 0;
  float xv = 0, yv = 0, zv = 0;
  double xc, yc, zc;
  int errCode = NO_ERROR;
  LTA *lta;

  if (transform->type == LINEAR_VOX_TO_VOX) {
    MATRIX *m;
    VECTOR *v_src, *v_dst;

    v_src = VectorAlloc(4, MATRIX_REAL);
    v_dst = VectorAlloc(4, MATRIX_REAL);
    *MATRIX_RELT(v_src, 4, 1) = 1.0;
    *MATRIX_RELT(v_dst, 4, 1) = 1.0;
    m = GCAgetPriorToSourceVoxelMatrix(gca, mri, transform);
    V3_X(v_src) = xp;
    V3_Y(v_src) = yp;
    V3_Z(v_src) = zp;
    MatrixMultiply(m, v_src, v_dst);
    xc = V3_X(v_dst);
    yc = V3_Y(v_dst);
    zc = V3_Z(v_dst);
    MatrixFree(&m);
    VectorFree(&v_src);
    VectorFree(&v_dst);
    width = mri->width;
    height = mri->height;
    depth = mri->depth;
    if (xc < 0)
      errCode = ERROR_BADPARM;
    else if (yc < 0)
      errCode = ERROR_BADPARM;
    else if (zc < 0)
      errCode = ERROR_BADPARM;
    else if (xc > (width - 1))
      errCode = ERROR_BADPARM;
    else if (yc > (height - 1))
      errCode = ERROR_BADPARM;
    else if (zc > (depth - 1))
      errCode = ERROR_BADPARM;
    *pxv = xc;
    *pyv = yc;
    *pzv = zc;
    return errCode;
  }

  // go to the template voxel position
  GCApriorToVoxelReal(gca, gca->mri_tal__, xp, yp, zp, &xt, &yt, &zt);
  // got the point in gca->mri_tal__ position
  if (transform->type != MORPH_3D_TYPE) {
    if (transform->type == LINEAR_VOX_TO_VOX)  // from src to talairach volume
    {
      width = mri->width;
      height = mri->height;
      depth = mri->depth;
      lta = (LTA *)transform->xform;
      // get the talairach to orig
      TransformWithMatrix(lta->inv_xforms[0].m_L, xt, yt, zt, &xc, &yc, &zc);
      // TransformSampleInverse(transform, xt, yt, zt, &xc, &yc, &zc);
      if (xc < 0) {
        errCode = ERROR_BADPARM;
      }
      else if (yc < 0) {
        errCode = ERROR_BADPARM;
      }
      else if (zc < 0) {
        errCode = ERROR_BADPARM;
      }
      else if (xc > (width - 1)) {
        errCode = ERROR_BADPARM;
      }
      else if (yc > (height - 1)) {
        errCode = ERROR_BADPARM;
      }
      else if (zc > (depth - 1)) {
        errCode = ERROR_BADPARM;
      }
      xv = xc;
      yv = yc;
      zv = zc;
    }
    else
      ErrorExit(ERROR_BADPARM, "GCApriorToSourceVoxelFloat: needs vox-to-vox transform");
  }
  else  // go directly from template to source
  {
    TransformSampleInverse(transform, xt, yt, zt, &xv, &yv, &zv);
  }
  *pxv = xv;
  *pyv = yv;
  *pzv = zv;
  return errCode;
}

int GCAnodeToSourceVoxelFloat(
    GCA *gca, MRI *mri, TRANSFORM *transform, int xn, int yn, int zn, float *pxv, float *pyv, float *pzv)
{
  int width, height, depth;
  double xt, yt, zt;
  float xv, yv, zv;
  double xc, yc, zc;
  int errCode = NO_ERROR;
  LTA *lta;
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  // get template voxel position
  GCAnodeToVoxelReal(gca, gca->mri_tal__, xn, yn, zn, &xt, &yt, &zt);
  if (transform->type != MORPH_3D_TYPE) {
    lta = (LTA *)transform->xform;
    // get the talairach to orig
    TransformWithMatrix(lta->inv_xforms[0].m_L, xt, yt, zt, &xc, &yc, &zc);
    // TransformSampleInverse(transform, xt, yt, zt, &xc, &yc, &zc);
    if (xc < 0) {
      errCode = ERROR_BADPARM;
    }
    else if (yc < 0) {
      errCode = ERROR_BADPARM;
    }
    else if (zc < 0) {
      errCode = ERROR_BADPARM;
    }
    else if (xc > (width - 1)) {
      errCode = ERROR_BADPARM;
    }
    else if (yc > (height - 1)) {
      errCode = ERROR_BADPARM;
    }
    else if (zc > (depth - 1)) {
      errCode = ERROR_BADPARM;
    }
    xv = xc;
    yv = yc;
    zv = zc;
  }
  else  // template to source
  {
    TransformSampleInverse(transform, xt, yt, zt, &xv, &yv, &zv);
  }
  *pxv = xv;
  *pyv = yv;
  *pzv = zv;
  return errCode;
}

int GCAnodeToSourceVoxel(GCA *gca, MRI *mri, TRANSFORM *transform, int xn, int yn, int zn, int *pxv, int *pyv, int *pzv)
{
  float xf, yf, zf;
  int errCode = NO_ERROR;
  errCode = GCAnodeToSourceVoxelFloat(gca, mri, transform, xn, yn, zn, &xf, &yf, &zf);
  if (xf < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (yf < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (zf < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (xf > (mri->width - 1)) {
    errCode = ERROR_BADPARM;
  }
  else if (yf > (mri->height - 1)) {
    errCode = ERROR_BADPARM;
  }
  else if (zf > (mri->depth - 1)) {
    errCode = ERROR_BADPARM;
  }

  *pxv = nint(xf);
  *pyv = nint(yf);
  *pzv = nint(zf);
  return errCode;
}

int GCApriorToSourceVoxel(GCA *gca,
                          const MRI *mri,
                          TRANSFORM *transform,
                          const int xp,
                          const int yp,
                          const int zp,
                          int *pxv,
                          int *pyv,
                          int *pzv)
{
  float xf, yf, zf;
  int errCode = NO_ERROR;

  errCode = GCApriorToSourceVoxelFloat(gca, mri, transform, xp, yp, zp, &xf, &yf, &zf);
  if (xf < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (yf < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (zf < 0) {
    errCode = ERROR_BADPARM;
  }
  else if (xf > (mri->width - 1)) {
    errCode = ERROR_BADPARM;
  }
  else if (yf > (mri->height - 1)) {
    errCode = ERROR_BADPARM;
  }
  else if (zf > (mri->depth - 1)) {
    errCode = ERROR_BADPARM;
  }

  *pxv = nint(xf);
  *pyv = nint(yf);
  *pzv = nint(zf);
  return errCode;
}

GCA *GCAalloc(int ninputs, float prior_spacing, float node_spacing, int width, int height, int depth, int flags)
{
  int max_labels;

  if (flags & GCA_NO_GCS) {
    max_labels = 0;
  }
  else {
    max_labels = DEFAULT_MAX_LABELS_PER_GCAN;
  }
  return (gcaAllocMax(ninputs, prior_spacing, node_spacing, width, height, depth, max_labels, flags));
}

GCA *gcaAllocMax(
    int ninputs, float prior_spacing, float node_spacing, int width, int height, int depth, int max_labels, int flags)
{
  GCA *gca;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  int x, y, z;

  gca = (GCA *)calloc(1, sizeof(GCA));
  if (!gca) {
    ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate struct");
  }

  gca->ninputs = ninputs;
  gca->prior_spacing = prior_spacing;
  gca->node_spacing = node_spacing;
  gca->type = GCA_UNKNOWN;  // mark it as unknown

  // setup default direction cosines
  gca->x_r = -1;
  gca->y_r = 0;
  gca->z_r = 0;
  gca->c_r = 0;
  gca->x_a = 0;
  gca->y_a = 0;
  gca->z_a = 1;
  gca->c_a = 0;
  gca->x_s = 0;
  gca->y_s = -1;
  gca->z_s = 0;
  gca->c_s = 0;
  gca->xsize = 1;
  gca->ysize = 1;
  gca->zsize = 1;
  //
  gca->width = width;
  gca->height = height;
  gca->depth = depth;

  /* ceil gives crazy results, I don't know why */
  gca->node_width = (int)(((float)width / node_spacing) + .99);
  gca->node_height = (int)((float)height / node_spacing + .99);
  gca->node_depth = (int)(((float)depth / node_spacing) + .99);
  gca->prior_width = (int)(((float)width / prior_spacing) + .99);
  gca->prior_height = (int)((float)height / prior_spacing + .99);
  gca->prior_depth = (int)(((float)depth / prior_spacing) + .99);
  gca->flags = flags;
  if (max_labels >= 0) {
    gca->nodes = (GCA_NODE ***)calloc(gca->node_width, sizeof(GCA_NODE **));
    if (!gca->nodes) {
      ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate nodes");
    }

    // setting vlaues gca->nodes volume
    for (x = 0; x < gca->node_width; x++) {
      gca->nodes[x] = (GCA_NODE **)calloc(gca->node_height, sizeof(GCA_NODE *));
      if (!gca->nodes[x]) {
        ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate %dth **", x);
      }

      for (y = 0; y < gca->node_height; y++) {
        gca->nodes[x][y] = (GCA_NODE *)calloc(gca->node_depth, sizeof(GCA_NODE));
        if (!gca->nodes[x][y]) ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate %d,%dth *", x, y);
        for (z = 0; z < gca->node_depth; z++) {
          gcan = &gca->nodes[x][y][z];
          // set max_labels
          gcan->max_labels = max_labels;

          if (max_labels > 0) {
            /* allocate new ones */
            gcan->gcs = alloc_gcs(gcan->max_labels, flags, gca->ninputs);
            if (!gcan->gcs) ErrorExit(ERROR_NOMEMORY, "GCANalloc: couldn't allocate gcs to %d", gcan->max_labels);
            // allocate label storage up to max_labels
            gcan->labels = (unsigned short *)calloc(gcan->max_labels, sizeof(unsigned short));
            if (!gcan->labels) ErrorExit(ERROR_NOMEMORY, "GCANalloc: couldn't allocate labels to %d", gcan->max_labels);
          }
        }
      }
    }

    gca->priors = (GCA_PRIOR ***)calloc(gca->prior_width, sizeof(GCA_PRIOR **));
    if (!gca->priors) {
      ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate priors");
    }
    // setting values to gca->prior volume
    for (x = 0; x < gca->prior_width; x++) {
      gca->priors[x] = (GCA_PRIOR **)calloc(gca->prior_height, sizeof(GCA_PRIOR *));
      if (!gca->priors[x]) {
        ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate %dth **", x);
      }

      for (y = 0; y < gca->prior_height; y++) {
        gca->priors[x][y] = (GCA_PRIOR *)calloc(gca->prior_depth, sizeof(GCA_PRIOR));
        if (!gca->priors[x][y]) ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate %d,%dth *", x, y);
        for (z = 0; z < gca->prior_depth; z++) {
          gcap = &gca->priors[x][y][z];
          if (gcap == NULL) {
            continue;
          }
          gcap->max_labels = max_labels;

          if (max_labels > 0) {
            /* allocate new ones */
            // allocate label space
            gcap->labels = (unsigned short *)calloc(max_labels, sizeof(unsigned short));
            if (!gcap->labels) ErrorExit(ERROR_NOMEMORY, "GCANalloc: couldn't allocate labels to %d", gcap->max_labels);
            // create prior space
            gcap->priors = (float *)calloc(max_labels, sizeof(float));
            if (!gcap->priors) ErrorExit(ERROR_NOMEMORY, "GCANalloc: couldn't allocate priors to %d", max_labels);
          }
        }
      }
    }
  }
  // setup
  gca->mri_node__ = 0;
  gca->mri_prior__ = 0;
  gca->mri_tal__ = 0;
  // initialize
  gca->node_i_to_r__ = NULL;
  gca->node_r_to_i__ = NULL;
  gca->prior_i_to_r__ = NULL;
  gca->prior_r_to_i__ = NULL;
  gca->tal_i_to_r__ = gca->tal_r_to_i__ = 0;

  GCAsetup(gca);

  return (gca);
}

int GCAfree(GCA **pgca)
{
  GCA *gca;
  int x, y, z;

  gca = *pgca;
  *pgca = NULL;

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        GCANfree(&gca->nodes[x][y][z], gca->ninputs);
      }
      free(gca->nodes[x][y]);
    }
    free(gca->nodes[x]);
  }
  free(gca->nodes);

  for (x = 0; x < gca->prior_width; x++) {
    for (y = 0; y < gca->prior_height; y++) {
      for (z = 0; z < gca->prior_depth; z++) {
        free(gca->priors[x][y][z].labels);
        free(gca->priors[x][y][z].priors);
      }
      free(gca->priors[x][y]);
    }
    free(gca->priors[x]);
  }

  free(gca->priors);
  GCAcleanup(gca);

  free(gca);

  return (NO_ERROR);
}

int GCANfree(GCA_NODE *gcan, int ninputs)
{
  if (gcan->nlabels) {
    free(gcan->labels);
    free_gcs(gcan->gcs, gcan->nlabels, ninputs);
  }
  return (NO_ERROR);
}

int GCAPfree(GCA_PRIOR *gcap)
{
  if (gcap->nlabels) {
    free(gcap->labels);
    free(gcap->priors);
  }
  return (NO_ERROR);
}

void PrintInfoOnLabels(GCA *gca, int label, int xn, int yn, int zn, int xp, int yp, int zp, int x, int y, int z)
{
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  GC1D *gc;
  int i;
  // using node to find the label
  gc = GCAfindGC(gca, xn, yn, zn, label);
  if (gc) {
    gcan = &gca->nodes[xn][yn][zn];
    fprintf(stdout, "\n Node (%3d, %3d, %3d) pos (%3d, %3d, %3d) label=%d, labels:", xn, yn, zn, x, y, z, label);
    for (i = 0; i < gcan->nlabels; ++i) {
      fprintf(stdout, "%4d ", gcan->labels[i]);
    }
    fprintf(stdout, "\n");
    gcap = &gca->priors[xp][yp][zp];
    if (gcap == NULL) {
      return;
    }
    fprintf(stdout, "Prior (%3d, %3d, %3d) pos (%3d, %3d, %3d) label=%d\n", xp, yp, zp, x, y, z, label);
    fprintf(stdout, "prior label histogram  (label):");
    for (i = 0; i < gcap->nlabels; ++i) {
      fprintf(stdout, "%4d ", gcap->labels[i]);
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "                     :(counts):");
    for (i = 0; i < gcap->nlabels; ++i) {
      fprintf(stdout, "%4.f ", gcap->priors[i]);
    }
    fprintf(stdout, "\n");
    fprintf(stdout, "mean: ");
    for (i = 0; i < gca->ninputs; i++) {
      fprintf(stdout, "%2.1f ", gc->means[i] / gc->ntraining);
    }
    fprintf(stdout, "\n");
  }
  fflush(stdout);
}

int GCAtrainCovariances(GCA *gca, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform)
{
  int x, y, z, width, height, depth, label, xn, yn, zn;
  float vals[MAX_GCA_INPUTS];
  int xp, yp, zp, holes_filled = 0;
  MRI *mri_mapped;

  /* convert transform to voxel coordinates */

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  mri_mapped = MRIalloc(gca->node_width, gca->node_height, gca->node_depth, MRI_UCHAR);
  width = mri_labels->width;
  height = mri_labels->height;
  depth = mri_labels->depth;
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        // get the segmented value
        label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));
        // get all input volume values at this point
        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

        // src -> talairach -> node
        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          if (!GCAsourceVoxelToPrior(gca, mri_inputs, transform, x, y, z, &xp, &yp, &zp)) {
            ///////////////// debug code ////////////////////////////
            if ((xp == Gxp && yp == Gyp && zp == Gzp) && (Ggca_label < 0 || Ggca_label == label)) {
              printf(
                  "src (%d, %d, %d), "
                  "prior (%d, %d, %d), "
                  "node (%d, %d, %d), label = %d\n",
                  x,
                  y,
                  z,
                  xp,
                  yp,
                  zp,
                  xn,
                  yn,
                  zn,
                  label);
            }
            ////////////////////////////////////////////////////////
            // update the value
            MRIsetVoxVal(mri_mapped, xn, yn, zn, 0, 1);
            GCAupdateNodeCovariance(gca, mri_inputs, xn, yn, zn, vals, label);

            //////////////debug code ////////////////////////////
            if (xn == Gxn && yn == Gyn && zn == Gzn) {
              fprintf(stdout, "Train Covariance\n");
              PrintInfoOnLabels(gca, label, xn, yn, zn, xp, yp, zp, x, y, z);
            }
            if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z && (label == Ggca_label || Ggca_label < 0)) {
              GC1D *gc;
              int i, nsamples;
              MATRIX *m;
              gc = GCAfindGC(gca, xn, yn, zn, label);
              if (gc) {
                nsamples = gc->ntraining - gc->n_just_priors;
                /* for no-intensity training */
                if (nsamples < 1) {
                  nsamples = 1;
                }
                printf("voxel(%d,%d,%d) = ", x, y, z);
                for (i = 0; i < gca->ninputs; i++) {
                  printf("%d ", nint(vals[i]));
                }

                printf(
                    " --> node(%d,%d,%d), "
                    "label %s (%d), mean ",
                    xn,
                    yn,
                    zn,
                    cma_label_to_name(label),
                    label);
                for (i = 0; i < gca->ninputs; i++) {
                  printf("%2.1f ", gc->means[i]);
                }
                printf("\ncovariances (det=%f):\n",
                       covariance_determinant(gc, gca->ninputs) / pow((double)nsamples, (double)gca->ninputs));
                m = load_covariance_matrix(gc, NULL, gca->ninputs);
                MatrixScalarMul(m, 1.0 / (nsamples), m);
                MatrixPrint(stdout, m);
                MatrixFree(&m);
              }
            }
            /////////////////////////////////////////////////
          }
        }  // if (!GCA...)
      }
    }
  }

  for (holes_filled = xn = 0; xn < gca->node_width; xn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (zn = 0; zn < gca->node_depth; zn++) {
        if (xn == Gx && yn == Gy && zn == Gz) {
          DiagBreak();
        }
        if (MRIgetVoxVal(mri_mapped, xn, yn, zn, 0) > 0) {
          continue;
        }
        if (!GCAnodeToSourceVoxel(gca, mri_inputs, transform, xn, yn, zn, &x, &y, &z)) {
          GC1D *gc;
          label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));
          gc = GCAfindGC(gca, xn, yn, zn, label);
          if (gc == NULL)  // label doesn't exist at this position
          {
            if (xn == Gx && yn == Gy && zn == Gz)
              printf("filling node hole at (%d, %d, %d) %s\n", xn, yn, zn, cma_label_to_name(label));
          }
          else {
            load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
            GCAupdateNodeCovariance(gca, mri_inputs, xn, yn, zn, vals, label);
          }
          holes_filled++;
        }
      }
    }
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_mapped, "m.mgz");
  }
  if (holes_filled > 0) {
    printf("%d prior holes filled\n", holes_filled);
  }
  MRIfree(&mri_mapped);
  fflush(stdout);
  return (NO_ERROR);
}
#include <unistd.h>

int GCAtrain(GCA *gca, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform, GCA *gca_prune, int noint)
{
  int x, y, z, width, height, depth, label, xn, yn, zn, holes_filled,
      /*i, xnbr, ynbr, znbr, xn_nbr, yn_nbr, zn_nbr,*/ xp, yp, zp;
  float vals[MAX_GCA_INPUTS];
  static int first_time = 1;
  FILE *logfp = NULL;
  static int logging = 0;
  GCA_PRIOR *gcap;
  GCA_NODE *gcan;
  MRI *mri_mapped;

  gca->total_training++;
  mri_mapped = MRIalloc(gca->prior_width, gca->prior_height, gca->prior_depth, MRI_UCHAR);
  if (first_time) {
    first_time = 0;
    logging = getenv("GCA_LOG") != NULL;
    if (logging) {
      printf("logging image intensities to GCA log file 'gca*.log'\n");
      for (label = 0; label <= MAX_CMA_LABEL; label++) {
        char fname[STRLEN];
        sprintf(fname, "gca%d.log", label);
        if (FileExists(fname)) {
          unlink(fname);
        }
      }
    }
  }

  /* convert transform to voxel coordinates */

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  // segmented volume
  width = mri_labels->width;
  height = mri_labels->height;
  depth = mri_labels->depth;
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        /// debugging /////////////////////////////////////
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }

        ///////////////////////////////////////////////////

        // get the segmented voxel label
        label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));
        if (label > gca->max_label) gca->max_label = label;
        // get all the values to vals[] in inputs at this point
        // mri_inputs are T1, PD etc.
        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
        // segmented volume->talairach volume->node
        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn))
          if (!GCAsourceVoxelToPrior(gca, mri_inputs, transform, x, y, z, &xp, &yp, &zp)) {
            // got node point (xn, yn. zn) and
            // prior point (xp, yp, zp) for
            // this label volume point (x, y, z)
            // update the value at this prior point
            if (xp == Gxp && yp == Gyp && zp == Gzp) {
              DiagBreak();
              if (Ggca_label < 0 || Ggca_label == label) {
                printf(
                    "src (%d, %d, %d), "
                    "prior (%d, %d, %d), "
                    "node (%d, %d, %d), label = %d\n",
                    x,
                    y,
                    z,
                    xp,
                    yp,
                    zp,
                    xn,
                    yn,
                    zn,
                    label);
              }
            }
            MRIsetVoxVal(mri_mapped, xp, yp, zp, 0, 1);
            GCAupdatePrior(gca, mri_inputs, xp, yp, zp, label);
            if ((GCAupdateNode(gca, mri_inputs, xn, yn, zn, vals, label, gca_prune, noint) == NO_ERROR) &&
                !(gca->flags & GCA_NO_MRF))
              //                 node        label point
              GCAupdateNodeGibbsPriors(gca, mri_labels, xn, yn, zn, x, y, z, label);

            /// debugging code //////////////////////////////////////
            {
              int n;
              GC1D *gc;
              gcap = &gca->priors[xp][yp][zp];
              if (gcap != NULL)
                for (n = 0; n < gcap->nlabels; n++) {
                  gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[n]);
                  if (gc == NULL) {
                    printf(
                        "(%d, %d, %d): gcap[%d][%d][%d]->labels[%d] ="
                        " %s - node (%d, %d, %d): no gc!\n",
                        x,
                        y,
                        z,
                        xp,
                        yp,
                        zp,
                        n,
                        cma_label_to_name(gcap->labels[n]),
                        xn,
                        yn,
                        zn);
                    DiagBreak();
                  }
                }
            }
            if (xn == Gxn && yn == Gyn && zn == Gzn) {
              PrintInfoOnLabels(gca, label, xn, yn, zn, xp, yp, zp, x, y, z);
            }
            if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z && (label == Ggca_label || (Ggca_label < 0)) &&
                (Ggca_nbr_label < 0)) {
              GC1D *gc;
              int i;
              gc = GCAfindGC(gca, xn, yn, zn, label);
              if (gc) {
                if (logging) {
                  char fname[STRLEN];
                  sprintf(fname, "gca%d.log", label);
                  logfp = fopen(fname, "a");
                }
                printf("voxel(%d,%d,%d) = ", x, y, z);
                for (i = 0; i < gca->ninputs; i++) {
                  printf("%2.1f ", (vals[i]));
                  if (logging) {
                    fprintf(logfp, "%2.1f ", (vals[i]));
                  }
                }

                printf(
                    " --> node(%d,%d,%d), "
                    "label %s (%d), mean ",
                    xn,
                    yn,
                    zn,
                    cma_label_to_name(label),
                    label);
                for (i = 0; i < gca->ninputs; i++) {
                  printf("%2.1f ", gc->means[i] / gc->ntraining);
                }
                printf("\n");
                gcan = &gca->nodes[xn][yn][zn];
                printf("   node labels:");
                for (i = 0; i < gcan->nlabels; ++i) {
                  printf("%d ", gcan->labels[i]);
                }
                printf("\n");
                printf(" --> prior (%d,%d,%d)\n", xp, yp, zp);
                gcap = &gca->priors[xp][yp][zp];
                if (gcap == NULL) {
                  continue;
                }
                printf("   prior labels:");
                for (i = 0; i < gcap->nlabels; ++i) {
                  printf("%d ", gcap->labels[i]);
                }
                printf("\n");
                if (logging) {
                  fprintf(logfp, "\n");
                  fclose(logfp);
                }
              }
              ///////////////////////////////////////////
            }
          }
      }
      if (gca->flags & GCA_NO_MRF) {
        continue;
      }
    }
  }

  for (holes_filled = xp = 0; xp < gca->prior_width; xp++) {
    for (yp = 0; yp < gca->prior_height; yp++) {
      for (zp = 0; zp < gca->prior_depth; zp++) {
        if (xp == Gxp && yp == Gyp && zp == Gzp) {
          DiagBreak();
        }
        if (MRIgetVoxVal(mri_mapped, xp, yp, zp, 0) > 0) {
          continue;
        }
        if (!GCApriorToSourceVoxel(gca, mri_inputs, transform, xp, yp, zp, &x, &y, &z)) {
          GC1D *gc;
          label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));
          GCAupdatePrior(gca, mri_inputs, xp, yp, zp, label);
          gcap = &gca->priors[xp][yp][zp];
          GCApriorToNode(gca, xp, yp, zp, &xn, &yn, &zn);
          gc = GCAfindGC(gca, xn, yn, zn, label);
          if (gc == NULL)  // label doesn't exist at this position
          {
            if (gcap->total_training > 1) {
              DiagBreak();
            }
            load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
            if (xp == Gxp && yp == Gyp && zp == Gzp)
              printf(
                  "filling prior hole at (%d, %d, %d) %s,"
                  " node = (%d, %d, %d)\n",
                  xp,
                  yp,
                  zp,
                  cma_label_to_name(label),
                  xn,
                  yn,
                  zn);
            if ((GCAupdateNode(gca, mri_inputs, xn, yn, zn, vals, label, gca_prune, noint) == NO_ERROR) &&
                !(gca->flags & GCA_NO_MRF))
              GCAupdateNodeGibbsPriors(gca, mri_labels, xn, yn, zn, x, y, z, label);
          }
          holes_filled++;
        }
      }
    }
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_mapped, "m.mgz");
  }
  if (holes_filled > 0) {
    printf("%d prior holes filled\n", holes_filled);
  }
  MRIfree(&mri_mapped);
  return (NO_ERROR);
}

int GCAwrite(GCA *gca, const char *fname)
{
  znzFile file;
  int x, y, z, n, i, j;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  GC1D *gc;
  int gzipped = 0;

  if (strstr(fname, ".gcz")) {
    gzipped = 1;
  }

  file = znzopen(fname, "wb", gzipped);
  if (znz_isnull(file)) {
    errno = 0;
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCAwrite(%s): could not open file", fname));
  }

  znzwriteFloat(GCA_INT_VERSION, file);
  znzwriteFloat(gca->prior_spacing, file);
  znzwriteFloat(gca->node_spacing, file);
  znzwriteInt(gca->prior_width, file);
  znzwriteInt(gca->prior_height, file);
  znzwriteInt(gca->prior_depth, file);
  znzwriteInt(gca->node_width, file);
  znzwriteInt(gca->node_height, file);
  znzwriteInt(gca->node_depth, file);
  znzwriteInt(gca->ninputs, file);
  znzwriteInt(gca->flags, file);

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == 139 && y == 103 && z == 139)
        /* wm should be pallidum */
        {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        znzwriteInt(gcan->nlabels, file);
        znzwriteInt(gcan->total_training, file);
        for (n = 0; n < gcan->nlabels; n++) {
          int r, c;
          gc = &gcan->gcs[n];
          znzwriteInt(gcan->labels[n], file);
          for (r = 0; r < gca->ninputs; r++) {
            znzwriteFloat(gc->means[r], file);
          }
          for (r = i = 0; r < gca->ninputs; r++)
            for (c = r; c < gca->ninputs; c++, i++) {
              znzwriteFloat(gc->covars[i], file);
            }

          if (gca->flags & GCA_NO_MRF) {
            continue;
          }
          for (i = 0; i < GIBBS_NEIGHBORS; i++) {
            znzwriteInt(gc->nlabels[i], file);
            for (j = 0; j < gc->nlabels[i]; j++) {
              znzwriteInt((int)gc->labels[i][j], file);
              znzwriteFloat(gc->label_priors[i][j], file);
            }
          }
        }
      }
    }
  }

  for (x = 0; x < gca->prior_width; x++) {
    for (y = 0; y < gca->prior_height; y++) {
      for (z = 0; z < gca->prior_depth; z++) {
        if (x == 139 && y == 103 && z == 139)
        /* wm should be pallidum */
        {
          DiagBreak();
        }
        gcap = &gca->priors[x][y][z];
        if (gcap == NULL) {
          continue;
        }
        znzwriteInt(gcap->nlabels, file);
        znzwriteInt(gcap->total_training, file);
        for (n = 0; n < gcap->nlabels; n++) {
          znzwriteInt((int)gcap->labels[n], file);
          znzwriteFloat(gcap->priors[n], file);
        }
      }
    }
  }

  // if (gca->type == GCA_FLASH || gca->type == GCA_PARAM)
  // always write gca->type
  {
    int n;

    znzwriteInt(FILE_TAG, file); /* beginning of tagged section */

    /* all tags are format: <int: tag> <int: num> <parm> <parm> .... */
    znzwriteInt(TAG_GCA_TYPE, file);
    znzwriteInt(1, file);
    znzwriteInt(gca->type, file);

    if (gca->type == GCA_FLASH) {
      znzwriteInt(TAG_PARAMETERS, file);
      znzwriteInt(3, file); /* currently only storing 3 parameters */
      for (n = 0; n < gca->ninputs; n++) {
        znzwriteFloat(gca->TRs[n], file);
        znzwriteFloat(gca->FAs[n], file);
        znzwriteFloat(gca->TEs[n], file);
      }
    }
  }

  if (gca->ct) {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      printf("writing colortable into GCA file...\n");
    }
    znzwriteInt(TAG_GCA_COLORTABLE, file);
    znzCTABwriteIntoBinary(gca->ct, file);
  }

  // write direction cosine information
  znzwriteInt(TAG_GCA_DIRCOS, file);
  znzwriteFloat(gca->x_r, file);
  znzwriteFloat(gca->x_a, file);
  znzwriteFloat(gca->x_s, file);
  znzwriteFloat(gca->y_r, file);
  znzwriteFloat(gca->y_a, file);
  znzwriteFloat(gca->y_s, file);
  znzwriteFloat(gca->z_r, file);
  znzwriteFloat(gca->z_a, file);
  znzwriteFloat(gca->z_s, file);
  znzwriteFloat(gca->c_r, file);
  znzwriteFloat(gca->c_a, file);
  znzwriteFloat(gca->c_s, file);
  znzwriteInt(gca->width, file);
  znzwriteInt(gca->height, file);
  znzwriteInt(gca->depth, file);
  znzwriteFloat(gca->xsize, file);
  znzwriteFloat(gca->ysize, file);
  znzwriteFloat(gca->zsize, file);

  znzclose(file);

  return (NO_ERROR);
}

GCA *GCAread(const char *fname)
{
  znzFile file;
  int x, y, z, n, i, j;
  GCA *gca;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  GC1D *gc;
  float version, node_spacing, prior_spacing;
  int node_width, node_height, node_depth, ninputs, flags;
  // int prior_width, prior_height, prior_depth;
  int tag;
  int gzipped = 0;
  int tempZNZ;

  if (strstr(fname, ".gcz")) {
    gzipped = 1;
  }

  file = znzopen(fname, "rb", gzipped);
  if (znz_isnull(file)) {
    ErrorReturn(NULL, (ERROR_BADPARM, "GCAread(%s): could not open file", fname));
  }

  if (!znzreadFloatEx(&version, file)) {
    ErrorReturn(NULL, (ERROR_BADPARM, "GCAread(%s): could not read file", fname));
  }

  if (version < GCA_UCHAR_VERSION) {
    node_spacing = znzreadFloat(file);
    node_width = znzreadInt(file);
    node_height = znzreadInt(file);
    node_depth = znzreadInt(file);
    ninputs = znzreadInt(file);
    if (version == 3.0) {
      flags = znzreadInt(file);
    }
    else {
      flags = 0;
    }

    gca = gcaAllocMax(ninputs,
                      node_spacing,
                      node_spacing,
                      node_spacing * node_width,
                      node_spacing * node_height,
                      node_spacing * node_depth,
                      0,
                      flags);
    if (!gca) {
      ErrorReturn(NULL, (Gerror, NULL));
    }

    for (x = 0; x < gca->node_width; x++) {
      for (y = 0; y < gca->node_height; y++) {
        for (z = 0; z < gca->node_depth; z++) {
          if (x == 28 && y == 39 && z == 39) {
            DiagBreak();
          }
          gcan = &gca->nodes[x][y][z];
          gcan->nlabels = znzreadInt(file);
          gcan->total_training = znzreadInt(file);
          if (gcan->nlabels) {
            gcan->labels = (unsigned short *)calloc(gcan->nlabels, sizeof(unsigned short));
            if (!gcan->labels)
              ErrorExit(ERROR_NOMEMORY,
                        "GCAread(%s): could not allocate %d "
                        "labels @ (%d,%d,%d)",
                        fname,
                        gcan->nlabels,
                        x,
                        y,
                        z);
            gcan->gcs = alloc_gcs(gcan->nlabels, flags, gca->ninputs);
            if (!gcan->gcs)
              ErrorExit(ERROR_NOMEMORY,
                        "GCAread(%s); could not allocated %d gcs "
                        "@ (%d,%d,%d)",
                        fname,
                        gcan->nlabels,
                        x,
                        y,
                        z);
          }
          else  // no labels assigned to this node
          {
            gcan->labels = 0;
            gcan->gcs = 0;
          }
          for (n = 0; n < gcan->nlabels; n++) {
            int r, c;
            gc = &gcan->gcs[n];
            znzread1(&tempZNZ, file);
            gcan->labels[n] = (unsigned short)tempZNZ;
            if (gcan->labels[n] > gca->max_label) gca->max_label = gcan->labels[n];
            for (r = 0; r < gca->ninputs; r++) {
              gc->means[r] = znzreadFloat(file);
            }
            for (i = r = 0; r < gca->ninputs; r++)
              for (c = r; c < gca->ninputs; c++, i++) {
                gc->covars[i] = znzreadFloat(file);
              }
            if (gca->flags & GCA_NO_MRF) {
              continue;
            }
            for (i = 0; i < GIBBS_NEIGHBORS; i++) {
              gc->nlabels[i] = znzreadInt(file);

              /* allocate new ones */
              gc->label_priors[i] = (float *)calloc(gc->nlabels[i], sizeof(float));
              if (!gc->label_priors[i])
                ErrorExit(ERROR_NOMEMORY,
                          "GCAread(%s): "
                          "couldn't expand gcs to %d",
                          fname,
                          gc->nlabels);
              gc->labels[i] = (unsigned short *)calloc(gc->nlabels[i], sizeof(unsigned short));
              if (!gc->labels)
                ErrorExit(ERROR_NOMEMORY,
                          "GCAread(%s): couldn't expand "
                          "labels to %d",
                          fname,
                          gc->nlabels[i]);
              for (j = 0; j < gc->nlabels[i]; j++) {
                gc->labels[i][j] = (unsigned short)znzreadInt(file);
                gc->label_priors[i][j] = znzreadFloat(file);
              }
            }
          }
        }
      }
    }
  }
  else /* current version - stores priors at different
                          resolution than densities */
  {
    if (!FEQUAL(version, GCA_UCHAR_VERSION) && !FEQUAL(version, GCA_INT_VERSION)) {
      // fclose(fp) ;
      // myclose(fp);
      znzclose(file);
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "GCAread(%s), version #%2.1f found, "
                   "%2.1f expected",
                   fname,
                   version,
                   GCA_INT_VERSION));
    }
    prior_spacing = znzreadFloat(file);
    node_spacing = znzreadFloat(file);
    // prior_width =
    znzreadInt(file);
    // prior_height =
    znzreadInt(file);
    // prior_depth =
    znzreadInt(file);
    node_width = znzreadInt(file);
    node_height = znzreadInt(file);
    node_depth = znzreadInt(file);
    ninputs = znzreadInt(file);
    flags = znzreadInt(file);

    gca = gcaAllocMax(ninputs,
                      prior_spacing,
                      node_spacing,
                      node_spacing * node_width,
                      node_spacing * node_height,
                      node_spacing * node_depth,
                      0,
                      flags);
    if (!gca) {
      ErrorReturn(NULL, (Gdiag, NULL));
    }

    for (x = 0; x < gca->node_width; x++) {
      for (y = 0; y < gca->node_height; y++) {
        for (z = 0; z < gca->node_depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          gcan = &gca->nodes[x][y][z];
          gcan->nlabels = znzreadInt(file);
          gcan->total_training = znzreadInt(file);
          if (gcan->nlabels) {
            gcan->labels = (unsigned short *)calloc(gcan->nlabels, sizeof(unsigned short));
            if (!gcan->labels)
              ErrorExit(ERROR_NOMEMORY,
                        "GCAread(%s): could not "
                        "allocate %d "
                        "labels @ (%d,%d,%d)",
                        fname,
                        gcan->nlabels,
                        x,
                        y,
                        z);
            gcan->gcs = alloc_gcs(gcan->nlabels, flags, gca->ninputs);
            if (!gcan->gcs)
              ErrorExit(ERROR_NOMEMORY,
                        "GCAread(%s); could not allocated %d gcs "
                        "@ (%d,%d,%d)",
                        fname,
                        gcan->nlabels,
                        x,
                        y,
                        z);
          }
          else  // no labels at this node
          {
            gcan->labels = 0;
            gcan->gcs = 0;
          }
          for (n = 0; n < gcan->nlabels; n++) {
            int r, c;

            gc = &gcan->gcs[n];

            if (version == GCA_UCHAR_VERSION) {
              // gcan->labels[n] = (unsigned short)znzread1(file) ;
              znzread1(&tempZNZ, file);
              gcan->labels[n] = (unsigned short)tempZNZ;
            }
            else {
              gcan->labels[n] = (unsigned short)znzreadInt(file);
            }

            for (r = 0; r < gca->ninputs; r++) {
              gc->means[r] = znzreadFloat(file);
            }
            for (i = r = 0; r < gca->ninputs; r++)
              for (c = r; c < gca->ninputs; c++, i++) {
                gc->covars[i] = znzreadFloat(file);
              }
            if (gca->flags & GCA_NO_MRF) {
              continue;
            }
            for (i = 0; i < GIBBS_NEIGHBORS; i++) {
              gc->nlabels[i] = znzreadInt(file);

              /* allocate new ones */
              gc->label_priors[i] = (float *)calloc(gc->nlabels[i], sizeof(float));
              if (!gc->label_priors[i])
                ErrorExit(ERROR_NOMEMORY,
                          "GCAread(%s): "
                          "couldn't expand gcs to %d",
                          fname,
                          gc->nlabels);
              gc->labels[i] = (unsigned short *)calloc(gc->nlabels[i], sizeof(unsigned short));
              if (!gc->labels)
                ErrorExit(ERROR_NOMEMORY,
                          "GCAread(%s): couldn't expand "
                          "labels to %d",
                          fname,
                          gc->nlabels[i]);
              for (j = 0; j < gc->nlabels[i]; j++) {
                gc->labels[i][j] = (unsigned short)znzreadInt(file);
                gc->label_priors[i][j] = znzreadFloat(file);
              }
            }
          }
        }
      }
    }

    for (x = 0; x < gca->prior_width; x++) {
      for (y = 0; y < gca->prior_height; y++) {
        for (z = 0; z < gca->prior_depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          gcap = &gca->priors[x][y][z];
          if (gcap == NULL) {
            continue;
          }
          gcap->nlabels = znzreadInt(file);
          gcap->total_training = znzreadInt(file);
          if (gcap->nlabels) {
            gcap->labels = (unsigned short *)calloc(gcap->nlabels, sizeof(unsigned short));
            if (!gcap->labels)
              ErrorExit(ERROR_NOMEMORY,
                        "GCAread(%s): could not "
                        "allocate %d "
                        "labels @ (%d,%d,%d)",
                        fname,
                        gcap->nlabels,
                        x,
                        y,
                        z);
            gcap->priors = (float *)calloc(gcap->nlabels, sizeof(float));
            if (!gcap->priors)
              ErrorExit(ERROR_NOMEMORY,
                        "GCAread(%s): could "
                        "not allocate %d "
                        "priors @ (%d,%d,%d)",
                        fname,
                        gcap->nlabels,
                        x,
                        y,
                        z);
          }
          else  // no labels assigned to this priors
          {
            gcap->labels = 0;
            gcap->priors = 0;
          }
          for (n = 0; n < gcap->nlabels; n++) {
            if (version == GCA_UCHAR_VERSION) {
              znzread1(&tempZNZ, file);
              gcap->labels[n] = (unsigned short)tempZNZ;
            }
            else {
              gcap->labels[n] = (unsigned short)znzreadInt(file);
            }
            if (gcap->labels[n] > gca->max_label) gca->max_label = gcap->labels[n];
            gcap->priors[n] = znzreadFloat(file);
          }
        }
      }
    }
  }

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        int xp, yp, zp;

        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        if (gcaNodeToPrior(gca, x, y, z, &xp, &yp, &zp) == NO_ERROR) {
          gcap = &gca->priors[xp][yp][zp];
          if (gcap == NULL) {
            continue;
          }
          for (n = 0; n < gcan->nlabels; n++) {
            gc = &gcan->gcs[n];
            gc->ntraining = gcan->total_training * getPrior(gcap, gcan->labels[n]);
          }
        }
      }
    }
  }

  while (znzreadIntEx(&tag, file)) {
    int n, nparms;

    if (tag == FILE_TAG) /* beginning of tagged section */
    {
      while (znzreadIntEx(&tag, file)) {
        /* all tags are format:
           <int: tag> <int: num> <parm> <parm> .... */
        switch (tag) {
          case TAG_GCA_COLORTABLE:
            /* We have a color table, read it with CTABreadFromBinary. If it
               fails, it will print its own error message. */
            fprintf(stdout, "reading colortable from GCA file...\n");
            gca->ct = znzCTABreadFromBinary(file);
            if (NULL != gca->ct)
              fprintf(stdout, "colortable with %d entries read (originally %s)\n", gca->ct->nentries, gca->ct->fname);
            break;
          case TAG_GCA_TYPE:
            znzreadInt(file); /* skip num=1 */
            gca->type = znzreadInt(file);
            if (DIAG_VERBOSE_ON) switch (gca->type) {
                case GCA_NORMAL:
                  printf("setting gca type = Normal gca type\n");
                  break;
                case GCA_PARAM:
                  printf("setting gca type = T1/PD gca type\n");
                  break;
                case GCA_FLASH:
                  printf("setting gca type = FLASH gca type\n");
                  break;
                default:
                  printf("setting gca type = Unknown\n");
                  gca->type = GCA_UNKNOWN;
                  break;
              }
            break;
          case TAG_PARAMETERS:
            nparms = znzreadInt(file);
            /* how many MR parameters are stored */
            printf("reading %d MR parameters out of GCA header...\n", nparms);
            for (n = 0; n < gca->ninputs; n++) {
              gca->TRs[n] = znzreadFloat(file);
              gca->FAs[n] = znzreadFloat(file);
              gca->TEs[n] = znzreadFloat(file);
              printf(
                  "input %d: TR=%2.1f msec, FA=%2.1f deg, "
                  "TE=%2.1f msec\n",
                  n,
                  gca->TRs[n],
                  DEGREES(gca->FAs[n]),
                  gca->TEs[n]);
            }
            break;
          case TAG_GCA_DIRCOS:
            gca->x_r = znzreadFloat(file);
            gca->x_a = znzreadFloat(file);
            gca->x_s = znzreadFloat(file);
            gca->y_r = znzreadFloat(file);
            gca->y_a = znzreadFloat(file);
            gca->y_s = znzreadFloat(file);
            gca->z_r = znzreadFloat(file);
            gca->z_a = znzreadFloat(file);
            gca->z_s = znzreadFloat(file);
            gca->c_r = znzreadFloat(file);
            gca->c_a = znzreadFloat(file);
            gca->c_s = znzreadFloat(file);
            gca->width = znzreadInt(file);
            gca->height = znzreadInt(file);
            gca->depth = znzreadInt(file);
            gca->xsize = znzreadFloat(file);
            gca->ysize = znzreadFloat(file);
            gca->zsize = znzreadFloat(file);

            if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
              printf("Direction cosines read:\n");
              printf(" x_r = % .4f, y_r = % .4f, z_r = % .4f\n", gca->x_r, gca->y_r, gca->z_r);
              printf(" x_a = % .4f, y_a = % .4f, z_a = % .4f\n", gca->x_a, gca->y_a, gca->z_a);
              printf(" x_s = % .4f, y_s = % .4f, z_s = % .4f\n", gca->x_s, gca->y_s, gca->z_s);
              printf(" c_r = % .4f, c_a = % .4f, c_s = % .4f\n", gca->c_r, gca->c_a, gca->c_s);
            }
            break;
          default:
            ErrorPrintf(ERROR_BADFILE, "GCAread(%s): unknown tag %x\n", fname, tag);
            break;
        }
      }
    }
  }

  GCAsetup(gca);

  znzclose(file);

  return (gca);
}

static int GCAupdatePrior(GCA *gca, MRI *mri, int xn, int yn, int zn, int label)
{
  int n;
  GCA_PRIOR *gcap;

  if (label >= MAX_CMA_LABEL)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "GCAupdatePrior(%d, %d, %d, %d): label out of range", xn, yn, zn, label));

  if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
    DiagBreak();
  }

  gcap = &gca->priors[xn][yn][zn];
  if (gcap == NULL) {
    return -1;
  }
  // find the label histogram index n
  for (n = 0; n < gcap->nlabels; n++) {
    if (gcap->labels[n] == label) {
      break;
    }
  }
  // if index is beyond what we have, then
  if (n >= gcap->nlabels) /* have to allocate a new classifier */
  {
    if (n >= gcap->max_labels) {
      int old_max_labels;
      unsigned short *old_labels;
      float *old_priors;

      old_max_labels = gcap->max_labels;
      gcap->max_labels += 2;
      old_labels = gcap->labels;
      old_priors = gcap->priors;

      /* allocate new ones */
      gcap->priors = (float *)calloc(gcap->max_labels, sizeof(float));

      if (!gcap->priors) ErrorExit(ERROR_NOMEMORY, "GCANupdatePriors: couldn't expand priors to %d", gcap->max_labels);
      gcap->labels = (unsigned short *)calloc(gcap->max_labels, sizeof(unsigned short));
      if (!gcap->labels) ErrorExit(ERROR_NOMEMORY, "GCANupdatePriors: couldn't expand labels to %d", gcap->max_labels);

      /* copy the old ones over */
      memmove(gcap->priors, old_priors, old_max_labels * sizeof(float));
      memmove(gcap->labels, old_labels, old_max_labels * sizeof(unsigned short));

      /* free the old ones */
      free(old_priors);
      free(old_labels);
    }
    // add one
    gcap->nlabels++;
  }

  /* these will be updated when training is complete */
  // give the value at n
  gcap->priors[n] += 1.0f;  // increment counter
  gcap->total_training++;   // increment training counter
  gcap->labels[n] = label;  // save the label value

  return (NO_ERROR);
}

static int GCAupdateNodeCovariance(GCA *gca, MRI *mri, int xn, int yn, int zn, float *vals, int label)
{
  int n, r, c, v;
  GCA_NODE *gcan;
  GC1D *gc;

  if (label >= MAX_CMA_LABEL)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "GCAupdateNodeCovariance(%d, %d, %d, %d): label out of range", xn, yn, zn, label));

  ///////////// debug ////////////////////////////////
  if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z && label == Ggca_label) {
    DiagBreak();
  }
  if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
    DiagBreak();
  }
  ////////////////////////////////////////////////////

  gcan = &gca->nodes[xn][yn][zn];

  for (n = 0; n < gcan->nlabels; n++) {
    if (gcan->labels[n] == label) {
      break;
    }
  }
  if (n >= gcan->nlabels)
    ErrorExit(ERROR_BADPARM, "GCAupdateNodeCovariance(%d, %d, %d, %d): could not find label", xn, yn, zn, label);

  gc = &gcan->gcs[n];

  /* remove means from vals */
  for (r = 0; r < gca->ninputs; r++) {
    vals[r] -= gc->means[r];
  }

  /* these will be normalized when training is complete */
  for (v = r = 0; r < gca->ninputs; r++) {
    for (c = r; c < gca->ninputs; c++, v++) {
      gc->covars[v] += vals[r] * vals[c];
    }
  }

  /* put means back in vals so caller isn't confused */
  for (r = 0; r < gca->ninputs; r++) {
    vals[r] += gc->means[r];
  }

  return (NO_ERROR);
}

static int GCAupdateNode(GCA *gca, MRI *mri, int xn, int yn, int zn, float *vals, int label, GCA *gca_prune, int noint)
{
  int n, i;
  GCA_NODE *gcan;
  GC1D *gc;

  if (label >= MAX_CMA_LABEL)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "GCAupdateNode(%d, %d, %d, %d): label out of range", xn, yn, zn, label));

  if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z && label == Ggca_label) {
    DiagBreak();
  }
  if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
    DiagBreak();
  }

  // if non-zero label and gca_prune is there
  if (label > 0 && gca_prune != NULL && !noint) {
    GCA_NODE *gcan_prune;
    GC1D *gc_prune;
    int nprune;

    gcan_prune = &gca_prune->nodes[xn][yn][zn];

    for (nprune = 0; nprune < gcan_prune->nlabels; nprune++) {
      if (gcan_prune->labels[nprune] == label) {
        break;
      }
    }
    if (nprune >= gcan_prune->nlabels)
      ErrorPrintf(ERROR_BADPARM,
                  "WARNING: pruning GCA at (%d,%d,%d) doesn't "
                  "contain label %d",
                  xn,
                  yn,
                  zn,
                  label);
    gc_prune = &gcan_prune->gcs[nprune];

    if (sqrt(GCAmahDist(gc_prune, vals, gca->ninputs)) > 2)
    /* more than 2 stds from mean */
    {
      total_pruned++;
      return (ERROR_BAD_PARM);
    }
  }

  // get one at this point
  gcan = &gca->nodes[xn][yn][zn];

  // look for this label in the array
  for (n = 0; n < gcan->nlabels; n++) {
    if (gcan->labels[n] == label) {
      break;
    }
  }
  // if not found
  if (n >= gcan->nlabels) /* have to allocate a new classifier */
  {
    if (n >= gcan->max_labels) {
      int old_max_labels;
      unsigned short *old_labels;
      GC1D *old_gcs;

      old_max_labels = gcan->max_labels;
      gcan->max_labels += 2;
      old_labels = gcan->labels;
      old_gcs = gcan->gcs;

/* allocate new ones */
      gcan->gcs = alloc_gcs(gcan->max_labels, gca->flags, gca->ninputs);

      if (!gcan->gcs) ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand gcs to %d", gcan->max_labels);
      gcan->labels = (unsigned short *)calloc(gcan->max_labels, sizeof(unsigned short));
      if (!gcan->labels) ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand labels to %d", gcan->max_labels);

/* copy the old ones over */
      copy_gcs(old_max_labels, old_gcs, gcan->gcs, gca->ninputs);
      memmove(gcan->labels, old_labels, old_max_labels * sizeof(unsigned short));

      /* free the old ones */
      free(old_gcs);
      free(old_labels);
    }
    gcan->nlabels++;
  }

  gc = &gcan->gcs[n];

  /* these will be updated when training is complete */
  if (noint) {
    gc->n_just_priors++;
  }
  else {
    // vals[] array is the values of inputs at this point
    for (i = 0; i < gca->ninputs; i++) {
      gc->means[i] += vals[i];
    }
    // get the mean valu (note it is not divided by ninputs!)
    /*gc->var += val*val ; */
  }
  if (gc->n_just_priors >= gc->ntraining) {
    DiagBreak();
  }
  gcan->total_training++;

  gcan->labels[n] = label;
  gc->ntraining++;

  return (NO_ERROR);
}

int GCAcompleteMeanTraining(GCA *gca)
{
  int x, y, z, n, total_nodes, total_gcs, i, j, holes_filled, total_brain_gcs, total_brain_nodes, r;
  float nsamples;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  GC1D *gc;

  total_nodes = gca->node_width * gca->node_height * gca->node_depth;
  total_brain_nodes = total_gcs = total_brain_gcs = 0;
  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        total_gcs += gcan->nlabels;
        if (gcan->nlabels > 1 || !IS_UNKNOWN(gcan->labels[0])) {
          total_brain_gcs += gcan->nlabels;
          total_brain_nodes++;
        }

        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        for (n = 0; n < gcan->nlabels; n++) {
          gc = &gcan->gcs[n];
          nsamples = gc->ntraining;
          if ((gca->flags & GCA_NO_MRF) == 0) {
            for (i = 0; i < GIBBS_NEIGHBORS; i++) {
              for (j = 0; j < gc->nlabels[i]; j++) {
                gc->label_priors[i][j] /= (float)nsamples;
                check_finite(
                    "GCAcompleteMeanTraining: "
                    "label_priors",
                    gc->label_priors[i][j]);
              }
            }
          }

          nsamples -= gc->n_just_priors;
          /* for no-intensity training */
          if (nsamples > 0) {
            for (r = 0; r < gca->ninputs; r++) {
              gc->means[r] /= nsamples;
              check_finite("GCAcompleteMeanTraining: mean", gc->means[r]);
            }
          }
          else {
            int r;
            for (r = 0; r < gca->ninputs; r++) {
              gc->means[r] = -1; /* mark it for later processing */
            }
          }
        }
      }
    }
  }

  for (x = 0; x < gca->prior_width; x++) {
    for (y = 0; y < gca->prior_height; y++) {
      for (z = 0; z < gca->prior_depth; z++) {
        gcap = &gca->priors[x][y][z];
        if (gcap == NULL) {
          continue;
        }
        for (n = 0; n < gcap->nlabels; n++) {
          gcap->priors[n] /= (float)gcap->total_training;
          check_finite("GCAcompleteMeanTraining: priors", gcap->priors[n]);
        }
      }
    }
  }
  printf("filling holes in the GCA...\n");
  for (holes_filled = x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          gc = &gcan->gcs[n];
          nsamples = gc->ntraining - gc->n_just_priors;
          if (nsamples <= 0) {
            GC1D *gc_nbr;
            int r;
            // int i;

            gc_nbr = GCAfindClosestValidGC(gca, x, y, z, gcan->labels[n], 0);
            if (!gc_nbr) {
              ErrorPrintf(ERROR_BADPARM,
                          "gca(%d,%d,%d,%d) - could not "
                          "find valid nbr label "
                          "%s (%d)",
                          x,
                          y,
                          z,
                          n,
                          cma_label_to_name(gcan->labels[n]),
                          gcan->labels[n]);
              continue;
            }
            holes_filled++;
            // for (i = r = 0; r < gca->ninputs; r++) {
            for (r = 0; r < gca->ninputs; r++) {
              gc->means[r] = gc_nbr->means[r];
            }
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
              printf("filling hole @ (%d, %d, %d)\n", x, y, z);
            }
          }
        }
      }
    }
  }

  printf("%d classifiers: %2.1f per node, %2.2f in brain (%d holes filled)\n",
         total_gcs,
         (float)total_gcs / (float)total_nodes,
         (float)total_brain_gcs / (float)total_brain_nodes,
         holes_filled);
  if (total_pruned > 0) {
    printf("%d samples pruned during training\n", total_pruned);
    total_pruned = 0;
  }
  return (NO_ERROR);
}

int GCAcompleteCovarianceTraining(GCA *gca)
{
  int x, y, z, n, r, c, v, nregs = 0, nsamples;
  GCA_NODE *gcan;
  GC1D *gc;

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];

        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        for (n = 0; n < gcan->nlabels; n++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
            DiagBreak();
          }
          gc = &gcan->gcs[n];
          nsamples = gc->ntraining - gc->n_just_priors;
          /* for no-intensity training */
          if (nsamples > 0) {
            for (r = v = 0; r < gca->ninputs; r++) {
              check_finite("GCAcompleteCovarianceTraining: mean", gc->means[r]);
              if (nsamples > 1) {
                for (c = r; c < gca->ninputs; c++, v++) {
                  gc->covars[v] /= (float)(nsamples - 1);
                  check_finite(
                      "GCAcompleteCovarianceTraining:"
                      " covar",
                      gc->covars[v]);
                  if (r == c)
                  /* diagonal should be positive definite */
                  {
                    if (gc->covars[v] < -0.1) {
                      DiagBreak();
                    }
                  }
                }
              }
            }
          }
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (gcan->labels[n] == Ggca_label || Ggca_label < 0)) {
            MATRIX *m;
            double det;

            det = covariance_determinant(gc, gca->ninputs);
            printf(
                "final covariance matrix for %s "
                "(nsamples=%d, det=%f):\n",
                cma_label_to_name(gcan->labels[n]),
                nsamples,
                det);
            m = load_covariance_matrix(gc, NULL, gca->ninputs);
            MatrixPrint(stdout, m);
            MatrixFree(&m);
            fflush(stdout);
          }
        }
      }
    }
  }

// holes_filled = 0;

  GCAfixSingularCovarianceMatrices(gca);
  /* shouldn't need the code that follows, but haven't made sure yet */

  /* find and fix singular covariance matrices */
  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          MATRIX *m_cov, *m_cov_inv;
          double det;

          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
            DiagBreak();
          }
          gc = &gcan->gcs[n];
          m_cov = load_covariance_matrix(gc, NULL, gca->ninputs);
          m_cov_inv = MatrixInverse(m_cov, NULL);
          det = covariance_determinant(gc, gca->ninputs);
          if (det <= 0) {
            MatrixFree(&m_cov_inv);
            m_cov_inv = NULL;
          }
          if (m_cov_inv == NULL) /* singular */
          {
            MATRIX *m_I;
            m_I = MatrixIdentity(gca->ninputs, NULL);
            MatrixScalarMul(m_I, MIN_VAR, m_I);
            MatrixAdd(m_cov, m_I, m_cov);
            m_cov_inv = MatrixInverse(m_cov, NULL);
            if (m_cov_inv == NULL)
              ErrorExit(ERROR_BADPARM,
                        "GCAcompleteCovarianceTraining: cannot "
                        "regularize singular covariance matrix at "
                        "(%d,%d,%d):%d",
                        x,
                        y,
                        z,
                        n);
            det = covariance_determinant(gc, gca->ninputs);
            if (det < 0) {
              DiagBreak();
            }
            for (v = r = 0; r < gca->ninputs; r++)
              for (c = r; c < gca->ninputs; c++, v++) {
                gc->covars[v] = *MATRIX_RELT(m_cov, r + 1, c + 1);
              }
            MatrixFree(&m_I);
            nregs++;
          }

          MatrixFree(&m_cov_inv);
          MatrixFree(&m_cov);
        }
      }
    }
  }
  gcaCheck(gca);

  if (total_pruned > 0) {
    printf("%d samples pruned during training\n", total_pruned);
    total_pruned = 0;
  }
  return (NO_ERROR);
}

MRI *GCAlabel(MRI *mri_inputs, GCA *gca, MRI *mri_dst, TRANSFORM *transform)
{
  int x, width, height, depth, num_pv, use_partial_volume_stuff;
#if INTERP_PRIOR
  float prior;
#endif

  use_partial_volume_stuff = (getenv("USE_PARTIAL_VOLUME_STUFF") != NULL);
  if (use_partial_volume_stuff) {
    printf("using partial volume calculations in labeling...\n");
  }

  // labeled volume has the same property of the inputs
  if (!mri_dst) {
    mri_dst = MRIalloc(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_INT);
    if (!mri_dst) {
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst");
    }
    MRIcopyHeader(mri_inputs, mri_dst);
  }

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  num_pv = 0;
  // if 0
  // ifdef HAVE_OPENMP
  // pragma omp parallel for if_ROMP(experimental) reduction(+: num_pv)
  // endif
  for (x = 0; x < width; x++) {
    int y, z, n, label, xn, yn, zn;
    // int max_n;
    float vals[MAX_GCA_INPUTS], max_p, p;
    GCA_NODE *gcan;
    GCA_PRIOR *gcap;
    GC1D *gc;
    // GC1D *max_gc;

    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }

        if (x == width / 2 && y == height / 2 && z == depth / 2) {
          DiagBreak();
        }

        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

          gcan = &gca->nodes[xn][yn][zn];
          gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
          if (gcap == NULL) {
            continue;
          }
          label = 0;
          // max_n = -1;
          // max_gc = NULL;
          max_p = 2 * GIBBS_NEIGHBORS * BIG_AND_NEGATIVE;
          // going through gcap labels
          for (n = 0; n < gcap->nlabels; n++) {
            gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[n]);
            if (gc == NULL) {
              gc = GCAfindClosestValidGC(gca, xn, yn, zn, gcap->labels[n], 0);
            }
            if (gc == NULL) {
              MRIsetVoxVal(mri_dst, x, y, z, 0, 0);  // unknown
              continue;
            }
#if INTERP_PRIOR
            prior = gcaComputePrior(gca, mri_inputs, transform, x, y, z, gcap->labels[n]);
            p = gcaComputeLogDensity(gc, vals, gca->ninputs, prior, gcap->labels[n]);
#else
            p = gcaComputeLogDensity(gc, vals, gca->ninputs, gcap->priors[n], gcap->labels[n]);
#endif
            // look for largest p
            if (p > max_p) {
              max_p = p;
              label = gcap->labels[n];
              // max_n = n;
              // max_gc = gc;
            }
          }

          if (use_partial_volume_stuff)
          //////////// start of partial volume stuff
          {
            int n1, l1, l2, max_l1, max_l2, max_n1, max_n2;
            double max_p_pv;

            max_p_pv = -10000;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
              DiagBreak();
            }
            max_l1 = label;
            max_l2 = max_n1 = max_n2 = 0;
            for (n = 0; n < gcap->nlabels; n++)
              for (n1 = n + 1; n1 < gcap->nlabels; n1++) {
                l1 = gcap->labels[n];
                l2 = gcap->labels[n1];
                p = compute_partial_volume_log_posterior(gca, gcan, gcap, vals, l1, l2);
                if (p > max_p_pv) {
                  max_l1 = l1;
                  max_l2 = l2;
                  max_p_pv = p;
                  max_n1 = n;
                  max_n2 = n1;
                }
                if (p > max_p && l1 != label && l2 != label) {
                  DiagBreak();
                }
              }

            /* not the label picked before - change it */
            if (max_p_pv > max_p && max_l1 != label && max_l2 != label) {
              double p1, p2;

              gc = GCAfindGC(gca, xn, yn, zn, max_l1);
              p1 = gcaComputeLogDensity(gc, vals, gca->ninputs, gcap->priors[max_n1], max_l1);
              gc = GCAfindGC(gca, xn, yn, zn, max_l2);
              p2 = gcaComputeLogDensity(gc, vals, gca->ninputs, gcap->priors[max_n2], max_l2);
              num_pv++;
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                printf(
                    "label @ %d, %d, %d: partial volume "
                    "from %s to %s\n",
                    x,
                    y,
                    z,
                    cma_label_to_name(label),
                    cma_label_to_name(p1 > p2 ? max_l1 : max_l2));
              label = p1 > p2 ? max_l1 : max_l2;
              DiagBreak();
            }
          }
          //////////// end of partial volume stuff

          // found the label
          ///////////////////////// debug code /////////////////////
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            int i;
            printf("(%d, %d, %d): inputs=", x, y, z);
            for (i = 0; i < gca->ninputs; i++) {
              printf("%2.1f ", vals[i]);
            }

            printf(
                "\nMAP (no MRF) label %s (%d), log(p)=%2.2e, "
                "node (%d, %d, %d)\n",
                cma_label_to_name(label),
                label,
                max_p,
                xn,
                yn,
                zn);
            dump_gcan(gca, gcan, stdout, 1, gcap);
          }
          /////////////////////////////////////////////
          // set the value
          MRIsetVoxVal(mri_dst, x, y, z, 0, label);
        }
        else {
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0);  // unknown
        }
      }  // z loop
    }    // y loop
  }      // x loop

  return (mri_dst);
}

MRI *GCAlabelProbabilities(MRI *mri_inputs, GCA *gca, MRI *mri_dst, TRANSFORM *transform)
{
  int x, width, height, depth;

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR);
    if (!mri_dst) {
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst");
    }
    MRIcopyHeader(mri_inputs, mri_dst);
  }

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  //#ifdef HAVE_OPENMP
  //#pragma omp parallel for if_ROMP(experimental)
  //#endif
  for (x = 0; x < width; x++) {
    int y, z, xn, yn, zn, n;
    // int label;
    GCA_NODE *gcan;
    GCA_PRIOR *gcap;
    // GC1D *gc;
    double max_p, p, total_p;
    float vals[MAX_GCA_INPUTS];

    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        ///////////////////////////////////////

        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          gcan = &gca->nodes[xn][yn][zn];
          gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
          if (gcap == NULL || gcap->nlabels <= 0) {
            continue;
          }
          // label = 0;
          max_p = 2 * GIBBS_NEIGHBORS * BIG_AND_NEGATIVE;
          // go through labels and find the one with max probability
          for (total_p = 0.0, n = 0; n < gcan->nlabels; n++) {
            // gc = &gcan->gcs[n];

            /* compute 1-d Mahalanobis distance */
            p = GCAcomputePosteriorDensity(gcap, gcan, n, -1, vals, gca->ninputs, xn, yn, zn, gca);
            if (p > max_p) {
              max_p = p;
              // label = gcan->labels[n];
            }
            total_p += p;
          }
          max_p = 255.0 * max_p / total_p;
          if (max_p > 255) {
            max_p = 255;
          }
          MRIsetVoxVal(mri_dst, x, y, z, 0, (BUFTYPE)max_p);
        }
        else {
          MRIsetVoxVal(mri_dst, x, y, z, 0, 255);  // 0;
        }
      }
    }
  }

  return (mri_dst);
}

MRI *GCAcomputeTopNProbabilities(
    MRI *mri_inputs, GCA *gca, MRI *mri_labels, MRI *mri_dst, TRANSFORM *transform, int max_labels)
{
  int x, y, z, width, height, depth, label, xn, yn, zn, n;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  GC1D *gc;
  double p, total_p;
  float vals[MAX_GCA_INPUTS];
  LABEL_PROB label_probs[MAX_LABELS_PER_GCAN];

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  if (!mri_dst) {
    mri_dst = MRIallocSequence(width, height, depth, MRI_FLOAT, 2 * max_labels);
    if (!mri_dst) ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst");
    MRIcopyHeader(mri_inputs, mri_dst);
  }

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == Gx && y == Gy && z == Gz) 
	  DiagBreak();
        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));

          gcan = &gca->nodes[xn][yn][zn];
          gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
          if (gcap == NULL) continue;

          for (total_p = 0.0, n = 0; n < gcap->nlabels; n++) {
            gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[n]);
            if (gc == NULL) {
              label_probs[n].prob = 0;
              label_probs[n].label = -1;
              continue;
            }
            MRIsetVoxVal(mri_labels, x, y, z, 0, gcap->labels[n]);  // change it for markov calculation

            p = gcaGibbsLogPosterior(gca, mri_labels, vals, gcap->labels[n], x, y, z, gcap, gcan, transform);
            p = label_probs[n].prob = exp(p);
            label_probs[n].label = gcap->labels[n];
            label_probs[n].index = n;
            label_probs[n].prob = p;
            total_p += p;
          }
          if (total_p > 0)
            for (n = 0; n < gcap->nlabels; n++) label_probs[n].prob /= -total_p;  // make it negative for sorting

          qsort(label_probs, gcap->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);

          MRIsetVoxVal(mri_labels, x, y, z, 0, label);  // change it back to original label
          for (n = 0; n < max_labels; n++) {
            if (n < gcap->nlabels) {
              MRIsetVoxVal(mri_dst, x, y, z, 2 * n, label_probs[n].label);
              MRIsetVoxVal(mri_dst, x, y, z, 2 * n + 1, -label_probs[n].prob);  // remove the negative from before
            }
            else {
              MRIsetVoxVal(mri_dst, x, y, z, 2 * n, -1);
              MRIsetVoxVal(mri_dst, x, y, z, 2 * n + 1, 0);
            }
          }
        }
      }
    }
  }

  return (mri_dst);
}

MRI *GCAcomputeProbabilities(MRI *mri_inputs, GCA *gca, MRI *mri_labels, MRI *mri_dst, TRANSFORM *transform)
{
  int x, y, z, width, height, depth, label, xn, yn, zn, n;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  // GC1D *gc;
  double label_p, p, total_p;
  float vals[MAX_GCA_INPUTS];

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR);
    if (!mri_dst) {
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst");
    }
    MRIcopyHeader(mri_inputs, mri_dst);
  }

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));

          gcan = &gca->nodes[xn][yn][zn];
          gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
          if (gcap == NULL) {
            continue;
          }
          for (label_p = total_p = 0.0, n = 0; n < gcan->nlabels; n++) {
            // gc = &gcan->gcs[n];

            p = GCAcomputePosteriorDensity(gcap, gcan, n, -1, vals, gca->ninputs, xn, yn, zn, gca);
            if (label == gcan->labels[n]) {
              label_p = p;
            }
            total_p += p;
          }
          label_p = 255.0 * label_p / total_p;
          if (label_p > 255) {
            label_p = 255;
          }
          MRIsetVoxVal(mri_dst, x, y, z, 0, (BUFTYPE)label_p);
        }
      }
    }
  }

  return (mri_dst);
}

#define STARTING_T 500

MRI *GCAannealUnlikelyVoxels(
    MRI *mri_inputs, GCA *gca, MRI *mri_dst, TRANSFORM *transform, int max_iter, MRI *mri_fixed, double prior_factor)
{
  int x, y, z, width, depth, height, *x_indices, *y_indices, *z_indices, nindices, index, iter, nchanged, xn, yn, zn, n,
      nbad, old_label;
  double log_posterior, T, delta_E, p, rn, new_posterior, old_posterior, total_posterior;
  GCA_NODE *gcan;
  MRI *mri_bad;

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  mri_bad = MRIalloc(width, height, depth, MRI_INT);
  MRIcopyHeader(mri_inputs, mri_bad);

  for (nindices = x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (mri_fixed && MRIgetVoxVal(mri_fixed, x, y, z, 0) > 0) {
          continue;
        }

        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          log_posterior = GCAvoxelGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, prior_factor);
          gcan = &gca->nodes[xn][yn][zn];
          if (log_posterior < log(1.0f / (float)gcan->total_training)) {
            MRIsetVoxVal(mri_bad, x, y, z, 0, 1);
            nindices++;
          }
        }
      }
    }
  }

  printf("annealing %d voxels...\n", nindices);
  x_indices = (int *)calloc(nindices, sizeof(int));
  y_indices = (int *)calloc(nindices, sizeof(int));
  z_indices = (int *)calloc(nindices, sizeof(int));
  for (index = x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (MRIgetVoxVal(mri_bad, x, y, z, 0) > 0) {
          x_indices[index] = x;
          y_indices[index] = y;
          z_indices[index] = z;
          index++;
        }
      }
    }
  }

  MRIfree(&mri_bad);
  T = STARTING_T;
  iter = 0;
  do {
    total_posterior = 0.0;
    for (nbad = nchanged = index = 0; index < nindices; index++) {
      x = x_indices[index];
      y = y_indices[index];
      z = z_indices[index];
      if (x == 155 && y == 126 && z == 128) {
        DiagBreak();
      }

      /* find the node associated with this coordinate and classify */
      if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
        gcan = &gca->nodes[xn][yn][zn];

        if (gcan->nlabels == 1) {
          continue;
        }
        n = (int)randomNumber(0.0, (double)gcan->nlabels - 0.0001);
        if (gcan->labels[n] == nint(MRIgetVoxVal(mri_dst, x, y, z, 0))) {
          continue;
        }
        old_posterior = GCAnbhdGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, prior_factor);
        old_label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0));
        MRIsetVoxVal(mri_dst, x, y, z, 0, gcan->labels[n]);
        new_posterior = GCAnbhdGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, prior_factor);
        delta_E = new_posterior - old_posterior;
        p = exp(delta_E / T);
        rn = randomNumber(0.0, 1.0);

        if (p > rn) {
          if (new_posterior < log(1.0f / (float)gcan->total_training)) {
            nbad++;
          }
          nchanged++;
          total_posterior += new_posterior;
        }
        else {
          total_posterior += old_posterior;
          if (old_posterior < log(1.0f / (float)gcan->total_training)) {
            nbad++;
          }
          MRIsetVoxVal(mri_dst, x, y, z, 0, old_label);
        }
      }
    }
    T = T * 0.99;
    fprintf(stdout,
            "%03d: T = %2.2f, nchanged %d, nbad = %d, ll=%2.2f\n",
            iter,
            T,
            nchanged,
            nbad,
            total_posterior / (double)nindices);
    if (!nchanged) {
      break;
    }
  } while (iter++ < max_iter);

  free(x_indices);
  free(y_indices);
  free(z_indices);
  fflush(stdout);
  return (mri_dst);
}

GCA *GCAreduce(GCA *gca_src)
{
  /* have to update to include priors at different resolution than nodes */
  ErrorReturn(NULL, (ERROR_UNSUPPORTED, "GCAreduce: not currently supported"););
}

MRI *GCAclassify(MRI *mri_inputs, GCA *gca, MRI *mri_dst, TRANSFORM *transform, int max_labels)
{
  int x, y, z, width, height, depth, xn, yn, zn, n;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  // GC1D *gc;
  // float max_p
  float p, total_p;
  LABEL_PROB label_probs[1000];
  float vals[MAX_GCA_INPUTS];

  if (max_labels > MAX_LABELS_PER_GCAN || max_labels <= 0) {
    max_labels = MAX_LABELS_PER_GCAN;
  }

  if (!mri_dst) {
    int width = mri_inputs->width, height = mri_inputs->height, depth = mri_inputs->depth;

    mri_dst = MRIallocSequence(width, height, depth, MRI_FLOAT, 2 * max_labels);
    if (!mri_dst) {
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst");
    }
    MRIcopyHeader(mri_inputs, mri_dst);
  }

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == 67 && y == 87 && z == 114) {
          DiagBreak();
        }

        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

        /* find the node associated with this coordinate and classify */
        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          gcan = &gca->nodes[xn][yn][zn];
          gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
          if (gcap == NULL) {
            continue;
          }
          // max_p = 2 * GIBBS_NEIGHBORS * BIG_AND_NEGATIVE;
          for (total_p = 0.0, n = 0; n < gcan->nlabels; n++) {
            // gc = &gcan->gcs[n];

            p = GCAcomputePosteriorDensity(gcap, gcan, n, -1, vals, gca->ninputs, xn, yn, zn, gca);
            total_p += p;
            label_probs[n].prob = p;
            label_probs[n].label = gcan->labels[n];
          }
          /* now sort the labels and probabilities */
          qsort(label_probs, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);

          for (n = 0; n < max_labels; n++) {
            if (n < gcan->nlabels) {
              MRIseq_vox(mri_dst, x, y, z, n * 2) = label_probs[n].label;
              MRIseq_vox(mri_dst, x, y, z, n * 2 + 1) = (BUFTYPE)nint(255.0 * label_probs[n].prob / total_p);
            }
            else {
              MRIseq_vox(mri_dst, x, y, z, n * 2) = 255;
              MRIseq_vox(mri_dst, x, y, z, n * 2 + 1) = 0;
            }
          }
        }
      }
    }
  }

  return (mri_dst);
}

int compare_sort_probabilities(const void *plp1, const void *plp2)
{
  LABEL_PROB *lp1, *lp2;

  lp1 = (LABEL_PROB *)plp1;
  lp2 = (LABEL_PROB *)plp2;

  if (lp1->prob > lp2->prob) {
    return (1);
  }
  else if (lp1->prob < lp2->prob) {
    return (-1);
  }

  return (0);
}

int GCAremoveOutlyingSamples(
    GCA *gca, GCA_SAMPLE *gcas, MRI *mri_inputs, TRANSFORM *transform, int nsamples, float nsigma)
{
  int x, y, z, i, nremoved, xp, yp, zp;
  // int width, height, depth;
  double dist;
  float vals[MAX_GCA_INPUTS];
  /*  GC1D       *gc ;*/

  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  // width = mri_inputs->width;
  // height = mri_inputs->height;
  // depth = mri_inputs->depth;
  TransformInvert(transform, mri_inputs);
  for (nremoved = i = 0; i < nsamples; i++) {
    if (i == Gdiag_no) {
      DiagBreak();
    }
    if (Gdiag_no == gcas[i].label) {
      DiagBreak();
    }
    if (gcas[i].label <= 0) {
      continue;
    }

    xp = gcas[i].xp;
    yp = gcas[i].yp;
    zp = gcas[i].zp;
    if (!GCApriorToSourceVoxel(gca, mri_inputs, transform, xp, yp, zp, &x, &y, &z)) {
      load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

      if (xp == Gxp && yp == Gyp && zp == Gzp) {
        DiagBreak();
      }


      dist = sqrt(GCAsampleMahDist(&gcas[i], vals, gca->ninputs));
      if (dist >= nsigma) {
        nremoved++;
        gcas[i].log_p = BIG_AND_NEGATIVE;
        gcas[i].label = 0;
      }
    }
  }

  if (DIAG_VERBOSE_ON) {
    printf("%d outlying samples removed...\n", nremoved);
  }

  return (nremoved);
}

float GCAnormalizedLogSampleProbability(
    GCA *gca, GCA_SAMPLE *gcas, MRI *mri_inputs, TRANSFORM *transform, int nsamples, double clamp)
{
  int x, y, z, xn, yn, zn, i, n, xp, yp, zp;
  // int width, height, depth;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  GC1D *gc;
  double total_log_p, log_p, norm_log_p;
  float vals[MAX_GCA_INPUTS];

  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  // width = mri_inputs->width;
  // height = mri_inputs->height;
  // depth = mri_inputs->depth;
  TransformInvert(transform, mri_inputs);
  for (total_log_p = 0.0, i = 0; i < nsamples; i++) {
    if (i == Gdiag_no) {
      DiagBreak();
    }
    if (Gdiag_no == gcas[i].label) {
      DiagBreak();
    }

    xp = gcas[i].xp;
    yp = gcas[i].yp;
    zp = gcas[i].zp;
    // do the processing for valid points only
    if (!GCApriorToSourceVoxel(gca, mri_inputs, transform, xp, yp, zp, &x, &y, &z))
      if (!GCApriorToNode(gca, xp, yp, zp, &xn, &yn, &zn)) {
        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
        // use node and prior values
        gcan = &gca->nodes[xn][yn][zn];
        gcap = &gca->priors[xp][yp][zp];
        if (gcap == NULL) {
          continue;
        }
        if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
          DiagBreak();
        }

        for (norm_log_p = 0.0f, n = 0; n < gcan->nlabels; n++) {
          gc = &gcan->gcs[n];
          norm_log_p += GCAcomputeConditionalDensity(gc, vals, gca->ninputs, gcan->labels[n]);
        }
        norm_log_p = log(norm_log_p);
        gc = GCAfindPriorGC(gca, xp, yp, zp, gcas[i].label);
        log_p = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, gcas[i].label);
        if (log_p < -clamp) {
          log_p = -clamp;
        }
        log_p = log(log_p) + gcas_getPriorLog(gcas[i]);
        log_p -= norm_log_p;
        total_log_p += log_p;
        gcas[i].log_p = log_p;

        if (!check_finite("1", total_log_p)) {
          fprintf(stdout, "total log p not finite at (%d, %d, %d)\n", x, y, z);
          DiagBreak();
        }
      }
  }
  fflush(stdout);
  return ((float)total_log_p);
}

float GCAcomputeLabelIntensityVariance(GCA *gca, GCA_SAMPLE *gcas, MRI *mri_inputs, TRANSFORM *transform, int nsamples)
{
  int i, tid, nthreads, label_counts[MAX_CMA_LABELS], label, max_label, nlabels;
  double label_means[MAX_CMA_LABELS], label_vars[MAX_CMA_LABELS], total_var;
  MATRIX *m_prior2source_voxel = NULL;
  VECTOR *v_src[_MAX_FS_THREADS], *v_dst[_MAX_FS_THREADS];

  memset(label_means, 0, sizeof(label_means));
  memset(label_vars, 0, sizeof(label_vars));
  memset(label_counts, 0, sizeof(label_counts));

  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  // store inverse transformation .. forward:input->gca template,
  // inv: gca template->input
  TransformInvert(transform, mri_inputs);

// go through all sample points

/* This loop used to be openmp'ed but there was something in it that
 was not thread-safe and caused unstable/nonrepeatable behavior with
 multimodal inputs. Removing did not seem to slow it down much.
*/
#ifdef HAVE_OPENMP
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif
  for (tid = 0; tid < nthreads; tid++) {
    v_src[tid] = VectorAlloc(4, MATRIX_REAL);
    v_dst[tid] = VectorAlloc(4, MATRIX_REAL);
    *MATRIX_RELT(v_src[tid], 4, 1) = 1.0;
    *MATRIX_RELT(v_dst[tid], 4, 1) = 1.0;
  }

  if (transform->type == MORPH_3D_TYPE)
    m_prior2source_voxel = NULL;
  else
    m_prior2source_voxel = GCAgetPriorToSourceVoxelMatrix(gca, mri_inputs, transform);

  max_label = 0;
  total_var = 0.0;
  // ifdef HAVE_OPENMP
  // pragma omp parallel for if_ROMP(experimental) firstprivate(gcas, tid, m_prior2source_voxel) reduction(max:max_label)
  // endif
  for (i = 0; i < nsamples; i++) {
    int x, y, z, xp, yp, zp;
    float vals[MAX_GCA_INPUTS];

#ifdef HAVE_OPENMP
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    /////////////////// diag code /////////////////////////////
    if (i == Gdiag_no) DiagBreak();
    if (Gdiag_no == gcas[i].label) DiagBreak();
    if (i == Gdiag_no || (gcas[i].xp == Gxp && gcas[i].yp == Gyp && gcas[i].zp == Gzp)) DiagBreak();
    ///////////////////////////////////////////////////////////

    // get prior coordinates in source volume
    if (transform->type == MORPH_3D_TYPE) {
      xp = gcas[i].xp;
      yp = gcas[i].yp;
      zp = gcas[i].zp;
      GCApriorToSourceVoxel(gca, mri_inputs, transform, xp, yp, zp, &x, &y, &z);
    }
    else {
      V3_X(v_src[tid]) = xp = gcas[i].xp;
      V3_Y(v_src[tid]) = yp = gcas[i].yp;
      V3_Z(v_src[tid]) = zp = gcas[i].zp;
      MatrixMultiply(m_prior2source_voxel, v_src[tid], v_dst[tid]);
      x = nint(V3_X(v_dst[tid]));
      y = nint(V3_Y(v_dst[tid]));
      z = nint(V3_Z(v_dst[tid]));
    }
    label = gcas[i].label;
    if (label > max_label) max_label = label;
    label_counts[label]++;
    if (MRIindexNotInVolume(mri_inputs, x, y, z) == 0) {
      // if it is inside the source voxel
      if (x == Gx && y == Gy && z == Gz) DiagBreak();

      // (x,y,z) is the source voxel position
      gcas[i].x = x;
      gcas[i].y = y;
      gcas[i].z = z;

      // get values from all inputs
      load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
      if (label > 0 && vals[0] < 20)  // real tissue shouldn't land where there is no data
      {
        total_var += 1000 * 1000;
      }
      label_means[label] += vals[0];
      label_vars[label] += (vals[0] * vals[0]);

      if (FZERO(vals[0]) && gcas[i].label == Gdiag_no) DiagBreak();
    }
    else {                                    // outside the voxel
      label_means[gcas[i].label] += 1000000;  // BIG
      label_vars[gcas[i].label] += 1000000;   // BIG
    }
  }
  fflush(stdout);

  for (nlabels = label = 0; label <= max_label; label++) {
    if (label_counts[label] > 0) {
      nlabels++;
      label_means[label] /= label_counts[label];
      label_vars[label] = label_vars[label] / label_counts[label] - label_means[label] * label_means[label];
      total_var += label_vars[label];
    }
  }

  for (tid = 0; tid < nthreads; tid++) {
    VectorFree(&v_src[tid]);
    VectorFree(&v_dst[tid]);
  }

  if (m_prior2source_voxel) MatrixFree(&m_prior2source_voxel);
  total_var = sqrt(total_var / nlabels);
  return ((float)-total_var);
}

float GCAcomputeLogSampleProbability(
    GCA * const gca, GCA_SAMPLE * const gcas, MRI * const mri_inputs, TRANSFORM * const transform, int const nsamples, double const clamp)
{

  //  double     outside_log_p = 0.;

  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  // store inverse transformation .. forward:input->gca template,
  // inv: gca template->input
  TransformInvert(transform, mri_inputs);

  // go through all sample points
  double total_log_p = 0.0;

/* This loop used to be openmp'ed but there was something in it that
 was not thread-safe and caused unstable/nonrepeatable behavior with
 multimodal inputs. Removing did not seem to slow it down much.
*/
  int const nthreads =
#ifdef HAVE_OPENMP
    omp_get_max_threads();
#else
    1;
#endif

  MATRIX * m_prior2source_voxel_nonconst = 
    (transform->type == MORPH_3D_TYPE)
    ? NULL
    : GCAgetPriorToSourceVoxelMatrix(gca, mri_inputs, transform);
  MATRIX const * const m_prior2source_voxel = m_prior2source_voxel_nonconst;    // just to make sure not changed below

  VECTOR *v_src[_MAX_FS_THREADS], *v_dst[_MAX_FS_THREADS];
  { int tid;
    for (tid = 0; tid < nthreads; tid++) {
      v_src[tid] = VectorAlloc(4, MATRIX_REAL);
      v_dst[tid] = VectorAlloc(4, MATRIX_REAL);
      *MATRIX_RELT(v_src[tid], 4, 1) = 1.0;
      *MATRIX_RELT(v_dst[tid], 4, 1) = 1.0;
    }
  }
  


#ifdef BEVIN_GCACOMPUTELOGSAMPLEPROBABILITY_REPRODUCIBLE

  #define ROMP_VARIABLE       i
  #define ROMP_LO             0
  #define ROMP_HI             nsamples
    
  #define ROMP_SUMREDUCTION0  total_log_p
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define total_log_p  ROMP_PARTIALSUM(0)
    
#else
  int i;

  ROMP_PF_begin     // important in mri_em_register
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+ : total_log_p)
#endif
  for (i = 0; i < nsamples; i++) {
    ROMP_PFLB_begin

#endif
    
    int x, y, z, xp, yp, zp;
    double log_p;
    float vals[MAX_GCA_INPUTS];

    int const tid =
#ifdef HAVE_OPENMP
      omp_get_thread_num();
#else
      0;
#endif

    /////////////////// diag code /////////////////////////////
    if (i == Gdiag_no) DiagBreak();
    if (Gdiag_no == gcas[i].label) DiagBreak();
    if (i == Gdiag_no || (gcas[i].xp == Gxp && gcas[i].yp == Gyp && gcas[i].zp == Gzp)) DiagBreak();
    ///////////////////////////////////////////////////////////

    // get prior coordinates in source volume
    if (transform->type == MORPH_3D_TYPE) {
      xp = gcas[i].xp;
      yp = gcas[i].yp;
      zp = gcas[i].zp;
      GCApriorToSourceVoxel(gca, mri_inputs, transform, xp, yp, zp, &x, &y, &z);
    }
    else {
      V3_X(v_src[tid]) = xp = gcas[i].xp;
      V3_Y(v_src[tid]) = yp = gcas[i].yp;
      V3_Z(v_src[tid]) = zp = gcas[i].zp;
      MatrixMultiply(m_prior2source_voxel, v_src[tid], v_dst[tid]);
      x = nint(V3_X(v_dst[tid]));
      y = nint(V3_Y(v_dst[tid]));
      z = nint(V3_Z(v_dst[tid]));
    }
    if (MRIindexNotInVolume(mri_inputs, x, y, z) == 0) {
      // if it is inside the source voxel
      if (x == Gx && y == Gy && z == Gz) DiagBreak();

      // (x,y,z) is the source voxel position
      gcas[i].x = x;
      gcas[i].y = y;
      gcas[i].z = z;

      // get values from all inputs

#ifdef FASTER_MRI_EM_REGISTER
      if (gca->ninputs > 1) 
        load_vals_xyzInt(mri_inputs, x, y, z, vals, gca->ninputs);
      else
#endif
        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

#ifdef FASTER_MRI_EM_REGISTER
      if (gca->ninputs == 1)
        log_p = gcaComputeSampleLogDensity_1_input(&gcas[i], vals[0]);
      else 
#endif
        log_p = gcaComputeSampleLogDensity(&gcas[i], vals, gca->ninputs);

      if (FZERO(vals[0]) && gcas[i].label == Gdiag_no) {
        if (fabs(log_p) < 5) DiagBreak();
        DiagBreak();
      }
      if (log_p < -clamp) log_p = -clamp;
    }
    else {               // outside the voxel
      log_p = -1000000;  // BIG_AND_NEGATIVE;
      // log(VERY_UNLIKELY); // BIG_AND_NEGATIVE;
      //      outside_log_p += log_p;
    }
    gcas[i].log_p = log_p;
    total_log_p += log_p;

#ifdef BEVIN_GCACOMPUTELOGSAMPLEPROBABILITY_REPRODUCIBLE

    #undef total_log_p
  #include "romp_for_end.h"
  
#else

    ROMP_PFLB_end
  }
  ROMP_PF_end

#endif

  fflush(stdout);

  { int tid;
    for (tid = 0; tid < nthreads; tid++) {
      VectorFree(&v_src[tid]);
      VectorFree(&v_dst[tid]);
  } }
  
  if (m_prior2source_voxel_nonconst) MatrixFree(&m_prior2source_voxel_nonconst);
  return ((float)total_log_p / nsamples);
}

float GCAcomputeLogSampleProbabilityLongitudinal(
    GCA *gca, GCA_SAMPLE *gcas, MRI *mri_inputs, TRANSFORM *transform, int nsamples, double clamp)
{
  int x, y, z, i, xp, yp, zp;
  // int width, height, depth;
  float vals[MAX_GCA_INPUTS];
  double total_log_p, log_p;
  int countOutside = 0, frame;
  double outside_log_p = 0.;
  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  // width = mri_inputs->width;
  // height = mri_inputs->height;
  // depth = mri_inputs->depth;
  // store inverse transformation .. forward:input->gca template,
  // inv: gca template->input
  TransformInvert(transform, mri_inputs);

  // go through all sample points
  for (total_log_p = 0.0, i = 0; i < nsamples; i++) {
    /////////////////// diag code /////////////////////////////
    if (i == Gdiag_no) {
      DiagBreak();
    }
    if (Gdiag_no == gcas[i].label) {
      DiagBreak();
    }
    if (i == Gdiag_no || (gcas[i].xp == Gxp && gcas[i].yp == Gyp && gcas[i].zp == Gzp)) {
      DiagBreak();
    }
    ///////////////////////////////////////////////////////////

    // get prior coordinates
    xp = gcas[i].xp;
    yp = gcas[i].yp;
    zp = gcas[i].zp;
    // if it is inside the source voxel
    if (!GCApriorToSourceVoxel(gca, mri_inputs, transform, xp, yp, zp, &x, &y, &z)) {
      if (x == Gx && y == Gy && z == Gz) {
        DiagBreak();
      }

      // (x,y,z) is the source voxel position
      gcas[i].x = x;
      gcas[i].y = y;
      gcas[i].z = z;
      // get values from all inputs
      load_vals(mri_inputs, x, y, z, vals, mri_inputs->nframes);
      for (frame = 0; frame < mri_inputs->nframes; frame++) {
        log_p = gcaComputeSampleLogDensity(&gcas[i], &vals[frame], 1);
        if (FZERO(vals[0]) && gcas[i].label == Gdiag_no) {
          if (fabs(log_p) < 5) {
            DiagBreak();
          }
          DiagBreak();
        }
        if (log_p < -clamp) {
          log_p = -clamp;
        }
        total_log_p += log_p;
        gcas[i].log_p = log_p;

        if (!check_finite("2", total_log_p)) {
          fprintf(stdout, "total log p not finite at (%d, %d, %d)\n", x, y, z);
          DiagBreak();
        }
      }
    }
    else  // outside the voxel
    {
      log_p = -1000000;  // BIG_AND_NEGATIVE;
      // log(VERY_UNLIKELY); // BIG_AND_NEGATIVE;
      total_log_p += log_p;
      gcas[i].log_p = log_p;
      outside_log_p += log_p;
      countOutside++;
    }
  }

#ifndef __OPTIMIZE__
#endif
  fflush(stdout);

  return ((float)total_log_p);
}

float GCAcomputeLogSampleProbabilityUsingCoords(
    GCA *gca, GCA_SAMPLE *gcas, MRI *mri_inputs, TRANSFORM *transform, int nsamples, double clamp)
{
  int x, y, z, xp, yp, zp, xn, yn, zn, i;
  // int width, height, depth;
  float vals[MAX_GCA_INPUTS];
  double total_log_p, log_p;

  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  // width = mri_inputs->width;
  // height = mri_inputs->height;
  // depth = mri_inputs->depth;
  for (total_log_p = 0.0, i = 0; i < nsamples; i++) {
    if (i == Gdiag_no) {
      DiagBreak();
    }
    if (Gdiag_no == gcas[i].label) {
      DiagBreak();
    }

    xp = gcas[i].xp;
    yp = gcas[i].yp;
    zp = gcas[i].zp;
    // do the processing only for valid points
    if (!GCApriorToNode(gca, xp, yp, zp, &xn, &yn, &zn)) {
      x = gcas[i].x;
      y = gcas[i].y;
      z = gcas[i].z;
      load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

      log_p = gcaComputeSampleLogDensity(&gcas[i], vals, gca->ninputs);
      if (log_p < -clamp) {
        log_p = -clamp;
      }
      total_log_p += log_p;
      gcas[i].log_p = log_p;

      if (!check_finite("3", total_log_p)) {
        DiagBreak();
        fprintf(stdout, "total log p not finite at (%d, %d, %d)\n", x, y, z);
      }
    }
    else {
      log_p = log(VERY_UNLIKELY);  // BIG_AND_NEGATIVE;
      total_log_p += log_p;
      gcas[i].log_p = log_p;
    }
  }
  fflush(stdout);

  return ((float)total_log_p);
}

/*
  compute the probability of the image given the transform and the class
  stats.
*/
float GCAcomputeLogImageProbability(GCA *gca, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform)
{
  int x, y, z, width, height, depth, xn, yn, zn, n, label;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  GC1D *gc;
  double total_log_p;
  float vals[MAX_GCA_INPUTS];

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  for (total_log_p = 0.0, x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == 85 && y == 89 && z == 135) {
          DiagBreak();
        }

        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
        label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));

        /* find the node associated with this coordinate and classify */
        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          gcan = &gca->nodes[xn][yn][zn];
          gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
          if (gcap == NULL) {
            continue;
          }
          for (n = 0; n < gcan->nlabels; n++) {
            if (gcan->labels[n] == label) {
              break;
            }
          }
          if (n < gcan->nlabels) {
            gc = &gcan->gcs[n];

            total_log_p += gcaComputeLogDensity(gc, vals, gca->ninputs, getPrior(gcap, label), label);
            if (!check_finite("4", total_log_p)) {
              DiagBreak();
              fprintf(stdout,
                      "total log p not finite at (%d, %d, %d)"
                      " n = %d,\n",
                      x,
                      y,
                      z,
                      n);
            }
          }
        }
      }
    }
  }
  fflush(stdout);

  return ((float)total_log_p);
}


#define MIN_SPACING 8

int GCAtransformSamples(GCA *gca_src, GCA *gca_dst, GCA_SAMPLE *gcas, int nsamples)
{
  int scale, i, label, xd, yd, zd, xs, ys, zs, n, xk, yk, zk, xd0, yd0, zd0, min_y_i, xdn, ydn, zdn;
  float max_p, vscale, min_y, prior;
  // float min_v;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  // GC1D *gc;
  MRI *mri_found;

  vscale = 1;
  mri_found =
      MRIalloc(gca_dst->node_width * vscale, gca_dst->node_height * vscale, gca_dst->node_depth * vscale, MRI_UCHAR);
  // change the voxel size
  mri_found->xsize = gca_dst->node_spacing * vscale;
  mri_found->ysize = gca_dst->node_spacing * vscale;
  mri_found->zsize = gca_dst->node_spacing * vscale;

  // copy direction cosines
  GCAcopyDCToMRI(gca_dst, mri_found);

  scale = gca_src->prior_spacing / gca_dst->prior_spacing;
  min_y = 10000;
  min_y_i = -1;
  for (i = 0; i < nsamples; i++) {
    label = gcas[i].label;

    xs = gcas[i].xp;
    ys = gcas[i].yp;
    zs = gcas[i].zp;
    xd0 = (int)(xs * scale);
    yd0 = (int)(ys * scale);
    zd0 = (int)(zs * scale);
    max_p = -1.0; /* find appropriate label with highest prior in dst */
    // min_v = 10000000.0f;
    for (xk = -scale / 2; xk <= scale / 2; xk++) {
      xd = MIN(MAX(0, xd0 + xk), gca_dst->prior_width - 1);
      for (yk = -scale / 2; yk <= scale / 2; yk++) {
        yd = MIN(MAX(0, yd0 + yk), gca_dst->prior_height - 1);
        for (zk = -scale / 2; zk <= scale / 2; zk++) {
          zd = MIN(MAX(0, zd0 + zk), gca_dst->prior_height - 1);
          if (MRIgetVoxVal(mri_found, (int)(xd * vscale), (int)(yd * vscale), (int)(zd * vscale), 0)) {
            continue;
          }
          if (!GCApriorToNode(gca_dst, xd, yd, zd, &xdn, &ydn, &zdn)) {
            gcan = &gca_dst->nodes[xdn][ydn][zdn];
            gcap = &gca_dst->priors[xd][yd][zd];
            if (gcap == NULL || gcap->nlabels <= 0) {
              continue;
            }
            for (n = 0; n < gcan->nlabels; n++) {
              // gc = &gcan->gcs[n];
              prior = getPrior(gcap, gcan->labels[n]);
              if (gcan->labels[n] == label && (prior > max_p
                                               /*|| FEQUAL(prior,max_p) && gc->var < min_v)*/)) {
                /*              min_v = gc->var ;*/
                max_p = prior;
                gcas[i].xp = xd;
                gcas[i].yp = yd;
                gcas[i].zp = zd;
              }
            }
          }
        }
      }
    }
    if (max_p < 0) {
      fprintf(stdout, "WARNING: label %d not found at (%d,%d,%d)\n", label, gcas[i].xp, gcas[i].yp, gcas[i].zp);
      DiagBreak();
    }
    MRIsetVoxVal(mri_found, (int)(gcas[i].xp * vscale), (int)(gcas[i].yp * vscale), (int)(gcas[i].zp * vscale), 0, 1);
    if (gcas[i].yp < min_y) {
      min_y = gcas[i].yp;
      min_y_i = i;
    }
  }

  i = min_y_i;
  fprintf(stdout, "min_y = (%d, %d, %d) at i=%d, label=%d\n", gcas[i].xp, gcas[i].yp, gcas[i].zp, i, gcas[i].label);
  MRIfree(&mri_found);
  fflush(stdout);
  return (NO_ERROR);
}

/* don't use a label with prior less than this */
#define MIN_MAX_PRIORS 0.5

static int exclude_classes[] = {
    0 /*1, 6, 21, 22, 23, 24, 25, 30, 57, 61, 62, 63*/
};

#define TILES 3
#define MAX_PCT .1

static int compare_gca_samples(const void *pc1, const void *pc2);

static int compare_gca_samples(const void *pgcas1, const void *pgcas2)
{
  GCA_SAMPLE *gcas1, *gcas2;

  gcas1 = (GCA_SAMPLE *)pgcas1;
  gcas2 = (GCA_SAMPLE *)pgcas2;

  /*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (FEQUAL(gcas_getPrior(*gcas1), gcas_getPrior(*gcas2))) {

    /* check size of determinants */
    return (1);
  }

  if (gcas_getPrior(*gcas1) > gcas_getPrior(*gcas2)) {
    return (-1);
  }
  else if (gcas_getPrior(*gcas1) < gcas_getPrior(*gcas2)) {
    return (1);
  }

  return (0);
}
#define MAX_SPACING 16 /* mm */

GCA_SAMPLE *GCAfindStableSamplesByLabel(GCA *gca, int nsamples, float min_prior)
{
  GCA_SAMPLE *gcas, *ordered_labels[MAX_DIFFERENT_LABELS], *gcas2;
  GCA_PRIOR *gcap;
  GC1D *gc;
  int found[MAX_DIFFERENT_LABELS], x, y, z, width, height, depth, n, label, nfound, samples_added[MAX_DIFFERENT_LABELS];
  float histo[MAX_DIFFERENT_LABELS], spacing, scale;
  int volume[TILES][TILES], max_class, total, extra, total_found, label_counts[MAX_DIFFERENT_LABELS],
      current_index[MAX_DIFFERENT_LABELS], ordered_label_counts[MAX_DIFFERENT_LABELS], i, index, *x_indices, *y_indices,
      *z_indices, nindices;

  memset(histo, 0, sizeof(histo));
  memset(label_counts, 0, sizeof(label_counts));
  memset(samples_added, 0, sizeof(samples_added));
  memset(current_index, 0, sizeof(current_index));
  memset(ordered_label_counts, 0, sizeof(ordered_label_counts));
  memset(found, 0, sizeof(found));
  memset(volume, 0, sizeof(volume));
  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE));
  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;

  /* compute the max priors and min variances for each class */
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        gcap = &gca->priors[x][y][z];
        if (gcap == NULL) {
          continue;
        }
        for (n = 0; n < gcap->nlabels; n++) {
          label = gcap->labels[n];
          if (label == Gdiag_no) {
            DiagBreak();
          }
          histo[label] += (float)nint(gcap->total_training * gcap->priors[n]);
          label_counts[label]++;
        }
      }
    }
  }

  for (max_class = MAX_DIFFERENT_LABELS; max_class >= 0; max_class--)
    if (histo[max_class] > 0) {
      break;
    }
  fprintf(stdout, "max class = %d\n", max_class);

  /* count total # of samples */
  for (total = 0, n = 1; n <= max_class; n++) {
    total += histo[n];
  }

  fprintf(stdout, "%d total training samples found\n", total);
  fflush(stdout);

  /* turn histogram into samples/class */
  for (n = 1; n <= max_class; n++) {
    histo[n] = (float)nint(0.25 + histo[n] * (float)nsamples / total);
  }

  for (n = 0; (unsigned)n < sizeof(exclude_classes) / sizeof(exclude_classes[0]); n++) {
    histo[exclude_classes[n]] = 0;
  }

  /* crop max # per class */
  for (n = 1; n <= max_class; n++)
    if (histo[n] > nint(MAX_PCT * nsamples)) {
      histo[n] = nint(MAX_PCT * nsamples);
    }

  for (extra = 0, n = 1; n <= max_class; n++) {
    extra += histo[n];
  }
  extra = nsamples - extra; /* from rounding */
  printf("%d extra samples for redistribution...\n", extra);

  /* first add to classes with only one sample */
  for (n = 1; extra > 0 && n <= max_class; n++) {
    if (histo[n] == 1) {
      histo[n]++;
      extra--;
    }
  }

  while (extra > 0) /* add 1 to each class */
  {
    for (n = 1; extra > 0 && n <= max_class; n++) {
      if (histo[n] >= 1) {
        histo[n]++;
        extra--;
      }
    }
  }

  {
    FILE *fp;
    fp = fopen("classes.dat", "w");
    for (n = 1; n <= max_class; n++)
      if (!FZERO(histo[n])) {
        fprintf(fp, "%d  %2.1f\n", n, histo[n]);
      }
    fclose(fp);
  }

  /* allocate arrays that will be used for sorting */
  for (n = 1; n <= max_class; n++)
    if (histo[n] > 0) {
      ordered_labels[n] = (GCA_SAMPLE *)calloc(label_counts[n], sizeof(GCA_SAMPLE));
      if (!ordered_labels[n])
        ErrorExit(ERROR_NO_MEMORY,
                  "GCAfindStableSamplesByLabel: "
                  "could not allocate %d samples "
                  "for label %d",
                  label_counts[n],
                  n);
    }

  nindices = width * height * depth;
  x_indices = (int *)calloc(nindices, sizeof(int));
  y_indices = (int *)calloc(nindices, sizeof(int));
  z_indices = (int *)calloc(nindices, sizeof(int));

  for (i = 0; i < nindices; i++) {
    x_indices[i] = i % width;
    y_indices[i] = (i / width) % height;
    z_indices[i] = (i / (width * height)) % depth;
  }
  for (i = 0; i < nindices; i++) {
    int tmp;
    index = (int)randomNumber(0.0, (double)(nindices - 0.0001));

    tmp = x_indices[index];
    x_indices[index] = x_indices[i];
    x_indices[i] = tmp;

    tmp = y_indices[index];
    y_indices[index] = y_indices[i];
    y_indices[i] = tmp;

    tmp = z_indices[index];
    z_indices[index] = z_indices[i];
    z_indices[i] = tmp;
  }

  for (index = 0; index < nindices; index++) {
    x = x_indices[index];
    y = y_indices[index];
    z = z_indices[index];
    gcap = &gca->priors[x][y][z];
    if (gcap == NULL) {
      continue;
    }
    for (n = 0; n < gcap->nlabels; n++) {
      label = gcap->labels[n];
      gc = GCAfindPriorGC(gca, x, y, z, label);
      if (histo[label] > 0) {
        i = ordered_label_counts[label];
        ordered_label_counts[label]++;
        ordered_labels[label][i].means = (float *)calloc(gca->ninputs, sizeof(float));
        ordered_labels[label][i].covars = (float *)calloc((gca->ninputs * (gca->ninputs + 1)) / 2, sizeof(float));
        if (!ordered_labels[label][i].means || !ordered_labels[label][i].covars)
          ErrorExit(ERROR_NOMEMORY,
                    "could not allocate mean (%d) and "
                    "covariance (%d) matrices",
                    gca->ninputs,
                    gca->ninputs * (gca->ninputs + 1) / 2);
        ordered_labels[label][i].xp = x;
        ordered_labels[label][i].yp = y;
        ordered_labels[label][i].zp = z;
        ordered_labels[label][i].label = label;
        gcas_setPrior(ordered_labels[label][i], getPrior(gcap, label));
        if (!gc) {
          int r, c, v;
          for (v = r = 0; r < gca->ninputs; r++) {
            ordered_labels[label][i].means[r] = 0.0;
            for (c = r; c < gca->ninputs; c++, v++) {
              if (c == r) {
                ordered_labels[label][i].covars[v] = 1.0;
              }
              else {
                ordered_labels[label][i].covars[v] = 0;
              }
            }
          }
        }
        else {
          int r, c, v;
          for (v = r = 0; r < gca->ninputs; r++) {
            ordered_labels[label][i].means[r] = gc->means[r];
            for (c = r; c < gca->ninputs; c++, v++) {
              ordered_labels[label][i].covars[v] = gc->covars[v];
            }
          }
        }
      }
    }
  }

  free(x_indices);
  free(y_indices);
  free(z_indices);

  /* sort each list of labels by prior, then by variance */
  for (n = 1; n <= max_class; n++)
    if (histo[n] > 0) {
      qsort(ordered_labels[n], ordered_label_counts[n], sizeof(GCA_SAMPLE), compare_gca_samples);
    }

  total_found = nfound = 0;

  for (spacing = MAX_SPACING; spacing >= gca->node_spacing; spacing /= 2) {
    MRI *mri_found;
    int xv, yv, zv;

    for (n = 1; n <= max_class; n++) {
      current_index[n] = 0;
    }

    scale = gca->prior_spacing / spacing;
    mri_found = MRIalloc(width * scale, height * scale, depth * scale, MRI_UCHAR);

    // change the size
    mri_found->xsize = gca->prior_spacing * scale;
    mri_found->ysize = gca->prior_spacing * scale;
    mri_found->zsize = gca->prior_spacing * scale;

    GCAcopyDCToMRI(gca, mri_found);

    for (i = 0; i < total_found; i++) {
      xv = gcas[i].xp * scale;
      yv = gcas[i].yp * scale;
      zv = gcas[i].zp * scale;
      MRIsetVoxVal(mri_found, xv, yv, zv, 0, 1);
    }

    do {
      nfound = 0;
      for (n = 1; n <= max_class; n++) {
        if (samples_added[n] < histo[n])
        /* add another sample from this class*/
        {
          while ((ordered_labels[n][current_index[n]].label < 0) && current_index[n] < ordered_label_counts[n]) {
            current_index[n]++;
          }
          if (current_index[n] < ordered_label_counts[n]) {
            gcas2 = &ordered_labels[n][current_index[n]];
            current_index[n]++;
            xv = gcas2->xp * scale;
            yv = gcas2->yp * scale;
            zv = gcas2->zp * scale;

            if (n == Gdiag_no) {
              DiagBreak();
            }

            if (!(int)MRIgetVoxVal(mri_found, xv, yv, zv, 0))
            /* none at this location yet */
            {
              if (gcas2->label == Gdiag_no) {
                DiagBreak();
              }
              samples_added[n]++;
              gcas[total_found + nfound++] = *gcas2;
              MRIsetVoxVal(mri_found, xv, yv, zv, 0, gcas2->label);
              gcas2->label = -1;
              if (nfound + total_found >= nsamples) {
                break;
              }
            }
          }
        }
      }

      total_found += nfound;
    } while ((nfound > 0) && total_found < nsamples);
    MRIfree(&mri_found);
  }

  if (total_found != nsamples) {
    ErrorPrintf(ERROR_BADPARM, "could only find %d samples!\n", total_found);
  }
  return (gcas);
}

GCA_SAMPLE *GCAfindContrastSamples(GCA *gca, int *pnsamples, int min_spacing, float min_prior)
{
  ErrorReturn(NULL, (ERROR_UNSUPPORTED, "GCAfindContrastSamples: no longer supported"));
}

GCA_SAMPLE *GCAfindStableSamples(GCA *gca,
                                 int *pnsamples,
                                 int min_spacing,
                                 float min_prior,
                                 int *exclude_list,
                                 int unknown_nbr_spacing,
                                 int vent_spacing)
{
  GCA_SAMPLE *gcas;
  int xi, yi, zi, width, height, depth, label, nfound;
  int label_counts[MAX_DIFFERENT_LABELS], best_label, nzeros, r;
  float max_prior, total_means[MAX_GCA_INPUTS];
  // float best_mean_dist;
  float priors[MAX_DIFFERENT_LABELS];
  float means[MAX_DIFFERENT_LABELS][MAX_GCA_INPUTS];
  float vars[MAX_DIFFERENT_LABELS], max_priors[MAX_DIFFERENT_LABELS];
  float prior_stride, x, y, z, min_unknown[MAX_GCA_INPUTS];
  float max_unknown[MAX_GCA_INPUTS], prior_factor;
  MRI *mri_filled, *mri_vent_dist = NULL;

  if (gca->max_label >= MAX_DIFFERENT_LABELS)
    ErrorExit(
        ERROR_UNSUPPORTED, "GCAfindStableSamples: max label %d too large (%d)\n", gca->max_label, MAX_DIFFERENT_LABELS);

#define MIN_UNKNOWN_DIST 2

  if (vent_spacing > 0) {
    MRI *mri_labels = GCAbuildMostLikelyLabelVolume(gca, NULL);
    MRI *mri_vent;
    mri_vent = MRIclone(mri_labels, NULL);
    MRIcopyLabel(mri_labels, mri_vent, Left_Lateral_Ventricle);
    MRIcopyLabel(mri_labels, mri_vent, Right_Lateral_Ventricle);
    //    MRIcopyLabel(mri_labels, mri_vent, Left_Inf_Lat_Vent) ;
    //    MRIcopyLabel(mri_labels, mri_vent, Right_Inf_Lat_Vent) ;
    MRIbinarize(mri_vent, mri_vent, 1, 0, 1);

    //    MRIwrite(mri_vent, "v.mgz") ;
    //    MRIwrite(mri_labels, "l.mgz") ;
    mri_vent_dist = MRIdistanceTransform(mri_vent, NULL, 1, mri_vent->width, DTRANS_MODE_SIGNED, NULL);
    //    MRIwrite(mri_vent_dist, "d.mgz") ;
    MRIfree(&mri_vent);
    MRIfree(&mri_labels);
  }
  gcaFindMaxPriors(gca, max_priors);
  gcaFindIntensityBounds(gca, min_unknown, max_unknown);
  printf("bounding unknown intensity as < ");
  for (r = 0; r < gca->ninputs; r++) {
    printf("%2.1f ", min_unknown[r]);
  }
  printf("or > ");
  for (r = 0; r < gca->ninputs; r++) {
    printf("%2.1f ", max_unknown[r]);
  }
  printf("\n");

  memset(label_counts, 0, sizeof(label_counts));

  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;

  // samples size allocated is prior voxel points
  gcas = (GCA_SAMPLE *)calloc(width * height * depth, sizeof(GCA_SAMPLE));
  if (!gcas)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate %dK x %d samples\n",
              "GCAfindStableSamples",
              width * height * depth / 1024,
              sizeof(GCA_SAMPLE));

  // create a full size volume prior_width*prior_spacing = image width
  mri_filled = MRIalloc(width * gca->prior_spacing, height * gca->prior_spacing, depth * gca->prior_spacing, MRI_UCHAR);
  // use the mri_prior_ header

  mri_filled->xsize = gca->xsize;
  mri_filled->ysize = gca->ysize;
  mri_filled->zsize = gca->zsize;

  GCAcopyDCToMRI(gca, mri_filled);

  prior_stride = (float)min_spacing / (float)gca->prior_spacing;
  if (prior_stride >= 2.5) {
    prior_factor = .5;
  }
  else if (prior_stride >= 1.5) {
    prior_factor = .75;
  }
  else {
    prior_factor = .9;
  }

  for (r = 0; r < gca->ninputs; r++) {
    total_means[r] = 0.0;
  }

  /* compute the max priors and min variances for each class */
  for (nzeros = nfound = x = 0; x < width; x += prior_stride) {
    fflush(stdout);  // nicknote: this prevents a segfault on Linux PowerPC
    // when -O2 optimization is used w/ gcc 3.3.3

    xi = nint(x);
    for (y = 0; y < height; y += prior_stride) {
      fflush(stdout);  // nicknote: this prevents a segfault on Linux PowerPC
      // when -O2 optimization is used w/ gcc 3.3.3

      yi = nint(y);
      for (z = 0; z < depth; z += prior_stride) {
        zi = nint(z);
        if (xi == Gx && yi == Gy && zi == Gz) {
          DiagBreak();
        }
        if (abs(x - 31) <= prior_stride && abs(y - 22) <= prior_stride && abs(z - 36) <= prior_stride) {
          DiagBreak();
        }
        gcaRegionStats(gca, x, y, z, prior_stride / 2, priors, vars, means);

        ///////////// diag code /////////////////////////////////////
        if (priors[4] > .5 * min_prior ||  /* left lat ven */
            priors[5] > .5 * min_prior ||  /* inf left lat ven */
            priors[14] > .5 * min_prior || /* 3rd ven */
            priors[15] > .5 * min_prior || /* 4th ven */
            priors[43] > .5 * min_prior || priors[44] > .5 * min_prior) {
          DiagBreak();
        }
        /////////////////////////////////////////////////////////////

        // best_mean_dist = 0.0;
        max_prior = -1;

        if ((different_nbr_max_labels(gca, x, y, z, unknown_nbr_spacing, 0) > 0) &&
            (exclude_list && exclude_list[0] == 0) && (priors[0] >= min_prior) && (priors[0] >= .5 * max_priors[0])) {
          best_label = 0;
        }
        else {
          best_label = -1;
        }

        for (label = 1; label < MAX_DIFFERENT_LABELS; label++) {
          /* ignore it if:
             1. It's prior is less than min_prior and
             it's prior is less than the max for this class.
             2. It's prior is less than 1/2 of min prior regardless.
          */
          if (((priors[label] < min_prior) && (priors[label] < prior_factor * max_priors[label])) ||
              (priors[label] < .5 * min_prior) || (exclude_list && exclude_list[label] > 0)) {
            continue;
          }

          if ((IS_WM(label) || IS_THALAMUS(label) || IS_CAUDATE(label)) && mri_vent_dist &&
              (MRIgetVoxVal(mri_vent_dist, x, y, z, 0) < vent_spacing))
            continue;

          if ((best_label == 0) || (priors[label] > max_prior) ||
              (FEQUAL(priors[label], max_prior) && label_counts[best_label] > label_counts[label])) {
            best_label = label;
            max_prior = priors[label];
          }
        }
        if (nfound > 0) {
          double p = randomNumber(0, 1.0);
          if (p < ((double)label_counts[best_label] / nfound)) {
            continue;
          }
        }
        if (best_label >= 0) {
          if (best_label == 0) {
            if (gcaFindBestSample(gca, x, y, z, best_label, prior_stride / 2, &gcas[nfound]) == NO_ERROR)
            {
              int xv, yv, zv;

              if (((means[0] <= min_unknown) || (means[0] >= max_unknown))) /* disabled check */
              {
                if (!GCApriorToVoxel(gca, mri_filled, x, y, z, &xv, &yv, &zv)) {
                  if (MRIgetVoxVal(mri_filled, xv, yv, zv, 0) == 0) {
                    mriFillRegion(mri_filled, xv, yv, zv, 0, 1, MIN_UNKNOWN_DIST);
                    /* MRIvox(mri_filled, xv, yv, zv) = 1 ;*/
                    nzeros++;
                    label_counts[best_label]++;
                    nfound++;
                  }
                }
              }
              else {
                DiagBreak();
              }
            }
          }
          else {
            /* find best gc with this label */
            if (gcaFindBestSample(gca, x, y, z, best_label, prior_stride / 2, &gcas[nfound]) == NO_ERROR) {
              for (r = 0; r < gca->ninputs; r++) {
                total_means[r] += means[best_label][r];
              }
              label_counts[best_label]++;
              nfound++;
            }
          }
        }
      }
    }
  }

  fprintf(stdout, "total sample mean = %2.1f (%d zeros)\n", total_means[0] / ((float)nfound - nzeros), nzeros);

  if (getenv("GCA_WRITE_CLASS")) {
    int n;
    FILE *fp;

    fp = fopen("classes.dat", "w");
    if (fp) {
      for (n = 0; n < MAX_DIFFERENT_LABELS; n++)
        if (label_counts[n] > 0) {
          fprintf(fp, "%d  %d\n", n, label_counts[n]);
        }
      fclose(fp);
    }
  }

  *pnsamples = nfound;
  fflush(stdout);

  if (mri_vent_dist) MRIfree(&mri_vent_dist);

  MRIfree(&mri_filled);
  return (gcas);
}
GCA_SAMPLE *GCAfindExteriorSamples(
    GCA *gca, int *pnsamples, int min_spacing, float min_prior, int unknown_nbr_spacing, int use_ventricles)
{
  GCA_SAMPLE *gcas;
  int xi, yi, zi, width, height, depth, nfound;
  int nzeros, r, best_label, label;
  float prior_stride, x, y, z, tmp[MAX_GCA_INPUTS], max_prior;
  float min_unknown[MAX_GCA_INPUTS], max_unknown[MAX_GCA_INPUTS], priors[MAX_DIFFERENT_LABELS], var[MAX_GCA_INPUTS],
      means[MAX_DIFFERENT_LABELS][MAX_GCA_INPUTS], vars[MAX_DIFFERENT_LABELS];
  MRI *mri_filled;

  memset(priors, 0, sizeof(priors));

  // make sure the unknowns we are selecting are outside
  // the range of normal tissue
  // intensities (i.e. > min_unknown || < max_unknown)
  GCAlabelMean(gca, Left_Lateral_Ventricle, max_unknown);
  GCAlabelMean(gca, Right_Lateral_Ventricle, tmp);
  for (r = 0; r < gca->ninputs; r++) {
    max_unknown[r] = (max_unknown[r] + tmp[r]) / 2;
  }

  GCAlabelVar(gca, Left_Lateral_Ventricle, var);
  GCAlabelVar(gca, Right_Lateral_Ventricle, tmp);
  for (r = 0; r < gca->ninputs; r++) {
    var[r] = (var[r] + tmp[r]) / 2;
  }
  printf("bounding unknown intensity as < ");
  for (r = 0; r < gca->ninputs; r++) {
    max_unknown[r] = max_unknown[r];
    if (max_unknown[r] < 10) {
      max_unknown[r] = 10;
    }
    printf("%2.1f ", max_unknown[r]);
  }
  printf("or > ");
  GCAlabelMean(gca, Left_Cerebral_White_Matter, min_unknown);
  GCAlabelMean(gca, Left_Cerebral_White_Matter, tmp);
  GCAlabelMean(gca, Left_Cerebral_Cortex, min_unknown);
  GCAlabelMean(gca, Left_Cerebral_Cortex, tmp);
  for (r = 0; r < gca->ninputs; r++) {
    min_unknown[r] = (min_unknown[r] + tmp[r]) / 2;
  }
  GCAlabelVar(gca, Left_Cerebral_White_Matter, var);
  GCAlabelVar(gca, Left_Cerebral_White_Matter, tmp);
  GCAlabelVar(gca, Left_Cerebral_Cortex, var);
  GCAlabelVar(gca, Left_Cerebral_Cortex, tmp);
  for (r = 0; r < gca->ninputs; r++) {
    var[r] = (var[r] + tmp[r]) / 2;
  }

  for (r = 0; r < gca->ninputs; r++) {
    min_unknown[r] = min_unknown[r] + sqrt(var[r]);
    printf("%2.1f ", min_unknown[r]);
  }
  printf("\n");

  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;

  // samples size allocated is prior voxel points
  gcas = (GCA_SAMPLE *)calloc(width * height * depth, sizeof(GCA_SAMPLE));
  if (!gcas)
    ErrorExit(ERROR_NOMEMORY,
              "%s: could not allocate %dK x %d samples\n",
              "GCAfindStableSamples",
              width * height * depth / 1024,
              sizeof(GCA_SAMPLE));

  // create a full size volume prior_width*prior_spacing = image width
  mri_filled = MRIalloc(width * gca->prior_spacing, height * gca->prior_spacing, depth * gca->prior_spacing, MRI_UCHAR);
  // use the mri_prior_ header

  mri_filled->xsize = gca->xsize;
  mri_filled->ysize = gca->ysize;
  mri_filled->zsize = gca->zsize;

  GCAcopyDCToMRI(gca, mri_filled);

  prior_stride = (float)min_spacing / (float)gca->prior_spacing;

  prior_stride = 1.0;
  /* compute the max priors and min variances for each class */
  for (nzeros = nfound = x = 0; x < width; x += prior_stride) {
    fflush(stdout);  // nicknote: this prevents a segfault on Linux PowerPC
    // when -O2 optimization is used w/ gcc 3.3.3

    xi = nint(x);
    for (y = 0; y < height; y += prior_stride) {
      fflush(stdout);  // nicknote: this prevents a segfault on Linux PowerPC
      // when -O2 optimization is used w/ gcc 3.3.3
      yi = nint(y);
      for (z = 0; z < depth; z += prior_stride) {
        zi = nint(z);
        if (xi == Gx && yi == Gy && zi == Gz) {
          DiagBreak();
        }
        gcaRegionStats(gca, x, y, z, prior_stride / 2, priors, vars, means);

        max_prior = 0;
        best_label = 0;
        for (label = 0; label < MAX_DIFFERENT_LABELS; label++) {
          if (priors[label] > max_prior) {
            best_label = label;
            max_prior = priors[label];
          }
        }
        if (best_label != Unknown && (!use_ventricles || !IS_LAT_VENT(best_label))) {
          continue;
        }
        if (priors[best_label] < min_prior) {
          continue;
        }
        if (best_label == 0 && (different_nbr_max_labels(gca, x, y, z, unknown_nbr_spacing, 0) == 0)) {
          continue;
        }

        if (gcaFindBestSample(gca, x, y, z, best_label, prior_stride / 2, &gcas[nfound]) == NO_ERROR) {
          int xv, yv, zv;

          if (best_label != 0 || ((means[0][0] <= max_unknown[0]) || (means[0][0] >= min_unknown[0]))) {
            if (!GCApriorToVoxel(gca, mri_filled, x, y, z, &xv, &yv, &zv)) {
              //              if (MRIvox(mri_filled, xv, yv, zv) == 0)
              {
                mriFillRegion(mri_filled, xv, yv, zv, 0, 1, MIN_UNKNOWN_DIST);
                /* MRIvox(mri_filled, xv, yv, zv) = 1 ;*/
                nzeros++;
                nfound++;
              }
            }
            else {
              DiagBreak();
            }
          }
        }
      }
    }
  }

  *pnsamples = nfound;
  fflush(stdout);

  MRIfree(&mri_filled);
  return (gcas);
}

//////////////////////////////////////////////////////////////////
// write segmented volume for sampled points
// the values are set only at prior_spacing/size
int GCAwriteSamples(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, int nsamples, const char *fname)
{
  int n, xv, yv, zv;
  MRI *mri_dst;

  mri_dst = MRIalloc(mri->width, mri->height, mri->depth, MRI_UCHAR);
  // add header information
  MRIcopyHeader(mri, mri_dst);
  GCAcopyDCToMRI(gca, mri_dst);

  // for each sample
  for (n = 0; n < nsamples; n++) {
    // use the gcas to get prior value to get the volume position
    if (!GCApriorToVoxel(gca, mri_dst, gcas[n].xp, gcas[n].yp, gcas[n].zp, &xv, &yv, &zv)) {
      if (gcas[n].label > 0)  // if label non-zero (!unknown)
      {
        MRIsetVoxVal(mri_dst, xv, yv, zv, 0, gcas[n].label);
      }
      else  // unknown label changed to
      {
        MRIsetVoxVal(mri_dst, xv, yv, zv, 0, Left_undetermined);
      }
      /* Left undetermined - to make it visible */

      /////////////// diagnostics //////////////////////////////////////
      if (DIAG_VERBOSE_ON && Gdiag & DIAG_SHOW)
        printf(
            "label %d: (%d, %d, %d) <-- (%d, %d, %d)\n", gcas[n].label, gcas[n].xp, gcas[n].yp, gcas[n].zp, xv, yv, zv);
      /////////////////////////////////////////////////////////////////
    }
  }
  MRIwrite(mri_dst, fname);
  MRIfree(&mri_dst);
  return (NO_ERROR);
}
/*--------------------------------------------------------------------------*/
/*
  \fn MRI *GCAtoMLLabel(GCA *gca, MRI *mri)
  \brief Creates a seg volume with the most likely a priori label in the GCA
 */
MRI *GCAtoMLLabel(GCA *gca, MRI *mri)
{
  int width, height, depth, x;

  width = gca->node_width;
  height = gca->node_height;
  depth = gca->node_depth;

  if (!mri) {
    mri = MRIallocSequence(width, height, depth, MRI_INT, 1);
    mri->xsize = gca->xsize * gca->node_spacing;
    mri->ysize = gca->ysize * gca->node_spacing;
    mri->zsize = gca->zsize * gca->node_spacing;
  }
  // in order to create the gca volume, the volume must have the same
  // direction cosines
  GCAcopyDCToMRI(gca, mri);

  for (x = 0; x < width; x++) {
    int nmax, y, z, xp, yp, zp, n, xn, yn, zn;
    double pmax;
    GC1D *gc;
    GCA_PRIOR *gcap;
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp)) continue;
        if (GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn)) continue;
        gcap = &gca->priors[xp][yp][zp];
        if (gcap == NULL) continue;
        pmax = 0;
        nmax = -1;
        for (n = 0; n < gcap->nlabels; n++) {
          gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[n]);
          if (!gc) continue;
          if (pmax < gcap->priors[n]) {
            pmax = gcap->priors[n];
            nmax = n;
          }
        }
        if (nmax == -1) continue;
        MRIsetVoxVal(mri, x, y, z, 0, gcap->labels[nmax]);
      }
    }
  }
  return (mri);
}
/*---------------------------------------------------------------------------------*/
/*
  \fn MRI *GCAtoAPrioriMax(GCA *gca, MRI *mri)
  \brief Creates a volume with the probability of the most likely a priori
  label in the GCA
 */
MRI *GCAtoAPrioriMax(GCA *gca, MRI *mri)
{
  int width, height, depth, x;

  width = gca->node_width;
  height = gca->node_height;
  depth = gca->node_depth;

  if (!mri) {
    mri = MRIallocSequence(width, height, depth, MRI_FLOAT, 1);
    mri->xsize = gca->xsize * gca->node_spacing;
    mri->ysize = gca->ysize * gca->node_spacing;
    mri->zsize = gca->zsize * gca->node_spacing;
  }
  // in order to create the gca volume, the volume must have the same
  // direction cosines
  GCAcopyDCToMRI(gca, mri);

  for (x = 0; x < width; x++) {
    int nmax, y, z, xp, yp, zp, n, xn, yn, zn;
    double pmax;
    GC1D *gc;
    GCA_PRIOR *gcap;
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        MRIsetVoxVal(mri, x, y, z, 0, 0);
        if (GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp)) continue;
        if (GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn)) continue;
        gcap = &gca->priors[xp][yp][zp];
        if (gcap == NULL) continue;
        pmax = 0;
        nmax = -1;
        for (n = 0; n < gcap->nlabels; n++) {
          gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[n]);
          if (!gc) continue;
          if (pmax < gcap->priors[n]) {
            pmax = gcap->priors[n];
            nmax = n;
          }
        }
        if (nmax == -1) continue;
        if (gcap->labels[nmax] == 0) continue;
        MRIsetVoxVal(mri, x, y, z, 0, pmax);
      }
    }
  }

  return (mri);
}

/*--------------------------------------------------------------------------------*/
MRI *GCAmri(GCA *gca, MRI *mri)
{
  int frame, width, height, depth, x, y, z, xp, yp, zp, n, xn, yn, zn;
  float val;
  GC1D *gc;
  GCA_PRIOR *gcap;

  if (!mri) {
    mri = MRIallocSequence(gca->node_width, gca->node_height, gca->node_depth, MRI_UCHAR, gca->ninputs);
    mri->xsize = gca->xsize * gca->node_spacing;
    mri->ysize = gca->ysize * gca->node_spacing;
    mri->zsize = gca->zsize * gca->node_spacing;
  }
  // in order to create the gca volume,
  // the volume must have the same direction cosines
  GCAcopyDCToMRI(gca, mri);

  width = mri->width;
  height = mri->height;
  depth = mri->depth;

  for (frame = 0; frame < gca->ninputs; frame++) {
    for (x = 0; x < width; x++) {
      for (y = 0; y < height; y++) {
        for (z = 0; z < depth; z++) {
          if (!GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp))
            if (!GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn)) {
              gcap = &gca->priors[xp][yp][zp];
              if (gcap == NULL) {
                continue;
              }
              for (val = 0.0, n = 0; n < gcap->nlabels; n++) {
                gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[n]);
                if (gc) {
                  val += gc->means[frame] * gcap->priors[n];
                }
              }
              MRIsetVoxVal(mri, x, y, z, frame, val);
            }
        }
      }
    }
  }
  return (mri);
}

MRI *GCAlabelMri(GCA *gca, MRI *mri, int label, TRANSFORM *transform)
{
  int width, height, depth, x, y, z, xn, yn, zn, n, frame;
  float val;
  GC1D *gc = NULL;
  GCA_NODE *gcan;
  MRI *mri_norm;

  if (!mri) {
    mri = MRIallocSequence(gca->node_width, gca->node_height, gca->node_depth, MRI_UCHAR, gca->ninputs);

    mri->xsize = gca->xsize * gca->node_spacing;
    mri->ysize = gca->ysize * gca->node_spacing;
    mri->zsize = gca->zsize * gca->node_spacing;
    GCAcopyDCToMRI(gca, mri);
    // mri=NULL, then use the gca volume
    // mri!=NULL, then use the volume as it is to write the label
  }
  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  mri_norm = MRIalloc(width, height, depth, MRI_SHORT);
  MRIcopyHeader(mri, mri_norm);

  for (frame = 0; frame < gca->ninputs; frame++) {
    MRIclear(mri_norm);
    for (x = 0; x < width; x++) {
      for (y = 0; y < height; y++) {
        for (z = 0; z < depth; z++) {
          if (x == width / 2 && y == height / 2 && z == depth / 2) {
            DiagBreak();
          }
          if (x == 19 && y == 14 && z == 15) {
            DiagBreak();
          }

          if (!GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn)) {
            gcan = &gca->nodes[xn][yn][zn];
            if (gcan->nlabels < 1) {
              continue;
            }
            for (n = 0; n < gcan->nlabels; n++) {
              gc = &gcan->gcs[n];
              if (gcan->labels[n] == label) {
                break;
              }
            }
            if (n >= gcan->nlabels) {
              continue;
            }
            val = gc->means[frame];
            switch (mri->type) {
              default:
                ErrorReturn(NULL,
                            (ERROR_UNSUPPORTED,
                             "GCAlabelMri: unsupported image"
                             " type %d",
                             mri->type));
              case MRI_UCHAR:
                MRIseq_vox(mri, x, y, z, frame) = (unsigned short)val;
                break;
              case MRI_SHORT:
                MRISseq_vox(mri, x, y, z, frame) = (short)val;
                break;
              case MRI_FLOAT:
                MRIFseq_vox(mri, x, y, z, frame) = (float)val;
                break;
            }
            MRISvox(mri_norm, x, y, z) += 1;
          }
        }
      }
    }
    for (x = 0; x < width; x++) {
      for (y = 0; y < height; y++) {
        for (z = 0; z < depth; z++) {
          if (MRISvox(mri_norm, x, y, z) > 0)
            MRISseq_vox(mri, x, y, z, frame) =
                nint((float)MRISseq_vox(mri, x, y, z, frame) / (float)MRISvox(mri_norm, x, y, z));
        }
      }
    }
  }
  MRIfree(&mri_norm);
  return (mri);
}

static int different_nbr_max_labels(GCA *gca, int x, int y, int z, int wsize, int label)
{
  int xk, yk, zk, xi, yi, zi, nbrs, n, width, height, depth, max_label;
  GCA_PRIOR *gcap;
  double max_p;

  // look at prior volume
  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;
  for (nbrs = 0, xk = -wsize; xk <= wsize; xk++) {
    xi = x + xk;
    if (xi < 0 || xi >= width) {
      continue;
    }

    for (yk = -wsize; yk <= wsize; yk++) {
      yi = y + yk;
      if (yi < 0 || yi >= height) {
        continue;
      }

      for (zk = -wsize; zk <= wsize; zk++) {
        zi = z + zk;
        if (zi < 0 || zi >= depth) {
          continue;
        }

        gcap = &gca->priors[xi][yi][zi];
        if (gcap == NULL) {
          continue;
        }
        if (gcap->nlabels == 0)  // no priors nor labels exist
        {
          continue;
        }

        max_p = gcap->priors[0];
        max_label = gcap->labels[0];
        for (n = 1; n < gcap->nlabels; n++)
          if (gcap->priors[n] >= max_p) {
            max_label = gcap->labels[n];
            max_p = gcap->priors[n];
          }
        if (max_label != label)
        // if the label is different from the one given
        {
          nbrs++;  // count them
        }
      }
    }
  }
  return (nbrs);
}
static int gcaRegionStats(GCA *gca,
                          int x,
                          int y,
                          int z,
                          int wsize,
                          float *priors,
                          float *vars,
                          float means[MAX_DIFFERENT_LABELS][MAX_GCA_INPUTS])
{
  int xk, yk, zk, xi, yi, zi, width, height, depth, label, n, total_training, r, l;
  GCA_PRIOR *gcap;
  GC1D *gc;
  float dof;

  if (x == 28 && y == 24 && z == 36) {
    DiagBreak();
  }

  memset(priors, 0, MAX_DIFFERENT_LABELS * sizeof(float));
  memset(vars, 0, MAX_DIFFERENT_LABELS * sizeof(float));
  for (l = 0; l < MAX_DIFFERENT_LABELS; l++)
    for (r = 0; r < gca->ninputs; r++) {
      means[l][r] = 0;
    }

  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;
  total_training = 0;
  for (xk = -wsize; xk <= wsize; xk++) {
    xi = x + xk;
    if (xi < 0 || xi >= width) {
      continue;
    }

    for (yk = -wsize; yk <= wsize; yk++) {
      yi = y + yk;
      if (yi < 0 || yi >= height) {
        continue;
      }

      for (zk = -wsize; zk <= wsize; zk++) {
        zi = z + zk;
        if (zi < 0 || zi >= depth) {
          continue;
        }

        gcap = &gca->priors[xi][yi][zi];
        if (gcap == NULL) {
          continue;
        }
        total_training += gcap->total_training;

        for (n = 0; n < gcap->nlabels; n++) {
          label = gcap->labels[n];
          if (label == Gdiag_no) {
            DiagBreak();
          }
          if (label > MAX_DIFFERENT_LABELS) {
            printf("ERROR: label = %d > %d n = %d, xyzi = %d %d %d\n", label, MAX_DIFFERENT_LABELS, n, xi, yi, zi);
            fflush(stdout);
          }
          gc = GCAfindPriorGC(gca, xi, yi, zi, label);
          if (gc) {
            dof = (gcap->priors[n] * (float)gcap->total_training);

            vars[label] += dof * covariance_determinant(gc, gca->ninputs);
            for (r = 0; r < gca->ninputs; r++) {
              means[label][r] += dof * gc->means[r];
            }
            priors[label] += dof;
          }
        }
      }
    }
  }

  for (label = 0; label < MAX_DIFFERENT_LABELS; label++) {
    dof = priors[label];
    if (dof > 0) {
      vars[label] /= dof;
      priors[label] /= total_training;
      for (r = 0; r < gca->ninputs; r++) {
        means[label][r] /= dof;
      }
    }
  }
  return (NO_ERROR);
}
static int gcaFindBestSample(GCA *gca, int x, int y, int z, int label, int wsize, GCA_SAMPLE *gcas)
{
  int xk, yk, zk, xi, yi, zi, width, height, depth, n, best_x, best_y, best_z, xn, yn, zn, r, c, v;
  // int best_n;
  GCA_PRIOR *gcap;
  GC1D *gc, *best_gc;
  float max_prior, prior;
  // float min_var;

  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;
  max_prior = 0.0;
  // min_var = 10000.0f;
  best_gc = NULL;
  best_x = best_y = best_z = -1;
  // best_n = -1;
  for (xk = -wsize; xk <= wsize; xk++) {
    xi = x + xk;
    if (xi < 0 || xi >= width) {
      continue;
    }

    for (yk = -wsize; yk <= wsize; yk++) {
      yi = y + yk;
      if (yi < 0 || yi >= height) {
        continue;
      }

      for (zk = -wsize; zk <= wsize; zk++) {
        zi = z + zk;
        if (zi < 0 || zi >= depth) {
          continue;
        }

        gcap = &gca->priors[xi][yi][zi];
        if (gcap == NULL) {
          continue;
        }
        for (n = 0; n < gcap->nlabels; n++) {
          if (gcap->labels[n] != label) {
            continue;
          }
          DiagBreak();
          prior = gcap->priors[n];
          if (!GCApriorToNode(gca, xi, yi, zi, &xn, &yn, &zn)) {
            gc = GCAfindGC(gca, xn, yn, zn, label);
            if (!gc) {
              continue;
            }

            if (prior > max_prior
            ) {
              max_prior = prior; /*min_var = gc->var ;*/
              best_gc = gc;
              best_x = xi;
              best_y = yi;
              best_z = zi;
              // best_n = n;
            }
          }
        }
      }
    }
  }

  if (best_x < 0) {
    ErrorPrintf(ERROR_BADPARM, "could not find GC1D for label %d at (%d,%d,%d)\n", label, x, y, z);
    return (ERROR_BADPARM);
  }
  if (best_x == 145 / 4 && best_y == 89 / 4 && best_z == 125 / 4) {
    DiagBreak();
  }
  gcas->xp = best_x;
  gcas->yp = best_y;
  gcas->zp = best_z;
  gcas->label = label;
  gcas->means = (float *)calloc(gca->ninputs, sizeof(float));
  gcas->covars = (float *)calloc((gca->ninputs * (gca->ninputs + 1)) / 2, sizeof(float));
  if (!gcas->means || !gcas->covars)
    ErrorExit(ERROR_NOMEMORY,
              "GCAfindStableSamples: could not allocate "
              "mean (%d) and covariance (%d) matrices",
              gca->ninputs,
              gca->ninputs * (gca->ninputs + 1) / 2);
  for (r = v = 0; r < gca->ninputs; r++) {
    gcas->means[r] = best_gc->means[r];
    for (c = r; c < gca->ninputs; c++, v++) {
      gcas->covars[v] = best_gc->covars[v];
    }
  }
  gcas_setPrior(*gcas, max_prior);

  return (NO_ERROR);
}

// same as GCAtransformAndWriteSamples except here we write
// the -logp into the volume instead of the seg label
int GCAtransformAndWriteSamplePvals(
    GCA *gca, MRI *mri, GCA_SAMPLE *gcas, int nsamples, const char *fname, TRANSFORM *transform)
{
  int i, n, xv, yv, zv;
  MRI *mri_dst;
  float vals[MAX_GCA_INPUTS];
  double pval, pval_total, p;
  GCA_PRIOR *gcap;

  mri_dst = MRIalloc(mri->width, mri->height, mri->depth, MRI_FLOAT);
  MRIcopyHeader(mri, mri_dst);
  TransformSetMRIVolGeomToDst(transform, mri_dst);

  TransformInvert(transform, mri);
  for (n = 0; n < nsamples; n++) {
    if (gcas[n].label == Gdiag_no) {
      DiagBreak();
    }
    if (!GCApriorToSourceVoxel(gca, mri_dst, transform, gcas[n].xp, gcas[n].yp, gcas[n].zp, &xv, &yv, &zv)) {
      gcap = &gca->priors[gcas[n].xp][gcas[n].yp][gcas[n].zp];
      load_vals(mri, xv, yv, zv, vals, gca->ninputs);
      pval = gcaComputeSampleConditionalDensity(&gcas[n], vals, gca->ninputs, gcas[n].label);
      pval *= getPrior(gcap, gcas[n].label);
      pval_total = 0.0;
      for (i = 0; i < gcap->nlabels; i++) {
        p = gcaComputeSampleConditionalDensity(&gcas[n], vals, gca->ninputs, gcap->labels[i]);
        p *= gcap->priors[i];
        pval_total += p;
      }
      if (!DZERO(pval_total)) {
        pval /= pval_total;
      }
      mriFillRegion(mri_dst, xv, yv, zv, 0, pval, 0);

      if (gcas[n].x == Gx && gcas[n].y == Gy && gcas[n].z == Gz) {
        DiagBreak();
      }

      gcas[n].x = xv;
      gcas[n].y = yv;
      gcas[n].z = zv;

      //////////// diag /////////////////////////
      if (gcas[n].x == Gx && gcas[n].y == Gy && gcas[n].z == Gz) {
        DiagBreak();
      }
      if (DIAG_VERBOSE_ON && Gdiag & DIAG_SHOW)
        printf(
            "label %d: (%d, %d, %d) <-- (%d, %d, %d)\n", gcas[n].label, gcas[n].xp, gcas[n].yp, gcas[n].zp, xv, yv, zv);
      //////////////////////////////////////////////
    }
  }
  fprintf(stdout, "writing sample pvals to %s\n", fname);
  MRIwrite(mri_dst, fname);
  MRIfree(&mri_dst);
  fflush(stdout);

  return (NO_ERROR);
}
// same as GCAtransformAndWriteSamples except here we write
// the mean into the volume instead of the seg label
int GCAtransformAndWriteSampleMeans(
    GCA *gca, MRI *mri, GCA_SAMPLE *gcas, int nsamples, const char *fname, TRANSFORM *transform)
{
  int n, xv, yv, zv, f;
  MRI *mri_dst;
  // GCA_PRIOR *gcap;

  mri_dst = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, gca->ninputs);
  MRIcopyHeader(mri, mri_dst);

  TransformInvert(transform, mri);
  for (n = 0; n < nsamples; n++) {
    if (gcas[n].label == Gdiag_no) {
      DiagBreak();
    }
    if (!GCApriorToSourceVoxel(gca, mri_dst, transform, gcas[n].xp, gcas[n].yp, gcas[n].zp, &xv, &yv, &zv)) {
      // gcap = &gca->priors[gcas[n].xp][gcas[n].yp][gcas[n].zp];
      for (f = 0; f < gca->ninputs; f++) {
        mriFillRegion(mri_dst, xv, yv, zv, f, gcas[n].means[f], 0);
      }

      if (gcas[n].x == Gx && gcas[n].y == Gy && gcas[n].z == Gz) {
        DiagBreak();
      }

      gcas[n].x = xv;
      gcas[n].y = yv;
      gcas[n].z = zv;

      //////////// diag /////////////////////////
      if (gcas[n].x == Gx && gcas[n].y == Gy && gcas[n].z == Gz) {
        DiagBreak();
      }
      if (DIAG_VERBOSE_ON && Gdiag & DIAG_SHOW)
        printf(
            "label %d: (%d, %d, %d) <-- (%d, %d, %d)\n", gcas[n].label, gcas[n].xp, gcas[n].yp, gcas[n].zp, xv, yv, zv);
      //////////////////////////////////////////////
    }
  }
  fprintf(stdout, "writing sample pvals to %s\n", fname);
  MRIwrite(mri_dst, fname);
  MRIfree(&mri_dst);
  fflush(stdout);

  return (NO_ERROR);
}
// create a volume mapped to the current mri volume
int GCAtransformAndWriteSamples(
    GCA *gca, MRI *mri, GCA_SAMPLE *gcas, int nsamples, const char *fname, TRANSFORM *transform)
{
  int n, xv, yv, zv, label;
  MRI *mri_dst;

  mri_dst = MRIalloc(mri->width, mri->height, mri->depth, MRI_UCHAR);
  MRIcopyHeader(mri, mri_dst);  // TransformSetMRIVolGeomToDst(transform, mri_dst) ;

  TransformInvert(transform, mri);
  for (n = 0; n < nsamples; n++) {
    if (gcas[n].label == Gdiag_no) {
      DiagBreak();
    }
    if (!GCApriorToSourceVoxel(gca, mri_dst, transform, gcas[n].xp, gcas[n].yp, gcas[n].zp, &xv, &yv, &zv)) {
      if (gcas[n].label > 0) {
        label = gcas[n].label;
      }
      else if (gcas[n].label == 0) {
        label = 29; /* Left undetermined - visible */
      }
      else {
        label = 0; /* Left undetermined - visible */
      }

      mriFillRegion(mri_dst, xv, yv, zv, 0, label, 0);

      if (gcas[n].x == Gx && gcas[n].y == Gy && gcas[n].z == Gz) {
        DiagBreak();
      }

      gcas[n].x = xv;
      gcas[n].y = yv;
      gcas[n].z = zv;

      //////////// diag /////////////////////////
      if (gcas[n].x == Gx && gcas[n].y == Gy && gcas[n].z == Gz) {
        DiagBreak();
      }
      if (DIAG_VERBOSE_ON && Gdiag & DIAG_SHOW)
        printf(
            "label %d: (%d, %d, %d) <-- (%d, %d, %d)\n", gcas[n].label, gcas[n].xp, gcas[n].yp, gcas[n].zp, xv, yv, zv);
      //////////////////////////////////////////////
    }
  }
  fprintf(stdout, "writing samples to %s\n", fname);
  MRIwrite(mri_dst, fname);
  MRIfree(&mri_dst);
  fflush(stdout);

  return (NO_ERROR);
}

static int mriFillRegion(MRI *mri, int x, int y, int z, int frame, float fill_val, int whalf)
{
  int xi, xk, yi, yk, zi, zk;

  for (xk = -whalf; xk <= whalf; xk++) {
    xi = mri->xi[x + xk];
    for (yk = -whalf; yk <= whalf; yk++) {
      yi = mri->yi[y + yk];
      for (zk = -whalf; zk <= whalf; zk++) {
        zi = mri->zi[z + zk];
        MRIsetVoxVal(mri, xi, yi, zi, frame, fill_val);
      }
    }
  }
  return (NO_ERROR);
}

static int gcaFindIntensityBounds(GCA *gca, float *pmin, float *pmax)
{
  int width, depth, height, x, y, z, n, label, r;
  GCA_NODE *gcan;
  GC1D *gc;
  float mn[MAX_GCA_INPUTS], mx[MAX_GCA_INPUTS], offset;

  for (r = 0; r < gca->ninputs; r++) {
    mn[r] = 100000.0f;
    mx[r] = -mn[r];
  }
  width = gca->node_width;
  height = gca->node_height;
  depth = gca->node_depth;

  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          gc = &gcan->gcs[n];
          label = gcan->labels[n];
          if (label == 0) /* don't include unknowns */
          {
            continue;
          }
          if (get_node_prior(gca, label, x, y, z) < 0.1) {
            continue;
          }
          /* exclude unlikely stuff (errors in labeling) */
          offset = 0.5 * pow(covariance_determinant(gc, gca->ninputs), 1.0 / (double)gca->ninputs);
          for (r = 0; r < gca->ninputs; r++) {
            if (gc->means[r] + offset > mx[r]) {
              mx[r] = gc->means[r] + offset;
            }
            if (gc->means[r] < mn[r]) {
              mn[r] = gc->means[r];
            }
          }
        }
      }
    }
  }
  for (r = 0; r < gca->ninputs; r++) {
    pmin[r] = mn[r];
    pmax[r] = mx[r];
  }
  return (NO_ERROR);
}
static int gcaFindMaxPriors(GCA *gca, float *max_priors)
{
  int width, depth, height, x, y, z, n, label;
  GCA_PRIOR *gcap;

  memset(max_priors, 0, MAX_DIFFERENT_LABELS * sizeof(float));
  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;

  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        gcap = &gca->priors[x][y][z];
        if (gcap == NULL || gcap->nlabels <= 0) {
          continue;
        }
        for (n = 0; n < gcap->nlabels; n++) {
          label = gcap->labels[n];
          if (gcap->priors[n] > max_priors[label]) {
            if (label == Ggca_label) {
              DiagBreak();
            }
            if ((gcap->priors[n] > 0.89) && (label == Ggca_label)) {
              DiagBreak();
            }
            max_priors[label] = gcap->priors[n];
          }
        }
      }
    }
  }
  return (NO_ERROR);
}


//                                 segmented volume here
static int GCAupdateNodeGibbsPriors(GCA *gca, MRI *mri, int xn, int yn, int zn, int xl, int yl, int zl, int label)
{
  int n, i, xnbr, ynbr, znbr, nbr_label;
  GCA_NODE *gcan;
  GC1D *gc;

  gcan = &gca->nodes[xn][yn][zn];

  // look for this label
  for (n = 0; n < gcan->nlabels; n++) {
    if (gcan->labels[n] == label) {
      break;
    }
  }
  // not found
  if (n >= gcan->nlabels) /* have to allocate a new classifier */
    ErrorExit(ERROR_BADPARM, "gca(%d, %d, %d): could not find label %d", xn, yn, zn, label);

  gc = &gcan->gcs[n];
  for (i = 0; i < GIBBS_NEIGHBORHOOD; i++) {
    /* coordinates of neighboring point */
    xnbr = mri->xi[xl + xnbr_offset[i]];  // xnbr_offset 1, -1, 0,  0, 0,  0
    ynbr = mri->yi[yl + ynbr_offset[i]];  // ynbr_offset 0,  0, 1, -1, 0,  0
    znbr = mri->zi[zl + znbr_offset[i]];  // znbr_offset 0,  0, 0,  0, 1, -1

    // get the label from the neighbor
    nbr_label = nint(MRIgetVoxVal(mri, xnbr, ynbr, znbr, 0));
    if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z && label == Ggca_label && (xnbr_offset[i] == 1) &&
        (nbr_label == Ggca_nbr_label)) {
      printf(
          "(%d, %d, %d) --> node (%d, %d, %d), "
          "label %s (%d), nbr (%d, %d, %d) = %s (%d)\n",
          xl,
          yl,
          zl,
          xn,
          yn,
          zn,
          cma_label_to_name(label),
          label,
          xnbr,
          ynbr,
          znbr,
          cma_label_to_name(nbr_label),
          nbr_label);
    }

    /* now see if this label exists already as a nbr */
    for (n = 0; n < gc->nlabels[i]; n++)
      if (gc->labels[i][n] == nbr_label) {
        break;
      }

    if (n >= gc->nlabels[i]) /* not there - reallocate stuff */
    {
      unsigned short *old_labels;
      float *old_label_priors;

      old_labels = gc->labels[i];
      old_label_priors = gc->label_priors[i];

      /* allocate new ones */
      gc->label_priors[i] = (float *)calloc(gc->nlabels[i] + 1, sizeof(float));
      if (!gc->label_priors[i])
        ErrorExit(ERROR_NOMEMORY,
                  "GCAupdateNodeGibbsPriors: "
                  "couldn't expand gcs to %d",
                  gc->nlabels[i] + 1);
      gc->labels[i] = (unsigned short *)calloc(gc->nlabels[i] + 1, sizeof(unsigned short));
      if (!gc->labels[i]) ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand labels to %d", gc->nlabels[i] + 1);

      if (gc->nlabels[i] > 0) /* copy the old ones over */
      {
        memmove(gc->label_priors[i], old_label_priors, gc->nlabels[i] * sizeof(float));
        memmove(gc->labels[i], old_labels, gc->nlabels[i] * sizeof(unsigned short));

        /* free the old ones */
        free(old_label_priors);
        free(old_labels);
      }
      gc->labels[i][gc->nlabels[i]++] = nbr_label;
    }
    gc->label_priors[i][n] += 1.0f;  // counter increment at this label
  }

  return (NO_ERROR);
}


MRI *GCAanneal(MRI *mri_inputs, GCA *gca, MRI *mri_dst, TRANSFORM *transform, int max_iter, double prior_factor)
{
  int x, y, z, width, height, depth, label, iter, xn, yn, zn, n, nchanged, index, nindices, old_label;
  // int  val;
  short *x_indices, *y_indices, *z_indices;
  GCA_NODE *gcan;
  double ll, lcma = 0.0, old_posterior, new_posterior, min_posterior;
  MRI *mri_changed, *mri_probs;

  printf("performing simulated annealing...\n");

  nindices = mri_dst->width * mri_dst->height * mri_dst->depth;
  x_indices = (short *)calloc(nindices, sizeof(short));
  y_indices = (short *)calloc(nindices, sizeof(short));
  z_indices = (short *)calloc(nindices, sizeof(short));
  if (!x_indices || !y_indices || !z_indices)
    ErrorExit(ERROR_NOMEMORY,
              "GCAanneal: "
              "could not allocate index set");

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  iter = 0;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_UCHAR);
    if (!mri_dst) {
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst");
    }
    MRIcopyHeader(mri_inputs, mri_dst);
  }

  mri_changed = MRIclone(mri_dst, NULL);

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  for (x = 0; x < width; x++)
    for (y = 0; y < height; y++)
      for (z = 0; z < depth; z++) {
        MRIsetVoxVal(mri_changed, x, y, z, 0, 1);
      }

  old_posterior = GCAgibbsImageLogPosterior(gca, mri_dst, mri_inputs, transform, prior_factor);
  old_posterior /= (double)(width * depth * height);
  fflush(stdout);

  do {
    if (iter == 0) {
      mri_probs = GCAlabelProbabilities(mri_inputs, gca, mri_dst, transform);
      MRIorderIndices(mri_probs, x_indices, y_indices, z_indices);
      MRIfree(&mri_probs);
    }
    else
      MRIcomputeVoxelPermutation(mri_inputs, x_indices, y_indices, z_indices);

    nchanged = 0;
    for (index = 0; index < nindices; index++) {
      x = x_indices[index];
      y = y_indices[index];
      z = z_indices[index];
      if (x == 155 && y == 126 && z == 128) {
        DiagBreak();
      }

      if ((int)MRIgetVoxVal(mri_changed, x, y, z, 0) == 0) {
        continue;
      }

      if (x == 63 && y == 107 && z == 120) {
        DiagBreak();
      }

      // val =
      MRIgetVoxVal(mri_inputs, x, y, z, 0);

      /* find the node associated with this coordinate and classify */
      if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
        gcan = &gca->nodes[xn][yn][zn];

        label = old_label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0));
        min_posterior = GCAnbhdGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, prior_factor);

        for (n = 0; n < gcan->nlabels; n++) {
          if (gcan->labels[n] == old_label) {
            continue;
          }
          MRIsetVoxVal(mri_dst, x, y, z, 0, gcan->labels[n]);
          new_posterior = GCAnbhdGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, prior_factor);
          if (new_posterior > min_posterior) {
            min_posterior = new_posterior;
            label = gcan->labels[n];
          }
        }
        if (label != old_label) {
          nchanged++;
          MRIsetVoxVal(mri_changed, x, y, z, 0, 1);
        }
        else {
          MRIsetVoxVal(mri_changed, x, y, z, 0, 0);
        }
        MRIsetVoxVal(mri_dst, x, y, z, 0, label);
      }
    }  // index loop
    if (nchanged > 10000) {
      ll = GCAgibbsImageLogPosterior(gca, mri_dst, mri_inputs, transform, prior_factor);
      ll /= (double)(width * depth * height);
      if (!FZERO(lcma))
        printf("pass %d: %d changed. image ll: %2.3f (CMA=%2.3f)\n", iter + 1, nchanged, ll, lcma);
      else
        printf("pass %d: %d changed. image ll: %2.3f\n", iter + 1, nchanged, ll);
    }
    else {
      printf("pass %d: %d changed.\n", iter + 1, nchanged);
    }
    MRIdilate(mri_changed, mri_changed);
  } while (nchanged > 0 && iter++ < max_iter);

  free(x_indices);
  free(y_indices);
  free(z_indices);
  MRIfree(&mri_changed);

  return (mri_dst);
}

char *gca_write_fname = NULL;
int gca_write_iterations = 0;


MRI *GCAreclassifyUsingGibbsPriors(MRI *mri_inputs,
                                   GCA *gca,
                                   MRI *mri_dst,
                                   TRANSFORM *transform,
                                   int max_iter,
                                   MRI *mri_fixed,
                                   int restart,
                                   void (*update_func)(MRI *),
                                   double min_prior_factor,
                                   double max_prior_factor)
{
  int x, y, z, width, height, depth, iter, nchanged, min_changed, index, nindices, fixed;
  short *x_indices, *y_indices, *z_indices;
  double prior_factor, old_posterior, lcma = 0.0;
  MRI *mri_changed, *mri_probs /*, *mri_zero */;

  prior_factor = min_prior_factor;
  // fixed is the label fixed volume, e.g. wm
  fixed = (mri_fixed != NULL);


  nindices = mri_dst->width * mri_dst->height * mri_dst->depth;
  x_indices = (short *)calloc(nindices, sizeof(short));
  y_indices = (short *)calloc(nindices, sizeof(short));
  z_indices = (short *)calloc(nindices, sizeof(short));
  if (!x_indices || !y_indices || !z_indices)
    ErrorExit(ERROR_NOMEMORY,
              "GCAreclassifyUsingGibbsPriors: "
              "could not allocate index set");

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  iter = 0;
  if (!mri_dst) {
    mri_dst = MRIalloc(width, height, depth, MRI_INT);
    if (!mri_dst) {
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst");
    }
    MRIcopyHeader(mri_inputs, mri_dst);
  }

  mri_changed = MRIclone(mri_dst, NULL);

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  // mark changed location
  // ifdef HAVE_OPENMP
  // pragma omp parallel for if_ROMP(experimental)
  // endif
  for (x = 0; x < width; x++) {
    int y, z;
    for (y = 0; y < height; y++)
      for (z = 0; z < depth; z++) {
        if (restart && mri_fixed) {
          // mark only fixed voxels
          if (MRIgetVoxVal(mri_fixed, x, y, z, 0)) {
            MRIsetVoxVal(mri_changed, x, y, z, 0, 1);
          }
        }
        else
        // everything is marked
        {
          MRIsetVoxVal(mri_changed, x, y, z, 0, 1);
        }
      }
  }
  if (restart && mri_fixed) {
    MRIdilate(mri_changed, mri_changed);
  }

  if (!restart) {
    // calculate the statistics
    old_posterior = GCAgibbsImageLogPosterior(gca, mri_dst, mri_inputs, transform, prior_factor);
    // get the per voxel value
    old_posterior /= (double)(width * depth * height);
  }
  else {
    old_posterior = 0;
  }


  prior_factor = min_prior_factor;
  do {
    if (restart) {
      for (index = x = 0; x < width; x++)
        for (y = 0; y < height; y++)
          for (z = 0; z < depth; z++) {
            // not fixed voxel, but changed
            if ((int)MRIgetVoxVal(mri_fixed, x, y, z, 0) == 0 && ((int)MRIgetVoxVal(mri_changed, x, y, z, 0) > 0)) {
              x_indices[index] = x;
              y_indices[index] = y;
              z_indices[index] = z;
              index++;
            }
          }
      nindices = index;
      mri_probs = NULL;
    }
    else if (iter == 0) {
      /*      char  fname[STRLEN], *cp ;*/
      /*      int   nfixed ;*/

      if (gca_write_iterations) {
        char fname[STRLEN];
        sprintf(fname, "%s%03d.mgz", gca_write_fname, iter);
        printf("writing snapshot to %s\n", fname);
        MRIwrite(mri_dst, fname);
      }
      // probs has 0 to 255 values
      mri_probs = GCAlabelProbabilities(mri_inputs, gca, NULL, transform);
      // sorted according to ascending order of probs
      MRIorderIndices(mri_probs, x_indices, y_indices, z_indices);
      MRIfree(&mri_probs);
    }
    else
      // randomize the indices value ((0 -> width*height*depth)
      MRIcomputeVoxelPermutation(mri_inputs, x_indices, y_indices, z_indices);

    ///////////////  relabel with neighborhood likelihood  //////////////
    nchanged = 0;
    if (G_write_probs && !mri_probs) {
      mri_probs = MRIalloc(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT);
      MRIcopyHeader(mri_inputs, mri_probs);
    }

    //#ifdef HAVE_OPENMP
    // pragma omp parallel for if_ROMP(experimental) reduction(+: nchanged)
    //#endif
    for (index = 0; index < nindices; index++) {
      int x, y, z, n, label, old_label;
      GCA_PRIOR *gcap;
      double new_posterior, max_posterior;
      // float val;

      x = x_indices[index];
      y = y_indices[index];
      z = z_indices[index];
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z) DiagBreak();

      // if the label is fixed, don't do anything
      if (mri_fixed && MRIgetVoxVal(mri_fixed, x, y, z, 0)) continue;

      // if not marked, don't do anything
      if (MRIgetVoxVal(mri_changed, x, y, z, 0) == 0) continue;

      // get the grey value
      // val =
      MRIgetVoxVal(mri_inputs, x, y, z, 0);

      /* find the node associated with this coordinate and classify */
      gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
      // it is not in the right place
      if (gcap == NULL) continue;

      // only one label associated, don't do anything
      if (gcap->nlabels == 1) continue;

      // save the current label
      label = old_label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0));
      // calculate neighborhood likelihood
      max_posterior = GCAnbhdGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, prior_factor);

      // go through all labels at this point
      for (n = 0; n < gcap->nlabels; n++) {
        // skip the current label
        if (gcap->labels[n] == old_label) continue;

        // assign the new label
        MRIsetVoxVal(mri_dst, x, y, z, 0, gcap->labels[n]);
        // calculate neighborhood likelihood
        new_posterior = GCAnbhdGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, prior_factor);
        // if it is bigger than the old one, then replace the label
        // and change max_posterior
        if (new_posterior > max_posterior) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z &&
              (label == Ggca_label || old_label == Ggca_label || Ggca_label < 0))
            fprintf(stdout,
                    "NbhdGibbsLogLikelihood at (%d, %d, %d):"
                    " old = %d (ll=%.2f) new = %d (ll=%.2f)\n",
                    x,
                    y,
                    z,
                    old_label,
                    max_posterior,
                    gcap->labels[n],
                    new_posterior);

          max_posterior = new_posterior;
          label = gcap->labels[n];
        }
      }

      /*#ifndef __OPTIMIZE__*/
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z &&
          (label == Ggca_label || old_label == Ggca_label || Ggca_label < 0)) {
        int xn, yn, zn;
        GCA_NODE *gcan;

        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          gcan = &gca->nodes[xn][yn][zn];
          printf(
              "(%d, %d, %d): old label %s (%d), "
              "new label %s (%d) (log(p)=%2.3f)\n",
              x,
              y,
              z,
              cma_label_to_name(old_label),
              old_label,
              cma_label_to_name(label),
              label,
              max_posterior);
          dump_gcan(gca, gcan, stdout, 0, gcap);
          if (label == Right_Caudate) {
            DiagBreak();
          }
        }
      }
      /*#endif*/

      // if label changed
      if (label != old_label) {
        nchanged++;
        // mark it as changed
        MRIsetVoxVal(mri_changed, x, y, z, 0, 1);
      }
      else {
        MRIsetVoxVal(mri_changed, x, y, z, 0, 0);
      }
      // assign new label
      MRIsetVoxVal(mri_dst, x, y, z, 0, label);
      if (mri_probs) {
        MRIsetVoxVal(mri_probs, x, y, z, 0, -max_posterior);
      }
    }
    if (mri_probs) {
      char fname[STRLEN];

      sprintf(fname, "%s%03d.mgz", G_write_probs, iter);
      printf("writing probabilities to %s\n", fname);
      MRIwrite(mri_probs, fname);
      MRIfree(&mri_probs);
    }

    if (0)
    // reclassify connected compontents of unlikely labels as a whole
    {
      MRI_SEGMENTATION *mriseg;
      MRI *mri_impossible, *mri_impossible_label;
      int label, s, nchanged, iter, max_label;
      char fname[STRLEN];
      double old_posterior, new_posterior;

      max_label = GCAmaxLabel(gca);
      old_posterior = GCAgibbsImageLogPosterior(gca, mri_dst, mri_inputs, transform, prior_factor) /
                      (double)(width * depth * height);
      iter = 0;
      printf("%02d: ll %2.5f\n", iter, old_posterior);

      do {
        nchanged = 0;
        mri_impossible = GCAmarkImpossible(gca, mri_dst, NULL, transform);
        if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
          MRIwrite(mri_impossible, "imp.mgz");
        }
        mri_impossible_label = MRIclone(mri_impossible, NULL);
        if (gca_write_iterations > 0) {
          char fname[STRLEN];
          static int fno = 0;

          sprintf(fname, "%s_before%03d.mgz", gca_write_fname, fno + 1);
          fno++;
          printf("writing snapshot to %s\n", fname);
          MRIwrite(mri_dst, fname);
        }
        for (label = 1; label <= max_label; label++) {
          MRIclear(mri_impossible_label);

          // find voxels that have label and aren't possible
          if (MRIcopyLabeledVoxels(mri_impossible, mri_dst, mri_impossible_label, label) == 0) {
            continue;  // no voxels in label
          }

          sprintf(fname, "mimp_label%d.mgz", label);
          if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
            MRIwrite(mri_impossible_label, fname);
          }
          mriseg = MRIsegment(mri_impossible_label, 0.5, 256);

          if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
            printf(
                "%d low probability segments "
                "found for label %s (%d)...\n",
                mriseg->nsegments,
                cma_label_to_name(label),
                label);

          for (s = 0; s < mriseg->nsegments; s++) {
            if (gcaRelabelSegment(gca, transform, mri_inputs, mri_dst, &mriseg->segments[s]) > 0) {
              MRIsetSegmentValue(mri_changed, mriseg, s, 1);
              nchanged++;
            }
          }
          MRIsegmentFree(&mriseg);
        }
        MRIfree(&mri_impossible_label);
        MRIfree(&mri_impossible);
        new_posterior = GCAgibbsImageLogPosterior(gca, mri_dst, mri_inputs, transform, prior_factor) /
                        (double)(width * depth * height);
        printf("%02d: ll %2.5f (%d segments changed)\n", iter + 1, new_posterior, nchanged);
      } while ((nchanged > 0) && (iter++ < 10));
    }

    /////////////////////////////////////////////////////////////////////////
    // print info
    if (nchanged > 10000 && iter < 2 && !restart) {
      double ll;

      ll = GCAgibbsImageLogPosterior(gca, mri_dst, mri_inputs, transform, prior_factor);
      // get the average value
      ll /= (double)(width * depth * height);
      if (!FZERO(lcma))
        printf(
            "pass %d: %d changed. image ll: %2.3f "
            "(CMA=%2.3f), PF=%2.3f\n",
            iter + 1,
            nchanged,
            ll,
            lcma,
            prior_factor);
      else  // currently this is executed
        printf("pass %d: %d changed. image ll: %2.3f, PF=%2.3f\n", iter + 1, nchanged, ll, prior_factor);
    }
    else {
      printf("pass %d: %d changed.\n", iter + 1, nchanged);
    }

    // get the largest 6 neighbor values,
    // that is, originally 0 could become 1
    MRIdilate(mri_changed, mri_changed);
    // if unpdate_func is present, use it
    if (update_func) {
      (*update_func)(mri_dst);
    }

#define MIN_CHANGED 5000
    min_changed = restart ? 0 : MIN_CHANGED;
    if (nchanged <= min_changed || (restart && iter >= max_iter)) {
      if (restart && restart != GCA_RESTART_ONCE) {
        iter = 0;
      }

      if (!restart) /* examine whole volume next time */
      {
        for (x = 0; x < width; x++)
          for (y = 0; y < height; y++)
            for (z = 0; z < depth; z++) {
              MRIsetVoxVal(mri_changed, x, y, z, 0, 1);
            }
      }
      if (fixed && !restart) {
        printf("removing fixed flag...\n");
        if (mri_fixed) {
          MRIclear(mri_fixed);
        }
        fixed = 0;
      }
      else {
        prior_factor *= 2;
        if (prior_factor < max_prior_factor) fprintf(stdout, "setting PRIOR_FACTOR to %2.4f\n", prior_factor);
      }
      if (gca_write_iterations > 0) {
        char fname[STRLEN];

        sprintf(fname, "%s_iter%d.mgz", gca_write_fname, iter + 1);
        printf("writing snapshot to %s\n", fname);
        MRIwrite(mri_dst, fname);
      }
    }
    if ((gca_write_iterations > 0) && !(iter % gca_write_iterations)) {
      char fname[STRLEN];
      sprintf(fname, "%s%03d.mgz", gca_write_fname, iter + 1);
      printf("writing snapshot to %s\n", fname);
      MRIwrite(mri_dst, fname);
    }
  } while ((nchanged > min_changed || prior_factor < max_prior_factor) && (iter++ < max_iter));


  if (mri_probs) {
    MRIfree(&mri_probs);
  }

  free(x_indices);
  free(y_indices);
  free(z_indices);
  MRIfree(&mri_changed);

  return (mri_dst);
}


double GCAgibbsImageLogPosterior(GCA *gca, MRI *mri_labels, MRI *mri_inputs, TRANSFORM *transform, double prior_factor)
{
  int x, width, depth, height;
  double total_log_posterior;

  width = mri_labels->width;
  height = mri_labels->height;
  depth = mri_labels->depth;

  total_log_posterior = 0.0;
  // ifdef HAVE_OPENMP
  // pragma omp parallel for if_ROMP(experimental) reduction(+: total_log_posterior)
  // endif
  for (x = 0; x < width; x++) {
    int y, z;
    double log_posterior;
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        log_posterior = GCAvoxelGibbsLogPosterior(gca, mri_labels, mri_inputs, x, y, z, transform, prior_factor);
        if (check_finite("gcaGibbsImageLogposterior", log_posterior) == 0) {
          DiagBreak();
        }
        total_log_posterior += log_posterior;
        if (total_log_posterior > 1e10) {
          DiagBreak();
        }
      }
    }
  }
  return (total_log_posterior);
}
static double gcaGibbsImpossibleConfiguration(GCA *gca, MRI *mri_labels, int x, int y, int z, TRANSFORM *transform)
{
  int xn, yn, zn, xnbr, ynbr, znbr, nbr_label, label, i, j, n;
  GCA_NODE *gcan;
  GC1D *gc;

  label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));

  /* find the node associated with this coordinate and classify */
  if (!GCAsourceVoxelToNode(gca, mri_labels, transform, x, y, z, &xn, &yn, &zn)) {
    gcan = &gca->nodes[xn][yn][zn];

    for (n = 0; n < gcan->nlabels; n++) {
      if (gcan->labels[n] == label) {
        break;
      }
    }
    if (n >= gcan->nlabels) {
      return (1); /* never occurred */
    }

    gc = &gcan->gcs[n];

    for (i = 0; i < GIBBS_NEIGHBORS; i++) {
      xnbr = mri_labels->xi[x + xnbr_offset[i]];
      ynbr = mri_labels->yi[y + ynbr_offset[i]];
      znbr = mri_labels->zi[z + znbr_offset[i]];
      nbr_label = nint(MRIgetVoxVal(mri_labels, xnbr, ynbr, znbr, 0));
      for (j = 0; j < gc->nlabels[i]; j++) {
        if (nbr_label == gc->labels[i][j]) {
          break;
        }
      }
      if (j < gc->nlabels[i]) {
        if (FZERO(gc->label_priors[i][j])) {
          return (1); /* never occurred (and this never should) */
        }
      }
      else /* never occurred - make it unlikely */
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        if (FZERO(gc->label_priors[i][j])) {
          return (1); /* never occurred (and this never should) */
        }
      }
    }
  }
  return (0);
}

double GCAnbhdGibbsLogPosterior(
    GCA *gca, MRI *mri_labels, MRI *mri_inputs, int x, int y, int z, TRANSFORM *transform, double gibbs_coef)
{
  double total_log_posterior /*, log_posterior*/;
  /*  int    i, xnbr, ynbr, znbr ;*/

  total_log_posterior = GCAvoxelGibbsLogPosterior(gca, mri_labels, mri_inputs, x, y, z, transform, gibbs_coef);


  if (check_finite("GCAnbhdGibbsLogPosterior: final", total_log_posterior) == 0) {
    DiagBreak();
  }

  return (total_log_posterior);
}

double GCAwindowPosteriorLogProbability(
    GCA *gca, MRI *mri_labels, MRI *mri_inputs, TRANSFORM *transform, int x, int y, int z, int whalf)
{
  int xi, yi, zi, xk, yk, zk, nvox;
  double total_posterior, posterior;

  for (total_posterior = 0.0, nvox = 0, xk = -whalf; xk <= whalf; xk++) {
    xi = mri_labels->xi[x + xk];
    for (yk = -whalf; yk <= whalf; yk++) {
      yi = mri_labels->yi[y + yk];
      for (zk = -whalf; zk <= whalf; zk++) {
        zi = mri_labels->zi[z + zk];
        if (xi == Ggca_x && yi == Ggca_y && zi == Ggca_z) DiagBreak();
        posterior = GCAvoxelGibbsLogPosterior(gca, mri_labels, mri_inputs, xi, yi, zi, transform, 1);
        total_posterior += posterior;
        nvox++;
      }
    }
  }

  return (total_posterior / nvox);
}

double GCAvoxelGibbsLogPosterior(
    GCA *gca, MRI *mri_labels, MRI *mri_inputs, int x, int y, int z, TRANSFORM *transform, double gibbs_coef)
{
  double log_posterior /*, dist*/, nbr_prior;
  int xn, yn, zn, xnbr, ynbr, znbr, nbr_label, label, i, j, n;
  GCA_NODE *gcan = 0;
  GCA_PRIOR *gcap = 0;
  GC1D *gc = 0;
  float vals[MAX_GCA_INPUTS];
#if INTERP_PRIOR
  float prior;
#endif
  // float     tmp = 0;

  // signify error
  log_posterior = 0.;

  // get the grey value
  load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
  // get the label
  label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));
  // what happens with higher number > CMA_MAX?
  /* find the node associated with this coordinate and classify */
  if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
    gcan = &gca->nodes[xn][yn][zn];
    gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
    if (gcap == NULL || gcap->nlabels <= 0) {
      if (label == Unknown)  // okay for there to be an
                             // unknown label out of the fov
      {
        return (0.0);
      }
      else {
        return (10 * BIG_AND_NEGATIVE);
      }
    }

    ////////////////// debug code (this should not occur ) /////////
    if (label > MAX_CMA_LABEL) {
      printf(
          "\nGCAvoxelGibbsLogPosterior() is called "
          "with label %d at (%d, %d, %d)\n",
          label,
          x,
          y,
          z);
      printf("gcan = %p, gcap = %p\n", gcan, gcap);
      if (gcan) {
        printf("gcan->nlabels = %d, gcan->total_training = %d ", gcan->nlabels, gcan->total_training);
        printf("log(return) = %.2f\n", log(0.01f / ((float)gcan->total_training * GIBBS_NEIGHBORS)));
        printf("labels for this location\n");
        for (n = 0; n < gcan->nlabels; n++)
          printf("label=%s (%d); ", cma_label_to_name(gcan->labels[n]), gcan->labels[n]);
      }
    }
    /////////////////////////////////////////////////////////////////
    for (n = 0; n < gcan->nlabels; n++) {
      if (gcan->labels[n] == label) {
        break;
      }
    }
    // could not find the label, then
    if (n >= gcan->nlabels) {
      gc = GCAfindClosestValidGC(gca, xn, yn, zn, label, 0);
    }
    else {
      gc = &gcan->gcs[n];
    }
    if (gc == NULL) {
      // if (gcan->total_training > 0)
      // return(log(0.01f/((float)gcan->total_training*GIBBS_NEIGHBORS))) ;
      /* 10*GIBBS_NEIGHBORS*BIG_AND_NEGATIVE*/
      // else
      return (10 * BIG_AND_NEGATIVE);
      // return(log(VERY_UNLIKELY)) ;
    }

    /* compute 1-d Mahalanobis distance */
    log_posterior = GCAcomputeConditionalLogDensity(gc, vals, gca->ninputs, label);
    if (check_finite("GCAvoxelGibbsLogPosterior: conditional log density", log_posterior) == 0) {
      DiagBreak();
    }

    nbr_prior = 0.0;
    if (gc->nlabels == NULL) {
      nbr_prior += log(0.1f / (float)gcan->total_training);
    }
    else {
      for (i = 0; i < GIBBS_NEIGHBORS; i++) {
        xnbr = mri_labels->xi[x + xnbr_offset[i]];
        ynbr = mri_labels->yi[y + ynbr_offset[i]];
        znbr = mri_labels->zi[z + znbr_offset[i]];
        nbr_label = nint(MRIgetVoxVal(mri_labels, xnbr, ynbr, znbr, 0));
        for (j = 0; j < gc->nlabels[i]; j++) {
          if (nbr_label == gc->labels[i][j]) {
            break;
          }
        }
        if (j < gc->nlabels[i]) {
          if (!FZERO(gc->label_priors[i][j])) {
            nbr_prior += log(gc->label_priors[i][j]);
          }
          else {
            nbr_prior += log(0.1f / (float)gcan->total_training);
          }
          /*BIG_AND_NEGATIVE */
          check_finite("GCAvoxelGibbsLogPosterior: label_priors", nbr_prior);
        }
        else /* never occurred - make it unlikely */
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          nbr_prior += log(0.1f / (float)gcan->total_training);
          /*BIG_AND_NEGATIVE*/
        }
      }
    }
// added to the previous value
#if INTERP_PRIOR
    prior = gcaComputePrior(gca, mri_inputs, transform, x, y, z, label);
    log_posterior += (gibbs_coef * nbr_prior + log(prior));
#else
    log_posterior += (gibbs_coef * nbr_prior + log(getPrior(gcap, label)));
#endif
    if (check_finite("GCAvoxelGibbsLogPosterior: final", log_posterior) == 0) {
      DiagBreak();
    }
  }
  else {
    return (10 * BIG_AND_NEGATIVE);
    // return (log(VERY_UNLIKELY)) ;
  }

// just check
  return (log_posterior);
}
// the posterior of an image given a segmentation without any MRF
double GCAimagePosteriorLogProbability(GCA *gca, MRI *mri_labels, MRI *mri_inputs, TRANSFORM *transform)
{
  double log_posterior;
  int x, y, z, num = 0;

  TransformInvert(transform, mri_inputs);
  for (log_posterior = 0.0, x = 0; x < mri_labels->width; x++)
    for (y = 0; y < mri_labels->height; y++)
      for (z = 0; z < mri_labels->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if ((int)MRIgetVoxVal(mri_labels, x, y, z, 0) == 0) {
          continue;
        }
        log_posterior += GCAvoxelLogPosterior(gca, mri_labels, mri_inputs, x, y, z, transform);
        num++;
        if (!std::isfinite(log_posterior)) {
          DiagBreak();
        }
      }
  return (log_posterior / (float)num);
}

double GCAvoxelLogPosterior(GCA *gca, MRI *mri_labels, MRI *mri_inputs, int x, int y, int z, TRANSFORM *transform)
{
  double log_posterior /*, dist*/;
  int xn, yn, zn, label, n;
  GCA_NODE *gcan = 0;
  GCA_PRIOR *gcap = 0;
  GC1D *gc = 0;
  float vals[MAX_GCA_INPUTS];
#if INTERP_PRIOR
  float prior;
#endif
  // float     tmp = 0;

  // signify error
  log_posterior = 0.;

  load_vals(mri_inputs, x, y, z, vals, gca->ninputs);  // get the grey value
  label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));  // get the label
  // what happens with higher number > CMA_MAX?
  /* find the node associated with this coordinate and classify */
  if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
    gcan = &gca->nodes[xn][yn][zn];
    gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
    if (gcap == NULL || gcap->nlabels <= 0) {
      if (label == Unknown)  // okay for there to be an
      {
        return (0.0);  // unknown label out of the fov
      }
      else {
        return (10 * BIG_AND_NEGATIVE);
      }
    }

    ////////////////// debug code (this should not occur ) /////////
    if (label > MAX_CMA_LABEL) {
      printf(
          "\nGCAvoxelGibbsLogPosterior() is called "
          "with label %d at (%d, %d, %d)\n",
          label,
          x,
          y,
          z);
      printf("gcan = %p, gcap = %p\n", gcan, gcap);
      if (gcan) {
        printf("gcan->nlabels = %d, gcan->total_training = %d ", gcan->nlabels, gcan->total_training);
        printf("log(return) = %.2f\n", log(0.01f / ((float)gcan->total_training * GIBBS_NEIGHBORS)));
        printf("labels for this location\n");
        for (n = 0; n < gcan->nlabels; n++)
          printf("label=%s (%d); ", cma_label_to_name(gcan->labels[n]), gcan->labels[n]);
      }
    }
    /////////////////////////////////////////////////////////////////
    for (n = 0; n < gcan->nlabels; n++) {
      if (gcan->labels[n] == label) {
        break;
      }
    }
    // could not find the label, then
    if (n >= gcan->nlabels) {
      // if (gcan->total_training > 0)
      // return(log(0.01f/((float)gcan->total_training*GIBBS_NEIGHBORS))) ;
      /* 10*GIBBS_NEIGHBORS*BIG_AND_NEGATIVE*/
      // else
      return (10 * BIG_AND_NEGATIVE);
      // return(log(VERY_UNLIKELY)) ;
    }

    gc = &gcan->gcs[n];

    /* compute 1-d Mahalanobis distance */
    log_posterior = GCAcomputeConditionalLogDensity(gc, vals, gca->ninputs, gcan->labels[n]);
    if (check_finite("GCAvoxelGibbsLogPosterior: conditional log density", log_posterior) == 0) {
      DiagBreak();
    }

// added to the previous value
#if INTERP_PRIOR
    prior = gcaComputePrior(gca, mri_inputs, transform, x, y, z, label);
    log_posterior += (log(prior));
#else
    log_posterior += (log(getPrior(gcap, label)));
#endif
    if (check_finite("GCAvoxelGibbsLogPosterior: final", log_posterior) == 0) {
      DiagBreak();
    }
  }
  else {
    return (10 * BIG_AND_NEGATIVE);
    // return (log(VERY_UNLIKELY)) ;
  }

  return (log_posterior);
}

MRI *GCAbuildMostLikelyVolume(GCA *gca, MRI *mri)
{
  int x, y, z, xn, yn, zn, width, depth, height, n, xp, yp, zp, r;
  const GCA_NODE *gcan;
  const GCA_PRIOR *gcap;
  double max_prior;
  int max_label;
  GC1D *gc_max;

  if (!mri) {
    mri = MRIallocSequence(gca->prior_width, gca->prior_height, gca->prior_depth, MRI_FLOAT, gca->ninputs);
    // hey create gca volume and thus copies gca prior values
    mri->xsize = gca->xsize * gca->prior_spacing;
    mri->ysize = gca->ysize * gca->prior_spacing;
    mri->zsize = gca->zsize * gca->prior_spacing;
  }
  // most likely volume should agree with direction cosines
  GCAcopyDCToMRI(gca, mri);

  if (mri->nframes != gca->ninputs)
    ErrorExit(ERROR_BADPARM,
              "GCAbuildMostLikelyVolume: mri->frames "
              "(%d) does not match gca->ninputs (%d)",
              mri->nframes,
              gca->ninputs);

  // mri is prior if mri = NULL
  width = mri->width;
  depth = mri->depth;
  height = mri->height;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        // get node value
        if (GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn) == NO_ERROR) {
          // get prior value
          if (GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp) == NO_ERROR) {
            gcan = &gca->nodes[xn][yn][zn];
            gcap = &gca->priors[xp][yp][zp];
            if (gcap == NULL || gcap->nlabels <= 0) {
              continue;
            }
            // initialize
            max_prior = gcap->priors[0];
            max_label = gcap->labels[0];
            gc_max = NULL;
            // prior labels
            for (n = 1; n < gcap->nlabels; n++) {
              if (gcap->priors[n] >= max_prior) {
                max_prior = gcap->priors[n];
                max_label = gcap->labels[n];
              }
            }
            // get max_prior, max_label
            // go through node labels
            for (n = 0; n < gcan->nlabels; n++) {
              if (gcan->labels[n] == max_label) {
                gc_max = &gcan->gcs[n];
              }
            }

            if (!gc_max) {
              continue;
            }
            for (r = 0; r < gca->ninputs; r++) {
              MRIsetVoxVal(mri, x, y, z, r, gc_max->means[r]);
            }
          }
          else {
            for (r = 0; r < gca->ninputs; r++) {
              MRIsetVoxVal(mri, x, y, z, r, 0);
            }
          }
        }
        else {
          for (r = 0; r < gca->ninputs; r++) {
            MRIsetVoxVal(mri, x, y, z, r, 0);
          }
        }
      }
    }
  }

  return (mri);
}

MRI *GCAbuildMostLikelyVolumeFrame(GCA *gca, MRI *mri, int frame)
{
  int x, y, z, xn, yn, zn, width, depth, height, n, xp, yp, zp;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  double max_prior;
  int max_label;
  GC1D *gc_max;

  if (!mri) {
    mri = MRIallocSequence(gca->prior_width, gca->prior_height, gca->prior_depth, MRI_FLOAT, 1);
    mri->xsize = gca->xsize * gca->prior_spacing;
    mri->ysize = gca->ysize * gca->prior_spacing;
    mri->zsize = gca->zsize * gca->prior_spacing;
  }
  // gca volume direction cosines must be copied
  GCAcopyDCToMRI(gca, mri);

  width = mri->width;
  depth = mri->depth;
  height = mri->height;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if (!GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn))
          if (!GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp)) {
            gcan = &gca->nodes[xn][yn][zn];
            gcap = &gca->priors[xp][yp][zp];
            if (gcap == NULL || gcap->nlabels <= 0) {
              continue;
            }
            if (gcap->nlabels == 0) {
              continue;
            }
            max_prior = gcap->priors[0];
            max_label = gcap->labels[0];
            gc_max = NULL;
            for (n = 1; n < gcap->nlabels; n++) {
              if (gcap->priors[n] >= max_prior) {
                max_prior = gcap->priors[n];
                max_label = gcap->labels[n];
              }
            }
            for (n = 0; n < gcan->nlabels; n++) {
              if (gcan->labels[n] == max_label) {
                gc_max = &gcan->gcs[n];
              }
            }

            if (!gc_max) {
              continue;
            }
            MRIsetVoxVal(mri, x, y, z, 0, gc_max->means[frame]);
          }
      }
    }
  }

  return (mri);
}

GC1D *GCAfindPriorGC(const GCA *gca, int xp, int yp, int zp, int label)
{
  int xn, yn, zn;

  if (!GCApriorToNode(gca, xp, yp, zp, &xn, &yn, &zn)) {
    return (GCAfindGC(gca, xn, yn, zn, label));
  }
  else {
    return NULL;
  }
}

GC1D *GCAfindGC(const GCA *gca, int xn, int yn, int zn, int label)
{
  int n;
  GCA_NODE *gcan;

  gcan = &gca->nodes[xn][yn][zn];

  for (n = 0; n < gcan->nlabels; n++) {
    if (gcan->labels[n] == label) {
      return (&gcan->gcs[n]);
    }
  }

  return (NULL);
}

#include "mrisegment.h"
/* each segment must be at least this much of total to be retained */
#define MIN_SEG_PCT 0.15

static int gcaReclassifySegment(GCA *gca,
                                MRI *mri_inputs,
                                MRI *mri_labels,
                                MRI_SEGMENT *mseg,
                                int old_label,
                                TRANSFORM *transform,
                                int *exclude_labels,
                                int nexcluded);
static int gcaReclassifyVoxel(GCA *gca,
                              MRI *mri_inputs,
                              MRI *mri_labels,
                              int x,
                              int y,
                              int z,
                              int old_label,
                              TRANSFORM *transform,
                              int *exclude_list,
                              int nexcluded);
static int choroid_labels[] = {Left_choroid_plexus, Right_choroid_plexus};
#define NCHOROID (sizeof(choroid_labels) / sizeof(choroid_labels[0]))

MRI *GCAconstrainLabelTopology(GCA *gca,
                               MRI *mri_inputs,
                               MRI *mri_src,
                               MRI *mri_dst,
                               TRANSFORM *transform,
                               double vent_topo_dist,
                               double vent_topo_volume_thresh1,
                               double vent_topo_volume_thresh2)

{
  int i, nex, *ex, j, nvox, ndilate; /*, x, y, z, width, height, depth*/
  MRI_SEGMENTATION *mriseg;
  MRI *mri_dilated, *mri_in_main_segment;
  double voxel_volume;

  voxel_volume = mri_inputs->xsize * mri_inputs->ysize * mri_inputs->zsize;
  ndilate = nint(vent_topo_dist / ((mri_inputs->xsize + mri_inputs->ysize + mri_inputs->zsize) / 3));
  mri_dst = MRIcopy(mri_src, mri_dst);

  for (i = 1; i <= MAX_CMA_LABEL; i++) {
    if (!IS_BRAIN(i))  // label is not brain, ignore
    {
      continue;
    }
    // count number of label i in dst
    nvox = MRIvoxelsInLabel(mri_dst, i);
    // no such label, continue
    if (!nvox) {
      continue;
    }
    // if hypointensities, continue
    if (LABEL_WITH_NO_TOPOLOGY_CONSTRAINT(i)) {
      continue;
    }

    /*
      for the ventricles we want to get rid of segments that are not close to the
      main/biggest one, but being wary of small disconnections. Dilate the labels
      a few times to close small disconnections, then erase segments that are still
      not connected.
    */
    if (IS_LAT_VENT(i)) {
      int max_segno;

      mri_dilated = MRIdilateLabel(mri_src, NULL, i, ndilate);
      mriseg = MRIsegment(mri_dilated, (float)i, (float)i);
      max_segno = MRIfindMaxSegmentNumber(mriseg);
      mri_in_main_segment = MRIsegmentFill(mriseg, max_segno, NULL, 1);  // only retain biggest segment
      MRIfree(&mri_dilated);
      MRIsegmentFree(&mriseg);
      nex = NCHOROID;
      ex = choroid_labels;
    }
    else {
      mri_in_main_segment = NULL;
      nex = 0;
      ex = NULL;
    }

    /*    printf("label %03d: %d voxels\n", i, nvox) ;*/
    mriseg = MRIsegment(mri_src, (float)i, (float)i);
    if (!mriseg) {
      ErrorPrintf(Gerror,
                  "GCAconstrainLabelTopology: "
                  "label %s failed (%d vox)",
                  cma_label_to_name(i),
                  nvox);
      continue;
    }

    /*    printf("\t%d segments:\n", mriseg->nsegments) ;*/
    for (j = 0; j < mriseg->nsegments; j++) {
      // every voxel in the lateral ventricle label should be in the max segment found above, or just big enough and
      // keep it
      if (IS_LAT_VENT(i) && ((mriseg->segments[j].nvoxels * voxel_volume > vent_topo_volume_thresh2) ||
                             ((mriseg->segments[j].nvoxels * voxel_volume > vent_topo_volume_thresh1) &&
                              MRIgetVoxVal(mri_in_main_segment,
                                           mriseg->segments[j].voxels[0].x,
                                           mriseg->segments[j].voxels[0].y,
                                           mriseg->segments[j].voxels[0].z,
                                           0) > 0))) {
        if (mriseg->segments[j].nvoxels * voxel_volume > vent_topo_volume_thresh2)
          printf(" !!!!!!!!! ventricle segment %d with volume %2.0f above threshold %2.0f - not erasing !!!!!!!!!!\n",
                 j,
                 mriseg->segments[j].nvoxels * voxel_volume,
                 vent_topo_volume_thresh2);
        continue;
      }
      /* printf("\t\t%02d: %d voxels", j, mriseg->segments[j].nvoxels) ;*/
      if ((float)mriseg->segments[j].nvoxels / (float)nvox < MIN_SEG_PCT) {
        // printf(" - reclassifying...") ;
        gcaReclassifySegment(gca, mri_inputs, mri_dst, &mriseg->segments[j], i, transform, ex, nex);
      }
      // printf("\n") ;
    }
    MRIsegmentFree(&mriseg);
    if (mri_in_main_segment) {
      MRIfree(&mri_in_main_segment);
    }
  }


  return (mri_dst);
}

static int gcaReclassifySegment(GCA *gca,
                                MRI *mri_inputs,
                                MRI *mri_labels,
                                MRI_SEGMENT *mseg,
                                int old_label,
                                TRANSFORM *transform,
                                int *exclude_list,
                                int nexcluded)
{
  int i;

  for (i = 0; i < mseg->nvoxels; i++) {
    gcaReclassifyVoxel(gca,
                       mri_inputs,
                       mri_labels,
                       mseg->voxels[i].x,
                       mseg->voxels[i].y,
                       mseg->voxels[i].z,
                       old_label,
                       transform,
                       exclude_list,
                       nexcluded);
  }

  return (NO_ERROR);
}

static int gcaReclassifyVoxel(GCA *gca,
                              MRI *mri_inputs,
                              MRI *mri_labels,
                              int x,
                              int y,
                              int z,
                              int old_label,
                              TRANSFORM *transform,
                              int *exclude_list,
                              int nexcluded)
{
  int nbr_labels[MAX_CMA_LABELS], xi, yi, zi, xk, yk, zk, i, new_label, ex, skip;
  double max_p, p;

  memset(nbr_labels, 0, sizeof(int) * MAX_CMA_LABELS);
  // get 6 neighbors
  for (zk = -1; zk <= 1; zk++) {
    zi = mri_labels->zi[z + zk];
    for (yk = -1; yk <= 1; yk++) {
      yi = mri_labels->yi[y + yk];
      for (xk = -1; xk <= 1; xk++) {
        xi = mri_labels->xi[x + xk];
        // get the label histogram
        nbr_labels[nint(MRIgetVoxVal(mri_labels, xi, yi, zi, 0))]++;
      }
    }
  }
  new_label = 0;
  max_p = 10 * BIG_AND_NEGATIVE;
  nbr_labels[old_label] = 0;  // don't look at old_label

  // for (i = 0 ; i <= 255 ; i++) // should not go up to 255 but MAX_CMA_LABEL
  for (i = 0; i < MAX_CMA_LABEL; i++) {
    if (mri_labels->type == MRI_UCHAR && i >= 256) {
      break;
    }
    for (skip = ex = 0; ex < nexcluded; ex++)
      if (exclude_list[ex] == i) {
        skip = 1;
        break;
      }
    if (skip) {
      continue;
    }

    if (nbr_labels[i] > 0)  // if neighbors has this label, then
    {
      MRIsetVoxVal(mri_labels, x, y, z, 0, i);
      // set to the current label and see what happens
      p = GCAvoxelGibbsLogPosterior(gca, mri_labels, mri_inputs, x, y, z, transform, PRIOR_FACTOR);
      // debug ///////////////////////////////////////
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
        printf("gcaReclassifyVoxel: nbr_labels[%d] = %d\n", i, nbr_labels[i]);
        printf(
            "gcaReclassifyVoxel:(%d, %d, %d): Label = %s (%d), "
            "VoxelGibbsLogPosterior = %.2f, max_p = %.2f\n",
            x,
            y,
            z,
            cma_label_to_name(i),
            i,
            p,
            max_p);
        if (p >= max_p)
          printf(
              "   replacing max_p with this p and "
              "label from %s(%d) to %s(%d)\n",
              cma_label_to_name(old_label),
              old_label,
              cma_label_to_name(i),
              i);
      }
      ////////////////////////////////////////////////
      if (p >= max_p) {
        max_p = p;
        new_label = i;
      }
    }
  }
  MRIsetVoxVal(mri_labels, x, y, z, 0, new_label);
  return (NO_ERROR);
}
#define MAX_VENTRICLE_ITERATIONS 30
#define V_WSIZE (5)
#define V_WHALF (V_WSIZE - 1) / 2
#define V_THRESH ((V_WSIZE * V_WSIZE / 2))
MRI *GCAexpandVentricle(GCA *gca, MRI *mri_inputs, MRI *mri_src, MRI *mri_dst, TRANSFORM *transform, int target_label)
{
  int nchanged, x, y, z, xn, yn, zn, xk, yk, zk, xi, yi, zi, label, total_changed, i, j, count, found, xmin, xmax, ymin,
      ymax, zmin, zmax;
  // int width, height, depth;
  // GCA_NODE *gcan;
  GC1D *gc;
  float v_means[MAX_GCA_INPUTS], v_var, dist_ven, dist_label, vals[MAX_GCA_INPUTS];

  /* compute label mean and variance */

  GCAcomputeLabelStats(gca, target_label, &v_var, v_means);
  printf("ventricle intensity = ");
  for (i = 0; i < gca->ninputs; i++) {
    printf("%2.1f ", v_means[i]);
  }
  printf(" +- %2.1f\n", sqrt(v_var));

  if (mri_src != mri_dst) {
    mri_dst = MRIcopy(mri_src, mri_dst);
  }

  // width = mri_src->width;
  // height = mri_src->height;
  // depth = mri_src->depth;

  i = total_changed = 0;
  xmin = mri_dst->width;
  ymin = mri_dst->height;
  zmin = mri_dst->depth;
  xmax = ymax = zmax = -1;
  for (z = 0; z < mri_dst->depth; z++) {
    for (y = 0; y < mri_dst->height; y++) {
      for (x = 0; x < mri_dst->width; x++) {
        label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0));
        if (label == target_label) {
          if (x > xmax) {
            xmax = x;
          }
          if (y > ymax) {
            ymax = y;
          }
          if (z > zmax) {
            zmax = z;
          }
          if (x < xmin) {
            xmin = x;
          }
          if (y < ymin) {
            ymin = y;
          }
          if (z < zmin) {
            zmin = z;
          }
        }
      }
    }
  }
  // expand bounding box to allow for ventricle to expand
  xmin = mri_dst->xi[xmin - 1];
  ymin = mri_dst->yi[ymin - 1];
  zmin = mri_dst->zi[zmin - 1];
  xmax = mri_dst->xi[xmax + 1];
  ymax = mri_dst->yi[ymax + 1];
  zmax = mri_dst->zi[zmax + 1];

  do {
    nchanged = 0;
    for (z = zmin; z <= zmax; z++) {
      for (y = ymin; y <= ymax; y++) {
        for (x = xmin; x <= xmax; x++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0));
          if (label == target_label) {
            continue;
          }
          found = 0;

          for (yk = -1; found == 0 && yk <= 1; yk += 2) {
            yi = mri_dst->yi[y + yk];  // look superior/inferior

            // only change it if there is a "wall" of
            // ventricle superior or inferior
            count = MRIlabelsInPlanarNbhd(mri_dst, x, yi, z, V_WHALF, target_label, MRI_HORIZONTAL);
            if (count >= V_THRESH) {
              found = 1;
            }
          }
          for (xk = -1; found == 0 && xk <= 1; xk += 2) {
            xi = mri_dst->xi[x + xk];  // look superior/inferior

            // only change it if there is a "wall" of
            // ventricle superior or inferior
            count = MRIlabelsInPlanarNbhd(mri_dst, xi, y, z, V_WHALF, target_label, MRI_SAGITTAL);
            if (count >= V_THRESH) {
              found = 1;
            }
          }
          for (zk = -1; found == 0 && zk <= 1; zk += 2) {
            zi = mri_dst->zi[z + zk];  // look superior/inferior

            // only change it if there is a "wall" of
            // ventricle anterior or posterior
            count = MRIlabelsInPlanarNbhd(mri_dst, x, y, zi, V_WHALF, target_label, MRI_CORONAL);
            if (count >= V_THRESH) {
              found = 1;
            }
          }
          if (found == 0) {
            continue;
          }
          if (GCAsourceVoxelToNode(gca, mri_dst, transform, x, y, z, &xn, &yn, &zn) == NO_ERROR) {
            // gcan = &gca->nodes[xn][yn][zn];

            if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
              DiagBreak();
            }

            load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
            for (dist_ven = 0, j = 0; j < gca->ninputs; j++) {
              dist_ven += SQR(vals[j] - v_means[j]);
            }
            gc = GCAfindGC(gca, xn, yn, zn, label);
            if (gc == NULL) {
              dist_label = 10000;  // ???
            }
            else {
              dist_label = GCAmahDistIdentityCovariance(gc, vals, gca->ninputs);
            }
            gc = GCAfindGC(gca, xn, yn, zn, target_label);
            if (2 * dist_ven < dist_label)  // much more like ventricle
                                            // than anything else
            {
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
                int olabel = nint(MRIgetVoxVal(mri_dst, x, y, z, 0));
                printf(
                    "GCAexpandVentricle:voxel"
                    "(%d, %d, %d) changed from %s (%d) "
                    "to %s (%d), because current label d = %2.0f "
                    "and new lable d = %2.0f\n",
                    x,
                    y,
                    z,
                    cma_label_to_name(olabel),
                    olabel,
                    cma_label_to_name(target_label),
                    target_label,
                    dist_label,
                    dist_ven);
              }
              // change it to ventricle
              nchanged++;
              MRIsetVoxVal(mri_dst, x, y, z, 0, target_label);
              if (x <= xmin) {
                xmin = mri_dst->xi[x - 1];
              }
              if (y <= ymin) {
                ymin = mri_dst->yi[y - 1];
              }
              if (z <= zmin) {
                zmin = mri_dst->zi[z - 1];
              }
              if (x >= xmax) {
                xmax = mri_dst->xi[x + 1];
              }
              if (y >= ymax) {
                ymax = mri_dst->yi[y + 1];
              }
              if (z >= zmax) {
                zmax = mri_dst->zi[z + 1];
              }
            }
          }
        }
      }
    }

    total_changed += nchanged;
    if (++i > MAX_VENTRICLE_ITERATIONS) {
      break;
    }
  } while (nchanged > 0);

  printf("%d labels changed to %s\n", total_changed, cma_label_to_name(target_label));
  return (mri_dst);
}
#define MAX_CORTICAL_ITERATIONS 10
MRI *GCAexpandCortex(GCA *gca, MRI *mri_inputs, MRI *mri_src, MRI *mri_dst, TRANSFORM *transform)
{
  int nchanged, x, y, z, width, height, depth, xn, yn, zn, xi, yi, zi, xk, yk, zk, label, total_changed, i, wm_nbr,
      gray_nbr, left;
  // GCA_NODE *gcan;
  float wm_means[MAX_GCA_INPUTS], wm_var, gray_means[MAX_GCA_INPUTS], gray_var, ldist, wdist, gdist;
  MRI *mri_tmp;
  float vals[MAX_GCA_INPUTS];
  GC1D *gc_wm, *gc_gm, *gc_label;

  /* compute label mean and variance */

  GCAcomputeLabelStats(gca, Left_Cerebral_White_Matter, &wm_var, wm_means);
  GCAcomputeLabelStats(gca, Left_Cerebral_Cortex, &gray_var, gray_means);
  printf("cortex mean - gray ");
  for (i = 0; i < gca->ninputs; i++) {
    printf("%2.1f ", gray_means[i]);
  }
  printf("+- %2.1f, white ", sqrt(gray_var));
  for (i = 0; i < gca->ninputs; i++) {
    printf("%2.1f ", wm_means[i]);
  }
  printf("+- %2.1f\n", sqrt(wm_var));

  if (mri_src != mri_dst) {
    mri_dst = MRIcopy(mri_src, mri_dst);
  }

  mri_tmp = MRIcopy(mri_dst, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  i = total_changed = 0;
  do {
    nchanged = 0;
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }

          label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0));

          if (label == Unknown || label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex ||
              label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter) {
            if (label != Unknown)
              left = (label == Left_Cerebral_Cortex || label == Left_Cerebral_White_Matter);
            else {
              left = -1;
              for (zk = -1; left < 0 && zk <= 1; zk++) {
                for (yk = -1; left < 0 && yk <= 1; yk++) {
                  for (xk = -1; left < 0 && xk <= 1; xk++) {
                    xi = mri_src->xi[x + xk];
                    yi = mri_src->yi[y + yk];
                    zi = mri_src->zi[z + zk];
                    label = nint(MRIgetVoxVal(mri_dst, xi, yi, zi, 0));
                    if (label == Left_Cerebral_Cortex || label == Left_Cerebral_White_Matter) {
                      left = 1;
                    }
                    else if (label == Right_Cerebral_Cortex || label == Right_Cerebral_White_Matter) {
                      left = 0;
                    }
                  }
                }
              }
            }

            gray_nbr = wm_nbr = 0;
            gc_wm = GCAfindSourceGC(
                gca, mri_src, transform, x, y, z, left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter);
            gc_gm =
                GCAfindSourceGC(gca, mri_src, transform, x, y, z, left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex);
            if (gc_gm) gray_nbr = left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex;
            if (gc_wm) wm_nbr = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
            if (gc_gm == NULL || gc_wm == NULL)
              for (zk = -1; zk <= 1; zk++) {
                for (yk = -1; yk <= 1; yk++) {
                  for (xk = -1; xk <= 1; xk++) {
                    if (fabs(xk) + fabs(yk) + fabs(zk) > 1) {
                      continue;
                    }
                    xi = mri_src->xi[x + xk];
                    yi = mri_src->yi[y + yk];
                    zi = mri_src->zi[z + zk];
                    label = nint(MRIgetVoxVal(mri_dst, xi, yi, zi, 0));
                    if ((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) && (gc_gm == NULL)) {
                      gray_nbr = label;
                      gc_gm = GCAfindSourceGC(gca, mri_src, transform, xi, yi, zi, label);
                      if (!gc_gm) {
                        gray_nbr = 0;
                      }
                      /* shouldn't ever happen */
                    }
                    else if ((label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter) &&
                             (gc_wm == NULL)) {
                      wm_nbr = label;
                      gc_wm = GCAfindSourceGC(gca, mri_src, transform, xi, yi, zi, label);
                      if (!gc_wm) {
                        wm_nbr = 0;
                      }
                      /* shouldn't ever happen */
                    }
                  }
                }
              }
            if (!wm_nbr && !gray_nbr) {
              continue;
            }

            load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
            if (wm_nbr)
              wdist = sqrt(GCAmahDistIdentityCovariance(gc_wm, vals, gca->ninputs));
            else {
              wdist = 1e10;
            }
            if (gray_nbr)
              gdist = sqrt(GCAmahDistIdentityCovariance(gc_gm, vals, gca->ninputs));
            else {
              gdist = 1e10;
            }

            if (gc_wm == NULL || sqrt(GCAmahDist(gc_wm, vals, gca->ninputs)) > 1.5) {
              wdist = 1e10; /* hack - don't label unlikely white */
            }
            if (gc_gm == NULL || sqrt(GCAmahDist(gc_gm, vals, gca->ninputs)) > 1.5) {
              gdist = 1e10; /* hack - don't label unlikely gray */
            }

            if (!GCAsourceVoxelToNode(gca, mri_dst, transform, x, y, z, &xn, &yn, &zn)) {
              // gcan = &gca->nodes[xn][yn][zn];
              gc_label = GCAfindGC(gca, xn, yn, zn, label);
              if (gc_label)
                ldist = sqrt(GCAmahDistIdentityCovariance(gc_label, vals, gca->ninputs));
              else {
                ldist = 1e10;
              }
              ldist *= .75; /* bias towards retaining label */

              if ((wdist < gdist) && gc_wm)
              /* might change to wm */
              {
                if (wdist < ldist && GCAisPossible(gca, mri_inputs, wm_nbr, transform, x, y, z, 0)) {
                  if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
                    int olabel = nint(MRIgetVoxVal(mri_tmp, x, y, z, 0));
                    if (olabel != wm_nbr)
                      printf(
                          "GCAexpandCortex:voxel (%d, %d, %d)"
                          " changed from %s (%d) "
                          "to %s (%d), wdist=%.2f, gdist=%.2f,"
                          " ldist=%.2f\n",
                          x,
                          y,
                          z,
                          cma_label_to_name(olabel),
                          olabel,
                          cma_label_to_name(wm_nbr),
                          wm_nbr,
                          wdist,
                          gdist,
                          ldist);
                  }
                  nchanged++;
                  MRIsetVoxVal(mri_tmp, x, y, z, 0, wm_nbr);
                }
              }
              else if (gc_gm) /* might change to gm */
              {
                if (gdist < ldist && GCAisPossible(gca, mri_inputs, gray_nbr, transform, x, y, z, 0)) {
                  int olabel = nint(MRIgetVoxVal(mri_tmp, x, y, z, 0));
                  if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
                    if (olabel != gray_nbr)
                      printf(
                          "GCAexpandCortex:voxel "
                          "(%d, %d, %d) changed from %s (%d) "
                          "to %s (%d), wdist=%.2f, "
                          "gdist=%.2f, ldist=%.2f\n",
                          x,
                          y,
                          z,
                          cma_label_to_name(olabel),
                          olabel,
                          cma_label_to_name(gray_nbr),
                          gray_nbr,
                          wdist,
                          gdist,
                          ldist);
                  }
                  if (olabel != gray_nbr) {
                    nchanged++;
                  }

                  MRIsetVoxVal(mri_tmp, x, y, z, 0, gray_nbr);
                }
              }
            }
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst);
    total_changed += nchanged;
    if (++i > MAX_CORTICAL_ITERATIONS) {
      break;
    }
  } while (nchanged > 0);

  MRIfree(&mri_tmp);
  printf("%d labels changed to cortex...\n", total_changed);
  return (mri_dst);
}

MRI *GCAexpandLabelIntoWM(
    GCA *gca, MRI *mri_inputs, MRI *mri_src, MRI *mri_dst, TRANSFORM *transform, MRI *mri_fixed, int target_label)
{
  int nchanged, x, y, z, width, height, depth, xn, yn, zn, xi, yi, zi, xk, yk, zk, nbr_label, n, label, total_changed,
      i;
  GCA_NODE *gcan, *gcan_nbr;
  GC1D *gc_label, *gc_wm;
  MRI *mri_tmp;
  float vals[MAX_GCA_INPUTS];
  double prior;

  if (mri_src != mri_dst) {
    mri_dst = MRIcopy(mri_src, mri_dst);
  }

  mri_tmp = MRIcopy(mri_dst, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  i = total_changed = 0;
  do {
    nchanged = 0;
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          // get the label
          label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0));

          if (!GCAsourceVoxelToNode(gca, mri_dst, transform, x, y, z, &xn, &yn, &zn)) {
            gcan = &gca->nodes[xn][yn][zn];

            // if this label is the same as the target label,
            // then expand
            // into neighbors if conditions are satisfied.
            if (label == target_label) {
              gc_label = GCAfindGC(gca, xn, yn, zn, label);
              // find wm gaussian classifier
              gc_wm = NULL;
              for (n = 0; n < gcan->nlabels; n++) {
                if ((gcan->labels[n] == Left_Cerebral_White_Matter) || (gcan->labels[n] == Right_Cerebral_White_Matter))
                  gc_wm = GCAfindGC(gca, xn, yn, zn, gcan->labels[n]);
              }
              // look around the neighbors
              for (zk = -1; zk <= 1; zk++) {
                for (yk = -1; yk <= 1; yk++) {
                  for (xk = -1; xk <= 1; xk++) {
                    if (fabs(xk) + fabs(yk) + fabs(zk) > 1) {
                      continue;
                    }
                    xi = mri_src->xi[x + xk];
                    yi = mri_src->yi[y + yk];
                    zi = mri_src->zi[z + zk];
                    if (xi == Ggca_x && yi == Ggca_y && zi == Ggca_z) {
                      DiagBreak();
                    }
                    // get the neighbor label
                    nbr_label = nint(MRIgetVoxVal(mri_dst, xi, yi, zi, 0));
                    if ((nbr_label == Right_Cerebral_White_Matter) || (nbr_label == Left_Cerebral_White_Matter)) {
                      // if it is wm, then load
                      // grey values at this neighbor
                      load_vals(mri_inputs, xi, yi, zi, vals, gca->ninputs);
                      if (!GCAsourceVoxelToNode(gca, mri_dst, transform, xi, yi, zi, &xn, &yn, &zn)) {
                        // get the wm gc
                        gc_wm = GCAfindGC(gca, xn, yn, zn, nbr_label);
                        // get the target label gc
                        gc_label = GCAfindGC(gca, xn, yn, zn, target_label);
                        if (!gc_wm || !gc_label) {
                          continue;
                        }
                        // calculate distance
                        // to wm and target label.
                        // if wm is bigger
                        if (GCAmahDistIdentityCovariance(gc_wm, vals, gca->ninputs) >
                            GCAmahDistIdentityCovariance(gc_label, vals, gca->ninputs)) {
                          gcan_nbr = &gca->nodes[xn][yn][zn];
                          for (prior = 0.0f, n = 0; n < gcan_nbr->nlabels; n++) {
                            // look for the
                            // target label
                            // in this neighbor
                            if (gcan_nbr->labels[n] == target_label) {
                              prior = get_node_prior(gca, target_label, xn, yn, zn);
                              if (prior != 0) {
                                break;
                              }
                            }
                          }
// if not found,
// prior stays 0.0f
#define PRIOR_THRESH 0.01
                          if (prior >= PRIOR_THRESH)
                          /* target is possible */
                          {
                            if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
                              printf(
                                  "GCAexpandLabelIntoWM:voxel (%d, %d, %d) "
                                  "changed from %s (%d) "
                                  "to %s (%d), prior=%.2f\n",
                                  xi,
                                  yi,
                                  zi,
                                  cma_label_to_name(nbr_label),
                                  nbr_label,
                                  cma_label_to_name(target_label),
                                  target_label,
                                  prior);
                            }
                            nchanged++;
                            MRIsetVoxVal(mri_tmp, xi, yi, zi, 0, target_label);
                            /* MRIvox(mri_fixed,
                               xi, yi, zi) = 0 ;*/
                          }
                        }
                      }
                      ///////////////////
                    }
                  }
                }
              }
            }
          }
          /////////////////////////////////////////////
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst);
    total_changed += nchanged;
    if (++i >= 1) {
      break;
    }
  } while (nchanged > 0);

  MRIfree(&mri_tmp);
  printf("%d labels changed to %s\n", total_changed, cma_label_to_name(target_label));
  return (mri_dst);
}

int GCArankSamples(GCA *gca, GCA_SAMPLE *gcas, int nsamples, int *ordered_indices)
{
  LABEL_PROB *label_probs;
  int i;

  label_probs = (LABEL_PROB *)calloc(nsamples, sizeof(LABEL_PROB));
  for (i = 0; i < nsamples; i++) {
    label_probs[i].label = i;
    label_probs[i].prob = gcas[i].log_p;
  }

  /* now sort the samples by probability */
  qsort(label_probs, nsamples, sizeof(LABEL_PROB), compare_sort_probabilities);

  for (i = 0; i < nsamples; i++) {
    ordered_indices[i] = label_probs[nsamples - (i + 1)].label;
  }

  free(label_probs);
  return (NO_ERROR);
}

#include "mrinorm.h"
MRI *GCAnormalizeSamples(
    MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas, int nsamples, TRANSFORM *transform, const char *ctl_point_fname)
{
  MRI *mri_dst, *mri_ctrl, *mri_bias;
  int xv, yv, zv, n, x, y, z, width, height, depth, xn, yn, zn, num, total, input, T1_index = 0;
  float bias;
  double mean, sigma;
  double val;
  float gm_means[MAX_GCA_INPUTS], gray_white_CNR;

  if (nsamples == 0) /* only using control points from file */
  {
    float wm_means[MAX_GCA_INPUTS], tmp[MAX_GCA_INPUTS];
    // float max_wm;
    int r;

    GCAlabelMean(gca, Left_Cerebral_White_Matter, wm_means);
    //    GCAlabelMean(gca, Left_Cerebral_White_Matter, tmp) ;
    GCAlabelMean(gca, Right_Cerebral_White_Matter, tmp);

    GCAlabelMean(gca, Left_Cerebral_Cortex, gm_means);
    gray_white_CNR = wm_means[0] - gm_means[0];
    for (r = 0; r < gca->ninputs; r++) {
      wm_means[r] = (wm_means[r] + tmp[r]) / 2;
      if ((wm_means[r] - gm_means[r]) > gray_white_CNR) {
        T1_index = r;
        gray_white_CNR = (wm_means[r] - gm_means[r]);
      }
    }
    // max_wm = wm_means[T1_index];
    printf("using volume %d as most T1-weighted for normalization\n", T1_index);
  }
  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  mri_dst = MRIclone(mri_in, NULL);
  mri_ctrl = MRIalloc(width, height, depth, MRI_UCHAR);
  MRIcopyHeader(mri_in, mri_ctrl);
  mri_bias = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_SHORT);
  if (!mri_bias)
    ErrorExit(ERROR_NOMEMORY,
              "GCAnormalizeSamples: could not allocate "
              "(%d,%d,%d,2) bias image",
              mri_in->width,
              mri_in->height,
              mri_in->depth);
  MRIcopyHeader(mri_in, mri_bias);

#define MAX_BIAS 1250
#define NO_BIAS 1000
#define MIN_BIAS 750

  if (ctl_point_fname) {
    MRI3dUseFileControlPoints(mri_ctrl, ctl_point_fname);
    MRInormAddFileControlPoints(mri_ctrl, CONTROL_MARKED, NULL);
  }

  /* add control points from file */
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        MRISvox(mri_bias, x, y, z) = NO_BIAS; /* by default */
        if (MRIvox(mri_ctrl, x, y, z) != CONTROL_MARKED)
        /* not read from file */
        {
          continue;
        }

        if (nsamples == 0) /* only using file control points */
        {
          MRIsampleVolumeFrame(mri_in, x, y, z, T1_index, &val);
          bias = NO_BIAS * DEFAULT_DESIRED_WHITE_MATTER_VALUE / val;
          MRISvox(mri_bias, x, y, z) = (short)nint(bias);
        }
        else /* find atlas point this maps to */
        {
          int n, max_n;
          GC1D *gc;
          GCA_NODE *gcan;
          GCA_PRIOR *gcap;
          double max_p;

          if (!GCAsourceVoxelToNode(gca, mri_dst, transform, x, y, z, &xn, &yn, &zn)) {
            gcan = &gca->nodes[xn][yn][zn];
            gcap = getGCAP(gca, mri_dst, transform, x, y, z);
            if (gcap == NULL) {
              continue;
            }
            max_p = 0;
            for (max_n = -1, n = 0; n < gcan->nlabels; n++) {
              if ((0 == IS_WM(gcan->labels[n])) && (0 == IS_CEREBELLAR_WM(gcan->labels[n])) &&
                  (gcan->labels[n] != Brain_Stem)) {
                continue;
              }
              gc = &gcan->gcs[n];
              if (getPrior(gcap, gcan->labels[n]) >= max_p) {
                max_p = getPrior(gcap, gcan->labels[n]);
                max_n = n;
              }
            }
            if (max_n < 0)
            /* couldn't find any valid label at this location */
            {
              continue;
            }
            gc = &gcan->gcs[max_n];

            for (bias = 0.0, input = 0; input < gca->ninputs; input++) {
              MRIsampleVolumeFrame(mri_in, x, y, z, input, &val);
              if (FZERO(val)) {
                val = 1;
              }
              bias += (float)NO_BIAS * ((float)gc->means[input] / val);
            }
            bias /= (float)gca->ninputs;
            if (bias < 100 || bias > 5000) {
              DiagBreak();
            }
            if (bias < MIN_BIAS) {
              bias = MIN_BIAS;
            }
            if (bias > MAX_BIAS) {
              bias = MAX_BIAS;
            }

            MRISvox(mri_bias, x, y, z) = (short)nint(bias);
          }
          /////////////////////////////////////////////////////
        }
      }
    }
  }

  TransformInvert(transform, mri_in);
  for (n = 0; n < nsamples; n++) {
    if (gcas[n].xp == Gxp && gcas[n].yp == Gyp && gcas[n].zp == Gzp) {
      DiagBreak();
    }

    if (!GCApriorToSourceVoxel(gca, mri_dst, transform, gcas[n].xp, gcas[n].yp, gcas[n].zp, &xv, &yv, &zv)) {
      if (xv == 181 && yv == 146 && zv == 128) {
        DiagBreak();
      }
      if (xv == Ggca_x && yv == Ggca_y && zv == Ggca_z) {
        DiagBreak();
      }
      if (gcas[n].label == 29 || gcas[n].label == 61) {
        gcas[n].label = 0;
        DiagBreak();
      }
      if (gcas[n].label > 0) {
        MRIvox(mri_ctrl, xv, yv, zv) = CONTROL_MARKED;

        for (bias = 0.0, input = 0; input < gca->ninputs; input++) {
          MRIsampleVolumeFrame(mri_in, xv, yv, zv, input, &val);
          if (FZERO(val)) {
            val = 1;
          }
          bias += (float)NO_BIAS * ((float)gcas[n].means[input] / val);
        }
        bias /= (float)gca->ninputs;
        if (bias < 100 || bias > 5000) {
          DiagBreak();
        }

        MRISvox(mri_bias, xv, yv, zv) = (short)nint(bias);
      }
      else {
        MRIvox(mri_ctrl, xv, yv, zv) = CONTROL_NONE;
      }
    }
  }

  /* now check for and remove outliers */
  mean = sigma = 0.0;
  for (num = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED) {
          num++;
          bias = (double)MRISvox(mri_bias, x, y, z);
          mean += bias;
          sigma += (bias * bias);
        }
      }
    }
  }

  if (num > 0) {
    mean /= (double)num;
    sigma = sqrt(sigma / (double)num - mean * mean);
    printf("bias field = %2.3f +- %2.3f\n", mean / NO_BIAS, sigma / NO_BIAS);
  }

  /* now check for and remove outliers */
  for (total = num = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED) {
          bias = (double)MRISvox(mri_bias, x, y, z);
          total++;
          if (fabs(bias - mean) > 4 * sigma) {
            MRIvox(mri_ctrl, x, y, z) = CONTROL_NONE;
            num++;
            MRISvox(mri_bias, x, y, z) = NO_BIAS;
          }
        }
      }
    }
  }

  printf("%d of %d control points discarded\n", num, total);

  MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias);
/*  MRIwrite(mri_bias, "bias.mgz") ;*/
  {
    MRI *mri_kernel, *mri_smooth, *mri_down;
    float sigma = 16.0f;

    mri_down = MRIdownsample2(mri_bias, NULL);
    mri_kernel = MRIgaussian1d(sigma, 100);
    mri_smooth = MRIconvolveGaussian(mri_down, NULL, mri_kernel);
    MRIfree(&mri_bias);
    MRIfree(&mri_kernel);
    mri_bias = MRIupsample2(mri_smooth, NULL);
    sigma = 2.0f;
    MRIfree(&mri_down);
    MRIfree(&mri_smooth);
    mri_kernel = MRIgaussian1d(sigma, 100);
    mri_smooth = MRIconvolveGaussian(mri_bias, NULL, mri_kernel);
    MRIfree(&mri_bias);
    mri_bias = mri_smooth;
    MRIfree(&mri_kernel);
  }
  /*  MRIwrite(mri_bias, "smooth_bias.mgz") ;*/

  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        bias = (float)MRISvox(mri_bias, x, y, z) / NO_BIAS;
        if (bias < 0) {
          DiagBreak();
        }
        for (input = 0; input < gca->ninputs; input++) {
          MRIsampleVolumeFrame(mri_in, x, y, z, input, &val);
          val *= bias; /* corrected value */
          switch (mri_in->type) {
            case MRI_UCHAR:
              if (val < 0) {
                val = 0;
              }
              else if (val > 255) {
                val = 255;
              }
              MRIseq_vox(mri_dst, x, y, z, input) = (BUFTYPE)nint(val);
              break;
            case MRI_SHORT:
              MRISseq_vox(mri_dst, x, y, z, input) = (short)nint(val);
              break;
            case MRI_FLOAT:
              MRIFseq_vox(mri_dst, x, y, z, input) = val;
              break;
            default:
              ErrorReturn(NULL,
                          (ERROR_UNSUPPORTED,
                           "GCAnormalizeSamples: "
                           "unsupported input type %d",
                           mri_in->type));
              break;
          }
        }
      }
    }
  }

  MRIfree(&mri_bias);
  MRIfree(&mri_ctrl);
  return (mri_dst);
}

void GCAnormalizeSamplesOneChannel(MRI *mri_in,
                                   GCA *gca,
                                   GCA_SAMPLE *gcas,
                                   int nsamples,
                                   TRANSFORM *transform,
                                   char *ctl_point_fname,
                                   int input_index,
                                   double bias_sigma)
/* This function is added by xhan, trying to normalize a single channel */
{
  MRI *mri_dst, *mri_ctrl, *mri_bias;
  int xv, yv, zv, n, x, y, z, width, height, depth, xn, yn, zn, num, total, nctrl;
  float bias;
  double mean, sigma;
  double val;

  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  mri_dst = mri_in;
  mri_ctrl = MRIalloc(width, height, depth, MRI_UCHAR);
  MRIcopyHeader(mri_in, mri_ctrl);
  mri_bias = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_SHORT);
  if (!mri_bias)
    ErrorExit(ERROR_NOMEMORY,
              "GCAnormalizeSamples: could not allocate "
              "(%d,%d,%d,2) bias image",
              mri_in->width,
              mri_in->height,
              mri_in->depth);
  MRIcopyHeader(mri_in, mri_bias);

#define MAX_BIAS 1250
#define NO_BIAS 1000
#define MIN_BIAS 750

  if (ctl_point_fname) {
    MRI3dUseFileControlPoints(mri_ctrl, ctl_point_fname);
    nctrl = MRInormAddFileControlPoints(mri_ctrl, CONTROL_MARKED, NULL);
  }
  else {
    nctrl = 0;
  }

  /* add control points from file */
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();

        MRISvox(mri_bias, x, y, z) = NO_BIAS;            /* by default */
        if (MRIvox(mri_ctrl, x, y, z) != CONTROL_MARKED) /* not read from file */
          continue;

        if (nsamples == 0) /* only using file control points */
        {
          MRIsampleVolumeFrame(mri_in, x, y, z, input_index, &val);
          bias = NO_BIAS * DEFAULT_DESIRED_WHITE_MATTER_VALUE / val;
          MRISvox(mri_bias, x, y, z) = (short)nint(bias);
        }
        else /* find atlas point this maps to */
        {
          int n, max_n;
          GC1D *gc;
          GCA_NODE *gcan;
          GCA_PRIOR *gcap;
          double max_p;

          if (!GCAsourceVoxelToNode(gca, mri_dst, transform, x, y, z, &xn, &yn, &zn)) {
            gcan = &gca->nodes[xn][yn][zn];
            gcap = getGCAP(gca, mri_dst, transform, x, y, z);
            if (gcap == NULL) continue;

            max_p = 0;
            for (max_n = -1, n = 0; n < gcan->nlabels; n++) {
              if ((0 == IS_WM(gcan->labels[n])) && (0 == IS_CEREBELLAR_WM(gcan->labels[n])) &&
                  (gcan->labels[n] != Brain_Stem))
                continue;

              gc = &gcan->gcs[n];
              if (getPrior(gcap, gcan->labels[n]) >= max_p) {
                max_p = getPrior(gcap, gcan->labels[n]);
                max_n = n;
              }
            }
            if (max_n < 0) /* couldn't find any valid label at this location */
              continue;

            gc = &gcan->gcs[max_n];

            MRIsampleVolumeFrame(mri_in, x, y, z, input_index, &val);
            if (FZERO(val)) {
              val = 1;
            }
            bias = (float)NO_BIAS * ((float)gc->means[input_index] / val);

            if (bias < 100 || bias > 5000) DiagBreak();

            if (bias < MIN_BIAS) bias = MIN_BIAS;

            if (bias > MAX_BIAS) bias = MAX_BIAS;

            MRISvox(mri_bias, x, y, z) = (short)nint(bias);
          }
          /////////////////////////////////////////////////////
        }
      }
    }
  }

  TransformInvert(transform, mri_in);
  for (n = 0; n < nsamples; n++) {
    if (gcas[n].xp == Ggca_x && gcas[n].yp == Ggca_y && gcas[n].zp == Ggca_z) DiagBreak();

    xv = gcas[n].x;
    yv = gcas[n].y;
    zv = gcas[n].z;
    if (xv >= 0 && xv <= mri_dst->width - 1 && yv >= 0 && yv <= mri_dst->height - 1 && zv >= 0 &&
        zv <= mri_dst->depth - 1)
    {
      if (xv == 181 && yv == 146 && zv == 128) {
        DiagBreak();
      }
      if (xv == Ggca_x && yv == Ggca_y && zv == Ggca_z) {
        DiagBreak();
      }
      if (gcas[n].label == 29 || gcas[n].label == 61) {
        gcas[n].label = 0;
        DiagBreak();
      }
      if (gcas[n].label > 0) {
        MRIvox(mri_ctrl, xv, yv, zv) = CONTROL_MARKED;

        MRIsampleVolumeFrame(mri_in, xv, yv, zv, input_index, &val);
        if (FZERO(val)) {
          val = 1;
        }
        bias = (float)NO_BIAS * ((float)gcas[n].means[input_index] / val);

        if (bias < 100 || bias > 5000) {
          DiagBreak();
        }

        MRISvox(mri_bias, xv, yv, zv) = (short)nint(bias);
      }
      else {
        MRIvox(mri_ctrl, xv, yv, zv) = CONTROL_NONE;
      }
    }
  }

  /* now check for and remove outliers */
  mean = sigma = 0.0;
  for (num = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED) {
          num++;
          bias = (double)MRISvox(mri_bias, x, y, z);
          mean += bias;
          sigma += (bias * bias);
        }
      }
    }
  }

  if (num > 0) {
    mean /= (double)num;
    sigma = sqrt(sigma / (double)num - mean * mean);
    printf("bias field = %2.3f +- %2.3f\n", mean / NO_BIAS, sigma / NO_BIAS);
  }

  /* now check for and remove outliers */
  for (total = num = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED) {
          bias = (double)MRISvox(mri_bias, x, y, z);
          total++;
          if (fabs(bias - mean) > 4 * sigma) {
            MRIvox(mri_ctrl, x, y, z) = CONTROL_NONE;
            num++;
            MRISvox(mri_bias, x, y, z) = NO_BIAS;
          }
        }
      }
    }
  }

  printf("%d of %d control points discarded\n", num, total);

  MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias);
  /*  MRIwrite(mri_bias, "bias.mgz") ;*/
  if (nctrl == 0)  // smooth it quite a bit
  {
    MRI *mri_kernel, *mri_smooth, *mri_down;
    float sigma = 16.0f;

    mri_down = MRIdownsample2(mri_bias, NULL);
    mri_kernel = MRIgaussian1d(sigma, 100);
    mri_smooth = MRIconvolveGaussian(mri_down, NULL, mri_kernel);
    MRIfree(&mri_bias);
    MRIfree(&mri_kernel);
    mri_bias = MRIupsample2(mri_smooth, NULL);
    if (mri_bias->width != mri_in->width || mri_bias->height != mri_in->height || mri_bias->depth != mri_in->depth) {
      MRI *mri_tmp;

      mri_tmp = MRIcloneDifferentType(mri_in, mri_bias->type);
      MRIextractInto(mri_bias, mri_tmp, 0, 0, 0, mri_bias->width, mri_bias->height, mri_bias->depth, 0, 0, 0);
      MRIfree(&mri_bias);
      mri_bias = mri_tmp;
    }
    sigma = 2.0f;
    MRIfree(&mri_down);
    MRIfree(&mri_smooth);
    mri_kernel = MRIgaussian1d(sigma, 100);
    mri_smooth = MRIconvolveGaussian(mri_bias, NULL, mri_kernel);
    MRIfree(&mri_bias);
    mri_bias = mri_smooth;
    MRIfree(&mri_kernel);
  }
  else  // if manually specified control points, don't let them be overwhelmed
  {
    MRI *mri_kernel, *mri_smooth;
    mri_smooth = MRIsoapBubble(mri_bias, mri_ctrl, NULL, 10, -1);
    mri_kernel = MRIgaussian1d(bias_sigma, 100);
    MRIconvolveGaussian(mri_smooth, mri_bias, mri_kernel);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      MRIwrite(mri_bias, "smooth_bias.mgz");
    }
    MRIfree(&mri_smooth);
    MRIfree(&mri_kernel);
  }

  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();

        bias = (float)MRISvox(mri_bias, x, y, z) / NO_BIAS;
        if (bias < 0) DiagBreak();

        MRIsampleVolumeFrame(mri_in, x, y, z, input_index, &val);
        val *= bias; /* corrected value */
        switch (mri_in->type) {
          case MRI_UCHAR:
            if (val < 0)
              val = 0;
            else if (val > 255)
              val = 255;
          default:
            break;
        }
        MRIsetVoxVal(mri_dst, x, y, z, input_index, val);
      }
    }
  }

  MRIfree(&mri_bias);
  MRIfree(&mri_ctrl);
  return;
}

MRI *GCAnormalizeSamplesAllChannels(MRI *mri_in,
                                    GCA *gca,
                                    GCA_SAMPLE *gcas,
                                    int nsamples,
                                    TRANSFORM *transform,
                                    char *ctl_point_fname,
                                    double bias_sigma)
{
  MRI *mri_dst;
  int input;

  mri_dst = MRIcopy(mri_in, NULL);

  if (gca->ninputs == 1 && mri_dst->nframes > 1)  // normalize each frame
  {
    for (input = 0; input < mri_dst->nframes; input++) {
      MRI *mri_frame = MRIcopyFrame(mri_dst, NULL, input, 0);
      GCAnormalizeSamplesOneChannel(mri_frame, gca, gcas, nsamples, transform, ctl_point_fname, 0, bias_sigma);
      MRIcopyFrame(mri_frame, mri_dst, 0, input);
      MRIfree(&mri_frame);
    }
  }
  else {
    for (input = 0; input < gca->ninputs; input++) {
      GCAnormalizeSamplesOneChannel(mri_dst, gca, gcas, nsamples, transform, ctl_point_fname, input, bias_sigma);
    }
  }

  return (mri_dst);
}

MRI *GCAnormalizeSamplesT1PD(
    MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas, int nsamples, TRANSFORM *transform, char *ctl_point_fname)
{
  MRI *mri_dst, *mri_ctrl, *mri_bias;
  int xv, yv, zv, n, x, y, z, width, height, depth, out_val, xn, yn, zn, num, total;
  double val;
  float bias;
  double mean, sigma;

  mri_dst = MRIclone(mri_in, NULL);
  mri_ctrl = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_UCHAR);
  MRIcopyHeader(mri_in, mri_ctrl);
  mri_bias = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_SHORT);
  if (!mri_bias)
    ErrorExit(ERROR_NOMEMORY,
              "GCAnormalize: could not allocate (%d,%d,%d,2) bias image",
              mri_in->width,
              mri_in->height,
              mri_in->depth);
  MRIcopyHeader(mri_in, mri_bias);

#define MAX_BIAS 1250
#define NO_BIAS 1000
#define MIN_BIAS 750

  if (ctl_point_fname) {
    MRI3dUseFileControlPoints(mri_ctrl, ctl_point_fname);
    MRInormAddFileControlPoints(mri_ctrl, CONTROL_MARKED, NULL);
  }
  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;

  /* add control points from file */
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        MRISvox(mri_bias, x, y, z) = NO_BIAS; /* by default */
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED)
        /* read from file */
        {
          int n, max_n;
          GC1D *gc;
          GCA_NODE *gcan;
          GCA_PRIOR *gcap;
          double max_p;

          if (!GCAsourceVoxelToNode(gca, mri_dst, transform, x, y, z, &xn, &yn, &zn)) {
            gcan = &gca->nodes[xn][yn][zn];
            gcap = getGCAP(gca, mri_dst, transform, x, y, z);
            if (gcap == NULL) {
              continue;
            }
            max_p = 0;
            for (max_n = -1, n = 0; n < gcan->nlabels; n++) {
              if ((0 == IS_WM(gcan->labels[n])) && (0 == IS_CEREBELLAR_WM(gcan->labels[n]))) {
                continue;
              }
              gc = &gcan->gcs[n];
              if (getPrior(gcap, gcan->labels[n]) >= max_p) {
                max_p = getPrior(gcap, gcan->labels[n]);
                max_n = n;
              }
            }
            if (max_n < 0)
            /* couldn't find any valid label at this location */
            {
              continue;
            }
            gc = &gcan->gcs[max_n];

            MRIsampleVolumeFrameType(mri_in, x, y, z, 1, SAMPLE_NEAREST, &val);
            if (FZERO(val)) {
              val = 1;
            }
            bias = (float)NO_BIAS * ((float)gc->means[1] / val);
            if (bias < 100 || bias > 5000) {
              DiagBreak();
            }
            if (bias < MIN_BIAS) {
              bias = MIN_BIAS;
            }
            if (bias > MAX_BIAS) {
              bias = MAX_BIAS;
            }

            MRISvox(mri_bias, x, y, z) = (short)nint(bias);
          }
          ////////////////////////////////////////////////
        }
      }
    }
  }

  TransformInvert(transform, mri_in);
  for (n = 0; n < nsamples; n++) {
    if (gcas[n].xp == Ggca_x && gcas[n].yp == Ggca_y && gcas[n].zp == Ggca_z) {
      DiagBreak();
    }

    if (!GCApriorToSourceVoxel(gca, mri_dst, transform, gcas[n].xp, gcas[n].yp, gcas[n].zp, &xv, &yv, &zv)) {
      if (xv == 181 && yv == 146 && zv == 128) {
        DiagBreak();
      }
      if (xv == Ggca_x && yv == Ggca_y && zv == Ggca_z) {
        DiagBreak();
      }
      if (gcas[n].label == 29 || gcas[n].label == 61) {
        gcas[n].label = 0;
        DiagBreak();
      }
      if (gcas[n].label > 0) {
        MRIvox(mri_ctrl, xv, yv, zv) = CONTROL_MARKED;
        MRIsampleVolumeFrameType(mri_in, xv, yv, zv, 1, SAMPLE_NEAREST, &val);
        if (FZERO(val)) {
          val = 1;
        }
        bias = (float)NO_BIAS * ((float)gcas[n].means[1] / val);
        if (bias < 100 || bias > 5000) {
          DiagBreak();
        }

        MRISvox(mri_bias, xv, yv, zv) = (short)nint(bias);
      }
      else {
        MRIvox(mri_ctrl, xv, yv, zv) = CONTROL_NONE;
      }
    }
  }

  /* now check for and remove outliers */
  mean = sigma = 0.0;
  for (num = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED) {
          num++;
          bias = (double)MRISvox(mri_bias, x, y, z);
          mean += bias;
          sigma += (bias * bias);
        }
      }
    }
  }

  if (num > 0) {
    mean /= (double)num;
    sigma = sqrt(sigma / (double)num - mean * mean);
    printf("bias field = %2.3f +- %2.3f\n", mean / NO_BIAS, sigma / NO_BIAS);
  }

  /* now check for and remove outliers */
  for (total = num = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        if (MRIvox(mri_ctrl, x, y, z) == CONTROL_MARKED) {
          bias = (double)MRISvox(mri_bias, x, y, z);
          total++;
          if (fabs(bias - mean) > 4 * sigma) {
            MRIvox(mri_ctrl, x, y, z) = CONTROL_NONE;
            num++;
            MRISvox(mri_bias, x, y, z) = NO_BIAS;
          }
        }
      }
    }
  }

  printf("%d of %d control points discarded\n", num, total);

  MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias);
/*  MRIwrite(mri_bias, "bias.mgz") ;*/
  MRIsoapBubble(mri_bias, mri_ctrl, mri_bias, 10, -1);
  /*  MRIwrite(mri_bias, "smooth_bias.mgz") ;*/

  width = mri_in->width;
  height = mri_in->height;
  depth = mri_in->depth;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        bias = (float)MRISvox(mri_bias, x, y, z) / NO_BIAS;
        if (bias < 0) {
          DiagBreak();
        }
        MRIsampleVolumeFrameType(mri_in, x, y, z, 1, SAMPLE_NEAREST, &val);
        out_val = nint((float)val * bias);
        MRISseq_vox(mri_dst, x, y, z, 1) = (short)nint(out_val);
      }
    }
  }

  MRIcopyFrame(mri_in, mri_dst, 0, 0); /* copy over T1 frame */
  MRIfree(&mri_bias);
  MRIfree(&mri_ctrl);
  return (mri_dst);
}

float GCAlabelProbability(MRI *mri_src, GCA *gca, TRANSFORM *transform, float x, float y, float z, int label)
{
  int xn, yn, zn;
  // GCA_NODE *gcan;
  GC1D *gc;
  float plabel, vals[MAX_GCA_INPUTS];

  if (x == Gx && y == Gy && z == Gz) {
    DiagBreak();
  }

  if (!GCAsourceVoxelToNode(gca, mri_src, transform, x, y, z, &xn, &yn, &zn)) {
    // gcan = &gca->nodes[xn][yn][zn];

    load_vals(mri_src, x, y, z, vals, gca->ninputs);
    gc = GCAfindGC(gca, xn, yn, zn, label);
    if (gc == NULL || gc->ntraining == 0) {
      gc = GCAfindClosestValidGC(gca, xn, yn, zn, label, 0);
    }

    if (gc == NULL) {
      plabel = 0.0;
    }
    else {
      plabel = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, label);
    }
  }
  else {
    plabel = 0.0;
  }
  return (plabel);
}

static GCA_NODE *findSourceGCAN(GCA *gca, MRI *mri_src, TRANSFORM *transform, int x, int y, int z)
{
  int xn, yn, zn;
  GCA_NODE *gcan = NULL;

  if (!GCAsourceVoxelToNode(gca, mri_src, transform, x, y, z, &xn, &yn, &zn)) {
    gcan = &gca->nodes[xn][yn][zn];
  }
  return (gcan);
}
MRI *GCAmaxPosteriorBorders(
    GCA *gca, MRI *mri_inputs, MRI *mri_src, MRI *mri_dst, TRANSFORM *transform, int max_iter, float min_ratio)
{
  int nchanged, x, y, z, width, height, depth, label, total_changed, i;
  MRI *mri_tmp;

  if (mri_src != mri_dst) {
    mri_dst = MRIcopy(mri_src, mri_dst);
  }

  mri_tmp = MRIcopy(mri_dst, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (total_changed = i = 0; i < max_iter; i++) {
    nchanged = 0;
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          if (x == 99 && y == 129 && z == 127) {
            DiagBreak(); /* gray should be wm */
          }
          if (x == 98 && y == 124 && z == 127) {
            DiagBreak(); /* wm should be hippo */
          }

          if (borderVoxel(mri_dst, x, y, z)) {
            label = GCAmaxLikelihoodBorderLabel(gca, mri_inputs, mri_dst, transform, x, y, z, min_ratio);
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z &&
                (label == Ggca_label || nint(MRIgetVoxVal(mri_tmp, x, y, z, 0)) == Ggca_label || Ggca_label < 0)) {
              DiagBreak();

              if (label != nint(MRIgetVoxVal(mri_dst, x, y, z, 0)))
                printf(
                    "MLE (%d, %d, %d): old label %s (%d), "
                    "new label %s (%d)\n",
                    x,
                    y,
                    z,
                    cma_label_to_name(nint(MRIgetVoxVal(mri_tmp, x, y, z, 0))),
                    nint(MRIgetVoxVal(mri_tmp, x, y, z, 0)),
                    cma_label_to_name(label),
                    label);
            }
            if (label != nint(MRIgetVoxVal(mri_dst, x, y, z, 0))) {
              nchanged++;
              MRIsetVoxVal(mri_tmp, x, y, z, 0, label);
            }
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst);
    total_changed += nchanged;
    if (!nchanged) {
      break;
    }
  }

  MRIfree(&mri_tmp);
  printf("%d border labels changed to MLE ...\n", total_changed);
  return (mri_dst);
}
static int borderVoxel(MRI *mri, int x, int y, int z)
{
  int xi, yi, zi, xk, yk, zk, label;

  label = nint(MRIgetVoxVal(mri, x, y, z, 0));

  for (xk = -1; xk <= 1; xk++) {
    xi = mri->xi[x + xk];
    for (yk = -1; yk <= 1; yk++) {
      for (zk = -1; zk <= 1; zk++) {
        if (abs(xk) + abs(yk) + abs(zk) != 1) {
          continue;
        }
        yi = mri->yi[y + yk];
        zi = mri->zi[z + zk];
        if (nint(MRIgetVoxVal(mri, xi, yi, zi, 0)) != label) {
          return (1);
        }
      }
    }
  }
  return (0);
}

static int GCAmaxLikelihoodBorderLabel(
    GCA *gca, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform, int x, int y, int z, float min_ratio)
{
  float p, max_p, vals[MAX_GCA_INPUTS];
  int label, i, xi, yi, zi, best_label, orig_label, n;
  GCA_NODE *gcan;
  GC1D *gc;

  if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
    DiagBreak();
  }

  load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
  orig_label = best_label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));

  // current GCA_NODE at this point
  gcan = findSourceGCAN(gca, mri_inputs, transform, x, y, z);
  if (gcan == NULL) {
    return (orig_label);
  }

  // current label
  // look for classifier for this label and get the p value
  gc = NULL;
  for (n = 0; n < gcan->nlabels; n++)
    if (gcan->labels[n] == best_label) {
      gc = &gcan->gcs[n];
      break;
    }

  if (gc == NULL) {
    ErrorPrintf(ERROR_BADPARM,
                "GCAmaxLikelihoodBorderLabel(%d, %d, %d): "
                "couldn't find gc for label %d",
                x,
                y,
                z,
                best_label);
    max_p = 0.0;
  }
  else {
    max_p = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, best_label);
  }

  // look around neighbors ///////////////////////////////
  for (i = 0; i < GIBBS_NEIGHBORS; i++) {
    xi = mri_inputs->xi[x + xnbr_offset[i]];
    yi = mri_inputs->yi[y + ynbr_offset[i]];
    zi = mri_inputs->zi[z + znbr_offset[i]];
    gcan = findSourceGCAN(gca, mri_inputs, transform, xi, yi, zi);
    if (gcan == NULL) {
      continue;
    }
    // get the neighbor label
    label = nint(MRIgetVoxVal(mri_labels, xi, yi, zi, 0));
    gc = NULL;
    for (n = 0; n < gcan->nlabels; n++)
      if (gcan->labels[n] == label) {
        // get the classifier for this label at this location
        gc = &gcan->gcs[n];
        break;
      }
    if (gc == NULL) {
      continue; /* label can't occur here */
    }

    // set the label to this neighbor value
    MRIsetVoxVal(mri_labels, x, y, z, 0, label);
    // check if possible
    if (gcaGibbsImpossibleConfiguration(gca, mri_labels, x, y, z, transform)) {
      MRIsetVoxVal(mri_labels, x, y, z, 0, orig_label);  // added by xh
      continue;                                          // shouldn't put the original label back ???? -xh
    }
    // restore the old value
    MRIsetVoxVal(mri_labels, x, y, z, 0, orig_label);

    // calculate p for this neighbor label
    p = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, label);
    //
    if (((best_label == orig_label && p > min_ratio * max_p) ||
         // starting loop
         (best_label != orig_label && p >= max_p)) &&
        // later in the loop
        GCAisPossible(gca, mri_labels, label, transform, x, y, z, 0)) {
      max_p = p;
      best_label = label;
    }
  }

  /* test to make sure that it is not an impossible Gibbs configuration */
  if (best_label != nint(MRIgetVoxVal(mri_labels, x, y, z, 0))) {
    label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));
    MRIsetVoxVal(mri_labels, x, y, z, 0, best_label); /* test potential new label */
    if (gcaGibbsImpossibleConfiguration(gca, mri_labels, x, y, z, transform)) {
      best_label = label; /* revert back to old label */
    }
    MRIsetVoxVal(mri_labels, x, y, z, 0, label);
    /* caller will change it if needed */
  }
  return (best_label);
}

int GCAcomputeLabelStats(GCA *gca, int target_label, float *pvar, float *means)
{
  int x, y, z, n, r;
  double var, dof, total_dof;
  GC1D *gc;
  GCA_NODE *gcan;
  float fval;

  var = total_dof = 0.0;
  memset(means, 0, gca->ninputs * sizeof(float));
  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];

        for (n = 0; n < gcan->nlabels; n++) {
          if (gcan->labels[n] == target_label) {
            gc = &gcan->gcs[n];
            fval = get_node_prior(gca, target_label, x, y, z);
            if (fval != 0) {
              dof = get_node_prior(gca, target_label, x, y, z) * gcan->total_training;
              for (r = 0; r < gca->ninputs; r++) {
                means[r] += dof * gc->means[r];
              }
              var += dof * covariance_determinant(gc, gca->ninputs);
              total_dof += dof;
            }
          }
        }
      }
    }
  }

  if (total_dof > 0.0) {
    for (r = 0; r < gca->ninputs; r++) {
      means[r] /= total_dof;
    }
    var /= total_dof;
  }
  if (pvar) {
    *pvar = var;
  }
  return (NO_ERROR);
}
int GCAhistogramTissueStatistics(
    GCA *gca, MRI *mri_T1, MRI *mri_PD, MRI *mri_labeled, TRANSFORM *transform, const char *fname)
{
  int x, y, z, n, label, biggest_label, T1, PD, xp, yp, zp;
  GCA_NODE *gcan;
  GCA_TISSUE_PARMS *gca_tp;
  VECTOR *v_parc, *v_T1;
  static VECTOR *v_ras_cor = NULL, *v_ras_flash;
  FILE *fp;
  float xf, yf, zf;

  fp = fopen(fname, "w");
  if (!fp) ErrorExit(ERROR_NOFILE, "GCAhistogramTissueStatistics: could not open %s", fname);

  v_parc = VectorAlloc(4, MATRIX_REAL);
  v_T1 = VectorAlloc(4, MATRIX_REAL);
  *MATRIX_RELT(v_parc, 4, 1) = 1.0;
  *MATRIX_RELT(v_T1, 4, 1) = 1.0;

  /* first build a list of all labels that exist */
  for (biggest_label = x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_height; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          gca->tissue_parms[(int)gcan->labels[n]].label = gcan->labels[n];
          if (gcan->labels[n] > biggest_label) {
            biggest_label = gcan->labels[n];
          }
        }
      }
    }
  }

  for (label = 0; label <= biggest_label; label++) {
    if (gca->tissue_parms[label].label <= 0) {
      continue;
    }
    gca_tp = &gca->tissue_parms[label];

    for (z = 0; z < mri_T1->depth; z++) {
      V3_Z(v_parc) = z;
      for (y = 0; y < mri_T1->height; y++) {
        V3_Y(v_parc) = y;
        for (x = 0; x < mri_T1->width; x++) {
          if (nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) != label) {
            continue;
          }
          if (borderVoxel(mri_labeled, x, y, z)) {
            continue;
          }
          if (transform) {
            MATRIX *m_tmp;
            V3_X(v_parc) = x;

            TransformSample(transform, x * mri_T1->xsize, y * mri_T1->ysize, z * mri_T1->zsize, &xf, &yf, &zf);
            xp = nint(xf);
            yp = nint(zf);
            zp = nint(zf);
            V3_X(v_T1) = xp;
            V3_Y(v_T1) = yp;
            V3_Z(v_T1) = zp;

            m_tmp = MRIgetVoxelToRasXform(mri_labeled);
            v_ras_cor = MatrixMultiply(m_tmp, v_parc, v_ras_cor);
            MatrixFree(&m_tmp);
            m_tmp = MRIgetVoxelToRasXform(mri_T1);
            v_ras_flash = MatrixMultiply(m_tmp, v_T1, v_ras_flash);
            MatrixFree(&m_tmp);
            if (!x && !y && !z && 0) {
              MatrixPrint(stdout, v_ras_cor);
              MatrixPrint(stdout, v_ras_flash);
            }

            if ((xp < 0 || xp >= mri_T1->width) || (yp < 0 || yp >= mri_T1->height) ||
                (zp < 0 || zp >= mri_T1->depth)) {
              continue;
            }
          }
          else {
            xp = x;
            yp = y;
            zp = z;
          }

          T1 = MRISvox(mri_T1, xp, yp, zp);
          PD = MRISvox(mri_PD, xp, yp, zp);
          fprintf(fp, "%d %d %d\n", label, T1, PD);
          gca_tp->total_training++;
          gca_tp->T1_mean += T1;
          gca_tp->T1_var += T1 * T1;
          gca_tp->PD_mean += PD;
          gca_tp->PD_var += PD * PD;
        }
      }
    }
  }

  fclose(fp);
  return (NO_ERROR);
}

int GCAnormalizeTissueStatistics(GCA *gca)
{
  int n;
  double nsamples;
  GCA_TISSUE_PARMS *gca_tp;

  for (n = 0; n < MAX_GCA_LABELS; n++) {
    gca_tp = &gca->tissue_parms[n];
    if (gca_tp->total_training <= 0) {
      continue;
    }
    nsamples = gca_tp->total_training;
    gca_tp->T1_mean /= nsamples;
    gca_tp->PD_mean /= nsamples;
    gca_tp->T2_mean /= nsamples;
    gca_tp->T1_var = gca_tp->T1_var / nsamples - gca_tp->T1_mean * gca_tp->T1_mean;
    gca_tp->PD_var = gca_tp->PD_var / nsamples - gca_tp->PD_mean * gca_tp->PD_mean;
    printf("%s: T1=%4d +- %4d, PD=%4d +- %4d \n",
           cma_label_to_name(n),
           nint(gca_tp->T1_mean),
           nint(sqrt(gca_tp->T1_var)),
           nint(gca_tp->PD_mean),
           nint(sqrt(gca_tp->PD_var)));
  }

  return (NO_ERROR);
}

const char *cma_label_to_name(int label)
{
  static char name[STRLEN];

  if (label == Unknown) {
    return ("Unknown");
  }
  if (label == Left_Cerebral_Exterior) {
    return ("Left_Cerebral_Exterior");
  }
  if (label == Left_Cerebral_White_Matter) {
    return ("Left_Cerebral_White_Matter");
  }
  if (label == Left_Cerebral_Cortex) {
    return ("Left_Cerebral_Cortex");
  }
  if (label == SUSPICIOUS) {
    return ("SUSPICIOUS");
  }
  if (label == Left_Lateral_Ventricle) {
    return ("Left_Lateral_Ventricle");
  }
  if (label == Left_Inf_Lat_Vent) {
    return ("Left_Inf_Lat_Vent");
  }
  if (label == Left_Cerebellum_Exterior) {
    return ("Left_Cerebellum_Exterior");
  }
  if (label == Left_Cerebellum_White_Matter) {
    return ("Left_Cerebellum_White_Matter");
  }
  if (label == Left_Cerebellum_Cortex) {
    return ("Left_Cerebellum_Cortex");
  }
  if (label == Left_Thalamus) {
    return ("Left_Thalamus");
  }
  if (label == Left_Caudate) {
    return ("Left_Caudate");
  }
  if (label == Left_Putamen) {
    return ("Left_Putamen");
  }
  if (label == Left_Pallidum) {
    return ("Left_Pallidum");
  }
  if (label == Third_Ventricle) {
    return ("Third_Ventricle");
  }
  if (label == Fourth_Ventricle) {
    return ("Fourth_Ventricle");
  }
  if (label == Brain_Stem) {
    return ("Brain_Stem");
  }
  if (label == Left_Hippocampus) {
    return ("Left_Hippocampus");
  }
  if (label == Left_Amygdala) {
    return ("Left_Amygdala");
  }
  if (label == Left_Amygdala_Anterior) {
    return ("Left_Amygdala_Anterior");
  }
  if (label == Left_Insula) {
    return ("Left_Insula");
  }
  if (label == Left_Operculum) {
    return ("Left_Operculum");
  }
  if (label == Line_1) {
    return ("Line_1");
  }
  if (label == Line_2) {
    return ("Line_2");
  }
  if (label == Line_3) {
    return ("Line_3");
  }
  if (label == CSF) {
    return ("CSF");
  }
  if (label == Left_Lesion) {
    return ("Left_Lesion");
  }
  if (label == Left_Accumbens_area) {
    return ("Left_Accumbens_area");
  }
  if (label == Left_Substancia_Nigra) {
    return ("Left_Substancia_Nigra");
  }
  if (label == Left_VentralDC) {
    return ("Left_VentralDC");
  }
  if (label == Left_undetermined) {
    return ("Left_undetermined");
  }
  if (label == Left_vessel) {
    return ("Left_vessel");
  }
  if (label == Left_choroid_plexus) {
    return ("Left_choroid_plexus");
  }
  if (label == Left_F3orb) {
    return ("Left_F3orb");
  }
  if (label == Left_lOg) {
    return ("Left_lOg");
  }
  if (label == Left_aOg) {
    return ("Left_aOg");
  }
  if (label == Left_mOg) {
    return ("Left_mOg");
  }
  if (label == Left_pOg) {
    return ("Left_pOg");
  }
  if (label == Left_Stellate) {
    return ("Left_Stellate");
  }
  if (label == Left_Porg) {
    return ("Left_Porg");
  }
  if (label == Left_Aorg) {
    return ("Left_Aorg");
  }
  if (label == Right_Cerebral_Exterior) {
    return ("Right_Cerebral_Exterior");
  }
  if (label == Right_Cerebral_White_Matter) {
    return ("Right_Cerebral_White_Matter");
  }
  if (label == Right_Cerebral_Cortex) {
    return ("Right_Cerebral_Cortex");
  }
  if (label == Right_Lateral_Ventricle) {
    return ("Right_Lateral_Ventricle");
  }
  if (label == Right_Inf_Lat_Vent) {
    return ("Right_Inf_Lat_Vent");
  }
  if (label == Right_Cerebellum_Exterior) {
    return ("Right_Cerebellum_Exterior");
  }
  if (label == Right_Cerebellum_White_Matter) {
    return ("Right_Cerebellum_White_Matter");
  }
  if (label == Right_Cerebellum_Cortex) {
    return ("Right_Cerebellum_Cortex");
  }
  if (label == Right_Thalamus) {
    return ("Right_Thalamus");
  }
  if (label == Right_Caudate) {
    return ("Right_Caudate");
  }
  if (label == Right_Putamen) {
    return ("Right_Putamen");
  }
  if (label == Right_Pallidum) {
    return ("Right_Pallidum");
  }
  if (label == Right_Hippocampus) {
    return ("Right_Hippocampus");
  }
  if (label == Right_Amygdala) {
    return ("Right_Amygdala");
  }
  if (label == Right_Amygdala_Anterior) {
    return ("Right_Amygdala_Anterior");
  }
  if (label == Right_Insula) {
    return ("Right_Insula");
  }
  if (label == Right_Operculum) {
    return ("Right_Operculum");
  }
  if (label == Right_Lesion) {
    return ("Right_Lesion");
  }
  if (label == Right_Accumbens_area) {
    return ("Right_Accumbens_area");
  }
  if (label == Right_Substancia_Nigra) {
    return ("Right_Substancia_Nigra");
  }
  if (label == Right_VentralDC) {
    return ("Right_VentralDC");
  }
  if (label == Right_undetermined) {
    return ("Right_undetermined");
  }
  if (label == Right_vessel) {
    return ("Right_vessel");
  }
  if (label == Right_choroid_plexus) {
    return ("Right_choroid_plexus");
  }
  if (label == Right_F3orb) {
    return ("Right_F3orb");
  }
  if (label == Right_lOg) {
    return ("Right_lOg");
  }
  if (label == Right_aOg) {
    return ("Right_aOg");
  }
  if (label == Right_mOg) {
    return ("Right_mOg");
  }
  if (label == Right_pOg) {
    return ("Right_pOg");
  }
  if (label == Right_Stellate) {
    return ("Right_Stellate");
  }
  if (label == Right_Porg) {
    return ("Right_Porg");
  }
  if (label == Right_Aorg) {
    return ("Right_Aorg");
  }
  if (label == Bone) {
    return ("Bone");
  }
  if (label == Fat) {
    return ("Fat");
  }
  if (label == Bright_Unknown) {
    return ("Bright Unknown");
  }
  if (label == Dark_Unknown) {
    return ("Dark Unknown");
  }

  if (label == Left_Interior) {
    return ("Left_Interior");
  }
  if (label == Right_Interior) {
    return ("Right_Interior");
  }
  if (label == WM_hypointensities) {
    return ("WM_hypointensities");
  }
  if (label == Left_future_WMSA) {
    return ("Left_future_WMSA");
  }
  if (label == Right_future_WMSA) {
    return ("Right_future_WMSA");
  }
  if (label == Left_WM_hypointensities) {
    return ("Left_WM_hypointensities");
  }
  if (label == Right_WM_hypointensities) {
    return ("Right_WM_hypointensities");
  }
  if (label == non_WM_hypointensities) {
    return ("non_WM_hypointensities");
  }
  if (label == Left_non_WM_hypointensities) {
    return ("Left_non_WM_hypointensities");
  }
  if (label == Right_non_WM_hypointensities) {
    return ("Right_non_WM_hypointensities");
  }
  if (label == Fifth_Ventricle) {
    return ("Fifth_Ventricle");
  }
  if (label == Optic_Chiasm) {
    return ("Optic_Chiasm");
  }
  if (label == Cranium) {
    return ("Cranium");
  }
  if (label == Dura) {
    return ("Dura");
  }
  if (label == CSF_SA) {
    return ("CSF_SA");
  }
  if (label == Ear) {
    return ("Ear");
  }
  if (label == Muscle) {
    return ("Muscle");
  }
  if (label == Epidermis) {
    return ("Epidermis");
  }
  if (label == Conn_Tissue) {
    return ("Conn_Tissue");
  }
  if (label == SC_FAT_MUSCLE) {
    return ("SC-Fat/Muscle");
  }
  if (label == Fatty_Tissue) {
    return ("Fatty_Tissue");
  }
  if (label == Spinal_Cord) {
    return ("Spinal_Cord");
  }
  if (label == Soft_Tissue) {
    return ("Soft_Tissue");
  }
  if (label == Nerve) {
    return ("Nerve");
  }
  if (label == Bone) {
    return ("Bone");
  }
  if (label == Air) {
    return ("Air");
  }
  if (label == Orbit) {
    return ("Orbit");
  }
  if (label == Tongue) {
    return ("Tongue");
  }
  if (label == Nasal_Structures) {
    return ("Nasal_Structures");
  }
  if (label == Globe) {
    return ("Globe");
  }
  if (label == Teeth) {
    return ("Teeth");
  }

  if (label == alveus) {
    return ("alveus");
  }
  if (label == perforant_pathway) {
    return ("perforant_pathway");
  }
  if (label == parasubiculum) {
    return ("parasubiculum");
  }
  if (label == presubiculum) {
    return ("presubiculum");
  }
  if (label == subiculum) {
    return ("subiculum");
  }
  if (label == CA1) {
    return ("CA1");
  }
  if (label == CA2) {
    return ("CA2");
  }
  if (label == CA3) {
    return ("CA3");
  }
  if (label == CA4) {
    return ("CA4");
  }
  if (label == GC_DG) {
    return ("GC_DG");
  }
  if (label == HATA) {
    return ("HATA");
  }
  if (label == fimbria) {
    return ("fimbria");
  }
  if (label == lateral_ventricle) {
    return ("lateral_ventricle");
  }
  if (label == molecular_layer_HP) {
    return ("molecular_layer_HP");
  }
  if (label == hippocampal_fissure) {
    return ("hippocampal_fissure");
  }
  if (label == entorhinal_cortex) {
    return ("entorhinal_cortex");
  }
  if (label == molecular_layer_subiculum) {
    return ("molecular_layer_subiculum");
  }
  if (label == Amygdala) {
    return ("Amygdala");
  }
  if (label == Cerebral_White_Matter) {
    return ("Cerebral_White_Matter");
  }
  if (label == Cerebral_Cortex) {
    return ("Cerebral_Cortex");
  }
  if (label == Inf_Lat_Vent) {
    return ("Inf_Lat_Vent");
  }
  if (Fornix == label) {
    return ("Fornix");
  }

  if (label == BA17) {
    return ("BA17");
  }
  if (label == BA18) {
    return ("BA18");
  }
  if (label == BA44) {
    return ("BA44");
  }
  if (label == BA45) {
    return ("BA45");
  }
  if (label == BA4a) {
    return ("BA4a");
  }
  if (label == BA4p) {
    return ("BA4p");
  }
  if (label == BA6) {
    return ("BA6");
  }
  if (label == BA2) {
    return ("BA2");
  }
  if (label == BAun1) {
    return ("BAun1");
  }
  if (label == BAun2) {
    return ("BAun2");
  }
  if (label == right_CA2_3) {
    return ("right_CA2_3");
  }
  if (label == right_alveus) {
    return ("right_alveus");
  }
  if (label == right_CA1) {
    return ("right_CA1");
  }
  if (label == right_fimbria) {
    return ("right_fimbria");
  }
  if (label == right_presubiculum) {
    return ("right_presubiculum");
  }
  if (label == right_hippocampal_fissure) {
    return ("right_hippocampal_fissure");
  }
  if (label == right_CA4_DG) {
    return ("right_CA4_DG");
  }
  if (label == right_subiculum) {
    return ("right_subiculum");
  }
  if (label == left_CA2_3) {
    return ("left_CA2_3");
  }
  if (label == left_alveus) {
    return ("left_alveus");
  }
  if (label == left_CA1) {
    return ("left_CA1");
  }
  if (label == left_fimbria) {
    return ("left_fimbria");
  }
  if (label == left_presubiculum) {
    return ("left_presubiculum");
  }
  if (label == left_hippocampal_fissure) {
    return ("left_hippocampal_fissure");
  }
  if (label == left_CA4_DG) {
    return ("left_CA4_DG");
  }
  if (label == left_subiculum) {
    return ("left_subiculum");
  }

  if (label == CC_Posterior) {
    return ("CC_Posterior");
  }
  if (label == CC_Mid_Posterior) {
    return ("CC_Mid_Posterior");
  }
  if (label == CC_Central) {
    return ("CC_Central");
  }
  if (label == CC_Mid_Anterior) {
    return ("CC_Mid_Anterior");
  }
  if (label == CC_Anterior) {
    return ("CC_Anterior");
  }
  if (label == Corpus_Callosum) {
    return ("Corpus_Callosum");
  }

  if (label == left_fornix) {
    return ("left_fornix");
  }
  if (label == right_fornix) {
    return ("right_fornix");
  }
  if (label == Tumor) {
    return ("Tumor");
  }

  if (label == lh_cst) {
    return ("Left Corticospinal Tract");
  }
  if (label == rh_cst) {
    return ("Right Corticospinal Tract");
  }
  if (label == lh_ilf) {
    return ("Left Inferior Longitudinal Fasciculus");
  }
  if (label == rh_ilf) {
    return ("Right Inferior Longitudinal Fasciculus");
  }
  if (label == lh_unc) {
    return ("Left Uncinate Fasciculus");
  }
  if (label == rh_unc) {
    return ("Right Uncinate Fasciculus");
  }
  if (label == fmajor) {
    return ("Corpus Callosum Forceps Major");
  }
  if (label == fminor) {
    return ("Corpus Callosum Forceps Minor");
  }
  if (label == lh_atr) {
    return ("Left Anterior Thalamic Radiation");
  }
  if (label == rh_atr) {
    return ("Right Anterior Thalamic Radiation");
  }
  if (label == lh_ccg) {
    return ("Left Cingulum - Cingulate Gyrus");
  }
  if (label == rh_ccg) {
    return ("Right Cingulum - Cingulate Gyrus");
  }
  if (label == lh_cab) {
    return ("Left Cingulum - Angular Bundle");
  }
  if (label == rh_cab) {
    return ("Right Cingulum - Angular Bundle");
  }
  if (label == lh_slfp) {
    return ("Left Superior Longitudinal Fasciculus - Parietal");
  }
  if (label == rh_slfp) {
    return ("Right Superior Longitudinal Fasciculus - Parietal");
  }
  if (label == lh_slft) {
    return ("Left Superior Longitudinal Fasciculus - Temporal");
  }
  if (label == rh_slft) {
    return ("Right Superior Longitudinal Fasciculus - Temporal");
  }
  if (label == lh_ifof) {
    return ("Left Inferior Fronto-Occipital Fasciculus  ");
  }
  if (label == rh_ifof) {
    return ("Right Inferior Fronto-Occipital Fasciculus  ");
  }  
  if (label == lh_fornix) {
    return ("Left Fornix");
  }
  if (label == rh_fornix) {
    return ("Right Fornix");
  }
  if (label == Cbm_Left_I_IV) {
    return ("Cbm_Left_I_IV");
  }
  if (label == Cbm_Right_I_IV) {
    return ("Cbm_Right_I_IV");
  }
  if (label == Cbm_Left_V) {
    return ("Cbm_Left_V");
  }
  if (label == Cbm_Right_V) {
    return ("Cbm_Right_V");
  }
  if (label == Cbm_Left_VI) {
    return ("Cbm_Left_VI");
  }
  if (label == Cbm_Vermis_VI) {
    return ("Cbm_Vermis_VI");
  }
  if (label == Cbm_Right_VI) {
    return ("Cbm_Right_VI");
  }
  if (label == Cbm_Left_CrusI) {
    return ("Cbm_Left_CrusI");
  }
  if (label == Cbm_Vermis_CrusI) {
    return ("Cbm_Vermis_CrusI");
  }
  if (label == Cbm_Right_CrusI) {
    return ("Cbm_Right_CrusI");
  }
  if (label == Cbm_Left_CrusII) {
    return ("Cbm_Left_CrusII");
  }
  if (label == Cbm_Vermis_CrusII) {
    return ("Cbm_Vermis_CrusII");
  }
  if (label == Cbm_Right_CrusII) {
    return ("Cbm_Right_CrusII");
  }
  if (label == Cbm_Left_VIIb) {
    return ("Cbm_Left_VIIb");
  }
  if (label == Cbm_Vermis_VIIb) {
    return ("Cbm_Vermis_VIIb");
  }
  if (label == Cbm_Right_VIIb) {
    return ("Cbm_Right_VIIb");
  }
  if (label == Cbm_Left_VIIIa) {
    return ("Cbm_Left_VIIIa");
  }
  if (label == Cbm_Vermis_VIIIa) {
    return ("Cbm_Vermis_VIIIa");
  }
  if (label == Cbm_Right_VIIIa) {
    return ("Cbm_Right_VIIIa");
  }
  if (label == Cbm_Left_VIIIb) {
    return ("Cbm_Left_VIIIb");
  }
  if (label == Cbm_Vermis_VIIIb) {
    return ("Cbm_Vermis_VIIIb");
  }
  if (label == Cbm_Right_VIIIb) {
    return ("Cbm_Right_VIIIb");
  }
  if (label == Cbm_Left_IX) {
    return ("Cbm_Left_IX");
  }
  if (label == Cbm_Vermis_IX) {
    return ("Cbm_Vermis_IX");
  }
  if (label == Cbm_Right_IX) {
    return ("Cbm_Right_IX");
  }
  if (label == Cbm_Left_X) {
    return ("Cbm_Left_X");
  }
  if (label == Cbm_Vermis_X) {
    return ("Cbm_Vermis_X");
  }
  if (label == Cbm_Right_X) {
    return ("Cbm_Right_X");
  }

  if (label == ctx_lh_unknown) return ("ctx_lh_unknown");
  if (label == ctx_lh_bankssts) return ("ctx_lh_bankssts");
  if (label == ctx_lh_caudalanteriorcingulate) return ("ctx_lh_caudalanteriorcingulate");
  if (label == ctx_lh_caudalmiddlefrontal) return ("ctx_lh_caudalmiddlefrontal");
  if (label == ctx_lh_corpuscallosum) return ("ctx_lh_corpuscallosum");
  if (label == ctx_lh_cuneus) return ("ctx_lh_cuneus");
  if (label == ctx_lh_entorhinal) return ("ctx_lh_entorhinal");
  if (label == ctx_lh_fusiform) return ("ctx_lh_fusiform");
  if (label == ctx_lh_inferiorparietal) return ("ctx_lh_inferiorparietal");
  if (label == ctx_lh_inferiortemporal) return ("ctx_lh_inferiortemporal");
  if (label == ctx_lh_isthmuscingulate) return ("ctx_lh_isthmuscingulate");
  if (label == ctx_lh_lateraloccipital) return ("ctx_lh_lateraloccipital");
  if (label == ctx_lh_lateralorbitofrontal) return ("ctx_lh_lateralorbitofrontal");
  if (label == ctx_lh_lingual) return ("ctx_lh_lingual");
  if (label == ctx_lh_medialorbitofrontal) return ("ctx_lh_medialorbitofrontal");
  if (label == ctx_lh_middletemporal) return ("ctx_lh_middletemporal");
  if (label == ctx_lh_parahippocampal) return ("ctx_lh_parahippocampal");
  if (label == ctx_lh_paracentral) return ("ctx_lh_paracentral");
  if (label == ctx_lh_parsopercularis) return ("ctx_lh_parsopercularis");
  if (label == ctx_lh_parsorbitalis) return ("ctx_lh_parsorbitalis");
  if (label == ctx_lh_parstriangularis) return ("ctx_lh_parstriangularis");
  if (label == ctx_lh_pericalcarine) return ("ctx_lh_pericalcarine");
  if (label == ctx_lh_postcentral) return ("ctx_lh_postcentral");
  if (label == ctx_lh_posteriorcingulate) return ("ctx_lh_posteriorcingulate");
  if (label == ctx_lh_precentral) return ("ctx_lh_precentral");
  if (label == ctx_lh_precuneus) return ("ctx_lh_precuneus");
  if (label == ctx_lh_rostralanteriorcingulate) return ("ctx_lh_rostralanteriorcingulate");
  if (label == ctx_lh_rostralmiddlefrontal) return ("ctx_lh_rostralmiddlefrontal");
  if (label == ctx_lh_superiorfrontal) return ("ctx_lh_superiorfrontal");
  if (label == ctx_lh_superiorparietal) return ("ctx_lh_superiorparietal");
  if (label == ctx_lh_superiortemporal) return ("ctx_lh_superiortemporal");
  if (label == ctx_lh_supramarginal) return ("ctx_lh_supramarginal");
  if (label == ctx_lh_frontalpole) return ("ctx_lh_frontalpole");
  if (label == ctx_lh_temporalpole) return ("ctx_lh_temporalpole");
  if (label == ctx_lh_transversetemporal) return ("ctx_lh_transversetemporal");
  if (label == ctx_lh_insula) return ("ctx_lh_insula");
  if (label == ctx_rh_unknown) return ("ctx_rh_unknown");
  if (label == ctx_rh_bankssts) return ("ctx_rh_bankssts");
  if (label == ctx_rh_caudalanteriorcingulate) return ("ctx_rh_caudalanteriorcingulate");
  if (label == ctx_rh_caudalmiddlefrontal) return ("ctx_rh_caudalmiddlefrontal");
  if (label == ctx_rh_corpuscallosum) return ("ctx_rh_corpuscallosum");
  if (label == ctx_rh_cuneus) return ("ctx_rh_cuneus");
  if (label == ctx_rh_entorhinal) return ("ctx_rh_entorhinal");
  if (label == ctx_rh_fusiform) return ("ctx_rh_fusiform");
  if (label == ctx_rh_inferiorparietal) return ("ctx_rh_inferiorparietal");
  if (label == ctx_rh_inferiortemporal) return ("ctx_rh_inferiortemporal");
  if (label == ctx_rh_isthmuscingulate) return ("ctx_rh_isthmuscingulate");
  if (label == ctx_rh_lateraloccipital) return ("ctx_rh_lateraloccipital");
  if (label == ctx_rh_lateralorbitofrontal) return ("ctx_rh_lateralorbitofrontal");
  if (label == ctx_rh_lingual) return ("ctx_rh_lingual");
  if (label == ctx_rh_medialorbitofrontal) return ("ctx_rh_medialorbitofrontal");
  if (label == ctx_rh_middletemporal) return ("ctx_rh_middletemporal");
  if (label == ctx_rh_parahippocampal) return ("ctx_rh_parahippocampal");
  if (label == ctx_rh_paracentral) return ("ctx_rh_paracentral");
  if (label == ctx_rh_parsopercularis) return ("ctx_rh_parsopercularis");
  if (label == ctx_rh_parsorbitalis) return ("ctx_rh_parsorbitalis");
  if (label == ctx_rh_parstriangularis) return ("ctx_rh_parstriangularis");
  if (label == ctx_rh_pericalcarine) return ("ctx_rh_pericalcarine");
  if (label == ctx_rh_postcentral) return ("ctx_rh_postcentral");
  if (label == ctx_rh_posteriorcingulate) return ("ctx_rh_posteriorcingulate");
  if (label == ctx_rh_precentral) return ("ctx_rh_precentral");
  if (label == ctx_rh_precuneus) return ("ctx_rh_precuneus");
  if (label == ctx_rh_rostralanteriorcingulate) return ("ctx_rh_rostralanteriorcingulate");
  if (label == ctx_rh_rostralmiddlefrontal) return ("ctx_rh_rostralmiddlefrontal");
  if (label == ctx_rh_superiorfrontal) return ("ctx_rh_superiorfrontal");
  if (label == ctx_rh_superiorparietal) return ("ctx_rh_superiorparietal");
  if (label == ctx_rh_superiortemporal) return ("ctx_rh_superiortemporal");
  if (label == ctx_rh_supramarginal) return ("ctx_rh_supramarginal");
  if (label == ctx_rh_frontalpole) return ("ctx_rh_frontalpole");
  if (label == ctx_rh_temporalpole) return ("ctx_rh_temporalpole");
  if (label == ctx_rh_transversetemporal) return ("ctx_rh_transversetemporal");
  if (label == ctx_rh_insula) return ("ctx_rh_insula");

  return (name);
}
MRI *GCArelabel_cortical_gray_and_white(GCA *gca, MRI *mri_inputs, MRI *mri_src, MRI *mri_dst, TRANSFORM *transform)
{
  int nchanged, x, y, z, width, height, depth, total_changed, label, xn, yn, zn, left, new_wm, new_gray;
  MRI *mri_tmp;
  // GCA_NODE *gcan;
  GC1D *gc_gray, *gc_wm;
  float vals[MAX_GCA_INPUTS], gray_dist, wm_dist;
  float grayPrior;
  float wmPrior;
  int xp, yp, zp;

  GCA_PRIOR *gcap = 0;

  if (mri_src != mri_dst) {
    mri_dst = MRIcopy(mri_src, mri_dst);
  }

  mri_tmp = MRIcopy(mri_dst, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  total_changed = new_wm = new_gray = 0;
  do {
    nchanged = 0;
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }

          label = nint(MRIgetVoxVal(mri_src, x, y, z, 0));
          if (label != Left_Cerebral_Cortex && label != Left_Cerebral_White_Matter && label != Right_Cerebral_Cortex &&
              label != Right_Cerebral_White_Matter) {
            continue;
          }

          load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
          if (!GCAsourceVoxelToNode(gca, mri_dst, transform, x, y, z, &xn, &yn, &zn)) {
            // gcan = &gca->nodes[xn][yn][zn];
            // label is left cortex or wm
            if (label == Left_Cerebral_Cortex || label == Left_Cerebral_White_Matter) {
              left = 1;
              gc_gray = GCAfindGC(gca, xn, yn, zn, Left_Cerebral_Cortex);
              gc_wm = GCAfindGC(gca, xn, yn, zn, Left_Cerebral_White_Matter);
            }
            // label is right cortex or wm
            else {
              gc_gray = GCAfindGC(gca, xn, yn, zn, Right_Cerebral_Cortex);
              gc_wm = GCAfindGC(gca, xn, yn, zn, Right_Cerebral_White_Matter);
              left = 0;
            }
            // can't be found
            if (!gc_wm || !gc_gray) {
              continue;
            }
            // calculate Mahalanobis distance
            gray_dist = sqrt(GCAmahDistIdentityCovariance(gc_gray, vals, gca->ninputs));
            wm_dist = sqrt(GCAmahDistIdentityCovariance(gc_wm, vals, gca->ninputs));

            // get the prior coordinate
            if (!GCAsourceVoxelToPrior(gca, mri_dst, transform, x, y, z, &xp, &yp, &zp)) {
              gcap = &gca->priors[xp][yp][zp];
              if (gcap == NULL) {
                continue;
              }
              // labeling gray but check Prior
              if (left) {
                grayPrior = getPrior(gcap, Left_Cerebral_Cortex);
                wmPrior = getPrior(gcap, Left_Cerebral_White_Matter);
              }
              else {
                grayPrior = getPrior(gcap, Right_Cerebral_Cortex);
                wmPrior = getPrior(gcap, Right_Cerebral_White_Matter);
              }
            }
            else {
              grayPrior = -1;
              wmPrior = -1;
            }
            // if grey < wm
            if (gray_dist < wm_dist && grayPrior > 0 && wmPrior > 0)
            // if ((5*gray_dist) < wm_dist &&
            // grayPrior > 0 && wmPrior > 0)
            {
              // if prior is not high, then you
              // can change white->gray
              if (wmPrior < 0.9)
                // label cortex
                label = left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex;
              // else
              // printf("(%d,%d,%d) had gray_dist(%.2f) "
              //"< wm_dist (%.2f),
              // but grayPrior(%.2f) < wmPrior(%.2f)\n",
              // x, y, z, gray_dist, wm_dist, grayPrior, wmPrior);
            }
            else {
              // if prior is not high,
              // then you can change gray->white
              if (grayPrior < 0.9 && grayPrior > 0 && wmPrior > 0)
                // label wm
                // if ((5*wm_dist < gray_dist)
                // && grayPrior > 0 && wmPrior > 0)
                label = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              // else
              // printf("(%d,%d,%d) had gray_dist(%.2f) > "
              //"wm_dist (%.2f), but grayPrior(%.2f)
              // > wmPrior(%.2f)\n",
              //    x, y, z, gray_dist, wm_dist,
              // grayPrior, wmPrior);
            }
            // if label changed from the current one
            // and it is possible
            if (label != nint(MRIgetVoxVal(mri_dst, x, y, z, 0)) &&
                GCAisPossible(gca, mri_dst, label, transform, x, y, z, 0)) {
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (label == Ggca_label || Ggca_label < 0)) {
                int in;
                VECTOR *v_means = 0;

                printf(
                    "relabel_cortical_gray_and_white: "
                    "voxel(%d,%d,%d), inputs=",
                    x,
                    y,
                    z);
                for (in = 0; in < gca->ninputs; in++) {
                  printf("%2.1f ", vals[in]);
                }
                printf("label=%s (%d)\n", cma_label_to_name(label), label);
                printf(" gray_dist = %.2f, wm_dist = %.2f\n", gray_dist, wm_dist);
                v_means = load_mean_vector(gc_gray, v_means, gca->ninputs);
                printf("v_means for gray = ");
                MatrixPrint(stdout, v_means);
                v_means = load_mean_vector(gc_wm, v_means, gca->ninputs);
                printf("v_means for wm = ");
                MatrixPrint(stdout, v_means);
                VectorFree(&v_means);
              }
              // count changed label
              if (label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) {
                new_gray++;
              }
              else {
                new_wm++;
              }
              nchanged++;
              MRIsetVoxVal(mri_tmp, x, y, z, 0, label);
            }
          }
          //////////////////////////////////
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst);
    total_changed += nchanged;
  } while (nchanged > 0);

  MRIfree(&mri_tmp);
  printf(
      "%d gm and wm labels changed (%%%2.0f to gray, %%%2.0f to "
      "white out of all changed labels)\n",
      total_changed,
      100.0f * (float)new_gray / total_changed,
      100.0f * (float)new_wm / total_changed);
  return (mri_dst);
}
int GCAdump(GCA *gca, MRI *mri, int x, int y, int z, TRANSFORM *transform, FILE *fp, int verbose)
{
  int xn, yn, zn, xp, yp, zp;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;

  if (!GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn)) {
    if (!GCAsourceVoxelToPrior(gca, mri, transform, x, y, z, &xp, &yp, &zp)) {
      printf(
          "\nGCA node at voxel (%d, %d, %d) --> node "
          "(%d, %d, %d), prior (%d, %d, %d)\n",
          x,
          y,
          z,
          xn,
          yn,
          zn,
          xp,
          yp,
          zp);
      gcan = &gca->nodes[xn][yn][zn];
      gcap = getGCAP(gca, mri, transform, x, y, z);
      if (gcap == NULL) {
        printf("\nGCAdump: prior point is outside.\n");
      }
      dump_gcan(gca, gcan, fp, verbose, gcap);
    }
    else {
      printf("\nGCAdump: prior point is outside.\n");
    }
  }
  else {
    printf("\nGCAdump: node point is outside.\n");
  }
  return (NO_ERROR);
}


static GCA_SAMPLE *gcaExtractThresholdedRegionLabelAsSamples(GCA *gca,
                                                             MRI *mri_labeled,
                                                             TRANSFORM *transform,
                                                             int *pnsamples,
                                                             int label,
                                                             int xp,
                                                             int yp,
                                                             int zp,
                                                             int wsize,
                                                             float pthresh)
{
  int i, x, y, z, xi, yi, zi, xk, yk, zk, whalf, r, c, v;
  // int  width, height, depth;
  int nsamples = 0;
  GCA_SAMPLE *gcas;
  GCA_PRIOR *gcap;
  GC1D *gc;
  float prior;

  // width = mri_labeled->width;
  // height = mri_labeled->height;
  // depth = mri_labeled->depth;
  whalf = (wsize - 1) / 2;
  gcap = &gca->priors[xp][yp][zp];
  if (gcap == NULL) {
    return NULL;
  }
  TransformInvert(transform, mri_labeled);
  if (!GCApriorToSourceVoxel(gca, mri_labeled, transform, xp, yp, zp, &x, &y, &z)) {
    for (nsamples = 0, zk = -whalf; zk <= whalf; zk++) {
      zi = mri_labeled->zi[z + zk];
      for (yk = -whalf; yk <= whalf; yk++) {
        yi = mri_labeled->yi[y + yk];
        for (xk = -whalf; xk <= whalf; xk++) {
          xi = mri_labeled->xi[x + xk];
          if (nint(MRIgetVoxVal(mri_labeled, xi, yi, zi, 0)) != label) {
            continue;
          }
          if (xi == Ggca_x && yi == Ggca_y && zi == Ggca_z) {
            DiagBreak();
          }
          nsamples++;
        }
      }
    }
  }

  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE));
  if (!gcas)
    ErrorExit(ERROR_NOMEMORY, "gcaExtractLabelAsSamples(%d): could not allocate %d samples\n", label, nsamples);

  /* go through region again and fill in samples */
  for (i = 0, zk = -whalf; zk <= whalf; zk++) {
    zi = mri_labeled->zi[z + zk];
    for (yk = -whalf; yk <= whalf; yk++) {
      yi = mri_labeled->yi[y + yk];
      for (xk = -whalf; xk <= whalf; xk++) {
        xi = mri_labeled->xi[x + xk];
        if (nint(MRIgetVoxVal(mri_labeled, xi, yi, zi, 0)) != label) {
          continue;
        }
        if (xi == Ggca_x && yi == Ggca_y && zi == Ggca_z) {
          DiagBreak();
        }

        gcap = getGCAP(gca, mri_labeled, transform, xi, yi, zi);
        if (gcap == NULL) {
          continue;
        }
        if (!GCAsourceVoxelToPrior(gca, mri_labeled, transform, xi, yi, zi, &xp, &yp, &zp)) {
          gc = GCAfindPriorGC(gca, xp, yp, zp, label);
          if (gcap) {
            prior = getPrior(gcap, label);
          }
          else {
            prior = 0;
          }
          if (!gc || !gcap || prior < pthresh) {
            nsamples--; /* shouldn't happen */
            continue;
          }

          gcas[i].label = label;
          gcas[i].xp = xp;
          gcas[i].yp = yp;
          gcas[i].zp = zp;
          gcas[i].x = xi;
          gcas[i].y = yi;
          gcas[i].z = zi;
          gcas[i].means = (float *)calloc(gca->ninputs, sizeof(float));
          gcas[i].covars = (float *)calloc((gca->ninputs * (gca->ninputs + 1)) / 2, sizeof(float));
          if (!gcas[i].means || !gcas[i].covars)
            ErrorExit(ERROR_NOMEMORY,
                      "GCArenormalizeAdapative: could not allocate "
                      "mean (%d) and covariance (%d) matrices",
                      gca->ninputs,
                      gca->ninputs * (gca->ninputs + 1) / 2);
          for (r = v = 0; r < gca->ninputs; r++) {
            gcas[i].means[r] = gc->means[r];
            for (c = r; c < gca->ninputs; c++) {
              gcas[i].covars[v] = gc->covars[v];
            }
          }
          gcas_setPrior(gcas[i],prior);
          if (FZERO(prior)) {
            DiagBreak();
          }
          i++;
        }
      }
    }
  }

  *pnsamples = nsamples;
  return (gcas);
}

static GCA_SAMPLE *gcaExtractLabelAsSamples(GCA *gca, MRI *mri_labeled, TRANSFORM *transform, int *pnsamples, int label)
{
  int i, nsamples, width, height, depth, x, y, z, xp, yp, zp, n, r, c, v;
  GCA_SAMPLE *gcas;
  GCA_PRIOR *gcap;
  GC1D *gc;

  width = mri_labeled->width;
  height = mri_labeled->height;
  depth = mri_labeled->depth;

  for (nsamples = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) != label) {
          continue;
        }
        nsamples++;
      }
    }
  }

  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE));
  if (!gcas)
    ErrorExit(ERROR_NOMEMORY, "gcaExtractLabelAsSamples(%d): could not allocate %d samples\n", label, nsamples);

  for (i = z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) != label) {
          continue;
        }

        if (!GCAsourceVoxelToPrior(gca, mri_labeled, transform, x, y, z, &xp, &yp, &zp)) {
          gcap = &gca->priors[xp][yp][zp];
          if (gcap == NULL) {
            continue;
          }
          for (n = 0; n < gcap->nlabels; n++) {
            if (gcap->labels[n] == label) {
              break;
            }
          }

          gc = GCAfindPriorGC(gca, xp, yp, zp, label);
          if (n >= gcap->nlabels || !gc) {
            nsamples--; /* doesn't exist at this location */
            continue;   /* ?? */
          }
          gcas[i].label = label;
          gcas[i].xp = xp;
          gcas[i].yp = yp;
          gcas[i].zp = zp;
          gcas[i].x = x;
          gcas[i].y = y;
          gcas[i].z = z;
          for (r = v = 0; r < gca->ninputs; r++) {
            gcas[i].means[r] = gc->means[r];
            for (c = r; c < gca->ninputs; c++) {
              gcas[i].covars[v] = gc->covars[v];
            }
          }
          gcas_setPrior(gcas[i], getPrior(gcap, label));
          if (FZERO(gcas_getPrior(gcas[i]))) {
            DiagBreak();
          }
          i++;  // this i is never used??????
        }       // !GCAsourceVoxelToPrior
      }
    }
  }

  *pnsamples = nsamples;
  return (gcas);
}

#define MIN_MEAN_SAMPLES 10
#define SAMPLE_PCT 0.20

int GCArenormalize(MRI *mri_in, MRI *mri_labeled, GCA *gca, TRANSFORM *transform)
{
  int x, y, z, n, label, biggest_label, nsamples, val, *ordered_indices, i, index, width, height, depth;
  float mean, var, *means, *stds, *gca_means;
  GCA_NODE *gcan;
  GCA_SAMPLE *gcas;
  GC1D *gc;

  if (gca->ninputs > 1) ErrorExit(ERROR_UNSUPPORTED, "GCArenormalize: can only renormalize scalars");

  /* first build a list of all labels that exist */
  for (nsamples = biggest_label = x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_height; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          if (gcan->labels[n] > biggest_label) {
            biggest_label = gcan->labels[n];
          }
        }
      }
    }
  }

  gca_means = (float *)calloc(biggest_label + 1, sizeof(float));
  means = (float *)calloc(biggest_label + 1, sizeof(float));
  stds = (float *)calloc(biggest_label + 1, sizeof(float));
  if (!gca_means || !means || !stds)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d vector", Progname, biggest_label + 1);

  /* do unknown labels separately */
  for (mean = var = 0.0, nsamples = x = 0; x < mri_in->width; x++) {
    for (y = 0; y < mri_in->height; y++) {
      for (z = 0; z < mri_in->height; z++) {
        if (nint(MRIgetVoxVal(mri_labeled, x, y, z, 0)) == Unknown) {
          nsamples++;
          val = MRIgetVoxVal(mri_in, x, y, z, 0);
          mean += val;
          var += (val * val);
        }
      }
    }
  }

  if (DIAG_VERBOSE_ON) {
    HISTOGRAM *histo;
    char fname[STRLEN];

    label = Right_Hippocampus;
    histo = MRIhistogramLabel(mri_in, mri_labeled, label, 256);
    sprintf(fname, "%s.plt", cma_label_to_name(label));
    HISTOplot(histo, fname);
    HISTOfree(&histo);
  }

  label = Unknown;
  if (!FZERO(nsamples)) {
    mean /= nsamples;
    means[label] = mean;
    stds[label] = sqrt(var / nsamples - mean * mean);
    GCAlabelMean(gca, label, &gca_means[label]);
    printf(
        "scaling label %s by %2.2f (%2.2f / %2.2f) "
        "(%d samples, std=%2.1f)\n",
        cma_label_to_name(label),
        means[label] / gca_means[label],
        means[label],
        gca_means[label],
        nsamples,
        stds[label]);
  }
  else {
    gca_means[label] = stds[label] = means[label] = 1.0;
  }

  gca_means[label] = stds[label] = means[label] = 1.0;

  for (label = 1; label <= biggest_label; label++) {
    gcas = gcaExtractLabelAsSamples(gca, mri_labeled, transform, &nsamples, label);
    if (!nsamples) {
      continue;
    }
    if (nsamples < MIN_MEAN_SAMPLES)
    /* not enough sample to estimate mean */
    {
      free(gcas);
      continue;
    }
    ordered_indices = (int *)calloc(nsamples, sizeof(int));
    GCAcomputeLogSampleProbability(gca, gcas, mri_in, transform, nsamples, DEFAULT_CLAMP);
    GCArankSamples(gca, gcas, nsamples, ordered_indices);

    if (nint(nsamples * SAMPLE_PCT) < MIN_MEAN_SAMPLES) {
      nsamples = MIN_MEAN_SAMPLES;
    }
    else {
      nsamples = nint(SAMPLE_PCT * (float)nsamples);
    }

    /* compute mean and variance of image intensities in this label */
    for (var = mean = 0.0f, i = 0; i < nsamples; i++) {
      index = ordered_indices[i];
      val = MRIgetVoxVal(mri_in, gcas[index].x, gcas[index].y, gcas[index].z, 0);
      mean += val;
      var += val * val;
    }
    mean /= (float)nsamples;
    var = var / nsamples - mean * mean;
    var = sqrt(var);
    means[label] = mean;
    stds[label] = var;
    GCAlabelMean(gca, label, &gca_means[label]);
    free(gcas);
    free(ordered_indices);
    printf(
        "scaling label %s by %2.2f (%2.2f / %2.2f) "
        "(%d samples, std=%2.1f)\n",
        cma_label_to_name(label),
        means[label] / gca_means[label],
        means[label],
        gca_means[label],
        nsamples,
        stds[label]);
    if (FZERO(gca_means[label])) {
      DiagBreak();
    }
  }

  width = gca->node_width;
  height = gca->node_height;
  depth = gca->node_depth;
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          label = gcan->labels[n];
          gc = &gcan->gcs[n];
          mean = gc->means[0] * means[label] / gca_means[label];
          gc->means[0] = mean;
        }
      }
    }
  }

  free(means);
  free(stds);
  free(gca_means);
  return (NO_ERROR);
}
int GCArenormalizeAdaptive(MRI *mri_in, MRI *mri_labeled, GCA *gca, TRANSFORM *transform, int wsize, float pthresh)
{
  int x, y, z, n, label, xp, yp, zp, peak, orig_wsize, frame;
  GCA_NODE *gcan;
  GCA_SAMPLE *gcas;
  GC1D *gc;
  HISTOGRAM *histo, *hsmooth;
  float fmin, fmax;
  int nsamples = 0;

  orig_wsize = wsize;
  MRIvalRange(mri_in, &fmin, &fmax);
  histo = HISTOalloc((int)(fmax - fmin + 1));
  hsmooth = HISTOalloc((int)(fmax - fmin + 1));

  /* go through GCA renormalizing each entry */
  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          label = gcan->labels[n];
          if (label == Unknown) {
            continue;
          }
          if (gcaNodeToPrior(gca, x, y, z, &xp, &yp, &zp) == NO_ERROR) {
            gc = &gcan->gcs[n];

#define MIN_SAMPLES 20
            wsize = orig_wsize;
            if (label == Ggca_label) {
              DiagBreak();
            }
            do {
              gcas = gcaExtractThresholdedRegionLabelAsSamples(
                  gca, mri_labeled, transform, &nsamples, label, xp, yp, zp, wsize, pthresh);

              if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
                DiagBreak();
              }
              wsize += 2;
              if (gcas && nsamples < MIN_SAMPLES) {
                GCAfreeSamples(&gcas, nsamples);
              }
            } while ((nsamples < MIN_SAMPLES) && (wsize < 2 * orig_wsize));

            if (nsamples < MIN_SAMPLES)
            /* couldn't find any in this nbhd */
            {
              continue;
            }

            for (frame = 0; frame < gca->ninputs; frame++) {
              gcaHistogramSamples(gca, gcas, mri_in, transform, nsamples, histo, frame);
              HISTOsmooth(histo, hsmooth, 2);
              if (IS_WM(label) && gca->ninputs == 1)
                peak = HISTOfindLastPeakRelative(hsmooth, HISTO_WINDOW_SIZE, .3);
              else if (IS_LAT_VENT(label) && gca->ninputs == 1)
                peak = HISTOfindFirstPeakRelative(hsmooth, HISTO_WINDOW_SIZE, .3);
              else
                peak = HISTOfindHighestPeakInRegion(hsmooth, 0, (fmax - fmin) * hsmooth->bin_size);
              if (peak < 0) {
                continue;
              }
              gc->means[frame] = (float)peak;
            }
            GCAfreeSamples(&gcas, nsamples);
          }
        }
      }
    }
  }

  HISTOfree(&histo);
  HISTOfree(&hsmooth);
  return (NO_ERROR);
}

#define MEAN_RESOLUTION 8

int GCArenormalizeLabels(MRI *mri_in, MRI *mri_labeled, GCA *gca, TRANSFORM *transform)
{
  int x, y, z, n, label, biggest_label, nsamples, xv, yv, zv, *ordered_indices, i, index, width, height, depth;
  float mean, var, val;
  GCA_NODE *gcan;
  GCA_SAMPLE *gcas;
  GC1D *gc;
  MRI *mri_means, *mri_control, *mri_tmp;
  char fname[STRLEN];

  if (gca->ninputs > 1) ErrorExit(ERROR_UNSUPPORTED, "GCArenormalizeLabels: can only renormalize scalars");

  /* first build a list of all labels that exist */
  for (biggest_label = x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_height; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          if (gcan->labels[n] > biggest_label) {
            biggest_label = gcan->labels[n];
          }
        }
      }
    }
  }

  for (label = 1; label <= biggest_label; label++) {
    gcas = gcaExtractLabelAsSamples(gca, mri_labeled, transform, &nsamples, label);
    if (!nsamples) {
      continue;
    }
    if (nsamples < MIN_SAMPLES / SAMPLE_PCT) {
      free(gcas);
      continue;
    }
    ordered_indices = (int *)calloc(nsamples, sizeof(int));
    GCAcomputeLogSampleProbability(gca, gcas, mri_in, transform, nsamples, DEFAULT_CLAMP);
    GCArankSamples(gca, gcas, nsamples, ordered_indices);

    if (nint(nsamples * SAMPLE_PCT) < MIN_SAMPLES) {
      nsamples = MIN_SAMPLES;
    }
    else {
      nsamples = nint(SAMPLE_PCT * (float)nsamples);
    }

    for (var = mean = 0.0f, i = 0; i < nsamples; i++) {
      index = ordered_indices[i];
      val = MRIgetVoxVal(mri_in, gcas[index].x, gcas[index].y, gcas[index].z, 0);
      mean += val;
      var += val * val;
    }
    mean /= (float)nsamples;
    var = var / nsamples - mean * mean;
    var = sqrt(var);
    printf("label %s: using %d samples to estimate mean = %2.1f +- %2.1f\n",
           cma_label_to_name(label),
           nsamples,
           mean,
           var);

    width = nint(mri_in->width / MEAN_RESOLUTION);
    height = nint(mri_in->height / MEAN_RESOLUTION);
    depth = nint(mri_in->depth / MEAN_RESOLUTION);
    mri_means = MRIalloc(width, height, depth, MRI_FLOAT);
    MRIcopyHeader(mri_in, mri_means);
    mri_control = MRIalloc(width, height, depth, MRI_SHORT);
    MRIcopyHeader(mri_in, mri_control);
    MRIsetResolution(mri_means, MEAN_RESOLUTION, MEAN_RESOLUTION, MEAN_RESOLUTION);
    MRIsetResolution(mri_control, MEAN_RESOLUTION, MEAN_RESOLUTION, MEAN_RESOLUTION);

    GCAlabelMri(gca, mri_means, label, transform);
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE) {
      sprintf(fname, "%s_label.mgz", cma_label_to_name(label));
      MRIwrite(mri_means, fname);
    }

    for (mean = 0.0f, n = z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (MRIFvox(mri_means, x, y, z) > 1) {
            n++;
            mean += MRIFvox(mri_means, x, y, z);
          }
        }
      }
    }
    mean /= (float)n;
    printf("mean GCA value %2.1f\n", mean);
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (MRIFvox(mri_means, x, y, z) < 1) {
            MRIFvox(mri_means, x, y, z) = mean;
          }
        }
      }
    }

    TransformInvert(transform, mri_in);
    for (i = 0; i < nsamples; i++) {
      index = ordered_indices[i];
      if (!GCApriorToSourceVoxel(
              gca, mri_means, transform, gcas[index].xp, gcas[index].yp, gcas[index].zp, &x, &y, &z)) {
        val = MRIgetVoxVal(mri_in, gcas[index].x, gcas[index].y, gcas[index].z, 0);
        if (x == 19 && y == 14 && z == 15) {
          DiagBreak();
        }
        if (MRISvox(mri_control, x, y, z) == 0) {
          MRIFvox(mri_means, x, y, z) = val;
        }
        else {
          MRIFvox(mri_means, x, y, z) += val;
        }
        MRISvox(mri_control, x, y, z)++;
      }
    }

    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == 19 && y == 14 && z == 15) {
            DiagBreak();
          }
          if (MRISvox(mri_control, x, y, z) > 0) {
            MRIFvox(mri_means, x, y, z) =
                nint((float)MRIFvox(mri_means, x, y, z) / (float)MRISvox(mri_control, x, y, z));
            MRISvox(mri_control, x, y, z) = 1;
          }
        }
      }
    }

    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE) {
      sprintf(fname, "%s_means.mgz", cma_label_to_name(label));
      MRIwrite(mri_means, fname);
    }

    mri_tmp = MRIalloc(width, height, depth, MRI_SHORT);
    MRIcopy(mri_means, mri_tmp);
    MRIfree(&mri_means);
    mri_means = mri_tmp;

    mri_tmp = MRIalloc(width, height, depth, MRI_UCHAR);
    MRIcopy(mri_control, mri_tmp);
    MRIfree(&mri_control);
    mri_control = mri_tmp;

    MRIsoapBubble(mri_means, mri_control, mri_means, 10, -1);
    MRIclear(mri_control);
    MRIsoapBubble(mri_means, mri_control, mri_means, 1, -1);

    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE) {
      sprintf(fname, "%s_soap.mgz", cma_label_to_name(label));
      MRIwrite(mri_means, fname);

      sprintf(fname, "%s_control.mgz", cma_label_to_name(label));
      MRIwrite(mri_control, fname);
    }

    width = gca->node_width;
    height = gca->node_height;
    depth = gca->node_depth;
    for (z = 0; z < depth; z++) {
      for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
          if (x == Gx && y == Gy && z == Gz) {
            DiagBreak();
          }
          gcan = &gca->nodes[x][y][z];
          for (n = 0; n < gcan->nlabels; n++) {
            if (gcan->labels[n] == label) {
              if (x == Gx && y == Gy && z == Gz) {
                DiagBreak();
              }
              gc = &gcan->gcs[n];
              if (GCAnodeToVoxel(gca, mri_means, x, y, z, &xv, &yv, &zv) == NO_ERROR)
                gc->means[0] = (float)MRISvox(mri_means, xv, yv, zv);
              break;
            }
          }
        }
      }
    }

    MRIfree(&mri_means);
    MRIfree(&mri_control);
    free(gcas);
    free(ordered_indices);
  }

  return (NO_ERROR);
}


/*
  labels and intensities have pairs of modes to update
  the gca with.
*/
int GCArenormalizeIntensities(GCA *gca, int *labels, float *intensities, int num)
{
  float *scales, mode;
  int xn, yn, zn, n, i, label;
  GCA_NODE *gcan;
  GC1D *gc;

  if (gca->ninputs > 1) ErrorExit(ERROR_UNSUPPORTED, "GCArenormalizeIntensities: can only renormalize scalars");

  scales = (float *)calloc(num, sizeof(float));
  if (scales == NULL) ErrorExit(ERROR_NOMEMORY, "GCArenormalizeIntensities: could not alloc %d elt array");

  for (i = 0; i < num; i++) {
    if (intensities[i] <= 0) {
      scales[i] = 1;
      continue;
    }
    GCAlabelMode(gca, labels[i], &mode);
    if (mode <= 0) GCAlabelMean(gca, Left_Lateral_Ventricle, &mode);

    if (mode <= 0)
      scales[i] = 1.0;
    else {
      scales[i] = intensities[i] / mode;
      printf("rescaling %s from %2.0f --> %2.0f\n", cma_label_to_name(labels[i]), mode, mode * scales[i]);
    }
  }

  /* now go through each label and scale it by these ratios */
  for (zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
          DiagBreak();
        }
        for (n = 0; n < gcan->nlabels; n++) {
          label = gcan->labels[n];
          if (label == Gdiag_no) {
            DiagBreak();
          }

          /* find index in lookup table for this label */
          for (i = 0; i < num; i++)
            if (label == labels[i]) {
              break;
            }
          if (i >= num) {
            continue;
          }

          gc = &gcan->gcs[n];
          if (gc->means[0] < 1 && gcan->labels[n] == 0)  // too close to 0 to scale
          {
            gc->means[0] = intensities[i] / 3;
          }
          else
	  {
            if ((xn == Gx && yn == Gy && zn == Gz) &&
                (Ggca_label == gcan->labels[i] || Ggca_label < 0))
	      printf("scaling gc for %s at (%d %d %d) from %2.1f --> %2.1f\n",
		     cma_label_to_name(label), xn, yn, zn, gc->means[0], scales[i] * gc->means[0]);
            gc->means[0] = scales[i] * gc->means[0];
	  }
	}
      }
    }
  }

  free(scales);
  return (NO_ERROR);
}
int GCAunifyVariance(GCA *gca)
{
  int xn, yn, zn, n, r, c, v;
  GCA_NODE *gcan;
  GC1D *gc;

  for (zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          gc = &gcan->gcs[n];
          for (r = v = 0; r < gca->ninputs; r++) {
            for (c = r; c < gca->ninputs; c++, v++) {
              if (r == c) {
                gc->covars[v] = MIN_VAR;
              }
              else {
                gc->covars[v] = 0.0;
              }
            }
          }
        }
      }
    }
  }

  return (NO_ERROR);
}

int GCAlabelMode(GCA *gca, int label, float *modes)
{
  int xn, yn, zn, n, r;
  GCA_NODE *gcan;
  GC1D *gc;
  float prior;
  HISTOGRAM *h;
  int b;

  h = HISTOalloc(256);
  for (b = 0; b < h->nbins; b++) {
    h->bins[b] = b;
  }

  memset(modes, 0, gca->ninputs * sizeof(float));
  for (zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          /* find index in lookup table for this label */
          if (gcan->labels[n] != label) {
            continue;
          }
          gc = &gcan->gcs[n];
          prior = get_node_prior(gca, label, xn, yn, zn);
          if (prior != 0) {
            for (r = 0; r < gca->ninputs; r++) {
              b = nint(gc->means[r]);
              h->counts[b] += prior;
              if (!std::isfinite(gc->means[r])) {
                DiagBreak();
              }
            }
          }
        }
      }
    }
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "gca_label%d.plt", label);
    HISTOplot(h, fname);
  }

  for (r = 0; r < gca->ninputs; r++) {
    b = HISTOfindHighestPeakInRegion(h, 0, h->nbins);
    modes[r] = h->bins[b];
  }
  return (NO_ERROR);
}
int GCAclassMode(GCA *gca, int classnum, float *modes)
{
  int xn, yn, zn, n, r;
  GCA_NODE *gcan;
  GC1D *gc;
  float prior;
  HISTOGRAM *h;
  int b, label;

  h = HISTOalloc(256);
  for (b = 0; b < h->nbins; b++) {
    h->bins[b] = b;
  }

  memset(modes, 0, gca->ninputs * sizeof(float));
  for (zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          label = gcan->labels[n];
          switch (classnum)  // check to make sure it is
                          // the specified class
          {
            case WM_CLASS:
              if (IS_WHITE_CLASS(gcan->labels[n]) == 0) {
                continue;
              }
              break;
            case GM_CLASS:
              if (IS_GRAY_CLASS(gcan->labels[n]) == 0) {
                continue;
              }
              break;
            case CSF_CLASS:
              if (IS_CSF_CLASS(gcan->labels[n]) == 0) {
                continue;
              }
              break;
            default:
              break;
          }
          prior = get_node_prior(gca, label, xn, yn, zn);
          gc = GCAfindGC(gca, xn, yn, zn, label);
          if (gc == NULL) {
            continue;
          }
          if (prior != 0) {
            for (r = 0; r < gca->ninputs; r++) {
              b = nint(gc->means[r]);
              h->counts[b] += prior;
              if (!std::isfinite(gc->means[r])) {
                DiagBreak();
              }
            }
          }
        }
      }
    }
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "gca_label%d.plt", classnum);
    HISTOplot(h, fname);
  }

  for (r = 0; r < gca->ninputs; r++) {
    b = HISTOfindHighestPeakInRegion(h, 0, h->nbins);
    modes[r] = h->bins[b];
  }
  return (NO_ERROR);
}

int GCAlabelMean(GCA *gca, int label, float *means)
{
  int xn, yn, zn, n, r;
  GCA_NODE *gcan;
  GC1D *gc;
  double wt;
  float prior;

  /* compute overall white matter mean to use as anchor for rescaling */
  memset(means, 0, gca->ninputs * sizeof(float));
  for (wt = 0.0, zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          /* find index in lookup table for this label */
          if (gcan->labels[n] != label) {
            continue;
          }
          gc = &gcan->gcs[n];
          prior = get_node_prior(gca, label, xn, yn, zn);
          if (prior != 0) {
            wt += prior;
            for (r = 0; r < gca->ninputs; r++) {
              means[r] += gc->means[r] * prior;
              if (!std::isfinite(gc->means[r])) {
                DiagBreak();
              }
            }
          }
        }
      }
    }
  }
  if (FZERO(wt)) {
    return (NO_ERROR);
  }
  for (r = 0; r < gca->ninputs; r++) {
    means[r] /= wt;
  }
  return (NO_ERROR);
}
int GCAlabelVar(GCA *gca, int label, float *vars)
{
  int xn, yn, zn, n, c, r, v;
  GCA_NODE *gcan;
  GC1D *gc;
  double wt;
  float prior;

  /* compute overall white matter mean to use as anchor for rescaling */
  memset(vars, 0, gca->ninputs * sizeof(float));
  for (wt = 0.0, zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          /* find index in lookup table for this label */
          if (gcan->labels[n] != label) {
            continue;
          }
          gc = &gcan->gcs[n];
          prior = get_node_prior(gca, label, xn, yn, zn);
          if (prior != 0) {
            wt += prior;
            for (r = v = 0; r < gca->ninputs; r++)
              for (c = r; c < gca->ninputs; c++, v++) {
                vars[r] += gc->covars[v] * prior;
                if (!std::isfinite(vars[r])) {
                  DiagBreak();
                }
              }
          }
        }
      }
    }
  }
  if (FZERO(wt)) {
    return (NO_ERROR);
  }
  for (r = 0; r < gca->ninputs; r++) {
    vars[r] /= wt;
  }
  return (NO_ERROR);
}
int GCAclassMean(GCA *gca, int classnum, float *means)
{
  int xn, yn, zn, n, r, label;
  GCA_NODE *gcan;
  GC1D *gc;
  double wt;
  float prior;

  /* compute overall white matter mean to use as anchor for rescaling */
  memset(means, 0, gca->ninputs * sizeof(float));
  for (wt = 0.0, zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          label = gcan->labels[n];
          switch (classnum)  // check to make sure it is
                          // the specified class
          {
            case WM_CLASS:
              if (IS_WHITE_CLASS(gcan->labels[n]) == 0) {
                continue;
              }
              break;
            case GM_CLASS:
              if (IS_GRAY_CLASS(gcan->labels[n]) == 0) {
                continue;
              }
              break;
            case CSF_CLASS:
              if (IS_CSF_CLASS(gcan->labels[n]) == 0) {
                continue;
              }
              break;
            default:
              break;
          }
          /* find index in lookup table for this label */
          gc = &gcan->gcs[n];
          prior = get_node_prior(gca, label, xn, yn, zn);
          if (prior != 0) {
            wt += prior;
            for (r = 0; r < gca->ninputs; r++) {
              means[r] += gc->means[r] * prior;
              if (!std::isfinite(gc->means[r])) {
                DiagBreak();
              }
            }
          }
        }
      }
    }
  }
  if (FZERO(wt)) {
    return (NO_ERROR);
  }
  for (r = 0; r < gca->ninputs; r++) {
    means[r] /= wt;
  }
  return (NO_ERROR);
}

int GCAregularizeConditionalDensities(GCA *gca, float smooth)
{
  int xn, yn, zn, n, i, label, max_label, r;
  GCA_NODE *gcan;
  GC1D *gc;
  double **means, *wts;
  float prior;

  means = (double **)calloc(gca->ninputs, sizeof(double *));
  wts = (double *)calloc(MAX_GCA_LABELS, sizeof(double));
  if (!means || !wts)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d wt and mean vectors", Progname, MAX_GCA_LABELS);
  for (r = 0; r < gca->ninputs; r++) {
    means[r] = (double *)calloc(MAX_GCA_LABELS, sizeof(double));
    if (!means[r]) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d wt and mean vectors", Progname, MAX_GCA_LABELS);
  }

  /* compute overall mean for each class */
  for (max_label = 1, zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          /* find index in lookup table for this label */
          label = gcan->labels[n];

          gc = &gcan->gcs[n];
          prior = get_node_prior(gca, label, xn, yn, zn);
          if (prior != 0) {
            wts[label] += prior;
            for (r = 0; r < gca->ninputs; r++) {
              means[r][label] += gc->means[r] * prior;
            }
            if (label > max_label) {
              max_label = label;
            }
          }
        }
      }
    }
  }

  for (i = 0; i <= max_label; i++) {
    if (!FZERO(wts[i])) {
      for (r = 0; r < gca->ninputs; r++) {
        means[r][i] /= wts[i];
      }
    }
  }

  /* now impose regularization */
  for (zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          /* find index in lookup table for this label */
          label = gcan->labels[n];
          if (label <= 0) {
            continue;
          }

          gc = &gcan->gcs[n];
          for (r = 0; r < gca->ninputs; r++) gc->means[r] = means[r][label] * smooth + gc->means[r] * (1.0f - smooth);
        }
      }
    }
  }

  for (r = 0; r < gca->ninputs; r++) {
    free(means[r]);
  }
  free(wts);
  free(means);
  return (NO_ERROR);
}
int GCArenormalizeToFlash(GCA *gca, char *tissue_parms_fname, MRI *mri)
{
  FILE *fp;
  char *cp, line[STRLEN];
  int labels[MAX_GCA_LABELS], nlabels;
  float intensities[MAX_GCA_LABELS], TR, alpha, T1, PD;

  TR = mri->tr;
  alpha = mri->flip_angle;

  fp = fopen(tissue_parms_fname, "r");
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE, "GCArenormalizeToFlash: could not open tissue parms file %s", tissue_parms_fname));

  cp = fgetl(line, STRLEN - 1, fp);
  nlabels = 0;
  while (cp) {
    if (sscanf(cp, "%d %f %f", &labels[nlabels], &T1, &PD) != 3)
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "GCArenormalizeToFlash: "
                   "could not parse %dth line %s in %s",
                   nlabels + 1,
                   cp,
                   tissue_parms_fname));

    intensities[nlabels] = FLASHforwardModel(T1, PD, TR, alpha, 3);

    nlabels++;
    cp = fgetl(line, STRLEN - 1, fp);
  }
  fclose(fp);
  GCArenormalizeIntensities(gca, labels, intensities, nlabels);
  return (NO_ERROR);
}

int GCAmeanFilterConditionalDensities(GCA *gca, float navgs)
{
  /* won't work for covariances */
  int xn, yn, zn, xn1, yn1, zn1, n, label, max_label, niter, f;
  GCA_NODE *gcan;
  GC1D *gc;
  double means[MAX_GCA_INPUTS], wt;
  MRI *mri_means;
  float prior;

  mri_means = MRIallocSequence(gca->node_width, gca->node_height, gca->node_depth, MRI_FLOAT, gca->ninputs);

  mri_means->xsize = gca->node_spacing;
  mri_means->ysize = gca->node_spacing;
  mri_means->zsize = gca->node_spacing;

  GCAcopyDCToMRI(gca, mri_means);

  /* compute overall mean for each class */
  for (max_label = 1, zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        if (xn == Gx && yn == Gy && zn == Gz) {
          DiagBreak();
        }
        for (n = 0; n < gcan->nlabels; n++) {
          /* find index in lookup table for this label */
          label = gcan->labels[n];

          gc = &gcan->gcs[n];
          if (label > max_label) {
            max_label = label;
          }
        }
      }
    }
  }

  /* now impose regularization */
  for (niter = 0; niter < navgs; niter++) {
    for (label = 0; label <= max_label; label++) {
      if (label == Gdiag_no) {
        DiagBreak();
      }
      if (IS_UNKNOWN(label)) {
        continue;
      }
      for (zn = 0; zn < gca->node_depth; zn++) {
        for (yn = 0; yn < gca->node_height; yn++) {
          for (xn = 0; xn < gca->node_width; xn++) {
            if (xn == Gx && yn == Gy && zn == Gz) {
              DiagBreak();
            }
            wt = 0.0;
            for (f = 0; f < gca->ninputs; f++) {
              means[f] = 0;
            }
            for (xn1 = xn - 1; xn1 <= xn + 1; xn1++) {
              if (xn1 < 0 || xn1 >= gca->node_width) {
                continue;
              }
              for (yn1 = yn - 1; yn1 <= yn + 1; yn1++) {
                if (yn1 < 0 || yn1 >= gca->node_height) {
                  continue;
                }
                for (zn1 = zn - 1; zn1 <= zn + 1; zn1++) {
                  if (zn1 < 0 || zn1 >= gca->node_depth) {
                    continue;
                  }
                  if (xn1 == xn && yn1 == yn && zn1 == zn) {
                    DiagBreak();
                  }
                  gcan = &gca->nodes[xn1][yn1][zn1];
                  for (n = 0; n < gcan->nlabels; n++) {
                    /* find index in lookup table
                       for this label */
                    if (gcan->labels[n] != label) {
                      continue;
                    }

                    gc = &gcan->gcs[n];
                    prior = get_node_prior(gca, label, xn1, yn1, zn1);
                    if (prior != 0) {
                      for (f = 0; f < gca->ninputs; f++) {
                        means[f] += prior * gc->means[f];
                        if (!std::isfinite(means[f] / wt)) {
                          DiagBreak();
                        }
                      }
                      wt += prior;
                      break;
                    }
                  }
                }
              }
            }
            if (FZERO(wt)) {
              continue; /* label didn't occur here */
            }
            for (f = 0; f < gca->ninputs; f++) {
              if (!std::isfinite(means[f] / wt)) {
                DiagBreak();
              }
              MRIFseq_vox(mri_means, xn, yn, zn, f) = means[f] / wt;
              if (!std::isfinite(MRIFseq_vox(mri_means, xn, yn, zn, f))) {
                DiagBreak();
              }
            }
          }
        }
      }

      for (zn = 0; zn < gca->node_depth; zn++) {
        for (yn = 0; yn < gca->node_height; yn++) {
          for (xn = 0; xn < gca->node_width; xn++) {
            gcan = &gca->nodes[xn][yn][zn];
            for (n = 0; n < gcan->nlabels; n++) {
              if (gcan->labels[n] != label) {
                continue;
              }
              gc = &gcan->gcs[n];
              if (xn == Gx && yn == Gy && zn == Gz) {
                DiagBreak();
              }
              for (f = 0; f < gca->ninputs; f++) {
                gc->means[f] = MRIFseq_vox(mri_means, xn, yn, zn, f);
                if (!std::isfinite(gc->means[f])) {
                  DiagBreak();
                }
              }
              break;
            }
          }
        }
      }
    }
  }

  MRIfree(&mri_means);
  return (NO_ERROR);
}

#define OFFSET_SIZE 25
double BOX_SIZE = 60; /* mm */
double HALF_BOX = (60 / 2);
int GCAhistoScaleImageIntensities(GCA *gca, MRI *mri, int noskull)
{
  float x0, y0, z0, fmin, fmax, min_real_val;
  int mri_peak, r, max_T1_weighted_image = 0, min_real_bin, peak, label, left_wm, right_wm;
  float wm_means[MAX_GCA_INPUTS], tmp[MAX_GCA_INPUTS], scales[MAX_GCA_INPUTS];
  // float max_wm /*, scale*/;
  HISTOGRAM *h_mri, *h_smooth, *h_gca;
  MRI_REGION box;
  MRI *mri_frame, *mri_mask, *mri_tmp, *mri_seg;
  float gm_means[MAX_GCA_INPUTS], gray_white_CNR;

  if ((gca->flags & GCA_NO_LH) == 0) GCAlabelMean(gca, Left_Cerebral_White_Matter, wm_means);
  // GCAlabelMean(gca, Left_Cerebral_White_Matter, tmp) ;
  if ((gca->flags & GCA_NO_RH) == 0) GCAlabelMean(gca, Right_Cerebral_White_Matter, tmp);

// the following rule some times fails for BIRN data, so changed it
// on 08-01-05 by xhan
  if ((gca->flags & GCA_NO_LH) == 0) GCAlabelMean(gca, Left_Cerebral_Cortex, gm_means);
  for (r = 0; r < gca->ninputs; r++) {
    // use modes instead of means
    if ((gca->flags & GCA_NO_LH) == 0)
      h_gca = gcaGetLabelHistogram(gca, Left_Cerebral_White_Matter, r, 0);
    else
      h_gca = gcaGetLabelHistogram(gca, Right_Cerebral_White_Matter, r, 0);
    peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
    printf("resetting wm mean[%d]: %2.0f --> %2.0f\n", r, wm_means[r], h_gca->bins[peak]);
    wm_means[r] = h_gca->bins[peak];
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      HISTOplot(h_gca, "wm.plt");
    }
    HISTOfree(&h_gca);

    if ((gca->flags & GCA_NO_LH) == 0)
      h_gca = gcaGetLabelHistogram(gca, Left_Cerebral_Cortex, r, 0);
    else
      h_gca = gcaGetLabelHistogram(gca, Right_Cerebral_Cortex, r, 0);
    peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
    gm_means[r] = h_gca->bins[peak];
    printf("resetting gm mean[%d]: %2.0f --> %2.0f\n", r, gm_means[r], h_gca->bins[peak]);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      HISTOplot(h_gca, "gm.plt");
    }
    HISTOfree(&h_gca);
  }
  gray_white_CNR = wm_means[0] - gm_means[0];
  for (r = 0; r < gca->ninputs; r++) {
    //                wm_means[r] = (wm_means[r] + tmp[r]) / 2 ;
    if ((wm_means[r] - gm_means[r]) > gray_white_CNR) {
      max_T1_weighted_image = r;
      gray_white_CNR = (wm_means[r] - gm_means[r]);
    }
  }
  printf("input volume #%d is the most T1-like\n", max_T1_weighted_image + 1);

  // max_wm = wm_means[max_T1_weighted_image];

  mri_frame = MRIcopyFrame(mri, NULL, max_T1_weighted_image, 0);
  mri_tmp = MRImean(mri_frame, NULL, 5);
  MRIfree(&mri_frame);
  mri_frame = mri_tmp;
  MRIvalRange(mri_frame, &fmin, &fmax);
  h_mri = MRIhistogram(mri_frame, nint(fmax - fmin + 1));
  HISTOclearZeroBin(h_mri); /* ignore background */
  h_smooth = HISTOsmooth(h_mri, NULL, 2);
  mri_peak = HISTOfindFirstPeak(h_smooth, 5, .1);
  min_real_bin = HISTOfindEndOfPeak(h_smooth, mri_peak, .25);
#define MAX_PEAK_BINS 40
  if (min_real_bin - mri_peak > MAX_PEAK_BINS) {
    min_real_bin = HISTOfindNextValley(h_smooth, mri_peak);
  }
  if (min_real_bin - mri_peak > MAX_PEAK_BINS) {
    min_real_bin = mri_peak + MAX_PEAK_BINS - 1;
  }
  min_real_val = h_smooth->bins[min_real_bin];
  if (noskull) {
    if (min_real_val > 0.25 * fmax) /* for skull-stripped images */
    {
      min_real_val = 0.25 * fmax;
    }
  }
  else {
    if (min_real_val > 0.5 * fmax) /* for non skull-stripped images */
    {
      min_real_val = 0.5 * fmax;
    }
  }
  printf("using real data threshold=%2.1f\n", min_real_val);

  MRIfindApproximateSkullBoundingBox(mri_frame, min_real_val, &box);
  printf("skull bounding box = (%d, %d, %d) --> (%d, %d, %d)\n",
         box.x,
         box.y,
         box.z,
         box.x + box.dx - 1,
         box.y + box.dy - 1,
         box.z + box.dz - 1);
  mri_mask = MRIbinarize(mri_frame, NULL, min_real_val, 0, 1);
  MRIfree(&mri_frame);
  HISTOfree(&h_mri);
  HISTOfree(&h_smooth);

  // why divided by 3 for x and y?? mistake or experience?? -xh
  // experience - don't want to be near midline (BRF)
  mri_seg = GCAbuildMostLikelyLabelVolume(gca, NULL) ;
  left_wm = MRIvoxelsInLabel(mri_seg, Left_Cerebral_White_Matter) ;
  right_wm = MRIvoxelsInLabel(mri_seg, Right_Cerebral_White_Matter) ;
  label = (left_wm > 2*right_wm) ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
  printf("finding center of left hemi white matter\n") ;
  MRIfree(&mri_seg) ;

  y0 = box.y + box.dy / 3;
  z0 = box.z + box.dz / 2;
  if (label == Left_Cerebral_White_Matter)
    x0 = box.x + 2*box.dx / 3;
  else
    x0 = box.x + box.dx / 3;
  printf("using (%.0f, %.0f, %.0f) as brain centroid of %s...\n", x0, y0, z0, cma_label_to_name(label));
  box.dx /= 4;
  box.x = x0 - box.dx / 2;
  box.dy /= 4;
  box.y = y0 - box.dy / 2;
  box.dz /= 4;
  box.z = z0 - box.dz / 2;
  for (r = 0; r < gca->ninputs; r++) {
    mri_frame = MRIcopyFrame(mri, NULL, r, 0);
    MRImask(mri_frame, mri_mask, mri_frame, 0, 0); /* remove stuff that is
                                                      background or csf */

    printf(
        "mean wm in atlas = %2.0f, using box (%d,%d,%d) --> (%d, %d,%d) "
        "to find MRI wm\n",
        wm_means[r],
        box.x,
        box.y,
        box.z,
        box.x + box.dx - 1,
        box.y + box.dy - 1,
        box.z + box.dz - 1);

    h_mri = MRIhistogramRegion(mri_frame, 0, NULL, &box);
    if (h_mri == NULL) return (Gerror);
    if (gca->ninputs == 1) {
      HISTOclearBins(h_mri, h_mri, 0, min_real_val);
    }
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      HISTOplot(h_mri, "mri.histo");
    }
    HISTOclearZeroBin(h_mri);
    if (gca->ninputs == 1 || r == max_T1_weighted_image) /* assume it is T1-weighted */
    {
#define MIN_CONFORMED_WM_VAL 50   // assume wm greater than this
#define MAX_CONFORMED_WM_VAL 240  // assume wm than this
      if (mriConformed(mri))      // use strong priors on where wm should be
      {
        HISTOclearBins(h_mri, h_mri, 0, MIN_CONFORMED_WM_VAL);
        HISTOclearBins(h_mri, h_mri, MAX_CONFORMED_WM_VAL + 1, 255);
      }
      mri_peak = HISTOfindLastPeak(h_mri, 2 * HISTO_WINDOW_SIZE, MIN_HISTO_PCT);
    }
    else {
      mri_peak = HISTOfindHighestPeakInRegion(h_mri, 1, h_mri->nbins);
    }
    if ((h_mri->nbins <= mri_peak) || (mri_peak < 0)) {
      printf(
          "WARNING: gca.c::GCAhistoScaleImageIntensities: "
          "h_mri->nbins=%d, mri_peak=%d\n",
          h_mri->nbins,
          mri_peak);
      mri_peak = 0;
    }
    else {
      mri_peak = h_mri->bins[mri_peak];
    }
    printf("before smoothing, mri peak at %d\n", mri_peak);
    h_smooth = HISTOsmooth(h_mri, NULL, 2);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      HISTOplot(h_smooth, "mri_smooth.histo");
    }
    /* assume it is the right-most peak of image
       is supposed to be T1-weighted */
    if (gca->ninputs == 1 && (gca->type == GCA_UNKNOWN || gca->type == GCA_NORMAL ||
                              (gca->type == GCA_FLASH && (DEGREES(mri->flip_angle) > 15))))
      mri_peak = HISTOfindLastPeak(h_smooth, HISTO_WINDOW_SIZE, MIN_HISTO_PCT);
    else {
      mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 1, h_mri->nbins);
    }
    {
      double mn, std, std_thresh = 10;
      if (getenv("FS_HISTO_STD_THRESH"))
      {
	char *cp = getenv("FS_HISTO_STD_THRESH") ;
	std_thresh = atof(cp) ;
	printf("FS_HISTO_STD_THRESH found in the environment resetting from 10 to %2.1f\n", std_thresh) ;
      }
      if (mri->xsize < .9) std_thresh *= 2;
      HISTOrobustGaussianFit(h_smooth, .3, &mn, &std);
      printf("robust fit to distribution - %2.0f +- %2.1f\n", mn, std);
      std_thresh *= (mn / wm_means[r]);
      if (std > std_thresh) {
        printf("distribution too broad for accurate scaling - disabling\n");
        mri_peak = HISTOfindBin(h_smooth, wm_means[r]);
      }
    }

    if ((h_mri->nbins <= mri_peak) || (mri_peak < 0)) {
      printf(
          "WARNING2: gca.c::GCAhistoScaleImageIntensities: "
          "h_mri->nbins=%d, mri_peak=%d\n",
          h_mri->nbins,
          mri_peak);
      mri_peak = 0;
    }
    else {
      mri_peak = h_smooth->bins[mri_peak];
    }

    printf(
        "after smoothing, mri peak at %d, scaling input intensities "
        "by %2.3f\n",
        mri_peak,
        wm_means[r] / mri_peak);
    if (mri_peak == 0) ErrorExit(ERROR_BADPARM, "GCAhistoScaleImageIntensities: could not find wm peak");
    scales[r] = wm_means[r] / mri_peak;

    MRIfree(&mri_frame);
    HISTOfree(&h_mri);
    HISTOfree(&h_smooth);
  }

  // scale each frame independently -xhan
  for (r = 0; r < gca->ninputs; r++) {
    printf("scaling channel %d by %g\n", r, scales[r]);
    MRIscalarMulFrame(mri, mri, scales[r], r);
  }
  MRIfree(&mri_mask);
  return (NO_ERROR);
}

/*
  each input frame is a different time point for the same subject.
*/
int GCAhistoScaleImageIntensitiesLongitudinal(GCA *gca, MRI *mri, int noskull)
{
  float x0, y0, z0, fmin, fmax, min_real_val, val;
  int mri_peak, r, max_T1_weighted_image = 0, min_real_bin, peak, frame;
  float wm_means[MAX_GCA_INPUTS], tmp[MAX_GCA_INPUTS], scales[MAX_GCA_INPUTS];
  // float max_wm /*, scale*/;
  HISTOGRAM *h_mri, *h_smooth, *h_gca;
  MRI_REGION box, tbox;
  MRI *mri_frame, *mri_mask, *mri_tmp;

  float gm_means[MAX_GCA_INPUTS], gray_white_CNR;

  GCAlabelMean(gca, Left_Cerebral_White_Matter, wm_means);
  GCAlabelMean(gca, Right_Cerebral_White_Matter, tmp);

  GCAlabelMean(gca, Left_Cerebral_Cortex, gm_means);
  for (r = 0; r < gca->ninputs; r++) {
    // use modes instead of means
    h_gca = gcaGetLabelHistogram(gca, Left_Cerebral_White_Matter, r, 0);
    peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
    printf("resetting wm mean[%d]: %2.0f --> %2.0f\n", r, wm_means[r], h_gca->bins[peak]);
    wm_means[r] = h_gca->bins[peak];
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      HISTOplot(h_gca, "wm.plt");
    }
    HISTOfree(&h_gca);

    h_gca = gcaGetLabelHistogram(gca, Left_Cerebral_Cortex, r, 0);
    peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
    gm_means[r] = h_gca->bins[peak];
    printf("resetting gm mean[%d]: %2.0f --> %2.0f\n", r, gm_means[r], h_gca->bins[peak]);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      HISTOplot(h_gca, "gm.plt");
    }
    HISTOfree(&h_gca);
  }
  gray_white_CNR = wm_means[0] - gm_means[0];
  for (r = 0; r < gca->ninputs; r++) {
    //                wm_means[r] = (wm_means[r] + tmp[r]) / 2 ;
    if ((wm_means[r] - gm_means[r]) > gray_white_CNR) {
      max_T1_weighted_image = r;
      gray_white_CNR = (wm_means[r] - gm_means[r]);
    }
  }
  printf("input volume #%d is the most T1-like\n", max_T1_weighted_image + 1);

  // max_wm = wm_means[max_T1_weighted_image];

  min_real_val = 0;
  box.x = box.y = box.z = box.dx = box.dy = box.dz = 0;
  for (frame = 0; frame < mri->nframes; frame++) {
    mri_tmp = MRIcopyFrame(mri, NULL, frame, 0);
    mri_frame = MRImean(mri_tmp, NULL, 5);
    MRIvalRange(mri_frame, &fmin, &fmax);
    h_mri = MRIhistogram(mri_frame, nint(fmax - fmin + 1));
    HISTOclearZeroBin(h_mri); /* ignore background */
    h_smooth = HISTOsmooth(h_mri, NULL, 2);
    mri_peak = HISTOfindFirstPeak(h_smooth, 5, .1);
    min_real_bin = HISTOfindEndOfPeak(h_smooth, mri_peak, .25);
#define MAX_PEAK_BINS 40
    if (min_real_bin - mri_peak > MAX_PEAK_BINS) {
      min_real_bin = HISTOfindNextValley(h_smooth, mri_peak);
    }
    if (min_real_bin - mri_peak > MAX_PEAK_BINS) {
      min_real_bin = mri_peak + MAX_PEAK_BINS - 1;
    }
    val = h_smooth->bins[min_real_bin];
    if (noskull) {
      if (val > 0.25 * fmax) /* for skull-stripped images */
      {
        val = 0.25 * fmax;
      }
    }
    else {
      if (val > 0.5 * fmax) /* for non skull-stripped images */
      {
        val = 0.5 * fmax;
      }
    }
    min_real_val += val;

    printf("using real data threshold=%2.1f\n", val);
    MRIfindApproximateSkullBoundingBox(mri_frame, val, &tbox);
    printf("skull bounding box = (%d, %d, %d) --> (%d, %d, %d)\n",
           tbox.x,
           tbox.y,
           tbox.z,
           tbox.x + tbox.dx - 1,
           tbox.y + tbox.dy - 1,
           tbox.z + tbox.dz - 1);
    box.x += tbox.x;
    box.y += tbox.y;
    box.z += tbox.z;
    box.dx += tbox.dx;
    box.dy += tbox.dy;
    box.dz += tbox.dz;
    MRIfree(&mri_frame);
    HISTOfree(&h_mri);
    HISTOfree(&h_smooth);
  }

  min_real_val /= mri->nframes;
  box.x /= mri->nframes;
  box.y /= mri->nframes;
  box.z /= mri->nframes;
  box.dx /= mri->nframes;
  box.dy /= mri->nframes;
  box.dz /= mri->nframes;
  // divide by 3 to avoid midline
  x0 = box.x + box.dx / 3;
  y0 = box.y + box.dy / 3;
  z0 = box.z + box.dz / 2;
  printf("using (%.0f, %.0f, %.0f) as brain centroid...\n", x0, y0, z0);
  box.dx /= 4;
  box.x = x0 - box.dx / 2;
  box.dy /= 4;
  box.y = y0 - box.dy / 2;
  box.dz /= 4;
  box.z = z0 - box.dz / 2;

  mri_tmp = MRIcopyFrame(mri, NULL, 0, 0);
  mri_frame = MRImean(mri_tmp, NULL, 5);
  mri_mask = MRIbinarize(mri_frame, NULL, min_real_val, 0, 1);
  MRIfree(&mri_frame);
  MRIfree(&mri_tmp);
  for (frame = 0; frame < mri->nframes; frame++) {
    mri_frame = MRIcopyFrame(mri, NULL, frame, 0);
    MRImask(mri_frame, mri_mask, mri_frame, 0, 0); /* remove stuff that is
                                                      background or csf */

    printf(
        "mean wm in atlas = %2.0f, using box (%d,%d,%d) --> (%d, %d,%d) "
        "to find MRI wm\n",
        wm_means[max_T1_weighted_image],
        box.x,
        box.y,
        box.z,
        box.x + box.dx - 1,
        box.y + box.dy - 1,
        box.z + box.dz - 1);

    h_mri = MRIhistogramRegion(mri_frame, 0, NULL, &box);
    if (gca->ninputs == 1) {
      HISTOclearBins(h_mri, h_mri, 0, min_real_val);
    }
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      HISTOplot(h_mri, "mri.histo");
    }
    HISTOclearZeroBin(h_mri);
    if (gca->ninputs == 1) /* assume it is T1-weighted */
    {
#define MIN_CONFORMED_WM_VAL 50   // assume wm greater than this
#define MAX_CONFORMED_WM_VAL 240  // assume wm than this
      if (mriConformed(mri))      // use strong priors on where wm should be
      {
        HISTOclearBins(h_mri, h_mri, 0, MIN_CONFORMED_WM_VAL);
        HISTOclearBins(h_mri, h_mri, MAX_CONFORMED_WM_VAL + 1, 255);
      }
      mri_peak = HISTOfindLastPeak(h_mri, 2 * HISTO_WINDOW_SIZE, MIN_HISTO_PCT);
    }
    else {
      mri_peak = HISTOfindHighestPeakInRegion(h_mri, 1, h_mri->nbins);
    }
    if ((h_mri->nbins <= mri_peak) || (mri_peak < 0)) {
      printf(
          "WARNING: gca.c::GCAhistoScaleImageIntensities: "
          "h_mri->nbins=%d, mri_peak=%d\n",
          h_mri->nbins,
          mri_peak);
      mri_peak = 0;
    }
    else {
      mri_peak = h_mri->bins[mri_peak];
    }
    printf("before smoothing, mri peak at %d\n", mri_peak);
    h_smooth = HISTOsmooth(h_mri, NULL, 2);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      HISTOplot(h_smooth, "mri_smooth.histo");
    }
    /* assume it is the right-most peak of image
       is supposed to be T1-weighted */
    if (gca->ninputs == 1 && (gca->type == GCA_UNKNOWN || gca->type == GCA_NORMAL ||
                              (gca->type == GCA_FLASH && (DEGREES(mri->flip_angle) > 15))))
      mri_peak = HISTOfindLastPeak(h_smooth, HISTO_WINDOW_SIZE, MIN_HISTO_PCT);
    else {
      mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 1, h_mri->nbins);
    }
    if ((h_mri->nbins <= mri_peak) || (mri_peak < 0)) {
      printf(
          "WARNING2: gca.c::GCAhistoScaleImageIntensities: "
          "h_mri->nbins=%d, mri_peak=%d\n",
          h_mri->nbins,
          mri_peak);
      mri_peak = 0;
    }
    else {
      mri_peak = h_smooth->bins[mri_peak];
    }
    printf(
        "after smoothing, mri peak at %d, scaling input intensities "
        "by %2.3f\n",
        mri_peak,
        wm_means[max_T1_weighted_image] / mri_peak);
    if (mri_peak == 0) ErrorExit(ERROR_BADPARM, "GCAhistoScaleImageIntensities: could not find wm peak");
    scales[frame] = wm_means[max_T1_weighted_image] / mri_peak;

    MRIfree(&mri_frame);
    HISTOfree(&h_mri);
    HISTOfree(&h_smooth);
  }

  // scale each frame independently -xhan
  for (r = 0; r < mri->nframes; r++) {
    printf("scaling channel %d by %g\n", r, scales[r]);
    MRIscalarMulFrame(mri, mri, scales[r], r);
  }
  MRIfree(&mri_mask);
  return (NO_ERROR);
}

int GCAhisto(GCA *gca, int nbins, int **pcounts)
{
  int *counts, x, y, z;
  GCA_NODE *gcan;

  *pcounts = counts = (int *)calloc(nbins + 1, sizeof(int));

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        if (gcan->nlabels == 0 || (gcan->nlabels == 1 && GCAfindGC(gca, x, y, z, Unknown) != NULL)) {
          continue;
        }
        counts[gcan->nlabels]++;
      }
    }
  }

  return (NO_ERROR);
}

int GCArenormalizeToExample(GCA *gca, MRI *mri_seg, MRI *mri_T1)
{
  float intensities[MAX_CMA_LABEL + 1];
  int x, y, z, label, labels[MAX_CMA_LABEL + 1], width, height, depth, counts[MAX_CMA_LABEL + 1];

  for (label = 0; label <= MAX_CMA_LABEL; label++) {
    labels[label] = label;
  }
  memset(intensities, 0, MAX_CMA_LABEL * sizeof(float));
  memset(counts, 0, MAX_CMA_LABEL * sizeof(int));

  width = mri_seg->width;
  height = mri_seg->height;
  depth = mri_seg->height;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        label = nint(MRIgetVoxVal(mri_seg, x, y, z, 0));
        if (label == Gdiag_no) {
          DiagBreak();
        }
        if (label > MAX_CMA_LABEL) {
          ErrorPrintf(ERROR_BADPARM, "GCArenormalizeToExample: bad label %d", label);
          continue;
        }
        intensities[label] += MRIgetVoxVal(mri_T1, x, y, z, 0);
        counts[label]++;
      }
    }
  }

  intensities[0] = 0;  // don't estimate unknowns
  for (label = 1; label <= MAX_CMA_LABEL; label++) {
    HISTOGRAM *h;

    if (counts[label] <= 10) {
      continue;
    }
    if (counts[label] > 20) {
      int peak;
      h = MRIhistogramLabel(mri_T1, mri_seg, label, 256);
      peak = HISTOfindHighestPeakInRegion(h, 0, h->nbins);
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        HISTOplot(h, "h.plt");
      }
      intensities[label] = h->bins[peak];
    }
    else {
      intensities[label] /= (float)counts[label];
    }
  }
  GCArenormalizeIntensities(gca, labels, intensities, MAX_CMA_LABEL);

  return (NO_ERROR);
}

static GC1D *findGCInWindow(GCA *gca, int x0, int y0, int z0, int label, int wsize)
{
  GC1D *gc, *gc_min;
  int x, y, z, n, whalf;
  double dist, min_dist;
  GCA_NODE *gcan;

  min_dist = gca->node_width + gca->node_height + gca->node_depth;
  gc_min = NULL;
  whalf = (wsize - 1) / 2;
  for (x = x0 - whalf; x <= x0 + whalf; x++) {
    if (x < 0 || x >= gca->node_width) {
      continue;
    }
    for (y = y0 - whalf; y <= y0 + whalf; y++) {
      if (y < 0 || y >= gca->node_height) {
        continue;
      }
      for (z = z0 - whalf; z <= z0 + whalf; z++) {
        if (z < 0 || z >= gca->node_depth) {
          continue;
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels && min_dist > 1; n++) {
          if (gcan->labels[n] != label) {
            continue;
          }
          gc = &gcan->gcs[n];
          dist = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          if (dist < min_dist) {
            gc_min = gc;
            min_dist = dist;
          }
        }
      }
    }
  }
  return (gc_min);
}

GC1D *GCAfindClosestValidGC(GCA *gca, int x0, int y0, int z0, int label, int check_var)
{
  GC1D *gc, *gc_min;
  int x, y, z, n, wsize;
  double dist, min_dist, det;
  GCA_NODE *gcan;
  static MATRIX *m_cov_inv = NULL;

  min_dist = gca->node_width + gca->node_height + gca->node_depth;
  wsize = 1;
  gc_min = NULL;
  do {
    for (x = x0 - wsize; x <= x0 + wsize; x++) {
      if (x < 0 || x >= gca->node_width) {
        continue;
      }
      for (y = y0 - wsize; y <= y0 + wsize; y++) {
        if (y < 0 || y >= gca->node_height) {
          continue;
        }
        for (z = z0 - wsize; z <= z0 + wsize; z++) {
          if (z < 0 || z >= gca->node_depth) {
            continue;
          }
          dist = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          if (dist < wsize - 2) {
            continue;  // already checked
          }

          gcan = &gca->nodes[x][y][z];
          for (n = 0; n < gcan->nlabels; n++) {
            if (gcan->labels[n] != label) {
              continue;
            }
            gc = &gcan->gcs[n];
            det = covariance_determinant(gc, gca->ninputs);
            m_cov_inv = load_inverse_covariance_matrix(gc, m_cov_inv, gca->ninputs);
            if (m_cov_inv == NULL) {
              det = -1;
            }
            else {
              det *= 1;  // MatrixFree(&m_cov_inv) ;
            }
            if ((check_var && det <= MIN_DET) || gc->ntraining == 0) {
              continue;
            }
            if (dist < min_dist) {
              gc_min = gc;
              min_dist = dist;
              if (FEQUAL(dist, 1.0)) {
                return (gc_min);
              }
            }
          }
        }
      }
    }
    wsize += 2;  // search next ring
    if (wsize > MAX(MAX(gca->node_width, gca->node_height), gca->node_depth)) {
      break;
    }
  } while (gc_min == NULL);

  if (gc_min) /* found one in immediate nbhd */
  {
    return (gc_min);
  }

  /* couldn't find one close - search everywhere */
  for (x = 0; x < gca->node_width && min_dist > 1; x++) {
    for (y = 0; y < gca->node_height && min_dist > 1; y++) {
      for (z = 0; z < gca->node_depth && min_dist > 1; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels && min_dist > 1; n++) {
          if (gcan->labels[n] != label) {
            continue;
          }
          gc = &gcan->gcs[n];
          det = covariance_determinant(gc, gca->ninputs);
          if ((check_var && det <= 0) || gc->ntraining == 0) {
            continue;
          }
          dist = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          if (dist < min_dist) {
            gc_min = gc;
            min_dist = dist;
          }
        }
      }
    }
  }
  if (gc_min == NULL) {
    DiagBreak();
  }
  return (gc_min);
}

static int gcaCheck(GCA *gca)
{
  int x, y, z, ret = NO_ERROR, n, r, c, v;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap ;
  GC1D *gc;

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          gc = &gcan->gcs[n];
          for (v = r = 0; r < gca->ninputs; r++) {
            if (gcan->labels[n] > 1 && gc->means[0] < 1 && gc->ntraining > 10) 
	      DiagBreak();
            for (c = r; c < gca->ninputs; c++, v++)
              if (!check_finite("gcaCheckN1", gc->means[r]) || !check_finite("gcaCheckN2", gc->covars[v])) {
                ret = ERROR_BADPARM;
                DiagBreak();
              }
          }
        }
      }
    }
  }
  for (x = 0; x < gca->prior_width; x++) {
    for (y = 0; y < gca->prior_height; y++) {
      for (z = 0; z < gca->prior_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcap = &gca->priors[x][y][z];
	check_finite("gcaCheck1", gcap->nlabels) ;
	check_finite("gcaCheck2", gcap->max_labels) ;
	check_finite("gcaCheck3", gcap->total_training) ;
        for (n = 0; n < gcap->nlabels; n++) 
	{
	  if (check_finite("gcaCheck4", gcap->priors[n]) == 0)
	  {
	    ret = ERROR_BADPARM;
	    DiagBreak();
          }
        }
      }
    }
  }
  return (ret);
}
int GCAcomputeVoxelLikelihoods(
    GCA *gca, MRI *mri_in, int x, int y, int z, TRANSFORM *transform, int *labels, double *likelihoods)
{
  GCA_NODE *gcan;
  int xn, yn, zn, n;
  float vals[MAX_GCA_INPUTS];

  if (!GCAsourceVoxelToNode(gca, mri_in, transform, x, y, z, &xn, &yn, &zn)) {
    load_vals(mri_in, x, y, z, vals, gca->ninputs);
    gcan = &gca->nodes[xn][yn][zn];
    for (n = 0; n < gcan->nlabels; n++) {
      labels[n] = gcan->labels[n];
      likelihoods[n] = GCAcomputeConditionalDensity(&gcan->gcs[n], vals, gca->ninputs, labels[n]);
    }
    return (gcan->nlabels);
  }
  else {
    return 0;
  }
}

int GCAmaxLikelihoodLabel(const GCA_NODE *gcan, float *vals, int ninputs, float *plikelihood)
{
  double p, max_p;
  int n, best_label;

  if (gcan->nlabels == 0) {
    return (ERROR_BADPARM);
  }

  best_label = gcan->labels[0];
  max_p = GCAcomputeConditionalDensity(&gcan->gcs[0], vals, ninputs, gcan->labels[0]);

  for (n = 1; n < gcan->nlabels; n++) {
    p = GCAcomputeConditionalDensity(&gcan->gcs[n], vals, ninputs, gcan->labels[n]);
    if (p >= max_p) {
      max_p = p;
      best_label = gcan->labels[n];
    }
  }
  if (plikelihood) {
    *plikelihood = max_p;
  }
  return (best_label);
}
double GCAcomputeConditionalDensity(GC1D *gc, float *vals, int ninputs, int label)
{
  double p, dist;

/* compute 1-d Mahalanobis distance */
  {
    dist = GCAmahDist(gc, vals, ninputs);
    p = (1.0 / (pow(2 * M_PI, ninputs / 2.0) * sqrt(covariance_determinant(gc, ninputs)))) * exp(-0.5 * dist);
  }
  return (p);
}
double GCAcomputeNormalizedConditionalDensity(GCA *gca, int xp, int yp, int zp, float *vals, int label)
{
  GCA_PRIOR *gcap;
  GC1D *gc;
  int n;
  double p, total, plabel;

  gcap = &gca->priors[xp][yp][zp];
  if (gcap == NULL) {
    return (1.0);
  }
  if (gcap->nlabels <= 1) {
    return (1.0);
  }

  for (plabel = total = 0.0, n = 0; n < gcap->nlabels; n++) {
    gc = GCAfindPriorGC(gca, xp, yp, zp, gcap->labels[n]);
    if (!gc) {
      continue;
    }
    p = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, gcap->labels[n]);
    if (gcap->labels[n] == label) {
      plabel = p;
    }
    total += p;
  }

  if (!FZERO(total)) {
    plabel /= total;
  }
  return (plabel);
}

static double gcaComputeLogDensity(GC1D *gc, float *vals, int ninputs, float prior, int label)
{
  double log_p;

  log_p = GCAcomputeConditionalLogDensity(gc, vals, ninputs, label);
  log_p += log(prior);
  return (log_p);
}

double GCAcomputeConditionalLogDensity(const GC1D *gc, float *vals, int ninputs, int label)
{
  double log_p, det;

/* compute 1-d Mahalanobis distance */
  {
    det = covariance_determinant(gc, ninputs);
    log_p = -log(sqrt(det)) - .5 * GCAmahDist(gc, vals, ninputs);
  }
  return (log_p);
}

GCA_SAMPLE *GCAfindAllSamples(GCA *gca, int *pnsamples, int *exclude_list, int unknown_nbr_spacing)
{
  GCA_SAMPLE *gcas;
  GCA_PRIOR *gcap;
  GC1D *gc;
  int x, y, z, width, height, depth, n, i, max_label, nsamples, r, c, v;
  // int label;
  // int max_n;
  float max_p, prior;
  int badcount = 0;

  // going through gca->prior volume
  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;

  /////////////////////////////////////////////////////////////////////////
  // just because of C, you have to repeat the loop twice :-(.....
  // first to find the size needed, second time it gets the info
  ////////////////////////////////////////////////////////////////////////
  /* find the # of places in the atlas that the highest prior class is
     other than unknown.
  */
  //////////// first loop of finding counts //////////////////////////////
  for (nsamples = x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        gcap = &gca->priors[x][y][z];
        if (gcap == NULL) {
          continue;
        }
        max_p = 0;
        // max_n = -1;
        max_label = 0;
        // if no label, ignore
        if (gcap->nlabels == 0) {
          continue;
        }

        //////////////// debug code /////////////////////////////////
        if (x * gca->prior_spacing == Gx && y * gca->prior_spacing == Gy && z * gca->prior_spacing == Gz) {
          DiagBreak();
        }
        /////////////////////////////////////////////////////////////

        // go through the labels to find the max value
        // for prior probablity
        // and its label
        for (n = 0; n < gcap->nlabels; n++) {
          // label = gcap->labels[n];
          prior = gcap->priors[n];
          if (prior >= max_p) {
            // max_n = n;
            max_p = prior;
            max_label = gcap->labels[n];
          }
        }
        // if this label is in the exclude list, continue
        if (exclude_list && exclude_list[max_label] > 0) {
          continue;
        }
        // if the label is unknown and no +/-1
        // neighbor has different label, continue
        if (IS_UNKNOWN(max_label) && (different_nbr_max_labels(gca, x, y, z, unknown_nbr_spacing, 0) == 0)) {
          continue;
        }

        nsamples++;
      }
    }
  }
  ///////////////////////////////////////////////////////////////////////

  // allocate nsamples worth GCA_SAMPLE
  // samples are non-unknowns and unknown with neighbor not unknown
  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE));

  badcount = 0;
  // get the values for gcas
  // repeat the same thing again so that we can add info
  ////////////////// second loop of adding info ////////////////////////
  for (i = x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        gcap = &gca->priors[x][y][z];
        if (gcap == NULL) {
          continue;
        }
        max_p = 0;
        // max_n = -1;
        max_label = 0;
        // no label, ignore
        if (gcap->nlabels == 0) {
          continue;
        }

        for (n = 0; n < gcap->nlabels; n++) {
          // label = gcap->labels[n];
          prior = gcap->priors[n];
          if (prior >= max_p) {
            // max_n = n;
            max_p = prior;
            max_label = gcap->labels[n];
          }
        }
        if (exclude_list && exclude_list[max_label] > 0) {
          continue;
        }
        if (IS_UNKNOWN(max_label) && (different_nbr_max_labels(gca, x, y, z, unknown_nbr_spacing, 0) == 0)) {
          continue;
        }

        // store prior coordinates
        gcas[i].xp = x;
        gcas[i].yp = y;
        gcas[i].zp = z;
        //////////////////////////////////////
        // store talarached coordinate
        // talarached and prior has the same direction cosines
        // thus prior->talarached is only the scale difference
        // Note that we pick only one point in prior_spacing, though.
        gcas[i].x = x * gca->prior_spacing / gca->xsize;
        gcas[i].y = y * gca->prior_spacing / gca->ysize;
        gcas[i].z = z * gca->prior_spacing / gca->zsize;
        //////////////////////////////////////
        gcas[i].label = max_label;
        gcas_setPrior(gcas[i], max_p);
        gcas[i].means = (float *)calloc(gca->ninputs, sizeof(float));
        gcas[i].covars = (float *)calloc((gca->ninputs * (gca->ninputs + 1)) / 2, sizeof(float));
        if (!gcas[i].means || !gcas[i].covars)
          ErrorExit(ERROR_NOMEMORY,
                    "GCAfindAllSamples: could not allocate mean "
                    "(%d) and covariance (%d) matrices",
                    gca->ninputs,
                    gca->ninputs * (gca->ninputs + 1) / 2);
        gc = GCAfindPriorGC(gca, x, y, z, max_label);
        if (gc == NULL) {
          int xn, yn, zn;
          GCApriorToNode(gca, x, y, z, &xn, &yn, &zn);
          gc = GCAfindClosestValidGC(gca, xn, yn, zn, max_label, 0);
        }
        if (gc)  // found
        {
          for (v = r = 0; r < gca->ninputs; r++) {
            gcas[i].means[r] = gc->means[r];
            for (c = r; c < gca->ninputs; c++, v++) {
              gcas[i].covars[v] = gc->covars[v];
            }
          }
        }
        else  // not found
        {
          badcount++;
          for (v = r = 0; r < gca->ninputs; r++) {
            gcas[i].means[r] = 0.0;
            for (c = r; c < gca->ninputs; c++, v++) {
              if (c == r) {
                gcas[i].covars[v] = 1.0;
              }
              else {
                gcas[i].covars[v] = 0.0;
              }
            }
          }
        }
        gcas[i].log_p = 0;  // initialize
        i++;
      }
    }
  }
  if (badcount > 0) {
    fprintf(stdout, "**************************************\n");
    fprintf(stdout, "those with gc cannot be found = %d\n", badcount);
    fprintf(stdout, "**************************************\n");
  }
  *pnsamples = nsamples;
  return (gcas);
}
int GCAcomputeMAPlabelAtLocation(GCA *gca, int xp, int yp, int zp, float *vals, int *pmax_n, float *plog_p)
{
  GCA_PRIOR *gcap;
  GC1D *gc;
  int n, max_n, max_label;
  float log_p, max_log_p;

  gcap = &gca->priors[xp][yp][zp];
  if (gcap == NULL) {
    return (Unknown);
  }
  if (gcap->nlabels == 0) {
    if (plog_p) {
      *plog_p = 0.0;
    }
    if (pmax_n) {
      *pmax_n = -1;
    }
    return (Unknown);
  }

  // start from label[0]
  max_label = gcap->labels[0];
  max_n = 0;
  gc = GCAfindPriorGC(gca, xp, yp, zp, gcap->labels[0]);
  if (gc)
    max_log_p = gcaComputeLogDensity(gc, vals, gca->ninputs, gcap->priors[0], max_label);
  else {
    max_log_p = -100000;
  }
  // go through all the labels here
  for (n = 1; n < gcap->nlabels; n++) {
    gc = GCAfindPriorGC(gca, xp, yp, zp, gcap->labels[n]);
    if (!gc) {
      continue;
    }
    // calculate log_p
    log_p = gcaComputeLogDensity(gc, vals, gca->ninputs, gcap->priors[n], gcap->labels[n]);
    if (log_p > max_log_p) {
      max_log_p = log_p;
      max_n = n;
      max_label = gcap->labels[n];
    }
  }

  if (plog_p) {
    *plog_p = max_log_p;
  }
  if (pmax_n) {
    *pmax_n = max_n;
  }
  return (max_label);  // return most probable label
}
int GCAcomputeMLElabelAtLocation(GCA *gca, int xp, int yp, int zp, float *vals, int *pmax_n, float *plog_p)
{
  GCA_PRIOR *gcap;
  GC1D *gc;
  int n, max_n, max_label;
  float log_p, max_log_p;

  gcap = &gca->priors[xp][yp][zp];
  if (gcap == NULL || gcap->nlabels <= 0) {
    return (Unknown);
  }
  if (gcap->nlabels == 0) {
    if (plog_p) {
      *plog_p = 0.0;
    }
    if (pmax_n) {
      *pmax_n = -1;
    }
    return (Unknown);
  }

  max_label = gcap->labels[0];
  max_n = 0;
  gc = GCAfindPriorGC(gca, xp, yp, zp, gcap->labels[0]);
  if (gc)
    max_log_p = GCAcomputeConditionalLogDensity(gc, vals, gca->ninputs, max_label);
  else {
    max_log_p = -100000;
  }
  for (n = 1; n < gcap->nlabels; n++) {
    gc = GCAfindPriorGC(gca, xp, yp, zp, gcap->labels[n]);
    if (!gc) {
      continue;
    }
    log_p = GCAcomputeConditionalLogDensity(gc, vals, gca->ninputs, gcap->labels[n]);
    if (log_p > max_log_p) {
      max_log_p = log_p;
      max_n = n;
      max_label = gcap->labels[n];
    }
  }

  if (plog_p) {
    *plog_p = max_log_p;
  }
  if (pmax_n) {
    *pmax_n = max_n;
  }
  return (max_label);
}
GC1D *alloc_gcs(int nlabels, int flags, int ninputs)
{
  GC1D *gcs;
  int i;

  //  printf( "%s: nlabels=%i ninputs=%i\n",
  //    __FUNCTION__, nlabels, ninputs );

  gcs = (GC1D *)calloc(nlabels, sizeof(GC1D));
  if (gcs == NULL) ErrorExit(ERROR_NOMEMORY, "alloc_gcs(%d, %x): could not allocated %d gcs", nlabels, flags, nlabels);
  for (i = 0; i < nlabels; i++) {
    gcs[i].means = (float *)calloc(ninputs, sizeof(float));
    gcs[i].covars = (float *)calloc((ninputs * (ninputs + 1)) / 2, sizeof(float));
    if (!gcs[i].means || !gcs[i].covars)
      ErrorExit(ERROR_NOMEMORY,
                "could not allocate mean (%d) "
                "and covariance (%d) matrices",
                ninputs,
                ninputs * (ninputs + 1) / 2);
  }
  if (flags & GCA_NO_MRF) {
    return (gcs);
  }

  for (i = 0; i < nlabels; i++) {
    gcs[i].nlabels = (short *)calloc(GIBBS_NEIGHBORHOOD, sizeof(short));
    gcs[i].labels = (unsigned short **)calloc(GIBBS_NEIGHBORHOOD, sizeof(unsigned short *));
    gcs[i].label_priors = (float **)calloc(GIBBS_NEIGHBORHOOD, sizeof(float *));
    if (!gcs[i].nlabels || !gcs[i].labels || !gcs[i].label_priors)
      ErrorExit(ERROR_NOMEMORY, "alloc_gcs(%d, %x): could not allocated %d gcs(%d)", nlabels, flags, nlabels, i);
  }
  return (gcs);
}

int free_gcs(GC1D *gcs, int nlabels, int ninputs)
{
  int i, j;

  for (i = 0; i < nlabels; i++) {
    if (gcs[i].means) {
      free(gcs[i].means);
    }
    if (gcs[i].covars) {
      free(gcs[i].covars);
    }
    if (gcs[i].nlabels) /* gibbs stuff allocated */
    {
      for (j = 0; j < GIBBS_NEIGHBORHOOD; j++) {
        if (gcs[i].labels[j]) {
          free(gcs[i].labels[j]);
        }
        if (gcs[i].label_priors[j]) {
          free(gcs[i].label_priors[j]);
        }
      }
      free(gcs[i].nlabels);
      free(gcs[i].labels);
      free(gcs[i].label_priors);
    }
  }

  free(gcs);
  return (NO_ERROR);
}

int copy_gcs(int nlabels, GC1D *gcs_src, GC1D *gcs_dst, int ninputs)
{
  int i, j, k, r, c, v;

  for (i = 0; i < nlabels; i++) {
    for (v = r = 0; r < ninputs; r++) {
      gcs_dst[i].means[r] = gcs_src[i].means[r];
      for (c = r; c < ninputs; c++, v++) {
        gcs_dst[i].covars[v] = gcs_src[i].covars[v];
      }
    }
    gcs_dst[i].ntraining = gcs_src[i].ntraining;
    if (gcs_dst[i].nlabels == NULL) /* NO_MRF flag must be set */
    {
      continue;
    }
    for (j = 0; j < GIBBS_NEIGHBORHOOD; j++) {
      gcs_dst[i].nlabels[j] = gcs_src[i].nlabels[j];
      if (!gcs_dst[i].label_priors[j]) {
        gcs_dst[i].label_priors[j] = (float *)calloc(gcs_src[i].nlabels[j], sizeof(float));
        gcs_dst[i].labels[j] = (unsigned short *)calloc(gcs_src[i].nlabels[j], sizeof(unsigned short));
        if (!gcs_dst[i].label_priors[j] || !gcs_dst[i].labels[j])
          ErrorExit(ERROR_NOMEMORY, "copy_gcs(%d): i=%d, j=%d, could not gibbs\n", nlabels, i, j);
      }
      for (k = 0; k < gcs_src[i].nlabels[j]; k++) {
        gcs_dst[i].label_priors[j][k] = gcs_src[i].label_priors[j][k];
        gcs_dst[i].labels[j][k] = gcs_src[i].labels[j][k];
      }
    }
  }
  return (NO_ERROR);
}

int GCAfreeGibbs(GCA *gca)
{
  int x, y, z, n, i;
  GCA_NODE *gcan;
  GC1D *gc;

  if (gca->flags & GCA_NO_MRF) {
    return (NO_ERROR); /* already done */
  }

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          gc = &gcan->gcs[n];
          for (i = 0; i < GIBBS_NEIGHBORS; i++) {
            free(gc->label_priors[i]);
            free(gc->labels[i]);
            gc->label_priors[i] = NULL;
            gc->labels[i] = NULL;
          }
          free(gc->nlabels);
          free(gc->labels);
          free(gc->label_priors);
          gc->nlabels = NULL;
          gc->labels = NULL;
          gc->label_priors = NULL;
        }
      }
    }
  }
  gca->flags |= GCA_NO_MRF;
  return (NO_ERROR);
}

int GCAcomputeSampleCoords(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, int nsamples, TRANSFORM *transform)
{
  int n, x, y, z;

  TransformInvert(transform, mri);
  if (transform->type == LINEAR_VOX_TO_VOX) {
    LTA *lta;
    lta = (LTA *)transform->xform;
    fprintf(stdout, "INFO: compute sample coordinates transform\n");
    MatrixPrint(stdout, lta->xforms[0].m_L);
  }
  fprintf(stdout, "INFO: transform used\n");

  for (n = 0; n < nsamples; n++) {
    if (gcas[n].label == Gdiag_no) {
      DiagBreak();
    }
    if (!GCApriorToSourceVoxel(gca, mri, transform, gcas[n].xp, gcas[n].yp, gcas[n].zp, &x, &y, &z)) {
      gcas[n].x = x;
      gcas[n].y = y;
      gcas[n].z = z;
      if (DIAG_VERBOSE_ON && Gdiag & DIAG_SHOW)
        printf("label %d: (%d, %d, %d) <-- (%d, %d, %d)\n", gcas[n].label, gcas[n].xp, gcas[n].yp, gcas[n].zp, x, y, z);
    }
  }
  return (NO_ERROR);
}

static HISTOGRAM *gcaHistogramSamples(
    GCA *gca, GCA_SAMPLE *gcas, MRI *mri, TRANSFORM *transform, int nsamples, HISTOGRAM *histo, int frame)
{
  int i;
  double mean, var;
  double val;
  float fmin, fmax;

  if (!histo) {
    MRIvalRangeFrame(mri, &fmin, &fmax, frame);
    histo = HISTOalloc(nint(fmax - fmin + 1));
  }
  else {
    HISTOclear(histo, histo);
  }

  for (mean = var = 0.0, i = 0; i < nsamples; i++) {
    MRIsampleVolumeFrame(mri, gcas[i].x, gcas[i].y, gcas[i].z, frame, &val);
    histo->counts[(int)val]++;
    mean += val;
    var += (val * val);
  }

  mean /= (double)nsamples;
  var = var / (double)nsamples - mean * mean;
  return (histo);
}
int GCArenormalizeFromAtlas(GCA *gca, GCA *gca_template)
{
  int xs, ys, zs, xt, yt, zt, ns, label, v;
  float scale;
  GCA_NODE *gcan;
  GC1D *gc, *gct;

  scale = (float)gca_template->node_width / (float)gca->node_width;

  for (xs = 0; xs < gca->node_width; xs++) {
    xt = (int)(scale * xs);
    for (ys = 0; ys < gca->node_width; ys++) {
      yt = (int)(scale * ys);
      for (zs = 0; zs < gca->node_width; zs++) {
        zt = (int)(scale * zs);
        if (xs == Ggca_x && ys == Ggca_y && zs == Ggca_z) {
          DiagBreak();
        }
        if (xt == Ggca_x && yt == Ggca_y && zt == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[xs][ys][zs];
        for (ns = 0; ns < gcan->nlabels; ns++) {
          label = gcan->labels[ns];
          gc = &gcan->gcs[ns];
          gct = GCAfindGC(gca_template, xt, yt, zt, label);
          if (gct == NULL) {
            continue; /* label not in template GCA */
          }
          for (v = 0; v < gca->ninputs; v++) {
            gc->means[v] = gct->means[v];
          }
        }
      }
    }
  }
  return (NO_ERROR);
}

void load_vals(const MRI *mri_inputs, float x, float y, float z, float *vals, int ninputs)
{
  // go through all inputs and get values from inputs
  int n;
  for (n = 0; n < ninputs; n++) {
    // trilinear value at float x, y, z
    double val;
    MRIsampleVolumeFrame(mri_inputs, x, y, z, n, &val);
    vals[n] = val;
  }
}

#ifdef FASTER_MRI_EM_REGISTER
static void load_vals_xyzInt(const MRI *mri_inputs, int x, int y, int z, float *vals, int ninputs)
{
  MRIsampleVolumeFrame_xyzInt_nRange_floats(mri_inputs, x, y, z, 0, ninputs, vals);
}
#endif

double GCAmahDist(const GC1D *gc, const float *vals, const int ninputs)
{
  static VECTOR *v_means = NULL, *v_vals = NULL;
  static MATRIX *m_cov = NULL, *m_cov_inv;
  int i;
  double dsq;

  if (ninputs == 1) {
    float v;
    v = vals[0] - gc->means[0];
    dsq = v * v / gc->covars[0];
    return (dsq);
  }
  // printf("In GCAMahDist...ninputs = %d\n", ninputs);
  if (v_vals && ninputs != v_vals->rows) {
    VectorFree(&v_vals);
  }
  if (v_means && ninputs != v_means->rows) {
    VectorFree(&v_means);
  }
  if (m_cov && (ninputs != m_cov->rows || ninputs != m_cov->cols)) {
    MatrixFree(&m_cov);
    MatrixFree(&m_cov_inv);
  }
  v_means = load_mean_vector(gc, v_means, ninputs);
  //  MatrixPrint(stdout,v_means); //lz
  m_cov = load_covariance_matrix(gc, m_cov, ninputs);
  //  MatrixPrint(stdout,m_cov); //lz
  if (v_vals == NULL) {
    v_vals = VectorClone(v_means);
  }
  for (i = 0; i < ninputs; i++) {
    VECTOR_ELT(v_vals, i + 1) = vals[i];
  }

  VectorSubtract(v_means, v_vals, v_vals); /* v_vals now has mean removed */
  // MatrixPrint(stdout,v_vals); //lz
  m_cov_inv = MatrixInverse(m_cov, m_cov_inv);
  // MatrixPrint(stdout,m_cov_inv); //lz
  if (!m_cov_inv) {
    ErrorExit(ERROR_BADPARM, "singular covariance matrix!");
  }
  MatrixSVDInverse(m_cov, m_cov_inv);

  MatrixMultiply(m_cov_inv, v_vals, v_means);
  /* v_means is now inverse(cov) * v_vals */
  dsq = VectorDot(v_vals, v_means);

  return (dsq);
}
double GCAmahDistIdentityCovariance(GC1D *gc, float *vals, int ninputs)
{
  static VECTOR *v_means = NULL, *v_vals = NULL;
  int i;
  double dsq;

  if (v_vals && ninputs != v_vals->rows) {
    VectorFree(&v_vals);
  }
  if (v_means && ninputs != v_means->rows) {
    VectorFree(&v_means);
  }
  v_means = load_mean_vector(gc, v_means, ninputs);
  if (v_vals == NULL) {
    v_vals = VectorClone(v_means);
  }
  for (i = 0; i < ninputs; i++) {
    VECTOR_ELT(v_vals, i + 1) = vals[i];
  }

  VectorSubtract(v_means, v_vals, v_vals); /* v_vals now has mean removed */
  dsq = VectorDot(v_vals, v_vals);

  return (dsq);
}

double GCAcomputePosteriorDensity(GCA_PRIOR *gcap,
                                  GCA_NODE *gcan,
                                  int node_n,
                                  int prior_n,
                                  float *vals,
                                  int ninputs,
                                  int xn,
                                  int yn,
                                  int zn,
                                  GCA *gca)
{
  GC1D *gc;
  double p;

  if (node_n >= 0) {
    gc = &gcan->gcs[node_n];

    p = GCAcomputeConditionalDensity(gc, vals, ninputs, gcan->labels[node_n]);
    p *= getPrior(gcap, gcan->labels[node_n]);
  }
  else  // n refers to gcap index
  {
    gc = gcanGetGC(gcan, gcap->labels[prior_n]);

    if (gc == NULL) {
      if (xn >= 0) {
        gc = GCAfindClosestValidGC(gca, xn, yn, zn, gcap->priors[prior_n], 0);
      }
      if (gc == NULL) {
        return ((gcap->total_training == 0 ? 1e-10 : 1.0 / (10 * gcap->total_training)));
      }
    }
    p = GCAcomputeConditionalDensity(gc, vals, ninputs, gcap->labels[prior_n]);
    p *= getPrior(gcap, gcap->labels[prior_n]);
  }
  return (p);
}

VECTOR *load_mean_vector(const GC1D *gc, VECTOR *v_means, const int ninputs)
{
  int n;

  if (v_means == NULL) {
    v_means = VectorAlloc(ninputs, MATRIX_REAL);
  }

  for (n = 0; n < ninputs; n++) {
    VECTOR_ELT(v_means, n + 1) = gc->means[n];
  }

  return (v_means);
}

MATRIX *load_covariance_matrix(const GC1D *gc, MATRIX *m_cov, const int ninputs)
{
  int n, m, i;
  // printf("In load cov matrix; ninpupts = %d \n", ninputs);
  if (m_cov == NULL) {
    m_cov = MatrixAlloc(ninputs, ninputs, MATRIX_REAL);
  }
  // MatrixPrint(stdout,m_cov);

  for (i = n = 0; n < ninputs; n++)
    for (m = n; m < ninputs; m++, i++) {
      *MATRIX_RELT(m_cov, n + 1, m + 1) = gc->covars[i];
      if (m != n) {
        *MATRIX_RELT(m_cov, m + 1, n + 1) = gc->covars[i];
      }
    }

  // MatrixPrint(stdout,m_cov);
  return (m_cov);
}

int set_mean_vector(GC1D *gc, VECTOR *v_means, int ninputs)
{
  int n;

  for (n = 0; n < ninputs; n++) {
    gc->means[n] = VECTOR_ELT(v_means, n + 1);
  }

  return (NO_ERROR);
}

int set_covariance_matrix(GC1D *gc, MATRIX *m_cov, int ninputs)
{
  int n, m, i;

  for (i = n = 0; n < ninputs; n++)
    for (m = n; m < ninputs; m++, i++) {
      gc->covars[i] = *MATRIX_RELT(m_cov, n + 1, m + 1);
      if (m != n) {
        gc->covars[i] = *MATRIX_RELT(m_cov, m + 1, n + 1);
      }
    }

  return (NO_ERROR);
}

MATRIX *load_inverse_covariance_matrix(GC1D *gc, MATRIX *m_inv_cov, int ninputs)
{
  static MATRIX *m_cov[_MAX_FS_THREADS] = {NULL};
#ifdef HAVE_OPENMP
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif

  if (m_cov[tid] && (m_cov[tid]->rows != ninputs || m_cov[tid]->cols != ninputs)) {
    MatrixFree(&m_cov[tid]);
  }
  m_cov[tid] = load_covariance_matrix(gc, m_cov[tid], ninputs);
  m_inv_cov = MatrixInverse(m_cov[tid], m_inv_cov);
  return (m_inv_cov);
}

double covariance_determinant(const GC1D *gc, const int ninputs)
{
  double det;
  int tid;
  static MATRIX *m_cov[_MAX_FS_THREADS];

  if (ninputs == 1) {
    return (gc->covars[0]);
  }
#ifdef HAVE_OPENMP
  tid = omp_get_thread_num();
#else
  tid = 0;
#endif
  if (m_cov[tid] && (m_cov[tid]->rows != ninputs || m_cov[tid]->cols != ninputs)) {
    MatrixFree(&m_cov[tid]);
  }
  m_cov[tid] = load_covariance_matrix(gc, m_cov[tid], ninputs);
  det = MatrixDeterminant(m_cov[tid]);
  return (det);
}

static double gcaComputeSampleLogDensity(GCA_SAMPLE *gcas, float *vals, int ninputs)
{
  double log_p;

  log_p = gcaComputeSampleConditionalLogDensity(gcas, vals, ninputs, gcas->label);
  log_p += gcas_getPriorLog(*gcas);
  return (log_p);
}


#ifdef FASTER_MRI_EM_REGISTER
static double gcaComputeSampleLogDensity_1_input(GCA_SAMPLE *gcas, float val)
{
  double log_p;

  log_p = gcaComputeSampleConditionalLogDensity_1_input(gcas, val, gcas->label);
  log_p += gcas_getPriorLog(*gcas);
  return (log_p);
}
#endif

static double gcaComputeSampleConditionalLogDensity(GCA_SAMPLE *gcas, float *vals, int ninputs, int label)
{
  double log_p, det;

  {
    det = sample_covariance_determinant(gcas, ninputs);
    log_p = -log(sqrt(det)) - .5 * GCAsampleMahDist(gcas, vals, ninputs);
  }
  return (log_p);
}
#ifdef FASTER_MRI_EM_REGISTER
static double gcaComputeSampleConditionalLogDensity_1_input(GCA_SAMPLE *gcas, float val, int label)
{
  double log_p, det;
  det = gcas->covars[0];
  log_p = -log(sqrt(det)) - .5 * GCAsampleMahDist_1_input(gcas, val);
  return (log_p);
}
#endif

static VECTOR *load_sample_mean_vector(GCA_SAMPLE *gcas, VECTOR *v_means, int ninputs)
{
  int n;

  if (v_means == NULL) {
    v_means = VectorAlloc(ninputs, MATRIX_REAL);
  }

  for (n = 0; n < ninputs; n++) {
    VECTOR_ELT(v_means, n + 1) = gcas->means[n];
  }

  return (v_means);
}
static MATRIX *load_sample_covariance_matrix(GCA_SAMPLE *gcas, MATRIX *m_cov, int ninputs)
{
  int n, m, i;

  if (m_cov == NULL) {
    m_cov = MatrixAlloc(ninputs, ninputs, MATRIX_REAL);
  }

  for (i = n = 0; n < ninputs; n++)
    for (m = n; m < ninputs; m++, i++) {
      *MATRIX_RELT(m_cov, n + 1, m + 1) = gcas->covars[i];
    }

  return (m_cov);
}

static double sample_covariance_determinant(GCA_SAMPLE *gcas, int ninputs)
{
  double det;
  static MATRIX *m_cov = NULL;

  if (ninputs == 1) {
    return (gcas->covars[0]);
  }
  if (m_cov && (m_cov->rows != ninputs || m_cov->cols != ninputs)) {
    MatrixFree(&m_cov);
  }

  m_cov = load_sample_covariance_matrix(gcas, m_cov, ninputs);
  det = MatrixDeterminant(m_cov);
  return (det);
}

#ifdef FASTER_MRI_EM_REGISTER
static double GCAsampleMahDist_1_input(GCA_SAMPLE *gcas, float val) 
{
    float v = val - gcas->means[0];
    return v * v / gcas->covars[0];
}
#endif

double GCAsampleMahDist(GCA_SAMPLE *gcas, float *vals, int ninputs)
{
  static VECTOR *v_means = NULL, *v_vals = NULL;
  static MATRIX *m_cov = NULL, *m_cov_inv;
  int i;
  double dsq;

  if (ninputs == 1) {
    float v;
    v = vals[0] - gcas->means[0];
    dsq = v * v / gcas->covars[0];
    return (dsq);
  }

  if (v_vals && ninputs != v_vals->rows) {
    VectorFree(&v_vals);
  }
  if (v_means && ninputs != v_means->rows) {
    VectorFree(&v_means);
  }
  if (m_cov && (ninputs != m_cov->rows || ninputs != m_cov->cols)) {
    MatrixFree(&m_cov);
    MatrixFree(&m_cov_inv);
  }

  v_means = load_sample_mean_vector(gcas, v_means, ninputs);
  m_cov = load_sample_covariance_matrix(gcas, m_cov, ninputs);

  if (v_vals == NULL) {
    v_vals = VectorClone(v_means);
  }
  for (i = 0; i < ninputs; i++) {
    VECTOR_ELT(v_vals, i + 1) = vals[i];
  }

  VectorSubtract(v_means, v_vals, v_vals); /* v_vals now has mean removed */
  m_cov_inv = MatrixInverse(m_cov, m_cov_inv);
  if (!m_cov_inv) {
    ErrorExit(ERROR_BADPARM, "singular covariance matrix!");
  }

  MatrixMultiply(m_cov_inv, v_vals, v_means);
  /* v_means is now inverse(cov) * v_vals */
  dsq = VectorDot(v_vals, v_means);

  return (dsq);
}
static double gcaComputeSampleConditionalDensity(GCA_SAMPLE *gcas, float *vals, int ninputs, int label)
{
  double p, dist;

/* compute 1-d Mahalanobis distance */
  {
    dist = GCAsampleMahDist(gcas, vals, ninputs);
    p = 1 / sqrt(sample_covariance_determinant(gcas, ninputs)) * exp(-dist);
  }
  return (p);
}
int GCAlabelExists(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z, int label)
{
  return (GCAisPossible(gca, mri, label, transform, x, y, z, 0));
}

GC1D *GCAfindSourceGC(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z, int label)
{
  int xn, yn, zn;
  GC1D *gc = NULL;

  if (!GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn)) {
    gc = GCAfindGC(gca, xn, yn, zn, label);
  }
  return (gc);
}
int GCAfreeSamples(GCA_SAMPLE **pgcas, int nsamples)
{
  GCA_SAMPLE *gcas;
  int i;

  gcas = *pgcas;
  *pgcas = NULL;

  if (gcas == NULL) {
    return (NO_ERROR);
  }

  for (i = 0; i < nsamples; i++) {
    if (gcas[i].means) {
      free(gcas[i].means);
    }
    if (gcas[i].covars) {
      free(gcas[i].covars);
    }
  }
  free(gcas);
  return (NO_ERROR);
}
int GCAnormalizePD(GCA *gca, MRI *mri_inputs, TRANSFORM *transform)
{
  double mri_PD = 0.0, gca_PD = 0.0;
  int n, x, y, z, brain_node, nvals, label;
  GCA_PRIOR *gcap;
  GC1D *gc;
  double val;

  if (mri_inputs->nframes != 2 || gca->ninputs != 2)
    ErrorReturn(
        ERROR_BADPARM,
        (ERROR_BADPARM, "GCAnormalizePD: invalid input size (%d) or gca size (%d)", mri_inputs->nframes, gca->ninputs));

  for (nvals = x = 0; x < mri_inputs->width; x++) {
    for (y = 0; y < mri_inputs->height; y++) {
      for (z = 0; z < mri_inputs->depth; z++) {
        gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
        if (gcap == NULL || gcap->nlabels <= 0) {
          continue;
        }
        for (label = brain_node = 0, n = 0; !brain_node && n < gcap->nlabels; n++)
          if (IS_WM(gcap->labels[n]) && gcap->priors[n] > 0.5) {
            brain_node = 1;
            label = gcap->labels[n];
          }
        if (!brain_node) {
          continue; /* no valid data here */
        }
        gc = GCAfindSourceGC(gca, mri_inputs, transform, x, y, z, label);
        if (gc) {
          MRIsampleVolumeFrameType(mri_inputs, x, y, z, 1, SAMPLE_NEAREST, &val);
          nvals++;
          mri_PD += val;
          gca_PD += gc->means[1];
        }
      }
    }
  }
  mri_PD /= (double)nvals;
  gca_PD /= (double)nvals;
  printf("mean PD in volume %2.1f, mean GCA PD %2.1f - scaling by %2.3f\n", mri_PD, gca_PD, gca_PD / mri_PD);
  MRIscalarMulFrame(mri_inputs, mri_inputs, gca_PD / mri_PD, 1);

  return (NO_ERROR);
}
GCA *GCAcreateWeightedFlashGCAfromParameterGCA(
    GCA *gca_T1PD, double *TR, double *fa, double *TE, int nflash, double *wts, double lambda)
{
  GCA *gca_flash;
  GCA_PRIOR *gcap_src, *gcap_dst;
  GCA_NODE *gcan_src, *gcan_dst;
  GC1D *gc_src, *gc_dst;
  MATRIX *m_jacobian, *m_cov_src = NULL, *m_cov_dst = NULL, *m_jacobian_T = NULL, *m_tmp = NULL;
  VECTOR *v_wts, *v_wts_T;
  int n, x, y, z, i, j, v, label_count = 0;
  double T1, PD, label_means[MAX_GCA_INPUTS];

  if (gca_T1PD->ninputs != 2)
    ErrorExit(ERROR_BADPARM,
              "GCAcreateWeightedFlashGCAfromParameterGCA: "
              "input gca must be T1/PD (ninputs=%d, should be 2",
              gca_T1PD->ninputs);
  gca_flash = GCAalloc(nflash,
                       gca_T1PD->prior_spacing,
                       gca_T1PD->node_spacing,
                       gca_T1PD->node_width * gca_T1PD->node_spacing,
                       gca_T1PD->node_height * gca_T1PD->node_spacing,
                       gca_T1PD->node_depth * gca_T1PD->node_spacing,
                       GCA_NO_FLAGS);

  m_jacobian = MatrixAlloc(gca_flash->ninputs, gca_T1PD->ninputs, MATRIX_REAL);
  v_wts = VectorAlloc(nflash, MATRIX_REAL);
  for (i = 1; i <= nflash; i++) {
    VECTOR_ELT(v_wts, i) = wts[i - 1];
  }
  v_wts_T = VectorTranspose(v_wts, NULL);

  /* first copy over priors */
  for (x = 0; x < gca_flash->prior_width; x++) {
    for (y = 0; y < gca_flash->prior_height; y++) {
      for (z = 0; z < gca_flash->prior_depth; z++) {
        gcap_src = &gca_T1PD->priors[x][y][z];
        if (gcap_src == NULL) {
          continue;
        }
        gcap_dst = &gca_flash->priors[x][y][z];
        if (gcap_dst == NULL) {
          continue;
        }
        gcap_dst->nlabels = gcap_src->nlabels;
        if (gcap_dst == NULL) {
          continue;
        }
        if (gcap_src->nlabels > gcap_dst->max_labels) {
          free(gcap_dst->priors);
          free(gcap_dst->labels);

          gcap_dst->labels = (unsigned short *)calloc(gcap_src->nlabels, sizeof(unsigned short));
          if (!gcap_dst->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateWeightedFlashGCAfromParameterGCA:"
                      " couldn't allocate %d labels",
                      gcap_src->nlabels);

          gcap_dst->priors = (float *)calloc(gcap_src->nlabels, sizeof(float));
          if (!gcap_dst->priors)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateWeightedFlashGCAfromParameterGCA:"
                      " couldn't allocate %d priors",
                      gcap_src->nlabels);
          gcap_dst->max_labels = gcap_dst->nlabels;
        }
        gcap_dst->total_training = gcap_src->total_training;
        for (n = 0; n < gcap_src->nlabels; n++) {
          gcap_dst->labels[n] = gcap_src->labels[n];
          gcap_dst->priors[n] = gcap_src->priors[n];
        }
      }
    }
  }

  /* now copy over classifiers and Markov stuff, using Jacobian to */
  /* map to new image space */
  for (x = 0; x < gca_flash->node_width; x++) {
    for (y = 0; y < gca_flash->node_height; y++) {
      for (z = 0; z < gca_flash->node_depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        gcan_src = &gca_T1PD->nodes[x][y][z];
        gcan_dst = &gca_flash->nodes[x][y][z];
        gcan_dst->nlabels = gcan_src->nlabels;
        gcan_dst->total_training = gcan_src->total_training;
        if (gcan_src->nlabels > gcan_dst->max_labels) {
          free(gcan_dst->labels);
          free_gcs(gcan_dst->gcs, gcan_dst->max_labels, gca_flash->ninputs);

          gcan_dst->labels = (unsigned short *)calloc(gcan_src->nlabels, sizeof(unsigned short));
          if (!gcan_dst->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateWeightedFlashGCAfromParameterGCA:"
                      " couldn't allocate %d labels",
                      gcan_src->nlabels);

          gcan_dst->gcs = alloc_gcs(gcan_src->nlabels, GCA_NO_FLAGS, nflash);
          gcan_dst->max_labels = gcan_dst->nlabels;
        }
        for (n = 0; n < gcan_src->nlabels; n++) {
          gcan_dst->labels[n] = gcan_src->labels[n];
          gc_src = &gcan_src->gcs[n];
          gc_dst = &gcan_dst->gcs[n];
          gc_dst->ntraining = gc_src->ntraining;
          gc_dst->n_just_priors = gc_src->n_just_priors;
          gc_dst->regularized = gc_src->regularized;
          for (i = 0; i < GIBBS_NEIGHBORS; i++) {
            gc_dst->nlabels[i] = gc_src->nlabels[i];
            gc_dst->label_priors[i] = (float *)calloc(gc_src->nlabels[i], sizeof(float));
            if (!gc_dst->label_priors[i])
              ErrorExit(ERROR_NOMEMORY,
                        "GCAcreateWeightedFlashGCAfromParameterGCA: "
                        "to %d",
                        gc_src->nlabels[i]);
            gc_dst->labels[i] = (unsigned short *)calloc(gc_src->nlabels[i], sizeof(unsigned short));
            if (!gc_dst->labels)
              ErrorExit(ERROR_NOMEMORY,
                        "GCAcreateWeightedFlashGCAfromParameterGCA:"
                        " to %d",
                        gc_src->nlabels[i]);
            for (j = 0; j < gc_src->nlabels[i]; j++) {
              gc_dst->label_priors[i][j] = gc_src->label_priors[i][j];
              gc_dst->labels[i][j] = gc_src->labels[i][j];
            }
          }

          /* now map intensity and covariance info over */
          T1 = gc_src->means[0];
          PD = gc_src->means[1];
          if (Ggca_label == gcan_dst->labels[n]) {
            label_count++;
          }
          for (i = 0; i < gca_flash->ninputs; i++) {
            gc_dst->means[i] = FLASHforwardModel(T1, PD, TR[i], fa[i], TE[i]);
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label < 0 || Ggca_label == gcan_src->labels[n]))
              printf(
                  "gcan(%d, %d, %d) %s: image[%d] "
                  "(fa=%2.1f) predicted mean "
                  "(%2.1f,%2.1f) --> %2.1f\n",
                  x,
                  y,
                  z,
                  cma_label_to_name(gcan_dst->labels[n]),
                  i,
                  DEGREES(fa[i]),
                  T1,
                  PD,
                  gc_dst->means[i]);
            *MATRIX_RELT(m_jacobian, i + 1, 1) = dFlash_dT1(T1, PD, TR[i], fa[i], TE[i]);
            *MATRIX_RELT(m_jacobian, i + 1, 2) = dFlash_dPD(T1, PD, TR[i], fa[i], TE[i]);
            if (gcan_dst->labels[n] == Ggca_label) {
              label_means[i] += gc_dst->means[i];
            }
          }
#define MIN_T1 50
          if (T1 < MIN_T1) {
            gc_dst->regularized = 1;
            m_cov_dst = MatrixIdentity(gca_flash->ninputs, m_cov_dst);
          }
          else {
            m_cov_src = load_covariance_matrix(gc_src, m_cov_src, gca_T1PD->ninputs);
            m_jacobian_T = MatrixTranspose(m_jacobian, m_jacobian_T);
            m_tmp = MatrixMultiply(m_cov_src, m_jacobian_T, m_tmp);
            m_cov_dst = MatrixMultiply(m_jacobian, m_tmp, m_cov_dst);
          }
          for (v = i = 0; i < gca_flash->ninputs; i++) {
            for (j = i; j < gca_flash->ninputs; j++, v++) {
              gc_dst->covars[v] = *MATRIX_RELT(m_cov_dst, i + 1, j + 1);
            }
          }
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label < 0 || Ggca_label == gcan_src->labels[n])) {
            printf("predicted covariance matrix:\n");
            MatrixPrint(stdout, m_cov_dst);
          }
        }
      }
    }
  }

  if (Ggca_label >= 0) {
    printf("label %s (%d): means = ", cma_label_to_name(Ggca_label), Ggca_label);
    for (i = 0; i < gca_flash->ninputs; i++) {
      label_means[i] /= (double)label_count;
      printf("%2.1f ", label_means[i]);
    }
    printf("\n");
  }

  /* check and fix singular covariance matrixces */
  GCAregularizeCovarianceMatrices(gca_flash, lambda);
  GCAfixSingularCovarianceMatrices(gca_flash);
  MatrixFree(&m_jacobian);
  MatrixFree(&m_cov_src);
  MatrixFree(&m_cov_dst);
  MatrixFree(&m_jacobian_T);
  MatrixFree(&m_tmp);
  VectorFree(&v_wts);
  VectorFree(&v_wts_T);

  gca_flash->type = GCA_FLASH;
  GCAcopyDCToGCA(gca_T1PD, gca_flash);

  return (gca_flash);
}
GCA *GCAcreateFlashGCAfromParameterGCA(GCA *gca_T1PD, double *TR, double *fa, double *TE, int nflash, double lambda)
{
  GCA *gca_flash;
  GCA_PRIOR *gcap_src, *gcap_dst;
  GCA_NODE *gcan_src, *gcan_dst;
  GC1D *gc_src, *gc_dst;
  MATRIX *m_jacobian, *m_cov_src = NULL, *m_cov_dst = NULL, *m_jacobian_T = NULL, *m_tmp = NULL;
  int n, x, y, z, i, j, v, label_count = 0;
  double T1, PD, label_means[MAX_GCA_INPUTS];

  if (gca_T1PD->ninputs != 2)
    ErrorExit(ERROR_BADPARM,
              "GCAcreateFlashGCAfromParameterGCA: "
              "input gca must be T1/PD (ninputs=%d, should be 2",
              gca_T1PD->ninputs);
  // gca_flash will have gca->ninputs = nflash
  gca_flash = GCAalloc(nflash,
                       gca_T1PD->prior_spacing,
                       gca_T1PD->node_spacing,
                       gca_T1PD->node_width * gca_T1PD->node_spacing,
                       gca_T1PD->node_height * gca_T1PD->node_spacing,
                       gca_T1PD->node_depth * gca_T1PD->node_spacing,
                       GCA_NO_FLAGS);

  m_jacobian = MatrixAlloc(gca_flash->ninputs, gca_T1PD->ninputs, MATRIX_REAL);

  /* first copy over priors */
  for (x = 0; x < gca_flash->prior_width; x++) {
    for (y = 0; y < gca_flash->prior_height; y++) {
      for (z = 0; z < gca_flash->prior_depth; z++) {
        gcap_src = &gca_T1PD->priors[x][y][z];
        if (gcap_src == NULL) {
          continue;
        }
        gcap_dst = &gca_flash->priors[x][y][z];
        if (gcap_dst == NULL) {
          continue;
        }
        gcap_dst->nlabels = gcap_src->nlabels;
        if (gcap_src->nlabels > gcap_dst->max_labels) {
          free(gcap_dst->priors);
          free(gcap_dst->labels);

          gcap_dst->labels = (unsigned short *)calloc(gcap_src->nlabels, sizeof(unsigned short));
          if (!gcap_dst->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateFlashGCAfromParameterGCA: "
                      "couldn't allocate %d labels",
                      gcap_src->nlabels);

          gcap_dst->priors = (float *)calloc(gcap_src->nlabels, sizeof(float));
          if (!gcap_dst->priors)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateFlashGCAfromParameterGCA: "
                      "couldn't allocate %d priors",
                      gcap_src->nlabels);
          gcap_dst->max_labels = gcap_dst->nlabels;
        }
        gcap_dst->total_training = gcap_src->total_training;
        for (n = 0; n < gcap_src->nlabels; n++) {
          gcap_dst->labels[n] = gcap_src->labels[n];
          gcap_dst->priors[n] = gcap_src->priors[n];
        }
      }
    }
  }

  /* now copy over classifiers and Markov stuff, */
  /* using Jacobian to map to new image space */
  for (x = 0; x < gca_flash->node_width; x++) {
    for (y = 0; y < gca_flash->node_height; y++) {
      for (z = 0; z < gca_flash->node_depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        gcan_src = &gca_T1PD->nodes[x][y][z];
        gcan_dst = &gca_flash->nodes[x][y][z];
        gcan_dst->nlabels = gcan_src->nlabels;
        gcan_dst->total_training = gcan_src->total_training;
        if (gcan_src->nlabels > gcan_dst->max_labels) {
          free(gcan_dst->labels);
          free_gcs(gcan_dst->gcs, gcan_dst->max_labels, gca_flash->ninputs);

          gcan_dst->labels = (unsigned short *)calloc(gcan_src->nlabels, sizeof(unsigned short));
          if (!gcan_dst->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateFlashGCAfromParameterGCA: "
                      "couldn't allocate %d labels",
                      gcan_src->nlabels);

          gcan_dst->gcs = alloc_gcs(gcan_src->nlabels, GCA_NO_FLAGS, nflash);
          gcan_dst->max_labels = gcan_dst->nlabels;
        }
        for (n = 0; n < gcan_src->nlabels; n++) {
          gcan_dst->labels[n] = gcan_src->labels[n];
          gc_src = &gcan_src->gcs[n];
          gc_dst = &gcan_dst->gcs[n];
          gc_dst->ntraining = gc_src->ntraining;
          gc_dst->n_just_priors = gc_src->n_just_priors;
          for (i = 0; i < GIBBS_NEIGHBORS; i++) {
            gc_dst->nlabels[i] = gc_src->nlabels[i];
            gc_dst->label_priors[i] = (float *)calloc(gc_src->nlabels[i], sizeof(float));
            if (!gc_dst->label_priors[i])
              ErrorExit(ERROR_NOMEMORY, "GCAcreateFlashGCAfromParameterGCA: to %d", gc_src->nlabels[i]);
            gc_dst->labels[i] = (unsigned short *)calloc(gc_src->nlabels[i], sizeof(unsigned short));
            if (!gc_dst->labels)
              ErrorExit(ERROR_NOMEMORY, "GCAcreateFlashGCAfromParameterGCA: to %d", gc_src->nlabels[i]);
            for (j = 0; j < gc_src->nlabels[i]; j++) {
              gc_dst->label_priors[i][j] = gc_src->label_priors[i][j];
              gc_dst->labels[i][j] = gc_src->labels[i][j];
            }
          }

          /* now map intensity and covariance info over */
          T1 = gc_src->means[0];
          PD = gc_src->means[1];
          if (T1 < 0) {
            printf("WARN: ******************************************\n");
            printf("WARN: (%d, %d, %d) has T1 = %f < 0 and PD = %f\n", x, y, z, T1, PD);
            printf("WARN: nlabels = %d\n", gcan_src->nlabels);
            for (i = 0; i < gcan_src->nlabels; ++i) printf("WARN: %d: label = %d\n", i, gcan_src->labels[i]);
            printf("WARN: make T1 = 10 msec\n");
            printf("WARN: ******************************************\n");
          }
          T1 = MAX(T1, 10);
          PD = MAX(PD, 0);
          if (Ggca_label == gcan_dst->labels[n]) {
            label_count++;
          }
          for (i = 0; i < gca_flash->ninputs; i++) {
            gc_dst->means[i] = FLASHforwardModel(T1, PD, TR[i], fa[i], TE[i]);
            if (x == Gx && y == Gy && z == Gz && (Ggca_label < 0 || Ggca_label == gcan_src->labels[n]))
              printf(
                  "gcan(%d, %d, %d) %s: image[%d] "
                  "(fa=%2.1f) predicted mean "
                  "(%2.1f,%2.1f) --> %2.1f\n",
                  x,
                  y,
                  z,
                  cma_label_to_name(gcan_dst->labels[n]),
                  i,
                  DEGREES(fa[i]),
                  T1,
                  PD,
                  gc_dst->means[i]);
            *MATRIX_RELT(m_jacobian, i + 1, 1) = dFlash_dT1(T1, PD, TR[i], fa[i], TE[i]);
            *MATRIX_RELT(m_jacobian, i + 1, 2) = dFlash_dPD(T1, PD, TR[i], fa[i], TE[i]);
            if (gcan_dst->labels[n] == Ggca_label) {
              label_means[i] += gc_dst->means[i];
            }
          }
#define MIN_T1 50
          if (T1 < MIN_T1) {
            m_cov_dst = MatrixIdentity(gca_flash->ninputs, m_cov_dst);
          }
          else {
            m_cov_src = load_covariance_matrix(gc_src, m_cov_src, gca_T1PD->ninputs);
            m_jacobian_T = MatrixTranspose(m_jacobian, m_jacobian_T);
            m_tmp = MatrixMultiply(m_cov_src, m_jacobian_T, m_tmp);
            m_cov_dst = MatrixMultiply(m_jacobian, m_tmp, m_cov_dst);
          }
          for (v = i = 0; i < gca_flash->ninputs; i++) {
            for (j = i; j < gca_flash->ninputs; j++, v++) {
              gc_dst->covars[v] = *MATRIX_RELT(m_cov_dst, i + 1, j + 1);
            }
          }
          if (x == Gx && y == Gy && z == Gz && (Ggca_label < 0 || Ggca_label == gcan_src->labels[n])) {
            printf("predicted covariance matrix:\n");
            MatrixPrint(stdout, m_cov_dst);
          }
        }
      }
    }
  }

  if (Ggca_label >= 0) {
    printf("label %s (%d): means = ", cma_label_to_name(Ggca_label), Ggca_label);
    for (i = 0; i < gca_flash->ninputs; i++) {
      label_means[i] /= (double)label_count;
      printf("%2.1f ", label_means[i]);
    }
    printf("\n");
  }

  /* check and fix singular covariance matrixces */
  GCAregularizeCovarianceMatrices(gca_flash, lambda);
  GCAfixSingularCovarianceMatrices(gca_flash);
  MatrixFree(&m_jacobian);
  MatrixFree(&m_cov_src);
  MatrixFree(&m_cov_dst);
  MatrixFree(&m_jacobian_T);
  MatrixFree(&m_tmp);

  gca_flash->type = GCA_FLASH;
  GCAcopyDCToGCA(gca_T1PD, gca_flash);

  return (gca_flash);
}

int GCAfixSingularCovarianceMatrices(GCA *gca)
{
  int x, y, z, fixed = 0, i, r, c, n, num, nparams, regularized = 0;
  GCA_NODE *gcan;
  GC1D *gc;
  double det, vars[MAX_GCA_INPUTS], min_det;
  MATRIX *m_cov_inv, *m_cov = NULL;

  nparams = (gca->ninputs * (gca->ninputs + 1)) / 2 + gca->ninputs;
  /* covariance matrix and means */

  memset(vars, 0, sizeof(vars));

  if (gca->total_training <= 1 && gca->prior_spacing <= 1 &&
      gca->node_spacing <= 1)  // degenerate case - can't estimate vars
  {
    for (num = 0, x = 0; x < gca->node_width; x++) {
      for (y = 0; y < gca->node_height; y++) {
        for (z = 0; z < gca->node_depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          gcan = &gca->nodes[x][y][z];
          for (n = 0; n < gcan->nlabels; n++) {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
              DiagBreak();
            }
            gc = &gcan->gcs[n];
            for (r = 0; r < gca->ninputs; r++) {
              vars[r] += SQR((gc->means[r] * 0.1));
            }
            num++;
          }
        }
      }
    }
  }
  else {
    for (num = 0, x = 0; x < gca->node_width; x++) {
      for (y = 0; y < gca->node_height; y++) {
        for (z = 0; z < gca->node_depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          gcan = &gca->nodes[x][y][z];
          for (n = 0; n < gcan->nlabels; n++) {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
              DiagBreak();
            }
            gc = &gcan->gcs[n];
            det = covariance_determinant(gc, gca->ninputs);
            if ((gc->ntraining == 0 && det > 1) || (gc->ntraining * gca->ninputs > 2.5 * nparams))
            /* enough to estimate parameters */
            {
              m_cov = load_covariance_matrix(gc, m_cov, gca->ninputs);
              for (r = 0; r < gca->ninputs; r++) {
                vars[r] += *MATRIX_RELT(m_cov, r + 1, r + 1);
              }
              num++;
            }
          }
        }
      }
    }
  }
  if (m_cov) {
    MatrixFree(&m_cov);
  }
  if (num >= 1) {
    printf("average std = ");
    for (min_det = 1.0, r = 0; r < gca->ninputs; r++) {
      vars[r] /= (float)num;
      printf("%2.1f ", sqrt(vars[r]));
      min_det *= vars[r];
    }
    min_det = min_det / pow(10.0, gca->ninputs);
    printf("  using min determinant for regularization = %2.1f\n", min_det);
  }
  else {
    min_det = MIN_DET;
  }

  for (regularized = fixed = x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
            DiagBreak();
          }
          gc = &gcan->gcs[n];
          det = covariance_determinant(gc, gca->ninputs);
          m_cov_inv = load_inverse_covariance_matrix(gc, NULL, gca->ninputs);

          if (det <= 0 || m_cov_inv == NULL) {
            fixed++;
            gc->regularized = 1;
            for (i = r = 0; r < gca->ninputs; r++) {
              for (c = r; c < gca->ninputs; c++, i++) {
                if (r == c) {
                  gc->covars[i] += vars[r];
                }
                /* mean of other variances at
                   this location */
              }
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
                MATRIX *m;
                printf(
                    "fixing singular covariance matrix for %s "
                    "@ (%d, %d, %d):\n",
                    cma_label_to_name(gcan->labels[n]),
                    x,
                    y,
                    z);
                m = load_covariance_matrix(gc, NULL, gca->ninputs);
                MatrixPrint(stdout, m);
                MatrixFree(&m);
              }
            }
          }
          else /* not singular - check if it is ill-conditioned */
          {
            if (gc->regularized == 0) {
              DiagBreak();
            }
            if ((gc->ntraining * gca->ninputs < 2 * nparams && det < 0.1) || (det < min_det)) {
              gc->regularized = 1;
              regularized++;
              for (i = r = 0; r < gca->ninputs; r++) {
                for (c = r; c < gca->ninputs; c++, i++) {
                  if (r == c) {
                    gc->covars[i] += vars[r];
                  }
                  /* mean of overall variance
                     in this image */
                }
              }
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
                MATRIX *m;
                printf(
                    "fixing ill-conditioned covariance "
                    "matrix for %s @ (%d, %d, %d):\n",
                    cma_label_to_name(gcan->labels[n]),
                    x,
                    y,
                    z);
                m = load_covariance_matrix(gc, NULL, gca->ninputs);
                MatrixPrint(stdout, m);
                MatrixFree(&m);
              }
            }
          }

          if (m_cov_inv) {
            MatrixFree(&m_cov_inv);
          }
          det = covariance_determinant(gc, gca->ninputs);
          m_cov_inv = load_inverse_covariance_matrix(gc, NULL, gca->ninputs);
          if (det <= min_det || m_cov_inv == NULL) {
            printf(
                "warning: regularization of node (%d, %d, %d) "
                "label %s failed\n",
                x,
                y,
                z,
                cma_label_to_name(gcan->labels[n]));
            DiagBreak();
          }
          if (m_cov_inv) {
            MatrixFree(&m_cov_inv);
          }
        }
      }
    }
  }

  printf(
      "%d singular and %d ill-conditioned covariance"
      " matrices regularized\n",
      fixed,
      regularized);
  return (NO_ERROR);
}
int GCAregularizeCovarianceMatrices(GCA *gca, double lambda)
{
  int x, y, z, r, n, num, nparams;
  GCA_NODE *gcan;
  GC1D *gc;
  double det, vars[MAX_GCA_INPUTS], min_det;
  MATRIX *m_cov = NULL;

  nparams = (gca->ninputs * (gca->ninputs + 1)) / 2 + gca->ninputs;
  /* covariance matrix and means */

  memset(vars, 0, sizeof(vars));

  for (num = 0, x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
            DiagBreak();
          }
          gc = &gcan->gcs[n];
          det = covariance_determinant(gc, gca->ninputs);
          if ((gc->ntraining == 0 && det > 1) || (gc->ntraining * gca->ninputs > 2.5 * nparams))
          /* enough to estimate parameters */
          {
            m_cov = load_covariance_matrix(gc, m_cov, gca->ninputs);
            for (r = 0; r < gca->ninputs; r++) {
              vars[r] += *MATRIX_RELT(m_cov, r + 1, r + 1);
            }
            num++;
          }
        }
      }
    }
  }

  if (num >= 1) {
    printf("average std = ");
    for (min_det = 1.0, r = 0; r < gca->ninputs; r++) {
      vars[r] /= (float)num;
      printf("%2.1f ", sqrt(vars[r]));
      min_det *= vars[r];
    }
    min_det = min_det / pow(10.0, gca->ninputs);
    printf("  using min determinant for regularization = %2.1f\n", min_det);
  }
  else {
    min_det = MIN_DET;
  }

  /* discard all previous stuff and just regularize by adding a */
  /* fixed constant independent of variance */
  for (r = 0; r < gca->ninputs; r++) {
    vars[r] = lambda;
  }

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          int r, c, v;
          MATRIX *m;

          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
            DiagBreak();
          }
          gc = &gcan->gcs[n];

          m = load_covariance_matrix(gc, NULL, gca->ninputs);
          for (r = v = 0; r < gca->ninputs; r++)
            for (c = r; c < gca->ninputs; c++, v++) {
              if (r == c) {
                gc->covars[v] += vars[r];
              }
            }
          MatrixFree(&m);
          m = load_covariance_matrix(gc, NULL, gca->ninputs);
          v = 0;
        }
      }
    }
  }

  return (NO_ERROR);
}
int GCAsetFlashParameters(GCA *gca, double *TRs, double *FAs, double *TEs)
{
  int n;

  printf("setting GCA flash parameters to:\n");
  for (n = 0; n < gca->ninputs; n++) {
    gca->TRs[n] = TRs[n];
    gca->FAs[n] = FAs[n];
    gca->TEs[n] = TEs[n];
    printf("TR=%2.1f msec, FA=%2.1f deg, TE=%2.1f msec\n", TRs[n], DEGREES(FAs[n]), TEs[n]);
  }
  gca->type = GCA_FLASH;  // mark gca type
  return (NO_ERROR);
}

GCA *GCAcreateFlashGCAfromFlashGCA(GCA *gca_flash_src, double *TR, double *fa, double *TE, int nflash)
{
  GCA *gca_flash_dst;
  GCA_PRIOR *gcap_src, *gcap_dst;
  GCA_NODE *gcan_src, *gcan_dst;
  GC1D *gc_src, *gc_dst;
  MATRIX *m_jac_dst, *m_jac_src, *m_cov_src, *m_cov_dst, *m_pinv_src, *m_cov_T1PD;
  int n, x, y, z, i, j, v, label_count = 0;
  double T1, PD, label_means[MAX_GCA_INPUTS];

  m_cov_T1PD = m_cov_dst = m_pinv_src = m_cov_src = NULL;
  FlashBuildLookupTables(gca_flash_src->ninputs, gca_flash_src->TRs, gca_flash_src->FAs, gca_flash_src->TEs);

  gca_flash_dst = GCAalloc(nflash,
                           gca_flash_src->prior_spacing,
                           gca_flash_src->node_spacing,
                           gca_flash_src->node_width * gca_flash_src->node_spacing,
                           gca_flash_src->node_height * gca_flash_src->node_spacing,
                           gca_flash_src->node_depth * gca_flash_src->node_spacing,
                           GCA_NO_FLAGS);

  for (n = 0; n < nflash; n++) {
    gca_flash_dst->TRs[n] = TR[n];
    gca_flash_dst->FAs[n] = fa[n];
    gca_flash_dst->TEs[n] = TE[n];
  }

  m_jac_src = MatrixAlloc(gca_flash_src->ninputs, 2, MATRIX_REAL);
  /* T1/PD --> src jacobian */
  m_jac_dst = MatrixAlloc(gca_flash_dst->ninputs, 2, MATRIX_REAL);
  /* T1/PD --> dst jacobian */

  /* first copy over priors */
  for (x = 0; x < gca_flash_dst->prior_width; x++) {
    for (y = 0; y < gca_flash_dst->prior_height; y++) {
      for (z = 0; z < gca_flash_dst->prior_depth; z++) {
        gcap_src = &gca_flash_src->priors[x][y][z];
        if (gcap_src == NULL) {
          continue;
        }
        gcap_dst = &gca_flash_dst->priors[x][y][z];
        if (gcap_dst == NULL) {
          continue;
        }
        gcap_dst->nlabels = gcap_src->nlabels;
        if (gcap_src->nlabels > gcap_dst->max_labels) {
          free(gcap_dst->priors);
          free(gcap_dst->labels);

          gcap_dst->labels = (unsigned short *)calloc(gcap_src->nlabels, sizeof(unsigned short));
          if (!gcap_dst->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateFlashGCAfromFlashGCA: "
                      "couldn't allocate %d labels",
                      gcap_src->nlabels);

          gcap_dst->priors = (float *)calloc(gcap_src->nlabels, sizeof(float));
          if (!gcap_dst->priors)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateFlashGCAfromFlashGCA: "
                      "couldn't allocate %d priors",
                      gcap_src->nlabels);
          gcap_dst->max_labels = gcap_dst->nlabels;
        }
        gcap_dst->total_training = gcap_src->total_training;
        for (n = 0; n < gcap_src->nlabels; n++) {
          gcap_dst->labels[n] = gcap_src->labels[n];
          gcap_dst->priors[n] = gcap_src->priors[n];
        }
      }
    }
  }

  /* now copy over classifiers and Markov stuff, */
  /* using Jacobian to map to new image space */
  for (x = 0; x < gca_flash_dst->node_width; x++) {
    for (y = 0; y < gca_flash_dst->node_height; y++) {
      for (z = 0; z < gca_flash_dst->node_depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        gcan_src = &gca_flash_src->nodes[x][y][z];
        gcan_dst = &gca_flash_dst->nodes[x][y][z];
        gcan_dst->nlabels = gcan_src->nlabels;
        gcan_dst->total_training = gcan_src->total_training;
        if (gcan_src->nlabels > gcan_dst->max_labels) {
          free(gcan_dst->labels);
          for (n = 0; n < gcan_dst->max_labels; n++) {
            gc_dst = &gcan_dst->gcs[n];
            for (i = 0; i < GIBBS_NEIGHBORS; i++) {
              if (gc_dst->label_priors[i]) {
                free(gc_dst->label_priors[i]);
              }
              if (gc_dst->labels[i]) {
                free(gc_dst->labels[i]);
              }
            }
            if (gc_dst->nlabels) {
              free(gc_dst->nlabels);
            }
            if (gc_dst->labels) {
              free(gc_dst->labels);
            }
            if (gc_dst->label_priors) {
              free(gc_dst->label_priors);
            }
          }

          gcan_dst->labels = (unsigned short *)calloc(gcan_src->nlabels, sizeof(unsigned short));
          if (!gcan_dst->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAcreateFlashGCAfromFlashGCA: "
                      "couldn't allocate %d labels",
                      gcan_src->nlabels);

          gcan_dst->gcs = alloc_gcs(gcan_src->nlabels, GCA_NO_FLAGS, nflash);
        }
        for (n = 0; n < gcan_src->nlabels; n++) {
          gcan_dst->labels[n] = gcan_src->labels[n];
          gc_src = &gcan_src->gcs[n];
          gc_dst = &gcan_dst->gcs[n];
          gc_dst->ntraining = gc_src->ntraining;
          gc_dst->n_just_priors = gc_src->n_just_priors;
          for (i = 0; i < GIBBS_NEIGHBORS; i++) {
            gc_dst->nlabels[i] = gc_src->nlabels[i];
            gc_dst->label_priors[i] = (float *)calloc(gc_src->nlabels[i], sizeof(float));
            if (!gc_dst->label_priors[i])
              ErrorExit(ERROR_NOMEMORY, "GCAcreateFlashGCAfromFlashGCA: to %d", gc_src->nlabels[i]);
            gc_dst->labels[i] = (unsigned short *)calloc(gc_src->nlabels[i], sizeof(unsigned short));
            if (!gc_dst->labels) ErrorExit(ERROR_NOMEMORY, "GCAcreateFlashGCAfromFlashGCA: to %d", gc_src->nlabels[i]);
            for (j = 0; j < gc_src->nlabels[i]; j++) {
              gc_dst->label_priors[i][j] = gc_src->label_priors[i][j];
              gc_dst->labels[i][j] = gc_src->labels[i][j];
            }
          }

          /* now map intensity and covariance info over */
          compute_T1_PD(gca_flash_src->ninputs,
                        gc_src->means,
                        gca_flash_src->TRs,
                        gca_flash_src->FAs,
                        gca_flash_src->TEs,
                        &T1,
                        &PD);
          T1 = MAX(T1, 10);
          PD = MAX(PD, 0);
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label < 0 || Ggca_label == gcan_src->labels[n])) {
            int j;
            printf("gcan(%d, %d, %d) %s: means=[", x, y, z, cma_label_to_name(gcan_src->labels[n]));
            for (j = 0; j < gca_flash_src->ninputs; j++) {
              printf(" %2.1f", gc_src->means[j]);
            }
            printf("] T1/PD=%2.1f/%2.1f\n", T1, PD);
          }
          if (Ggca_label == gcan_dst->labels[n]) {
            label_count++;
          }
          for (i = 0; i < gca_flash_dst->ninputs; i++) {
            gc_dst->means[i] = FLASHforwardModel(T1, PD, TR[i], fa[i], TE[i]);
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label < 0 || Ggca_label == gcan_src->labels[n]))
              printf(
                  "gcan(%d, %d, %d) %s: image[%d] "
                  "(fa=%2.1f) predicted mean "
                  "(%2.1f,%2.1f) --> %2.1f\n",
                  x,
                  y,
                  z,
                  cma_label_to_name(gcan_dst->labels[n]),
                  i,
                  DEGREES(fa[i]),
                  T1,
                  PD,
                  gc_dst->means[i]);
            *MATRIX_RELT(m_jac_dst, i + 1, 1) =
                dFlash_dT1(T1, PD, gca_flash_dst->TRs[i], gca_flash_dst->FAs[i], gca_flash_dst->TEs[i]);
            *MATRIX_RELT(m_jac_dst, i + 1, 2) =
                dFlash_dPD(T1, PD, gca_flash_dst->TRs[i], gca_flash_dst->FAs[i], gca_flash_dst->TEs[i]);

            //  *MATRIX_RELT(m_jac_src, i+1, 1) =
            //  dFlash_dT1(T1, PD, gca_flash_src->TRs[i],
            // gca_flash_src->FAs[i], gca_flash_src->TEs[i]) ;
            //  *MATRIX_RELT(m_jac_src, i+1, 2) =
            //  dFlash_dPD(T1, PD, gca_flash_src->TRs[i],
            // gca_flash_src->FAs[i], gca_flash_src->TEs[i]) ;
            if (gcan_dst->labels[n] == Ggca_label) {
              label_means[i] += gc_dst->means[i];
            }
          }
          // this allows src and dst have different ninputs
          for (i = 0; i < gca_flash_src->ninputs; i++) {
            *MATRIX_RELT(m_jac_src, i + 1, 1) =
                dFlash_dT1(T1, PD, gca_flash_src->TRs[i], gca_flash_src->FAs[i], gca_flash_src->TEs[i]);
            *MATRIX_RELT(m_jac_src, i + 1, 2) =
                dFlash_dPD(T1, PD, gca_flash_src->TRs[i], gca_flash_src->FAs[i], gca_flash_src->TEs[i]);
          }

#define MIN_T1 50
          if (x == Gx && y == Gy && z == Gz) {
            DiagBreak();
          }
          if (T1 < MIN_T1)
            m_cov_dst = MatrixIdentity(gca_flash_dst->ninputs, m_cov_dst);
          else {
            m_cov_src = load_covariance_matrix(gc_src, m_cov_src, gca_flash_src->ninputs);
            m_pinv_src = MatrixPseudoInverse(m_jac_src, m_pinv_src);
            m_cov_T1PD = MatrixSimilarityTransform(m_cov_src, m_pinv_src, m_cov_T1PD);
            m_cov_dst = MatrixSimilarityTransform(m_cov_T1PD, m_jac_dst, m_cov_dst);
          }
          for (v = i = 0; i < gca_flash_dst->ninputs; i++) {
            for (j = i; j < gca_flash_dst->ninputs; j++, v++) {
              gc_dst->covars[v] = *MATRIX_RELT(m_cov_dst, i + 1, j + 1);
            }
          }
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label < 0 || Ggca_label == gcan_src->labels[n])) {
            printf("predicted covariance matrix:\n");
            MatrixPrint(stdout, m_cov_dst);
          }
        }
      }
    }
  }

  if (Ggca_label >= 0) {
    printf("label %s (%d): means = ", cma_label_to_name(Ggca_label), Ggca_label);
    for (i = 0; i < gca_flash_dst->ninputs; i++) {
      label_means[i] /= (double)label_count;
      printf("%2.1f ", label_means[i]);
    }
    printf("\n");
  }

  /* check and fix singular covariance matrixces */
  GCAfixSingularCovarianceMatrices(gca_flash_dst);
  MatrixFree(&m_jac_dst);
  MatrixFree(&m_jac_src);
  if (m_cov_src) {
    MatrixFree(&m_cov_src);
  }
  MatrixFree(&m_cov_dst);
  if (m_pinv_src) {
    MatrixFree(&m_pinv_src);
  }
  if (m_cov_T1PD) {
    MatrixFree(&m_cov_T1PD);
  }

  gca_flash_dst->type = GCA_FLASH;
  GCAcopyDCToGCA(gca_flash_src, gca_flash_dst);
  return (gca_flash_dst);
}

int GCAnormalizeMeans(GCA *gca, float target)
{
  int x, y, z, frame, n;
  double norm;
  double val;
  GCA_NODE *gcan;
  GC1D *gc;

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          gc = &gcan->gcs[n];
          for (frame = 0, norm = 0; frame < gca->ninputs; frame++) {
            val = gc->means[frame];
            norm += (val * val);
          }
          norm = sqrt(norm) / target;
          if (FZERO(norm)) {
            norm = 1;
          }
          for (frame = 0; frame < gca->ninputs; frame++) {
            gc->means[frame] /= norm;
          }
        }
      }
    }
  }

  return (NO_ERROR);
}
int GCAregularizeCovariance(GCA *gca, float regularize)
{
  int x, y, z, i, r, c, n, num, nparams;
  GCA_NODE *gcan;
  GC1D *gc;
  double det, vars[MAX_GCA_INPUTS];
  MATRIX *m_cov = NULL;

  nparams = (gca->ninputs * (gca->ninputs + 1)) / 2 + gca->ninputs;
  /* covariance matrix and means */
  memset(vars, 0, sizeof(vars));
  for (num = 0, x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
            DiagBreak();
          }
          gc = &gcan->gcs[n];
          det = covariance_determinant(gc, gca->ninputs);
          if ((gc->ntraining == 0 && det > 1) || (gc->ntraining * gca->ninputs > 2.5 * nparams))
          /* enough to estimate parameters */
          {
            m_cov = load_covariance_matrix(gc, m_cov, gca->ninputs);
            for (r = 0; r < gca->ninputs; r++) {
              vars[r] += *MATRIX_RELT(m_cov, r + 1, r + 1);
            }
            num++;
          }
        }
      }
    }
  }
  if (m_cov) {
    MatrixFree(&m_cov);
  }
  if (num >= 1) {
    for (r = 0; r < gca->ninputs; r++) {
      vars[r] /= (float)num;
      printf("average std[%d] = %2.1f\n", r, sqrt(vars[r]));
    }
  }
  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
            DiagBreak();
          }
          gc = &gcan->gcs[n];
          for (i = r = 0; r < gca->ninputs; r++) {
            for (c = r; c < gca->ninputs; c++, i++) {
              gc->covars[i] = (1 - regularize) * gc->covars[i];
              if (r == c) {
                gc->covars[i] += regularize * vars[r];
              }
            }
          }
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z && (Ggca_label == gcan->labels[n] || Ggca_label < 0)) {
            MATRIX *m;
            printf(
                "regularizing covariance matrix "
                "for %s @ (%d, %d, %d):\n",
                cma_label_to_name(gcan->labels[n]),
                x,
                y,
                z);
            m = load_covariance_matrix(gc, NULL, gca->ninputs);
            MatrixPrint(stdout, m);
            MatrixFree(&m);
          }
        }
      }
    }
  }
  return (NO_ERROR);
}

MATRIX *GCAlabelCovariance(GCA *gca, int label, MATRIX *m_total)
{
  int xn, yn, zn, n;
  GCA_NODE *gcan;
  GC1D *gc;
  double wt;
  float prior;
  MATRIX *m_cov = NULL;

  /* compute overall white matter mean to use as anchor for rescaling */
  for (wt = 0.0, zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          /* find index in lookup table for this label */
          if (gcan->labels[n] != label) {
            continue;
          }
          gc = &gcan->gcs[n];
          prior = get_node_prior(gca, label, xn, yn, zn);
          if (prior != 0) {
            wt += prior;
            m_cov = load_covariance_matrix(gc, m_cov, gca->ninputs);
            MatrixScalarMul(m_cov, prior, m_cov);
            if (m_total == NULL) {
              m_total = MatrixCopy(m_cov, NULL);
            }
            else {
              MatrixAdd(m_cov, m_total, m_total);
            }
          }
        }
      }
    }
  }

  if (m_total == NULL) {
    return (NULL);
  }

  if (FZERO(wt)) {
    wt = 1;
  }
  MatrixScalarMul(m_total, 1 / wt, m_total);
  MatrixFree(&m_cov);
  return (m_total);
}

int GCAlabelMeanFromImage(GCA *gca, TRANSFORM *transform, MRI *mri, int label, float *means)
{
  int xn, yn, zn, n, r, xv, yv, zv;
  GCA_NODE *gcan;
  GC1D *gc;
  double wt, val;
  float prior;
  double MIN_MEAN_PRIOR = 0.5;

  /* compute overall white matter mean to use as anchor for rescaling */
  memset(means, 0, gca->ninputs * sizeof(float));
  for (wt = 0.0, zn = 0; zn < gca->node_depth; zn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        if (GCAnodeToSourceVoxel(gca, mri, transform, xn, yn, zn, &xv, &yv, &zv) == NO_ERROR) {
          for (n = 0; n < gcan->nlabels; n++) {
            /* find index in lookup table for this label */
            if (gcan->labels[n] != label) {
              continue;
            }
            gc = &gcan->gcs[n];
            prior = get_node_prior(gca, label, xn, yn, zn);
            if (prior < MIN_MEAN_PRIOR) {
              continue;
            }
            wt += prior;
            for (r = 0; r < gca->ninputs; r++) {
              MRIsampleVolumeFrame(mri, xv, yv, zv, r, &val);
              means[r] += val * prior;
              if (!std::isfinite(gc->means[r])) {
                DiagBreak();
              }
            }
          }
        }
      }
    }
  }
  for (r = 0; r < gca->ninputs; r++) {
    means[r] /= wt;
  }
  return (NO_ERROR);
}

/* don't try to estimate cortex directly - too hard to
   get alignment. We'll estimate it from other gm classes
   but how about just use initial linear registration?
*/
static int align_labels[] = {Left_Lateral_Ventricle,
                             Right_Lateral_Ventricle,
                             Right_Pallidum,
                             Left_Pallidum,
                             Right_Hippocampus,
                             Left_Hippocampus,
                             Right_Cerebral_White_Matter,
                             Left_Cerebral_White_Matter,
                             Left_Cerebral_Cortex,
                             Right_Cerebral_Cortex,
                             Right_Caudate,
                             Left_Caudate,
                             Left_Cerebellum_Cortex,
                             Right_Cerebellum_Cortex,
                             Left_Cerebellum_White_Matter,
                             Right_Cerebellum_White_Matter,
                             Left_Amygdala,
                             Right_Amygdala,
                             Left_Thalamus,
                             Right_Thalamus,
                             Left_Putamen,
                             Right_Putamen,
                             Brain_Stem,
                             Right_VentralDC,
                             Left_VentralDC,
                             Third_Ventricle,
                             Fourth_Ventricle};

#define NALIGN_LABELS (sizeof(align_labels) / sizeof(align_labels[0]))
#define BORDER_SIZE 2
#define WM_BORDER_SIZE 5

static int gm_labels[] = {
    Left_Hippocampus,
    Right_Hippocampus,
    Left_Amygdala,
    Right_Amygdala,
    //    Left_Caudate,
    //    Right_Caudate,
    Left_Cerebral_Cortex,
    Right_Cerebral_Cortex
    //          Left_Putamen,
    //  Right_Putamen
};

#define NGM_LABELS (sizeof(gm_labels) / sizeof(gm_labels[0]))

static int wm_labels[] = {Left_Cerebral_White_Matter, Right_Cerebral_White_Matter};
#define NWM_LABELS (sizeof(wm_labels) / sizeof(wm_labels[0]))

static int csf_labels[] = {Left_Lateral_Ventricle, Right_Lateral_Ventricle, Third_Ventricle, Fourth_Ventricle};

#define NCSF_LABELS (sizeof(csf_labels) / sizeof(csf_labels[0]))


static int lh_labels[] = {
    Left_Lateral_Ventricle,
    Left_Cerebral_White_Matter,
    Left_Hippocampus,
    Left_Cerebral_Cortex,
    Left_Caudate,
    Left_Cerebellum_Cortex,
    Left_Cerebellum_White_Matter,
    Left_Amygdala,
    Left_Thalamus,
    Left_Putamen,
    Left_Pallidum,
    Left_VentralDC,
};
static int rh_labels[] = {
    Right_Lateral_Ventricle,
    Right_Cerebral_White_Matter,
    Right_Hippocampus,
    Right_Cerebral_Cortex,
    Right_Caudate,
    Right_Cerebellum_Cortex,
    Right_Cerebellum_White_Matter,
    Right_Amygdala,
    Right_Thalamus,
    Right_Putamen,
    Right_Pallidum,
    Right_VentralDC,
};

#define NHEMI_LABELS (sizeof(rh_labels) / sizeof(rh_labels[0]))
int GCAmapRenormalizeWithAlignment(GCA *gca,
                                   MRI *mri,
                                   TRANSFORM *transform,
                                   FILE *logfp,
                                   const char *base_name,
                                   LTA **plta,
                                   int handle_expanded_ventricles)
{
  return (GCAcomputeRenormalizationWithAlignment(
      gca, mri, transform, logfp, base_name, plta, handle_expanded_ventricles, NULL, NULL, NULL, NULL));
}

int GCAmapRenormalizeWithAlignmentLongitudinal(GCA *gca,
                                               MRI *mri,
                                               TRANSFORM *transform,
                                               FILE *logfp,
                                               const char *base_name,
                                               LTA **plta,
                                               int handle_expanded_ventricles)
{
  return (GCAcomputeRenormalizationWithAlignmentLongitudinal(
      gca, mri, transform, logfp, base_name, plta, handle_expanded_ventricles, NULL, NULL, NULL, NULL));
}

/* Sequential routine adjusts the former scales,offsets,peaks */
int GCAseqRenormalizeWithAlignment(GCA *gca,
                                   MRI *mri,
                                   TRANSFORM *transform,
                                   FILE *logfp,
                                   const char *base_name,
                                   LTA **plta,
                                   int handle_expanded_ventricles,
                                   float *old_label_scales,
                                   float *old_label_offsets,
                                   float *old_label_peaks,
                                   int *old_label_computed)
{
  if (old_label_scales == NULL || old_label_offsets == NULL || old_label_peaks == NULL || old_label_computed == NULL)
    ErrorExit(ERROR_BADPARM,
              "%s ERROR:  GCA seqRenormalize old_label needs to "
              "be set for sequential call.\n",
              Progname);

  float plabel_scales[MAX_CMA_LABELS], plabel_offsets[MAX_CMA_LABELS], plabel_peaks[MAX_CMA_LABELS];
  int plabel_computed[MAX_CMA_LABELS];
  int ret = GCAcomputeRenormalizationWithAlignment(gca,
                                                   mri,
                                                   transform,
                                                   logfp,
                                                   base_name,
                                                   plta,
                                                   handle_expanded_ventricles,
                                                   plabel_scales,
                                                   plabel_offsets,
                                                   plabel_peaks,
                                                   plabel_computed);

  // merge plabel with old_label results
  int l;

  for (l = 0; l < MAX_CMA_LABELS; l++) {
    if (plabel_computed[l] != 0)  // we computed this label
    {
      if (old_label_computed[l] != 0)  // merge offsets and scales
      {
        // printf("%d   old: %f   new: %f   final: ",l,old_label_scales[l], plabel_scales[l]);
        old_label_offsets[l] = plabel_scales[l] * old_label_offsets[l] + plabel_offsets[l];
        old_label_scales[l] *= plabel_scales[l];
      }
      else  // no old label info -> coyp offsets and scales
      {
        old_label_offsets[l] = plabel_offsets[l];
        old_label_scales[l] = plabel_scales[l];
      }
      // in any case
      old_label_peaks[l] = plabel_peaks[l];
      old_label_computed[l] = plabel_computed[l];
      // printf("%f \n",old_label_scales[l]);
    }
    // else: if not computed in this call, just keep old_label info
  }

  if (logfp) {
    fprintf(logfp, " Labels after sequential adjustment:\n");
    for (l = 0; l < MAX_CMA_LABELS; l++)
      if (old_label_computed[l] != 0)
        fprintf(logfp,
                "label %s: scaling by %2.2f  + %2.1f to %2.0f\n",
                cma_label_to_name(l),
                old_label_scales[l],
                old_label_offsets[l],
                old_label_peaks[l]);
    fflush(logfp);
  }
  if (DIAG_VERBOSE_ON) {
    FILE *fp;
    fp = fopen("norm_offset.plt", "w");
    for (l = 0; l < MAX_CMA_LABELS; l++)
      if (old_label_computed[l] != 0) {
        fprintf(fp, "%d %f %f\n", l, old_label_scales[l], old_label_offsets[l]);
      }
    fclose(fp);
  }

  if (base_name) {
    FILE *fp;
    char fname[STRLEN];
    sprintf(fname, "%s.label_intensities.txt", base_name);
    printf("saving sequentially combined intensity scales to %s\n", fname);
    fp = fopen(fname, "w");
    if (fp == NULL) ErrorExit(ERROR_NOFILE, "%s: could not open intensity tracking file %s", Progname, fname);

    for (l = 0; l < MAX_CMA_LABELS; l++)
      if (old_label_computed[l] != 0)
        fprintf(fp,
                "%d %s %2.2f %2.1f %2.0f\n",
                l,
                cma_label_to_name(l),
                old_label_scales[l],
                old_label_offsets[l],
                old_label_peaks[l]);

    fflush(fp);
    fclose(fp);
    if (getenv("EXIT_AFTER_INT") != NULL) exit(0);
  }
  return ret;
}

static int mtl_labels[] = {Left_Amygdala, Right_Amygdala, Left_Hippocampus, Right_Hippocampus};
#define MTL_LABELS (sizeof(mtl_labels) / sizeof(mtl_labels[0]))

MRI *MRImarkPossibleBorders(MRI *mri, MRI *mri_borders, float border_thresh)
{
  int x, y, z, xk, yk, zk, xi, yi, zi;
  float val, val0;

  if (mri_borders == NULL) {
    mri_borders = MRIclone(mri, NULL);
  }

  for (x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }

        val0 = MRIgetVoxVal(mri, x, y, z, 0);
        for (xk = -1; xk <= 1; xk++) {
          xi = mri->xi[x + xk];
          for (yk = -1; yk <= 1; yk++) {
            yi = mri->yi[y + yk];
            for (zk = -1; zk <= 1; zk++) {
              if (fabs(xk) + fabs(yk) + fabs(zk) > 1) {
                continue;  // only 6-connected nbrs
              }
              zi = mri->zi[z + zk];
              val = MRIgetVoxVal(mri, xi, yi, zi, 0);
              if (fabs(val - val0) >= border_thresh) {
                MRIsetVoxVal(mri_borders, x, y, z, 0, 1);
                xk = yk = zk = 2;  // force breaking out of loop
                break;
              }
            }
          }
        }
      }
  return (mri_borders);
}
typedef struct
{
  int nsubjects;
  int nlabels;
  int ntypes;  // how many different acquisition types are there
  int *labels;
  int *types;                  // indicated the acquisition sequence for each subject
  double *similarities;        // correlation values
  double *avg_similarities;    // correlation values avrgd by acquisition type
  HISTOGRAM ***histos;         // one for each subject for each label
  HISTOGRAM **subject_histos;  // overall histo (not by label)
  HISTOGRAM **mri_histos;      // histos found in the volume being examined
  HISTOGRAM *mri_histo;
  char **subjects;
} LABEL_HISTOS, LH;

static LABEL_HISTOS *lhAlloc(int nsubjects, int nlabels, int ntypes)
{
  LABEL_HISTOS *lh;
  int s, n;

  lh = (LABEL_HISTOS *)calloc(1, sizeof(LABEL_HISTOS));

  lh->nsubjects = nsubjects;
  lh->nlabels = nlabels;
  lh->labels = (int *)calloc(nsubjects, sizeof(int));
  lh->ntypes = ntypes;
  lh->avg_similarities = (double *)calloc(ntypes, sizeof(double));
  lh->histos = (HISTOGRAM ***)calloc(nsubjects, sizeof(HISTOGRAM **));
  lh->subject_histos = (HISTOGRAM **)calloc(nsubjects, sizeof(HISTOGRAM *));
  lh->subjects = (char **)calloc(nsubjects, sizeof(char *));
  lh->types = (int *)calloc(nsubjects, sizeof(int));
  lh->similarities = (double *)calloc(nsubjects, sizeof(double));

  lh->mri_histos = (HISTOGRAM **)calloc(nlabels, sizeof(HISTOGRAM *));
  for (n = 0; n < lh->nlabels; n++) {
    lh->mri_histos[n] = HISTOinit(NULL, 256, 0, 255);
  }
  lh->mri_histo = HISTOinit(NULL, 256, 0, 255);
  for (s = 0; s < nsubjects; s++) {
    lh->subject_histos[s] = HISTOinit(NULL, 256, 0, 255);
    lh->histos[s] = (HISTOGRAM **)calloc(nlabels, sizeof(HISTOGRAM *));
    for (n = 0; n < lh->nlabels; n++) {
      lh->histos[s][n] = HISTOinit(NULL, 256, 0, 255);
    }
    lh->subjects[s] = (char *)calloc(100, sizeof(char));
  }

  return (lh);
}

static int lhComputeMRIHistos(LH *lh, MRI *mri_norm, int nerode)
{
  MRI *mri_label, *mri_eroded;
  //  int   n ;

  mri_label = MRIclone(mri_norm, NULL);
  mri_eroded = MRIclone(mri_norm, NULL);
  MRIbinarize(mri_norm, mri_label, 1, 0, 1);
  HISTOfree(&lh->mri_histo);
  lh->mri_histo = MRIhistogramLabel(mri_norm, mri_label, 1, 256);
  MRIfree(&mri_label);
  MRIfree(&mri_eroded);
  return (NO_ERROR);
}

static double lhComputeSimilarity(LH *lh, int s)
{
  double total_sim;
  total_sim = HISTOcorrelate(lh->mri_histo, lh->subject_histos[s]);

  return (total_sim);
}

#define MAX_SUBJECTS 10000
static int lhFindMostSimilarAcquisitionType(LH *lh, MRI *mri_norm)
{
  int max_ind, s, nsubjects[MAX_SUBJECTS], n;
  double max_sim, sim;

  memset(nsubjects, 0, sizeof(nsubjects));
  lhComputeMRIHistos(lh, mri_norm, 1);
  max_sim = lhComputeSimilarity(lh, max_ind = 0);
  printf("subject %s (%d): sim = %f\n", lh->subjects[0], 0, max_sim);
  for (s = 1; s < lh->nsubjects; s++) {
    sim = lhComputeSimilarity(lh, s);
    printf("type %d, subject %s, (%d): sim = %f\n", lh->types[s], lh->subjects[s], s, sim);
    lh->avg_similarities[lh->types[s]] += sim;
    nsubjects[lh->types[s]]++;
    if (sim > max_sim) {
      max_sim = sim;
      max_ind = s;
    }
  }
  max_ind = 0;
  max_sim = 0;
  for (n = 0; n < lh->ntypes; n++) {
    lh->avg_similarities[n] /= nsubjects[n];
    printf("type %d (%d subjects): avg corr = %2.3f\n", n, nsubjects[n], lh->avg_similarities[n]);
    if (lh->avg_similarities[n] > max_sim) {
      max_sim = lh->avg_similarities[n];
      max_ind = n;
    }
  }
  printf("best matching acquisition = %d\n", max_ind);
  return (max_ind);
}

static int lhFree(LH **plh)
{
  LH *lh;
  int s, n;

  lh = *plh;
  *plh = NULL;
  free(lh->labels);
  for (s = 0; s < lh->nsubjects; s++) {
    free(lh->subjects[s]);
    for (n = 0; n < lh->nlabels; n++) {
      HISTOfree(&lh->histos[s][n]);
    }
    free(lh->histos[s]);
  }
  HISTOfree(&lh->mri_histo);
  for (n = 0; n < lh->nlabels; n++) {
    HISTOfree(&lh->mri_histos[n]);
  }
  free(lh->mri_histos);
  free(lh->histos);
  free(lh->subjects);
  free(lh);
  return (NO_ERROR);
}

static LABEL_HISTOS *lhRead(const char *fname)
{
  FILE *fp;
  int nsubjects, nlabels, s, n, cnt, i, label, ntypes;
  LABEL_HISTOS *lh;
  char *cp, line[2 * STRLEN];

  fp = fopen(fname, "r");
  if (fp == NULL) {
    ErrorExit(ERROR_NOFILE, "lhRead: could not open file %s", fname);
  }

  cp = fgetl(line, STRLEN - 1, fp);
  sscanf(cp, "nsubjects = %d", &nsubjects);
  cp = fgetl(line, STRLEN - 1, fp);
  sscanf(cp, "ntypes = %d", &ntypes);
  cp = fgetl(line, STRLEN - 1, fp);
  sscanf(cp, "nlabels = %d", &nlabels);
  printf("reading label histos for %d subjects with %d labels each\n", nsubjects, nlabels);

  lh = lhAlloc(nsubjects, nlabels, ntypes);
  cp = fgetl(line, STRLEN - 1, fp);
  cp = strtok(cp, " ");
  for (n = 0; n < lh->nlabels; n++) {
    sscanf(cp, "%d", &lh->labels[n]);
    cp = strtok(NULL, " ");
  }

  for (s = 0; s < lh->nsubjects; s++) {
    cp = fgetl(line, STRLEN - 1, fp);
    sscanf(cp, "subject = %s type = %d", lh->subjects[s], &lh->types[s]);
    //    printf("reading subject %d: %s, acquisition type = %d\n", s, lh->subjects[s], lh->types[s]) ;

    // read in global histogram
    cp = fgetl(line, 2 * STRLEN - 1, fp);
    for (i = 0; i < 256; i++) {
      sscanf(cp, "%d", &cnt);
      lh->subject_histos[s]->counts[i] = cnt;
      lh->subject_histos[s]->bins[i] = i;
      cp = strchr(cp, ' ');
      if (cp == NULL) {
        break;
      }
      cp++;  // past space
    }
    for (n = 0; n < nlabels; n++) {
      if (s == 3 && n == 35) {
        DiagBreak();
      }
      cp = fgetl(line, STRLEN - 1, fp);
      sscanf(cp, "label = %d", &label);
      cp = fgetl(line, 2 * STRLEN - 1, fp);
      for (i = 0; i < 256; i++) {
        if (i == 110 && s == 0 && n == 0) {
          DiagBreak();
        }
        sscanf(cp, "%d", &cnt);
        lh->histos[s][n]->counts[i] = cnt;
        lh->histos[s][n]->bins[i] = i;
        cp = strchr(cp, ' ');
        if (cp == NULL) {
          break;
        }
        cp++;  // skip space
      }
      HISTOmakePDF(lh->histos[s][n], lh->histos[s][n]);
    }
  }

  fclose(fp);
  return (lh);
}

int GCAmapRenormalizeWithHistograms(GCA *gca,
                                    MRI *mri,
                                    TRANSFORM *transform,
                                    FILE *logfp,
                                    const char *base_name,
                                    float *plabel_scales,
                                    float *plabel_offsets,
                                    float *plabel_peaks,
                                    int *plabel_computed)
{
  int l, frame, j, gca_peak, acq_type, s, num;
  float label_scales[MAX_CMA_LABELS], label_peaks[MAX_CMA_LABELS], label_offsets[MAX_CMA_LABELS];
  // float lower_thresh, upper_thresh;
  float scale, offset, avg_scale, avg_offset;
  LH *lh;
  int equiv_class[MAX_CMA_LABELS];
  HISTOGRAM *h_gca;

  if (FileExists("label_scales.dat")) {
    lh = lhRead("label_scales.dat");
  }
  else {
    char *cp = getenv("FREESURFER_HOME"), fname[STRLEN];

    sprintf(fname, "%s/average/label_scales.dat", cp);
    lh = lhRead(fname);
  }

  acq_type = lhFindMostSimilarAcquisitionType(lh, mri);
  if (plabel_scales == NULL) {
    plabel_scales = label_scales;
  }
  if (plabel_offsets == NULL) {
    plabel_offsets = label_offsets;
  }
  if (plabel_peaks == NULL) {
    plabel_peaks = label_peaks;
  }

  memset(label_peaks, 0, sizeof(label_peaks));

  set_equilavent_classes(equiv_class);

  printf("renormalizing by  histogram matching....\n");
  for (frame = 0; frame < mri->nframes; frame++) {
    for (l = 0; l < lh->nlabels; l++) {
      if (l == Gdiag_no) {
        DiagBreak();
      }
      label_scales[l] = 1.0;
      label_offsets[l] = 0.0;
    }

    printf("renormalizing input #%d\n", frame);
    for (j = 0; j < lh->nlabels; j++) {
      l = lh->labels[j];
      if (l == Gdiag_no) {
        DiagBreak();
      }

      h_gca = gcaGetLabelHistogram(gca, l, 0, 0);
      HISTOmakePDF(h_gca, h_gca);
      for (avg_scale = avg_offset = 0.0, num = s = 0; s < lh->nsubjects; s++) {
        if (lh->types[s] != acq_type) {
          continue;
        }
        num++;
        HISTOfindLinearFit(h_gca, lh->histos[s][j], .025, 4, 0, 0, &scale, &offset);
        avg_scale += scale;
        avg_offset += offset;
      }

      avg_scale /= num;
      avg_offset /= num;
      label_offsets[l] = avg_offset;
      label_scales[l] = avg_scale;

      gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
      if (gca_peak < 0) {
        // fprintf(stderr,
        //      "INFO: GCAmapRenormalizeWithAlignment: "
        //      "gca peak(=%d) < 0\n",gca_peak);
        // fflush(stderr);
        continue;
      }
      if (gca_peak >= h_gca->nbins) {
        fprintf(stderr,
                "ERROR: GCAmapRenormalizeWithAlignment: "
                "gca peak(=%d) >= h_gca->nbins(=%d)\n",
                gca_peak,
                h_gca->nbins);
        fflush(stderr);
        exit(1);
      }
      if (gca_peak >= 0) {
        printf("%s: gca peak = %2.5f (%2.0f), scaling to %2.1f\n",
               cma_label_to_name(l),
               h_gca->counts[gca_peak],
               h_gca->bins[gca_peak],
               h_gca->bins[gca_peak] * label_scales[l] + label_offsets[l]);
      }

      label_peaks[l] = h_gca->bins[gca_peak];
      fflush(stdout);

      // switch (l) {
      //   case Brain_Stem:
      //   case Left_VentralDC:
      //   case Right_VentralDC:
      //     lower_thresh = 80;
      //     upper_thresh = 110;
      //     break;
      //   case Left_Caudate:
      //   case Right_Caudate:
      //     lower_thresh = 50;
      //     upper_thresh = 100;
      //     break;
      //   case Left_Cerebral_Cortex:
      //   case Right_Cerebral_Cortex:
      //     lower_thresh = 40;
      //     upper_thresh = 95;
      //     break;
      //   case Left_Pallidum:
      //   case Right_Pallidum:
      //     lower_thresh = 75;
      //     upper_thresh = 135;
      //     break;
      //   case Left_Thalamus:
      //   case Right_Thalamus:
      //     lower_thresh = 75;
      //     upper_thresh = 120;
      //     break;
      //   case Left_Amygdala:
      //   case Right_Amygdala:
      //     lower_thresh = 50;
      //     upper_thresh = 90;
      //     break;
      //   case Left_Cerebral_White_Matter:
      //   case Right_Cerebral_White_Matter:
      //     lower_thresh = 90;
      //     upper_thresh = 130;
      //     break;
      //   case Left_Putamen:
      //   case Right_Putamen:
      //     lower_thresh = 60;
      //     upper_thresh = 107;
      //     break;
      //   case Left_Lateral_Ventricle:
      //   case Right_Lateral_Ventricle:
      //   case Third_Ventricle:
      //   case Fourth_Ventricle:
      //   case CSF:
      //     lower_thresh = 0;
      //     upper_thresh = 55;
      //     break;
      //   case Left_Inf_Lat_Vent:
      //   case Right_Inf_Lat_Vent:
      //     lower_thresh = 0;
      //     upper_thresh = 65;
      //     break;
      //   default:
      //     lower_thresh = 0;
      //     upper_thresh = 256;
      //     break;
      // }
      fflush(stdout);
    }
  }

  if (logfp) {
    for (l = 0; l < MAX_CMA_LABELS; l++)
      fprintf(logfp,
              "label %s: scaling by %2.2f  + %2.1f to %2.0f\n",
              cma_label_to_name(l),
              label_scales[l],
              label_offsets[l],
              label_peaks[l]);
    fflush(logfp);
  }
  if (DIAG_VERBOSE_ON) {
    FILE *fp;
    fp = fopen("norm_offset.plt", "w");
    for (l = 0; l < MAX_CMA_LABELS; l++) {
      fprintf(fp, "%d %f %f\n", l, label_scales[l], label_offsets[l]);
    }
    fclose(fp);
  }

  if (base_name) {
    FILE *fp;
    char fname[STRLEN];
    sprintf(fname, "%s.label_intensities.txt", base_name);
    printf("saving intensity scales to %s\n", fname);
    fp = fopen(fname, "w");
    if (fp == NULL) {
      ErrorExit(ERROR_NOFILE, "%s: could not open intensity tracking file %s", Progname, fname);
    }

    for (l = 0; l < MAX_CMA_LABELS; l++)
      fprintf(
          fp, "%d %s %2.2f %2.1f %2.0f\n", l, cma_label_to_name(l), label_scales[l], label_offsets[l], label_peaks[l]);

    fflush(fp);
    fclose(fp);
    if (getenv("EXIT_AFTER_INT") != NULL) {
      exit(0);
    }
  }
  gcaCheck(gca);
  for (frame = 0; frame < mri->nframes; frame++) {
    GCAapplyRenormalization(gca, label_scales, label_offsets, frame);
  }
  gcaCheck(gca);
  lhFree(&lh);
  return (NO_ERROR);
}
int GCAcomputeRenormalizationWithAlignment(GCA *gca,
                                           MRI *mri,
                                           TRANSFORM *transform,
                                           FILE *logfp,
                                           const char *base_name,
                                           LTA **plta,
                                           int handle_expanded_ventricles,
                                           float *plabel_scales,
                                           float *plabel_offsets,
                                           float *plabel_peaks,
                                           int *plabel_computed)
{
  HISTOGRAM *h_mri, *h_gca, *h_mtl = NULL, *h_caudate;
  unsigned int k, j;
  int l, nbins, i, x, y, z, num, frame, bin, n, computed[MAX_CMA_LABELS], b, label, border = BORDER_SIZE,
                                                                                          gca_peak, mri_peak;
  float fmin, fmax, label_scales[MAX_CMA_LABELS], overlap, mean_gm_scale, mean_wm_scale, mean_csf_scale,
      label_peaks[MAX_CMA_LABELS], label_offsets[MAX_CMA_LABELS], mean_wm_offset, mean_csf_offset, mean_gm_offset,
      lower_thresh, upper_thresh;
  double val /*, scale*/;
  MRI *mri_seg = NULL, *mri_aligned, *mri_labels = NULL, *mri_borders;
  char fname[STRLEN];
  MATRIX *m_L, *m_by_label[MAX_CMA_LABELS];
  LTA *lta;

  if (plabel_scales == NULL) {
    plabel_scales = label_scales;
  }
  if (plabel_offsets == NULL) {
    plabel_offsets = label_offsets;
  }
  if (plabel_peaks == NULL) {
    plabel_peaks = label_peaks;
  }
  if (plabel_computed == NULL) {
    plabel_computed = computed;
  }

  double det = -1;
  float peak_threshold = 0.03;
  float overlap_threshold = 0.001;
  int equiv_class[MAX_CMA_LABELS];

  memset(label_peaks, 0, sizeof(label_peaks));
  if (transform->type == MORPH_3D_TYPE) {
    peak_threshold = 0.01;
    overlap_threshold = -1.0;  // at mri_ca_label stage;
    // trust the registration more
  }

  set_equilavent_classes(equiv_class);

  printf("renormalizing by structure alignment....\n");
  if (plta) {
    lta = *plta;
  }
  else {
    lta = NULL;
  }
  mri_borders = MRImarkPossibleBorders(mri, NULL, 15);
  for (frame = 0; frame < mri->nframes; frame++) {
    for (l = 0; l < MAX_CMA_LABELS; l++) {
      if (l == Gdiag_no) {
        DiagBreak();
      }
      label_scales[l] = 1.0;
      label_offsets[l] = 0.0;
      computed[l] = 0;
      m_by_label[l] = NULL;  // not estimated yet
    }

    printf("renormalizing input #%d\n", frame);
    MRIvalRangeFrame(mri, &fmin, &fmax, frame);
    nbins = 256;
    h_mri = HISTOalloc(nbins);
    h_mtl = HISTOalloc(nbins);
    h_caudate = HISTOalloc(nbins);
    for (j = 0; j < NALIGN_LABELS; j++) {
      l = align_labels[j];
      if (l == Gdiag_no) {
        DiagBreak();
      }

      mri_seg = MRIclone(mri, mri_seg);
      mri_labels = MRIclone(mri, mri_labels);

      /* include 2 voxel border to get context around structure.
        e.g. hippo is made easier to find by wm inferior and ventricle
        posterior.
      */
      if (transform->type != MORPH_3D_TYPE) {
        if (IS_HIPPO(l) || IS_AMYGDALA(l) || IS_CAUDATE(l) || IS_PUTAMEN(l) || IS_PALLIDUM(l)) {
          border = BORDER_SIZE + 1;  // need more context for hippo
        }
        else if (IS_GM(l)) {
          border = 0;
        }
        else if (IS_WHITE_CLASS(l)) {
          border = WM_BORDER_SIZE;
        }
        else {
          border = BORDER_SIZE;
        }
        GCAbuildMostLikelyVolumeForStructure(gca, mri_seg, l, border, transform, mri_labels);
        for (x = 0; x < mri_labels->width; x++) {
          for (y = 0; y < mri_labels->height; y++) {
            for (z = 0; z < mri_labels->depth; z++) {
              if (x == Gx && y == Gy && z == Gz) {
                DiagBreak();
              }
              label = MRIgetVoxVal(mri_labels, x, y, z, 0);
              if (computed[label] == 0) {
                continue;
              }
              val = MRIgetVoxVal(mri_seg, x, y, z, frame);
              val = val * label_scales[label] + label_offsets[label];
              MRIsetVoxVal(mri_seg, x, y, z, frame, val);
            }
          }
        }

        /* ventricle at the posterior part of hippo
          frequently makes local minima
          in alignment energy functional - remove them.
        */
        if (l == Left_Hippocampus || l == Right_Hippocampus) {
          for (x = 0; x < mri_labels->width; x++) {
            for (y = 0; y < mri_labels->height; y++) {
              for (z = 0; z < mri_labels->depth; z++) {
                if (x == Gx && y == Gy && z == Gz) {
                  DiagBreak();
                }
                label = MRIgetVoxVal(mri_labels, x, y, z, 0);
                if (IS_LAT_VENT(label) || IS_INF_LAT_VENT(label)) {
                  MRIsetVoxVal(mri_seg, x, y, z, frame, 0);
                }
              }
            }
          }
        }
        if (l == Left_Cerebral_White_Matter || l == Right_Cerebral_White_Matter) {
          MRI *mri_tmp, *mri_border;

          // create a volume that is the wm eroded 3 times,
          // plus the border voxels
          mri_tmp = MRIclone(mri_seg, NULL);
          GCAbuildMostLikelyVolumeForStructure(gca, mri_tmp, l, BORDER_SIZE, transform, NULL);
          mri_border = MRIsubtract(mri_seg, mri_tmp, NULL);
          // just outside

          // erode just the interior 4 times to get to high prob regions
          MRIerode(mri_tmp, mri_tmp);
          MRIerode(mri_tmp, mri_tmp);
          MRIerode(mri_tmp, mri_tmp);
          MRIerode(mri_tmp, mri_tmp);
          MRIerode(mri_tmp, mri_tmp);
          MRIerode(mri_tmp, mri_tmp);            // two more to remove border
          MRIadd(mri_tmp, mri_border, mri_seg);  // border + interior
          MRIfree(&mri_tmp);
          MRIfree(&mri_border);
        }

        if (Gdiag & DIAG_WRITE) {
          sprintf(fname, "%s_label%d.mgz", base_name, l);
          MRIwrite(mri_seg, fname);
        }

        if (lta)  // try to find a previously computed one
        {
          for (n = 0; n < lta->num_xforms; n++)
            if (lta->xforms[n].label == l) {
              break;
            }
          if (n >= lta->num_xforms) {
            n = -1;  // indicate no xform found
          }
        }
        else  // no transform specified by caller
        {
          n = -1;
        }
        if (n < 0)  // no transform - compute one
        {
          // float  evalues[4] ;
          // MATRIX *m_evectors ;

          printf("aligning %s\n", cma_label_to_name(l));
          m_L = MRIgetVoxelToVoxelXform(mri_seg, mri);
          if (!IS_GM(l)) {
            // will use alignment of WM for GM
            //  MRIpowellAlignImages(mri_seg, mri,
            // m_L, &scale_factor, NULL,NULL, NULL, 0) ;
            if ((l == Left_Lateral_Ventricle || l == Right_Lateral_Ventricle) && (transform->type != MORPH_3D_TYPE) &&
                (handle_expanded_ventricles == 1)) {
              char label_base_name[STRLEN];
              sprintf(label_base_name, "%s_label%d", base_name, l);
              initialize_ventricle_alignment(mri_seg, mri, m_L, label_base_name, border, l);
            }
            MRIfaridAlignImages(mri_seg, mri, m_L);
          }
          else {
            // assume that cortical gm goes as wm
            if ((l == Left_Cerebral_Cortex) && computed[Left_Cerebral_White_Matter] != 0) {
              if (m_by_label[Left_Cerebral_White_Matter]) m_L = MatrixCopy(m_by_label[Left_Cerebral_White_Matter], m_L);
            }
            if ((l == Right_Cerebral_Cortex) && computed[Right_Cerebral_White_Matter] != 0) {
              if (m_by_label[Right_Cerebral_White_Matter])
                m_L = MatrixCopy(m_by_label[Right_Cerebral_White_Matter], m_L);
            }
          }

          det = MatrixDeterminant(m_L);
          if (det > 4 && det < 8 && (l == Left_Lateral_Ventricle || l == Right_Lateral_Ventricle)) {
            det = 1;  // allow large determinants for the ventricles
          }


          printf("det = %2.3f, M=\n", det);
          MatrixPrint(stdout, m_L);
          if (det < 0.25 || det > 4) {
            printf("invalid transform detected (det=%2.4f) \n", det);
            det = -1;  // mark it as invalid for later
            MatrixFree(&m_L);
            m_L = MRIgetVoxelToVoxelXform(mri_seg, mri);
          }
        }
        else  // use previously computed transform
        {
          m_L = MatrixCopy(lta->xforms[n].m_L, NULL);
        }

        if (l == Gdiag_no) {
          DiagBreak();
        }

        if (Gdiag & DIAG_WRITE) {
          sprintf(fname, "%s_label%d_after.mgz", base_name, l);
          mri_aligned = MRIlinearTransform(mri_seg, NULL, m_L);
          MRIwrite(mri_aligned, fname);
          MRIfree(&mri_aligned);
        }

        if (l == Left_Cerebral_White_Matter || l == Right_Cerebral_White_Matter) {
          // wm so big it's hard to localize with a linear xform
          GCAbuildMostLikelyVolumeForStructure(gca, mri_seg, l, 0, transform, NULL);
          MRIerode(mri_seg, mri_seg);
          MRIerode(mri_seg, mri_seg);
        }
        else {
          /* put ventricles back in for erosion to remove
            (otherwise a bunch of hippo
            gets removed */
          if (l == Left_Hippocampus || l == Right_Hippocampus) {
            for (x = 0; x < mri_labels->width; x++) {
              for (y = 0; y < mri_labels->height; y++) {
                for (z = 0; z < mri_labels->depth; z++) {
                  if (x == Gx && y == Gy && z == Gz) {
                    DiagBreak();
                  }
                  label = MRIgetVoxVal(mri_labels, x, y, z, 0);
                  if (IS_LAT_VENT(label) || IS_INF_LAT_VENT(label)) {
                    MRIsetVoxVal(mri_seg, x, y, z, frame, 128);
                  }
                }
              }
            }
          }
          for (b = 0; b < border; b++) {
            MRIerode(mri_seg, mri_seg);  // get rid of outside border
          }
          MRIerode(mri_seg, mri_seg);  // get rid of inside border
        }

        mri_aligned = MRIlinearTransform(mri_seg, NULL, m_L);
      }
      else  // 3d morph already done - don't bother aligning
      {
        m_L = NULL;
        if (l == Left_Cerebral_White_Matter || l == Right_Cerebral_White_Matter) {
          // wm so big it's hard to localize with a linear xform
          GCAbuildMostLikelyVolumeForStructure(gca, mri_seg, l, 0, transform, NULL);
          MRIerode(mri_seg, mri_seg);
        }
        else {
          GCAbuildMostLikelyVolumeForStructure(gca, mri_seg, l, 0, transform, NULL);
        }
        mri_aligned = MRIerode(mri_seg, NULL);
      }

      MRIbinarize(mri_aligned, mri_aligned, 1, 0, 128);
      if (Gdiag & DIAG_WRITE) {
        sprintf(fname, "%s_label%d_eroded.mgz", base_name, l);
        MRIwrite(mri_aligned, fname);
      }
      if (l == Gdiag_no) {
        DiagBreak();
      }
      HISTOclear(h_mri, h_mri);
      h_mri->bin_size = (fmax - fmin) / 255.0;
      if (h_mri->bin_size < 1 && (mri->type == MRI_UCHAR || mri->type == MRI_SHORT)) {
        h_mri->bin_size = 1;
      }
      for (i = 0; i < nbins; i++) {
        h_mri->bins[i] = (i + 1) * h_mri->bin_size;
      }

      for (num = x = 0; x < mri_aligned->width; x++) {
        for (y = 0; y < mri_aligned->height; y++) {
          for (z = 0; z < mri_aligned->depth; z++) {
            if (x == Gx && y == Gy && z == Gz) {
              DiagBreak();
            }
            MRIsampleVolumeFrame(mri_aligned, x, y, z, frame, &val);
            if (DZERO(val))  // not in this structure
            {
              continue;
            }
            MRIsampleVolumeFrame(mri, x, y, z, frame, &val);

            if (FZERO(val))  // skull stripped
            {
              continue;
            }
            if (MRIgetVoxVal(mri_borders, x, y, z, 0) > 0)  // avoid border voxels - matches erosion
            {
              continue;
            }
            bin = nint((val - fmin) / h_mri->bin_size);
            if (bin >= h_mri->nbins) {
              bin = h_mri->nbins - 1;
            }
            else if (bin < 0) {
              bin = 0;
            }

            h_mri->counts[bin]++;
            num++;
          }
        }
      }
      if (l == Gdiag_no) {
        DiagBreak();
      }
      MRIfree(&mri_aligned);

      h_gca = gcaGetLabelHistogram(gca, l, 0, 0);
      gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
      HISTOmakePDF(h_gca, h_gca);

      if (gca_peak < 0) {
        // fprintf(stderr,
        //      "INFO: GCAmapRenormalizeWithAlignment: "
        //      "gca peak(=%d) < 0\n",gca_peak);
        // fflush(stderr);
        continue;
      }
      if (gca_peak >= h_gca->nbins) {
        fprintf(stderr,
                "ERROR: GCAmapRenormalizeWithAlignment: "
                "gca peak(=%d) >= h_gca->nbins(=%d)\n",
                gca_peak,
                h_gca->nbins);
        fflush(stderr);
        exit(1);
      }
      if (gca_peak >= 0) {
        printf("gca peak = %2.5f (%2.0f)\n", h_gca->counts[gca_peak], h_gca->bins[gca_peak]);
      }

      label_peaks[l] = h_gca->bins[gca_peak];
      fflush(stdout);

      switch (l) {
        case Brain_Stem:
        case Left_VentralDC:
        case Right_VentralDC:
          lower_thresh = 80;
          upper_thresh = 110;
          break;
        case Left_Caudate:
        case Right_Caudate:
          lower_thresh = 50;
          upper_thresh = 100;
          break;
        case Left_Cerebral_Cortex:
        case Right_Cerebral_Cortex:
          lower_thresh = 40;
          upper_thresh = 95;
          break;
        case Left_Pallidum:
        case Right_Pallidum:
          lower_thresh = 75;
          upper_thresh = 135;
          break;
        case Left_Thalamus:
        case Right_Thalamus:
          lower_thresh = 75;
          upper_thresh = 120;
          break;
        case Left_Hippocampus:
        case Right_Hippocampus:
        case Left_Amygdala:
        case Right_Amygdala:
          lower_thresh = 50;
          upper_thresh = 90;
          break;
        case Left_Cerebral_White_Matter:
        case Right_Cerebral_White_Matter:
          lower_thresh = 90;
          upper_thresh = 130;
          break;
        case Left_Putamen:
        case Right_Putamen:
          lower_thresh = 60;
          upper_thresh = 107;
          break;
        case Left_Lateral_Ventricle:
        case Right_Lateral_Ventricle:
        case Third_Ventricle:
        case Fourth_Ventricle:
        case CSF:
          lower_thresh = 0;
          upper_thresh = 55;
          break;
        case Left_Inf_Lat_Vent:
        case Right_Inf_Lat_Vent:
          lower_thresh = 0;
          upper_thresh = 65;
          break;
        default:
          lower_thresh = 0;
          upper_thresh = 256;
          break;
      }
      HISTOclearBins(h_mri, h_mri, 0, lower_thresh - 1);
      HISTOclearBins(h_mri, h_mri, upper_thresh + 1, 255);
      mri_peak = HISTOfindHighestPeakInRegion(h_mri, 0, h_mri->nbins);
      HISTOfillHoles(h_mri);
      if (l == Left_Caudate || l == Right_Caudate) {
        h_caudate->bin_size = 1;
        h_caudate->min = 0;
        h_caudate->max = 255;
        for (i = 0; i < nbins; i++) {
          h_caudate->bins[i] = (i + 1) * h_caudate->bin_size;
        }
        HISTOadd(h_mri, h_caudate, h_caudate);
      }

      if (l == Left_Hippocampus || l == Right_Hippocampus || l == Left_Amygdala || l == Right_Amygdala) {
        h_mtl->bin_size = 1;
        h_mtl->min = 0;
        h_mtl->max = 255;
        for (i = 0; i < nbins; i++) {
          h_mtl->bins[i] = (i + 1) * h_mtl->bin_size;
        }
        HISTOadd(h_mri, h_mtl, h_mtl);
      }
      HISTOmakePDF(h_mri, h_mri);
      if (mri_peak >= 0) printf("mri peak = %2.5f (%2.0f)\n", h_mri->counts[mri_peak], h_mri->bins[mri_peak]);
      fflush(stdout);

      if (IS_CSF(l) && h_mri->bins[mri_peak] > 55) {
        printf("CSF peak too bright - rejecting\n");
        continue;
      }
      if (h_mri->counts[mri_peak] < peak_threshold || num <= 50)
      /* not enough to reliably estimate density */
      {
        if (h_mri->counts[mri_peak] < peak_threshold)
          printf(
              "uniform distribution in MR - "
              "rejecting arbitrary fit\n");
        if (m_L) {
          MatrixFree(&m_L);
        }
        continue;
      }
      if (m_L) {
        if (plta && (!IS_GM(l)))  // GM will be copied from WM later
        {
          m_by_label[l] = m_L;  // store if for assembling an LTA later
        }
        else {
          MatrixFree(&m_L);
        }
      }

      if (Gdiag & DIAG_WRITE) {
        sprintf(fname, "%s_label%d_mri.plt", base_name, l);
        HISTOplot(h_mri, fname);
        sprintf(fname, "%s_label%d_gca.plt", base_name, l);
        HISTOplot(h_gca, fname);
        DiagBreak();
      }
      overlap = HISTOthreshSum(h_mri, h_gca, .025);

      // if (overlap > 0.01)
      //      if (overlap > 0.001)
      if (IS_LAT_VENT(l) || overlap > overlap_threshold) {
        //                        if (l == Gdiag_no)
        //  HISTOfindLinearFit(h_gca, h_mri, .025, 10, -75, 75,
        // &label_scales[l],  &label_offsets[l]) ;
        //              HISTOfindLinearFit(h_gca, h_mri, .025,
        // 4, -125, 125, &label_scales[l], &label_offsets[l]) ;
        HISTOfindLinearFit(h_gca, h_mri, .025, 4, 0, 0, &label_scales[l], &label_offsets[l]);

        val = h_gca->bins[gca_peak] * label_scales[l] + label_offsets[l];
        if ((val < lower_thresh || val > upper_thresh) ||
            (h_mri->bins[mri_peak] < lower_thresh || h_mri->bins[mri_peak] > upper_thresh)) {
          //          if (transform->type != MORPH_3D_TYPE)
          {
            printf(
                "%s: unreasonable value (%2.1f/%2.1f), "
                "not in range [%2.0f, %2.0f] - rejecting\n",
                cma_label_to_name(l),
                val,
                h_mri->bins[mri_peak],
                lower_thresh,
                upper_thresh);
            label_scales[l] = 1.0;
            label_offsets[l] = 1.0;
            continue;
          }
        }

        // only allow certain labels to be used for initializing the 3d morph
        // (which is what happens when computed[l] = 2)
        computed[l] = (det > 0) ? 2 : 1;
        printf(
            "%s (%d): linear fit = %2.2f x + %2.1f "
            "(%d voxels, overlap=%2.3f)\n",
            cma_label_to_name(l),
            l,
            label_scales[l],
            label_offsets[l],
            num,
            overlap);

        // note that the following range need be changed
        // if both scale and offset are allowed' 1/1.5 = 0.67
        if (IS_LAT_VENT(l)) {
          if (label_scales[l] < 0.4) {
            label_scales[l] = 0.4;
          }
          else if (label_scales[l] > 1.5) {
            label_scales[l] = 1.5;
          }
        }
        if ((label_scales[l] < 0.67 || (label_scales[l] > 1.5)) && !IS_LAT_VENT(l)) {
          /*
           if(IS_CSF(l)){
           if(label_scales[l] < 0.67) label_scales[l] = 0.67;
           else if(label_scales[l] > 1.5) label_scales[l] = 1.5;
           } else
          */
          {
            // scaling is unreliable, ignore it
            computed[l] = 0;
            m_by_label[l] = NULL;
          }
        }
        label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];

        if (logfp) {
          fprintf(logfp,
                  "%s (%d): linear fit = %2.2f x + "
                  "%2.1f (%d voxels, peak = %2.0f), gca=%2.1f\n",
                  cma_label_to_name(l),
                  l,
                  label_scales[l],
                  label_offsets[l],
                  num,
                  val,
                  label_peaks[l]);
          fflush(logfp);
        }
        fprintf(stdout,
                "%s (%d): linear fit = %2.2f x + "
                "%2.1f (%d voxels, peak = %2.0f), gca=%2.1f\n",
                cma_label_to_name(l),
                l,
                label_scales[l],
                label_offsets[l],
                num,
                val,
                label_peaks[l]);
        fflush(stdout);
        {
          HISTOlinearScale(h_gca, h_gca, label_scales[l], label_offsets[l]);
          if (Gdiag & DIAG_WRITE) {
            sprintf(fname, "%s_label%d_gca_scaled.plt", base_name, l);
            HISTOplot(h_gca, fname);
          }
        }
      }
      else {
        printf("overlap = %g, overlap_threshold = %g\n", overlap, overlap_threshold);
        printf("insufficient overlap %2.4f in histograms - rejecting\n", overlap);
      }

      if (l == Gdiag_no) {
        DiagBreak();
      }
      if (l > 100) {
        break;
      }
    }

    if (getenv("FS_USE_HISTO_COMBOS") != NULL) {
      unsigned int j;
      HISTOmakePDF(h_caudate, h_caudate);
      mri_peak = HISTOfindHighestPeakInRegion(h_caudate, 0, h_caudate->nbins);
      for (j = 0; j <= 1; j++) {
        if (j == 0) {
          l = Left_Caudate;
        }
        else {
          l = Right_Caudate;
        }
        h_gca = gcaGetLabelHistogram(gca, l, 0, 1);
        gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
        HISTOmakePDF(h_gca, h_gca);
        label_peaks[l] = h_gca->bins[gca_peak];
        HISTOfindLinearFit(h_gca, h_caudate, .025, 4, 0, 0, &label_scales[l], &label_offsets[l]);
        val = h_gca->bins[gca_peak] * label_scales[l] + label_offsets[l];
        lower_thresh = 50;
        upper_thresh = 100;
        if ((val < lower_thresh || val > upper_thresh) ||
            (h_caudate->bins[mri_peak] < lower_thresh || h_caudate->bins[mri_peak] > upper_thresh)) {
          printf(
              "%s: unreasonable value (%2.1f/%2.1f), "
              "not in range [%2.0f, %2.0f] - rejecting\n",
              cma_label_to_name(l),
              val,
              h_caudate->bins[mri_peak],
              lower_thresh,
              upper_thresh);
          label_scales[l] = 1.0;
          label_offsets[l] = 1.0;
        }
        else {
          label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];

          if (logfp) {
            fprintf(logfp,
                    "%s (%d): linear fit = %2.2f x + "
                    "%2.1f (peak = %2.0f), gca=%2.1f\n",
                    cma_label_to_name(l),
                    l,
                    label_scales[l],
                    label_offsets[l],
                    val,
                    label_peaks[l]);
            fflush(logfp);
          }
          fprintf(stdout,
                  "%s (%d): linear fit = %2.2f x + "
                  "%2.1f (peak = %2.0f), gca=%2.1f\n",
                  cma_label_to_name(l),
                  l,
                  label_scales[l],
                  label_offsets[l],
                  val,
                  label_peaks[l]);
          fflush(stdout);
          computed[l] = 1;
        }
      }

      HISTOmakePDF(h_mtl, h_mtl);
      mri_peak = HISTOfindHighestPeakInRegion(h_mtl, 0, h_mtl->nbins);
      for (j = 0; j < MTL_LABELS; j++) {
        l = mtl_labels[j];
        h_gca = gcaGetLabelHistogram(gca, l, 0, 1);
        gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
        HISTOmakePDF(h_gca, h_gca);
        label_peaks[l] = h_gca->bins[gca_peak];
        HISTOfindLinearFit(h_gca, h_mtl, .025, 4, 0, 0, &label_scales[l], &label_offsets[l]);
        val = h_gca->bins[gca_peak] * label_scales[l] + label_offsets[l];
        lower_thresh = 50;
        upper_thresh = 90;
        if ((val < lower_thresh || val > upper_thresh) ||
            (h_mtl->bins[mri_peak] < lower_thresh || h_mtl->bins[mri_peak] > upper_thresh)) {
          printf(
              "%s: unreasonable value (%2.1f/%2.1f), "
              "not in range [%2.0f, %2.0f] - rejecting\n",
              cma_label_to_name(l),
              val,
              h_mtl->bins[mri_peak],
              lower_thresh,
              upper_thresh);
          label_scales[l] = 1.0;
          label_offsets[l] = 1.0;
        }
        else {
          label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];

          if (logfp) {
            fprintf(logfp,
                    "%s (%d): linear fit = %2.2f x + "
                    "%2.1f (peak = %2.0f), gca=%2.1f\n",
                    cma_label_to_name(l),
                    l,
                    label_scales[l],
                    label_offsets[l],
                    val,
                    label_peaks[l]);
            fflush(logfp);
          }
          fprintf(stdout,
                  "%s (%d): linear fit = %2.2f x + "
                  "%2.1f (peak = %2.0f), gca=%2.1f\n",
                  cma_label_to_name(l),
                  l,
                  label_scales[l],
                  label_offsets[l],
                  val,
                  label_peaks[l]);
          fflush(stdout);
          computed[l] = 1;
        }
      }
    }

    HISTOfree(&h_gca);
    HISTOfree(&h_mri);
    HISTOfree(&h_mtl);
    HISTOfree(&h_caudate);

    // make sure non-computed labels don't scale
    for (l = 0; l < MAX_CMA_LABELS; l++) {
      if (computed[l] == 0) {
        label_scales[l] = 1.0;
        label_offsets[l] = 0.0;
        h_gca = gcaGetLabelHistogram(gca, l, 0, 1);
        gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
        HISTOmakePDF(h_gca, h_gca);

        if (gca_peak < 0) {
          // fprintf(stderr,
          //      "INFO: GCAmapRenormalizeWithAlignment: "
          //      "gca peak(=%d) < 0\n",gca_peak);
          // fflush(stderr);
          continue;
        }
        if (gca_peak >= h_gca->nbins) {
          fprintf(stderr,
                  "ERROR: GCAmapRenormalizeWithAlignment: "
                  "gca peak(=%d) >= h_gca->nbins(=%d)\n",
                  gca_peak,
                  h_gca->nbins);
          fflush(stderr);
          exit(1);
        }
        if (gca_peak >= 0) {
          printf("gca peak %s = %2.5f (%2.0f)\n", cma_label_to_name(l), h_gca->counts[gca_peak], h_gca->bins[gca_peak]);
        }

        label_peaks[l] = h_gca->bins[gca_peak];
        fflush(stdout);
      }
    }

    if (DIAG_VERBOSE_ON) {
      FILE *fp;
      float scale, offset;
      fp = fopen("norm_offset.plt", "r");
      if (fp != NULL) {
        for (l = 0; l < MAX_CMA_LABELS; l++) {
          if (fscanf(fp, "%d %f %f", &l, &scale, &offset) != 3) {
	      fprintf(stderr,"%s:%d Did not fscanf 3 items\n", __FILE__, __LINE__);
	      exit(1);
	  }
          label_scales[l] = scale;
          label_offsets[l] = offset;
          computed[l] = 1;
        }
        fclose(fp);
      }
    }
    fprintf(stdout, "not using caudate to estimate GM means\n");
    for (k = 0; k < NHEMI_LABELS; k++) {
      int lhl, rhl;
      if (computed[lh_labels[k]] && !computed[rh_labels[k]]) {
        lhl = lh_labels[k];
        rhl = rh_labels[k];
        label_scales[rhl] = label_scales[lhl];
        label_offsets[rhl] = label_offsets[lhl];
        label_peaks[rhl] = label_peaks[lhl];
        computed[rhl] = 1;
        fprintf(stdout,
                "setting label %s based on %s = %2.2f x + %2.0f: %2.0f\n",
                cma_label_to_name(rhl),
                cma_label_to_name(lhl),
                label_scales[rhl],
                label_offsets[rhl],
                label_peaks[rhl]);
      }
      else if (computed[rh_labels[k]] && !computed[lh_labels[k]]) {
        lhl = lh_labels[k];
        rhl = rh_labels[k];
        label_scales[lhl] = label_scales[rhl];
        label_offsets[lhl] = label_offsets[rhl];
        label_peaks[lhl] = label_peaks[rhl];
        computed[lhl] = 1;
        fprintf(stdout,
                "setting label %s based on %s = %2.2f x + %2.0f: %2.0f\n",
                cma_label_to_name(lhl),
                cma_label_to_name(rhl),
                label_scales[lhl],
                label_offsets[lhl],
                label_peaks[lhl]);
      }
    }

    num = 0;
    mean_gm_scale = 0;
    mean_gm_offset = 0;
    for (k = 0; k < NGM_LABELS; k++) {
      label = gm_labels[k];
      if (computed[label]) {
        mean_gm_scale += label_scales[label];
        mean_gm_offset += label_offsets[label];
        num++;
      }
    }
    if (num == 0) {
      mean_gm_scale = 1;
      mean_gm_offset = 0;
    }
    else {
      mean_gm_scale /= (float)num;
      mean_gm_offset /= (float)num;
    }

    num = 0;
    mean_wm_scale = 0;
    mean_wm_offset = 0;
    for (k = 0; k < NWM_LABELS; k++) {
      label = wm_labels[k];
      if (computed[label]) {
        mean_wm_scale += label_scales[label];
        mean_wm_offset += label_offsets[label];
        num++;
      }
    }
    if (num == 0) {
      mean_wm_scale = 1;
      mean_wm_offset = 0;
    }
    else {
      mean_wm_scale /= (float)num;
      mean_wm_offset /= (float)num;
    }

    num = 0;
    mean_csf_scale = 0;
    mean_csf_offset = 0;
    for (k = 0; k < NCSF_LABELS; k++) {
      label = csf_labels[k];
      if (computed[label]) {
        mean_csf_scale += label_scales[label];
        mean_csf_offset += label_offsets[label];
        num++;
      }
    }
    if (num == 0) {
      mean_csf_scale = 1;
      mean_csf_offset = 0;
    }
    else {
      mean_csf_scale /= (float)num;
      mean_csf_offset /= (float)num;
    }

    printf("estimating mean gm scale to be %2.2f x + %2.1f\n", mean_gm_scale, mean_gm_offset);
    printf("estimating mean wm scale to be %2.2f x + %2.1f\n", mean_wm_scale, mean_wm_offset);
    printf("estimating mean csf scale to be %2.2f x + %2.1f\n", mean_csf_scale, mean_csf_offset);

    // assume that cortical gm goes as wm
    if (computed[Left_Cerebral_Cortex] == 0 && computed[Left_Cerebral_White_Matter] != 0) {
      if (m_by_label[Left_Cerebral_White_Matter])
        m_by_label[Left_Cerebral_Cortex] = MatrixCopy(m_by_label[Left_Cerebral_White_Matter], NULL);
      label_scales[Left_Cerebral_Cortex] = mean_gm_scale;
      label_offsets[Left_Cerebral_Cortex] = mean_gm_offset;
      computed[Left_Cerebral_Cortex] = 1;
      l = Left_Cerebral_Cortex;
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    }
    if (computed[Left_Cerebellum_Cortex] == 0) {
      label_scales[Left_Cerebellum_Cortex] = mean_gm_scale;
      label_offsets[Left_Cerebellum_Cortex] = mean_gm_offset;
      computed[Left_Cerebellum_Cortex] = 1;
      l = Left_Cerebellum_Cortex;
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
      printf("setting left cbm cortex = %2.2f x + %2.2f\n", mean_gm_scale, mean_gm_offset);
    }
    if (computed[Right_Cerebellum_Cortex] == 0) {
      label_scales[Right_Cerebellum_Cortex] = mean_gm_scale;
      label_offsets[Right_Cerebellum_Cortex] = mean_gm_offset;
      computed[Right_Cerebellum_Cortex] = 1;
      l = Right_Cerebellum_Cortex;
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
      printf("setting right cbm cortex = %2.2f x + %2.2f\n", mean_gm_scale, mean_gm_offset);
    }
    if (computed[Right_Cerebral_Cortex] == 0 && computed[Right_Cerebral_White_Matter] != 0) {
      if (m_by_label[Right_Cerebral_White_Matter])
        m_by_label[Right_Cerebral_Cortex] = MatrixCopy(m_by_label[Right_Cerebral_White_Matter], NULL);
      label_scales[Right_Cerebral_Cortex] = mean_gm_scale;
      label_offsets[Right_Cerebral_Cortex] = mean_gm_offset;
      l = Right_Cerebral_Cortex;
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
      computed[Right_Cerebral_Cortex] = 1;
    }

    // lock some labels scaling to others that have been estimated
    if (computed[Left_Caudate]) {
      label_offsets[Left_Accumbens_area] = label_offsets[Left_Caudate];
      label_scales[Left_Accumbens_area] = label_scales[Left_Caudate];
      l = Left_Accumbens_area;
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
      computed[Left_Accumbens_area] = 1;
    }
    if (computed[Right_Caudate]) {
      label_offsets[Right_Accumbens_area] = label_offsets[Right_Caudate];
      label_scales[Right_Accumbens_area] = label_scales[Right_Caudate];
      l = Right_Accumbens_area;
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
      computed[Right_Accumbens_area] = 1;
    }
    if (computed[Left_Inf_Lat_Vent] == 0 && computed[Left_Hippocampus] != 0) {
      label_scales[Left_Inf_Lat_Vent] = label_scales[Left_Hippocampus];
      label_offsets[Left_Inf_Lat_Vent] = label_offsets[Left_Hippocampus];
      l = Left_Inf_Lat_Vent;
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
      computed[Left_Inf_Lat_Vent] = 1;
    }
    if (computed[Right_Inf_Lat_Vent] == 0 && computed[Right_Hippocampus] != 0) {
      label_scales[Right_Inf_Lat_Vent] = label_scales[Right_Hippocampus];
      label_offsets[Right_Inf_Lat_Vent] = label_offsets[Right_Hippocampus];
      l = Right_Inf_Lat_Vent;
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
      computed[Right_Inf_Lat_Vent] = 1;
    }

    label_scales[CSF] = mean_csf_scale;
    label_scales[Fifth_Ventricle] = mean_csf_scale;
    label_offsets[CSF] = mean_csf_offset;
    label_offsets[Fifth_Ventricle] = mean_csf_offset;
    computed[CSF] = computed[Fifth_Ventricle] = 1;
    l = CSF;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    l = Fifth_Ventricle;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];

    // set the scale and offset for the rest; added by xhan
    for (l = 0; l < MAX_CMA_LABELS; l++) {
      if (l == Gdiag_no) {
        DiagBreak();
      }
      if (computed[l] > 0) {
        continue;
      }
      if (equiv_class[l] == 1) {
        label_scales[l] = mean_csf_scale;
        label_offsets[l] = mean_csf_offset;
        computed[l] = 1;
      }
      else if (equiv_class[l] == 2) {
        label_scales[l] = mean_wm_scale;
        label_offsets[l] = mean_wm_offset;
        computed[l] = 1;
      }
      else if (equiv_class[l] == 3) {
        label_scales[l] = mean_gm_scale;
        label_offsets[l] = mean_gm_offset;
        computed[l] = 1;
      }
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    }

    // make sure labels are self-consistent
    for (l = 0; l < MAX_CMA_LABELS; l++) {
      double peak, scale;
      switch (l) {
        case Left_Pallidum:
        case Right_Pallidum:
          if ((label_peaks[l] >= .97 * label_peaks[Left_Cerebral_White_Matter]) ||
              (label_peaks[l] >= .97 * label_peaks[Right_Cerebral_White_Matter])) {
            // don't let pallidum be as bright as wm
            peak = 0.97 * (label_peaks[Left_Cerebral_White_Matter] + label_peaks[Right_Cerebral_White_Matter]) / 2;
            scale = peak / label_peaks[l];
            printf(
                "%s too bright - rescaling by %2.3f "
                "(from %2.3f) to %2.1f (was %2.1f)\n",
                cma_label_to_name(l),
                scale,
                label_scales[l],
                peak,
                label_peaks[l]);
            label_scales[l] *= scale;
            label_peaks[l] = peak;
          }
          break;
        case Left_Putamen:
        case Right_Putamen:
          if ((label_peaks[l] >= .9 * label_peaks[Left_Cerebral_White_Matter]) ||
              (label_peaks[l] >= .9 * label_peaks[Right_Cerebral_White_Matter])) {
            // don't let putamen be as bright as wm
            peak = 0.9 * (label_peaks[Left_Cerebral_White_Matter] + label_peaks[Right_Cerebral_White_Matter]) / 2;
            scale = peak / label_peaks[l];
            printf(
                "%s too bright - rescaling by %2.3f "
                "(from %2.3f) to %2.1f (was %2.1f)\n",
                cma_label_to_name(l),
                scale,
                label_scales[l],
                peak,
                label_peaks[l]);
            label_scales[l] *= scale;
            label_peaks[l] = peak;
          }
          break;
      }
    }

// now use precomputed eigenstructure of intensity
// volumes to adjust corrections

    if (logfp) {
      for (l = 0; l < MAX_CMA_LABELS; l++)
        if (computed[l] != 0)
          fprintf(logfp,
                  "label %s: scaling by %2.2f  + %2.1f to %2.0f\n",
                  cma_label_to_name(l),
                  label_scales[l],
                  label_offsets[l],
                  label_peaks[l]);
      fflush(logfp);
    }
    if (DIAG_VERBOSE_ON) {
      FILE *fp;
      fp = fopen("norm_offset.plt", "w");
      for (l = 0; l < MAX_CMA_LABELS; l++)
        if (computed[l] != 0) {
          fprintf(fp, "%d %f %f\n", l, label_scales[l], label_offsets[l]);
        }
      fclose(fp);
    }

    if (base_name) {
      FILE *fp;
      char fname[STRLEN];
      sprintf(fname, "%s.label_intensities.txt", base_name);
      printf("saving intensity scales to %s\n", fname);
      fp = fopen(fname, "w");
      if (fp == NULL) {
        ErrorExit(ERROR_NOFILE, "%s: could not open intensity tracking file %s", Progname, fname);
      }

      for (l = 0; l < MAX_CMA_LABELS; l++)
        if (computed[l] != 0)
          fprintf(fp,
                  "%d %s %2.2f %2.1f %2.0f\n",
                  l,
                  cma_label_to_name(l),
                  label_scales[l],
                  label_offsets[l],
                  label_peaks[l]);

      fflush(fp);
      fclose(fp);
      if (getenv("EXIT_AFTER_INT") != NULL) {
        exit(0);
      }
    }
    gcaCheck(gca);
    GCAapplyRenormalization(gca, label_scales, label_offsets, frame);
    gcaCheck(gca);
  }
  if (plta)  // return linear transform array to caller
  {
    int i;

    // count # of xforms
    for (i = l = 0; l < MAX_CMA_LABELS; l++) {
      if (m_by_label[l] != NULL && computed[l] == 2) {
        i++;
      }
    }

    if (i > 0)  // should always be true
    {
      *plta = lta = LTAalloc(i, mri);
      for (i = l = 0; l < MAX_CMA_LABELS; l++) {
        if (m_by_label[l] != NULL && computed[l] == 2) {
          MatrixCopy(m_by_label[l], lta->xforms[i].m_L);
          MatrixFree(&m_by_label[l]);
          lta->xforms[i].label = l;
          i++;
        }
      }
    }
    printf("%d transforms computed\n", i);
  }

  if (label_scales != plabel_scales)  // copy to user-supplied space
  {
    memmove(plabel_scales, label_scales, MAX_CMA_LABELS * sizeof(label_scales[0]));
  }
  if (label_offsets != plabel_offsets)  // copy to user-supplied space
  {
    memmove(plabel_offsets, label_offsets, MAX_CMA_LABELS * sizeof(label_offsets[0]));
  }
  if (label_peaks != plabel_peaks)  // copy to user-supplied space
  {
    memmove(plabel_peaks, label_peaks, MAX_CMA_LABELS * sizeof(label_peaks[0]));
  }
  if (computed != plabel_computed)  // copy to user-supplied space
  {
    memmove(plabel_computed, computed, MAX_CMA_LABELS * sizeof(computed[0]));
  }

  if (mri_seg) {
    MRIfree(&mri_seg);
  }
  return (NO_ERROR);
}
/*
  !!!!!!!!!!!!!! Note the difference in the longitudinal version is that every frame is
  the same modality but a different time point of the same subject.
*/
int GCAcomputeRenormalizationWithAlignmentLongitudinal(GCA *gca,
                                                       MRI *mri,
                                                       TRANSFORM *transform,
                                                       FILE *logfp,
                                                       const char *base_name,
                                                       LTA **plta,
                                                       int handle_expanded_ventricles,
                                                       float *plabel_scales,
                                                       float *plabel_offsets,
                                                       float *plabel_peaks,
                                                       int *plabel_computed)
{
  HISTOGRAM *h_mri, *h_gca, *h_mtl = NULL, *h_caudate;
  unsigned int k, j;
  int l, nbins, i, x, y, z, num, frame, bin, computed[MAX_CMA_LABELS], label, gca_peak, mri_peak;
  float fmin, fmax, label_scales[MAX_CMA_LABELS], overlap, mean_gm_scale, mean_wm_scale, mean_csf_scale,
      label_peaks[MAX_CMA_LABELS], label_offsets[MAX_CMA_LABELS], mean_wm_offset, mean_csf_offset, mean_gm_offset,
      lower_thresh, upper_thresh;
  double val /*, scale*/;
  MRI *mri_seg = NULL, *mri_aligned, *mri_labels = NULL, *mri_borders;
  char fname[STRLEN];
  MATRIX *m_L, *m_by_label[MAX_CMA_LABELS];
  LTA *lta;
  double det = -1;
  float peak_threshold = 0.03;
  float overlap_threshold = 0.001;
  int equiv_class[MAX_CMA_LABELS];

  if (plabel_scales == NULL) plabel_scales = label_scales;

  if (plabel_offsets == NULL) plabel_offsets = label_offsets;

  if (plabel_peaks == NULL) plabel_peaks = label_peaks;

  if (plabel_computed == NULL) plabel_computed = computed;

  memset(label_peaks, 0, sizeof(label_peaks));
  if (transform->type == MORPH_3D_TYPE) {
    peak_threshold = 0.01;
    overlap_threshold = -1.0;  // at mri_ca_label stage;
    // trust the registration more
  }

  set_equilavent_classes(equiv_class);

  printf("renormalizing by structure alignment....\n");
  if (plta)
    lta = *plta;
  else
    lta = NULL;

  mri_borders = MRImarkPossibleBorders(mri, NULL, 15);
  for (l = 0; l < MAX_CMA_LABELS; l++) {
    if (l == Gdiag_no) DiagBreak();

    label_scales[l] = 1.0;
    label_offsets[l] = 0.0;
    computed[l] = 0;
    m_by_label[l] = NULL;  // not estimated yet
  }

  MRIvalRange(mri, &fmin, &fmax);
  nbins = 256;
  h_mri = HISTOalloc(nbins);
  h_mtl = HISTOalloc(nbins);
  h_caudate = HISTOalloc(nbins);
  for (j = 0; j < NALIGN_LABELS; j++) {
    l = align_labels[j];
    if (l == Gdiag_no) DiagBreak();

    mri_seg = MRIalloc(mri->width, mri->height, mri->depth, MRI_UCHAR);
    MRIcopyHeader(mri, mri_seg);
    mri_labels = MRIclone(mri_seg, mri_labels);

    /* include 2 voxel border to get context around structure.
       e.g. hippo is made easier to find by wm inferior and ventricle
       posterior.
    */
    if (transform->type != MORPH_3D_TYPE) ErrorExit(ERROR_UNSUPPORTED, "Long renorm must use 3d morph\n");

    m_L = NULL;
    if (l == Left_Cerebral_White_Matter || l == Right_Cerebral_White_Matter) {
      // wm so big it's hard to localize with a linear xform
      GCAbuildMostLikelyVolumeForStructure(gca, mri_seg, l, 0, transform, NULL);
      MRIerode(mri_seg, mri_seg);
    }
    else {
      GCAbuildMostLikelyVolumeForStructure(gca, mri_seg, l, 0, transform, NULL);
    }
    mri_aligned = MRIerode(mri_seg, NULL);

    MRIbinarize(mri_aligned, mri_aligned, 1, 0, 128);
    if (Gdiag & DIAG_WRITE) {
      sprintf(fname, "%s_label%d_eroded.mgz", base_name, l);
      MRIwrite(mri_aligned, fname);
    }
    if (l == Gdiag_no) DiagBreak();

    HISTOclear(h_mri, h_mri);
    h_mri->bin_size = (fmax - fmin) / 255.0;
    if (h_mri->bin_size < 1 && (mri->type == MRI_UCHAR || mri->type == MRI_SHORT)) h_mri->bin_size = 1;

    for (i = 0; i < nbins; i++) h_mri->bins[i] = (i + 1) * h_mri->bin_size;

    for (num = x = 0; x < mri_aligned->width; x++) {
      for (y = 0; y < mri_aligned->height; y++) {
        for (z = 0; z < mri_aligned->depth; z++) {
          if (x == Gx && y == Gy && z == Gz) DiagBreak();

          MRIsampleVolume(mri_aligned, x, y, z, &val);
          if (DZERO(val))  // not in this structure
            continue;

          for (frame = 0; frame < mri->nframes; frame++) {
            MRIsampleVolumeFrame(mri, x, y, z, frame, &val);

            if (FZERO(val))  // skull stripped
              continue;

            if (MRIgetVoxVal(mri_borders, x, y, z, 0) > 0)  // avoid border voxels - matches erosion
              continue;

            bin = nint((val - fmin) / h_mri->bin_size);
            if (bin >= h_mri->nbins)
              bin = h_mri->nbins - 1;
            else if (bin < 0)
              bin = 0;

            h_mri->counts[bin]++;
            num++;
          }
        }
      }
    }
    if (l == Gdiag_no) DiagBreak();

    MRIfree(&mri_aligned);

    h_gca = gcaGetLabelHistogram(gca, l, 0, 0);
    gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
    HISTOmakePDF(h_gca, h_gca);

    if (gca_peak < 0) {
      // fprintf(stderr,
      //      "INFO: GCAmapRenormalizeWithAlignment: "
      //      "gca peak(=%d) < 0\n",gca_peak);
      // fflush(stderr);
      continue;
    }
    if (gca_peak >= h_gca->nbins) {
      fprintf(stderr,
              "ERROR: GCAmapRenormalizeWithAlignment: "
              "gca peak(=%d) >= h_gca->nbins(=%d)\n",
              gca_peak,
              h_gca->nbins);
      fflush(stderr);
      exit(1);
    }
    if (gca_peak >= 0) printf("gca peak = %2.5f (%2.0f)\n", h_gca->counts[gca_peak], h_gca->bins[gca_peak]);

    label_peaks[l] = h_gca->bins[gca_peak];
    fflush(stdout);

    switch (l) {
      case Brain_Stem:
      case Left_VentralDC:
      case Right_VentralDC:
        lower_thresh = 80;
        upper_thresh = 110;
        break;
      case Left_Caudate:
      case Right_Caudate:
        lower_thresh = 50;
        upper_thresh = 100;
        break;
      case Left_Cerebral_Cortex:
      case Right_Cerebral_Cortex:
        lower_thresh = 40;
        upper_thresh = 95;
        break;
      case Left_Pallidum:
      case Right_Pallidum:
        lower_thresh = 75;
        upper_thresh = 135;
        break;
      case Left_Thalamus:
      case Right_Thalamus:
        lower_thresh = 75;
        upper_thresh = 120;
        break;
      case Left_Hippocampus:
      case Right_Hippocampus:
      case Left_Amygdala:
      case Right_Amygdala:
        lower_thresh = 50;
        upper_thresh = 90;
        break;
      case Left_Cerebral_White_Matter:
      case Right_Cerebral_White_Matter:
        lower_thresh = 90;
        upper_thresh = 130;
        break;
      case Left_Putamen:
      case Right_Putamen:
        lower_thresh = 60;
        upper_thresh = 107;
        break;
      case Left_Lateral_Ventricle:
      case Right_Lateral_Ventricle:
      case Third_Ventricle:
      case Fourth_Ventricle:
      case CSF:
        lower_thresh = 0;
        upper_thresh = 55;
        break;
      case Left_Inf_Lat_Vent:
      case Right_Inf_Lat_Vent:
        lower_thresh = 0;
        upper_thresh = 65;
        break;
      default:
        lower_thresh = 0;
        upper_thresh = 256;
        break;
    }
    HISTOclearBins(h_mri, h_mri, 0, lower_thresh - 1);
    HISTOclearBins(h_mri, h_mri, upper_thresh + 1, 255);
    mri_peak = HISTOfindHighestPeakInRegion(h_mri, 0, h_mri->nbins);
    HISTOfillHoles(h_mri);
    if (l == Left_Caudate || l == Right_Caudate) {
      h_caudate->bin_size = 1;
      h_caudate->min = 0;
      h_caudate->max = 255;
      for (i = 0; i < nbins; i++) h_caudate->bins[i] = (i + 1) * h_caudate->bin_size;

      HISTOadd(h_mri, h_caudate, h_caudate);
    }

    if (l == Left_Hippocampus || l == Right_Hippocampus || l == Left_Amygdala || l == Right_Amygdala) {
      h_mtl->bin_size = 1;
      h_mtl->min = 0;
      h_mtl->max = 255;
      for (i = 0; i < nbins; i++) h_mtl->bins[i] = (i + 1) * h_mtl->bin_size;

      HISTOadd(h_mri, h_mtl, h_mtl);
    }
    HISTOmakePDF(h_mri, h_mri);
    if (mri_peak >= 0) printf("mri peak = %2.5f (%2.0f)\n", h_mri->counts[mri_peak], h_mri->bins[mri_peak]);
    fflush(stdout);

    if (IS_CSF(l) && h_mri->bins[mri_peak] > 55) {
      printf("CSF peak too bright - rejecting\n");
      continue;
    }
    if (h_mri->counts[mri_peak] < peak_threshold || num <= 50)
    /* not enough to reliably estimate density */
    {
      if (h_mri->counts[mri_peak] < peak_threshold) printf("uniform distribution in MR - rejecting arbitrary fit\n");
      if (m_L) MatrixFree(&m_L);

      continue;
    }
    if (m_L) {
      if (plta && (!IS_GM(l)))  // GM will be copied from WM later
        m_by_label[l] = m_L;    // store if for assembling an LTA later
      else
        MatrixFree(&m_L);
    }

    if (Gdiag & DIAG_WRITE) {
      sprintf(fname, "%s_label%d_mri.plt", base_name, l);
      HISTOplot(h_mri, fname);
      sprintf(fname, "%s_label%d_gca.plt", base_name, l);
      HISTOplot(h_gca, fname);
      DiagBreak();
    }
    overlap = HISTOthreshSum(h_mri, h_gca, .025);

    // if (overlap > 0.01)
    //      if (overlap > 0.001)
    if (IS_LAT_VENT(l) || overlap > overlap_threshold) {
      //                        if (l == Gdiag_no)
      //  HISTOfindLinearFit(h_gca, h_mri, .025, 10, -75, 75,
      // &label_scales[l],  &label_offsets[l]) ;
      //              HISTOfindLinearFit(h_gca, h_mri, .025,
      // 4, -125, 125, &label_scales[l], &label_offsets[l]) ;
      HISTOfindLinearFit(h_gca, h_mri, .025, 4, 0, 0, &label_scales[l], &label_offsets[l]);

      val = h_gca->bins[gca_peak] * label_scales[l] + label_offsets[l];
      if ((val < lower_thresh || val > upper_thresh) ||
          (h_mri->bins[mri_peak] < lower_thresh || h_mri->bins[mri_peak] > upper_thresh)) {
        //          if (transform->type != MORPH_3D_TYPE)
        {
          printf(
              "%s: unreasonable value (%2.1f/%2.1f), "
              "not in range [%2.0f, %2.0f] - rejecting\n",
              cma_label_to_name(l),
              val,
              h_mri->bins[mri_peak],
              lower_thresh,
              upper_thresh);
          label_scales[l] = 1.0;
          label_offsets[l] = 1.0;
          continue;
        }
      }

      // only allow certain labels to be used for initializing the 3d morph
      // (which is what happens when computed[l] = 2)
      computed[l] = (det > 0) ? 2 : 1;
      printf(
          "%s (%d): linear fit = %2.2f x + %2.1f "
          "(%d voxels, overlap=%2.3f)\n",
          cma_label_to_name(l),
          l,
          label_scales[l],
          label_offsets[l],
          num,
          overlap);

      // note that the following range need be changed
      // if both scale and offset are allowed' 1/1.5 = 0.67
      if (IS_LAT_VENT(l)) {
        if (label_scales[l] < 0.4)
          label_scales[l] = 0.4;
        else if (label_scales[l] > 1.5)
          label_scales[l] = 1.5;
      }
      if ((label_scales[l] < 0.67 || (label_scales[l] > 1.5)) && !IS_LAT_VENT(l)) {
        /*
          if(IS_CSF(l)){
          if(label_scales[l] < 0.67) label_scales[l] = 0.67;
          else if(label_scales[l] > 1.5) label_scales[l] = 1.5;
          } else
        */
        {
          // scaling is unreliable, ignore it
          computed[l] = 0;
          m_by_label[l] = NULL;
        }
      }
      label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];

      if (logfp) {
        fprintf(logfp,
                "%s (%d): linear fit = %2.2f x + "
                "%2.1f (%d voxels, peak = %2.0f), gca=%2.1f\n",
                cma_label_to_name(l),
                l,
                label_scales[l],
                label_offsets[l],
                num,
                val,
                label_peaks[l]);
        fflush(logfp);
      }
      fprintf(stdout,
              "%s (%d): linear fit = %2.2f x + "
              "%2.1f (%d voxels, peak = %2.0f), gca=%2.1f\n",
              cma_label_to_name(l),
              l,
              label_scales[l],
              label_offsets[l],
              num,
              val,
              label_peaks[l]);
      fflush(stdout);
      {
        HISTOlinearScale(h_gca, h_gca, label_scales[l], label_offsets[l]);
        if (Gdiag & DIAG_WRITE) {
          sprintf(fname, "%s_label%d_gca_scaled.plt", base_name, l);
          HISTOplot(h_gca, fname);
        }
      }
    }
    else {
      printf("overlap = %g, overlap_threshold = %g\n", overlap, overlap_threshold);
      printf("insufficient overlap %2.4f in histograms - rejecting\n", overlap);
    }

    if (l == Gdiag_no) DiagBreak();

    if (l > 100) break;
  }

  if (getenv("FS_USE_HISTO_COMBOS") != NULL) {
    unsigned int j;
    HISTOmakePDF(h_caudate, h_caudate);
    mri_peak = HISTOfindHighestPeakInRegion(h_caudate, 0, h_caudate->nbins);
    for (j = 0; j <= 1; j++) {
      if (j == 0)
        l = Left_Caudate;
      else
        l = Right_Caudate;

      h_gca = gcaGetLabelHistogram(gca, l, 0, 1);
      gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
      HISTOmakePDF(h_gca, h_gca);
      label_peaks[l] = h_gca->bins[gca_peak];
      HISTOfindLinearFit(h_gca, h_caudate, .025, 4, 0, 0, &label_scales[l], &label_offsets[l]);
      val = h_gca->bins[gca_peak] * label_scales[l] + label_offsets[l];
      lower_thresh = 50;
      upper_thresh = 100;
      if ((val < lower_thresh || val > upper_thresh) ||
          (h_caudate->bins[mri_peak] < lower_thresh || h_caudate->bins[mri_peak] > upper_thresh)) {
        printf(
            "%s: unreasonable value (%2.1f/%2.1f), "
            "not in range [%2.0f, %2.0f] - rejecting\n",
            cma_label_to_name(l),
            val,
            h_caudate->bins[mri_peak],
            lower_thresh,
            upper_thresh);
        label_scales[l] = 1.0;
        label_offsets[l] = 1.0;
      }
      else {
        label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];

        if (logfp) {
          fprintf(logfp,
                  "%s (%d): linear fit = %2.2f x + "
                  "%2.1f (peak = %2.0f), gca=%2.1f\n",
                  cma_label_to_name(l),
                  l,
                  label_scales[l],
                  label_offsets[l],
                  val,
                  label_peaks[l]);
          fflush(logfp);
        }
        fprintf(stdout,
                "%s (%d): linear fit = %2.2f x + "
                "%2.1f (peak = %2.0f), gca=%2.1f\n",
                cma_label_to_name(l),
                l,
                label_scales[l],
                label_offsets[l],
                val,
                label_peaks[l]);
        fflush(stdout);
        computed[l] = 1;
      }
    }

    HISTOmakePDF(h_mtl, h_mtl);
    mri_peak = HISTOfindHighestPeakInRegion(h_mtl, 0, h_mtl->nbins);
    for (j = 0; j < MTL_LABELS; j++) {
      l = mtl_labels[j];
      h_gca = gcaGetLabelHistogram(gca, l, 0, 1);
      gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
      HISTOmakePDF(h_gca, h_gca);
      label_peaks[l] = h_gca->bins[gca_peak];
      HISTOfindLinearFit(h_gca, h_mtl, .025, 4, 0, 0, &label_scales[l], &label_offsets[l]);
      val = h_gca->bins[gca_peak] * label_scales[l] + label_offsets[l];
      lower_thresh = 50;
      upper_thresh = 90;
      if ((val < lower_thresh || val > upper_thresh) ||
          (h_mtl->bins[mri_peak] < lower_thresh || h_mtl->bins[mri_peak] > upper_thresh)) {
        printf(
            "%s: unreasonable value (%2.1f/%2.1f), "
            "not in range [%2.0f, %2.0f] - rejecting\n",
            cma_label_to_name(l),
            val,
            h_mtl->bins[mri_peak],
            lower_thresh,
            upper_thresh);
        label_scales[l] = 1.0;
        label_offsets[l] = 1.0;
      }
      else {
        label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];

        if (logfp) {
          fprintf(logfp,
                  "%s (%d): linear fit = %2.2f x + "
                  "%2.1f (peak = %2.0f), gca=%2.1f\n",
                  cma_label_to_name(l),
                  l,
                  label_scales[l],
                  label_offsets[l],
                  val,
                  label_peaks[l]);
          fflush(logfp);
        }
        fprintf(stdout,
                "%s (%d): linear fit = %2.2f x + "
                "%2.1f (peak = %2.0f), gca=%2.1f\n",
                cma_label_to_name(l),
                l,
                label_scales[l],
                label_offsets[l],
                val,
                label_peaks[l]);
        fflush(stdout);
        computed[l] = 1;
      }
    }
  }

  HISTOfree(&h_gca);
  HISTOfree(&h_mri);
  HISTOfree(&h_mtl);
  HISTOfree(&h_caudate);

  // make sure non-computed labels don't scale
  for (l = 0; l < MAX_CMA_LABELS; l++) {
    if (computed[l] == 0) {
      label_scales[l] = 1.0;
      label_offsets[l] = 0.0;
      h_gca = gcaGetLabelHistogram(gca, l, 0, 1);
      gca_peak = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
      HISTOmakePDF(h_gca, h_gca);

      if (gca_peak < 0) {
        // fprintf(stderr,
        //      "INFO: GCAmapRenormalizeWithAlignment: "
        //      "gca peak(=%d) < 0\n",gca_peak);
        // fflush(stderr);
        continue;
      }
      if (gca_peak >= h_gca->nbins) {
        fprintf(stderr,
                "ERROR: GCAmapRenormalizeWithAlignment: "
                "gca peak(=%d) >= h_gca->nbins(=%d)\n",
                gca_peak,
                h_gca->nbins);
        fflush(stderr);
        exit(1);
      }
      if (gca_peak >= 0) {
        printf("gca peak %s = %2.5f (%2.0f)\n", cma_label_to_name(l), h_gca->counts[gca_peak], h_gca->bins[gca_peak]);
      }

      label_peaks[l] = h_gca->bins[gca_peak];
      fflush(stdout);
    }
  }

  if (DIAG_VERBOSE_ON) {
    FILE *fp;
    float scale, offset;
    fp = fopen("norm_offset.plt", "r");
    if (fp != NULL) {
      for (l = 0; l < MAX_CMA_LABELS; l++) {
        if (fscanf(fp, "%d %f %f", &l, &scale, &offset) != 3) {
           fprintf(stderr,"%s:%d Did not fscanf 3 items\n", __FILE__, __LINE__);
	   exit(1);
	}
        label_scales[l] = scale;
        label_offsets[l] = offset;
        computed[l] = 1;
      }
      fclose(fp);
    }
  }
  fprintf(stdout, "not using caudate to estimate GM means\n");
  for (k = 0; k < NHEMI_LABELS; k++) {
    int lhl, rhl;
    if (computed[lh_labels[k]] && !computed[rh_labels[k]]) {
      lhl = lh_labels[k];
      rhl = rh_labels[k];
      label_scales[rhl] = label_scales[lhl];
      label_offsets[rhl] = label_offsets[lhl];
      label_peaks[rhl] = label_peaks[lhl];
      computed[rhl] = 1;
      fprintf(stdout,
              "setting label %s based on %s = %2.2f x + %2.0f: %2.0f\n",
              cma_label_to_name(rhl),
              cma_label_to_name(lhl),
              label_scales[rhl],
              label_offsets[rhl],
              label_peaks[rhl]);
    }
    else if (computed[rh_labels[k]] && !computed[lh_labels[k]]) {
      lhl = lh_labels[k];
      rhl = rh_labels[k];
      label_scales[lhl] = label_scales[rhl];
      label_offsets[lhl] = label_offsets[rhl];
      label_peaks[lhl] = label_peaks[rhl];
      computed[lhl] = 1;
      fprintf(stdout,
              "setting label %s based on %s = %2.2f x + %2.0f: %2.0f\n",
              cma_label_to_name(lhl),
              cma_label_to_name(rhl),
              label_scales[lhl],
              label_offsets[lhl],
              label_peaks[lhl]);
    }
  }

  num = 0;
  mean_gm_scale = 0;
  mean_gm_offset = 0;
  for (k = 0; k < NGM_LABELS; k++) {
    label = gm_labels[k];
    if (computed[label]) {
      mean_gm_scale += label_scales[label];
      mean_gm_offset += label_offsets[label];
      num++;
    }
  }
  if (num == 0) {
    mean_gm_scale = 1;
    mean_gm_offset = 0;
  }
  else {
    mean_gm_scale /= (float)num;
    mean_gm_offset /= (float)num;
  }

  num = 0;
  mean_wm_scale = 0;
  mean_wm_offset = 0;
  for (k = 0; k < NWM_LABELS; k++) {
    label = wm_labels[k];
    if (computed[label]) {
      mean_wm_scale += label_scales[label];
      mean_wm_offset += label_offsets[label];
      num++;
    }
  }
  if (num == 0) {
    mean_wm_scale = 1;
    mean_wm_offset = 0;
  }
  else {
    mean_wm_scale /= (float)num;
    mean_wm_offset /= (float)num;
  }

  num = 0;
  mean_csf_scale = 0;
  mean_csf_offset = 0;
  for (k = 0; k < NCSF_LABELS; k++) {
    label = csf_labels[k];
    if (computed[label]) {
      mean_csf_scale += label_scales[label];
      mean_csf_offset += label_offsets[label];
      num++;
    }
  }
  if (num == 0) {
    mean_csf_scale = 1;
    mean_csf_offset = 0;
  }
  else {
    mean_csf_scale /= (float)num;
    mean_csf_offset /= (float)num;
  }

  printf("estimating mean gm scale to be %2.2f x + %2.1f\n", mean_gm_scale, mean_gm_offset);
  printf("estimating mean wm scale to be %2.2f x + %2.1f\n", mean_wm_scale, mean_wm_offset);
  printf("estimating mean csf scale to be %2.2f x + %2.1f\n", mean_csf_scale, mean_csf_offset);

  // assume that cortical gm goes as wm
  if (computed[Left_Cerebral_Cortex] == 0 && computed[Left_Cerebral_White_Matter] != 0) {
    if (m_by_label[Left_Cerebral_White_Matter])
      m_by_label[Left_Cerebral_Cortex] = MatrixCopy(m_by_label[Left_Cerebral_White_Matter], NULL);
    label_scales[Left_Cerebral_Cortex] = mean_gm_scale;
    label_offsets[Left_Cerebral_Cortex] = mean_gm_offset;
    computed[Left_Cerebral_Cortex] = 1;
    l = Left_Cerebral_Cortex;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
  }
  if (computed[Left_Cerebellum_Cortex] == 0) {
    label_scales[Left_Cerebellum_Cortex] = mean_gm_scale;
    label_offsets[Left_Cerebellum_Cortex] = mean_gm_offset;
    computed[Left_Cerebellum_Cortex] = 1;
    l = Left_Cerebellum_Cortex;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    printf("setting left cbm cortex = %2.2f x + %2.2f\n", mean_gm_scale, mean_gm_offset);
  }
  if (computed[Right_Cerebellum_Cortex] == 0) {
    label_scales[Right_Cerebellum_Cortex] = mean_gm_scale;
    label_offsets[Right_Cerebellum_Cortex] = mean_gm_offset;
    computed[Right_Cerebellum_Cortex] = 1;
    l = Right_Cerebellum_Cortex;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    printf("setting right cbm cortex = %2.2f x + %2.2f\n", mean_gm_scale, mean_gm_offset);
  }
  if (computed[Right_Cerebral_Cortex] == 0 && computed[Right_Cerebral_White_Matter] != 0) {
    if (m_by_label[Right_Cerebral_White_Matter])
      m_by_label[Right_Cerebral_Cortex] = MatrixCopy(m_by_label[Right_Cerebral_White_Matter], NULL);
    label_scales[Right_Cerebral_Cortex] = mean_gm_scale;
    label_offsets[Right_Cerebral_Cortex] = mean_gm_offset;
    l = Right_Cerebral_Cortex;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    computed[Right_Cerebral_Cortex] = 1;
  }

  // lock some labels scaling to others that have been estimated
  if (computed[Left_Caudate]) {
    label_offsets[Left_Accumbens_area] = label_offsets[Left_Caudate];
    label_scales[Left_Accumbens_area] = label_scales[Left_Caudate];
    l = Left_Accumbens_area;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    computed[Left_Accumbens_area] = 1;
  }
  if (computed[Right_Caudate]) {
    label_offsets[Right_Accumbens_area] = label_offsets[Right_Caudate];
    label_scales[Right_Accumbens_area] = label_scales[Right_Caudate];
    l = Right_Accumbens_area;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    computed[Right_Accumbens_area] = 1;
  }
  if (computed[Left_Inf_Lat_Vent] == 0 && computed[Left_Hippocampus] != 0) {
    label_scales[Left_Inf_Lat_Vent] = label_scales[Left_Hippocampus];
    label_offsets[Left_Inf_Lat_Vent] = label_offsets[Left_Hippocampus];
    l = Left_Inf_Lat_Vent;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    computed[Left_Inf_Lat_Vent] = 1;
  }
  if (computed[Right_Inf_Lat_Vent] == 0 && computed[Right_Hippocampus] != 0) {
    label_scales[Right_Inf_Lat_Vent] = label_scales[Right_Hippocampus];
    label_offsets[Right_Inf_Lat_Vent] = label_offsets[Right_Hippocampus];
    l = Right_Inf_Lat_Vent;
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
    computed[Right_Inf_Lat_Vent] = 1;
  }

  label_scales[CSF] = mean_csf_scale;
  label_scales[Fifth_Ventricle] = mean_csf_scale;
  label_offsets[CSF] = mean_csf_offset;
  label_offsets[Fifth_Ventricle] = mean_csf_offset;
  computed[CSF] = computed[Fifth_Ventricle] = 1;
  l = CSF;
  label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
  l = Fifth_Ventricle;
  label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];

  // set the scale and offset for the rest; added by xhan
  for (l = 0; l < MAX_CMA_LABELS; l++) {
    if (l == Gdiag_no) DiagBreak();

    if (computed[l] > 0) continue;

    if (equiv_class[l] == 1) {
      label_scales[l] = mean_csf_scale;
      label_offsets[l] = mean_csf_offset;
      computed[l] = 1;
    }
    else if (equiv_class[l] == 2) {
      label_scales[l] = mean_wm_scale;
      label_offsets[l] = mean_wm_offset;
      computed[l] = 1;
    }
    else if (equiv_class[l] == 3) {
      label_scales[l] = mean_gm_scale;
      label_offsets[l] = mean_gm_offset;
      computed[l] = 1;
    }
    label_peaks[l] = label_peaks[l] * label_scales[l] + label_offsets[l];
  }

  // make sure labels are self-consistent
  for (l = 0; l < MAX_CMA_LABELS; l++) {
    double peak, scale;
    switch (l) {
      case Left_Pallidum:
      case Right_Pallidum:
        if ((label_peaks[l] >= .97 * label_peaks[Left_Cerebral_White_Matter]) ||
            (label_peaks[l] >= .97 * label_peaks[Right_Cerebral_White_Matter])) {
          // don't let pallidum be as bright as wm
          peak = 0.97 * (label_peaks[Left_Cerebral_White_Matter] + label_peaks[Right_Cerebral_White_Matter]) / 2;
          scale = peak / label_peaks[l];
          printf(
              "%s too bright - rescaling by %2.3f "
              "(from %2.3f) to %2.1f (was %2.1f)\n",
              cma_label_to_name(l),
              scale,
              label_scales[l],
              peak,
              label_peaks[l]);
          label_scales[l] *= scale;
          label_peaks[l] = peak;
        }
        break;
      case Left_Putamen:
      case Right_Putamen:
        if ((label_peaks[l] >= .9 * label_peaks[Left_Cerebral_White_Matter]) ||
            (label_peaks[l] >= .9 * label_peaks[Right_Cerebral_White_Matter])) {
          // don't let putamen be as bright as wm
          peak = 0.9 * (label_peaks[Left_Cerebral_White_Matter] + label_peaks[Right_Cerebral_White_Matter]) / 2;
          scale = peak / label_peaks[l];
          printf(
              "%s too bright - rescaling by %2.3f "
              "(from %2.3f) to %2.1f (was %2.1f)\n",
              cma_label_to_name(l),
              scale,
              label_scales[l],
              peak,
              label_peaks[l]);
          label_scales[l] *= scale;
          label_peaks[l] = peak;
        }
        break;
    }
  }

// now use precomputed eigenstructure of intensity
// volumes to adjust corrections

  if (logfp) {
    for (l = 0; l < MAX_CMA_LABELS; l++)
      if (computed[l] != 0)
        fprintf(logfp,
                "label %s: scaling by %2.2f  + %2.1f to %2.0f\n",
                cma_label_to_name(l),
                label_scales[l],
                label_offsets[l],
                label_peaks[l]);
    fflush(logfp);
  }
  if (DIAG_VERBOSE_ON) {
    FILE *fp;
    fp = fopen("norm_offset.plt", "w");
    for (l = 0; l < MAX_CMA_LABELS; l++)
      if (computed[l] != 0) {
        fprintf(fp, "%d %f %f\n", l, label_scales[l], label_offsets[l]);
      }
    fclose(fp);
  }

  if (base_name) {
    FILE *fp;
    char fname[STRLEN];
    sprintf(fname, "%s.label_intensities.txt", base_name);
    printf("saving intensity scales to %s\n", fname);
    fp = fopen(fname, "w");
    if (fp == NULL) {
      ErrorExit(ERROR_NOFILE, "%s: could not open intensity tracking file %s", Progname, fname);
    }

    for (l = 0; l < MAX_CMA_LABELS; l++)
      if (computed[l] != 0)
        fprintf(fp,
                "%d %s %2.2f %2.1f %2.0f\n",
                l,
                cma_label_to_name(l),
                label_scales[l],
                label_offsets[l],
                label_peaks[l]);

    fflush(fp);
    fclose(fp);
    if (getenv("EXIT_AFTER_INT") != NULL) exit(0);
  }
  gcaCheck(gca);
  GCAapplyRenormalization(gca, label_scales, label_offsets, 0);
  gcaCheck(gca);

  if (plta)  // return linear transform array to caller
  {
    int i;

    // count # of xforms
    for (i = l = 0; l < MAX_CMA_LABELS; l++) {
      if (m_by_label[l] != NULL && computed[l] == 2) {
        i++;
      }
    }

    if (i > 0)  // should always be true
    {
      *plta = lta = LTAalloc(i, mri);
      for (i = l = 0; l < MAX_CMA_LABELS; l++) {
        if (m_by_label[l] != NULL && computed[l] == 2) {
          MatrixCopy(m_by_label[l], lta->xforms[i].m_L);
          MatrixFree(&m_by_label[l]);
          lta->xforms[i].label = l;
          i++;
        }
      }
    }
    printf("%d transforms computed\n", i);
  }

  if (label_scales != plabel_scales)  // copy to user-supplied space
  {
    memmove(plabel_scales, label_scales, MAX_CMA_LABELS * sizeof(label_scales[0]));
  }
  if (label_offsets != plabel_offsets)  // copy to user-supplied space
  {
    memmove(plabel_offsets, label_offsets, MAX_CMA_LABELS * sizeof(label_offsets[0]));
  }
  if (label_peaks != plabel_peaks)  // copy to user-supplied space
  {
    memmove(plabel_peaks, label_peaks, MAX_CMA_LABELS * sizeof(label_peaks[0]));
  }
  if (computed != plabel_computed)  // copy to user-supplied space
  {
    memmove(plabel_computed, computed, MAX_CMA_LABELS * sizeof(computed[0]));
  }

  if (mri_seg) MRIfree(&mri_seg);

  return (NO_ERROR);
}

int GCAapplyRenormalization(GCA *gca, float *label_scales, float *label_offsets, int frame)
{
  int xn, yn, zn, l, i;
  GCA_NODE *gcan;
  GC1D *gc;

  for (xn = 0; xn < gca->node_width; xn++) {
    double means_before[MAX_GCA_LABELS], means_after[MAX_GCA_LABELS];
    // double scales[MAX_GCA_LABELS];
    double delta_i, delta_j;
    int xp, yp, zp;
    // int labels[MAX_GCA_LABELS];
    int niter;
    LABEL_PROB ranks_before[MAX_GCA_LABELS], ranks_after[MAX_GCA_LABELS];

    for (yn = 0; yn < gca->node_height; yn++) {
      for (zn = 0; zn < gca->node_depth; zn++) {
        if (xn == Gx && yn == Gy && zn == Gz) {
          DiagBreak();
        }
        gcan = &gca->nodes[xn][yn][zn];
        if (gcan->nlabels <= 0) {
          continue;
        }

        for (i = 0; i < gcan->nlabels; i++) {
          gc = &gcan->gcs[i];
          l = gcan->labels[i];
          // labels[i] = l;
          // scales[i] = label_scales[l];
          means_before[i] = gc->means[frame];
          ranks_before[i].label = l;
          ranks_before[i].prob = means_before[i];
          ranks_before[i].index = i;
        }
        qsort(ranks_before, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
        niter = 0;
        for (i = 0; i < gcan->nlabels; i++) {
          gc = &gcan->gcs[i];
          l = gcan->labels[i];
          means_after[i] = means_before[i] * label_scales[l] + label_offsets[l];
          if (means_after[i] < 0) {
            means_after[i] = 0;
          }
          ranks_after[i].label = l;
          ranks_after[i].prob = means_after[i];
          ranks_after[i].index = i;
        }
        qsort(ranks_after, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
        for (i = 0; i < gcan->nlabels; i++) {
          if (ranks_before[i].label != ranks_after[i].label) {
            double pi, pj, lambda;
            int j, ind_j, ind_i;
            GCA_PRIOR *gcap;

            /* two have swapped position - put them */
            /* back in the right order */
            for (j = 0; j < gcan->nlabels; j++)
              if (ranks_after[j].label == ranks_before[i].label) {
                break;
              }
            if (j >= gcan->nlabels) {
              DiagBreak();
              continue;
            }
            gcaNodeToPrior(gca, xn, yn, zn, &xp, &yp, &zp);
            gcap = &gca->priors[xp][yp][zp];
            pi = getPrior(gcap, ranks_after[i].label);
            pj = getPrior(gcap, ranks_after[j].label);
            if (FZERO(pi) && FZERO(pj)) {
              break;  // both labels will never happen
            }
            lambda = pi / (pi + pj);
            ind_j = ranks_after[j].index;
            ind_i = ranks_after[i].index;
            delta_j = (means_after[ind_j] - means_after[ind_i]) * lambda;
            delta_i = (means_after[ind_i] - means_after[ind_j]) * (1 - lambda);

            if ((fabs(delta_j) < 1) && (fabs(delta_i) < 1)) {
              // this will move one mean to the
              // other side of the other
              if ((fabs(delta_j) > fabs(delta_i)) && !FZERO(delta_j)) {
                delta_j /= fabs(delta_j);  // make it +-1
              }
              else if (!FZERO(delta_i)) {
                delta_i /= fabs(delta_i);  // make it +-1
              }
            }
            if (!std::isfinite(delta_i) || !std::isfinite(delta_j)) {
              DiagBreak();
              break;
            }
            ranks_after[j].prob = means_after[ind_j] = means_after[ind_j] - delta_j;
            ranks_after[i].prob = means_after[ind_i] = means_after[ind_i] - delta_i;
            if ((xn == Gx && yn == Gy && zn == Gz) && (ranks_after[i].label == gcan->labels[i] ||
                                                       ranks_after[j].label == gcan->labels[j] || Ggca_label < 0)) {
              printf(
                  "ordering of labels %s and %s changed, "
                  "modifying means by %2.0f (%2.1f) "
                  "and %2.0f (%2.1f)\n",
                  cma_label_to_name(ranks_after[i].label),
                  cma_label_to_name(ranks_after[j].label),
                  means_after[i],
                  delta_i,
                  means_after[j],
                  delta_i);
            }

            qsort(ranks_after, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
            i = -1; /* start loop over */
            if (niter++ > 9) {
              DiagBreak();
              break;
            }
            continue;
          }
        }

        for (i = 0; i < gcan->nlabels; i++) {
          if (FZERO(label_scales[gcan->labels[i]])) {
            continue;
          }
          gc = &gcan->gcs[i];
          if ((xn == Gx && yn == Gy && zn == Gz) && (Ggca_label == gcan->labels[i] || Ggca_label < 0)) {
            printf(
                "scaling gc for label %s at "
                "(%d, %d, %d) from %2.1f to %2.1f\n",
                cma_label_to_name(gcan->labels[i]),
                xn,
                yn,
                zn,
                means_before[i],
                means_after[i]);
            DiagBreak();
          }
          gc->means[frame] = means_after[i];
          check_finite("after rescaling", gc->means[frame]);
        }
      }
    }
  }
  gcaCheck(gca);
  return (NO_ERROR);
}

static float pthresh = 0.5;
int GCAmapRenormalize(GCA *gca, MRI *mri, TRANSFORM *transform)
{
  HISTOGRAM *h, *hsmooth;
  int l, xp, yp, zp, nbins, i, x, y, z, xn, yn, zn, num, frame, bin;
  float fmin, fmax, prior, label_scales[MAX_CMA_LABELS], label_modes[MAX_CMA_LABELS], modes[MAX_GCA_INPUTS], std, peak,
      smooth_peak;
  double val /*, scale*/;
  GCA_PRIOR *gcap;
  GCA_NODE *gcan;
  GC1D *gc;
  MATRIX *m_cov;
  MRI *mri_fsamples = NULL;  // diag volume

/* for each class, build a histogram of values
   (weighted by priors) to determine
   p(I|u,c) p(c).
*/


  for (frame = 0; frame < mri->nframes; frame++) {
    printf("renormalizing input #%d\n", frame);
    MRIvalRangeFrame(mri, &fmin, &fmax, frame);
    nbins = 256;
    h = HISTOalloc(nbins);

    hsmooth = HISTOcopy(h, NULL);
    for (l = 0; l <= MAX_CMA_LABELS; l++) /* don't do Unknown class */
    {
      label_scales[l] = 1; /* mark it as unusable */
      GCAlabelMode(gca, l, modes);
      m_cov = GCAlabelCovariance(gca, l, NULL);
      if (m_cov == NULL) {
        continue;
      }
      std = 4 * sqrt(*MATRIX_RELT(m_cov, frame + 1, frame + 1));
      MatrixFree(&m_cov);
      label_modes[l] = modes[frame];
      if (IS_UNKNOWN(l) || IS_INF_LAT_VENT(l)) {
        continue;
      }

      printf("%s (%d): mode = %2.2f +- %2.1f\n", cma_label_to_name(l), l, label_modes[l], std);
      if (l == Gdiag_no) {
        mri_fsamples = MRIclone(mri, NULL);
        DiagBreak();
      }
      if (FZERO(label_modes[l])) {
        continue;
      }
      HISTOclear(h, h);
      h->bin_size = (fmax - fmin) / 255.0;
      if (h->bin_size < 1 && (mri->type == MRI_UCHAR || mri->type == MRI_SHORT)) {
        h->bin_size = 1;
      }
      for (i = 0; i < nbins; i++) {
        h->bins[i] = (i + 1) * h->bin_size;
      }

      for (num = xp = 0; xp < gca->prior_width; xp++) {
        for (yp = 0; yp < gca->prior_height; yp++) {
          for (zp = 0; zp < gca->prior_depth; zp++) {
            if (xp == Gxp && yp == Gyp && zp == Gzp) {
              DiagBreak();
            }
            gcap = &gca->priors[xp][yp][zp];
            if (gcap == NULL) {
              continue;
            }
            prior = getPrior(gcap, l);
            if (prior < pthresh) {
              continue;
            }
            if (!GCApriorToSourceVoxel(gca, mri, transform, xp, yp, zp, &x, &y, &z)) {
              MRIsampleVolumeFrame(mri, x, y, z, frame, &val);
              if (FZERO(val))  // skull stripped
              {
                continue;
              }
              bin = nint((val - fmin) / h->bin_size);
              if (bin >= h->nbins) {
                bin = h->nbins - 1;
              }
              else if (bin < 0) {
                bin = 0;
              }

              h->counts[bin] += prior;
              num++;
              if (mri_fsamples != NULL) {
                MRIsetVoxVal(mri_fsamples, x, y, z, 0, 128);
              }
            }
          }
        }
      }
      if (num <= 50) /* not enough to reliably estimate density */
      {
        continue;
      }
      HISTOfillHoles(h);
      if (l == Gdiag_no) {
        HISTOplot(h, "h.plt");
        if (mri_fsamples && (Gdiag & DIAG_WRITE)) {
          char fname[STRLEN];
          sprintf(fname, "fsamples%d.mgz", l);
          printf("writing fsamples for class %s to %s\n", cma_label_to_name(l), fname);
          MRIwrite(mri_fsamples, fname);
          MRIfree(&mri_fsamples);
        }
        DiagBreak();
      }
      HISTOsmooth(h, hsmooth, 1);
      if (l == Gdiag_no) {
        HISTOplot(hsmooth, "hs.plt");
      }
      peak = h->bins[HISTOfindHighestPeakInRegion(h, 0, h->nbins)];
      smooth_peak = hsmooth->bins[HISTOfindHighestPeakInRegion(hsmooth, 0, hsmooth->nbins)];

      label_scales[l] = (float)smooth_peak / label_modes[l];
      printf(
          "%s (%d): peak at %2.2f, smooth at %2.2f (%d voxels), "
          "scaling by %2.2f\n",
          cma_label_to_name(l),
          l,
          peak,
          smooth_peak,
          num,
          label_scales[l]);
      bin = nint((modes[frame] - fmin) / hsmooth->bin_size);
#ifdef WSIZE
#undef WSIZE
#endif
#define WSIZE 11
#define WHALF ((WSIZE - 1) / 2)
      bin = HISTOfindCurrentPeak(hsmooth, bin, WSIZE, .2);
      smooth_peak = hsmooth->bins[bin];
      if (bin < 0 || smooth_peak <= 0) {
        continue;
      }
      if (num < 200 && hsmooth->counts[bin] < 5)
      /* not very much data - check more */
      {
        int other_bin;
        other_bin = HISTOfindPreviousPeak(hsmooth, bin, WHALF);
        if (other_bin >= 0) {
          if ((hsmooth->counts[other_bin] / hsmooth->counts[bin]) > 0.9) {
            printf(
                "!!!!!!!!!additional peak detected "
                "at %2.1f (was %2.1f) - unreliable estimate...\n",
                hsmooth->bins[other_bin],
                hsmooth->bins[bin]);
            label_scales[l] = 1.0;
            continue;
          }
        }
        other_bin = HISTOfindNextPeak(hsmooth, bin, WHALF);
        if (other_bin >= 0) {
          if (hsmooth->counts[other_bin] / hsmooth->counts[bin] > 0.9) {
            printf(
                "!!!!!!!!!additional peak detected "
                "at %2.1f (was %2.1f) - unreliable estimate...\n",
                hsmooth->bins[other_bin],
                hsmooth->bins[bin]);
            label_scales[l] = 1.0;
            continue;
          }
        }
      }
      label_scales[l] = (float)smooth_peak / label_modes[l];
      if ((label_scales[l] < 0.5 || label_scales[l] > 1.5) && !IS_LAT_VENT(l)) {
        printf("!!!!!!! rejecting excessive scaling %2.2f\n", label_scales[l]);
        label_scales[l] = 1.0;
      }
      printf(
          "%s (%d): AFTER PRIOR: peak at %2.2f, smooth "
          "at %2.2f (%d voxels), scaling by %2.2f\n",
          cma_label_to_name(l),
          l,
          peak,
          smooth_peak,
          num,
          label_scales[l]);

      if (l == Gdiag_no) {
        DiagBreak();
      }
      if (l > 100) {
        break;
      }
    }

    l = Left_Inf_Lat_Vent;
    label_scales[Left_Inf_Lat_Vent] = (.25 + .75 * label_scales[Left_Lateral_Ventricle]);
    printf(
        "%s (%d): scaling by %2.2f = %2.1f "
        "(based on %2.2f for lateral ventricle)\n",
        cma_label_to_name(l),
        l,
        label_scales[Left_Inf_Lat_Vent],
        label_modes[Left_Inf_Lat_Vent] * label_scales[Left_Inf_Lat_Vent],
        label_scales[Left_Lateral_Ventricle]);
    l = Right_Inf_Lat_Vent;
    label_scales[Right_Inf_Lat_Vent] = (.25 + .75 * label_scales[Right_Lateral_Ventricle]);
    printf(
        "%s (%d): scaling by %2.2f = %2.1f "
        "(based on %2.2f for lateral ventricle)\n",
        cma_label_to_name(l),
        l,
        label_scales[Right_Inf_Lat_Vent],
        label_modes[Right_Inf_Lat_Vent] * label_scales[Right_Inf_Lat_Vent],
        label_scales[Right_Lateral_Ventricle]);

    for (xn = 0; xn < gca->node_width; xn++) {
      double means_before[MAX_GCA_LABELS], means_after[MAX_GCA_LABELS], scales[MAX_GCA_LABELS];
      // int labels[MAX_GCA_LABELS]
      int niter;
      LABEL_PROB ranks_before[MAX_GCA_LABELS], ranks_after[MAX_GCA_LABELS];

      for (yn = 0; yn < gca->node_height; yn++) {
        for (zn = 0; zn < gca->node_depth; zn++) {
          if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
            DiagBreak();
          }
          gcan = &gca->nodes[xn][yn][zn];
          if (gcan->nlabels <= 0) {
            continue;
          }

          for (i = 0; i < gcan->nlabels; i++) {
            gc = &gcan->gcs[i];
            l = gcan->labels[i];
            // labels[i] = l;
            scales[i] = label_scales[l];
            means_before[i] = gc->means[frame];
            ranks_before[i].label = l;
            ranks_before[i].prob = means_before[i];
            ranks_before[i].index = i;
          }
          qsort(ranks_before, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
          niter = 0;
          for (i = 0; i < gcan->nlabels; i++) {
            gc = &gcan->gcs[i];
            l = gcan->labels[i];
            means_after[i] = means_before[i] * scales[i];
            ranks_after[i].label = l;
            ranks_after[i].prob = means_after[i];
            ranks_after[i].index = i;
          }
          qsort(ranks_after, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
          for (i = 0; i < gcan->nlabels; i++) {
            if (ranks_before[i].label != ranks_after[i].label) {
              double pi, pj, lambda, delta_i, delta_j;
              int j, ind_j, ind_i;
              GCA_PRIOR *gcap;

              /* two have swapped position - put them */
              /* back in the right order */
              for (j = 0; j < gcan->nlabels; j++)
                if (ranks_after[j].label == ranks_before[i].label) {
                  break;
                }
              if (j >= gcan->nlabels) {
                DiagBreak();
                continue;
              }
              gcaNodeToPrior(gca, xn, yn, zn, &xp, &yp, &zp);
              gcap = &gca->priors[xp][yp][zp];
              pi = getPrior(gcap, ranks_after[i].label);
              pj = getPrior(gcap, ranks_after[j].label);
              if (FZERO(pi) && FZERO(pj)) {
                break;  // both labels will never happen
              }
              lambda = pi / (pi + pj);
              ind_j = ranks_after[j].index;
              ind_i = ranks_after[i].index;
              delta_j = (means_after[ind_j] - means_after[ind_i]) * lambda;
              delta_i = (means_after[ind_i] - means_after[ind_j]) * (1 - lambda);

              if ((fabs(delta_j) < 1) && (fabs(delta_i) < 1)) {
                // this will move one mean to the
                // other side of the other
                if ((fabs(delta_j) > fabs(delta_i)) && !FZERO(delta_j)) {
                  delta_j /= fabs(delta_j);  // make it +-1
                }
                else if (!FZERO(delta_i)) {
                  delta_i /= fabs(delta_i);  // make it +-1
                }
              }
              if (!std::isfinite(delta_i) || !std::isfinite(delta_j)) {
                DiagBreak();
                break;
              }
              ranks_after[j].prob = means_after[ind_j] = means_after[ind_j] - delta_j;
              ranks_after[i].prob = means_after[ind_i] = means_after[ind_i] - delta_i;
              if ((xn == Gx && yn == Gy && zn == Gz) && (ranks_after[i].label == gcan->labels[i] ||
                                                         ranks_after[j].label == gcan->labels[j] || Ggca_label < 0)) {
                printf(
                    "ordering of labels %s and %s changed, "
                    "modifying means by %2.0f (%2.1f) "
                    "and %2.0f (%2.1f)\n",
                    cma_label_to_name(ranks_after[i].label),
                    cma_label_to_name(ranks_after[j].label),
                    means_after[i],
                    delta_i,
                    means_after[j],
                    delta_i);
              }

              qsort(ranks_after, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
              i = -1; /* start loop over */
              if (niter++ > 9) {
                DiagBreak();
                break;
              }
              continue;
            }
          }

          for (i = 0; i < gcan->nlabels; i++) {
            if (FZERO(label_scales[gcan->labels[i]])) {
              continue;
            }
            gc = &gcan->gcs[i];
            if ((xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) && (Ggca_label == gcan->labels[i] || Ggca_label < 0)) {
              printf(
                  "scaling gc for label %s at "
                  "(%d, %d, %d) from %2.1f to %2.1f\n",
                  cma_label_to_name(gcan->labels[i]),
                  xn,
                  yn,
                  zn,
                  means_before[i],
                  means_after[i]);
              DiagBreak();
            }
            gc->means[frame] = means_after[i];
          }
        }
      }
    }
  }

  return (NO_ERROR);
}

int GCAmapRenormalizeByClass(GCA *gca, MRI *mri, TRANSFORM *transform)
{
  HISTOGRAM *h, *hsmooth;
  int l, nbins, i, x, y, z, max_p_label, xn, yn, zn, num, frame, bin, n, c, max_label;
  // int label;
  float fmin, fmax, prior, label_scales[MAX_CMA_LABELS], class_modes[NTISSUE_CLASSES], class_scales[NTISSUE_CLASSES],
      modes[MAX_GCA_INPUTS], peak, smooth_peak;
  double val;
  float vals[MAX_GCA_INPUTS];
  GCA_PRIOR *gcap;
  GCA_NODE *gcan;
  GC1D *gc;
  double p, max_p;

  /* for each class, build a histogram of values
     (weighted by priors) to determine
     p(I|u,c) p(c).
  */

  max_label = GCAmaxLabel(gca);


  for (frame = 0; frame < mri->nframes; frame++) {
    printf("renormalizing input #%d\n", frame);
    MRIvalRangeFrame(mri, &fmin, &fmax, frame);
    nbins = 256;
    h = HISTOalloc(nbins);

    hsmooth = HISTOcopy(h, NULL);
    for (c = 0; c < NTISSUE_CLASSES; c++) /* don't do Unknown class */
    {
      class_scales[c] = 1; /* mark it as unusable */
      GCAclassMode(gca, c, modes);
      class_modes[c] = modes[frame];
      printf("%s (%d): mode = %2.2f\n", c == CSF_CLASS ? "CSF" : c == GM_CLASS ? "GM" : "WM", c, class_modes[c]);
      if (c == Gdiag_no) {
        DiagBreak();
      }
      if (FZERO(class_modes[c])) {
        continue;
      }
      HISTOclear(h, h);
      h->bin_size = (fmax - fmin) / 255.0;
      if (h->bin_size < 1 && (mri->type == MRI_UCHAR || mri->type == MRI_SHORT)) {
        h->bin_size = 1;
      }
      for (i = 0; i < nbins; i++) {
        h->bins[i] = (i + 1) * h->bin_size;
      }

      for (num = x = 0; x < mri->width; x++) {
        for (y = 0; y < mri->height; y++) {
          for (z = 0; z < mri->depth; z++) {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
              DiagBreak();
            }
            if (GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn) != NO_ERROR) {
              continue;
            }
            gcan = &gca->nodes[xn][yn][zn];
            gcap = getGCAP(gca, mri, transform, x, y, z);
            if (gcap == NULL || gcan == NULL) {
              continue;
            }
            if (gcan->nlabels == 1 && IS_UNKNOWN(gcan->labels[0])) {
              continue;
            }
            load_vals(mri, x, y, z, vals, gca->ninputs);
            val = vals[frame];
            if (FZERO(val))  // skull stripped
            {
              continue;
            }

            // find class with highest posterior likilihood
            max_p = GCAcomputePosteriorDensity(gcap, gcan, 0, -1, vals, gca->ninputs, xn, yn, zn, gca);
            max_p_label = gcan->labels[0];

            for (n = 1; n < gcan->nlabels; n++) {
              // label = gcan->labels[n];

              p = GCAcomputePosteriorDensity(gcap, gcan, n, -1, vals, gca->ninputs, xn, yn, zn, gca);
              if (p > max_p) {
                max_p = p;
                max_p_label = gcan->labels[n];
              }
            }

            if (IS_CLASS(max_p_label, c) == 0)
            // a label in a different class
            {
              continue;
            }
            prior = getPrior(gcap, max_p_label);
            if (prior < pthresh) {
              continue;
            }

            bin = nint((val - fmin) / h->bin_size);
            if (bin >= h->nbins) {
              bin = h->nbins - 1;
            }
            else if (bin < 0) {
              bin = 0;
            }

            h->counts[bin] += prior;
            num++;
          }
        }
      }
      if (num <= 50) /* not enough to reliably estimate density */
      {
        continue;
      }
      HISTOfillHoles(h);
      if (c == Gdiag_no) {
        HISTOplot(h, "h.plt");
        DiagBreak();
      }
      HISTOsmooth(h, hsmooth, 1);
      if (c == Gdiag_no) {
        HISTOplot(hsmooth, "hs.plt");
      }
      peak = h->bins[HISTOfindHighestPeakInRegion(h, 0, h->nbins)];
      smooth_peak = hsmooth->bins[HISTOfindHighestPeakInRegion(hsmooth, 0, hsmooth->nbins)];

      class_scales[c] = (float)smooth_peak / class_modes[c];
      printf(
          "%s (%d): peak at %2.2f, smooth at %2.2f (%d voxels), "
          "scaling by %2.2f\n",
          c == CSF_CLASS ? "CSF" : c == GM_CLASS ? "GM" : "WM",
          c,
          peak,
          smooth_peak,
          num,
          class_scales[c]);
      if (c == Gdiag_no) {
        DiagBreak();
      }
    }

    for (l = 0; l <= max_label; l++) {
      GCAlabelMode(gca, l, modes);
      if (FZERO(modes[frame]))  // no real data
      {
        label_scales[l] = 1.0;
        continue;
      }

// computed from manually labeled images
//(gray matter coef, wm is 1-l)
#define L_THALAMUS 0.6
#define L_PALLIDUM 0.6
#define L_PUTAMEN 0.8
#define L_VENTRALDC 0.2  // mostly white

      if (IS_GRAY_CLASS(l)) {
        label_scales[l] = class_scales[GM_CLASS];
      }
      else if (IS_CSF_CLASS(l)) {
        label_scales[l] = class_scales[CSF_CLASS];
      }
      else if (IS_WHITE_CLASS(l)) {
        label_scales[l] = class_scales[WM_CLASS];
      }
      else
        switch (l) {
          case Right_VentralDC:
          case Left_VentralDC:
            label_scales[l] = L_VENTRALDC * class_scales[GM_CLASS] + (1 - L_VENTRALDC) * class_scales[WM_CLASS];
            break;
          case Left_Thalamus:
          case Right_Thalamus:
            //      GCAlabelMean(gca, Left_Thalamus, means) ;
            //      GCAlabelMean(gca, Right_Thalamus, rmeans) ;
            //      means[frame] = (means[frame] + rmeans[frame]) / 2 ;
            label_scales[l] = L_THALAMUS * class_scales[GM_CLASS] + (1 - L_THALAMUS) * class_scales[WM_CLASS];
            break;
          case Left_Pallidum:
          case Right_Pallidum:
            label_scales[l] = L_PALLIDUM * class_scales[GM_CLASS] + (1 - L_PALLIDUM) * class_scales[WM_CLASS];
            break;
          case Left_Putamen:
          case Right_Putamen:
            label_scales[l] = (L_PUTAMEN)*class_scales[GM_CLASS] + (1 - L_PUTAMEN) * class_scales[WM_CLASS];
            break;
          case Left_Inf_Lat_Vent:  // just CSF
          case Right_Inf_Lat_Vent:
            label_scales[l] = class_scales[CSF_CLASS];
            break;
          case Left_Hippocampus:  // just GM
          case Right_Hippocampus:
          case Left_Amygdala:
          case Right_Amygdala:
          case Left_Accumbens_area:
          case Right_Accumbens_area:
          case Left_Cerebellum_Cortex:
          case Left_Cerebellum_Exterior:
          case Right_Cerebellum_Cortex:
          case Right_Cerebellum_Exterior:
          case Right_Cerebral_Exterior:
          case Left_Cerebral_Exterior:
            label_scales[l] = class_scales[GM_CLASS];
            break;
          case Right_Cerebellum_White_Matter:  // just WM
          case Left_Cerebellum_White_Matter:
          case Brain_Stem:
            label_scales[l] = class_scales[WM_CLASS];
            break;
          default:
            label_scales[l] = 1.0;
        }
      printf("%s (%d): scaling by %2.2f = %2.1f (was %2.1f)\n",
             cma_label_to_name(l),
             l,
             label_scales[l],
             modes[frame] * label_scales[l],
             modes[frame]);
    }

    for (xn = 0; xn < gca->node_width; xn++) {
      double means_before[MAX_GCA_LABELS], means_after[MAX_GCA_LABELS], scales[MAX_GCA_LABELS];
      // int labels[MAX_GCA_LABELS];
      int niter;
      LABEL_PROB ranks_before[MAX_GCA_LABELS], ranks_after[MAX_GCA_LABELS];

      for (yn = 0; yn < gca->node_height; yn++) {
        for (zn = 0; zn < gca->node_depth; zn++) {
          if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
            DiagBreak();
          }
          gcan = &gca->nodes[xn][yn][zn];
          if (gcan->nlabels <= 0) {
            continue;
          }

          for (i = 0; i < gcan->nlabels; i++) {
            gc = &gcan->gcs[i];
            l = gcan->labels[i];
            // labels[i] = l;
            scales[i] = label_scales[l];
            means_before[i] = gc->means[frame];
            ranks_before[i].label = l;
            ranks_before[i].prob = means_before[i];
            ranks_before[i].index = i;
          }
          qsort(ranks_before, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
          niter = 0;
          for (i = 0; i < gcan->nlabels; i++) {
            gc = &gcan->gcs[i];
            l = gcan->labels[i];
            means_after[i] = means_before[i] * scales[i];
            ranks_after[i].label = l;
            ranks_after[i].prob = means_after[i];
            ranks_after[i].index = i;
          }
          qsort(ranks_after, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
          for (i = 0; i < gcan->nlabels; i++) {
            if (ranks_before[i].label != ranks_after[i].label) {
              double diff, avg;
              int j;

              /* two have swapped position - put them */
              /* back in the right order */
              for (j = 0; j < gcan->nlabels; j++)
                if (ranks_after[j].label == ranks_before[i].label) {
                  break;
                }
              diff = means_before[ranks_after[i].index] - means_before[ranks_before[i].index];
              avg = (means_after[ranks_after[i].index] + means_after[ranks_before[i].index]) / 2;
              ranks_after[i].prob = means_after[ranks_after[i].index] = avg + diff / 4;
              ranks_after[j].prob = means_after[ranks_after[j].index] = avg - diff / 4;
              qsort(ranks_after, gcan->nlabels, sizeof(LABEL_PROB), compare_sort_probabilities);
              i = -1; /* start loop over */
              if (niter++ > 9) {
                DiagBreak();
                break;
              }
              continue;
            }
          }

          for (i = 0; i < gcan->nlabels; i++) {
            if (FZERO(label_scales[gcan->labels[i]])) {
              continue;
            }
            gc = &gcan->gcs[i];
            if ((xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) && (Ggca_label == gcan->labels[i] || Ggca_label < 0)) {
              printf(
                  "scaling gc for label %s at "
                  "(%d, %d, %d) from %2.1f to %2.1f\n",
                  cma_label_to_name(gcan->labels[i]),
                  xn,
                  yn,
                  zn,
                  means_before[i],
                  means_after[i]);
              DiagBreak();
            }
            gc->means[frame] = means_after[i];
          }
        }
      }
    }
  }

  return (NO_ERROR);
}

#define NLABELS 4

static int labels[NLABELS] = {Dura, Bone, SC_FAT_MUSCLE, CSF_SA};
MRI *GCArelabelNonbrain(GCA *gca, MRI *mri_inputs, MRI *mri_src, MRI *mri_dst, TRANSFORM *transform)
{
  int x, y, z, xn, yn, zn, width, height, depth, label, i, total_changed = 0, n, nchanged;
  int max_i = 0;
  GC1D *gcs[NLABELS];
  double pvals[NLABELS], max_p;
  // GCA_NODE *gcan;
  float vals[MAX_GCA_INPUTS], means[NLABELS][MAX_GCA_INPUTS], vars[NLABELS][MAX_GCA_INPUTS], dist;
  MRI *mri_tmp;

  if (mri_src != mri_dst) {
    mri_dst = MRIcopy(mri_src, mri_dst);
  }

  mri_tmp = MRIcopy(mri_dst, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (i = 0; i < NLABELS; i++) {
    GCAcomputeLabelStats(gca, labels[i], vars[i], means[i]);
  }

  /* replace Epidermis with SC_FAT */
  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        label = nint(MRIgetVoxVal(mri_src, x, y, z, 0));
        if (label == Epidermis) {
          label = SC_FAT_MUSCLE;
        }
        if (label == Cranium) {
          label = Bone;
        }
        MRIsetVoxVal(mri_src, x, y, z, 0, label);
        MRIsetVoxVal(mri_tmp, x, y, z, 0, label);
      }
    }
  }

  do {
    nchanged = 0;
    for (x = 0; x < width; x++) {
      for (y = 0; y < height; y++) {
        for (z = 0; z < depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }

          label = nint(MRIgetVoxVal(mri_tmp, x, y, z, 0));

          if (label == SC_FAT_MUSCLE)
          /* check to see whether at borders of skull */
          {
            if (MRIneighborsInWindow(mri_src, x, y, z, 3, Unknown) > 1) {
              continue;
            }
          }

          if (label != Dura && label != Bone && label != SC_FAT_MUSCLE && label != CSF_SA) {
            continue;
          }
          if (MRIneighborsInWindow(mri_src, x, y, z, 3, label) >= 24) {
            continue; /* in the body of the label - ignore */
          }
          if (MRIneighborsInWindow(mri_src, x, y, z, 5, label) >= 100) {
            continue; /* in the body of the label - ignore */
          }
          load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
          if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
            // gcan = &gca->nodes[xn][yn][zn];
            max_i = -1;
            max_p = -1;
            for (i = 0; i < NLABELS; i++) {
              gcs[i] = GCAfindGC(gca, xn, yn, zn, labels[i]);
              if (gcs[i] == NULL) gcs[i] = findGCInWindow(gca, xn, yn, zn, labels[i], 3);
              if (gcs[i] == NULL) {
                DiagBreak();
                continue;
              }
              for (dist = 0, n = 0; n < gca->ninputs; n++) {
                dist += SQR(vals[n] - means[i][n]);
              }
              pvals[i] = exp(-dist);
              if (pvals[i] >= max_p) {
                max_p = pvals[i];
                max_i = i;
              }
            }
          }
          ///////////////
          if ((labels[max_i] == Dura) && (MRIneighborsInWindow(mri_tmp, x, y, z, 5, Bone) +
                                              MRIneighborsInWindow(mri_tmp, x, y, z, 5, SC_FAT_MUSCLE) >=
                                          120)) {
            continue;
          }
          if (labels[max_i] != label && ((x == Ggca_x && y == Ggca_y && z == Ggca_z))) {
            printf(
                "GCArelabelNonbrain: changing label at "
                "(%d, %d, %d) from %s (%d) to %s (%d)\n",
                x,
                y,
                z,
                cma_label_to_name(label),
                label,
                cma_label_to_name(labels[max_i]),
                labels[max_i]);
          }
          MRIsetVoxVal(mri_tmp, x, y, z, 0, labels[max_i]);
          if (labels[max_i] != label) {
            nchanged++;
          }
        }
      }
    }
    total_changed += nchanged;
    break;
  } while (nchanged > 0);

  /* look for dura between bone and skin -
     it is partial volumed of these two */
  do {
    nchanged = 0;
    for (x = 0; x < width; x++) {
      for (y = 0; y < height; y++) {
        for (z = 0; z < depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }

          label = MRIgetVoxVal(mri_tmp, x, y, z, 0);

          if (label != Dura) {
            continue;
          }
          if (((MRIneighborsInWindow(mri_tmp, x, y, z, 3, Bone) < 7) ||
               (MRIneighborsInWindow(mri_tmp, x, y, z, 3, SC_FAT_MUSCLE) < 7)) &&
              ((MRIneighborsInWindow(mri_tmp, x, y, z, 5, Bone) < 30) ||
               (MRIneighborsInWindow(mri_tmp, x, y, z, 5, SC_FAT_MUSCLE) < 30))) {
            continue;
          }
          load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
          if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
            // gcan = &gca->nodes[xn][yn][zn];
            max_i = -1;
            max_p = -1;
            for (i = 0; i < NLABELS; i++) {
              if (labels[i] != SC_FAT_MUSCLE && labels[i] != Bone) {
                continue;
              }
              gcs[i] = GCAfindGC(gca, xn, yn, zn, labels[i]);
              if (gcs[i] == NULL) gcs[i] = findGCInWindow(gca, xn, yn, zn, labels[i], 3);
              if (gcs[i] == NULL) {
                DiagBreak();
                continue;
              }
              for (dist = 0, n = 0; n < gca->ninputs; n++) {
                dist += SQR(vals[n] - means[i][n]);
              }
              pvals[i] = exp(-dist);
              if (pvals[i] >= max_p) {
                max_p = pvals[i];
                max_i = i;
              }
            }
          }
          /////////////////////////
          if (labels[max_i] != label && ((x == Ggca_x && y == Ggca_y && z == Ggca_z))) {
            printf(
                "GCArelabelNonbrain: changing label at "
                "(%d, %d, %d) from %s (%d) to %s (%d)\n",
                x,
                y,
                z,
                cma_label_to_name(label),
                label,
                cma_label_to_name(labels[max_i]),
                labels[max_i]);
          }
          MRIsetVoxVal(mri_tmp, x, y, z, 0, labels[max_i]);
          if (labels[max_i] != label) {
            nchanged++;
          }
        }
      }
    }
    total_changed += nchanged;
    printf("%d voxels labeled dura changed to bone or skin...\n", nchanged);
    if (nchanged < 10) {
      break;
    }
  } while (nchanged > 0);

  /* change dura and cortex labels that are
     in the middle of tons of SC_FAT to SC_FAT */
  do {
    nchanged = 0;
    for (x = 0; x < width; x++) {
      for (y = 0; y < height; y++) {
        for (z = 0; z < depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }

          label = nint(MRIgetVoxVal(mri_tmp, x, y, z, 0));

          if (label != Dura && !IS_CORTEX(label)) {
            continue;
          }
          if (MRIneighborsInWindow(mri_tmp, x, y, z, 5, SC_FAT_MUSCLE) < 100) {
            continue;
          }
          label = SC_FAT_MUSCLE;
          if ((label != nint(MRIgetVoxVal(mri_tmp, x, y, z, 0))) && ((x == Ggca_x && y == Ggca_y && z == Ggca_z))) {
            printf(
                "GCArelabelNonbrain: changing label at "
                "(%d, %d, %d) from %s (%d) to %s (%d)\n",
                x,
                y,
                z,
                cma_label_to_name(nint(MRIgetVoxVal(mri_tmp, x, y, z, 0))),
                nint(MRIgetVoxVal(mri_tmp, x, y, z, 0)),
                cma_label_to_name(label),
                label);
          }
          if (label != nint(MRIgetVoxVal(mri_tmp, x, y, z, 0))) {
            nchanged++;
          }
          MRIsetVoxVal(mri_tmp, x, y, z, 0, label);
        }
      }
    }
    total_changed += nchanged;
    printf("%d voxels labeled dura and/or cortex changed to skin...\n", nchanged);
    if (nchanged < 10) {
      break;
    }
  } while (nchanged > 0);

  MRIcopy(mri_tmp, mri_dst);
  MRIfree(&mri_tmp);

  printf("%d voxels changed in MRIrelabelNonbrain\n", total_changed);
  return (mri_dst);
}

int GCAreplaceLabels(GCA *gca, int in_label, int out_label)
{
  int x, y, z, n;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++)
          if (gcan->labels[n] == in_label) {
            gcan->labels[n] = out_label;
          }
      }
    }
  }

  for (x = 0; x < gca->prior_width; x++) {
    for (y = 0; y < gca->prior_height; y++) {
      for (z = 0; z < gca->prior_depth; z++) {
        gcap = &gca->priors[x][y][z];
        if (gcap == NULL) {
          continue;
        }
        for (n = 0; n < gcap->nlabels; n++)
          if (gcap->labels[n] == in_label) {
            gcap->labels[n] = out_label;
          }
      }
    }
  }

  return (NO_ERROR);
}

int GCAreplaceRightWithLeft(GCA *gca)
{
  GCAreplaceLabels(gca, Right_Cerebral_Exterior, Left_Cerebral_Exterior);
  GCAreplaceLabels(gca, Right_Cerebral_White_Matter, Left_Cerebral_White_Matter);
  GCAreplaceLabels(gca, Right_Cerebral_Cortex, Left_Cerebral_Cortex);
  GCAreplaceLabels(gca, Right_Lateral_Ventricle, Left_Lateral_Ventricle);
  GCAreplaceLabels(gca, Right_Inf_Lat_Vent, Left_Inf_Lat_Vent);
  GCAreplaceLabels(gca, Right_Cerebellum_Exterior, Left_Cerebellum_Exterior);
  GCAreplaceLabels(gca, Right_Cerebellum_White_Matter, Left_Cerebellum_White_Matter);
  GCAreplaceLabels(gca, Right_Cerebellum_Cortex, Left_Cerebellum_Cortex);
  GCAreplaceLabels(gca, Right_Thalamus, Left_Thalamus);
  GCAreplaceLabels(gca, Right_Caudate, Left_Caudate);
  GCAreplaceLabels(gca, Right_Putamen, Left_Putamen);
  GCAreplaceLabels(gca, Right_Pallidum, Left_Pallidum);
  GCAreplaceLabels(gca, Right_Hippocampus, Left_Hippocampus);
  GCAreplaceLabels(gca, Right_Amygdala, Left_Amygdala);
  GCAreplaceLabels(gca, Right_Insula, Left_Insula);
  GCAreplaceLabels(gca, Right_Operculum, Left_Operculum);
  GCAreplaceLabels(gca, Right_Lesion, Left_Lesion);
  GCAreplaceLabels(gca, Right_Accumbens_area, Left_Accumbens_area);
  GCAreplaceLabels(gca, Right_Substancia_Nigra, Left_Substancia_Nigra);
  GCAreplaceLabels(gca, Right_VentralDC, Left_VentralDC);
  GCAreplaceLabels(gca, Right_undetermined, Left_undetermined);
  return (NO_ERROR);
}

GCA_NODE *GCAbuildRegionalGCAN(const GCA *gca, int xn, int yn, int zn, int wsize)
{
  GCA_NODE *gcan, *gcan_nbr;
  GCA_PRIOR *gcap;
  int n, xi, yi, zi, xk, yk, zk, nlabels = 0, whalf, xp, yp, zp, label, total_training[MAX_CMA_LABELS + 1],
                                 used[MAX_CMA_LABELS + 1];
  float label_priors[MAX_CMA_LABEL + 1], p;
  MATRIX *m_cov[MAX_CMA_LABEL + 1], *m_tmp = NULL;
  VECTOR *v_means[MAX_CMA_LABEL + 1], *v_tmp = NULL;

  gcan = (GCA_NODE *)calloc(1, sizeof(GCA_NODE));

  memset(label_priors, 0, sizeof(label_priors));
  memset(total_training, 0, sizeof(total_training));
  memset(used, 0, sizeof(used));
  nlabels = 0;
  whalf = (wsize - 1) / 2;
  for (xk = -whalf; xk <= whalf; xk++) {
    xi = xn + xk;
    if (xi < 0 || xi >= gca->node_width) {
      continue;
    }
    for (yk = -whalf; yk <= whalf; yk++) {
      yi = yn + yk;
      if (yi < 0 || yi >= gca->node_height) {
        continue;
      }
      for (zk = -whalf; zk <= whalf; zk++) {
        zi = zn + zk;
        if (zi < 0 || zi >= gca->node_depth) {
          continue;
        }
        gcan_nbr = &gca->nodes[xi][yi][zi];
        if (gcaNodeToPrior(gca, xi, yi, zi, &xp, &yp, &zp) != NO_ERROR) {
          continue;
        }
        gcap = &gca->priors[xp][yp][zp];
        for (n = 0; n < gcan_nbr->nlabels; n++) {
          label = gcan_nbr->labels[n];
          if (used[label] == 0) /* first time for this label */
          {
            used[label] = 1;
            m_cov[label] = load_covariance_matrix(&gcan_nbr->gcs[n], NULL, gca->ninputs);
            v_means[label] = load_mean_vector(&gcan_nbr->gcs[n], NULL, gca->ninputs);
            MatrixClear(m_cov[label]);
            MatrixClear(v_means[label]);
            nlabels++;
          }
          total_training[label] += gcan_nbr->gcs[n].ntraining;
          p = getPrior(gcap, label);
          label_priors[label] += p;

          m_tmp = load_covariance_matrix(&gcan_nbr->gcs[n], m_tmp, gca->ninputs);
          MatrixScalarMul(m_tmp, p, m_tmp);
          MatrixAdd(m_tmp, m_cov[label], m_cov[label]);

          v_tmp = load_mean_vector(&gcan_nbr->gcs[n], v_tmp, gca->ninputs);
          MatrixScalarMul(v_tmp, p, v_tmp);
          MatrixAdd(v_tmp, v_means[label], v_means[label]);
        }
        gcan->total_training += gcan_nbr->total_training;
      }
    }
  }

  gcan->nlabels = gcan->max_labels = nlabels;
  gcan->gcs = alloc_gcs(nlabels, GCA_NO_MRF, gca->ninputs);
  gcan->labels = (unsigned short *)calloc(nlabels, sizeof(unsigned short));

  for (nlabels = 0, n = 0; n <= MAX_CMA_LABELS; n++) {
    if (used[n] > 0) {
      gcan->labels[nlabels] = n;
      MatrixScalarMul(m_cov[n], 1.0 / label_priors[n], m_cov[n]);
      MatrixScalarMul(v_means[n], 1.0 / label_priors[n], v_means[n]);
      set_mean_vector(&gcan->gcs[nlabels], v_means[n], gca->ninputs);
      set_covariance_matrix(&gcan->gcs[nlabels], m_cov[n], gca->ninputs);
      MatrixFree(&m_cov[n]);
      VectorFree(&v_means[n]);
      gcan->gcs[nlabels].ntraining = total_training[n];
      nlabels++;
    }
  }

  MatrixFree(&m_tmp);
  VectorFree(&v_tmp);
  return (gcan);
}
GCA_PRIOR *GCAbuildRegionalGCAP(const GCA *gca, int xp, int yp, int zp, int wsize)
{
  GCA_PRIOR *gcap, *gcap_nbr;
  int n, xi, yi, zi, xk, yk, zk, nlabels = 0, whalf, label, total_training[MAX_CMA_LABELS + 1],
                                 used[MAX_CMA_LABELS + 1];
  float label_priors[MAX_CMA_LABEL + 1], p;
  double total_p;

  gcap = (GCA_PRIOR *)calloc(1, sizeof(GCA_PRIOR));

  memset(label_priors, 0, sizeof(label_priors));
  memset(total_training, 0, sizeof(total_training));
  memset(used, 0, sizeof(used));
  nlabels = 0;
  whalf = (wsize - 1) / 2;
  for (xk = -whalf; xk <= whalf; xk++) {
    xi = xp + xk;
    if (xi < 0 || xi >= gca->prior_width) continue;

    for (yk = -whalf; yk <= whalf; yk++) {
      yi = yp + yk;
      if (yi < 0 || yi >= gca->prior_height) continue;

      for (zk = -whalf; zk <= whalf; zk++) {
        zi = zp + zk;
        if (zi < 0 || zi >= gca->prior_depth) continue;

        gcap_nbr = &gca->priors[xi][yi][zi];
        for (n = 0; n < gcap_nbr->nlabels; n++) {
          label = gcap_nbr->labels[n];
          if (used[label] == 0) /* first time for this label */
          {
            used[label] = 1;
            nlabels++;
          }
          total_training[label] += gcap_nbr->total_training;
          p = getPrior(gcap, label);
          label_priors[label] += p;
        }
        gcap->total_training += gcap_nbr->total_training;
      }
    }
  }

  gcap->nlabels = gcap->max_labels = nlabels;
  gcap->labels = (unsigned short *)calloc(nlabels, sizeof(unsigned short));
  gcap->priors = (float *)calloc(nlabels, sizeof(float));

  for (total_p = 0., n = 0; n <= MAX_CMA_LABELS; n++) {
    if (used[n] > 0) {
      total_p += label_priors[n];
      gcap->total_training += total_training[n];
    }
  }

  for (nlabels = 0, n = 0; n <= MAX_CMA_LABELS; n++) {
    if (used[n] > 0) {
      gcap->labels[nlabels] = n;
      gcap->priors[nlabels] = label_priors[n] / total_p;
      nlabels++;
    }
  }

  return (gcap);
}

int GCAfreeRegionalGCAN(GCA_NODE **pgcan)
{
  GCA_NODE *gcan;

  gcan = *pgcan;
  *pgcan = NULL;
  free_gcs(gcan->gcs, GCA_NO_MRF, gcan->nlabels);
  free(gcan->labels);
  free(gcan);
  return (NO_ERROR);
}

GCA *GCAcompactify(GCA *gca)
{
  int width, height, depth;
  GCA_PRIOR *gcap = 0;
  GCA_NODE *gcan = 0;
  float *old_priors;
  unsigned short *old_labels;
  GC1D *old_gcs;
  int n, nmax;
  int i, j, k;
  double byteSaved = 0.;

  width = gca->prior_width;
  height = gca->prior_height;
  depth = gca->prior_depth;

  for (k = 0; k < depth; ++k)
    for (j = 0; j < height; ++j)
      for (i = 0; i < width; ++i) {
        gcap = &gca->priors[i][j][k];
        // typedef struct
        // {
        //   short nlabels ;
        //   short max_labels ;        modify
        //   unsigned short  *labels ;  modify
        //   float *priors ;           modify
        //   int   total_training ;
        // } GCA_PRIOR ;
        //
        if (gcap) {
          n = gcap->nlabels;
          nmax = gcap->max_labels;
          if (n < nmax) {
            // printf("prior has more than needed (%d,%d,%d)
            // nlabels=%d, max_labels=%d\n", i,j,k, n, nmax);
            old_priors = gcap->priors;
            old_labels = gcap->labels;
            gcap->priors = (float *)calloc(n, sizeof(float));
            if (!gcap->priors)
              ErrorExit(ERROR_NOMEMORY, "GCANupdatePriors: couldn't expand priors to %d", gcap->max_labels);
            gcap->labels = (unsigned short *)calloc(n, sizeof(unsigned short));
            if (!gcap->labels)
              ErrorExit(ERROR_NOMEMORY, "GCANupdatePriors: couldn't expand labels to %d", gcap->max_labels);
            /* copy the old ones over */
            memmove(gcap->priors, old_priors, n * sizeof(float));
            memmove(gcap->labels, old_labels, n * sizeof(unsigned short));

            /* free the old ones */
            free(old_priors);
            free(old_labels);
            gcap->max_labels = gcap->nlabels;

            byteSaved += (sizeof(float) + sizeof(unsigned short)) * (nmax - n);
          }
        }
      }

  width = gca->node_width;
  height = gca->node_height;
  depth = gca->node_depth;
  for (k = 0; k < depth; ++k)
    for (j = 0; j < height; ++j)
      for (i = 0; i < width; ++i) {
        gcan = &gca->nodes[i][j][k];
        if (gcan) {
          n = gcan->nlabels;
          nmax = gcan->max_labels;
          if (n < nmax) {
            // printf("node has more than needed (%d,%d,%d) "
            // "nlabels=%d, max_labels=%d\n", i,j,k, n, nmax);
            // typedef struct
            // {
            // int  nlabels ;
            // int  max_labels ;         modify
            // unsigned short *labels ;   modify
            // GC1D *gcs ;               modify
            // int  total_training ;
            /* total # of times this node was was accessed */
            // } GCA_NODE ;
            old_labels = gcan->labels;
            old_gcs = gcan->gcs;
            // only allocate what is needed
            gcan->gcs = alloc_gcs(n, gca->flags, gca->ninputs);
            if (!gcan->gcs) ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand gcs to %d", gcan->max_labels);
            // only allocate what is needed
            gcan->labels = (unsigned short *)calloc(n, sizeof(unsigned short));
            if (!gcan->labels)
              ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand labels to %d", gcan->max_labels);
            copy_gcs(n, old_gcs, gcan->gcs, gca->ninputs);
            memmove(gcan->labels, old_labels, n * sizeof(unsigned short));

            /* free the old ones */
            free(old_gcs);
            free(old_labels);
            gcan->max_labels = n;
            byteSaved += (sizeof(float) + sizeof(unsigned short)) * (nmax - n);
          }
        }
      }

  if (DIAG_VERBOSE_ON) {
    printf("GCAcompactify reduced the memory use by %.f bytes.\n", byteSaved);
  }

  return gca;
}

MRI *GCAreplaceImpossibleLabels(
    MRI *mri_inputs, GCA *gca, MRI *mri_in_labels, MRI *mri_out_labels, TRANSFORM *transform)
{
  int x, y, z, width, height, depth, label, xn, yn, zn, n, found, nchanged;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  float max_p, p, vals[MAX_GCA_INPUTS];

  mri_out_labels = MRIcopy(mri_in_labels, mri_out_labels);

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;
  for (nchanged = x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        label = nint(MRIgetVoxVal(mri_out_labels, x, y, z, 0));
        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          gcan = &gca->nodes[xn][yn][zn];
          gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
          if (gcap == NULL) {
            continue;
          }
          for (found = n = 0; n < gcap->nlabels; n++) {
            if (gcap->labels[n] == label) {
              found = 1;
              break;
            }
          }
          if (found) {
            continue;
          }
          load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
          nchanged++;
          max_p = GCAcomputePosteriorDensity(gcap, gcan, 0, -1, vals, gca->ninputs, xn, yn, zn, gca);
          label = gcap->labels[0];
          for (n = 1; n < gcap->nlabels; n++) {
            p = GCAcomputePosteriorDensity(gcap, gcan, n, -1, vals, gca->ninputs, xn, yn, zn, gca);
            if (p >= max_p) {
              max_p = p;
              label = gcap->labels[n];
            }
          }
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            printf("GCAreplaceImpossibleLabels: changing label at (%d, %d, %d) from %s (%d) to %s (%d)\n",
                   x,
                   y,
                   z,
                   cma_label_to_name(nint(MRIgetVoxVal(mri_out_labels, x, y, z, 0))),
                   nint(MRIgetVoxVal(mri_out_labels, x, y, z, 0)),
                   cma_label_to_name(label),
                   label);
          MRIsetVoxVal(mri_out_labels, x, y, z, 0, label);
        }
      }
    }
  }

  printf("%d impossible labels replaced...\n", nchanged);
  return (mri_out_labels);
}

double compute_partial_volume_log_posterior(GCA *gca, GCA_NODE *gcan, GCA_PRIOR *gcap, float *vals, int l1, int l2)
{
  GC1D *gc1, *gc2;
  int i;
  double p1, p2, p, alpha, u, v, p_alpha, dist;

  gc1 = gc2 = NULL;
  for (i = 0; i < gcan->nlabels; i++) {
    if (gcan->labels[i] == l1) {
      gc1 = &gcan->gcs[i];
    }
    else if (gcan->labels[i] == l2) {
      gc2 = &gcan->gcs[i];
    }
  }
  if (gc1 == NULL || gc2 == NULL) {
    return (VERY_UNLIKELY);
  }

  p1 = p2 = 0;
  for (i = 0; i < gcap->nlabels; i++) {
    if (gcap->labels[i] == l1) {
      p1 = gcap->priors[i];
    }
    else if (gcap->labels[i] == l2) {
      p2 = gcap->priors[i];
    }
  }

#define D_ALPHA 0.01
  for (p = alpha = 0.0; alpha <= 1.0; alpha += D_ALPHA) {
    u = alpha * gc1->means[0] + (1 - alpha) * gc2->means[0];
    v = alpha * gc1->covars[0] + (1 - alpha) * gc2->covars[0];
    dist = SQR(u - vals[0]) / v;
    p_alpha = 1.0 / sqrt(v) * exp(-0.5 * dist) * pow(p1, alpha) * pow(p2, 1 - alpha);
    p += (p_alpha * D_ALPHA);
  }
  return (log(p));
}

static int gcaRelabelSegment(GCA *gca, TRANSFORM *transform, MRI *mri_inputs, MRI *mri_dst, MRI_SEGMENT *mseg)
{
  int i, n, x, y, z, labels[MAX_CMA_LABEL + 1], label, max_label, old_label, debug = 0;
  double max_posterior, new_posterior;
  GCA_PRIOR *gcap;

  memset(labels, 0, sizeof(labels));

  /* build a list of all possible labels that could occur within this segment,
     and also compute current log likelihood */
  for (max_posterior = 0.0, old_label = max_label = i = 0; i < mseg->nvoxels; i++) {
    x = mseg->voxels[i].x;
    y = mseg->voxels[i].y;
    z = mseg->voxels[i].z;
    if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
      debug = 1;
    }
    max_posterior += GCAnbhdGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, PRIOR_FACTOR);

    if (i == 0) old_label = max_label = MRIgetVoxVal(mri_dst, x, y, z, 0);  // current label
    gcap = getGCAP(gca, mri_dst, transform, x, y, z);
    for (n = 0; n < gcap->nlabels; n++) {
      labels[gcap->labels[n]]++;
    }
    /* count # of possible labels in segment */
  }

  for (label = 0; label <= MAX_CMA_LABEL; label++) {
    if (labels[label] <= 0) /* not possible in this nbhd */
    {
      continue;
    }

    /* change all labels to the new one */
    for (i = 0; i < mseg->nvoxels; i++) {
      x = mseg->voxels[i].x;
      y = mseg->voxels[i].y;
      z = mseg->voxels[i].z;
      MRIsetVoxVal(mri_dst, x, y, z, 0, label);
    }
    for (new_posterior = 0.0, i = 0; i < mseg->nvoxels; i++) {
      x = mseg->voxels[i].x;
      y = mseg->voxels[i].y;
      z = mseg->voxels[i].z;
      new_posterior += GCAnbhdGibbsLogPosterior(gca, mri_dst, mri_inputs, x, y, z, transform, PRIOR_FACTOR);
    }
    if (new_posterior > max_posterior) {
      max_label = label;
      max_posterior = new_posterior;
    }
  }
  if ((old_label != max_label) && (debug || getenv("DEBUG_MRM")))
    printf(
        "changing segment at (%2.0f, %2.0f, %2.0f) "
        "from %s (%d) to %s (%d)\n",
        mseg->cx,
        mseg->cy,
        mseg->cz,
        cma_label_to_name(old_label),
        old_label,
        cma_label_to_name(max_label),
        max_label);
  for (i = 0; i < mseg->nvoxels; i++) {
    x = mseg->voxels[i].x;
    y = mseg->voxels[i].y;
    z = mseg->voxels[i].z;
    MRIsetVoxVal(mri_dst, x, y, z, 0, max_label);
  }
  return (old_label != max_label);
}

MRI *GCAmarkImpossible(GCA *gca, MRI *mri_labeled, MRI *mri_dst, TRANSFORM *transform)
{
  int x, y, z, label;

  if (mri_dst == NULL) {
    mri_dst = MRIclone(mri_labeled, NULL);
  }
  for (x = 0; x < mri_labeled->width; x++) {
    for (y = 0; y < mri_labeled->height; y++) {
      for (z = 0; z < mri_labeled->depth; z++) {
        label = MRIgetVoxVal(mri_labeled, x, y, z, 0);
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }
        if (GCAisPossible(gca, mri_labeled, label, transform, x, y, z, 1) == 0) {
          MRIsetVoxVal(mri_dst, x, y, z, 0, 1);
        }
      }
    }
  }
  return (mri_dst);
}

int GCAmaxLabel(GCA *gca)
{
  int x, y, z, max_label, n;
  GCA_PRIOR *gcap;

  max_label = 0;
  for (x = 0; x < gca->prior_width; x++) {
    for (y = 0; y < gca->prior_height; y++) {
      for (z = 0; z < gca->prior_depth; z++) {
        gcap = &gca->priors[x][y][z];
        for (n = 0; n < gcap->nlabels; n++) {
          if (gcap->labels[n] > max_label) {
            max_label = gcap->labels[n];
          }
        }
      }
    }
  }
  return (max_label);
}

MRI *GCAbuildMostLikelyVolumeForStructure(
    const GCA *gca, MRI *mri, int label, int border, TRANSFORM *transform, MRI *mri_labels)
{
  int z, width, depth, height;
  MRI *mri_tmp;
  MRI_SEGMENTATION *mriseg;
  int free_transform = 0;
#ifdef GCAbuildMostLikelyVolumeForStructure_TIMERS
  Timer tTot;
#endif

  if (transform == NULL) {
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL);
    free_transform = 1;
  }
  if (!mri) {
    mri = MRIallocSequence(gca->prior_width, gca->prior_height, gca->prior_depth, MRI_FLOAT, gca->ninputs);
    // hey create gca volume and thus copies gca prior values
    mri->xsize = gca->xsize * gca->prior_spacing;
    mri->ysize = gca->ysize * gca->prior_spacing;
    mri->zsize = gca->zsize * gca->prior_spacing;
  }
  // most likely volume should agree with direction cosines
  //  GCAcopyDCToMRI(gca, mri);

  if (mri->nframes != gca->ninputs)
    ErrorExit(ERROR_BADPARM,
              "GCAbuildMostLikelyVolume: mri->frames "
              "(%d) does not match gca->ninputs (%d)",
              mri->nframes,
              gca->ninputs);

  // mri is prior if mri = NULL
  width = mri->width;
  depth = mri->depth;
  height = mri->height;
  //#ifdef HAVE_OPENMP
  //#pragma omp parallel for if_ROMP(experimental)
  //#endif
  for (z = 0; z < depth; z++) {
    int x, y, xn, yn, zn, max_label, n, r, xp, yp, zp;
    double max_prior;
    const GCA_NODE *gcan;
    const GCA_PRIOR *gcap;
    const GC1D *gc_max;

    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        // get node value
        if (GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn) == NO_ERROR) {
          // get prior value
          if (GCAsourceVoxelToPrior(gca, mri, transform, x, y, z, &xp, &yp, &zp) == NO_ERROR) {
            gcan = &gca->nodes[xn][yn][zn];
            gcap = &gca->priors[xp][yp][zp];
            if (gcap == NULL || gcap->nlabels <= 0) {
              continue;
            }
            // initialize
            max_prior = gcap->priors[0];
            max_label = gcap->labels[0];
            gc_max = NULL;
            // prior labels
            for (n = 1; n < gcap->nlabels; n++) {
              if (gcap->priors[n] >= max_prior) {
                max_prior = gcap->priors[n];
                max_label = gcap->labels[n];
              }
            }
            if (max_label != label) {
              if (mri_labels) {
                MRIsetVoxVal(mri_labels, x, y, z, 0, 0);
              }
              for (r = 0; r < gca->ninputs; r++) {
                MRIsetVoxVal(mri, x, y, z, r, 0);
              }
              continue;
            }
            // get max_prior, max_label
            // go through node labels
            for (n = 0; n < gcan->nlabels; n++) {
              if (gcan->labels[n] == max_label) {
                gc_max = &gcan->gcs[n];
              }
            }

            if (!gc_max) {
              continue;
            }
            if (mri_labels) {
              MRIsetVoxVal(mri_labels, x, y, z, 0, label);
            }
            for (r = 0; r < gca->ninputs; r++) {
              MRIsetVoxVal(mri, x, y, z, r, gc_max->means[r]);
            }
          }
          else {
            for (r = 0; r < gca->ninputs; r++) {
              MRIsetVoxVal(mri, x, y, z, r, 0);
            }
          }
        }
        else {
          for (r = 0; r < gca->ninputs; r++) {
            MRIsetVoxVal(mri, x, y, z, r, 0);
          }
        }
      }
    }
  }

  mriseg = MRIsegment(mri, 1.0, 255.0);
  if (!mriseg) {
    ErrorPrintf(Gerror,
                "GCAmostLikelyVolumeForStructure: "
                "label %s segmentation failed",
                cma_label_to_name(label));
  }
  else if (mriseg->nsegments > 1)  // use largest segment
  {
  }
  MRIsegmentFree(&mriseg);

  // add voxels from labels on the border of this stuct
  if (border > 0) {
    mri_tmp = MRIcopy(mri, NULL);

    //#ifdef HAVE_OPENMP
    //#pragma omp parallel for if_ROMP(experimental)
    //#endif
    for (z = 0; z < depth; z++) {
      int y, xn, yn, zn, max_label, n, r, xp, yp, zp, x;
      double max_prior;
      const GCA_NODE *gcan;
      const GCA_PRIOR *gcap;
      const GC1D *gc_max;
      for (y = 0; y < height; y++)

      {
        for (x = 0; x < width; x++) {
          if (x == Gx && y == Gy && z == Gz) {
            DiagBreak();
          }
          if (MRIgetVoxVal(mri, x, y, z, 0) > 0) {
            continue;  // already filled in
          }

          if (MRIareNonzeroInNbhd(mri, (2 * border) + 1, x, y, z) == 0) {
            continue;
          }

          // get node value
          if (GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn) == NO_ERROR) {
            // get prior value
            if (GCAsourceVoxelToPrior(gca, mri, transform, x, y, z, &xp, &yp, &zp) == NO_ERROR) {
              gcan = &gca->nodes[xn][yn][zn];
              gcap = &gca->priors[xp][yp][zp];
              if (gcap == NULL || gcap->nlabels <= 0) {
                continue;
              }
              // initialize
              max_prior = gcap->priors[0];
              max_label = gcap->labels[0];
              gc_max = NULL;
              // prior labels
              for (n = 1; n < gcap->nlabels; n++) {
                if (gcap->priors[n] >= max_prior) {
                  max_prior = gcap->priors[n];
                  max_label = gcap->labels[n];
                }
              }
              // get max_prior, max_label
              // go through node labels
              for (n = 0; n < gcan->nlabels; n++) {
                if (gcan->labels[n] == max_label) {
                  gc_max = &gcan->gcs[n];
                }
              }

              if (!gc_max || max_prior < .25) {
                continue;
              }
              if (mri_labels) {
                MRIsetVoxVal(mri_labels, x, y, z, 0, max_label);
              }
              for (r = 0; r < gca->ninputs; r++) {
                MRIsetVoxVal(mri_tmp, x, y, z, r, gc_max->means[r]);
              }
            }
          }
        }
      }
    }

    MRIcopy(mri_tmp, mri);
    MRIfree(&mri_tmp);
  }
  if (free_transform) {
    TransformFree(&transform);
  }

#ifdef GCAbuildMostLikelyVolumeForStructure_TIMERS
  printf("%s: Complete in %d ms\n", __FUNCTION__, tTot.milliseconds());
#endif

  return (mri);
}

static HISTOGRAM *gcaGetLabelHistogram(GCA *gca, int label, int frame, int border)
{
  HISTOGRAM *h_gca;
  int zn, b;

  // build histogram of this label
  h_gca = HISTOalloc(256);
  for (b = 0; b < h_gca->nbins; b++) h_gca->bins[b] = b;

  // NJS: 8June2015 commented-out parallelizing this loop as testing found that
  // mri_em_register would sometimes give different results.  unknown as to why,
  // but perhaps h_gca->counts[b] += prior is not thread-safe.
  // parallelizing this loop only trimmed two seconds from em_reg.
  // ifdef HAVE_OPENMP
  // pragma omp parallel for if_ROMP(experimental)
  // endif
  for (zn = 0; zn < gca->node_depth; zn++) {
    GCA_NODE *gcan;
    GC1D *gc;
    float prior;
    int b, yn, xn, n;
    for (yn = 0; yn < gca->node_height; yn++) {
      for (xn = 0; xn < gca->node_width; xn++) {
        gcan = &gca->nodes[xn][yn][zn];
        for (n = 0; n < gcan->nlabels; n++) {
          /* find index in lookup table for this label */
          if (gcan->labels[n] != label) continue;

          gc = &gcan->gcs[n];
          prior = get_node_prior(gca, label, xn, yn, zn);
          if (prior != 0) {
            //  for (r = 0 ; r < gca->ninputs ; r++)
            {
              b = nint(gc->means[frame]);
              if (b < 0 || b >= h_gca->nbins) {
                DiagBreak();
                if (b < 0) b = 0;

                if (b >= h_gca->nbins) b = h_gca->nbins - 1;
              }
              h_gca->counts[b] += prior;
              if (!std::isfinite(gc->means[frame])) DiagBreak();
            }
          }
        }
      }
    }
  }
  return (h_gca);
}

#if INTERP_PRIOR
static float gcaComputePrior(GCA *gca, MRI *mri, TRANSFORM *transform, int x0, int y0, int z0, int label)
{
  double x, y, z, xmd, ymd, zmd, xpd, ypd, zpd, prior;
  int xm, ym, zm, xp, yp, zp;
  GCA_PRIOR *gcap;
  float total_prior;

  if (x0 == Ggca_x && y0 == Ggca_y && z0 == Ggca_z) {
    DiagBreak();
  }

  GCAsourceVoxelToPriorReal(gca, mri, transform, x0, y0, z0, &x, &y, &z);
  xm = MAX((int)x, 0);
  xp = MIN(gca->prior_width - 1, xm + 1);
  ym = MAX((int)y, 0);
  yp = MIN(gca->prior_height - 1, ym + 1);
  zm = MAX((int)z, 0);
  zp = MIN(gca->prior_depth - 1, zm + 1);

  xmd = x - (float)xm;
  ymd = y - (float)ym;
  zmd = z - (float)zm;
  xpd = (1.0f - xmd);
  ypd = (1.0f - ymd);
  zpd = (1.0f - zmd);

  gcap = &gca->priors[xp][yp][zp];
  prior = getPrior(gcap, label);
  total_prior = prior * xpd * ypd * zpd;

  gcap = &gca->priors[xp][yp][zm];
  prior = getPrior(gcap, label);
  total_prior += prior * xpd * ypd * zmd;

  gcap = &gca->priors[xp][ym][zp];
  prior = getPrior(gcap, label);
  total_prior += prior * xpd * ymd * zpd;

  gcap = &gca->priors[xp][ym][zm];
  prior = getPrior(gcap, label);
  total_prior += prior * xpd * ymd * zmd;

  gcap = &gca->priors[xm][yp][zp];
  prior = getPrior(gcap, label);
  total_prior += prior * xmd * ypd * zpd;

  gcap = &gca->priors[xm][yp][zm];
  prior = getPrior(gcap, label);
  total_prior += prior * xmd * ypd * zmd;

  gcap = &gca->priors[xm][ym][zp];
  prior = getPrior(gcap, label);
  total_prior += prior * xmd * ymd * zpd;

  gcap = &gca->priors[xm][ym][zm];
  prior = getPrior(gcap, label);
  total_prior += prior * xmd * ymd * zmd;

  return (total_prior);
}

#endif

static void set_equilavent_classes(int *equivalent_classes)
{
  int i;

  for (i = 0; i < MAX_CMA_LABELS; i++) {
    equivalent_classes[i] = i;
  }

  // CSF 1, GM 3, WM 2
  equivalent_classes[0] = 0;
  equivalent_classes[1] = 1;
  equivalent_classes[2] = 2;
  equivalent_classes[3] = 3;
  equivalent_classes[4] = 1;
  equivalent_classes[5] = 1;
  equivalent_classes[6] = 0;
  equivalent_classes[7] = 2;
  equivalent_classes[8] = 3;
  equivalent_classes[9] = 3;
  equivalent_classes[10] = 3;
  equivalent_classes[11] = 3;
  equivalent_classes[12] = 3;
  equivalent_classes[13] = 3;
  equivalent_classes[14] = 1;
  equivalent_classes[15] = 1;
  equivalent_classes[16] = 2;  // this is brain stem, really WM??
  equivalent_classes[17] = 3;
  equivalent_classes[18] = 3;
  equivalent_classes[19] = 0;  // Left_Insula
  equivalent_classes[20] = 0;  // Left_Operculum
  equivalent_classes[21] = 0;  // Line_1
  equivalent_classes[22] = 0;  // Line_2
  equivalent_classes[23] = 0;  // Line_3
  equivalent_classes[24] = 1;  // CSF
  equivalent_classes[25] = 0;  // Left_lesion
  equivalent_classes[26] = 3;  // left_accumbens_area
  equivalent_classes[27] = 0;  // Left_Substancia_Nigra
  equivalent_classes[28] = 0;  // Left_VentralDC;
  equivalent_classes[29] = 0;  // left_undetermined
  equivalent_classes[30] = 0;  // Left_vessel
  equivalent_classes[31] = 0;  // left_choroid_plexus
  equivalent_classes[32] = 0;  // lEFT_f3ORB
  equivalent_classes[33] = 0;  // Left_lOg
  equivalent_classes[34] = 0;  // Left_aOg
  equivalent_classes[35] = 0;  // Left_mOg
  equivalent_classes[36] = 0;  // Left_pOg
  equivalent_classes[37] = 0;  // Left_Stellate
  equivalent_classes[38] = 0;  // Left_Porg
  equivalent_classes[39] = 0;  // Left_Aorg
  equivalent_classes[40] = equivalent_classes[1];
  equivalent_classes[41] = equivalent_classes[2];
  equivalent_classes[42] = equivalent_classes[3];
  equivalent_classes[43] = equivalent_classes[4];
  equivalent_classes[44] = equivalent_classes[5];
  equivalent_classes[45] = equivalent_classes[6];
  equivalent_classes[46] = equivalent_classes[7];
  equivalent_classes[47] = equivalent_classes[8];
  equivalent_classes[48] = equivalent_classes[9];
  equivalent_classes[49] = equivalent_classes[10];
  equivalent_classes[50] = equivalent_classes[11];
  equivalent_classes[51] = equivalent_classes[12];
  equivalent_classes[52] = equivalent_classes[13];
  equivalent_classes[53] = equivalent_classes[17];
  equivalent_classes[54] = equivalent_classes[18];
  equivalent_classes[55] = equivalent_classes[19];
  equivalent_classes[56] = equivalent_classes[20];
  equivalent_classes[57] = equivalent_classes[25];
  equivalent_classes[58] = equivalent_classes[26];
  equivalent_classes[59] = equivalent_classes[27];
  equivalent_classes[60] = equivalent_classes[28];
  equivalent_classes[61] = equivalent_classes[29];
  equivalent_classes[62] = equivalent_classes[30];
  equivalent_classes[63] = equivalent_classes[31];  // choroid_plexus
  equivalent_classes[64] = equivalent_classes[32];
  equivalent_classes[65] = equivalent_classes[33];
  equivalent_classes[66] = equivalent_classes[34];
  equivalent_classes[67] = equivalent_classes[35];
  equivalent_classes[68] = equivalent_classes[36];
  equivalent_classes[69] = equivalent_classes[37];
  equivalent_classes[70] = equivalent_classes[38];
  equivalent_classes[71] = equivalent_classes[39];
  equivalent_classes[72] = 1;
  equivalent_classes[73] = 0;  // Left_Interior
  equivalent_classes[74] = equivalent_classes[73];
  equivalent_classes[75] = 1;
  equivalent_classes[76] = equivalent_classes[75];
  equivalent_classes[77] = 2;
  equivalent_classes[78] = 2;
  equivalent_classes[79] = equivalent_classes[78];
  equivalent_classes[80] = 2;
  equivalent_classes[81] = 2;
  equivalent_classes[82] = equivalent_classes[81];
  equivalent_classes[83] = 0;
  equivalent_classes[84] = equivalent_classes[83];
  equivalent_classes[186] = 2;
  equivalent_classes[187] = equivalent_classes[186];

  return;
}

MRI *GCAbuildMostLikelyLabelVolume(GCA *gca, MRI *mri)
{
  /* this function creates a label volume and will be used to register
     a subject's manual label to it, as a way to get linear registration
     from the subject to the gca for gca training */
  int x, y, z, xn, yn, zn, width, depth, height, n, xp, yp, zp;
  // GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  double max_prior;
  int max_label;

  if (mri == NULL) {
    mri = MRIallocSequence(gca->prior_width, gca->prior_height, gca->prior_depth, MRI_FLOAT, gca->ninputs);
    // hey create gca volume and thus copies gca prior values
    mri->xsize = gca->xsize * gca->prior_spacing;
    mri->ysize = gca->ysize * gca->prior_spacing;
    mri->zsize = gca->zsize * gca->prior_spacing;
    // most likely volume should agree with direction cosines
    GCAcopyDCToMRI(gca, mri);
  }

  width = mri->width;
  depth = mri->depth;
  height = mri->height;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        // get node value
        if (GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn) == NO_ERROR) {
          // get prior value
          if (GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp) == NO_ERROR) {
            // gcan = &gca->nodes[xn][yn][zn];
            gcap = &gca->priors[xp][yp][zp];
            if (gcap == NULL || gcap->nlabels <= 0) {
              continue;
            }
            // initialize
            max_prior = gcap->priors[0];
            max_label = gcap->labels[0];
            // prior labels
            for (n = 1; n < gcap->nlabels; n++) {
              if (gcap->priors[n] >= max_prior) {
                max_prior = gcap->priors[n];
                max_label = gcap->labels[n];
              }
            }
            MRIsetVoxVal(mri, x, y, z, 0, max_label);
          }
          else {
            MRIsetVoxVal(mri, x, y, z, 0, 0);
          }
        }
        else {
          MRIsetVoxVal(mri, x, y, z, 0, 0);
        }
      }
    }
  }

  return (mri);
}

MRI *GCAbuildMostLikelyLabelProbabilityVolume(GCA *gca)
{
  /* this function creates a label volume and will be used to register
     a subject's manual label to it, as a way to get linear registration
     from the subject to the gca for gca training */
  int x, y, z, xn, yn, zn, width, depth, height, n, xp, yp, zp;
  // GCA_NODE *gcan;
  GCA_PRIOR *gcap;
  MRI *mri;
  double max_prior;
  // int max_label;

  mri = MRIallocSequence(gca->prior_width, gca->prior_height, gca->prior_depth, MRI_FLOAT, gca->ninputs);
  // hey create gca volume and thus copies gca prior values
  mri->xsize = gca->xsize * gca->prior_spacing;
  mri->ysize = gca->ysize * gca->prior_spacing;
  mri->zsize = gca->zsize * gca->prior_spacing;
  // most likely volume should agree with direction cosines
  GCAcopyDCToMRI(gca, mri);

  width = mri->width;
  depth = mri->depth;
  height = mri->height;
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        // get node value
        if (GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn) == NO_ERROR) {
          // get prior value
          if (GCAvoxelToPrior(gca, mri, x, y, z, &xp, &yp, &zp) == NO_ERROR) {
            // gcan = &gca->nodes[xn][yn][zn];
            gcap = &gca->priors[xp][yp][zp];
            if (gcap == NULL || gcap->nlabels <= 0) {
              continue;
            }
            // initialize
            max_prior = gcap->priors[0];
            // max_label = gcap->labels[0];
            // prior labels
            for (n = 1; n < gcap->nlabels; n++) {
              if (gcap->priors[n] >= max_prior) {
                max_prior = gcap->priors[n];
                // max_label = gcap->labels[n];
              }
            }
            MRIsetVoxVal(mri, x, y, z, 0, max_prior);
          }
          else {
            MRIsetVoxVal(mri, x, y, z, 0, 0);
          }
        }
        else {
          MRIsetVoxVal(mri, x, y, z, 0, 0);
        }
      }
    }
  }

  return (mri);
}

int GCAcomputeLabelMeansAndCovariances(GCA *gca, int target_label, MATRIX **p_mcov, VECTOR **p_vmeans)
{
  int x, y, z, n, r;
  // double var;
  double dof, total_dof;
  GC1D *gc;
  GCA_NODE *gcan;
  float fval;
  MATRIX *m_cov = NULL, *m_cov_total;
  VECTOR *v_means;

  m_cov_total = MatrixAlloc(gca->ninputs, gca->ninputs, MATRIX_REAL);
  v_means = VectorAlloc(gca->ninputs, MATRIX_REAL);

  // var = 
  total_dof = 0.0;
  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];

        for (n = 0; n < gcan->nlabels; n++) {
          if (gcan->labels[n] == target_label) {
            gc = &gcan->gcs[n];
            fval = get_node_prior(gca, target_label, x, y, z);
            if (fval != 0) {
              dof = get_node_prior(gca, target_label, x, y, z) * gcan->total_training;
              for (r = 0; r < gca->ninputs; r++) {
                VECTOR_ELT(v_means, r + 1) += dof * gc->means[r];
              }
              m_cov = load_covariance_matrix(gc, m_cov, gca->ninputs);
              MatrixScalarMul(m_cov, dof, m_cov);
              MatrixAdd(m_cov, m_cov_total, m_cov_total);
              total_dof += dof;
            }
          }
        }
      }
    }
  }

  if (total_dof > 0.0) {
    MatrixScalarMul(m_cov_total, 1 / (double)total_dof, m_cov_total);
    VectorScalarMul(v_means, 1 / (double)total_dof, v_means);
  }

  *p_mcov = m_cov_total;
  *p_vmeans = v_means;
  MatrixFree(&m_cov);
  return (NO_ERROR);
}
#ifdef WSIZE
#undef WSIZE
#endif
#ifdef WHALF
#undef WHALF
#endif

#define WSIZE 5
#define WHALF ((WSIZE - 1) / 2)
static double compute_conditional_density(MATRIX *m_inv_cov, VECTOR *v_means, VECTOR *v_vals)
{
  double p, dist, det;
  int ninputs;

  ninputs = m_inv_cov->rows;

  det = MatrixDeterminant(m_inv_cov);
  dist = MatrixMahalanobisDistance(v_means, m_inv_cov, v_vals);
  p = (1.0 / (pow(2 * M_PI, ninputs / 2.0) * sqrt(1.0 / det))) * exp(-0.5 * dist);
  return (p);
}

static int load_val_vector(VECTOR *v_means, MRI *mri_inputs, int x, int y, int z)
{
  int n;

  for (n = 0; n < mri_inputs->nframes; n++) {
    VECTOR_ELT(v_means, n + 1) = MRIgetVoxVal(mri_inputs, x, y, z, n);
  }
  return (NO_ERROR);
}

MRI *GCAlabelWMandWMSAs(GCA *gca, MRI *mri_inputs, MRI *mri_src_labels, MRI *mri_dst_labels, TRANSFORM *transform)
{
  int h, wm_label, wmsa_label, x, y, z, label, nwm, nwmsa, nunknown, ngm, ncaudate, caudate_label, gm_label, n, found,
      i, NotWMSA;
  MATRIX *m_cov_wm, *m_cov_wmsa, *m_inv_cov_wmsa, *m_inv_cov_wm, *m_cov_un, *m_inv_cov_un;
  // MATRIX *m_I;
  VECTOR *v_mean_wm, *v_mean_wmsa, *v_vals, *v_dif_label, *v_dif_wmsa, *v_mean_caudate, *v_mean_un;
  double pwm, pwmsa, wmsa_dist, wm_dist, wm_mdist, wmsa_mdist;
  GCA_PRIOR *gcap;
  MRI *mri_tmp = NULL;

  mri_dst_labels = MRIcopy(mri_src_labels, mri_dst_labels);

  v_vals = VectorAlloc(mri_inputs->nframes, MATRIX_REAL);
  v_dif_label = VectorAlloc(mri_inputs->nframes, MATRIX_REAL);
  v_dif_wmsa = VectorAlloc(mri_inputs->nframes, MATRIX_REAL);
  // m_I = MatrixIdentity(mri_inputs->nframes, NULL);
  for (h = 0; h <= 1; h++) {
    if (h == 0)  // lh
    {
      wm_label = Left_Cerebral_White_Matter;
      wmsa_label = Left_WM_hypointensities;
      caudate_label = Left_Caudate;
      gm_label = Left_Cerebral_Cortex;
    }
    else {
      wm_label = Right_Cerebral_White_Matter;
      wmsa_label = Right_WM_hypointensities;
      caudate_label = Right_Caudate;
      gm_label = Right_Cerebral_Cortex;
    }

    GCAcomputeLabelMeansAndCovariances(gca, Unknown, &m_cov_un, &v_mean_un);
    GCAcomputeLabelMeansAndCovariances(gca, wm_label, &m_cov_wm, &v_mean_wm);
    GCAcomputeLabelMeansAndCovariances(gca, caudate_label, &m_cov_wm, &v_mean_caudate);
    GCAcomputeLabelMeansAndCovariances(gca, wmsa_label, &m_cov_wmsa, &v_mean_wmsa);
    m_inv_cov_wm = MatrixInverse(m_cov_wm, NULL);
    if (m_inv_cov_wm == NULL)
      ErrorExit(ERROR_BADPARM,
                "%s: could not compute inverse covariance for %s (%d)",
                Progname,
                cma_label_to_name(wm_label),
                wm_label);
    m_inv_cov_un = MatrixInverse(m_cov_un, NULL);
    if (m_inv_cov_un == NULL)
      ErrorExit(ERROR_BADPARM,
                "%s: could not compute inverse covariance for %s (%d)",
                Progname,
                cma_label_to_name(Unknown),
                Unknown);
    m_inv_cov_wmsa = MatrixInverse(m_cov_wmsa, NULL);
    if (m_inv_cov_wmsa == NULL)
      ErrorExit(ERROR_BADPARM,
                "%s: could not compute inverse covariance for %s (%d)",
                Progname,
                cma_label_to_name(wmsa_label),
                wmsa_label);

    // do max likelihood reclassification of possible wmsa voxels
    // if they are in a nbhd with likely labels
    for (x = 0; x < mri_inputs->width; x++) {
      for (y = 0; y < mri_inputs->height; y++) {
        for (z = 0; z < mri_inputs->depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          label = MRIgetVoxVal(mri_src_labels, x, y, z, 0);
          if (label != wm_label && label != wmsa_label && label != Unknown) {
            continue;
          }
          // only process it if it's in the body of the wm
          nwm = MRIlabelsInNbhd(mri_src_labels, x, y, z, WHALF, wm_label);
          nwmsa = MRIlabelsInNbhd(mri_src_labels, x, y, z, WHALF, wmsa_label);

          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            printf("(%d, %d, %d) - %s (nbrs = %d + %d = %2.2f%%)\n",
                   x,
                   y,
                   z,
                   cma_label_to_name(label),
                   nwm,
                   nwmsa,
                   (double)(nwm + nwmsa) * 100.0 / (WSIZE * WSIZE * WSIZE));
          if (label == Unknown) {
            // only unknowns that are close to wm
            if (nwm + nwmsa < 0.5 * WSIZE * WSIZE * WSIZE) {
              continue;
            }
            nunknown = MRIlabelsInNbhd(mri_src_labels, x, y, z, WHALF, Unknown);
            if (nwm + nwmsa + nunknown < 0.9 * WSIZE * WSIZE * WSIZE) {
              continue;
            }
          }
          else if (nwm + nwmsa < .9 * WSIZE * WSIZE * WSIZE)
          // somewhat arbitrary - the bulk of the nbhd
          {
            continue;
          }

          gcap = getGCAP(gca, mri_dst_labels, transform, x, y, z);
          for (found = n = 0; n < gcap->nlabels; n++)
            if ((IS_WHITE_CLASS(gcap->labels[n]) && gcap->priors[n] > 0.1) || IS_HYPO(gcap->labels[n])) {
              found = 1;
            }
          if (found == 0)  // no chance of wm or wmsa here
          {
            continue;
          }

          if (label == Unknown) {
            DiagBreak();
          }
          load_val_vector(v_vals, mri_inputs, x, y, z);
          pwm = compute_conditional_density(m_inv_cov_wm, v_mean_wm, v_vals);
          pwmsa = compute_conditional_density(m_inv_cov_wmsa, v_mean_wmsa, v_vals);
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) printf("         - pwm = %2.3e, pwmsa = %2.3e\n", pwm, pwmsa);
          if (label == wm_label && pwmsa > pwm) {
            wm_dist = VectorDistance(v_mean_wm, v_vals);
            wmsa_dist = VectorDistance(v_mean_wmsa, v_vals);
            wm_mdist = MatrixMahalanobisDistance(v_mean_wm, m_inv_cov_wm, v_vals);
            wmsa_mdist = MatrixMahalanobisDistance(v_mean_wmsa, m_inv_cov_wmsa, v_vals);
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              printf(
                  "         - wm_dist = %2.0f, "
                  "wmsa_dist = %2.0f, mdists = (%2.0f, %2.0f)\n",
                  wm_dist,
                  wmsa_dist,
                  wm_mdist,
                  wmsa_mdist);
            if ((wm_dist > wmsa_dist) && (wm_mdist > wmsa_mdist)) {
              VectorSubtract(v_vals, v_mean_wm, v_dif_label);
              VectorSubtract(v_vals, v_mean_wmsa, v_dif_wmsa);

              /* If the abs distance to wmsa is less than the abs dist
              to the label for all modes OR (the euclidan distance to
              WM is more than twice the dist to WMSA and the
              mahalanobis distance to WM is more than twice the dist
              to WMSA) then change the label to WMSA.
              */
              NotWMSA = 0;
              for (n = 0; n < mri_inputs->nframes; n++)
                if (fabs(VECTOR_ELT(v_dif_wmsa, n + 1)) >= fabs(VECTOR_ELT(v_dif_label, n + 1))) NotWMSA = 1;
              if (!NotWMSA || ((2 * wmsa_dist < wm_dist) && (2 * wmsa_mdist < wm_mdist))) label = wmsa_label;
            }
          }
          MRIsetVoxVal(mri_dst_labels, x, y, z, 0, label);
        }
      }
    }

    // now do 3 iterations of region growing
    for (i = 0; i < 3; i++) {
      mri_tmp = MRIcopy(mri_dst_labels, mri_tmp);
      for (x = 0; x < mri_inputs->width; x++) {
        for (y = 0; y < mri_inputs->height; y++) {
          for (z = 0; z < mri_inputs->depth; z++) {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
              DiagBreak();
            }
            label = MRIgetVoxVal(mri_dst_labels, x, y, z, 0);
            if (label != wm_label && label != Unknown && label != caudate_label) {
              continue;
            }
            load_val_vector(v_vals, mri_inputs, x, y, z);
            nwmsa = MRIlabelsInNbhd(mri_dst_labels, x, y, z, 1, wmsa_label);
            if (nwmsa < 1) {
              continue;
            }
            gcap = getGCAP(gca, mri_dst_labels, transform, x, y, z);
            for (found = n = 0; n < gcap->nlabels; n++)
              if ((IS_WHITE_CLASS(gcap->labels[n]) && gcap->priors[n] > 0.1) || IS_HYPO(gcap->labels[n])) {
                found = 1;
              }
            if (found == 0)  // no chance of wm or wmsa here
            {
              continue;
            }

// only process it if it's in the body of the wm
#undef WSIZE
#define WSIZE 5
#define WHALF ((WSIZE - 1) / 2)

            nwm = MRIlabelsInNbhd(mri_tmp, x, y, z, WHALF, wm_label);
            nwmsa = MRIlabelsInNbhd(mri_tmp, x, y, z, WHALF, wmsa_label);
            nunknown = MRIlabelsInNbhd(mri_tmp, x, y, z, WHALF, Unknown);
            ncaudate = MRIlabelsInNbhd(mri_tmp, x, y, z, WHALF, caudate_label);
            ngm = MRIlabelsInNbhd(mri_tmp, x, y, z, WHALF, gm_label);

            // took gm out for now
            if (ncaudate + nwm + nwmsa + nunknown < .9 * WSIZE * WSIZE * WSIZE) /* somewhat arbitrary -
                                                                                   the bulk of the nbhd */
            {
              continue;
            }
            ngm = MRIlabelsInNbhd(mri_tmp, x, y, z, 1, gm_label);
            if (ngm > 0)  // not if there are any gm nearest nbrs
            {
              continue;
            }

            if (nwm + nwmsa == 0) {
              continue;
            }
            wm_dist = VectorDistance(v_mean_wm, v_vals);
            wmsa_dist = VectorDistance(v_mean_wmsa, v_vals);
            VectorSubtract(v_vals, v_mean_wmsa, v_dif_wmsa);
            if (label == caudate_label) {
              VectorSubtract(v_vals, v_mean_caudate, v_dif_label);
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                printf(
                    "         - wm_dist = %2.0f, "
                    "wmsa_dist = %2.0f\n",
                    wm_dist,
                    wmsa_dist);

              NotWMSA = 0;
              for (n = 0; n < mri_inputs->nframes; n++)
                if (fabs(VECTOR_ELT(v_dif_wmsa, n + 1)) >= fabs(VECTOR_ELT(v_dif_label, n + 1))) NotWMSA = 1;
              if (!NotWMSA) label = wmsa_label;
            }
            else if (label == wm_label) {
              VectorSubtract(v_vals, v_mean_wm, v_dif_label);
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                printf(
                    "         - wm_dist = %2.0f, "
                    "wmsa_dist = %2.0f\n",
                    wm_dist,
                    wmsa_dist);
              NotWMSA = 0;
              for (n = 0; n < mri_inputs->nframes; n++)
                if (fabs(VECTOR_ELT(v_dif_wmsa, n + 1)) >= fabs(VECTOR_ELT(v_dif_label, n + 1))) NotWMSA = 1;
              if (!NotWMSA || (wmsa_dist * 3 < wm_dist)) label = wmsa_label;
            }
            else if (label == Unknown) {
              VectorSubtract(v_vals, v_mean_un, v_dif_label);
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                printf(
                    "         - wm_dist = %2.0f, "
                    "wmsa_dist = %2.0f\n",
                    wm_dist,
                    wmsa_dist);
              NotWMSA = 0;
              for (n = 0; n < mri_inputs->nframes; n++)
                if (fabs(VECTOR_ELT(v_dif_wmsa, n + 1)) >= fabs(VECTOR_ELT(v_dif_label, n + 1))) NotWMSA = 1;
              if (!NotWMSA) label = wmsa_label;
            }
            MRIsetVoxVal(mri_dst_labels, x, y, z, 0, label);
          }
        }
      }
    }
    MatrixFree(&m_cov_un);
    MatrixFree(&m_inv_cov_un);
    MatrixFree(&m_cov_wm);
    MatrixFree(&m_inv_cov_wm);
    MatrixFree(&m_cov_wmsa);
    MatrixFree(&m_inv_cov_wmsa);
    VectorFree(&v_mean_wm);
    VectorFree(&v_mean_wmsa);
    VectorFree(&v_mean_caudate);
    VectorFree(&v_mean_un);
  }
  VectorFree(&v_vals);
  VectorFree(&v_dif_label);
  VectorFree(&v_dif_wmsa);
  MRIfree(&mri_tmp);
  return (mri_dst_labels);
}
#ifdef WSIZE
#undef WSIZE
#endif
#define WSIZE 7
static MRI *expand_ventricle(MRI *mri_src, MRI *mri, MRI *mri_dst, float thresh, int nvox_thresh, int xdir)
{
  int nfilled, x, y, z, xmin, xmax, ymin, ymax, zmin, zmax, nvox, prev_nfilled, iter = 0;
  MRI_REGION box;
  MRI *mri_tmp = NULL;

  mri_dst = MRIcopy(mri_src, mri_dst);
  mri_tmp = MRIcopy(mri_dst, NULL);
  nfilled = 0;
  nvox_thresh = WSIZE * WSIZE * WSIZE - 2;
  do {
    MRIboundingBox(mri_dst, 0, &box);
    prev_nfilled = nfilled;
    nfilled = 0;

    if (xdir > 0) {
      xmin = MAX(0, box.x);
      xmax = MIN(mri->width - 1, box.x + box.dx);  // expand by 1
    }
    else {
      xmin = MAX(0, box.x - 1);
      xmax = MIN(mri->width - 1, box.x + box.dx - 1);
    }
    ymin = MAX(0, box.y - 1);
    ymax = MIN(mri->height - 1, box.y + box.dy);  // expand by 1
    zmin = MAX(0, box.z - 1);
    zmax = MIN(mri->depth - 1, box.z + box.dz);  // expand by 1
    for (x = xmin; x <= xmax; x++)
      for (y = ymin; y <= ymax; y++)
        for (z = zmin; z <= zmax; z++) {
          if (x == Gx && y == Gy && z == Gz) {
            DiagBreak();
          }
          if (MRIgetVoxVal(mri_dst, x, y, z, 0) > 0)  // already on
          {
            continue;
          }
          if (MRIcountNonzeroInNbhd(mri_dst, WSIZE, x, y, z) < (WSIZE * WSIZE)) {
            continue;  // doesn't nbr any enough voxels
          }
          nvox = MRIcountThreshInNbhd(mri, WSIZE, x, y, z, thresh);
          nvox = (WSIZE * WSIZE * WSIZE) - nvox;  // # below thresh instead of above it
          if (nvox < nvox_thresh) {
            continue;  // doesn't nbr any on voxels
          }
          nfilled++;
          MRIsetVoxVal(mri_tmp, x, y, z, 0, 1);
        }

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      char fname[STRLEN];
      sprintf(fname, "v%d.mgz", iter);
      MRIwrite(mri_tmp, fname);
    }
    MRIcopy(mri_tmp, mri_dst);
    iter++;
    /*    prev_nfilled = nfilled ; */  // disable check for now
  } while (nfilled > 0 && nfilled >= prev_nfilled);
  MRIfree(&mri_tmp);
  return (mri_dst);
}
static MRI *fill_ventricles(MRI *mri_seg, MRI *mri, int border, MRI *mri_dst, int label)
{
  int i, left;
  float fmin, fmax;

  for (i = 0; i < border + 2; i++) {
    if (i == 0) {
      mri_dst = MRIerode(mri_seg, mri_dst);
    }
    else {
      MRIerode(mri_dst, mri_dst);
    }
  }

  MRIbinarize(mri_dst, mri_dst, 1, 0, 1);
  MRIlabelValRange(mri, mri_dst, 1, &fmin, &fmax);

  left = label == Left_Lateral_Ventricle;
  expand_ventricle(mri_dst, mri, mri_dst, fmax, 26, left);

  return (mri_dst);
}
static int initialize_ventricle_alignment(MRI *mri_seg, MRI *mri, MATRIX *m_L, char *base_name, int border, int label)
{
  MRI *mri_vent, *mri_vent_dist, *mri_gca_vent, *mri_gca_vent_dist, *mri_mask;
  MRI_REGION box;

  printf(
      "Handling expanded ventricles... "
      "(gca::initialize_ventricle_alignment)\n");

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_seg, "s.mgz");
  }
  mri_gca_vent = MRIbinarize(mri_seg, NULL, 1, 0, 1);
  while (border-- > 0) {
    MRIerode(mri_gca_vent, mri_gca_vent);
  }

  mri_vent = fill_ventricles(mri_seg, mri, border, NULL, label);
  mri_vent_dist = MRIdistanceTransform(mri_vent, NULL, 1, mri_vent->width, DTRANS_MODE_SIGNED, NULL);
  mri_gca_vent_dist = MRIdistanceTransform(mri_gca_vent, NULL, 1, mri_gca_vent->width, DTRANS_MODE_SIGNED, NULL);
  MRIwrite(mri_gca_vent, "gv.mgz");
  MRIwrite(mri_vent, "v.mgz");
  MRIwrite(mri_gca_vent_dist, "gvd.mgz");
  MRIwrite(mri_vent_dist, "vd.mgz");
  mri_mask = MRIcopy(mri_gca_vent, NULL);
  MRIcopyLabel(mri_vent, mri_mask, 1);
  MRIboundingBox(mri_mask, 0, &box);
  box.x -= 1;
  box.y -= 1;
  box.z -= 1;
  box.dx += 2;
  box.dy += 2;
  box.dz += 2;
  MRIcropBoundingBox(mri_mask, &box);
  MRIfillBox(mri_mask, &box, 1.0);
  MRIpowellAlignImages(mri_gca_vent_dist, mri_vent_dist, m_L, NULL, NULL, mri_gca_vent, mri_vent, 1);

  MRIfree(&mri_gca_vent);
  MRIfree(&mri_vent);
  MRIfree(&mri_mask);
  MRIfree(&mri_gca_vent_dist);
  MRIfree(&mri_vent_dist);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_seg, "s.mgz");
  }
  MRIcomputeOptimalLinearXform(mri_seg,
                               mri,
                               m_L,
                               -RADIANS(20),
                               RADIANS(20),
                               .75,
                               1.5,  // allow for ventricular expansion
                               -10,
                               10,
                               3,
                               4,
                               3,
                               4,
                               base_name,
                               0);
  MRIcomputeOptimalLinearXform(mri_seg,
                               mri,
                               m_L,
                               -RADIANS(5),
                               RADIANS(5),
                               .8,
                               1.2,  // allow for ventricular expansion
                               -5,
                               5,
                               3,
                               3,
                               3,
                               3,
                               base_name,
                               0);
  return (NO_ERROR);
}

/*-------------------------------------------------------------------------
  GCAcolorTableCMA() - construct a color table from the unique entries
  in the GCA.  RGBs are random. Indices that are not represented have
  their entries in the ctab set to NULL.  Note that the names in cma.h
  do not match the names in FreeSurferLUT.txt exactly, so
  FreeSurferLUT.txt is loaded and the names are extracted rather than
  using those in cma.h.
  -------------------------------------------------------------------------*/
COLOR_TABLE *GCAcolorTableCMA(GCA *gca)
{
  int nl, n, c, r, s, nmax;
  int labelhitlist[1000];  // probably only need 256
  COLOR_TABLE *ct, *ct0;
  char ctabfile[2000];

  sprintf(ctabfile, "%s/FreeSurferColorLUT.txt", getenv("FREESURFER_HOME"));
  printf("GCAcolorTableCMA: using ctab %s\n", ctabfile);
  ct0 = CTABreadASCII(ctabfile);
  if (ct0 == NULL) {
    printf("ERROR: reading %s\n", ctabfile);
    exit(1);
  }

  // Init the hit
  for (n = 0; n < 1000; n++) {
    labelhitlist[n] = 0;
  }

  // Go thru each node
  for (c = 0; c < gca->node_width; c++) {
    for (r = 0; r < gca->node_height; r++) {
      for (s = 0; s < gca->node_depth; s++) {
        nl = gca->nodes[c][r][s].nlabels;
        // Go thru each label in the node
        for (n = 0; n < nl; n++) {
          // Get the index (labels[n] is an index, not string)
          // Mark the index as represented
          labelhitlist[gca->nodes[c][r][s].labels[n]] = 1;
        }  // n
      }    // s
    }      // r
  }        // c

  // determine the maximum index
  nmax = 0;
  for (n = 0; n < 1000; n++)
    if (labelhitlist[n]) {
      nmax = n;
    }

  ct = CTABalloc(nmax + 1);
  for (n = 0; n <= nmax; n++) {
    if (labelhitlist[n]) {
      // If this index is represented, then set up its
      // entry in the color table. The name is derived
      // from the CMA.
      ct->entries[n]->ri = nint(randomNumber(0, 255));
      ct->entries[n]->gi = nint(randomNumber(0, 255));
      ct->entries[n]->bi = nint(randomNumber(0, 255));
      ct->entries[n]->ai = 255;
      ct->entries[n]->rf = (float)ct->entries[n]->ri / 255.0;
      ct->entries[n]->gf = (float)ct->entries[n]->gi / 255.0;
      ct->entries[n]->bf = (float)ct->entries[n]->bi / 255.0;
      ct->entries[n]->af = (float)ct->entries[n]->ai / 255.0;
      // There is a mismatch between the names listed in
      // FreeSurferColorLUT.txt and the names in cma.h but the
      // indices are the same (I hope), so we need to
      // use the names from FreeSurferColorLUT.txt.
      sprintf(ct->entries[n]->name, "%s", ct0->entries[n]->name);
      // printf("%d %s %s\n", n,cma_label_to_name(n),ct->entries[n]->name);
    }
    else {
      // If this index is not represented, then free and NULL
      // its entry.
      free(ct->entries[n]);
      ct->entries[n] = NULL;
    }
  }
  CTABfree(&ct0);
  return (ct);
}
double GCAimageLogLikelihood(GCA *gca, MRI *mri_inputs, TRANSFORM *transform, int penalize_zero_brain, MRI *mri_orig)
{
  int x, y, z, xn, yn, zn, label, xp, yp, zp, num, nz, nout, n;
  GCA_PRIOR *gcap;
  GC1D *gc;
  double total_log_p, log_p = 0, prior, min_log_p;
  float vals[MAX_GCA_INPUTS], fmin, fmax;
  MRI *mri_posterior;

  if (DIAG_VERBOSE_ON)
    mri_posterior = MRIalloc(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT);
  else {
    mri_posterior = NULL;
  }

  if (Gx >= 0) {
    GCAsourceVoxelToPrior(gca, mri_inputs, transform, Gx, Gy, Gz, &Gxp, &Gyp, &Gzp);
    GCApriorToSourceVoxel(gca, mri_inputs, transform, Gxp, Gyp, Gzp, &x, &y, &z);
  }

  fmin = 100000;
  fmax = -fmin;
  for (min_log_p = 0, nz = nout = num = xp = 0; xp < gca->prior_width; xp++) {
    for (yp = 0; yp < gca->prior_height; yp++) {
      for (zp = 0; zp < gca->prior_depth; zp++) {
        if (xp == Gxp && yp == Gyp && zp == Gzp) {
          DiagBreak();
        }
        gcap = &gca->priors[xp][yp][zp];
        if (gcap == NULL || gcap->nlabels == 0) {
          continue;
        }
        min_log_p += gcap->total_training;
        label = gcapGetMaxPriorLabel(gcap, &prior);
        if (IS_BRAIN(label) && !IS_CSF(label)) {
          GCApriorToNode(gca, xp, yp, zp, &xn, &yn, &zn);
          gc = GCAfindGC(gca, xn, yn, zn, label);
          if (gc->means[0] > fmax) {
            fmax = gc->means[0];
          }
          if (gc->means[0] < fmin) {
            fmin = gc->means[0];
          }
        }
      }
    }
  }
  fprintf(stderr, "valid intensity range = [%2.0f, %2.0f]\n", fmin, fmax);

  min_log_p = log(0.01 / min_log_p);
  if (min_log_p > -50) {
    min_log_p = -50;  // to penalize enough (I know it's a hack)
  }
  for (total_log_p = 0, nz = num = xp = 0; xp < gca->prior_width; xp++) {
    for (yp = 0; yp < gca->prior_height; yp++) {
      for (zp = 0; zp < gca->prior_depth; zp++) {
        if (xp == Gxp && yp == Gyp && zp == Gzp) {
          DiagBreak();
        }
        gcap = &gca->priors[xp][yp][zp];
        if (gcap == NULL || gcap->nlabels == 0) {
          num++;
          nout++;
          total_log_p += min_log_p;
          continue;
        }
        if (GCApriorToSourceVoxel(gca, mri_inputs, transform, xp, yp, zp, &x, &y, &z) != NO_ERROR) {
          num++;
          nout++;
          total_log_p += min_log_p;
          continue;
        }
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        GCApriorToNode(gca, xp, yp, zp, &xn, &yn, &zn);
        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
        for (n = 0; n < gcap->nlabels; n++) {
          label = gcap->labels[n];
          prior = gcap->priors[n];
          gc = GCAfindGC(gca, xn, yn, zn, label);
          if (gc == NULL) {
            num++;
            nout++;
            total_log_p += min_log_p;
            continue;
          }
          //        if (IS_BRAIN(label) == 0 || FZERO(prior))
          if (FZERO(prior)) {
            continue;
          }
          num++;

          log_p = gcaComputeLogDensity(gc, vals, gca->ninputs, prior, label);
          if (penalize_zero_brain && FZERO(vals[0]) && IS_BRAIN(label) && !IS_CSF(label)) {
            if (mri_orig)  // use unstripped volume to
                           // see if this could have been brain
            {
              float ovals[MAX_GCA_INPUTS];
              load_vals(mri_orig, x, y, z, ovals, gca->ninputs);
              if (ovals[0] >= fmin && ovals[0] <= fmax)  // could have been brain
              {
                nz++;
                log_p += min_log_p;
              }
            }
            else {
              log_p += min_log_p;
              nz++;
            }
          }
        }
        if (mri_posterior) {
          MRIsetVoxVal(mri_posterior, x, y, z, 0, log_p);
        }
        if (!std::isfinite(log_p)) {
          DiagBreak();
        }
        total_log_p += log_p;
      }
    }
  }

  if (mri_posterior) {
    MRIwrite(mri_posterior, "ll.mgz");
  }

  if (penalize_zero_brain) {
    fprintf(stderr, "%d zero brain voxels\n", nz);
  }
  return (total_log_p / (double)num);
}

MRI *GCAcomputeLikelihoodImage(GCA *gca, MRI *mri_inputs, MRI *mri_labeled, TRANSFORM *transform)
{
  int x, y, z, label;
  double likelihood;
  float vals[MAX_GCA_INPUTS];
  MRI *mri_likelihood;
  GC1D *gc;

  mri_likelihood = MRIcloneDifferentType(mri_inputs, MRI_FLOAT);
  if (mri_likelihood == NULL) ErrorExit(ERROR_NOMEMORY, "%s: could allocate likelihood volume", Progname);

  for (x = 0; x < mri_inputs->width; x++)
    for (y = 0; y < mri_inputs->height; y++)
      for (z = 0; z < mri_inputs->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        label = MRIgetVoxVal(mri_labeled, x, y, z, 0);
        gc = GCAfindSourceGC(gca, mri_inputs, transform, x, y, z, label);
        if (gc == NULL) continue;  // outside of the domain of the classifiers - leave likelihood at 0
        load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
        likelihood = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, label);
        if (DZERO(likelihood)) likelihood = 1e-10;
        MRIsetVoxVal(mri_likelihood, x, y, z, 0, -log10(likelihood));
      }

  return (mri_likelihood);
}
static int gcaMaxPriorLabel(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z)
{
  int n, max_label;
  GCA_PRIOR *gcap;
  float maxp;

  gcap = getGCAP(gca, mri, transform, x, y, z);
  if (gcap->nlabels == 0) {
    return (0);
  }
  maxp = gcap->priors[0];
  max_label = gcap->labels[0];
  for (n = 1; n < gcap->nlabels; n++) {
    if (gcap->priors[n] > maxp) {
      maxp = gcap->priors[n];
      max_label = gcap->labels[n];
    }
  }

  return (max_label);
}
static int entropy_labels[] = {
    Left_Cerebral_White_Matter,
    Left_Cerebral_Cortex,
    Left_Cerebellum_White_Matter,
    Left_Cerebellum_Cortex,
    Left_Amygdala,
    Left_Hippocampus,
    Left_Thalamus,
    Left_Pallidum,
    Left_Caudate,
    Left_Putamen,
    Left_Lateral_Ventricle,
    Left_Inf_Lat_Vent,
    Left_VentralDC,
    Brain_Stem,
    Third_Ventricle,
    Fourth_Ventricle,
};

static int contra_entropy_labels[] = {
    Right_Cerebral_White_Matter,
    Right_Cerebral_Cortex,
    Right_Cerebellum_White_Matter,
    Right_Cerebellum_Cortex,
    Right_Amygdala,
    Right_Hippocampus,
    Right_Thalamus,
    Right_Pallidum,
    Right_Caudate,
    Right_Putamen,
    Right_Lateral_Ventricle,
    Right_Inf_Lat_Vent,
    Right_VentralDC,
    Brain_Stem,
    Third_Ventricle,
    Fourth_Ventricle,
};
#define NUM_ENTROPY_LABELS (sizeof(entropy_labels) / sizeof(entropy_labels[0]))
#define NUM_CONTRA_LABELS (sizeof(contra_entropy_labels) / sizeof(contra_entropy_labels[0]))

static int compute_posterior_scale_change(GCA *gca,
                                          MRI *mri,
                                          MRI *mri_aseg,
                                          TRANSFORM *transform,
                                          int *labels,
                                          int *contra_labels,
                                          float *scales,
                                          int nlabels,
                                          float step_size)
{
  float dk[NUM_ENTROPY_LABELS];
  double prior, plike, dist;
  int x, y, z, xn, yn, zn, l, num[NUM_ENTROPY_LABELS], ind;
  GCA_PRIOR *gcap;
  GC1D *gc;
  float vals[MAX_GCA_INPUTS];
  int label_indices[MAX_CMA_LABEL + 1], label;

  memset(dk, 0, sizeof(dk));
  memset(num, 0, sizeof(num));

  for (l = 0; l <= MAX_CMA_LABEL; l++) {
    label_indices[l] = -1;
  }
  for (l = 0; l < nlabels; l++) {
    label_indices[labels[l]] = l;
    label_indices[contra_labels[l]] = l;
  }

  for (x = 0; x < mri->width; x++) {
    for (y = 0; y < mri->height; y++) {
      for (z = 0; z < mri->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        load_vals(mri, x, y, z, vals, gca->ninputs);
        label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0);
        gcap = getGCAP(gca, mri, transform, x, y, z);
        ind = label_indices[label];
        if (gcap == NULL || gcap->nlabels <= 1 || ind < 0) {
          continue;
        }
        if (MRIlabelsInNbhd(mri_aseg, x, y, z, 1, label) < 26)  // avoid border voxels
        {
          continue;
        }
        if (!GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn)) {
          gc = GCAfindGC(gca, xn, yn, zn, label);
          if (gc == NULL) {
            continue;
          }
          plike = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, label);
          prior = getPrior(gcap, label);

          dist = vals[0] - (gc->means[0] * scales[ind]);
          if (FZERO(dist)) {
            continue;
          }
          num[ind]++;
          dk[ind] += plike * prior * dist / (fabs(dist));
          if (!std::isfinite(dk[ind])) {
            DiagBreak();
          }
        }
      }
    }
  }
  for (l = 0; l < nlabels; l++)
    if (num[l] > 0) {
      if (!std::isfinite(dk[l])) {
        DiagBreak();
      }
      dk[l] /= num[l];
      scales[l] += step_size * dk[l];
    }
  return (NO_ERROR);
}
int GCArenormalizeWithEntropyMinimization(GCA *gca, MRI *mri, TRANSFORM *transform, FILE *logfp)
{
  float scales[NUM_ENTROPY_LABELS], ll, last_posterior, peaks[NUM_ENTROPY_LABELS], contra_peaks[NUM_CONTRA_LABELS],
      pct_change, last_scales[NUM_ENTROPY_LABELS];
  unsigned int i;
  int done = 0, peak_bin;
  HISTOGRAM *h_gca;
  MRI *mri_aseg;

  for (i = 0; i < NUM_ENTROPY_LABELS; i++) {
    scales[i] = 1.0;
    h_gca = gcaGetLabelHistogram(gca, entropy_labels[i], 0, 0);
    peak_bin = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
    peaks[i] = h_gca->bins[peak_bin];
    HISTOfree(&h_gca);

    h_gca = gcaGetLabelHistogram(gca, contra_entropy_labels[i], 0, 0);
    peak_bin = HISTOfindHighestPeakInRegion(h_gca, 0, h_gca->nbins);
    contra_peaks[i] = h_gca->bins[peak_bin];
    HISTOfree(&h_gca);
  }

  mri_aseg = GCAlabel(mri, gca, NULL, transform);
  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    sprintf(fname, "seg%3.3d.mgz", 0);
    printf("writing current segmentation snapshot to %s\n", fname);
    MRIwrite(mri_aseg, fname);
  }
  last_posterior = ll = GCAimagePosteriorLogProbability(gca, mri_aseg, mri, transform);
  printf("%3.3d: ll = %2.7f\n", 0, ll);
  i = 0;
  do {
    memmove(last_scales, scales, sizeof(scales));
    compute_posterior_scale_change(
        gca, mri, mri_aseg, transform, entropy_labels, contra_entropy_labels, scales, NUM_ENTROPY_LABELS, .1);
    gcaScale(gca, entropy_labels, contra_entropy_labels, scales, NUM_ENTROPY_LABELS, 1);
    {
      for (unsigned int n = 0; n < NUM_ENTROPY_LABELS; n++)
        printf("scales[%s] = %2.3f, peak = %2.0f (rh=%2.0f)\n",
               cma_label_to_name(entropy_labels[n]),
               scales[n],
               peaks[n] * scales[n],
               contra_peaks[n] * scales[n]);
    }
    ll = GCAimagePosteriorLogProbability(gca, mri_aseg, mri, transform);
    gcaScale(gca, entropy_labels, contra_entropy_labels, scales, NUM_ENTROPY_LABELS, -1);

    pct_change = 100 * (last_posterior - ll) / last_posterior;
    if (pct_change < 0.01) {
      done = 1;
    }
    i++;
    printf("%3.3d: ll = %2.7f (%2.3f%%)\n", i, ll, pct_change);
    if (logfp) {
      fprintf(logfp, "%3.3d: ll = %2.7f (%2.3f%%)\n", i, ll, pct_change);
    }
    if (!((i + 1) % 1)) {
      printf("recomputing MAP labels...\n");
      gcaScale(gca, entropy_labels, contra_entropy_labels, scales, NUM_ENTROPY_LABELS, 1);
      GCAlabel(mri, gca, mri_aseg, transform);
      if (Gdiag & DIAG_WRITE && (!((i + 1) % 2))) {
        char fname[STRLEN];
        sprintf(fname, "seg%3.3d.mgz", i + 1);
        printf("writing current segmentation snapshot to %s\n", fname);
        MRIwrite(mri_aseg, fname);
      }
      ll = GCAimagePosteriorLogProbability(gca, mri_aseg, mri, transform);
      gcaScale(gca, entropy_labels, contra_entropy_labels, scales, NUM_ENTROPY_LABELS, -1);
    }
    if (last_posterior < ll) {
      memmove(scales, last_scales, sizeof(scales));
    }

    last_posterior = ll;
    if (i < 8) {
      done = 0;
    }
    if (logfp) {
      fflush(logfp);
    }
  } while (!done);

  for (i = 0; i < NUM_ENTROPY_LABELS; i++) {
    printf("scaling %s by %2.3f from %2.0f to %2.0f (rh=%2.0f)\n",
           cma_label_to_name(entropy_labels[i]),
           scales[i],
           peaks[i],
           peaks[i] * scales[i],
           contra_peaks[i] * scales[i]);
    if (logfp)
      fprintf(logfp,
              "scaling %s by %2.3f from %2.0f to %2.0f (rh=%2.0f)\n",
              cma_label_to_name(entropy_labels[i]),
              scales[i],
              peaks[i],
              peaks[i] * scales[i],
              contra_peaks[i] * scales[i]);
  }
  if (logfp) {
    fflush(logfp);
  }
  gcaScale(gca, entropy_labels, contra_entropy_labels, scales, NUM_ENTROPY_LABELS, 1);
  return (NO_ERROR);
}

double GCAcomputeMeanEntropy(GCA *gca, MRI *mri, TRANSFORM *transform)
{
  double entropy, entropy_total, p[MAX_LABELS_PER_GCAN], ptotal, max_p;
  int x, y, z, c, xn, yn, zn, num;
  // int max_c;
  GCA_PRIOR *gcap;
  GC1D *gc;
  float vals[MAX_GCA_INPUTS];

  num = 0;
  for (entropy_total = 0.0, x = 0; x < mri->width; x++) {
    for (y = 0; y < mri->height; y++) {
      for (z = 0; z < mri->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        gcap = getGCAP(gca, mri, transform, x, y, z);
        if (gcap == NULL || gcap->nlabels <= 1) {
          continue;
        }
        if (!GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn)) {
          load_vals(mri, x, y, z, vals, gca->ninputs);
          max_p = 0;
          // max_c = 0;
          for (ptotal = 0, c = 0; c < gcap->nlabels; c++) {
            gc = GCAfindGC(gca, xn, yn, zn, gcap->labels[c]);
            if (gc == NULL) {
              p[c] = 0;
              continue;
            }
            p[c] = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, gcap->labels[c]);
            ptotal += p[c];
            if (p[c] > max_p) {
              max_p = p[c];
              // max_c = c;
            }
          }
          if (FZERO(ptotal)) {
            continue;
          }
          for (entropy = 0.0, c = 0; c < gcap->nlabels; c++) {
            p[c] /= ptotal;
            if (DZERO(p[c])) {
              continue;
            }
            entropy += p[c] * log(p[c]);
          }
          num++;
          if (!std::isfinite(entropy)) {
            DiagBreak();
          }
          entropy_total += entropy;
        }
      }
    }
  }

  return (-entropy_total / num);
}
static int gcaScale(GCA *gca, int *labels, int *contra_labels, float *scales, int nlabels, int dir)
{
  int label_indices[MAX_CMA_LABEL + 1], x, y, z, n, label, ind, i;
  GC1D *gc;
  GCA_NODE *gcan;
  float scale;

  for (i = 0; i <= MAX_CMA_LABEL; i++) {
    label_indices[i] = -1;
  }
  for (i = 0; i < nlabels; i++) {
    label_indices[labels[i]] = i;
    label_indices[contra_labels[i]] = i;
  }

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          label = gcan->labels[n];
          ind = label_indices[label];
          if (ind < 0) {
            continue;
          }
          scale = scales[ind];
          gc = &gcan->gcs[n];
          for (i = 0; i < gca->ninputs; i++) {
            if (dir > 0) {
              gc->means[i] *= scale;
            }
            else {
              gc->means[i] /= scale;
            }
          }
        }
      }
    }
  }
  return (NO_ERROR);
}

int GCAgetMaxPriorLabelAtVoxel(GCA *gca, MRI *mri, int x, int y, int z, TRANSFORM *transform, double *p_prior)
{
  int xp, yp, zp;
  double xpf, ypf, zpf;

  if (!GCAsourceVoxelToPrior(gca, mri, transform, x, y, z, &xp, &yp, &zp)) {
    xpf = xp;
    ypf = yp;
    zpf = zp;
    xp = nint(xpf);
    yp = nint(ypf);
    zp = nint(zpf);
    return (GCAgetMaxPriorLabel(gca, xp, yp, zp, p_prior));
  }
  return (-1);
}

int GCAgetMaxPriorLabel(GCA *gca, int xp, int yp, int zp, double *p_prior)
{
  GCA_PRIOR *gcap;
  gcap = &gca->priors[xp][yp][zp];
  return (gcapGetMaxPriorLabel(gcap, p_prior));
}
static int gcapGetMaxPriorLabel(GCA_PRIOR *gcap, double *p_prior)
{
  int n, best_label = 0;
  float max_p;

  if (p_prior) {
    *p_prior = 0;
  }
  if (gcap == NULL) {
    return (0);
  }

  max_p = -1;
  for (n = 0; n < gcap->nlabels; n++) {
    if (gcap->priors[n] > max_p) {
      max_p = gcap->priors[n];
      best_label = gcap->labels[n];
    }
  }
  if (p_prior) {
    *p_prior = max_p;
  }
  return (best_label);
}

float GCAcomputeLabelPosterior(GCA *gca, TRANSFORM *transform, MRI *mri, float x, float y, float z, int label)
{
  float p, plabel, vals[MAX_GCA_INPUTS], ptotal;
  GCA_PRIOR *gcap;
  int n, olabel;

  load_vals(mri, x, y, z, vals, gca->ninputs);
  gcap = getGCAPfloat(gca, mri, transform, x, y, z);
  if (gcap == NULL) {
    return (0.0);
  }

  if (gcap->total_training > 0) {
    plabel = .01 / gcap->total_training;
  }
  else {
    plabel = 0.0;
  }
  for (ptotal = 0.0, n = 0; n < gcap->nlabels; n++) {
    olabel = gcap->labels[n];
    p = GCAlabelProbability(mri, gca, transform, x, y, z, olabel) * gcap->priors[n];
    if (olabel == label) {
      plabel = p;
    }

    ptotal += p;
  }

  if (!FZERO(ptotal)) {
    plabel /= ptotal;
  }
  return (plabel);
}

float GCAcomputeLabelLikelihood(GCA *gca, TRANSFORM *transform, MRI *mri, float x, float y, float z, int label)
{
  float p, plabel, vals[MAX_GCA_INPUTS], ptotal;
  GCA_PRIOR *gcap;
  int n, olabel, found = 0;

  gcap = getGCAPfloat(gca, mri, transform, x, y, z);
  if (gcap == NULL) {
    return (0.0);
  }

  if (gcap->total_training > 0) {
    plabel = .01 / gcap->total_training;
  }
  else {
    plabel = 0.0;
  }
  for (ptotal = 0.0, n = 0; n < gcap->nlabels; n++) {
    olabel = gcap->labels[n];
    if (olabel == 9)  // debugging!!! Should be removed
    {
      continue;
    }
    p = GCAlabelProbability(mri, gca, transform, x, y, z, olabel);
    if (olabel == label) {
      found = 1;
      plabel = p;
    }

    ptotal += p;
  }

  if (found == 0) {
    GC1D *gc;
    int xn, yn, zn;

    if (!GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn)) {
      gc = GCAfindClosestValidGC(gca, xn, yn, zn, label, 0);
      if (gc == NULL) {
        DiagBreak();
      }
      else {
        load_vals(mri, x, y, z, vals, gca->ninputs);
        plabel = GCAcomputeConditionalDensity(gc, vals, gca->ninputs, label);
        ptotal += plabel;
      }
    }
  }

  if (!FZERO(ptotal)) {
    plabel /= ptotal;
  }
  return (plabel);
}

int GCAstructureBoundingBox(GCA *gca, int label, MRI_REGION *box)
{
  int x, y, z, xmin, ymin, zmin, xmax, ymax, zmax, n;
  GCA_PRIOR *gcap;

  xmin = gca->prior_width;
  ymin = gca->prior_height;
  zmin = gca->prior_depth;
  xmax = ymax = zmax = -1;

  for (x = 0; x < gca->prior_width; x++)
    for (y = 0; y < gca->prior_height; y++)
      for (z = 0; z < gca->prior_depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        gcap = &gca->priors[x][y][z];
        for (n = 0; n < gcap->nlabels; n++)
          if (gcap->labels[n] == label && gcap->priors[n] > 0.05) {
            break;
          }
        if (n < gcap->nlabels)  // label exists at this location
        {
          if (y < 60) {
            DiagBreak();
          }
          if (z < 42) {
            DiagBreak();
          }
          xmin = MIN(xmin, x);
          ymin = MIN(ymin, y);
          zmin = MIN(zmin, z);
          xmax = MAX(xmax, x);
          ymax = MAX(ymax, y);
          zmax = MAX(zmax, z);
        }
      }
  xmin *= gca->prior_spacing;
  ymin *= gca->prior_spacing;
  zmin *= gca->prior_spacing;
  xmax *= gca->prior_spacing;
  ymax *= gca->prior_spacing;
  zmax *= gca->prior_spacing;
  box->x = xmin;
  box->y = ymin;
  box->z = zmin;
  box->dx = xmax - xmin + 1;
  box->dy = ymax - ymin + 1;
  box->dz = zmax - zmin + 1;
  return (NO_ERROR);
}
int GCArenormalizeClass(GCA *gca, int classnum, float scale_to_wm)
{
  float wm_mode, class_mode, scale;
  int x, y, z, n, same_class, r;
  GCA_NODE *gcan;

  GCAclassMode(gca, WM_CLASS, &wm_mode);
  GCAclassMode(gca, classnum, &class_mode);

  scale = scale_to_wm * wm_mode / class_mode;
  printf("current wm/csf ratio = %2.2f (%2.0f/%2.0f), scaling by %2.3f\n",
         wm_mode / class_mode,
         wm_mode,
         class_mode,
         scale);
  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          switch (classnum) {
            case CSF_CLASS:
              same_class = IS_FLUID(gcan->labels[n]);
              break;
            case GM_CLASS:
              same_class = IS_GRAY_MATTER(gcan->labels[n]);
              break;
            case WM_CLASS:
              same_class = IS_WHITE_MATTER(gcan->labels[n]);
              break;
            default:
              same_class = 0;
              break;
          }
          if (same_class) {
            for (r = 0; r < gca->ninputs; r++) {
              gcan->gcs[n].means[r] *= scale;
            }
          }
        }
      }
    }
  }

  return (NO_ERROR);
}

static int nbr_labels[GIBBS_NEIGHBORS][MAX_CMA_LABELS][MAX_CMA_LABELS];
static GCA_NODE *gcaBuildNbhdGCAN(GCA *gca, int x, int y, int z, float sigma)
{
  GCA_NODE *gcan, *gcan_total;
  int labels[MAX_CMA_LABELS], n, xi, yi, zi, xk, yk, zk, whalf, nlabels, n2, n3, n4, r, c, v, i;
  GC1D *gc, *gc_total;
  double norm, wt, two_sigma_sq, dsq, total_wt, gc_training[MAX_LABELS_PER_GCAN], total_training;
  static int max_label = -1;
  int wsize =   (int)nint(3 * 2 * sigma) + 1;

  sigma /= gca->node_spacing ;  // to node vox stride
  wsize =   1+2*(int)ceil(3 * sigma);
  total_training = 0.0;
  memset(gc_training, 0, sizeof(gc_training));

  two_sigma_sq = (2 * sigma * sigma);
  if (max_label < 0) {
    max_label = GCAmaxLabel(gca);
  }

  whalf = ceil((float)(wsize - 1) / (2.0)) ;
  memset(labels, 0, sizeof(labels));
  //  memset(nbr_labels, 0, sizeof(nbr_labels)) ;

  gcan_total = (GCA_NODE *)calloc(1, sizeof(GCA_NODE));
  if (gcan_total == NULL) {
    ErrorExit(ERROR_NOMEMORY, "gcaBuildNbhdGCAN: could not allocate node");
  }

  // first count # of different labels in nbhd
  for (total_wt = 0.0, xi = -whalf; xi <= whalf; xi++) {
    xk = x + xi;
    if (xk < 0 || xk >= gca->node_width) {
      continue;
    }
    for (yi = -whalf; yi <= whalf; yi++) {
      yk = y + yi;
      if (yk < 0 || yk >= gca->node_height) {
        continue;
      }
      for (zi = -whalf; zi <= whalf; zi++) {
        zk = z + zi;
        if (zk < 0 || zk >= gca->node_depth) {
          continue;
        }
        gcan = &gca->nodes[xk][yk][zk];
        dsq = xi * xi + yi * yi + zi * zi;
        wt = exp(-dsq / two_sigma_sq);
        if (FZERO(wt)) {
          continue;
        }
        total_wt += wt;
      }
    }
  }
  for (xi = -whalf; xi <= whalf; xi++) {
    xk = x + xi;
    if (xk < 0 || xk >= gca->node_width) {
      continue;
    }
    for (yi = -whalf; yi <= whalf; yi++) {
      yk = y + yi;
      if (yk < 0 || yk >= gca->node_height) {
        continue;
      }
      for (zi = -whalf; zi <= whalf; zi++) {
        zk = z + zi;
        if (zk < 0 || zk >= gca->node_depth) {
          continue;
        }
        gcan = &gca->nodes[xk][yk][zk];
        dsq = xi * xi + yi * yi + zi * zi;
        wt = exp(-dsq / two_sigma_sq) / total_wt;
        if (FZERO(wt)) {
          continue;
        }
        for (n = 0; n < gcan->nlabels; n++)
          if (gcan->gcs[n].ntraining > 0) {
            labels[gcan->labels[n]]++;
            if (!(gca->flags & GCA_NO_MRF)) {
              for (i = 0; i < GIBBS_NEIGHBORS; i++) {
                for (n2 = 0; n2 < gcan->gcs[n].nlabels[i]; n2++) {
                  nbr_labels[i][gcan->labels[n]][gcan->gcs[n].labels[i][n2]] = 1;
                }
              }
            }
          }
      }
    }
  }
  for (n = nlabels = 0; n <= max_label; n++) {
    if (labels[n] > 0) {
      nlabels++;
    }
  }
  if (nlabels > MAX_LABELS_PER_GCAN)
    ErrorExit(ERROR_NOMEMORY,
              "gcaBuildNbhdGCAN(%d, %d, %d, %d, %2.3f): too many nbrs %d > %d",
              x,
              y,
              z,
              wsize,
              sigma,
              nlabels,
              MAX_LABELS_PER_GCAN);

  gcan_total->labels = (unsigned short *)calloc(nlabels, sizeof(unsigned short));
  gcan_total->nlabels = nlabels;
  if (gcan_total->labels == NULL) {
    ErrorExit(ERROR_NOMEMORY, "gcaBuildNbhdGCAN: could not allocate node labels");
  }

  gcan_total->gcs = alloc_gcs(nlabels, gca->flags, gca->ninputs);
  if (gcan_total->gcs == NULL) {
    ErrorExit(ERROR_NOMEMORY, "gcaBuildNbhdGCAN: could not allocate node gcs");
  }
  for (n = nlabels = 0; n <= max_label; n++) {
    if (labels[n] > 0) {
      gc_total = &gcan_total->gcs[nlabels];
      gcan_total->labels[nlabels++] = n;
      if (!(gca->flags & GCA_NO_MRF)) {
        for (i = 0; i < GIBBS_NEIGHBORS; i++) {
          for (n2 = 0; n2 <= max_label; n2++) {
            if (nbr_labels[i][n][n2]) {
              gc_total->nlabels[i]++;
            }
          }
        }
      }
    }
  }

  // now allocate label_priors[i] and labels[i] since we know nlabels[i]
  for (n = 0; n < gcan_total->nlabels; n++) {
    gc_total = &gcan_total->gcs[n];
    if (!(gca->flags & GCA_NO_MRF)) {
      for (i = 0; i < GIBBS_NEIGHBORS; i++) {
        if (gc_total->nlabels[i] == 0) {
          DiagBreak();
        }
        else {
          gc_total->labels[i] = (unsigned short *)calloc(gc_total->nlabels[i], sizeof(short));
          if (gc_total->labels[i] == NULL)
            ErrorExit(ERROR_NOMEMORY, "gcaBuildNbhdGCAN: could not allocate %d:%d labels", gc_total->nlabels[i], i);

          gc_total->label_priors[i] = (float *)calloc(gc_total->nlabels[i], sizeof(float));
          if (gc_total->label_priors[i] == NULL)
            ErrorExit(
                ERROR_NOMEMORY, "gcaBuildNbhdGCAN: could not allocate %d:%d label priors", gc_total->nlabels[i], i);
        }
        for (n2 = nlabels = 0; n2 <= max_label; n2++) {
          if (nbr_labels[i][gcan_total->labels[n]][n2]) {
            if (nlabels >= gc_total->nlabels[i]) {
              DiagBreak();
            }
            else {
              gc_total->labels[i][nlabels++] = n2;  // label n2 is nbr in position i
            }
          }
        }
      }
    }
  }

  // accumulate statistics over the nbhd
  for (xi = -whalf; xi <= whalf; xi++) {
    xk = x + xi;
    if (xk < 0 || xk >= gca->node_width) {
      continue;
    }
    for (yi = -whalf; yi <= whalf; yi++) {
      yk = y + yi;
      if (yk < 0 || yk >= gca->node_height) {
        continue;
      }
      for (zi = -whalf; zi <= whalf; zi++) {
        zk = z + zi;
        if (zk < 0 || zk >= gca->node_depth) {
          continue;
        }
        gcan = &gca->nodes[xk][yk][zk];
        dsq = xi * xi + yi * yi + zi * zi;
        wt = exp(-dsq / two_sigma_sq) / total_wt;
        if (FZERO(wt)) {
          continue;
        }
        total_training += wt * gcan->total_training;
        for (n = 0; n < gcan->nlabels; n++) {
          for (n2 = 0; n2 < gcan_total->nlabels; n2++)
            if (gcan_total->labels[n2] == gcan->labels[n]) {
              break;  // found the corresponding index
            }
          if (n2 >= gcan_total->nlabels) {
            continue;  // should never happen
          }
          gc = &gcan->gcs[n];
          gc_total = &gcan_total->gcs[n2];
          if (gcan->labels[n] == Gdiag_no) {
            DiagBreak();
          }
          for (v = r = 0; r < gca->ninputs; r++) {
            gc_total->means[r] += gc->means[r] * gc->ntraining * wt;
            for (c = r; c < gca->ninputs; c++, v++) {
              gc_total->covars[v] += gc->covars[v] * gc->ntraining * wt;
            }
          }
          gc_total->ntraining += wt * gc->ntraining;
          gc_training[n2] += wt * gc->ntraining;
          if (!(gca->flags & GCA_NO_MRF)) {
            for (i = 0; i < GIBBS_NEIGHBORS; i++) {
              for (n3 = 0; n3 < gc->nlabels[i]; n3++) {
                // find the index of this label in the total struct
                for (n4 = 0; n4 < gc_total->nlabels[i]; n4++)
                  if (gc_total->labels[i][n4] == gc->labels[i][n3]) {
                    break;
                  }
                if (n4 < gc_total->nlabels[i]) {
                  gc_total->label_priors[i][n4] += wt * gc->label_priors[i][n3];
                }
              }
            }
          }
        }
      }
    }
  }

  gcan_total->total_training = nint(total_training);

  // now normalize everything
  for (n = 0; n < gcan_total->nlabels; n++) {
    gc_total = &gcan_total->gcs[n];
    gc_total->ntraining = MAX(1, nint(gc_training[n]));
    for (v = r = 0; r < gca->ninputs; r++) {
      gc_total->means[r] /= gc_training[n];
      if (gc_total->means[r] < 5 && gcan_total->labels[n] != Unknown) {
        DiagBreak();
      }
      for (c = r; c < gca->ninputs; c++, v++) {
        if (gc_total->ntraining == 0) {
          DiagBreak();
        }
        else {
          gc_total->covars[v] /= gc_training[n];
        }
      }
    }
  }

  // shouldn't be necessary - normalize them to sum to 1
  for (n = 0; n < gcan_total->nlabels; n++) {
    gc_total = &gcan_total->gcs[n];
    if (!(gca->flags & GCA_NO_MRF)) {
      for (i = 0; i < GIBBS_NEIGHBORS; i++) {
        for (n2 = 0, norm = 0; n2 < gc_total->nlabels[i]; n2++) {
          norm += gc_total->label_priors[i][n2];
        }
        for (n2 = 0; n2 < gc_total->nlabels[i]; n2++) {
          gc_total->label_priors[i][n2] /= norm;

          // to avoid *long* memset
          nbr_labels[i][gcan_total->labels[n]][gc_total->labels[i][n2]] = 0;
          if (FZERO(gc_total->label_priors[i][n2])) {
            DiagBreak();
          }
        }
      }
    }
  }

  if (DIAG_VERBOSE_ON || Gdiag_no == 999) {
    for (i = 0; i < GIBBS_NEIGHBORS; i++)
      for (n = 0; n < MAX_CMA_LABELS; n++)
        for (n2 = 0; n2 < MAX_CMA_LABELS; n2++)
          if (nbr_labels[i][n][n2] > 0) {
            DiagBreak();
          }
  }

  return (gcan_total);
}
static GCA_PRIOR *gcaBuildNbhdGCAP(GCA *gca, int x, int y, int z, float sigma, GCA *gca_smooth)
{
  GCA_PRIOR *gcap, *gcap_total;
  int labels[MAX_CMA_LABELS], n, xi, yi, zi, xk, yk, zk, whalf, nlabels, n2;
  double norm, wt, two_sigma_sq, dsq, total_wt, total_training;
  static int max_label = -1;
  int  wsize ;

  sigma /= gca->prior_spacing ;  // into prior strides
  wsize = 1+2*(int)ceil(3*sigma) ;  // 3 sigmas in each dir and odd
  total_training = 0.0;
  if (max_label < 0) {
    max_label = GCAmaxLabel(gca);
  }
  two_sigma_sq = (2 * sigma * sigma);

  whalf = (int)ceil((float)(wsize - 1) / (2.0));
  memset(labels, 0, sizeof(labels));

  gcap_total = (GCA_PRIOR *)calloc(1, sizeof(GCA_PRIOR));
  if (gcap_total == NULL) {
    ErrorExit(ERROR_NOMEMORY, "gcaBuildNbhdGCAP: could not allocate node");
  }

  // first count # of different labels in nbhd
  for (total_wt = 0.0, xi = -whalf; xi <= whalf; xi++) {
    xk = x + xi;
    if (xk < 0 || xk >= gca->prior_width) {
      continue;
    }
    for (yi = -whalf; yi <= whalf; yi++) {
      yk = y + yi;
      if (yk < 0 || yk >= gca->prior_height) {
        continue;
      }
      for (zi = -whalf; zi <= whalf; zi++) {
        zk = z + zi;
        if (zk < 0 || zk >= gca->prior_depth) {
          continue;
        }
        dsq = xi * xi + yi * yi + zi * zi;
        wt = exp(-dsq / two_sigma_sq);
        if (FZERO(wt)) {
          continue;
        }
        gcap = &gca->priors[xk][yk][zk];
        total_wt += wt;
      }
    }
  }
  for (xi = -whalf; xi <= whalf; xi++) {
    xk = x + xi;
    if (xk < 0 || xk >= gca->prior_width) {
      continue;
    }
    for (yi = -whalf; yi <= whalf; yi++) {
      yk = y + yi;
      if (yk < 0 || yk >= gca->prior_height) {
        continue;
      }
      for (zi = -whalf; zi <= whalf; zi++) {
        zk = z + zi;
        if (zk < 0 || zk >= gca->prior_depth) {
          continue;
        }
        dsq = xi * xi + yi * yi + zi * zi;
        wt = exp(-dsq / two_sigma_sq) / total_wt;
        if (FZERO(wt)) {
          continue;
        }
        gcap = &gca->priors[xk][yk][zk];
        for (n = 0; n < gcap->nlabels; n++) {
          int xn, yn, zn;

          GCApriorToNode(gca, x, y, z, &xn, &yn, &zn);
          if (gca_smooth) {
            int xn2, yn2, zn2;
            xn2 = xn * (gca->node_spacing / gca_smooth->node_spacing);
            yn2 = yn * (gca->node_spacing / gca_smooth->node_spacing);
            zn2 = zn * (gca->node_spacing / gca_smooth->node_spacing);
            if (GCAfindGC(gca_smooth, xn2, yn2, zn2, gcap->labels[n]) == NULL) {
              continue;
            }
          }
	  if (gcap->labels[n] == Gdiag_no && x == Gxp && y == Gyp && z == Gzp)
	    DiagBreak() ;
          if (1 || GCAfindGC(gca, xn, yn, zn, gcap->labels[n]) != NULL) { // DISABLED (BRF)
            labels[gcap->labels[n]]++;
          }
        }
      }
    }
  }
  for (n = nlabels = 0; n <= max_label; n++) {
    if (labels[n] > 0) {
      nlabels++;
    }
  }

  gcap_total->labels = (unsigned short *)calloc(nlabels, sizeof(unsigned short));
  if (gcap_total->labels == NULL) {
    ErrorExit(ERROR_NOMEMORY, "gcaBuildNbhdGCAP: could not allocate node labels");
  }

  gcap_total->priors = (float *)calloc(nlabels, sizeof(float));
  if (gcap_total->labels == NULL) {
    ErrorExit(ERROR_NOMEMORY, "gcaBuildNbhdGCAP: could not allocate node priors");
  }

  gcap_total->nlabels = nlabels;
  for (n = nlabels = 0; n <= max_label; n++) {
    if (labels[n] > 0) {
      gcap_total->labels[nlabels++] = n;
    }
  }

  // accumulate statistics over the nbhd
  norm = 0.0;
  for (xi = -whalf; xi <= whalf; xi++) {
    xk = x + xi;
    if (xk < 0 || xk >= gca->prior_width) {
      continue;
    }
    for (yi = -whalf; yi <= whalf; yi++) {
      yk = y + yi;
      if (yk < 0 || yk >= gca->prior_height) {
        continue;
      }
      for (zi = -whalf; zi <= whalf; zi++) {
        zk = z + zi;
        if (zk < 0 || zk >= gca->prior_depth) {
          continue;
        }
        dsq = xi * xi + yi * yi + zi * zi;
        wt = exp(-dsq / two_sigma_sq) / total_wt;
        if (FZERO(wt)) {
          continue;
        }
        gcap = &gca->priors[xk][yk][zk];
        for (n = 0; n < gcap->nlabels; n++) {
          int xn, yn, zn;

          GCApriorToNode(gca, x, y, z, &xn, &yn, &zn);
          if (gca_smooth) {
            int xn2, yn2, zn2;
            xn2 = xn * (gca->node_spacing / gca_smooth->node_spacing);
            yn2 = yn * (gca->node_spacing / gca_smooth->node_spacing);
            zn2 = zn * (gca->node_spacing / gca_smooth->node_spacing);
            if (GCAfindGC(gca_smooth, xn2, yn2, zn2, gcap->labels[n]) == NULL) {
              continue;
            }
          }
          if (0 && GCAfindGC(gca, xn, yn, zn, gcap->labels[n]) == NULL) {  // DISABLED (BRF)
            continue;
          }
          for (n2 = 0; n2 < gcap_total->nlabels; n2++)
            if (gcap_total->labels[n2] == gcap->labels[n]) {
              break;  // found the corresponding index
            }
          if (n2 >= gcap_total->nlabels) {
            continue;  // should never happen
          }
          gcap_total->priors[n2] += wt * gcap->priors[n];
        }
        total_training += gcap->total_training == 0 ? wt : wt * gcap->total_training;
      }
    }
  }

  gcap_total->total_training = nint(total_training);
  if (gcap_total->total_training == 0) {
    DiagBreak();
  }
  gcap_total->total_training = MAX(1, nint(total_training));

  // now normalize everything
  for (norm = 0.0, n = 0; n < gcap_total->nlabels; n++) {
    norm += gcap_total->priors[n];
  }
  for (n = 0; n < gcap_total->nlabels; n++) {
    gcap_total->priors[n] /= norm;  // to account for sigma weighting
    if (DZERO(gcap_total->priors[n])) {
      DiagBreak();
    }
  }
  return (gcap_total);
}

GC1D *gcanGetGC(GCA_NODE *gcan, int label)
{
  int n;
  GC1D *gc = NULL;

  for (n = 0; n < gcan->nlabels; n++)
    if (gcan->labels[n] == label) {
      gc = &gcan->gcs[n];
      break;
    }
  return (gc);
}

MRI *GCAreclassifyUnlikelyVoxels(GCA *gca,
                                 TRANSFORM *transform,
                                 MRI *mri_inputs,
                                 MRI *mri_aseg,
                                 MRI *mri_aseg_changed,
                                 float mah_dist_thresh,
                                 int wsize,
                                 float sigma)
{
  int x, y, z, max_label, nchanged, iter = 0;
  GCA_NODE *gcan_total;
  // GCA_NODE *gcan;
  GCA_PRIOR *gcap_total, *gcap;
  int xp, yp, zp, xn, yn, zn, width, height, depth, label, n;
  // int max_n;
  double p, max_p, d_old;
  // double d_new;
  float vals[MAX_GCA_INPUTS];
  GC1D *gc, *max_gc;
  MRI *mri_changed;

  mri_changed = MRIclone(mri_aseg, NULL);
  MRIsetValues(mri_changed, 1.0);
  printf("relabeling unlikely (>%2.0f sigmas from mean) voxels in %d mm window, with sigma = %2.0f\n",
         mah_dist_thresh,
         wsize,
         sigma);
  mri_aseg_changed = MRIcopy(mri_aseg, mri_aseg_changed);

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;

  do {
    nchanged = 0;
    for (x = 0; x < width; x++) {
      for (y = 0; y < height; y++) {
        for (z = 0; z < depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }

          if (nint(MRIgetVoxVal(mri_changed, x, y, z, 0)) == 0) {
            continue;
          }
          if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
            load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

            // gcan = &gca->nodes[xn][yn][zn];
            gcap = getGCAP(gca, mri_inputs, transform, x, y, z);
            if (gcap == NULL) {
              continue;
            }
            label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0);
            gc = GCAfindGC(gca, xn, yn, zn, label);

            if (gc == NULL) {
              continue;
            }
            d_old = sqrt(GCAmahDist(gc, vals, gca->ninputs));
            if (d_old > mah_dist_thresh) {
              GCAsourceVoxelToPrior(gca, mri_inputs, transform, x, y, z, &xp, &yp, &zp);
              if (label != Unknown && !IS_CORTEX(label)) {
                DiagBreak();
              }
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
                DiagBreak();
              }
              gcan_total = gcaBuildNbhdGCAN(gca, xn, yn, zn, sigma);
              gcap_total = gcaBuildNbhdGCAP(gca, xp, yp, zp, sigma, NULL);
              max_label = 0;
              // max_n = -1;
              max_p = 2 * GIBBS_NEIGHBORS * BIG_AND_NEGATIVE;
              max_gc = NULL;

              // going through gcap labels
              for (n = 0; n < gcap_total->nlabels; n++) {
                gc = gcanGetGC(gcan_total, gcap_total->labels[n]);
                if (gc == NULL) {
                  continue;
                }
                p = gcaComputeLogDensity(gc, vals, gca->ninputs, gcap_total->priors[n], gcap_total->labels[n]);
                p = gcaGibbsLogPosterior(
                    gca, mri_aseg_changed, vals, gcap_total->labels[n], x, y, z, gcap_total, gcan_total, transform);
                if (p > max_p) {
                  max_p = p;
                  max_label = gcap_total->labels[n];
                  // max_n = n;
                  max_gc = gc;
                }
              }

              MRIsetVoxVal(mri_changed, x, y, z, 0, 0);  // will change to 1 if changed
              if (max_gc) {
                // d_new = sqrt(GCAmahDist(max_gc, vals, gca->ninputs));
                // if (d_new < mah_dist_thresh/2)
                {
                  if (max_label != label) {
                    nchanged++;
                    MRIsetVoxVal(mri_changed, x, y, z, 0, 1);
                  }

                  if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                    printf(
                        "GCAreclassifyUnlikelyVoxels: "
                        "(%d, %d, %d) %s (%d) --> %s (%d)\n",
                        x,
                        y,
                        z,
                        cma_label_to_name(label),
                        label,
                        cma_label_to_name(max_label),
                        max_label);
                  MRIsetVoxVal(mri_aseg_changed, x, y, z, 0, max_label);
                }
              }
              GCANfree(gcan_total, gca->ninputs);
              free(gcan_total);
              GCAPfree(gcap_total);
              free(gcap_total);
            }
          }
        }
      }
    }
    MRIdilate(mri_changed, mri_changed);
    printf("iter %d: nchanged = %d\n", ++iter, nchanged);
  } while (nchanged > 50 && iter < 10);

  MRIfree(&mri_changed);
  return (mri_aseg_changed);
}
static double gcaGibbsLogPosterior(GCA *gca,
                                   MRI *mri_labels,
                                   float *vals,
                                   int label,
                                   int x,
                                   int y,
                                   int z,
                                   GCA_PRIOR *gcap,
                                   GCA_NODE *gcan,
                                   TRANSFORM *transform)
{
  double log_posterior, nbr_prior;
  int xnbr, ynbr, znbr, nbr_label, i, j, n;
  GC1D *gc = 0;
#if INTERP_PRIOR
  double prior;
#endif
  // signify error
  log_posterior = 0.;

  /////////////////////////////////////////////////////////////////
  for (n = 0; n < gcan->nlabels; n++) {
    if (gcan->labels[n] == label) {
      break;
    }
  }
  // could not find the label, then
  if (n >= gcan->nlabels) {
    // if (gcan->total_training > 0)
    // return(log(0.01f/((float)gcan->total_training*GIBBS_NEIGHBORS))) ;
    /* 10*GIBBS_NEIGHBORS*BIG_AND_NEGATIVE*/
    // else
    return (10 * BIG_AND_NEGATIVE);
    // return(log(VERY_UNLIKELY)) ;
  }

  gc = &gcan->gcs[n];

  /* compute 1-d Mahalanobis distance */
  log_posterior = GCAcomputeConditionalLogDensity(gc, vals, gca->ninputs, gcan->labels[n]);
  if (check_finite("GCAvoxelGibbsLogPosterior: conditional log density", log_posterior) == 0) {
    DiagBreak();
  }

  nbr_prior = 0.0;
  for (i = 0; i < GIBBS_NEIGHBORS; i++) {
    xnbr = mri_labels->xi[x + xnbr_offset[i]];
    ynbr = mri_labels->yi[y + ynbr_offset[i]];
    znbr = mri_labels->zi[z + znbr_offset[i]];
    nbr_label = nint(MRIgetVoxVal(mri_labels, xnbr, ynbr, znbr, 0));
    for (j = 0; j < gc->nlabels[i]; j++) {
      if (nbr_label == gc->labels[i][j]) {
        break;
      }
    }
    if (j < gc->nlabels[i]) {
      if (!FZERO(gc->label_priors[i][j])) {
        nbr_prior += log(gc->label_priors[i][j]);
      }
      else {
        nbr_prior += log(0.1f / (float)gcan->total_training);
      }
      /*BIG_AND_NEGATIVE */
      check_finite("gcaGibbsLogPosterior: label_priors", nbr_prior);
    }
    else /* never occurred - make it unlikely */
    {
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
        DiagBreak();
      }
      nbr_prior += log(0.1f / (float)gcan->total_training);
    }
  }
// added to the previous value
#if INTERP_PRIOR
  prior = gcaComputePrior(gca, mri_labels, transform, x, y, z, label);
  log_posterior += (nbr_prior + log(prior));
#else
  log_posterior += (nbr_prior + log(getPrior(gcap, label)));
#endif
  if (check_finite("GCAvoxelGibbsLogPosterior: final", log_posterior) == 0) {
    DiagBreak();
  }

  return (log_posterior);
}

MRI *GCAcomputeOptimalScale(GCA *gca,
                            TRANSFORM *transform,
                            MRI *mri_inputs,
                            MRI *mri_aseg,
                            MRI *mri_sigma,
                            int wsize,
                            double min_sigma,
                            double max_sigma,
                            double delta_sigma)
{
  int x, y, z;
  GCA_NODE *gcan_total;
  GCA_PRIOR *gcap_total;
  int xp, yp, zp, xn, yn, zn, width, height, depth, n;
  double p, best_p, best_sigma, total_p, scale, first_sigma, last_sigma;
  float vals[MAX_GCA_INPUTS];
  double sigma;

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;

  if (mri_sigma == NULL) {
    mri_sigma = MRIcloneDifferentType(mri_aseg, MRI_FLOAT);
  }

  for (x = 0; x < width; x++) {
    for (y = 0; y < height; y++) {
      for (z = 0; z < depth; z++) {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
          DiagBreak();
        }

        if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
          load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

          GCAsourceVoxelToPrior(gca, mri_inputs, transform, x, y, z, &xp, &yp, &zp);
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }
          best_p = -1;
          best_sigma = -1;
          first_sigma = min_sigma;
          last_sigma = max_sigma;
          for (scale = max_sigma / 2; scale >= 1; scale /= 2) {
            for (sigma = first_sigma; sigma <= last_sigma; sigma += scale * delta_sigma) {
              gcan_total = gcaBuildNbhdGCAN(gca, xn, yn, zn, sigma);
              gcap_total = gcaBuildNbhdGCAP(gca, xp, yp, zp, sigma, NULL);

              // marginalize over labels
              total_p = 0;
              for (n = 0; n < gcap_total->nlabels; n++) {
                p = gcaGibbsLogPosterior(
                    gca, mri_aseg, vals, gcap_total->labels[n], x, y, z, gcap_total, gcan_total, transform);
                total_p += exp(p);
              }

              if (total_p > best_p) {
                best_p = total_p;
                best_sigma = sigma;
              }

              GCANfree(gcan_total, gca->ninputs);
              free(gcan_total);
              GCAPfree(gcap_total);
              free(gcap_total);
            }
            first_sigma = MAX(min_sigma, best_sigma - scale * delta_sigma);
            last_sigma = MIN(max_sigma, best_sigma + scale * delta_sigma);
          }
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            printf("setting optimal scale at (%d, %d, %d) to %2.3f\n", x, y, z, best_sigma);
          }

          MRIsetVoxVal(mri_sigma, x, y, z, 0, best_sigma);
        }
      }
    }
  }

  return (mri_sigma);
}
MRI *GCAreclassifyVoxelsAtOptimalScale(
    GCA *gca, TRANSFORM *transform, MRI *mri_inputs, MRI *mri_aseg, MRI *mri_aseg_changed, MRI *mri_sigma, int wsize)
{
  int x, y, z, max_label, nchanged, iter = 0;
  GCA_NODE *gcan_total;
  GCA_PRIOR *gcap_total;
  int xp, yp, zp, xn, yn, zn, width, height, depth, label, n;
  // int max_n;
  double p, max_p, sigma;
  float vals[MAX_GCA_INPUTS];
  GC1D *gc, *max_gc;
  MRI *mri_changed;

  mri_changed = MRIclone(mri_aseg, NULL);
  MRIsetValues(mri_changed, 1.0);
  printf("relabeling voxels using optimal scale estimates\n");
  mri_aseg_changed = MRIcopy(mri_aseg, mri_aseg_changed);

  width = mri_inputs->width;
  height = mri_inputs->height;
  depth = mri_inputs->depth;

  do {
    nchanged = 0;
    for (x = 0; x < width; x++) {
      for (y = 0; y < height; y++) {
        for (z = 0; z < depth; z++) {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
            DiagBreak();
          }

          if (nint(MRIgetVoxVal(mri_changed, x, y, z, 0)) == 0) {
            continue;
          }
          if (!GCAsourceVoxelToNode(gca, mri_inputs, transform, x, y, z, &xn, &yn, &zn)) {
            load_vals(mri_inputs, x, y, z, vals, gca->ninputs);

            GCAsourceVoxelToPrior(gca, mri_inputs, transform, x, y, z, &xp, &yp, &zp);
            label = (int)MRIgetVoxVal(mri_aseg_changed, x, y, z, 0);
            if (label != Unknown && !IS_CORTEX(label)) {
              DiagBreak();
            }
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z) {
              DiagBreak();
            }
            sigma = MRIgetVoxVal(mri_sigma, x, y, z, 0);
            gcan_total = gcaBuildNbhdGCAN(gca, xn, yn, zn, sigma);
            gcap_total = gcaBuildNbhdGCAP(gca, xp, yp, zp, sigma, NULL);
            max_label = 0;
            // max_n = -1;
            max_p = 2 * GIBBS_NEIGHBORS * BIG_AND_NEGATIVE;
            max_gc = NULL;

            // going through gcap labels
            for (n = 0; n < gcap_total->nlabels; n++) {
              gc = gcanGetGC(gcan_total, gcap_total->labels[n]);
              if (gc == NULL) {
                continue;
              }
              p = gcaComputeLogDensity(gc, vals, gca->ninputs, gcap_total->priors[n], gcap_total->labels[n]);
              p = gcaGibbsLogPosterior(
                  gca, mri_aseg_changed, vals, gcap_total->labels[n], x, y, z, gcap_total, gcan_total, transform);
              if (p > max_p) {
                max_p = p;
                max_label = gcap_total->labels[n];
                // max_n = n;
                max_gc = gc;
              }
            }

            MRIsetVoxVal(mri_changed, x, y, z, 0, 0);  // will change to 1 if changed
            if (max_gc && (max_label != label)) {
              nchanged++;
              MRIsetVoxVal(mri_changed, x, y, z, 0, 1);

              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                printf(
                    "GCAreclassifyVoxelsAtOptimalScale: "
                    "(%d, %d, %d, s=%2.1f) %s (%d) --> %s (%d)\n",
                    x,
                    y,
                    z,
                    sigma,
                    cma_label_to_name(label),
                    label,
                    cma_label_to_name(max_label),
                    max_label);
              MRIsetVoxVal(mri_aseg_changed, x, y, z, 0, max_label);
            }
            GCANfree(gcan_total, gca->ninputs);
            free(gcan_total);
            GCAPfree(gcap_total);
            free(gcap_total);
          }
        }
      }
    }
    MRIdilate(mri_changed, mri_changed);
    printf("iter %d: nchanged = %d\n", ++iter, nchanged);
  } while (nchanged > 50 && iter < 10);

  MRIfree(&mri_changed);
  return (mri_aseg_changed);
}
int GCAcheck(GCA *gca)
{
  int error, xp, yp, zp, n, xn, yn, zn, label;
  GCA_PRIOR *gcap;
  GC1D *gc;

  error = NO_ERROR;
  for (xp = 0; xp < gca->prior_width; xp++)
    for (yp = 0; yp < gca->prior_height; yp++)
      for (zp = 0; zp < gca->prior_depth; zp++) {
        if (xp == Gxp && yp == Gyp && zp == Gzp) {
          DiagBreak();
        }
        gcap = &gca->priors[xp][yp][zp];
        GCApriorToNode(gca, xp, yp, zp, &xn, &yn, &zn);
        for (n = 0; n < gcap->nlabels; n++) {
          label = gcap->labels[n];
          gc = GCAfindGC(gca, xn, yn, zn, label);
          if (gc == NULL)
            printf("gcap(%d, %d, %d), labels[%d] = %s (%d) --> node (%d, %d, %d), no gc!\n",
                   xp,
                   yp,
                   zp,
                   n,
                   cma_label_to_name(label),
                   label,
                   xn,
                   yn,
                   zn);
          error = ERROR_BADPARM;
        }
      }

  return (error);
}
GCA *GCAcopy(GCA *gca_src, GCA *gca_dst)
{
  if (gca_dst == NULL)
    gca_dst = GCAalloc(gca_src->ninputs,
                       gca_src->prior_spacing,
                       gca_src->node_spacing,
                       gca_src->node_width * gca_src->node_spacing,
                       gca_src->node_height * gca_src->node_spacing,
                       gca_src->node_depth * gca_src->node_spacing,
                       gca_src->flags);

  gca_dst->ninputs = gca_src->ninputs;
  gca_dst->type = gca_src->type;

  gca_dst->xsize = gca_src->xsize;
  gca_dst->ysize = gca_src->ysize;
  gca_dst->zsize = gca_src->zsize;
  gca_dst->total_training = gca_src->total_training;
  GCAcopyDCToGCA(gca_src, gca_dst);

  memmove(gca_dst->tissue_parms, gca_src->tissue_parms, sizeof(gca_src->tissue_parms));
  memmove(gca_dst->TRs, gca_src->TRs, sizeof(gca_src->TRs));
  memmove(gca_dst->FAs, gca_src->FAs, sizeof(gca_src->FAs));
  memmove(gca_dst->TEs, gca_src->TEs, sizeof(gca_src->TEs));
  return (gca_dst);
}
GCA *GCAsmooth(GCA *gca, double sigma)
{
  GCA_NODE *gcan_total, *gcan;
  GCA_PRIOR *gcap_total, *gcap;
  int xp, yp, zp, xn, yn, zn, n;
  GCA *gca_smooth;

  gca_smooth = GCAcopy(gca, NULL);

  for (xn = 0; xn < gca->node_width; xn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (zn = 0; zn < gca->node_depth; zn++) {
        if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
          DiagBreak();
        }

        gcan = &gca_smooth->nodes[xn][yn][zn];
        gcan_total = gcaBuildNbhdGCAN(gca, xn, yn, zn, sigma);
        if (gcan_total->nlabels == 0) {
          gcan_total->nlabels = 1;
          gcan_total->labels[0] = 0;
          gcan_total->gcs = alloc_gcs(gcan_total->nlabels, gca->flags, gca->ninputs);
          gcan_total->gcs[0].covars[0] = 25;
        }
        if (gcan->max_labels < gcan_total->nlabels) {
          free(gcan->labels);
          gcan->labels = (unsigned short *)calloc(gcan_total->nlabels, sizeof(unsigned short));
          if (gcan->labels == NULL)
            ErrorExit(ERROR_NOMEMORY, "GCAsmooth(%2.2f) couldn't allocate %d label node", sigma, gcan_total->nlabels);
        }
        free_gcs(gcan->gcs, gcan->nlabels, gca->ninputs);
        gcan->gcs = alloc_gcs(gcan_total->nlabels, gca->flags, gca->ninputs);
        copy_gcs(gcan_total->nlabels, gcan_total->gcs, gcan->gcs, gca->ninputs);
        gcan->nlabels = gcan_total->nlabels;
        for (n = 0; n < gcan_total->nlabels; n++) {
          if (fabs(gcan_total->gcs[n].covars[0] - 25) > 5) {
            DiagBreak();
          }
          gcan->labels[n] = gcan_total->labels[n];
          if (!std::isfinite(gcan->gcs[n].means[0])) {
            DiagBreak();
          }
          if (gcan->gcs[n].covars == NULL) {
            DiagBreak();
          }
          if (fabs(gcan->gcs[n].covars[0] - 25) > 0.1) {
            DiagBreak();
          }
          gcan->gcs[n].covars[0] = 25;
        }
        gcan->total_training = gcan_total->total_training;
        GCANfree(gcan_total, gca->ninputs);
        free(gcan_total);
      }
    }
  }

  for (xp = 0; xp < gca->prior_width; xp++) {
    for (yp = 0; yp < gca->prior_height; yp++) {
      for (zp = 0; zp < gca->prior_depth; zp++) {
        if (xp == Ggca_x && yp == Ggca_y && zp == Ggca_z) {
          DiagBreak();
        }
        if (xp == Gx && yp == Gy && zp == Gz) {
          DiagBreak();
        }
        gcap_total = gcaBuildNbhdGCAP(gca, xp, yp, zp, sigma, gca_smooth);
        gcap = &gca_smooth->priors[xp][yp][zp];
        gcap->nlabels = gcap_total->nlabels;
        if (gcap_total->nlabels > gcap->max_labels) {
          free(gcap->labels);
          free(gcap->priors);
          gcap->labels = (unsigned short *)calloc(gcap->nlabels, sizeof(unsigned short));
          if (!gcap->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAsmooth(%2.2f): could not "
                      "allocate %d "
                      "labels @ (%d,%d,%d)",
                      sigma,
                      gcap->nlabels,
                      xp,
                      yp,
                      zp);
          gcap->priors = (float *)calloc(gcap->nlabels, sizeof(float));
          if (!gcap->priors)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAsmooth(%s): could "
                      "not allocate %d "
                      "priors @ (%d,%d,%d)",
                      sigma,
                      gcap->nlabels,
                      xp,
                      yp,
                      zp);
        }

        for (n = 0; n < gcap_total->nlabels; n++) {
          GC1D *gc;
          gcap->labels[n] = gcap_total->labels[n];
          gcap->priors[n] = gcap_total->priors[n];
          if (!std::isfinite(gcap->priors[n])) {
            DiagBreak();
          }
          GCApriorToNode(gca_smooth, xp, yp, zp, &xn, &yn, &zn);
          gcan = &gca_smooth->nodes[xn][yn][zn];
          gc = gcanGetGC(gcan, gcap->labels[n]);
          if (gc == NULL) {
            DiagBreak();
          }
        }
        if (gcap->nlabels == 0) {
          gcap->labels[0] = 0;
          gcap->priors[0] = 1.0;
          gcap->nlabels = 1;
        }

        GCAPfree(gcap_total);
        free(gcap_total);
      }
    }
  }
  gcaCheck(gca_smooth) ;
  return (gca_smooth);
}
GCA *GCAnodeDownsample2(GCA *gca)
{
  GCA_NODE *gcan_total, *gcan;
  GCA_PRIOR *gcap_total, *gcap;
  int xp, yp, zp, xn, yn, zn, wsize, n, xn2, yn2, zn2;
  GCA *gca_smooth;
  float sigma;

  sigma = gca->node_spacing / 2;
  wsize = ceil(4 * sigma) + 1;
  gca_smooth =
      GCAalloc(gca->ninputs, gca->prior_spacing, gca->node_spacing * 2, gca->width, gca->height, gca->depth, 0);
  GCAcopy(gca, gca_smooth);

  for (xn = 0; xn < gca->node_width; xn++) {
    for (yn = 0; yn < gca->node_height; yn++) {
      for (zn = 0; zn < gca->node_depth; zn++) {
        if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) {
          DiagBreak();
        }

        xn2 = xn / 2;
        yn2 = yn / 2;
        zn2 = zn / 2;
        gcan_total = gcaBuildNbhdGCAN(gca, xn, yn, zn, sigma);
        gcan = &gca_smooth->nodes[xn2][yn2][zn2];
        if (gcan_total->nlabels == 0) {
          gcan_total->nlabels = 1;
          gcan_total->labels[0] = 0;
          gcan_total->gcs = alloc_gcs(gcan_total->nlabels, gca->flags, gca->ninputs);
          gcan_total->gcs[0].covars[0] = 25;
        }
        if (gcan->max_labels < gcan_total->nlabels) {
          free(gcan->labels);
          gcan->labels = (unsigned short *)calloc(gcan_total->nlabels, sizeof(unsigned short));
          if (gcan->labels == NULL)
            ErrorExit(ERROR_NOMEMORY, "GCAsmooth(%2.2f) couldn't allocate %d label node", sigma, gcan_total->nlabels);
        }
        free_gcs(gcan->gcs, gcan->nlabels, gca->ninputs);
        gcan->gcs = alloc_gcs(gcan_total->nlabels, gca->flags, gca->ninputs);
        copy_gcs(gcan_total->nlabels, gcan_total->gcs, gcan->gcs, gca->ninputs);
        gcan->nlabels = gcan_total->nlabels;
        for (n = 0; n < gcan_total->nlabels; n++) {
          if (fabs(gcan_total->gcs[n].covars[0] - 25) > 5) {
            DiagBreak();
          }
          gcan->labels[n] = gcan_total->labels[n];
          if (!std::isfinite(gcan->gcs[n].means[0])) {
            DiagBreak();
          }
          if (gcan->gcs[n].covars == NULL) {
            DiagBreak();
          }
          if (fabs(gcan->gcs[n].covars[0] - 25) > 0.1) {
            DiagBreak();
          }
          gcan->gcs[n].covars[0] = 25;
        }
        gcan->total_training = gcan_total->total_training;
        GCANfree(gcan_total, gca->ninputs);
        free(gcan_total);
      }
    }
  }

  for (xp = 0; xp < gca->prior_width; xp++) {
    for (yp = 0; yp < gca->prior_height; yp++) {
      for (zp = 0; zp < gca->prior_depth; zp++) {
        if (xp == Ggca_x && yp == Ggca_y && zp == Ggca_z) {
          DiagBreak();
        }
        if (xp == Gx && yp == Gy && zp == Gz) {
          DiagBreak();
        }
        gcap_total = gcaBuildNbhdGCAP(gca, xp, yp, zp, 0, gca_smooth);
        gcap = &gca_smooth->priors[xp][yp][zp];
        gcap->nlabels = gcap_total->nlabels;
        if (gcap_total->nlabels > gcap->max_labels) {
          free(gcap->labels);
          free(gcap->priors);
          gcap->labels = (unsigned short *)calloc(gcap->nlabels, sizeof(unsigned short));
          if (!gcap->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAsmooth(%2.2f): could not "
                      "allocate %d "
                      "labels @ (%d,%d,%d)",
                      sigma,
                      gcap->nlabels,
                      xp,
                      yp,
                      zp);
          gcap->priors = (float *)calloc(gcap->nlabels, sizeof(float));
          if (!gcap->priors)
            ErrorExit(ERROR_NOMEMORY,
                      "GCAsmooth(%s): could "
                      "not allocate %d "
                      "priors @ (%d,%d,%d)",
                      sigma,
                      gcap->nlabels,
                      xp,
                      yp,
                      zp);
        }

        for (n = 0; n < gcap_total->nlabels; n++) {
          GC1D *gc;
          gcap->labels[n] = gcap_total->labels[n];
          gcap->priors[n] = gcap_total->priors[n];
          if (!std::isfinite(gcap->priors[n])) {
            DiagBreak();
          }
          GCApriorToNode(gca_smooth, xp, yp, zp, &xn, &yn, &zn);
          gcan = &gca_smooth->nodes[xn][yn][zn];
          gc = gcanGetGC(gcan, gcap->labels[n]);
          if (gc == NULL) {
            DiagBreak();
          }
        }
        if (gcap->nlabels == 0) {
          gcap->labels[0] = 0;
          gcap->priors[0] = 1.0;
          gcap->nlabels = 1;
        }

        GCAPfree(gcap_total);
        free(gcap_total);
      }
    }
  }
  return (gca_smooth);
}
int GCAsourceVoxelToPriorReal(
    GCA *gca, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv, double *pxp, double *pyp, double *pzp)
{
  float xt = 0, yt = 0, zt = 0;
  double xrt, yrt, zrt;
  // int retval;

  LTA *lta;
  if (transform->type != MORPH_3D_TYPE) {
    if (transform->type == LINEAR_VOX_TO_VOX) {
      lta = (LTA *)transform->xform;
      // transform point to talairach volume point
      TransformWithMatrix(lta->xforms[0].m_L, xv, yv, zv, &xrt, &yrt, &zrt);
      xt = xrt;
      yt = yrt;
      zt = zrt;
      // TransformSample(transform, xv, yv, zv, &xt, &yt, &zt) ;
    }
    else
      ErrorExit(ERROR_BADPARM, "GCAsourceVoxelToPrior: needs vox-to-vox transform");
  }
  else  // morph 3d type can go directly from source to template
  {
    // printf("GCAsourcevoxelToPriorReal: (xv,yv,zv) = (%d,%d,%d).\n", xv, yv, zv);
    TransformSample(transform, xv, yv, zv, &xt, &yt, &zt);
    // printf("GCAsourcevoxelToPriorReal: (xt,yt,zt) = (%f,%f,%f).\n", xt, yt, zt); //LZ
  }
  // get the position in gca from talairach volume
  // printf("GCAsourcevoxelToPriorReal: printing mri_tal\n");
  if (!gca->mri_tal__) {
    gca->mri_tal__ = MRIallocHeader(gca->width, gca->height, gca->depth, MRI_UCHAR, 1);
    gca->mri_tal__->xsize = gca->xsize;
    gca->mri_tal__->ysize = gca->ysize;
    gca->mri_tal__->zsize = gca->zsize;

    GCAcopyDCToMRI(gca, gca->mri_tal__);

    gca->tal_i_to_r__ = extract_i_to_r(gca->mri_tal__);
    gca->tal_r_to_i__ = extract_r_to_i(gca->mri_tal__);
  }

  // LZ : MRIprintStats(gca->mri_tal__, stdout);
  GCAvoxelToPriorReal(gca, gca->mri_tal__, xt, yt, zt, pxp, pyp, pzp);
  if (*pxp < 0 || *pyp < 0 || *pzp < 0 || *pxp >= gca->prior_width || *pyp >= gca->prior_height ||
      *pzp >= gca->prior_depth) {
    // retval = (ERROR_BADPARM);
  }
  if (*pxp < 0) {
    *pxp = 0;
  }
  if (*pyp < 0) {
    *pyp = 0;
  }
  if (*pzp < 0) {
    *pzp = 0;
  }
  if (*pxp >= gca->prior_width) {
    *pxp = gca->prior_width - 1;
  }
  if (*pyp >= gca->prior_height) {
    *pyp = gca->prior_height - 1;
  }
  if (*pzp >= gca->prior_depth) {
    *pzp = gca->prior_depth - 1;
  }
  return (NO_ERROR);
}

int GCAinsertLabels(GCA *gca,
                    MRI *mri,
                    TRANSFORM *transform,
                    int ninsertions,
                    int *insert_labels,
                    int *insert_intensities,
                    int insert_coords[MAX_INSERTIONS][3],
                    int *insert_whalf)
{
  int xn, yn, zn, xp, yp, zp, i, l, label, x, y, z, found, n, whalf;
  int xk, yk, zk, xi, yi, zi;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;

  for (l = 0; l < ninsertions; l++) {
    whalf = insert_whalf[l];
    label = insert_labels[l];
    x = insert_coords[l][0];
    y = insert_coords[l][1];
    z = insert_coords[l][2];
    printf("inserting label %d (%s) at (%d, %d, %d) with intensity = %d\n",
           label,
           cma_label_to_name(label),
           x,
           y,
           z,
           insert_intensities[l]);

    for (xk = -whalf; xk <= whalf; xk++)
      for (yk = -whalf; yk <= whalf; yk++)
        for (zk = -whalf; zk <= whalf; zk++) {
          xi = mri->xi[x + xk];
          yi = mri->yi[y + yk];
          zi = mri->zi[z + zk];

          found = 0;
          if (!GCAsourceVoxelToNode(gca, mri, transform, xi, yi, zi, &xn, &yn, &zn)) {
            gcan = &gca->nodes[xn][yn][zn];
            for (i = 0; i < gcan->nlabels; i++)
              if (gcan->labels[i] == label) {
                found = 1;
                break;
              }
            if (found == 0)  // insert it
            {
              GC1D *gc, *gcs;

              gcs = alloc_gcs(gcan->nlabels + 1, gca->flags, gca->ninputs);
              gc = &gcs[gcan->nlabels];
              printf("inserting label %s at node (%d, %d, %d)\n", cma_label_to_name(label), xn, yn, zn);
              gc->means[0] = insert_intensities[l];
              gc->covars[0] = MIN_VAR;
              if ((gca->flags & GCA_NO_MRF) == 0) {
                for (i = 0; i < GIBBS_NEIGHBORS; i++) {
                  gc->nlabels[i] = gcan->nlabels + 1;  // assume any of existing labels can occur plus this one
                  gc->labels[i] = (unsigned short *)calloc(gc->nlabels[i], sizeof(unsigned short));
                  if (!gc->labels) {
                    ErrorExit(ERROR_NOMEMORY, "GCAinsertLabels: couldn't expand labels to %d", gc->nlabels[i]);
                  }
                  gc->label_priors[i] = (float *)calloc(gc->nlabels[i], sizeof(float));
                  if (!gc->label_priors[i]) {
                    ErrorExit(ERROR_NOMEMORY, "GCAinsertLabel: couldn't expand gcs to %d", gc->nlabels);
                  }
                  for (n = 0; n < gcan->nlabels; n++) {
                    gc->labels[i][n] = gcan->labels[n];
                    gc->label_priors[i][n] = 1.0 / gc->nlabels[i];
                  }
                  gc->labels[i][n] = label;
                  gc->label_priors[i][n] = 1.0 / gc->nlabels[i];
                }
              }
              copy_gcs(gcan->nlabels, gcan->gcs, gcs, gca->ninputs);
              free_gcs(gcan->gcs, gcan->nlabels, gca->ninputs);
              gc->ntraining = gcan->total_training;  // arbitrary
              gcan->total_training *= 2;
              gcan->gcs = gcs;
              gcan->labels[gcan->nlabels++] = label;
            }
          }
          if (!GCAsourceVoxelToPrior(gca, mri, transform, xi, yi, zi, &xp, &yp, &zp)) {
            found = 0;
            gcap = &gca->priors[xp][yp][zp];
            for (i = 0; i < gcap->nlabels; i++)
              if (gcap->labels[i] == label) {
                found = 1;
                break;
              }
            if (found == 0) {
              gcap = &gca->priors[xp][yp][zp];
              printf("inserting label %s at prior (%d, %d, %d)\n", cma_label_to_name(label), xp, yp, zp);
              gcap->nlabels = 1;
              gcap->priors[0] = 1.0;
              gcap->labels[0] = label;
            }
          }
        }
  }
  return (NO_ERROR);
}

float GCAcomputeNumberOfGoodFittingSamples(
    GCA *gca, GCA_SAMPLE *gcas, MRI *mri_inputs, TRANSFORM *transform, int nsamples)

{
  int x, y, z, i, xp, yp, zp;
  // int width, height, depth;
  float vals[MAX_GCA_INPUTS];
  double total_log_p, log_p, dist;
  int countOutside = 0;
  double outside_log_p = 0.;
  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  // width = mri_inputs->width;
  // height = mri_inputs->height;
  // depth = mri_inputs->depth;
  // store inverse transformation .. forward:input->gca template,
  // inv: gca template->input
  TransformInvert(transform, mri_inputs);

  // go through all sample points
  for (total_log_p = 0.0, i = 0; i < nsamples; i++) {
    /////////////////// diag code /////////////////////////////
    if (i == Gdiag_no) {
      DiagBreak();
    }
    if (Gdiag_no == gcas[i].label) {
      DiagBreak();
    }
    if (i == Gdiag_no || (gcas[i].xp == Gxp && gcas[i].yp == Gyp && gcas[i].zp == Gzp)) {
      DiagBreak();
    }
    ///////////////////////////////////////////////////////////

    // get prior coordinates
    xp = gcas[i].xp;
    yp = gcas[i].yp;
    zp = gcas[i].zp;
    // if it is inside the source voxel
    if (!GCApriorToSourceVoxel(gca, mri_inputs, transform, xp, yp, zp, &x, &y, &z)) {
      if (x == Gx && y == Gy && z == Gz) {
        DiagBreak();
      }

      // (x,y,z) is the source voxel position
      gcas[i].x = x;
      gcas[i].y = y;
      gcas[i].z = z;
      // get values from all inputs
      load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
      log_p = gcaComputeSampleLogDensity(&gcas[i], vals, gca->ninputs);
      if (FZERO(vals[0]) && gcas[i].label == Gdiag_no) {
        if (fabs(log_p) < 5) {
          DiagBreak();
        }
        DiagBreak();
      }

      if (!FZERO(vals[0])) {
        DiagBreak();
      }
      if (gcas[i].label != Unknown) {
        DiagBreak();
      }
      if (i == Gdiag_no) {
        DiagBreak();
      }
      dist = sqrt(GCAsampleMahDist(&gcas[i], vals, gca->ninputs));
      if (dist > 1) {
        total_log_p -= 1.0;
      }
      if (dist > 2) {
        total_log_p -= 1.0;
      }
      if (dist > 3) {
        total_log_p -= 1.0;
      }
      if (dist > 4) {
        total_log_p -= 1.0;
      }
      gcas[i].log_p = log_p;

      if (!check_finite("2", total_log_p)) {
        fprintf(stdout, "total log p not finite at (%d, %d, %d)\n", x, y, z);
        DiagBreak();
      }
    }
    else  // outside the volume
    {
      log_p = -1000000;  // BIG_AND_NEGATIVE;
      // log(VERY_UNLIKELY); // BIG_AND_NEGATIVE;
      total_log_p--;
      gcas[i].log_p = log_p;
      outside_log_p += log_p;
      countOutside++;
    }
  }

#ifndef __OPTIMIZE__
#endif
  fflush(stdout);

  return ((float)total_log_p);
}

int GCAreadLabelIntensities(char *fname, float *label_scales, float *label_offsets)
{
  FILE *fp;
  int l;
  char *cp, line[STRLEN];
  float scale, offset;

  for (l = 0; l < MAX_CMA_LABELS; l++) {
    label_scales[l] = 1.0;
    label_offsets[l] = 0;
  }
  printf("reading intensity scales from %s\n", fname);
  fp = fopen(fname, "r");
  if (fp == NULL) {
    ErrorExit(ERROR_NOFILE, "%s: could not open intensity tracking file %s", Progname, fname);
  }

  cp = fgetl(line, STRLEN - 1, fp);
  while (cp) {
    sscanf(cp, "%d %*s %f %f %*f\n", &l, &scale, &offset);
    label_scales[l] = scale;
    label_offsets[l] = offset;
    printf("label %s (%d): %2.2f + %2.1f\n", cma_label_to_name(l), l, scale, offset);
    cp = fgetl(line, STRLEN - 1, fp);
  }

  fclose(fp);
  return (NO_ERROR);
}

int GCAremoveHemi(GCA *gca, int lh)
{
  int x, y, z, n, r;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          if ((lh && (IS_LH_CLASS(gcan->labels[n]))) || (!lh && (IS_RH_CLASS(gcan->labels[n])))) {
            gcan->labels[n] = Unknown;
            for (r = 0; r < gca->ninputs; r++) {
              gcan->gcs[n].means[r] = 0;
            }
          }
        }
      }
    }
  }

  for (x = 0; x < gca->prior_width; x++) {
    for (y = 0; y < gca->prior_height; y++) {
      for (z = 0; z < gca->prior_depth; z++) {
        int max_label;
        double max_prior;

        gcap = &gca->priors[x][y][z];
        max_prior = 0;
        max_label = 0;
        for (n = 0; n < gcap->nlabels; n++) {
          if (gcap->priors[n] > max_prior) {
            max_prior = gcap->priors[n];
            max_label = gcap->labels[n];
          }
        }
        // most likely class is in the hemi to be erased - erase everything
        if ((lh && (IS_LH_CLASS(max_label))) || (!lh && (IS_RH_CLASS(max_label)))) {
          int xn, yn, zn, r;

          GCApriorToNode(gca, x, y, z, &xn, &yn, &zn);
          gcan = &gca->nodes[xn][yn][zn];
          gcan->nlabels = 1;
          gcan->labels[0] = Unknown;
          for (r = 0; r < gca->ninputs; r++) {
            gcan->gcs[0].means[r] = 0;
          }
          gcap->nlabels = 1;
          gcap->labels[0] = Unknown;
        }
        else {
          for (n = 0; n < gcap->nlabels; n++) {
            if ((lh && (IS_LH_CLASS(gcap->labels[n]))) || (!lh && (IS_RH_CLASS(gcap->labels[n])))) {
              gcap->labels[n] = Unknown;
            }
          }
        }
      }
    }
  }

  if (lh)
    gca->flags |= GCA_NO_LH;
  else
    gca->flags |= GCA_NO_RH;
  return (NO_ERROR);
}

int GCAupdateDistributions(GCA *gca, MRI *mri, TRANSFORM *transform)
{
  int x, y, z, nvals, label, xv, yv, zv, n, node_label;
  float vals[MAX_GCA_INPUTS], *all_vals, med;
  GCA_NODE *gcan;

  all_vals = (float *)calloc(gca->node_width * gca->node_height * gca->node_depth, sizeof(float));
  for (label = 0; label <= MAX_CMA_LABEL; label++) {
    //    printf("updating means for label %s\n", cma_label_to_name(label)) ;
    nvals = 0;
    for (x = 0; x < gca->node_width; x++) {
      for (y = 0; y < gca->node_height; y++) {
        for (z = 0; z < gca->node_depth; z++) {
          gcan = &gca->nodes[x][y][z];
          if (GCAnodeToSourceVoxel(gca, mri, transform, x, y, z, &xv, &yv, &zv) == NO_ERROR) {
            node_label = gcaMaxPriorLabel(gca, mri, transform, xv, yv, zv);
            if (node_label == label) {
              load_vals(mri, xv, yv, zv, vals, gca->ninputs);
              all_vals[nvals++] = vals[0];
            }
          }
        }
      }
    }
    if (nvals == 0) {
      // printf("label %s not found, skipping\n", cma_label_to_name(label)) ;
      continue;
    }
    med = median(all_vals, nvals);
    printf("updating label %s to use median %f\n", cma_label_to_name(label), med);
    for (x = 0; x < gca->node_width; x++) {
      for (y = 0; y < gca->node_height; y++) {
        for (z = 0; z < gca->node_depth; z++) {
          gcan = &gca->nodes[x][y][z];
          if (GCAnodeToSourceVoxel(gca, mri, transform, x, y, z, &xv, &yv, &zv) == NO_ERROR) {
            for (n = 0; n < gcan->nlabels; n++) {
              if (gcan->labels[n] == label) {
                gcan->gcs[n].means[0] = med;
              }
            }
          }
        }
      }
    }
  }
  free(all_vals);
  return (NO_ERROR);
}

int GCAremoveLabel(GCA *gca, int label)
{
  int x, y, z, n, r;
  GCA_NODE *gcan;
  GCA_PRIOR *gcap;

  for (x = 0; x < gca->node_width; x++) {
    for (y = 0; y < gca->node_height; y++) {
      for (z = 0; z < gca->node_depth; z++) {
        gcan = &gca->nodes[x][y][z];
        for (n = 0; n < gcan->nlabels; n++) {
          if (gcan->labels[n] == label) {
            gcan->labels[n] = Unknown;
            for (r = 0; r < gca->ninputs; r++) {
              gcan->gcs[n].means[r] = 0;
            }
          }
        }
      }
    }
  }

  for (x = 0; x < gca->prior_width; x++) {
    for (y = 0; y < gca->prior_height; y++) {
      for (z = 0; z < gca->prior_depth; z++) {
        gcap = &gca->priors[x][y][z];
        for (n = 0; n < gcap->nlabels; n++) {
          if (gcap->labels[n] == label) {
            gcap->labels[n] = Unknown;
          }
        }
      }
    }
  }

  return (NO_ERROR);
}
int is_possible_wmsa(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z, int whalf)
{
  GCA_PRIOR *gcap;
  int xp, yp, zp, xi, yi, zi;
  float left_wmsa, right_wmsa, never_thresh, wmsa;

  for (xi = -whalf; xi <= whalf; xi++)
    for (yi = -whalf; yi <= whalf; yi++)
      for (zi = -whalf; zi <= whalf; zi++) {
        GCAsourceVoxelToPrior(gca, mri, transform, x + xi, y + yi, z + zi, &xp, &yp, &zp);
        gcap = getGCAP(gca, mri, transform, x, y, z);
        if (gcap == NULL) {
          continue;
        }
        never_thresh = .5 / gcap->total_training;

        left_wmsa = getPrior(gcap, Left_WM_hypointensities);
        right_wmsa = getPrior(gcap, Right_WM_hypointensities);
        wmsa = getPrior(gcap, WM_hypointensities);
        if (left_wmsa > never_thresh || right_wmsa > never_thresh || wmsa > never_thresh) {
          return (1);  // wmsa  occurred here
        }
      }
  return (0);  // wmsa has never occurred here
}

double cortex_prior(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z)
{
  GCA_PRIOR *gcap;
  int xp, yp, zp;
  float left_gm, right_gm;

  GCAsourceVoxelToPrior(gca, mri, transform, x, y, z, &xp, &yp, &zp);
  gcap = getGCAP(gca, mri, transform, x, y, z);
  if (gcap == NULL) {
    return (0.0);
  }
  left_gm = getPrior(gcap, Left_Cerebral_Cortex);
  right_gm = getPrior(gcap, Right_Cerebral_Cortex);
  return (left_gm + right_gm);
}

double wm_prior(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z)
{
  GCA_PRIOR *gcap;
  int xp, yp, zp;
  float left_wm, right_wm, wmsa;

  GCAsourceVoxelToPrior(gca, mri, transform, x, y, z, &xp, &yp, &zp);
  gcap = getGCAP(gca, mri, transform, x, y, z);
  if (gcap == NULL) {
    return (0.0);
  }
  wmsa = getPrior(gcap, WM_hypointensities);
  left_wm = getPrior(gcap, Left_Cerebral_White_Matter);
  right_wm = getPrior(gcap, Right_Cerebral_White_Matter);

  return (left_wm + right_wm + wmsa);
}
double csf_prior(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z)
{
  GCA_PRIOR *gcap;
  int xp, yp, zp, n;
  float prior;

  GCAsourceVoxelToPrior(gca, mri, transform, x, y, z, &xp, &yp, &zp);
  gcap = getGCAP(gca, mri, transform, x, y, z);
  if (gcap == NULL) {
    return (0.0);
  }
  for (prior = 0.0, n = 0; n < gcap->nlabels; n++)
    if (IS_CSF(gcap->labels[n])) {
      prior += gcap->priors[n];
    }

  return (prior);
}

double gm_prior(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z)
{
  GCA_PRIOR *gcap;
  int xp, yp, zp, n;
  float prior;

  GCAsourceVoxelToPrior(gca, mri, transform, x, y, z, &xp, &yp, &zp);
  gcap = getGCAP(gca, mri, transform, x, y, z);
  if (gcap == NULL) {
    return (0.0);
  }
  for (prior = 0.0, n = 0; n < gcap->nlabels; n++)
    if (IS_GRAY_MATTER(gcap->labels[n])) {
      prior += gcap->priors[n];
    }

  return (prior);
}

int GCAisLeftHemisphere(GCA *gca, MRI *mri, TRANSFORM *transform, int x, int y, int z)
{
  GCA_PRIOR *gcap;
  float lh_prior, rh_prior;
  int n;

  gcap = getGCAP(gca, mri, transform, x, y, z);
  for (lh_prior = rh_prior = 0.0, n = 0; n < gcap->nlabels; n++)
    if (IS_LH_CLASS(gcap->labels[n])) {
      lh_prior += gcap->priors[n];
    }
    else if (IS_RH_CLASS(gcap->labels[n])) {
      rh_prior += gcap->priors[n];
    }

  return (lh_prior > rh_prior);
}

/*
  \fn MRI *GCAsampleToVol(MRI *mri, GCA *gca, TRANSFORM *transform, MRI
  **seg, MRI *out) \brief Samples the intensities of the GCA into the
  space of the mri (which should be the conformed space of a subject,
  eg, orig.mgz). transform is the non-linear transform that maps
  between the conformed space and the atlas space (eg, talairach.m3z)
  and must have been inverted.  If seg is non-NULL then the most
  likely seg in the output space is returned. The output type will be
  float and will have as many frames as the gca has modes. Takes only
  a few min to run.

  MRI *orig, *gcaseg;
  orig = MRIread("orig.mgz");
  gca = GCAread("wmsa_new_eesmith.gca");
  transform = TransformRead("transforms/talairach.m3z");
  TransformInvert(transform, orig);
  gcavol = GCAsample(orig, gca, transform, &gcaseg, NULL);
  MRIwrite(gcavol,"gca.intensities.mgh");
  MRIwrite(gcaseg,"gcaseg.mgh");

  The output can be converted to uchar
    mri_convert gca.intensities.mgh  gca.mode.mgh --split -odt uchar --no_scale 1
  The result can be fed to mri_ca_normalize
    mri_ca_normalize -c cp.mgh -mask brainmask.mgz \
       gca.mode0000.mgh gca.mode0001.mgh gca.mode0002.mgh \
       wmsa_new_eesmith.gca talairach.m3z \
       gca.mode0000.can.mgh gca.mode0001.can.mgh gca.mode0002.can.mgh
  gca.mode000?.can.mgh should be pretty close to gca.mode0000.mgh. It is not
  exact, maybe because of the conversion to uchar (which mri_ca_normalize requires).
*/
MRI *GCAsampleToVol(MRI *mri, GCA *gca, TRANSFORM *transform, MRI **seg, MRI *out)
{
  int z;

  if (out == NULL) {
    out = MRIcloneBySpace(mri, MRI_FLOAT, gca->ninputs);
    if (out == NULL) return (NULL);
  }
  if (MRIdimMismatch(mri, out, 0)) {
    printf("ERROR: GCAsampleToVol(): output dimension mismatch\n");
    return (NULL);
  }
  if (seg != NULL) {
    if (*seg == NULL) {
      *seg = MRIcloneBySpace(mri, MRI_INT, 1);
      if (*seg == NULL) return (NULL);
    }
    if (MRIdimMismatch(mri, *seg, 0)) {
      printf("ERROR: GCAsampleToVol(): seg dimention mismatch\n");
      return (NULL);
    }
  }

  for (z = 0; z < mri->depth; z++) {
    int y, x, xn, yn, zn, n, max_n, err, f;
    GC1D *gc;
    GCA_NODE *gcan;
    GCA_PRIOR *gcap;
    double max_p;

    for (y = 0; y < mri->height; y++) {
      for (x = 0; x < mri->width; x++) {
        // Set to 0 in case it does not make it through to the end
        for (f = 0; f < gca->ninputs; f++) MRIsetVoxVal(out, x, y, z, f, 0);
        if (seg) MRIsetVoxVal(*seg, x, y, z, 0, 0);

        // Get the node for the xyz of this voxel
        err = GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn);
        if (err) continue;

        gcan = &gca->nodes[xn][yn][zn];
        gcap = getGCAP(gca, mri, transform, x, y, z);
        if (gcap == NULL) continue;

        // get structure with max prior prob
        max_p = 0;
        max_n = -1;
        for (n = 0; n < gcan->nlabels; n++) {
          if (getPrior(gcap, gcan->labels[n]) >= max_p) {
            max_p = getPrior(gcap, gcan->labels[n]);
            max_n = n;
          }
        }
        if (max_n < 0) continue; /* couldn't find any valid label at this location */

        // Set the output
        gc = &gcan->gcs[max_n];
        for (f = 0; f < gca->ninputs; f++) MRIsetVoxVal(out, x, y, z, f, gc->means[f]);
        if (seg) MRIsetVoxVal(*seg, x, y, z, 0, gcan->labels[max_n]);
      }
    }
  }

  return (out);
}

MRI *GCAsampleToVolWMSAprob(MRI *mri, GCA *gca, TRANSFORM *transform, MRI *out)
{
  int z;

  if (out == NULL) {
    out = MRIcloneBySpace(mri, MRI_FLOAT, 1);
    if (out == NULL) return (NULL);
  }
  if (MRIdimMismatch(mri, out, 0)) {
    printf("ERROR: GCAsampleToVol(): output dimension mismatch\n");
    return (NULL);
  }

  for (z = 0; z < mri->depth; z++) {
    int y, x, xn, yn, zn, n, wmsa_n, err;

    GCA_NODE *gcan;
    GCA_PRIOR *gcap;
    double wmsa_p;

    for (y = 0; y < mri->height; y++) {
      for (x = 0; x < mri->width; x++) {
        // Set to 0 in case it does not make it through to the end
        MRIsetVoxVal(out, x, y, z, 0, 0);

        // Get the node for the xyz of this voxel
        err = GCAsourceVoxelToNode(gca, mri, transform, x, y, z, &xn, &yn, &zn);
        if (err) continue;

        gcan = &gca->nodes[xn][yn][zn];
        gcap = getGCAP(gca, mri, transform, x, y, z);
        if (gcap == NULL) continue;

        // get priors for voxels with WMSA having highest prob
        wmsa_p = 0;
        wmsa_n = -1;
        for (n = 0; n < gcan->nlabels; n++) {
          if (IS_HYPO(gcan->labels[n])) {
            wmsa_p = getPrior(gcap, gcan->labels[n]);
            wmsa_n = n;
          }
        }
        if (wmsa_n < 0) continue; /* couldn't find any valid label at this location */

        // Set the output to be the prior at each WMSA voxel
        MRIsetVoxVal(out, x, y, z, 0, wmsa_p);
      }
    }
  }

  return (out);
}

// n is the half neighborhood size
// relabels a given voxel to WMSA if it was not already
// and if its mahalanobis distance to WMSAs is closer to that of WM
// as determined by the labeled voxels in a n-half neighborhood
int MRIwmsaHalo(MRI *mri_inputs, MRI *mri_labeled, int n)
{
  int h, label, nwmsa, ngm, x, y, z, wmsa_label, wm_label, gm_label, a, b, changed;
  // int nwm;
  double wm_mdist = 0, wmsa_mdist = 0, caudate_mdist = 0, vent_mdist = 0;
  VECTOR *v_vals = NULL, *wm_means, *wmsa_means, *caudate_means, *vent_means;
  MATRIX *wmsa_cov, *wm_cov, *caudate_cov, *vent_cov, *m_inv_cov_wmsa = NULL, *m_inv_cov_wm = NULL,
                                                      *m_inv_cov_caudate = NULL, *m_inv_cov_vent = NULL;
  MRI *label_copy = NULL;

  // Use ENTIRE wm
  MRIcomputeWMMeansandCovariances(mri_inputs, mri_labeled, &wm_cov, &wm_means);
  m_inv_cov_wm = MatrixAlloc(wm_cov->rows, wm_cov->cols, wm_cov->type);
  MatrixInverse(wm_cov, m_inv_cov_wm);

  MRIcomputeWMSAMeansandCovariances(mri_inputs, mri_labeled, &wmsa_cov, &wmsa_means);
  m_inv_cov_wmsa = MatrixAlloc(wmsa_cov->rows, wmsa_cov->cols, wmsa_cov->type);
  MatrixInverse(wmsa_cov, m_inv_cov_wmsa);

  MRIcomputeCaudateMeansandCovariances(mri_inputs, mri_labeled, &caudate_cov, &caudate_means);
  m_inv_cov_caudate = MatrixAlloc(caudate_cov->rows, caudate_cov->cols, caudate_cov->type);
  MatrixInverse(caudate_cov, m_inv_cov_caudate);

  MRIcomputeVentMeansandCovariances(mri_inputs, mri_labeled, &vent_cov, &vent_means);
  m_inv_cov_vent = MatrixAlloc(vent_cov->rows, vent_cov->cols, vent_cov->type);
  MatrixInverse(vent_cov, m_inv_cov_vent);

  changed = 1;
  // b = 1;
  // while(changed > 0)
  //{
  for (b = 0; b < 2; b++) {
    changed = 0;
    v_vals = VectorAlloc(mri_inputs->nframes, MATRIX_REAL);

    printf("copying label file...\n");
    label_copy = MRIcopy(mri_labeled, NULL);
    printf("Finished copying label file...\n");

    for (h = 0; h < 2; h++) {
      if (h == 0) {
        wmsa_label = Right_WM_hypointensities;
        wm_label = Right_Cerebral_White_Matter;
        gm_label = Right_Cerebral_Cortex;
      }
      else {
        wmsa_label = Left_WM_hypointensities;
        wm_label = Left_Cerebral_White_Matter;
        gm_label = Left_Cerebral_Cortex;
      }

      for (x = 0; x < mri_labeled->width; x++) {
        for (y = 0; y < mri_labeled->height; y++) {
          for (z = 0; z < mri_labeled->depth; z++) {
            label = MRIvox(label_copy, x, y, z);

            // nwm = 
            MRIlabelsInNbhd(label_copy, x, y, z, n, wm_label);
            nwmsa = MRIlabelsInNbhd(label_copy, x, y, z, n, wmsa_label);
            ngm = MRIlabelsInNbhd(label_copy, x, y, z, 1, gm_label);

            if (nwmsa < 1 || ngm > 0)
              // if (nwmsa < 1 || nwm < 1 || ngm > 0)
              // if (nwm < n || ngm > 0)
              continue;

            for (a = 0; a < mri_inputs->nframes; a++) {
              VECTOR_ELT(v_vals, a + 1) = MRIgetVoxVal(mri_inputs, x, y, z, a);
            }
            // printf("Computing means and covariances\n");
            // fflush(stdout);
            // MRIcomputeNbhdMeansandCovariances(mri_inputs, label_copy, wm_label, x, y, z, n, &wm_cov, &wm_means);
            // printf("Computing means and covariances\n");
            // fflush(stdout);
            // MRIcomputeNbhdMeansandCovariances(mri_inputs, label_copy, wmsa_label, x, y, z, n, &wmsa_cov,
            // &wmsa_means);
            // printf("finished\n");
            // fflush(stdout);

            // m_inv_cov_wm = MatrixAlloc(wm_cov->rows, wm_cov->cols, wm_cov->type) ;
            // m_inv_cov_wmsa = MatrixAlloc(wmsa_cov->rows, wmsa_cov->cols, wmsa_cov->type) ;

            // printf("Calculating WMSA inv matrix\n");
            // fflush(stdout);
            // MatrixInverse(wmsa_cov, m_inv_cov_wmsa) ;
            // printf("Done calculating WMSA inv matrix\n") ;
            // fflush(stdout);

            // printf("Calculating WM inv matrix\n");
            // fflush(stdout);
            // MatrixInverse(wm_cov, m_inv_cov_wm) ;
            // printf("Done calculating WM inv matrix\n") ;
            // fflush(stdout);

            wmsa_mdist = MatrixMahalanobisDistance(wmsa_means, m_inv_cov_wmsa, v_vals);

            if (x == Gx && y == Gy && z == Gz)
              printf("WM_mdist is %f and WMSA_mdist is %f for voxel %d, %d, %d\n", wm_mdist, wmsa_mdist, x, y, z);

            if (IS_WM(label)) {
              wm_mdist = MatrixMahalanobisDistance(wm_means, m_inv_cov_wm, v_vals);
              if ((VECTOR_ELT(v_vals, 1) < VECTOR_ELT(wm_means, 1)) &&
                  (VECTOR_ELT(v_vals, 2) > VECTOR_ELT(wm_means, 2)) &&
                  (VECTOR_ELT(v_vals, 3) > VECTOR_ELT(wm_means, 3)) &&
                  ((7 * wmsa_mdist < wm_mdist) || (3 * wmsa_mdist < wm_mdist && nwmsa > 3))) {
                MRIsetVoxVal(mri_labeled, x, y, z, 0, wmsa_label);
                changed++;
              }
            }
            if (IS_CAUDATE(label)) {
              caudate_mdist = MatrixMahalanobisDistance(caudate_means, m_inv_cov_caudate, v_vals);
              if ((VECTOR_ELT(v_vals, 1) < VECTOR_ELT(caudate_means, 1)) &&
                  (VECTOR_ELT(v_vals, 2) > VECTOR_ELT(caudate_means, 2)) &&
                  (VECTOR_ELT(v_vals, 3) > VECTOR_ELT(caudate_means, 3)) &&
                  ((7 * wmsa_mdist < caudate_mdist) || (2 * wmsa_mdist < caudate_mdist && nwmsa > 3))) {
                MRIsetVoxVal(mri_labeled, x, y, z, 0, wmsa_label);
                changed++;
              }
            }

            if (IS_VENTRICLE(label)) {
              vent_mdist = MatrixMahalanobisDistance(vent_means, m_inv_cov_vent, v_vals);
              printf("T1 mean for vent is %f\n", VECTOR_ELT(vent_means, 1));
              printf("PD mean for vent is %f\n", VECTOR_ELT(vent_means, 2));
              printf("T2 mean for vent is %f\n", VECTOR_ELT(vent_means, 3));
              printf("T1 val for voxel %d, %d, %d, is %f\n", x, y, z, VECTOR_ELT(v_vals, 1));
              // printf("Voxel %d, %d, %d is ventricle...\n", x, y, z);
              // printf("Ventricle m dist is: %f\n", vent_mdist) ;
              // printf("WMSA m dist is: %f\n", wmsa_mdist) ;
              if ((VECTOR_ELT(v_vals, 1) > 2 * VECTOR_ELT(vent_means, 1)) &&
                  ((7 * wmsa_mdist < vent_mdist) || (wmsa_mdist < vent_mdist && nwmsa > 3)))
              //  if( (VECTOR_ELT(v_vals,1) > 2*VECTOR_ELT(vent_means,1)) &&
              //   ((7*wmsa_mdist<vent_mdist)))
              {
                MRIsetVoxVal(mri_labeled, x, y, z, 0, wmsa_label);
                changed++;
              }
            }
          }
        }
      }
    }
    printf("Iteration %d, %d halo voxels changed to WMSA\n", b, changed);
    fflush(stdout);
    // b++;
    MRIfree(&label_copy);
  }
  MatrixFree(&wm_cov);
  MatrixFree(&wm_means);
  MatrixFree(&m_inv_cov_wm);
  MatrixFree(&wmsa_cov);
  MatrixFree(&wmsa_means);
  MatrixFree(&m_inv_cov_wmsa);
  return (NO_ERROR);
}

WMSA *WMSAalloc(int nreftissues)
{
  WMSA *wmsa;
  wmsa = (WMSA *)calloc(sizeof(WMSA), 1);
  wmsa->nreftissues = nreftissues;
  wmsa->reftissues = (int *)calloc(nreftissues, sizeof(int));
  wmsa->hardthresh = (double *)calloc(nreftissues, sizeof(double));
  wmsa->softthresh = (double *)calloc(nreftissues, sizeof(double));
  return (wmsa);
}

int MRIwmsaHalo2(WMSA *wmsa)
{
  int ok, i, j, k, a, b, c, d, x, y, z, changed, nwmsa1, nwmsa2, nwmsa, ngm1, ngm2, ngm, nmatches, label;
  int wmsa_labels[] = {78, 79};
  int gm_labels[] = {3, 42};
  double reftissue_dist = 0, wmsa_mdist;
  MRI *label_copy = NULL;
  MATRIX *v_vals;
  VECTOR *reftissue_means[100], *wmsa_means = NULL;
  MATRIX *reftissue_covs[100], *reftissue_inv_covs[100], *wmsa_covs = NULL, *wmsa_inv_covs = NULL;
  int nreftissues = wmsa->nreftissues;
  int nmodalities = wmsa->modalities->nframes;
  int wmsacomp[nreftissues][nmodalities];
  int llabel[1];

  // For each ref tissue, compute the means for each modality and covariance matrix across all modalities
  // 0 = dont use any neighbors
  printf("size %d\n", (int)sizeof(llabel));
  for (i = 0; i < nreftissues; i++) {
    llabel[0] = wmsa->reftissues[i];
    MRIcomputeLabelMeansandCovariances(
        wmsa->modalities, wmsa->seg, &reftissue_covs[i], &reftissue_means[i], &llabel[0], 1, 0);
    reftissue_inv_covs[i] = MatrixInverse(reftissue_covs[i], NULL);
  }

  // Means and Cov of WMSA voxels
  MRIcomputeLabelMeansandCovariances(wmsa->modalities, wmsa->seg, &wmsa_covs, &wmsa_means, wmsa_labels, 2, 0);
  wmsa_inv_covs = MatrixInverse(wmsa_covs, NULL);

  // For each modality, we need to figure out if the mean of ref tissues
  // are greater or less than the means of WMSA. 1 means greater.
  for (j = 0; j < nreftissues; j++) {
    for (k = 0; k < nmodalities; k++) {
      if (reftissue_means[j]->rptr[1][k + 1] < wmsa_means->rptr[1][k + 1])
        wmsacomp[j][k] = 0;
      else
        wmsacomp[j][k] = 1;
    }
  }

  // set up vector for voxel intensities for each modality
  v_vals = VectorAlloc(wmsa->modalities->nframes, MATRIX_REAL);

  // perform iterations of halo-growing (default should be 3)
  for (b = 0; b < wmsa->niters; b++) {
    changed = 0;

    // create a static copy of the label file so we can dynamically update the true label file
    label_copy = MRIcopy(wmsa->seg, NULL);

    for (x = 0; x < wmsa->seg->width; x++) {
      for (y = 0; y < wmsa->seg->height; y++) {
        for (z = 0; z < wmsa->seg->depth; z++) {
          // get segmentation label for this voxel
          label = MRIvox(label_copy, x, y, z);

          // If label is already WMSA skip it
          if ((label == wmsa_labels[0]) || (label == wmsa_labels[1])) continue;

          // Only continue with voxel if its label matches one of the user's seg inputs
          ok = 0;
          for (c = 0; c < nreftissues; c++) {
            if (label == wmsa->reftissues[c]) {
              ok = 1;
              break;
            }
          }
          if (!ok) continue;
          // To get here, the voxel must be one of the reference tissues
          // c is now the nth reference tissue

          // find number WMSA and GM neighbors
          nwmsa1 = MRIlabelsInNbhd(label_copy, x, y, z, wmsa->nbrwhalf, wmsa_labels[0]);  // lh
          nwmsa2 = MRIlabelsInNbhd(label_copy, x, y, z, wmsa->nbrwhalf, wmsa_labels[1]);  // rh
          nwmsa = nwmsa1 + nwmsa2;                                       // one of them will probably be 0
          ngm1 = MRIlabelsInNbhd(label_copy, x, y, z, 1, gm_labels[0]);  // lh
          ngm2 = MRIlabelsInNbhd(label_copy, x, y, z, 1, gm_labels[1]);  // rh
          ngm = ngm1 + ngm2;                                             // one of them will probably be 0

          // Ignore this voxel if it has 0 WMSA neighbors or a single GM neighbor
          if (nwmsa < 1 || ngm > 0) continue;

          // Extract voxel intensities for each modality
          for (a = 0; a < nmodalities; a++) v_vals->rptr[a + 1][1] = MRIgetVoxVal(wmsa->modalities, x, y, z, a);

          // If this vox is really a WMSA, then each modality at this voxel should have
          // the same intensity relationship with WMSA as the reference tissue.
          // Use wmsacomp matrix to check if voxel vals should be higher or lower than seg of interest
          // If wmsacomp=1, then
          nmatches = 0;
          for (d = 0; d < nmodalities; d++) {
            if (wmsacomp[c][d] == 0) {  // meanWMSA > meanRefTissue
              if (v_vals->rptr[d + 1][1] < reftissue_means[c]->rptr[d + 1][1]) nmatches = nmatches + 1;
            }
            else  // meanWMSA < meanRefTissue
                if (v_vals->rptr[d + 1][1] > reftissue_means[c]->rptr[d + 1][1])
              nmatches = nmatches + 1;

          }  // modality
          if (nmatches != nmodalities) continue;

          // Compute MD of this voxel from WMSA
          wmsa_mdist = MatrixMahalanobisDistance(wmsa_means, wmsa_inv_covs, v_vals);

          // Compute distance to ref tissue
          reftissue_dist = MatrixMahalanobisDistance(reftissue_means[c], reftissue_inv_covs[c], v_vals);

          // Impose distance contraints
          ok = 0;
          // Fraction of Dist to WMSA is less than distance to ref tissue
          if (wmsa->hardthresh[c] * wmsa_mdist < reftissue_dist) ok = 1;
          // Fraction of Dist to WMSA is less than distance to ref tissue and
          // there are enough WMSA neighbors (default should be 3)
          if (wmsa->softthresh[c] * wmsa_mdist < reftissue_dist && nwmsa > wmsa->nbrthresh) ok = 1;
          if (!ok) continue;

          if (nwmsa1 > nwmsa2) {  // left hemisphere
            MRIsetVoxVal(wmsa->seg, x, y, z, 0, wmsa_labels[0]);
            changed++;
          }
          else {  // right hemi
            MRIsetVoxVal(wmsa->seg, x, y, z, 0, wmsa_labels[1]);
            changed++;
          }
        }  // for z
      }    // for y
    }      // for x

    printf("Iteration %d, %d halo voxels changed to WMSA\n", b, changed);
    fflush(stdout);
    // if changed == 0, should we bail out?

    MRIfree(&label_copy);
  }  // for b

  for (i = 0; i < nreftissues; i++) {
    MatrixFree(&reftissue_covs[i]);
    MatrixFree(&reftissue_means[i]);
    MatrixFree(&reftissue_inv_covs[i]);
  }
  MatrixFree(&v_vals);
  MatrixFree(&wmsa_covs);
  MatrixFree(&wmsa_means);
  MatrixFree(&wmsa_inv_covs);
  return (NO_ERROR);
}
void GCAinitLabelsFromMRI(GCA *gca, MRI *mri_labels)
{
  int x, y, z, label;
  GCA_PRIOR *gcap;
  GCA_NODE *gcan;

  if (gca->width != mri_labels->width || gca->height != mri_labels->height || gca->depth != mri_labels->depth)
    ErrorExit(ERROR_BADPARM, "GCAinitLabelsFromMRI: GCA and MRI must have same dimensions");

  for (x = 0; x < gca->width; x++)
    for (y = 0; y < gca->height; y++)
      for (z = 0; z < gca->depth; z++) {
        label = nint(MRIgetVoxVal(mri_labels, x, y, z, 0));
        gcan = &gca->nodes[x][y][z];
        gcap = &gca->priors[x][y][z];
        gcap->labels[0] = label;
        gcap->priors[0] = 1.0;
        gcan->labels[0] = label;
        gcan->gcs = alloc_gcs(1, 0, 1);
      }
}

/*!
\fn unsigned short *GCAmergeLabelLists(unsigned short *labels1, int nlabels1, unsigned short *labels2, int nlabels2, int *pnlist)
\brief Merges two unsigned short integer lists into a single list. Good for working with GCA label lists.
*/
unsigned short *GCAmergeLabelLists(unsigned short *labels1, int nlabels1, unsigned short *labels2, int nlabels2, int *pnlist)
{
  int labellist[500],nlist,n,m,hit;
  unsigned short *ulist;

  // Count and make a list of the unique labels from both gcaps
  nlist = 0;
  for(n=0; n < nlabels1; n++){
    hit = 0;
    for(m=0; m < nlist; m++){
      if(labels1[n] == labellist[m]){
	hit = 1;
	break;
      }
    }
    if(!hit){
      labellist[nlist] = labels1[n];
      nlist++;
    }
  }
  for(n=0; n < nlabels2; n++){
    hit = 0;
    for(m=0; m < nlist; m++){
      if(labels2[n] == labellist[m]){
	hit = 2;
	break;
      }
    }
    if(!hit){
      labellist[nlist] = labels2[n];
      nlist++;
    }
  }

  *pnlist = nlist;
  ulist = (unsigned short*)calloc(sizeof(unsigned short),nlist);
  for(m=0; m < nlist; m++)
    ulist[m] = labellist[m];

  return(ulist);
}

/*!
\fn GCA_PRIOR *GCAPcopy(GCA_PRIOR *gcap, int symmetrize, GCA_PRIOR *gcapcopy)
\brief Make a copy of the GCA prior. If symmetrize, then the copy will have
label codes set to the code of the contralateral structural.
*/
GCA_PRIOR *GCAPcopy(GCA_PRIOR *gcap, int symmetrize, GCA_PRIOR *gcapcopy)
{
  int n;

  if(gcapcopy==NULL){
    // Alloc the new gcap to hold the merged
    gcapcopy = (GCA_PRIOR *) calloc(sizeof(GCA_PRIOR),1);
  } 
  else {
    if(gcapcopy->labels) free(gcapcopy->labels);
    if(gcapcopy->priors) free(gcapcopy->priors);
  }
  gcapcopy->nlabels = gcap->nlabels;
  gcapcopy->labels = (unsigned short*)calloc(sizeof(unsigned short),gcap->nlabels);
  gcapcopy->priors = (float*)calloc(sizeof(float),gcap->nlabels);

  for(n=0; n < gcap->nlabels; n++){
    if(symmetrize)
      gcapcopy->labels[n] = MRIasegContraLatLabel(gcap->labels[n]);
    else
      gcapcopy->labels[n] = gcap->labels[n];
    gcapcopy->priors[n] = gcap->priors[n];
  }
  return(gcapcopy);
}

int GCAPfree(GCA_PRIOR **pgcap)
{
  GCA_PRIOR *gcap = *pgcap;
  free(gcap->labels);
  free(gcap->priors);
  free(*pgcap);
  *pgcap = NULL;
  return(0);
}

/*!
\fn GCA_PRIOR *GCAPmerge(GCA_PRIOR *gcap1, GCA_PRIOR *gcap2, GCA_PRIOR *gcapm)
\brief Merge two GCA priors to create a new prior.
*/
GCA_PRIOR *GCAPmerge(GCA_PRIOR *gcap1, GCA_PRIOR *gcap2, GCA_PRIOR *gcapm)
{
  int n,m,nlabels;

  if(gcapm==NULL){
    // Alloc the new gcap to hold the merged
    gcapm = (GCA_PRIOR *) calloc(sizeof(GCA_PRIOR),1);
  } 
  else {
    if(gcapm->labels) free(gcapm->labels);
    if(gcapm->priors) free(gcapm->priors);
  }

  // Count and make a list of the unique labels from both gcaps
  gcapm->labels =  GCAmergeLabelLists(gcap1->labels, gcap1->nlabels, gcap2->labels, gcap2->nlabels, &nlabels);
  gcapm->priors = (float*)calloc(sizeof(float),nlabels);
  gcapm->nlabels = nlabels;
  gcapm->max_labels = gcapm->nlabels;
  gcapm->total_training = gcap1->total_training;

  // Assign the merged values. Divide by two to account for the merge
  // and assure that the sum(priors)=1
  for(m=0; m < gcapm->nlabels; m++){
    for(n=0; n < gcap1->nlabels; n++){    
      if(gcap1->labels[n] == gcapm->labels[m]){
	gcapm->priors[m] += gcap1->priors[n]/2.0;
	break;
      }
    }
    for(n=0; n < gcap2->nlabels; n++){    
      if(gcap2->labels[n] == gcapm->labels[m]){
	gcapm->priors[m] += gcap2->priors[n]/2.0;
	break;
      }
    }
  }

  return(gcapm);
}

int GCAPprint(FILE *fp, GCA_PRIOR *gcap)
{
  int n,contralabel;
  fprintf(fp,"nlabels %d\n",gcap->nlabels);
  fprintf(fp,"total_training %d\n",gcap->total_training);
  fprintf(fp,"max_labels %d\n",gcap->max_labels);
  for(n=0; n < gcap->nlabels; n++){
    // labelcode priorprob (note: sum(priorprob) should = 1)
    contralabel = MRIasegContraLatLabel(gcap->labels[n]);
    fprintf(fp,"  %d/%d  %f\n",gcap->labels[n],contralabel,gcap->priors[n]);
  }
  fflush(fp);
  return(0);
}

int GCANprint(FILE *fp, GCA_NODE *node, int ninputs)
{
  int n;
  fprintf(fp,"nlabels %d\n",node->nlabels);
  fprintf(fp,"total_training %d\n",node->total_training);
  fprintf(fp,"max_labels %d\n",node->max_labels);
  for(n=0; n < node->nlabels; n++){
    printf("nthLabel=%d LabelCode=%d (contra=%d)==========\n",n,node->labels[n],MRIasegContraLatLabel(node->labels[n]));
    GC1Dprint(fp, &(node->gcs[n]), ninputs);
  }
  return(0);
}

int GC1Dprint(FILE *fp, GC1D *gc1d, int ninputs)
{
  int r,c,v;
  fprintf(fp,"n_just_priors %d\n",gc1d->n_just_priors);
  fprintf(fp,"ntraining %d\n",gc1d->ntraining);
  fprintf(fp,"regularized %d\n",gc1d->regularized);
  v = 0;
  for (r = 0; r < ninputs; r++) {
    fprintf(fp,"  Mode %d   Mean %g  CoVar ",r,gc1d->means[r]);
    for (c = r; c < ninputs; c++) {
      fprintf(fp," %g ",gc1d->covars[v] );
      v++;
    }
    fprintf(fp,"\n");
  }

  for(r=0; r < GIBBS_NEIGHBORS; r++){
    fprintf(fp,"  Nbr %d  nlabels %2d ",r,gc1d->nlabels[r]);
    for(c=0; c < gc1d->nlabels[r]; c++){
      // Node: sum(gc1d->label_priors[r][c]) over c should = 1
      fprintf(fp," (%d,%g) ",gc1d->labels[r][c],gc1d->label_priors[r][c]);
    }
    fprintf(fp," \n");
  }
  fflush(fp);
  
  return(0);
}

/*!
\fn GC1D *GC1Dcopy(GC1D *gc, int ninputs, int symmetrize, GC1D *gccopy)
\brief Makes a copy of the GC1D. If symmetrize=1, then the GC is 
symmetrized along the way, meaning that labels are set to their
contralateral counterparts, and the 1st and 2nd MRF/GIBBS neighbors
are swapped. ninputs = gca->ninputs
*/
GC1D *GC1Dcopy(GC1D *gc, int ninputs, int symmetrize, GC1D *gccopy)
{
  int r,c,v,rr;

  if(gccopy == NULL) {
    gccopy = (GC1D*) calloc(sizeof(GC1D),1);
  }
  else {
    // Free stuff if needed
    if(gccopy->means)  free(gccopy->means);
    if(gccopy->covars) free(gccopy->covars);
    if(gccopy->nlabels) free(gccopy->nlabels);
    if(gccopy->label_priors){
      for(r=0; r < GIBBS_NEIGHBORS; r++){
	if(gccopy->label_priors[r]) free(gccopy->label_priors[r]);
      }
      free(gccopy->label_priors);
    }
    if(gccopy->labels){
      for(r=0; r < GIBBS_NEIGHBORS; r++){
	if(gccopy->labels[r]) free(gccopy->labels[r]);
      }
      free(gccopy->labels);
    }
  }
  gccopy->ntraining   = gc->ntraining;
  gccopy->regularized = gc->regularized;
  gccopy->n_just_priors = gc->n_just_priors;

  // Copy means and covars 
  gccopy->means  = (float*)calloc(sizeof(float),ninputs);
  gccopy->covars = (float*)calloc(sizeof(float),(ninputs*ninputs-ninputs)/2+ninputs);
  v = 0;
  for (r = 0; r < ninputs; r++) {
    gccopy->means[r] = gc->means[r];
    for (c = r; c < ninputs; c++) {
      gccopy->covars[r] = gc->covars[r];
      v++;
    }
  }

  // Copy MRF/Gibbs stuff. If GIBBS_NEIGHBORS != 6, then problems
  gccopy->nlabels = (short*)calloc(sizeof(short),GIBBS_NEIGHBORS);
  gccopy->labels = (unsigned short**)calloc(sizeof(unsigned short*),GIBBS_NEIGHBORS);
  gccopy->label_priors = (float**)calloc(sizeof(float*),GIBBS_NEIGHBORS);

  for(r=0; r < GIBBS_NEIGHBORS; r++){
    rr = r;
    if(symmetrize){
      // When symmetrizing, the neighbor order reverses for the column face neighbors 
      // It is assumed that the col nhbrs are the first two see xnbr_offset[] = {1, -1, 0, 0, 0, 0};
      if(r == 0) rr = 1;
      if(r == 1) rr = 0;
    }
    gccopy->nlabels[r] = gc->nlabels[rr];
    gccopy->labels[r] = (unsigned short *)calloc(sizeof(unsigned short),gccopy->nlabels[r]);
    for(c=0; c < gc->nlabels[rr]; c++){
      if(symmetrize){
	// When symmetrizing, the label has to be changed to the contra lateral 
	gccopy->labels[r][c] = MRIasegContraLatLabel(gc->labels[rr][c]);
      }
      else{
	gccopy->labels[r][c] = gc->labels[rr][c];
      }
    }

    // Copy the probabilties. 
    gccopy->label_priors[r] = (float *)calloc(sizeof(float),gccopy->nlabels[r]);
    for(c=0; c < gccopy->nlabels[r]; c++){
      gccopy->label_priors[r][c] = gc->label_priors[rr][c];
    }
  }

  return(gccopy);
}

/*!
\fn GC1D *GC1Dmerge(GC1D *gc1, GC1D *gc2, int ninputs, GC1D *gcm)
\brief Merges two GC1D to create a new GC. ninputs = gca->ninputs
*/
GC1D *GC1Dmerge(GC1D *gc1, GC1D *gc2, int ninputs, GC1D *gcm)
{
  int r,c,v,nlabels,n;

  if(gcm == NULL) gcm = (GC1D*) calloc(sizeof(GC1D),1);

  gcm->ntraining = gc1->ntraining + gc2->ntraining; //??
  gcm->regularized = gc1->regularized; // should check to make sure consistent
  gcm->n_just_priors = gc1->n_just_priors; // ???

  // Free stuff if needed
  if(gcm->means)  free(gcm->means);
  if(gcm->covars) free(gcm->covars);
  if(gcm->nlabels) free(gcm->nlabels);
  if(gcm->label_priors){
    for(r=0; r < GIBBS_NEIGHBORS; r++){
      if(gcm->label_priors[r]) free(gcm->label_priors[r]);
    }
    free(gcm->label_priors);
  }
  if(gcm->labels){
    for(r=0; r < GIBBS_NEIGHBORS; r++){
      if(gcm->labels[r]) free(gcm->labels[r]);
    }
    free(gcm->labels);
  }

  // Merge means and covars by averaging
  gcm->means  = (float*)calloc(sizeof(float),ninputs);
  gcm->covars = (float*)calloc(sizeof(float),(ninputs*ninputs-ninputs)/2+ninputs);
  v = 0;
  for (r = 0; r < ninputs; r++) {
    gcm->means[r] = (gc1->means[r]+gc2->means[r])/2.0;
    for (c = r; c < ninputs; c++) {
      gcm->covars[r] = (gc1->covars[r]+gc2->covars[r])/2.0;
      v++;
    }
  }

  // Merge MRF/Gibbs stuff. If GIBBS_NEIGHBORS != 6 then problems
  gcm->nlabels = (short*)calloc(sizeof(short),GIBBS_NEIGHBORS);
  gcm->labels = (unsigned short**)calloc(sizeof(unsigned short*),GIBBS_NEIGHBORS);
  gcm->label_priors = (float**)calloc(sizeof(float*),GIBBS_NEIGHBORS);

  for(r=0; r < GIBBS_NEIGHBORS; r++){
    // Merge the labels
    gcm->labels[r] = GCAmergeLabelLists(gc1->labels[r], gc1->nlabels[r], gc2->labels[r], gc2->nlabels[r], &nlabels);
    gcm->nlabels[r] = nlabels;
    gcm->label_priors[r] = (float*)calloc(sizeof(float),nlabels);

    // Merge the probabilties. Dividing by 2 assures that psum=1
    for(c=0; c < gcm->nlabels[r]; c++){
      for(n=0; n < gc1->nlabels[r]; n++){
	if(gcm->labels[r][c] == gc1->labels[r][n]){
	  gcm->label_priors[r][c] += gc1->label_priors[r][n]/2.0;
	  break;
	}
      }
      for(n=0; n < gc2->nlabels[r]; n++){
	if(gcm->labels[r][c] == gc2->labels[r][n]){
	  gcm->label_priors[r][c] += gc2->label_priors[r][n]/2.0;
	  break;
	}
      }
    }
  }

  return(gcm);
}

/*!
\fn GCA_NODE *GCANmerge(GCA_NODE *node1, GCA_NODE *node2, int ninputs, int symmetrize, GCA_NODE *nodem)
\brief Merges two nodes to create a new node. If symmetrize=1, then the second node is symmetrized
*/
GCA_NODE *GCANmerge(GCA_NODE *node1, GCA_NODE *node2, int ninputs, int symmetrize, GCA_NODE *nodem)
{
  int n, n1, n2, n1ok, n2ok;
  unsigned short *node2labels;
  GC1D *gc2;

  if(symmetrize){
    // Create a new array, convert the node2 labels to contralateral counterparts
    node2labels = (unsigned short*)calloc(sizeof(unsigned short),node2->nlabels);
    for(n=0; n < node2->nlabels; n++){
      node2labels[n] = MRIasegContraLatLabel(node2->labels[n]);
    }
  }
  else node2labels = node2->labels;

  if(nodem==NULL){
    // Alloc the merged load
    nodem = (GCA_NODE*)calloc(sizeof(GCA_NODE),1);
  } else {
    if(nodem->labels) free(nodem->labels);
    if(nodem->gcs) GC1Dfree(&nodem->gcs, ninputs);
  }
  // Get list of unique labels
  nodem->labels = GCAmergeLabelLists(node1->labels, node1->nlabels, node2labels, node2->nlabels, &nodem->nlabels);
  // Alloc the GC1D array, one for each label
  nodem->gcs = (GC1D*)calloc(sizeof(GC1D),nodem->nlabels);

  nodem->total_training = node1->total_training + node2->total_training;
  nodem->max_labels = nodem->nlabels;

  // Go through each label
  for(n=0; n < nodem->nlabels; n++){
    // Determine whether this label is represented in node1
    n1ok = 0;
    for(n1=0; n1 < node1->nlabels; n1++){
      if(nodem->labels[n] == node1->labels[n1]){
	n1ok = 1;
	break;
      }
    }
    // Determine whether this label is represented in node2
    n2ok = 0;
    for(n2=0; n2 < node2->nlabels; n2++){
      if(nodem->labels[n] == node2labels[n2]) {
	n2ok = 1;
	break;
      }
    }
    if(n1ok && n2ok){
      if(symmetrize)
	gc2 = GC1Dcopy(&node2->gcs[n2], ninputs, symmetrize, NULL);
      else
	gc2 = &node2->gcs[n2];
      GC1Dmerge(&node1->gcs[n1], gc2, ninputs, &nodem->gcs[n]);
      if(symmetrize) GC1Dfree(&gc2, ninputs);
      continue;
    }
    if(n1ok){ // This label not reprsented in node2, just copy GC from node1
      GC1Dcopy(&node1->gcs[n1], ninputs, 0, &nodem->gcs[n]);
      continue;
    }
    if(n2ok){ // This label not reprsented in node1, just copy GC from node2, sym if needed
      GC1Dcopy(&node2->gcs[n2], ninputs, symmetrize, &nodem->gcs[n]);
      continue;
    }
  }

  if(symmetrize) free(node2labels);

  return(nodem);
}

int GC1Dfree(GC1D **pgc, int ninputs)
{
  GC1D *gc = *pgc;
  int r;
  free(gc->means);
  free(gc->covars);
  free(gc->nlabels);
  for(r=0; r < GIBBS_NEIGHBORS; r++){
    free(gc->labels[r]);
    free(gc->label_priors[r]);
  }
  free(*pgc);
  *pgc = NULL;
  return(0);
}

/*!
\fn GCA *GCAsymmetrize(GCA *gca)
\brief Symmetrizes atlas by making sure that the statistics are symmetric. New
statistics are computed effectively by averaging the non-flipped with the flipped.
Takes into account the changing of segmentation labels to the corresponding 
contralateral code.
 */
GCA *GCAsymmetrize(GCA *gca)
{
  GCA *gcasym;
  GCA_NODE  *node,  *contranode, *symnode;
  GCA_PRIOR *prior, *contraprior, *symprior;
  int c,r,s,csym;

  gcasym = GCAalloc(gca->ninputs, gca->prior_spacing, gca->node_spacing, 
		    gca->width, gca->height, gca->depth, gca->flags) ;

  printf("width %d\n",gca->width);
  printf("node  dim %d\n",gca->node_width);
  printf("prior dim %d\n",gca->prior_width);

  gcasym->max_label = 0;
  printf("Syming nodes\n");
  for(c=0; c < gca->node_width; c++){
    csym = gca->node_width - c - 1;
    for(r=0; r < gca->node_height; r++){
      for(s=0; s < gca->node_depth; s++){
	node       = &(gca->nodes[c][r][s]);
	contranode = &(gca->nodes[csym][r][s]);
	symnode    = &(gcasym->nodes[c][r][s]);
	GCANmerge(node, contranode, gca->ninputs, 1, symnode);
	if(symnode->max_labels > gcasym->max_label) gcasym->max_label = symnode->max_labels;
      } // s
    } // r
  } //c 

  printf("Syming priors\n");
  for(c=0; c < gca->prior_width; c++){
    csym = gca->prior_width - c - 1;
    for(r=0; r < gca->prior_height; r++){
      for(s=0; s < gca->prior_depth; s++){
	prior       = &(gca->priors[c][r][s]);
	contraprior = GCAPcopy(&(gca->priors[csym][r][s]),1,NULL);
	symprior    = &(gcasym->priors[c][r][s]);
	GCAPmerge(prior, contraprior, symprior);
	if(0 && (c==28 || csym==28) && r==51 && s==51){
	  printf("prior %d/%d %d %d -------------------------------------\n",c,csym,r,s);
	  GCAPprint(stdout, prior);
	  printf("contra prior -------------------------------------\n");
	  GCAPprint(stdout, contraprior);
	  printf("sym prior -------------------------------------\n");
	  GCAPprint(stdout, symprior);
	}
	GCAPfree(contraprior);
	if(symprior->max_labels > gcasym->max_label) gcasym->max_label = symprior->max_labels;
      } // s
    } // r
  } //c 

  return(gcasym);
}

/*!
\fn int GCAisNotSymmetric(GCA *gca)
\brief Tests whether GCA is NOT symmetric. Returns 0 if it is symmetric.
 */
int GCAisNotSymmetric(GCA *gca)
{
  int n,m, c,r,s,csym, ok, k, kk, q,v ;
  GCA_NODE  *node, *symnode;
  GC1D *gcs, *symgcs;
  GCA_PRIOR *prior, *symprior;
  int labelcontra ;

  printf("Testing whether GCA is symmetric\n");
  printf("node  dim %d\n",gca->node_width);
  printf("prior dim %d\n",gca->prior_width);

  printf("Testing nodes\n");
  for(c=1; c < gca->node_width; c++){
    csym = gca->node_width - c - 1;
    for(r=1; r < gca->node_height; r++){
      for(s=1; s < gca->node_depth; s++){
	node = &(gca->nodes[c][r][s]);
	symnode = &(gca->nodes[csym][r][s]);

	// Test that the number of labels are the symmetric at this node
	if(node->nlabels != symnode->nlabels){
	  printf("Node (%d/%d,%d,%d) not sym nlabels %d %d\n",c,csym,r,s,node->nlabels,symnode->nlabels);
	  return(1);
	}

	for(n=0; n < node->nlabels; n++){
	  // Test that label identity is symmetric. The order might not be the same
	  ok = 0;
	  for(m=0; m < symnode->nlabels; m++){
	    labelcontra = MRIasegContraLatLabel(symnode->labels[m]);
	    if(node->labels[n] == labelcontra){
	      ok = 1;
	      break;
	    }
	  }
	  if(!ok){
	    printf("Node (%d/%d,%d,%d) cannot find %dth label %d in contra\n",c,csym,r,s,n,node->labels[n]);
	    return(2);
	  }

	  gcs = &node->gcs[n];
	  symgcs = &symnode->gcs[m];

	  // Test the means and covars
	  v = 0;
	  for(k = 0; k < gca->ninputs; k++){
	    if(gcs->means[k] != symgcs->means[k]){
	      printf("Node (%d/%d,%d,%d) not sym means %d   %lf %lf\n",c,csym,r,s,n,gcs->means[k],symgcs->means[k]);
	      return(3);
	    }
	    for(m=k; m < gca->ninputs; m++){
	      if(gcs->covars[v] != symgcs->covars[v]){
		printf("Node (%d,%d,%d) not sym covars %d   %lf %lf\n",c,r,s,v,gcs->covars[v],symgcs->covars[v]);
		return(4);
	      }
	      v++;
	    } // m
	  } // k

	  // Test MRF neighbors
	  for(k=0; k < GIBBS_NEIGHBORS; k++){
	    kk = k;
	    if(k==0) kk = 1;
	    if(k==1) kk = 0;
	    if(gcs->nlabels[k] != symgcs->nlabels[kk]){
	      printf("Node (%d/%d,%d,%d), label %d, MRF nbr %d,  not sym in number of labels %d %d\n",
		     c,csym,r,s,n,k,gcs->nlabels[k],symgcs->nlabels[kk]);
	      return(20);
	    }
	    // Now go through each label of this neighbor to find a match
	    ok = 0;
	    for(m=0; m < gcs->nlabels[k]; m++){
	      labelcontra = MRIasegContraLatLabel(gcs->labels[k][m]);
	      for(q=0; q < symgcs->nlabels[kk]; q++){
		if(symgcs->labels[kk][q] == labelcontra){
		  ok = 1;
		  break;
		}
	      }
	      if(!ok){
		printf("Node (%d/%d,%d,%d), label %d, MRF nbr %d, cannot find a match for label %d/%d\n",
		       c,csym,r,s,n,k,gcs->labels[k][m],labelcontra);
		printf("gca -------------------------------\n");
		GC1Dprint(stdout,gcs,gca->ninputs);
		printf("gca sym -------------------------------\n");
		GC1Dprint(stdout,symgcs,gca->ninputs);
		return(20);
	      }
	      // There is a match, check the prior probs
	      if(gcs->label_priors[k][m] != symgcs->label_priors[kk][q]){
		printf("Node (%d/%d,%d,%d), label %d, MRF nbr %d, probs dont match %f %f\n",
		       c,csym,r,s,n,k,gcs->label_priors[k][m],symgcs->label_priors[kk][q]);
		printf("gca -------------------------------\n");
		GC1Dprint(stdout,gcs,gca->ninputs);
		printf("gca sym -------------------------------\n");
		GC1Dprint(stdout,symgcs,gca->ninputs);
		return(21);
	      }
	    } // mth label of MRF neighbor
	  } // kth MRF neighbor
	} // nth label in the node
      } // s
    } // r
  } //c 

  printf("Testing priors\n");
  for(c=1; c < gca->prior_width; c++){
    csym = gca->prior_width - c - 1;
    for(r=1; r < gca->prior_height; r++){
      for(s=1; s < gca->prior_depth; s++){
	prior = &(gca->priors[c][r][s]);
	symprior = &(gca->priors[csym][r][s]);
	// Test whether the number of labels at this point is symmetric
	//if(prior->nlabels < 2 && symprior->nlabels < 2) continue;// not sure why this needed,edge
	if(prior->nlabels != symprior->nlabels){
	  printf("Prior ((%d,%d),%d,%d) not sym nlabels %d %d\n",c,csym,r,s,prior->nlabels,symprior->nlabels);
	  return(11);
	}
	// Test whether each label exists in the symmetric location
	for(n=0; n < prior->nlabels; n++){
	  ok = 0;
	  for(m=0; m < prior->nlabels; m++){
	    labelcontra = MRIasegContraLatLabel(symprior->labels[m]);
	    if(prior->labels[n] == labelcontra){
	      ok = 1;
	      break;
	    }
	  }
	  if(! ok){
	    printf("Prior (%d/%d,%d,%d) cannot find label %d in sym\n",c,csym,r,s,prior->labels[n]);
	    for(m=0; m < prior->nlabels; m++) printf(" %d",prior->labels[m]);
	    printf("\n");
	    for(m=0; m < prior->nlabels; m++) printf(" %d",symprior->labels[m]);
	    printf("\n");
	    return(12);
	  }
	  // Test whether the probability of this label at this point is symmetric
	  if(prior->priors[n] != symprior->priors[m]){
	    printf("Prior (%d/%d,%d,%d) not sym pval n=%d m=%d   %lf %lf\n",
		   c,csym,r,s,n,m,prior->priors[n],symprior->priors[n]);
	    return(13);
	  }
	}
      } // s
    } // r
  } //c 

  return(0);
}
