/**
 * @brief utilities for MRI classification using a variety of classifiers
 *
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diag.h"
#include "error.h"
#include "gclass.h"
#include "macros.h"
#include "mri.h"
#include "mriclass.h"
#include "proto.h"
#include "region.h"
#include "utils.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define POLV_WSIZE 5
#define DEBUG_POINT(x, y, z) (((x) == 35) && ((y) == 4) && ((z) == 31))

#define MAX_FILES 100

/*-----------------------------------------------------
                       STRUCTURES
-------------------------------------------------------*/

typedef struct
{
  FILE *fp;
  int npixels[MAX_INPUTS];
  MRIC *mric;
  int round;
} GET_INPUT_PARMS;

/*-----------------------------------------------------
                      GLOBAL DATA
-------------------------------------------------------*/

const char *class_names[GAUSSIAN_NCLASSES] = {"CSF", "GREY MATTER", "THIN STRANDS", "BORDER PIXELS", "WHITE MATTER", "BRIGHT MATTER"};

/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/

static int total_calls = 0, buffered = 0, total_computed = 0;

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int mricComputeGCStatistics(MRIC *mric, FILE *fp, int nfiles, int round);
static int mricTrainRBF(MRIC *mric, FILE *fp, int nfiles, int round);
static int mricRetrainRBF(MRIC *mric, FILE *fp, int nfiles, int round);
static int mricGetClassifierInput(VECTOR *v_inputs, int no, void *parm, int same_class, int *pclass);

static int mricGetClassObservation(GET_INPUT_PARMS *parms, VECTOR *v_obs, int obs_no, int classnum);
static int mricFillInputParms(MRIC *mric, FILE *fp, int round, GET_INPUT_PARMS *parms);

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Free an MRIC structure and all its members.
------------------------------------------------------*/
int MRICfree(MRIC **pmric)
{
  MRIC *mric;
  int round;

  mric = *pmric;
  *pmric = NULL;
  if (mric) {
    for (round = 0; round < mric->nrounds; round++) {
      switch (mric->type[round]) {
        case CLASSIFIER_RBF:
          RBFfree(&mric->classifier[round].rbf);
          break;
        case CLASSIFIER_GAUSSIAN:
          GCfree(&mric->classifier[round].gc);
          break;
        default:
          break;
      }
    }
    free(mric);
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIC *MRICalloc(int nrounds, int *types, int *features, void *parms)
{
  MRIC *mric;
  unsigned f;
  int ninputs, round;

  mric = (MRIC *)calloc(1, sizeof(MRIC));
  if (!mric) ErrorExit(ERROR_NO_MEMORY, "MRICalloc(%d): could not allocate struct", nrounds);

  mric->nrounds = nrounds;
  for (round = 0; round < nrounds; round++) {
    for (ninputs = 0, f = 0x001; f != MAX_FEATURE; f <<= 1)
      if (f & features[round]) ninputs++;

    if (ninputs < 1 || ninputs > MAX_INPUTS)
      ErrorReturn(NULL, (ERROR_BADPARM, "MRICalloc(%d): bad # of inputs %d", round, ninputs));

    mric->type[round] = types[round];
    mric->ninputs[round] = ninputs;
    mric->features[round] = features[round];
    switch (types[round]) {
      case CLASSIFIER_RBF: {
        RBF_PARMS *rbf_parms;

        if (parms) {
          rbf_parms = (RBF_PARMS *)parms;
          mric->classifier[round].rbf = RBFinit(ninputs, NCLASSES, rbf_parms->max_clusters, class_names);
        }
        break;
      }
      case CLASSIFIER_GAUSSIAN:
        mric->classifier[round].gc = GCalloc(GAUSSIAN_NCLASSES, ninputs, class_names);
        if (!mric->classifier[round].gc) {
          free(mric);
          return (NULL);
        }
        break;
      default:
        MRICfree(&mric);
        ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRICalloc: classifier %d not supported", types[round]));
        break;
    }
  }
  return (mric);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRICretrain(MRIC *mric, char *file_name)
{
  char line[300], *cp;
  FILE *fp;
  int nfiles, round;

  /* first figure out the total # of files */
  fp = fopen(file_name, "r");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRICtrain(%s): could not open file", file_name));

  nfiles = 0;
  while ((cp = fgetl(line, 299, fp)) != NULL) nfiles++;
  fprintf(stderr, "processing %d files\n", nfiles);

  for (round = 0; round < mric->nrounds; round++) {
    rewind(fp);
    switch (mric->type[round]) {
      case CLASSIFIER_RBF:
        MRICsetRegionSize(mric, 3, 3, 3); /* inputs are at random locations */
        mricRetrainRBF(mric, fp, nfiles, round);
        break;
      case CLASSIFIER_GAUSSIAN:
        break;
    }
  }
  fclose(fp);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRICtrain(MRIC *mric, char *file_name, char *prior_fname)
{
  char line[300], *cp;
  FILE *fp;
  int nfiles, round;

  if (prior_fname) {
    FileNameAbsolute(prior_fname, mric->prior_fname);
    mric->mri_priors = MRIread(prior_fname);
    if (!mric->mri_priors)
      ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRICtrain: could not load prior file '%s'", prior_fname));
  }
  else
    mric->mri_priors = NULL;

  /* first figure out the total # of files */
  fp = fopen(file_name, "r");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRICtrain(%s): could not open file", file_name));

  nfiles = 0;
  while ((cp = fgetl(line, 299, fp)) != NULL) nfiles++;
  fprintf(stderr, "processing %d files\n", nfiles);

  for (round = 0; round < mric->nrounds; round++) {
    rewind(fp);
    switch (mric->type[round]) {
      case CLASSIFIER_RBF:
        if (mricTrainRBF(mric, fp, nfiles, round) != NO_ERROR) return (Gerror);
        break;
      case CLASSIFIER_GAUSSIAN:
        mricComputeGCStatistics(mric, fp, nfiles, round);
        break;
    }
  }
  fclose(fp);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIC *MRICread(char *fname)
{
  MRIC *mric;
  int ninputs, type[MAX_ROUNDS], features[MAX_ROUNDS], round, nrounds;
  FILE *fp;
  char prior_fname[100], line[100], *cp;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NO_FILE, "MRICread(%s): could not open file", fname));

  prior_fname[0] = 0;
  cp = fgetl(line, 99, fp);
  if (sscanf(cp, "%d %s", &nrounds, prior_fname) < 1) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_BADFILE, "MRICread(%s): could not scan parms", fname));
  }
  for (round = 0; round < nrounds; round++) {
    cp = fgetl(line, 99, fp);
    if (sscanf(cp, "%d %d 0x%x", &type[round], &ninputs, &features[round]) < 3) {
      fclose(fp);
      ErrorReturn(NULL, (ERROR_BADFILE, "MRICread(%s): could not scan parms", fname));
    }
  }

  mric = MRICalloc(nrounds, type, features, NULL);
  if (*prior_fname) {
    mric->mri_priors = MRIread(prior_fname);
    strcpy(mric->prior_fname, prior_fname);
  }
  for (round = 0; round < nrounds; round++) {
    switch (type[round]) {
      case CLASSIFIER_RBF:
        mric->classifier[round].rbf = RBFreadFrom(fp);
        break;
      case CLASSIFIER_GAUSSIAN:
        GCasciiReadFrom(fp, mric->classifier[round].gc);
        break;
      default:
        break;
    }
  }

  fclose(fp);

  return (mric);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRIC *MRICquickRead(char *fname)
{
  MRIC *mric;
  int ninputs, type[MAX_ROUNDS], features[MAX_ROUNDS], round, nrounds;
  FILE *fp;
  char line[100], *cp, prior_fname[100];

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NO_FILE, "MRICread(%s): could not open file", fname));

  prior_fname[0] = 0;
  cp = fgetl(line, 99, fp);
  if (sscanf(cp, "%d %s", &nrounds, prior_fname) < 1) {
    fclose(fp);
    ErrorReturn(NULL, (ERROR_BADFILE, "MRICread(%s): could not scan parms", fname));
  }
  for (round = 0; round < nrounds; round++) {
    cp = fgetl(line, 99, fp);
    if (sscanf(cp, "%d %d 0x%x", &type[round], &ninputs, &features[round]) < 3) {
      fclose(fp);
      ErrorReturn(NULL, (ERROR_BADFILE, "MRICread(%s): could not scan parms", fname));
    }
  }

  mric = MRICalloc(nrounds, type, features, NULL);
  if (prior_fname[0]) strcpy(mric->prior_fname, prior_fname);

  for (round = 0; round < nrounds; round++) {
    switch (type[round]) {
      case CLASSIFIER_GAUSSIAN:
        GCasciiReadFrom(fp, mric->classifier[round].gc);
        break;
      default:
        break;
    }
  }

  fclose(fp);

  return (mric);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int MRICwrite(MRIC *mric, char *fname)
{
  FILE *fp;
  int round;

  fp = fopen(fname, "wb");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRICwrite(%s): could not open file", fname));

  fprintf(fp, "%d ", mric->nrounds);
  if (mric->mri_priors) fprintf(fp, " %s", mric->prior_fname);
  fprintf(fp, "\n");

  for (round = 0; round < mric->nrounds; round++)
    fprintf(fp, "%d %d 0x%x\n", mric->type[round], mric->ninputs[round], mric->features[round]);

  for (round = 0; round < mric->nrounds; round++) {
    switch (mric->type[round]) {
      case CLASSIFIER_RBF:
        RBFwriteInto(mric->classifier[round].rbf, fp);
        break;
      case CLASSIFIER_GAUSSIAN:
        GCasciiWriteInto(fp, mric->classifier[round].gc);
        break;
      default:
        break;
    }
  }
  fclose(fp);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
#define PRETTY_SURE .90f
#define CLASS_SCALE ((BUFTYPE)1)

MRI *MRICclassify(MRIC *mric, MRI *mri_src, MRI *mri_dst, float conf, MRI *mri_probs, MRI *mri_classes)
{
  MATRIX *m_priors;
  VECTOR *v_inputs;
  GCLASSIFY *gc;
  int x, y, z, classno, nclasses, xt, yt, zt, type, round, x1, y1, z1, x0, c, bg, max_white_class, max_non_white_class;
  // int width, depth, height;
  BUFTYPE *psrc, src, *pdst, *pclasses;
  float prob, *pprobs = NULL, total, min_output, fmin, fmax, white_prob, non_white_prob, max_white_prob,
              max_non_white_prob, p;
  double xrt, yrt, zrt;
  MRI *mri_priors, *mri_in;
  MRI_REGION bounding_box;
  RBF *rbf;

  max_white_class = max_non_white_class = -1;
#if 0
  MRIboundingBox(mri_src, DEFINITELY_BACKGROUND, &bounding_box) ;
#else
  MRIvalRange(mri_src, &fmin, &fmax);
  bg = nint(.2 * (fmax - fmin) + fmin);
  MRIboundingBox(mri_src, bg, &bounding_box);
#endif
  REGIONexpand(&bounding_box, &bounding_box, POLV_WSIZE / 2);
  MRIclipRegion(mri_src, &bounding_box, &bounding_box);
  REGIONcopy(&bounding_box, &mri_src->roi);
  x0 = bounding_box.x;
  x1 = bounding_box.x + bounding_box.dx - 1;
  y1 = bounding_box.y + bounding_box.dy - 1;
  z1 = bounding_box.z + bounding_box.dz - 1;
  if (conf < 0.0f || conf >= 1.0f) conf = PRETTY_SURE;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL);

  mri_priors = mric->mri_priors;
  if (mri_priors)
    m_priors = MatrixAlloc(mric->classifier[0].gc->nclasses, 1, MATRIX_REAL);
  else
    m_priors = NULL;

  // width = mri_src->width;
  // height = mri_src->height;
  // depth = mri_src->depth;

  mri_in = mri_src;
  for (round = 0; round < mric->nrounds; round++) {
    type = mric->type[round];
    gc = mric->classifier[round].gc;
    rbf = mric->classifier[round].rbf;
    nclasses = gc->nclasses;
    v_inputs = VectorAlloc(mric->ninputs[round], MATRIX_REAL);

    for (z = bounding_box.z; z <= z1; z++) {
      DiagHeartbeat((float)((z - bounding_box.z) + round * bounding_box.dz) / (float)(bounding_box.dz * mric->nrounds));
      for (y = bounding_box.y; y <= y1; y++) {
        psrc = &MRIvox(mri_src, x0, y, z);
        pdst = &MRIvox(mri_dst, x0, y, z);
        if (mri_probs)
          pprobs = &MRIFvox(mri_probs, x0, y, z);
        else
          pprobs = NULL;
        if (mri_classes)
          pclasses = &MRIvox(mri_classes, x0, y, z);
        else
          pclasses = NULL;
        for (x = x0; x <= x1; x++) {
          src = *psrc++;
          if (mri_priors) {
            MRIvoxelToVoxel(mri_src, mri_priors, (double)x, (double)y, (double)z, &xrt, &yrt, &zrt);
            xt = mri_priors->xi[nint(xrt)];
            yt = mri_priors->yi[nint(yrt)];
            zt = mri_priors->zi[nint(zrt)];
            for (classno = 0; classno < nclasses; classno++)
              m_priors->rptr[classno + 1][1] = MRIFseq_vox(mri_priors, xt, yt, zt, classno);
          }
          MRICcomputeInputs(mric, mri_in, x, y, z, v_inputs, mric->features[round]);

          switch (type) /* now classify this observation */
          {
            default:
              ErrorReturn(NULL, (ERROR_UNSUPPORTED, "MRICclassify: unsupported classifier type %d", type));
            case CLASSIFIER_GAUSSIAN:
              classno = GCclassify(gc, v_inputs, m_priors, &prob);
              break;
            case CLASSIFIER_RBF:
              classno = RBFclassify(rbf, v_inputs);
              min_output = RVECTOR_ELT(rbf->v_outputs, 1);
              for (c = 2; c <= rbf->noutputs; c++)
                if (RVECTOR_ELT(rbf->v_outputs, c) < min_output) min_output = RVECTOR_ELT(rbf->v_outputs, c);
              for (total = 0.0f, c = 1; c <= rbf->noutputs; c++) total += RVECTOR_ELT(rbf->v_outputs, c) - min_output;

              /* must sum all probs of white and not white */
              max_non_white_prob = max_white_prob = -1.0f;
              non_white_prob = white_prob = 0.0f;
#if 0
            if (DEBUG_POINT(x,y,z))
              DiagBreak() ;
#endif

              for (c = 0; c < rbf->noutputs; c++) {
                p = RVECTOR_ELT(rbf->v_outputs, c + 1);
                if (ISWHITE(c)) {
                  if (p > max_white_prob) {
                    max_white_prob = p;
                    max_white_class = c;
                  }
                  white_prob += p;
                }
                else {
                  if (p > max_non_white_prob) {
                    max_non_white_prob = p;
                    max_non_white_class = c;
                  }
                  non_white_prob += p;
                }
              }
              if (white_prob > non_white_prob) {
                classno = max_white_class;
                prob = white_prob;
              }
              else {
                classno = max_non_white_class;
                prob = non_white_prob;
              }

              if (!FZERO(total))
                prob = (RVECTOR_ELT(rbf->v_outputs, classno + 1) - min_output) / total;
              else
                prob = 0.0f;
              break;
          }

          if (pclasses) *pclasses++ = (BUFTYPE)classno * CLASS_SCALE;
          if (pprobs) *pprobs++ = prob;
          if (ISWHITE(classno) && (prob > conf)) {
            if (!src && mric->features[round] == FEATURE_CPOLV_MEDIAN5)
              *pdst++ = 255;
            else
              *pdst++ = src;
          }
          else
            *pdst++ = 0 /*src*/;
        }
      }
    }
    VectorFree(&v_inputs);
    mri_in = mri_classes; /* all subsequent rounds */
  }

  if (Gdiag & DIAG_SHOW) {
    fprintf(stderr, "total = %d, buffered = %d, computed = %d\n", total_calls, buffered, total_computed);
    fprintf(stderr,
            "efficiency = %2.3f%%\n",
            100.0f * (float)(mri_src->width * mri_src->height * mri_src->depth) / (float)total_computed);
  }
  if (m_priors) MatrixFree(&m_priors);

  return (mri_dst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
#define DEFAULT_REGION_SIZE 8
static int region_size = DEFAULT_REGION_SIZE;
static int region_height = 0;
static int region_width = 0;

int MRICsetRegionSize(MRIC *mric, int rwidth, int rheight, int rdepth)
{
  int old_rsize = region_size;
  region_size = rdepth;
  region_height = rheight;
  region_width = rwidth;
  return (old_rsize);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRICresetRegionSize(MRIC *mric)
{
  int old_rsize = region_size;
  region_size = DEFAULT_REGION_SIZE;
  region_height = region_width = 0;
  return (old_rsize);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRICcomputeInputs(MRIC *mric, MRI *mri, int x, int y, int z, VECTOR *v_inputs, int features)
{
  static int old_features = 0;
  static MRI_REGION region = {0, 0, 0, 0, 0, 0};
  static MRI *mri_prev = NULL, *mri_zscore3 = NULL, *mri_zscore5 = NULL, *mri_direction = NULL, *mri_mean3 = NULL,
             *mri_mean5 = NULL, *mri_cpolv = NULL, *mri_cpolv_mean3 = NULL, *mri_cpolv_mean5 = NULL,
             *mri_cpolv_median3 = NULL, *mri_cpolv_median5 = NULL, *mri_min3 = NULL, *mri_min5 = NULL, *mri_min7 = NULL,
             *mri_cpolv_zscore5 = NULL, *mri_cpolv_zscore7 = NULL, *mri_cpolv_curv5 = NULL, *mri_cpolv_curv7 = NULL,
             *mri_cpolv_order = NULL;

  int x0, y0, z0, index;
  char *cp;

  if (mri != mri_prev) /* reset counters */
    total_calls = buffered = total_computed = 0;

  if ((features & FEATURE_CPOLV) && (!mri_cpolv || !MRImatch(mri, mri_cpolv) || (mri != mri_prev))) {
    if (mri_cpolv) MRIfree(&mri_cpolv);
    if (Gdiag & DIAG_HEARTBEAT) fprintf(stderr, "computing cpolv...\n");
    cp = getenv("cpolv");
    if (cp && !mri_cpolv) {
      mri_cpolv = MRIread(cp);
      if (!MRImatch(mri, mri_cpolv)) MRIfree(&mri_cpolv);
    }

    if (!mri_cpolv) mri_cpolv = MRIcentralPlaneOfLeastVarianceNormal(mri, NULL, POLV_WSIZE);
    if (Gdiag & DIAG_WRITE) MRIwrite(mri_cpolv, "cpolv.mnc");
#if 0
    if (features & FEATURE_CPOLV_MEDIAN5)
    {
      if (mri_cpolv_median5)
        MRIfree(&mri_cpolv_median5) ;
      mri_cpolv_median5 = MRIpolvMedian(mri, NULL, mri_cpolv, 5) ;
    }
#endif
  }
  total_calls++;
  /*
     if the specified point is outside of the precomputed window,
     Update the window and compute a new set of input images.
     */
  if (!mri_prev || mri_prev != mri || (old_features != features) ||
      (REGIONinside(&region, x, y, z) == REGION_OUTSIDE)) {
    MRI *mri_std, *mri_mean, *mri_region, *mri_grad, *mri_big = NULL;
    MRI_REGION rbig;
    int x0, y0, z0; /* starting coords in mri_region */

    old_features = features;

    /* compute clip region with border */
    region.x = x;
    region.y = y;
    region.z = z;
    region.dx = region_width > 0 ? region_width : mri->width;
    region.dy = region_height > 0 ? region_height : mri->height;
    region.dz = region_size;
    MRIclipRegion(mri, &region, &region);
    total_computed += (region.dx * region.dy * region.dz);

    REGIONexpand(&region, &rbig, 2);
    MRIclipRegion(mri, &rbig, &rbig);
    mri_region = MRIextractRegion(mri, NULL, &rbig);
    x0 = region.x - rbig.x;
    y0 = region.y - rbig.y;
    z0 = region.z - rbig.z;

    mri_prev = mri;
    if (mri_cpolv_curv5) MRIfree(&mri_cpolv_curv5);
    if (mri_cpolv_curv7) MRIfree(&mri_cpolv_curv7);
    if (mri_cpolv_zscore5) MRIfree(&mri_cpolv_zscore5);
    if (mri_cpolv_order) MRIfree(&mri_cpolv_order);
    if (mri_cpolv_zscore7) MRIfree(&mri_cpolv_zscore7);

    if (mri_min3) MRIfree(&mri_min3);
    if (mri_min5) MRIfree(&mri_min5);
    if (mri_min7) MRIfree(&mri_min7);
    if (mri_zscore3) MRIfree(&mri_zscore3);
    if (mri_zscore5) MRIfree(&mri_zscore5);
    if (mri_mean3) MRIfree(&mri_mean3);
    if (mri_mean5) MRIfree(&mri_mean5);
    if (mri_direction) MRIfree(&mri_direction);
    if (mri_cpolv_mean3) MRIfree(&mri_cpolv_mean3);
    if (mri_cpolv_mean5) MRIfree(&mri_cpolv_mean5);
    if (mri_cpolv_median3) MRIfree(&mri_cpolv_median3);
#if 1
    if (mri_cpolv_median5) MRIfree(&mri_cpolv_median5);
#endif

    if (features & FEATURE_CPOLV_ORDER) {
      static int first = 1;
#define ORDER_THRESH 2
      mri_big = MRIpolvOrder(mri_region, NULL, mri_cpolv, 5, ORDER_THRESH);
      mri_cpolv_order = MRIextract(mri_big, NULL, x0, y0, z0, region.dx, region.dy, region.dz);
      if (first && (Gdiag & DIAG_WRITE)) {
        MRIwrite(mri_region, "region.mnc");
        MRIwrite(mri_big, "big.mnc");
        MRIwrite(mri_cpolv_order, "order.mnc");
        first = 0;
      }

      MRIfree(&mri_big);
    }

    if (features & FEATURE_CPOLV_ZSCORE5) {
      mri_big = MRIpolvZscore(mri_region, NULL, mri_cpolv, 5);
      mri_cpolv_zscore5 = MRIextract(mri_big, NULL, x0, y0, z0, region.dx, region.dy, region.dz);
      MRIfree(&mri_big);
    }

    if (features & FEATURE_CPOLV_ZSCORE7) {
      mri_big = MRIpolvZscore(mri_region, NULL, mri_cpolv, 7);
      mri_cpolv_zscore7 = MRIextract(mri_big, NULL, x0, y0, z0, region.dx, region.dy, region.dz);
      MRIfree(&mri_big);
    }

    if (features & FEATURE_CPOLV_CURV5) {
      mri_big = MRIpolvNormalCurvature(mri_region, NULL, mri_cpolv, 5);
      mri_cpolv_curv5 = MRIextract(mri_big, NULL, x0, y0, z0, region.dx, region.dy, region.dz);
      MRIfree(&mri_big);
    }

    if (features & FEATURE_CPOLV_CURV7) {
      mri_big = MRIpolvNormalCurvature(mri_region, NULL, mri_cpolv, 7);
      mri_cpolv_curv7 = MRIextract(mri_big, NULL, x0, y0, z0, region.dx, region.dy, region.dz);
      MRIfree(&mri_big);
    }

    if (features & FEATURE_ZSCORE3) {
      static int first = 1;
      mri_mean = MRImeanRegion(mri, NULL, 3, &region);
      mri_std = MRIstdRegion(mri, NULL, mri_mean, 3, &region);
      mri_zscore3 = MRIzScoreRegion(mri, NULL, mri_mean, mri_std, &region);
      if (Gdiag & DIAG_WRITE && first) {
        first = 0;
        MRIwrite(mri_mean, "mean.mnc");
        MRIwrite(mri_std, "std.mnc");
        MRIwrite(mri_zscore3, "zscore.mnc");
      }
      MRIfree(&mri_mean);
      MRIfree(&mri_std);
    }
    if (features & FEATURE_ZSCORE5) {
      static int first = 1;
      mri_mean = MRImeanRegion(mri, NULL, 5, &region);
      mri_std = MRIstdRegion(mri, NULL, mri_mean, 5, &region);
      mri_zscore5 = MRIzScoreRegion(mri, NULL, mri_mean, mri_std, &region);
      if (Gdiag & DIAG_WRITE && first) {
        first = 0;
        MRIwrite(mri_mean, "mean.mnc");
        MRIwrite(mri_std, "std.mnc");
        MRIwrite(mri_zscore5, "zscore.mnc");
      }
      MRIfree(&mri_mean);
      MRIfree(&mri_std);
    }
    if (features & FEATURE_DIRECTION) {
      static int first = 1;

      /* expand region by 2 voxels to avoid border effects */
      mri_grad = MRIsobel(mri_region, NULL, NULL);
      mri_big = MRIdirectionMap(mri_grad, NULL, 3);

      mri_direction = MRIextract(mri_big, NULL, x0, y0, z0, region.dx, region.dy, region.dz);
      if ((Gdiag & DIAG_WRITE) && first) {
        first = 0;
        MRIwrite(mri_direction, "dir.mnc");
        MRIwrite(mri_region, "region.mnc");
        MRIwrite(mri_grad, "grad.mnc");
        MRIwrite(mri_big, "tmp.mnc");
      }
      MRIfree(&mri_grad);
    }
    MRIfree(&mri_region);
    if (features & FEATURE_MEAN3) mri_mean3 = MRImeanRegion(mri, NULL, 3, &region);
    if (features & FEATURE_MEAN5) mri_mean5 = MRImeanRegion(mri, NULL, 5, &region);
    if (features & FEATURE_CPOLV_MEAN3) mri_cpolv_mean3 = MRIpolvMeanRegion(mri, NULL, mri_cpolv, 3, &region);
    if (features & FEATURE_CPOLV_MEAN5) mri_cpolv_mean5 = MRIpolvMeanRegion(mri, NULL, mri_cpolv, 5, &region);
    if (features & FEATURE_CPOLV_MEDIAN3) mri_cpolv_median3 = MRIpolvMedianRegion(mri, NULL, mri_cpolv, 3, &region);
#if 1
    if (features & FEATURE_CPOLV_MEDIAN5) {
      static int first = 1;

      mri_cpolv_median5 = MRIpolvMedianRegion(mri, NULL, mri_cpolv, 5, &region);
      if ((Gdiag & DIAG_WRITE) && first) {
        first = 0;
        MRIwrite(mri_cpolv_median5, "median5.mnc");
      }
    }
#endif
    if (features & FEATURE_MIN3) mri_min3 = MRIerodeRegion(mri, NULL, 3, &region);
    if (features & FEATURE_MIN5) mri_min5 = MRIerodeRegion(mri, NULL, 5, &region);
    if (features & FEATURE_MIN7) mri_min7 = MRIerodeRegion(mri, NULL, 7, &region);
    if (features & FEATURE_MIN5) {
      static int first = 1;

      if ((Gdiag & DIAG_WRITE) && first) {
        first = 0;
        MRIwrite(mri_min5, "min5.mnc");
      }
    }
  }
  else
    buffered++;

  /* x0,y0,z0 are coordinates in region based images (not input mri or polv) */
  x0 = x - region.x;
  y0 = y - region.y;
  z0 = z - region.z;
  index = 1; /* inputs are 1-based because of NRC */
  if (features & FEATURE_INTENSITY) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri, x, y, z);
  if (features & FEATURE_ZSCORE3) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_zscore3, x0, y0, z0);
  if (features & FEATURE_ZSCORE5) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_zscore5, x0, y0, z0);
  if (features & FEATURE_DIRECTION) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_direction, x0, y0, z0);
  if (features & FEATURE_MEAN3) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_mean3, x0, y0, z0);
  if (features & FEATURE_MEAN5) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_mean5, x0, y0, z0);
  if (features & FEATURE_CPOLV_MEAN3) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri_cpolv_mean3, x0, y0, z0);
  if (features & FEATURE_CPOLV_MEAN5) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri_cpolv_mean5, x0, y0, z0);
  if (features & FEATURE_CPOLV_MEDIAN3) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri_cpolv_median3, x0, y0, z0);
#if 1
  if (features & FEATURE_CPOLV_MEDIAN5) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri_cpolv_median5, x0, y0, z0);
#else
  if (features & FEATURE_CPOLV_MEDIAN5) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri_cpolv_median5, x, y, z);
#endif
  if (features & FEATURE_MIN3) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri_min3, x0, y0, z0);
  if (features & FEATURE_MIN5) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri_min5, x0, y0, z0);
  if (features & FEATURE_MIN7) VECTOR_ELT(v_inputs, index++) = (float)MRIvox(mri_min7, x0, y0, z0);
  if (features & FEATURE_X_POSITION) VECTOR_ELT(v_inputs, index++) = (float)x;
  if (features & FEATURE_Y_POSITION) VECTOR_ELT(v_inputs, index++) = (float)y;
  if (features & FEATURE_Z_POSITION) VECTOR_ELT(v_inputs, index++) = (float)z;
  if ((features & FEATURE_PRIORS) && mric->mri_priors) {
    double xrt, yrt, zrt;
    int xt, yt, zt;
    MRI *mri_priors;

    mri_priors = mric->mri_priors;
    MRIvoxelToVoxel(mri, mri_priors, (double)x, (double)y, (double)z, &xrt, &yrt, &zrt);
    xt = mri_priors->xi[nint(xrt)];
    yt = mri_priors->yi[nint(yrt)];
    zt = mri_priors->zi[nint(zrt)];
    VECTOR_ELT(v_inputs, index++) = MRIFseq_vox(mric->mri_priors, xt, yt, zt, WHITE_MATTER);
  }
  if (features & FEATURE_CPOLV_ZSCORE5) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_cpolv_zscore5, x0, y0, z0);
  if (features & FEATURE_CPOLV_ZSCORE7) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_cpolv_zscore7, x0, y0, z0);
  if (features & FEATURE_CPOLV_CURV5) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_cpolv_curv5, x0, y0, z0);
  if (features & FEATURE_CPOLV_CURV7) VECTOR_ELT(v_inputs, index++) = MRIFvox(mri_cpolv_curv7, x0, y0, z0);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Build an image of target values. Each point in the image
          contains the class # based on intensity thresholds either
          provided by the user, or defaults if parms are 0.
------------------------------------------------------*/
MRI *MRICbuildTargetImage(MRI *mri_src, MRI *mri_target, MRI *mri_wm, int lo_lim, int hi_lim)
{
  int x, y, z, width, height, depth;
  BUFTYPE *psrc, *pwm, *ptarget, src, wm, target, *pthin, thin;
  MRI *mri_thin;

  mri_thin = MRIfindThinWMStrands(mri_wm, NULL, 2);
  if (lo_lim <= 0) lo_lim = LO_LIM;
  if (hi_lim <= 0) hi_lim = HI_LIM;

  if (!mri_target) mri_target = MRIclone(mri_wm, NULL);

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      psrc = &MRIvox(mri_src, 0, y, z);
      ptarget = &MRIvox(mri_target, 0, y, z);
      pwm = &MRIvox(mri_wm, 0, y, z);
      pthin = &MRIvox(mri_thin, 0, y, z);
      for (x = 0; x < width; x++) {
        src = *psrc++;
        wm = *pwm++;
        thin = *pthin++;
        if (thin)
          target = THIN_STRAND;
        else if (wm) {
          int xi, yi, zi, xk, yk, zk;

          target = WHITE_MATTER;
          for (zk = z - 1; zk <= z + 1; zk++) {
            zi = mri_target->zi[zk];
            for (yk = y - 1; yk <= y + 1; yk++) {
              yi = mri_target->yi[yk];
              for (xk = x - 1; xk <= x + 1; xk++) {
                xi = mri_target->xi[xk];
                if (!MRIvox(mri_wm, xi, yi, zi)) target = BORDER_MATTER;
              }
            }
          }
        }
        else if (src > hi_lim)
          target = BRIGHT_MATTER;
        else if (src < lo_lim)
          target = BACKGROUND;
        else
          target = GREY_MATTER;
        *ptarget++ = target;
      }
    }
  }
  return (mri_target);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Update the prior estimates using a frequency count.
------------------------------------------------------*/
#define COUNT_IMAGE GAUSSIAN_NCLASSES

MRI *MRICupdatePriors(MRI *mri_target, MRI *mri_priors, int scale)
{
  int width, height, depth, x, y, z, classnum, w, h, d, xt, yt, zt;
  BUFTYPE *ptarget;
  double xrt, yrt, zrt;

  width = mri_target->width;
  height = mri_target->height;
  depth = mri_target->depth;
  w = (int)(((float)width / (float)scale) + 0.99f);
  h = (int)(((float)height / (float)scale) + 0.99f);
  d = (int)(((float)depth / (float)scale) + 0.99f);

  /*
     allocate one image for each class, plus one to keep track of the
     # of pixels mapped to that location.
     */
  if (!mri_priors) {
    mri_priors = MRIallocSequence(w, h, d, MRI_FLOAT, GAUSSIAN_NCLASSES + 1);
    MRIcopyHeader(mri_target, mri_priors);
    mri_priors->xsize *= scale;
    mri_priors->ysize *= scale;
    mri_priors->zsize *= scale;
    mri_priors->linear_transform = mri_priors->inverse_linear_transform = NULL;
  }

  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      ptarget = &MRIvox(mri_target, 0, y, z);
      for (x = 0; x < width; x++) {
        MRIvoxelToTalairachVoxel(mri_target, (double)x, (double)y, (double)z, &xrt, &yrt, &zrt);
        xt = mri_priors->xi[nint(xrt / scale)];
        yt = mri_priors->yi[nint(yrt / scale)];
        zt = mri_priors->zi[nint(zrt / scale)];
        classnum = *ptarget++;
        MRIFseq_vox(mri_priors, xt, yt, zt, classnum) = MRIFseq_vox(mri_priors, xt, yt, zt, classnum) + 1.0f;
        MRIFseq_vox(mri_priors, xt, yt, zt, COUNT_IMAGE) = MRIFseq_vox(mri_priors, xt, yt, zt, COUNT_IMAGE) + 1.0f;
      }
    }
  }
  return (mri_priors);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Normalize the priors by dividing by the total # of
           times each pixel has been mapped.
------------------------------------------------------*/
int MRInormalizePriors(MRI *mri_priors)
{
  int width, height, depth, x, y, z, classnum;
  float *pnorm, norm;

  width = mri_priors->width;
  height = mri_priors->height;
  depth = mri_priors->depth;

  /*
     allocate one image for each class, plus one to keep track of the
     # of pixels mapped to that location.
     */
  for (z = 0; z < depth; z++) {
    for (y = 0; y < height; y++) {
      pnorm = &MRIFseq_vox(mri_priors, 0, y, z, COUNT_IMAGE);
      for (x = 0; x < width; x++) {
        norm = *pnorm++;
        if (!FZERO(norm))
          for (classnum = 0; classnum < GAUSSIAN_NCLASSES; classnum ++) {
            MRIFseq_vox(mri_priors, x, y, z, classnum) /= norm;
          }
      }
    }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRICupdateStatistics(MRIC *mric, int round, MRI *mri_src, MRI *mri_wm, MRI_REGION *box)
{
  GCLASSIFY *gc;
  GCLASS *gcl;
  int x, y, z, classno, row, col, x1, y1, z1;
  BUFTYPE *ptarget;
  // BUFTYPE *psrc;
  // int width, height, depth, src, nclasses;
  float covariance;
  VECTOR *v_inputs;
  MRI *mri_target;

  mri_target = MRICbuildTargetImage(mri_src, NULL, mri_wm, LO_LIM, HI_LIM);
  v_inputs = VectorAlloc(mric->ninputs[round], MATRIX_REAL);

  // nclasses = mric->classifier[round].gc->nclasses;
  gc = mric->classifier[round].gc;

  // width = mri_src->width;
  // height = mri_src->height;
  // depth = mri_src->depth;

  x1 = box->x + box->dx - 1;
  y1 = box->y + box->dy - 1;
  z1 = box->z + box->dz - 1;
  for (z = box->z; z <= z1; z++) {
    DiagHeartbeat((float)((z - box->z) + round * box->dz) / (float)(box->dz * mric->nrounds));
    for (y = box->y; y <= y1; y++) {
      // psrc = &MRIvox(mri_src, box->x, y, z);
      ptarget = &MRIvox(mri_target, box->x, y, z);
      for (x = box->x; x <= x1; x++) {
        // psrc is not used for anything, so commented out /clarsen
        // src = *psrc++;
        classno = (int)*ptarget++;
        gcl = &gc->classes[classno];
        gcl->nobs++;

        MRICcomputeInputs(mric, mri_src, x, y, z, v_inputs, mric->features[round]);
        for (row = 1; row <= mric->ninputs[round]; row++) gcl->m_u->rptr[row][1] += VECTOR_ELT(v_inputs, row);
        for (row = 1; row <= gcl->m_covariance->rows; row++) {
          for (col = 1; col <= row; col++) {
            covariance = gcl->m_covariance->rptr[row][col] + VECTOR_ELT(v_inputs, row) * VECTOR_ELT(v_inputs, col);
            gcl->m_covariance->rptr[row][col] = covariance;
            /*            gcl->m_covariance->rptr[col][row] = covariance;*/
          }
        }
      }
    }
  }

  VectorFree(&v_inputs);
  MRIfree(&mri_target);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRICcomputeStatistics(MRIC *mric, int round)
{
  GCLASSIFY *gc;
  GCLASS *gcl;
  int classno, nclasses, nobs, row, col;
  float mean_a, mean_b, covariance;

  gc = mric->classifier[round].gc;
  nclasses = gc->nclasses;
  for (classno = 0; classno < nclasses; classno++) {
    gcl = &gc->classes[classno];
    nobs = gcl->nobs;
    if (nobs) {
      for (row = 1; row <= gcl->m_u->rows; row++) gcl->m_u->rptr[row][1] /= (float)nobs;
      for (row = 1; row <= gcl->m_covariance->rows; row++) {
        mean_a = gcl->m_u->rptr[row][1];
        for (col = 1; col <= row; col++) {
          mean_b = gcl->m_u->rptr[col][1];
          covariance = gcl->m_covariance->rptr[row][col] / nobs - mean_a * mean_b;
          gcl->m_covariance->rptr[row][col] = covariance;
          gcl->m_covariance->rptr[col][row] = covariance;
        }
        if (Gdiag & DIAG_SHOW)
          fprintf(stderr,
                  "mean %24.24s, feature %d = %-2.3f, var = %-2.3f\n",
                  class_names[classno],
                  row,
                  mean_a,
                  gcl->m_covariance->rptr[row][row]);
      }
    }
    GCinit(gc, classno);
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
const char *MRICclassName(MRIC *mric, int round, int classno)
{
  const char *class_name = "unknown";

  switch (mric->type[round]) {
    default:
    case CLASSIFIER_GAUSSIAN:
      class_name = class_names[classno];
      break;
  }

  return (class_name);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRICdump(FILE *fp, MRIC *mric)
{
  GCLASSIFY *gc;
  GCLASS *gcl;
  int classno, nclasses, fno, round;

  fprintf(
      stderr, "mric with %d round, prior file %s\n", mric->nrounds, mric->mri_priors ? mric->prior_fname : "(none)");
  for (round = 0; round < mric->nrounds; round++) {
    gc = mric->classifier[round].gc;
    nclasses = gc->nclasses;
    for (classno = 0; classno < nclasses; classno++) {
      gcl = &gc->classes[classno];
      fprintf(stderr, "  %s:\n", class_names[classno]);
      for (fno = 0; fno < mric->ninputs[round]; fno++) {
        fprintf(stderr,
                "    feature %10.10s, mean = %+2.3f, std = %+2.3f\n",
                MRICfeatureName(mric, round, fno),
                gcl->m_u->rptr[fno + 1][1],
                sqrt(gcl->m_covariance->rptr[fno + 1][fno + 1]));
      }
    }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
const char *MRICfeatureName(MRIC *mric, int round, int feature_number)
{
  unsigned f;
  int fno;

  /* first ninputs-1 correspond to inputs #s, rest to frames in priors */

  /* find bit which corresponds to this # */
  for (f = 0x001, fno = 0; f != MAX_FEATURE; f <<= 1)
    if ((f & mric->features[round]) && (fno++ == feature_number)) break;

  if (f & FEATURE_INTENSITY) return ("INTENSITY");
  if (f & FEATURE_ZSCORE3) return ("ZSCORE3");
  if (f & FEATURE_ZSCORE5) return ("ZSCORE5");
  if (f & FEATURE_CPOLV_ZSCORE5) return ("CPOLV ZSCORE5");
  if (f & FEATURE_CPOLV_ZSCORE7) return ("CPOLV ZSCORE7");
  if (f & FEATURE_CPOLV_CURV5) return ("CPOLV NORMAL CURVATURE 5");
  if (f & FEATURE_CPOLV_CURV7) return ("CPOLV NORMAL CURVATURE 7");
  if (f & FEATURE_DIRECTION) return ("DIRECTION");
  if (f & FEATURE_MEAN3) return ("MEAN3");
  if (f & FEATURE_MEAN5) return ("MEAN5");
  if (f & FEATURE_CPOLV_MEAN3) return ("CPOLV MEAN3");
  if (f & FEATURE_CPOLV_MEAN5) return ("CPOLV MEAN5");
  if (f & FEATURE_CPOLV_MEDIAN3) return ("CPOLV MEDIAN3");
  if (f & FEATURE_CPOLV_MEDIAN5) return ("CPOLV MEDIAN 5");
  if (f & FEATURE_MIN3) return ("MIN3");
  if (f & FEATURE_MIN5) return ("MIN5");
  if (f & FEATURE_MIN7) return ("MIN7");
  if (f & FEATURE_PRIORS) return ("PRIORS");
  if (f & FEATURE_POSITION) return ("POSITION");

  return ("unknown");
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
const char *MRICfeatureNumberToName(int feature_number)
{
  unsigned f;
  int fno;

  /* find bit which corresponds to this # */
  for (f = 0x001, fno = 0; f != MAX_FEATURE; f <<= 1)
    if (fno++ == feature_number) break;

  if (f & FEATURE_INTENSITY) return ("INTENSITY");
  if (f & FEATURE_ZSCORE3) return ("ZSCORE3");
  if (f & FEATURE_ZSCORE5) return ("ZSCORE5");
  if (f & FEATURE_CPOLV_ZSCORE5) return ("CPOLV ZSCORE5");
  if (f & FEATURE_CPOLV_ZSCORE7) return ("CPOLV ZSCORE7");
  if (f & FEATURE_CPOLV_CURV5) return ("CPOLV NORMAL CURVATURE 5");
  if (f & FEATURE_CPOLV_CURV7) return ("CPOLV NORMAL CURVATURE 7");
  if (f & FEATURE_DIRECTION) return ("DIRECTION");
  if (f & FEATURE_MEAN3) return ("MEAN3");
  if (f & FEATURE_MEAN5) return ("MEAN5");
  if (f & FEATURE_CPOLV_MEAN3) return ("CPOLV MEAN3");
  if (f & FEATURE_CPOLV_MEAN5) return ("CPOLV MEAN5");
  if (f & FEATURE_CPOLV_MEDIAN3) return ("CPOLV MEDIAN3");
  if (f & FEATURE_CPOLV_MEDIAN5) return ("CPOLV MEDIAN 5");
  if (f & FEATURE_MIN3) return ("MIN3");
  if (f & FEATURE_MIN5) return ("MIN5");
  if (f & FEATURE_MIN7) return ("MIN7");
  if (f & FEATURE_PRIORS) return ("PRIORS");
  if (f & FEATURE_POSITION) return ("POSITION");
#if 0
  if (f & FEATURE_X_POSITION)
    return("X POSITION") ;
  if (f & FEATURE_Y_POSITION)
    return("Y POSITION") ;
  if (f & FEATURE_Z_POSITION)
    return("Z POSITION") ;
#endif
  if (f & FEATURE_CPOLV_ORDER) return ("CPOLV ORDER");

  return ("unknown");
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRICfeatureNumberCode(int feature_number)
{
  unsigned f;
  int fno;

  /* find bit which corresponds to this # */
  for (f = 0x001, fno = 0; f != MAX_FEATURE; f <<= 1)
    if (fno++ == feature_number) break;

  if (f & FEATURE_POSITION) return (FEATURE_POSITION);

  return (f);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRICfeatureCode(MRIC *mric, int round, int feature_number)
{
  unsigned f;
  int fno;

  /* first ninputs-1 correspond to inputs #s, rest to frames in priors */

  /* find bit which corresponds to this # */
  for (f = 0x001, fno = 0; f != MAX_FEATURE; f <<= 1)
    if ((f & mric->features[round]) && (fno++ == feature_number)) break;

  return (f);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
int MRICfeatureNumber(MRIC *mric, int round, int feature_code)
{
  unsigned f;
  int fno;

  /* first ninputs-1 correspond to inputs #s, rest to frames in priors */

  /* find bit which corresponds to this # */
  for (f = 0x001, fno = 0; f != MAX_FEATURE; f <<= 1) {
    if (f & feature_code) return (fno);
    if (f & mric->features[round]) fno++;
  }

  return (-1);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          compute the statistics (means and covariances) for a
          Gaussian classifier.
------------------------------------------------------*/
static int mricComputeGCStatistics(MRIC *mric, FILE *fp, int nfiles, int round)
{
  char source_fname[100], target_fname[100], line[300], *cp;
  int fno;
  MRI *mri_src, *mri_target;
  MRI_REGION bounding_box;

  /* now calculate statistics */
  fprintf(stderr, "computing classifier statistics...\n");
  fno = 0;
  while ((cp = fgetl(line, 299, fp)) != NULL) {
    sscanf(cp, "%s %s", source_fname, target_fname);
    fprintf(stderr, "file[%d]: %s --> %s\n", fno, source_fname, target_fname);
    mri_src = MRIread(source_fname);
    if (!mri_src) {
      fprintf(stderr, "could not read MR image %s\n", source_fname);
      continue;
    }

    mri_target = MRIread(target_fname);
    if (!mri_target) {
      fprintf(stderr, "could not read MR image %s\n", target_fname);
      MRIfree(&mri_src);
      continue;
    }

    MRIboundingBox(mri_src, DEFINITELY_BACKGROUND, &bounding_box);
    REGIONexpand(&bounding_box, &bounding_box, POLV_WSIZE / 2);
    MRIclipRegion(mri_src, &bounding_box, &bounding_box);
    REGIONcopy(&bounding_box, &mri_src->roi);
    MRICupdateStatistics(mric, round, mri_src, mri_target, &bounding_box);

    MRIfree(&mri_src);
    MRIfree(&mri_target);
    fno++;
  }

  MRICcomputeStatistics(mric, round); /* divide by # of observations */
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Train a Radial Basis Function classifier.
------------------------------------------------------*/
static int mricTrainRBF(MRIC *mric, FILE *fp, int nfiles, int round)
{
  char *cp, source_fname[100], line[100];
  GET_INPUT_PARMS parms;
  int fno;
  MRI *mri;

  parms.fp = fp;
  parms.mric = mric;
  parms.round = round;

  fno = 0;
  rewind(fp);
  while ((cp = fgetl(line, 299, fp)) != NULL) {
    sscanf(cp, "%s %*s", source_fname);
    mri = MRIreadInfo(source_fname);
    parms.npixels[fno] = mri->width * mri->height * mri->depth;
    fno++;
    MRIfree(&mri);
  }
  rewind(fp);

  if (RBFtrain(mric->classifier[round].rbf, mricGetClassifierInput, &parms, 0.0f) != NO_ERROR) return (Gerror);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Train a Radial Basis Function classifier.
------------------------------------------------------*/
static int mricRetrainRBF(MRIC *mric, FILE *fp, int nfiles, int round)
{
  char *cp, source_fname[100], line[100];
  GET_INPUT_PARMS parms;
  int fno;
  MRI *mri;

  parms.fp = fp;
  parms.mric = mric;
  parms.round = round;

  fno = 0;
  rewind(fp);
  while ((cp = fgetl(line, 299, fp)) != NULL) {
    sscanf(cp, "%s %*s", source_fname);
    mri = MRIreadInfo(source_fname);
    parms.npixels[fno] = mri->width * mri->height * mri->depth;
    fno++;
    MRIfree(&mri);
  }
  rewind(fp);

  if (RBFretrain(mric->classifier[round].rbf, mricGetClassifierInput, &parms, 0.0f) != NO_ERROR) return (Gerror);

  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
           Fill the elements of a GET_INPUT_PARMS structure.
------------------------------------------------------*/
static int mricFillInputParms(MRIC *mric, FILE *fp, int round, GET_INPUT_PARMS *parms)
{
  char *cp, source_fname[100], line[100];
  int fno;
  MRI *mri;

  parms->fp = fp;
  parms->mric = mric;
  parms->round = round;

  fno = 0;
  rewind(fp);
  while ((cp = fgetl(line, 299, fp)) != NULL) {
    sscanf(cp, "%s %*s", source_fname);
    mri = MRIreadInfo(source_fname);
    parms->npixels[fno] = mri->width * mri->height * mri->depth;
    fno++;
    MRIfree(&mri);
  }
  rewind(fp);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          generate a set of inputs for a specific point in space for
          a specific image.
------------------------------------------------------*/
static int mricGetClassifierInput(VECTOR *v_inputs, int no, void *parm, int same_class, int *pclass)
{
  char *cp, source_fname[100], line[100], wm_fname[100];
  int ino, x, y, z, obs_no, width, height, depth, fno, round;
  GET_INPUT_PARMS *parms;
  MRI *mri_wm;
  MRIC *mric;
  FILE *fp;
  static MRI *mri_src = NULL, *mri_target = NULL;
  static int current_ino = -1;

  parms = (GET_INPUT_PARMS *)parm;
  round = parms->round;
  mric = parms->mric;
  fp = parms->fp;

  /* find appropriate image # */
  ino = 0;
  obs_no = no;
  while (obs_no > parms->npixels[ino]) {
    obs_no -= parms->npixels[ino];
    ino++;
  }

  /* now load the appropriate MR image */
  if (current_ino != ino) /* not the same as last call */
  {
    if (mri_src) /* free old images */
      MRIfree(&mri_src);

    current_ino = ino;
    rewind(fp);
    fno = 0;
    while ((cp = fgetl(line, 299, fp)) != NULL) {
      sscanf(cp, "%s %s", source_fname, wm_fname);
      if (fno++ == ino) break;
    }
    if (!cp) return (ERROR_NO_FILE);
    mri_src = MRIread(source_fname);
    if (!mri_src) {
      fclose(fp);
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "mricGetClassifierInput: could not "
                   "open src file %s",
                   source_fname));
    }
    mri_wm = MRIread(wm_fname);
    if (!mri_wm) {
      MRIfree(&mri_src);
      fclose(fp);
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "mricGetClassifierInput: could not "
                   "open wmfilter file %s",
                   wm_fname));
    }
#if 0
    MRIvox(mri_src, 60, 63, 63) = 200 ;
    MRIvox(mri_src, 61, 63, 63) = 201 ;
    MRIvox(mri_src, 62, 63, 63) = 199 ;
    MRIvox(mri_src, 63, 63, 63) = 198 ;
    MRIvox(mri_wm, 60, 63, 63) = 0 ;
    MRIvox(mri_wm, 61, 63, 63) = 0 ;
    MRIvox(mri_wm, 62, 63, 63) = 0 ;
    MRIvox(mri_wm, 63, 63, 63) = 0 ;
#endif
    mri_target = MRICbuildTargetImage(mri_src, mri_target, mri_wm, 0, 0);
    MRIfree(&mri_wm);
  }

  width = mri_src->width;
  height = mri_src->height;
  depth = mri_src->depth;
  z = obs_no / (width * height);
  if (z >= depth) return (ERROR_BADPARM);
  y = (obs_no - z * width * height) / width;
  x = obs_no % width;
  MRICcomputeInputs(mric, mri_src, x, y, z, v_inputs, mric->features[round]);
  *pclass = MRIvox(mri_target, x, y, z);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Examine the statistical distribution of features in
          the training set (means and variances). Useful for
          clustering analysis.
------------------------------------------------------*/
int MRICexamineTrainingSet(MRIC *mric, char *file_name, int round)
{
  char *cp, source_fname[100], line[100];
  GET_INPUT_PARMS parms;
  int fno;
  MRI *mri;
  FILE *fp;

  fp = fopen(file_name, "r");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRICexamineTrainingSet(%s): could not open file", file_name));
  parms.fp = fp;
  parms.mric = mric;
  parms.round = round;

  fno = 0;
  rewind(fp);
  while ((cp = fgetl(line, 299, fp)) != NULL) {
    sscanf(cp, "%s %*s", source_fname);
    mri = MRIreadInfo(source_fname);
    parms.npixels[fno] = mri->width * mri->height * mri->depth;
    fno++;
  }
  rewind(fp);

  parms.fp = fp;
  parms.mric = mric;
  parms.round = round;
  RBFexamineTrainingSet(mric->classifier[round].rbf, mricGetClassifierInput, &parms);

  return (NO_ERROR);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Build a pairwise scatter plot of the (normalized)
          features in the given mric.
------------------------------------------------------*/
#define SCATTER_ROUND 0

int MRICbuildScatterPlot(MRIC *mric, int classnum, MATRIX *m_scatter, char *training_file_name)
{
  int obs_no = 0, i, x, y, nbins;
  // int half_bins, bin_offset;
  VECTOR *v_obs;
#if 0
  static int       first = 0 ;
  float            *means, *stds, mean, std ;
#endif
  float v, z[2] = {0, 0}, *mins, *maxs, mn, mx;
  RBF *rbf;
  GET_INPUT_PARMS parms;
  FILE *fp;

  nbins = m_scatter->rows;
// half_bins = (nbins - 1) / 2;
// bin_offset = half_bins;
#if 0
  if (!first)
  {
    first = 1 ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "examining training set...") ;
    MRICexamineTrainingSet(mric, training_file_name, SCATTER_ROUND) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "done.\n") ;
  }
#endif

  fp = fopen(training_file_name, "r");
  if (!fp)
    ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MRICbuildScatterPlot(%s): could not open file", training_file_name));
  mricFillInputParms(mric, fp, SCATTER_ROUND, &parms);

  rbf = mric->classifier[SCATTER_ROUND].rbf;

  v_obs = VectorAlloc(rbf->ninputs, MATRIX_REAL);

  /* not really a good idea to be messing with the CS internal parameters,
     but it's much easier than the alternative, because all the other CS stuff
     is class-specific, but the means and stds for normalizing the inputs
     should be across all classes.
     */

  /* now fill in scatter plot */
  mins = rbf->min_inputs;
  maxs = rbf->max_inputs;
  while (mricGetClassObservation(&parms, v_obs, obs_no++, classnum) == NO_ERROR) {
    for (i = 0; i < rbf->ninputs; i++) {
      mn = mins[i];
      mx = maxs[i];
      v = VECTOR_ELT(v_obs, i + 1);
      z[i] = (v - mn) / (mx - mn); /* scale it to 0 --> 1 */
    }

    /* map 0 to half_bins, -MAX_SIGMA*sigma to 0, MAX_SIGMA*sigma to nbins */
    x = z[0] * (nbins - 1) + 1;
    y = z[1] * (nbins - 1) + 1;
    if (!FZERO(VECTOR_ELT(v_obs, 1)) || !FZERO(VECTOR_ELT(v_obs, 2))) m_scatter->rptr[x][y] += 1.0f;
  }

  VectorFree(&v_obs);
  return (NO_ERROR);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Find the next observation for a given class. Stolen
          from rbf.c, but both really need it so....
------------------------------------------------------*/
static int mricGetClassObservation(GET_INPUT_PARMS *parms, VECTOR *v_obs, int desired_class_obs_no, int classnum)
{
  int ret, classno, obs_no, class_obs_no;
  static int last_class = -1;
  static int last_class_obs = -1;
  static int last_obs = -1;

  desired_class_obs_no++; /* it is a count - not an index */
  classno = classnum;
  if (classno == last_class && desired_class_obs_no > last_class_obs) {
    /* start at one past previous observation, not at start */
    class_obs_no = last_class_obs;
    obs_no = last_obs + 1;
  }
  else
    class_obs_no = obs_no = 0;

  do {
    ret = mricGetClassifierInput(v_obs, obs_no++, parms, 1, &classno);
    if ((ret == NO_ERROR) && (classno == classnum)) class_obs_no++;
  } while ((ret == NO_ERROR) && (class_obs_no < desired_class_obs_no));

  if (ret == NO_ERROR) {
    last_class = classno;
    last_obs = obs_no - 1;
    last_class_obs = desired_class_obs_no;
  }
  else
    last_class = last_obs = last_class_obs = -1;

  return (ret);
}
