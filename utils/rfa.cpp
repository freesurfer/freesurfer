/**
 * @brief Creates the Random Forest Array (RFA) atlas from training set
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

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cma.h"
#include "diag.h"
#include "error.h"
#include "gca.h"
#include "gcamorph.h"
#include "macros.h"
#include "rfa.h"
#include "talairachex.h"
#include "transform.h"
#include "utils.h"

#define MAX_LABELS_PER_NODE 10

RFA *RFAalloc(RFA_PARMS *parms, int alloc_trees)
{
  int x, y, z;
  RF_NODE *node;
  RFA *rfa;

  rfa = (RFA *)calloc(1, sizeof(RFA));
  if (rfa == NULL) ErrorExit(ERROR_NOMEMORY, "RFAalloc: could not allocate rfa");

  rfa->spacing = parms->spacing;
  rfa->vg.width = parms->width / parms->spacing;
  rfa->vg.height = parms->height / parms->spacing;
  rfa->vg.depth = parms->depth / parms->spacing;
  rfa->width = rfa->vg.width;
  rfa->height = rfa->vg.height;
  rfa->depth = rfa->vg.depth;
  rfa->vg.xsize = parms->spacing;
  rfa->vg.ysize = parms->spacing;
  rfa->vg.zsize = parms->spacing;
  rfa->vg.valid = 0;
  rfa->ninputs = parms->nvols;
  rfa->wsize = parms->wsize;
  parms->training_index = 0;  // this will get incremented in RFAaddInput

  // +3 is for the 3 spatial coordinates
  rfa->nfeatures = (parms->wsize * parms->wsize * parms->wsize) * parms->nvols + 3;
  rfa->nodes = (RF_NODE ***)calloc(rfa->vg.width, sizeof(RF_NODE **));
  if (rfa->nodes == NULL) ErrorExit(ERROR_NOMEMORY, "RFAalloc: could not allocate node X");

  for (x = 0; x < rfa->vg.width; x++) {
    rfa->nodes[x] = (RF_NODE **)calloc(rfa->vg.height, sizeof(RF_NODE *));
    if (rfa->nodes[x] == NULL) ErrorExit(ERROR_NOMEMORY, "RFAalloc: could not allocate node Y");
    for (y = 0; y < rfa->vg.height; y++) {
      rfa->nodes[x][y] = (RF_NODE *)calloc(rfa->vg.height, sizeof(RF_NODE));
      if (rfa->nodes[x][y] == NULL) ErrorExit(ERROR_NOMEMORY, "RFAalloc: could not allocate node Z");
      for (z = 0; z < rfa->vg.depth; z++) {
        node = &rfa->nodes[x][y][z];
        node->label_priors = (float *)calloc(MAX_LABELS_PER_NODE, sizeof(node->label_priors[0]));
        if (node->label_priors == NULL) ErrorExit(ERROR_NOMEMORY, "RFAalloc: could not allocate node label_priors");
        node->labels = (unsigned short *)calloc(MAX_LABELS_PER_NODE, sizeof(node->labels[0]));
        if (node->labels == NULL) ErrorExit(ERROR_NOMEMORY, "RFAalloc: could not allocate node labels");
        node->max_labels = MAX_LABELS_PER_NODE;
        if (alloc_trees)
          node->rf = RFalloc(parms->ntrees, rfa->nfeatures, MAX_LABELS_PER_NODE, parms->max_depth, NULL, 10);
      }
    }
  }

  return (rfa);
}

int RFAtrain(RFA *rfa, MRI *mri_inputs, MRI *mri_labels, TRANSFORM *transform) { return (NO_ERROR); }

int RFAcompleteTraining(RFA *rfa, RFA_PARMS *parms)
{
  RFAupdateTraining(rfa, parms);
  return (NO_ERROR);
}

int RFAaddInput(RFA *rfa, MRI *mri_seg, MRI *mri_inputs, TRANSFORM *transform, RFA_PARMS *parms)
{
  int index;

  index = (parms->training_index % parms->training_size);
  parms->training_index++;
  parms->mri_segs[index] = mri_seg;
  parms->mri_inputs[index] = mri_inputs;
  parms->transforms[index] = transform;

  if ((parms->training_index + 1) % parms->training_size == 0) {
    int i;

    RFAupdateTraining(rfa, parms);

    for (i = 0; i < parms->training_size; i++) {
      MRIfree(&parms->mri_inputs[i]);
      MRIfree(&parms->mri_segs[i]);
      TransformFree(&parms->transforms[i]);
    }
  }

  return (NO_ERROR);
}

int RFAvoxelToNode(const RFA *rfa, double xt, double yt, double zt, double *px, double *py, double *pz)
{
  *px = xt / rfa->spacing;
  *py = yt / rfa->spacing;
  *pz = zt / rfa->spacing;
  return (NO_ERROR);
}

int RFAnodeToVoxel(const RFA *rfa, double xt, double yt, double zt, double *px, double *py, double *pz)
{
  *px = xt * rfa->spacing;
  *py = yt * rfa->spacing;
  *pz = zt * rfa->spacing;
  return (NO_ERROR);
}

int RFAsourceVoxelToNode(
    const RFA *rfa, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv, int *px, int *py, int *pz)
{
  float xt = 0, yt = 0, zt = 0;
  double xrt, yrt, zrt, xd, yd, zd;
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
      ErrorExit(ERROR_BADPARM, "RFAsourceVoxelToNode: needs vox-to-vox transform");
  }
  else  // morph 3d type can go directly from source to template
  {
    TransformSample(transform, xv, yv, zv, &xt, &yt, &zt);
  }
  RFAvoxelToNode(rfa, xt, yt, zt, &xd, &yd, &zd);
  *px = (int)(xd);
  *py = (int)(yd);
  *pz = (int)(zd);

  // if (*px < 0 || *py < 0 || *pz < 0 || *px >= rfa->width || *py >= rfa->height || *pz >= rfa->depth)
  //   retval = (ERROR_BADPARM);
  if (*px < 0) *px = 0;
  if (*py < 0) *py = 0;
  if (*pz < 0) *pz = 0;
  if (*px >= rfa->width) *px = rfa->width - 1;
  if (*py >= rfa->height) *py = rfa->height - 1;
  if (*pz >= rfa->depth) *pz = rfa->depth - 1;
  return (NO_ERROR);
}
int RFAsourceVoxelToAtlas(
    const RFA *rfa, MRI *mri, TRANSFORM *transform, int xv, int yv, int zv, double *px, double *py, double *pz)
{
  float xt, yt, zt;
  LTA *lta;

  if (transform->type != MORPH_3D_TYPE) {
    if (transform->type == LINEAR_VOX_TO_VOX) {
      lta = (LTA *)transform->xform;
      // transform point to talairach volume point
      TransformWithMatrix(lta->xforms[0].m_L, xv, yv, zv, px, py, pz);
      // TransformSample(transform, xv, yv, zv, &xt, &yt, &zt) ;
    }
    else
      ErrorExit(ERROR_BADPARM, "RFAsourceVoxelToNode: needs vox-to-vox transform");
  }
  else  // morph 3d type can go directly from source to template
  {
    TransformSample(transform, xv, yv, zv, &xt, &yt, &zt);
    *px = (double)xt;
    *py = (double)yt;
    *pz = (double)zt;
  }
#if 0
  if (*px < 0 || *py < 0 || *pz < 0 ||
      *px >= rfa->width ||
      *py >= rfa->height ||
      *pz >= rfa->depth)
    retval = (ERROR_BADPARM) ;
  if (*px < 0)
    *px= 0 ;
  if (*py < 0)
    *py= 0 ;
  if (*pz < 0)
    *pz= 0 ;
  if (*px >= rfa->width)
    *px = rfa->width-1 ;
  if (*py >= rfa->height)
    *py = rfa->height-1 ;
  if (*pz >= rfa->depth)
    *pz = rfa->depth-1 ;
#endif
  return (NO_ERROR);
}
int RFAnodeToSourceVoxel(
    RFA *rfa, MRI *mri, TRANSFORM *transform, int xn, int yn, int zn, double *pxv, double *pyv, double *pzv)
{
  int width, height, depth;
  double xt, yt, zt;
  double xv, yv, zv;
  double xc, yc, zc;
  float xf, yf, zf;
  int errCode = NO_ERROR;
  LTA *lta;

  width = mri->width;
  height = mri->height;
  depth = mri->depth;
  // get template voxel position
  RFAnodeToVoxel(rfa, xn, yn, zn, &xt, &yt, &zt);
  if (transform->type != MORPH_3D_TYPE) {
    lta = (LTA *)transform->xform;
    // get the talairach to orig
    TransformWithMatrix(lta->inv_xforms[0].m_L, xt, yt, zt, &xc, &yc, &zc);
    // TransformSampleInverse(transform, xt, yt, zt, &xc, &yc, &zc);
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
    xv = xc;
    yv = yc;
    zv = zc;
  }
  else  // template to source
  {
    TransformSampleInverse(transform, xt, yt, zt, &xf, &yf, &zf);
    xv = (double)xf;
    yv = (double)yf;
    zv = (double)zf;
  }
  *pxv = xv;
  *pyv = yv;
  *pzv = zv;
  return errCode;
}

int RFAupdateTraining(RFA *rfa, RFA_PARMS *parms)
{
  int ntraining_sets, ntraining, start, x, y, z, xn, yn, zn, i, *training_classes, *wmsa_permutation;
  int xk, yk, zk, xi, yi, zi, tno, label, index, total_correct, total_training, n, nwmsa, max_count, wmsa_count,
      wmsa_label = WM_hypointensities;
  double **training_data, xd, yd, zd, **wmsa_data, xatlas, yatlas, zatlas;
  MRI *mri_in, *mri_seg;
  TRANSFORM *transform;
  RF_NODE *node;

  x = y = z = 0;  // silly compiler warning
  ntraining_sets = parms->training_index % parms->training_size;
  start = parms->training_index - ntraining_sets;
  if (start < 0 || ntraining_sets <= 0) return (NO_ERROR);  // no more to do

  ntraining = ntraining_sets * rfa->spacing * rfa->spacing * rfa->spacing;
  training_data = (double **)calloc(2 * ntraining, sizeof(training_data[0]));
  training_classes = (int *)calloc(2 * ntraining, sizeof(training_classes[0]));
  if (training_data == NULL || training_classes == NULL)
    ErrorExit(ERROR_NOFILE, "RFAupdateTraining: could not allocate %d-length training buffers", ntraining);

  total_correct = total_training = 0;
  if (rfa->vg.valid == 0) {
    if (parms->transforms[0]->type == MORPH_3D_TYPE) {
      GCA_MORPH *gcam = (GCA_MORPH *)(parms->transforms[0]->xform);
      rfa->vg = *(&gcam->atlas);
    }
    else {
      LTA *lta = (LTA *)(parms->transforms[0]->xform);
      rfa->vg = *(&lta->xforms[0].dst);
    }
  }
  for (i = 0; i < ntraining; i++)  // allow for augmenting with other wmsa examples
  {
    training_data[i] = (double *)calloc(rfa->nfeatures, sizeof(double));
    if (training_data[i] == NULL)
      ErrorExit(ERROR_NOMEMORY, "RFAupdateTraining: could not allocate %d-len feature vector #%d", rfa->nfeatures, i);
  }

  // count total # of wmsas in training set
  for (nwmsa = i = 0; i < ntraining_sets; i++) {
    mri_in = parms->mri_inputs[i];
    mri_seg = parms->mri_segs[i];
    transform = parms->transforms[i];
    nwmsa += MRIvoxelsInLabel(mri_seg, WM_hypointensities);
    nwmsa += MRIvoxelsInLabel(mri_seg, Left_WM_hypointensities);
    nwmsa += MRIvoxelsInLabel(mri_seg, Right_WM_hypointensities);
  }
  printf("%d wmsa voxels found in training data\n", nwmsa);
  wmsa_data = (double **)calloc(nwmsa, sizeof(wmsa_data[0]));
  if (wmsa_data == NULL) ErrorExit(ERROR_NOFILE, "RFAupdateTraining: could not allocate %d-length wmsa buffers", nwmsa);
  wmsa_permutation = (int *)calloc(nwmsa, sizeof(wmsa_permutation[0]));
  if (wmsa_permutation == NULL)
    ErrorExit(ERROR_NOFILE, "RFAupdateTraining: could not allocate %d-length wmsa buffers", nwmsa);
  for (nwmsa = i = 0; i < ntraining_sets; i++) {
    mri_in = parms->mri_inputs[i];
    mri_seg = parms->mri_segs[i];
    transform = parms->transforms[i];
    for (x = 0; x < mri_in->width; x++)
      for (y = 0; y < mri_in->height; y++)
        for (z = 0; z < mri_in->depth; z++) {
          label = MRIgetVoxVal(mri_seg, x, y, z, 0);
          if (IS_HYPO(label)) {
            wmsa_data[nwmsa] = (double *)calloc(rfa->nfeatures, sizeof(double));
            if (wmsa_data[nwmsa] == NULL)
              ErrorExit(ERROR_NOMEMORY, "RFAtrain: could not allocate %dth wmsa feature vector", nwmsa);
            extract_feature(mri_in, rfa->wsize, x, y, z, wmsa_data[nwmsa++], -1, -1, -1);
          }
        }
  }
  compute_permutation(nwmsa, wmsa_permutation);

  // for every RFA node
  for (xn = 0; xn < rfa->width; xn++) {
    if ((xn % 10) == 0) printf("training node series %d of %d\n", xn, rfa->width);
    for (yn = 0; yn < rfa->height; yn++)
      for (zn = 0; zn < rfa->depth; zn++) {
        if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) DiagBreak();
        node = &rfa->nodes[xn][yn][zn];
        for (tno = i = 0; i < ntraining_sets; i++) {
          mri_in = parms->mri_inputs[i];
          mri_seg = parms->mri_segs[i];
          transform = parms->transforms[i];
          RFAnodeToSourceVoxel(rfa, mri_in, transform, xn, yn, zn, &xd, &yd, &zd);
          x = (int)(xd);
          y = (int)(yd);
          z = (int)(zd);
          for (xk = 0; xk < rfa->spacing; xk++) {
            xi = mri_in->xi[x + xk];
            for (yk = 0; yk < rfa->spacing; yk++) {
              yi = mri_in->yi[y + yk];
              for (zk = 0; zk < rfa->spacing; zk++) {
                zi = mri_in->zi[z + zk];
                RFAsourceVoxelToAtlas(rfa, mri_in, transform, xi, yi, zi, &xatlas, &yatlas, &zatlas);
                extract_feature(mri_in, rfa->wsize, xi, yi, zi, training_data[tno], xatlas, yatlas, zatlas);
                label = MRIgetVoxVal(mri_seg, xi, yi, zi, 0);
                for (index = 0; index < node->nlabels; index++)
                  if (node->labels[index] == label) break;
                if (index >= node->nlabels) {
                  if (index >= node->max_labels) {
                    float label_priors[100 * MAX_LABELS_PER_NODE];
                    unsigned short labels[100 * MAX_LABELS_PER_NODE];
                    memmove(label_priors, node->label_priors, node->nlabels * sizeof(node->label_priors[0]));
                    memmove(labels, node->labels, node->nlabels * sizeof(node->labels[0]));
                    free(node->label_priors);
                    free(node->labels);
                    node->max_labels += 5;
                    node->labels = (unsigned short *)calloc(node->max_labels, sizeof(node->labels[0]));
                    if (node->labels == NULL)
                      ErrorExit(ERROR_NOMEMORY,
                                "could not reallocate %d node labels @(%d %d %d)",
                                node->max_labels,
                                xn,
                                yn,
                                zn);
                    node->label_priors = (float *)calloc(node->max_labels, sizeof(node->label_priors[0]));
                    if (node->label_priors == NULL)
                      ErrorExit(ERROR_NOMEMORY,
                                "could not reallocate %d node priors @(%d %d %d)",
                                node->max_labels,
                                xn,
                                yn,
                                zn);
                    memmove(node->label_priors, label_priors, node->nlabels * sizeof(node->label_priors[0]));
                    memmove(node->labels, labels, node->nlabels * sizeof(node->labels[0]));
                  }
                  node->labels[node->nlabels++] = label;
                }
                if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z && (Ggca_label < 0 || Ggca_label == label)) {
                  int j;
                  printf("tdata[%d]: ", tno);
                  for (j = 0; j < rfa->nfeatures; j++) printf("%.0f ", training_data[tno][j]);
                  printf(" - %s (%d)\n", cma_label_to_name(label), label);
                }
                training_classes[tno] = index;
                node->label_priors[index]++;
                tno++;
              }
            }
          }
        }
        for (max_count = wmsa_count = n = 0; n < node->nlabels; n++) {
          if (IS_HYPO(node->labels[n])) {
            wmsa_label = n;
            wmsa_count += node->label_priors[n];
          }
          if (node->label_priors[n] > max_count) max_count = node->label_priors[n];
          node->label_priors[n] /= ntraining;
        }
        if (node->nlabels > 1) {
          int start_wmsa_index, wmsa_index, ind;
          int correct, extra = 0;

          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) DiagBreak();
          if (wmsa_count > 0) DiagBreak();
          if (wmsa_count > 0 && max_count > 5 * wmsa_count)  // augment training set
          {
            extra = (max_count / 2 - wmsa_count);
            extra = 0;  // disabled!!
            start_wmsa_index = (int)randomNumber(0.0, (double)((nwmsa - extra) - .00001));
            for (wmsa_index = start_wmsa_index; wmsa_index < start_wmsa_index + extra; wmsa_index++) {
              ind = ntraining + wmsa_index - start_wmsa_index;
              if (ind < 0 || ind >= 2 * ntraining) DiagBreak();
              training_data[ind] = wmsa_data[wmsa_permutation[wmsa_index]];
              training_classes[ind] = wmsa_label;
            }
          }
          RFsetNumberOfClasses(node->rf, node->nlabels);
          node->rf->min_step_size = .5;
          node->rf->max_steps = 256;
          RFtrain(node->rf,
                  parms->feature_fraction,
                  parms->training_fraction,
                  training_classes,
                  training_data,
                  ntraining + extra);
          correct = RFcomputeOutOfBagCorrect(node->rf, training_classes, training_data, ntraining + extra);
          total_training += (ntraining + extra);
          total_correct += correct;
        }
      }
  }

  printf("out of bag accuracy: %d of %d = %2.2f%%\n",
         total_correct,
         total_training,
         100.0 * total_correct / total_training);
  for (i = 0; i < ntraining; i++) free(training_data[i]);
  for (i = 0; i < nwmsa; i++) free(wmsa_data[i]);

  free(training_data);
  free(training_classes);
  free(wmsa_permutation);
  free(wmsa_data);
  return (NO_ERROR);
}
int RFAwrite(RFA *rfa, char *fname)
{
  FILE *fp;
  int xn, yn, zn, n;
  RF_NODE *node;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "RFAwrite(%s): could not open file", fname));
  fprintf(fp,
          "%d %d %d %d %d %d %d %2.1f\n",
          rfa->total_training,
          rfa->wsize,
          rfa->nfeatures,
          rfa->wsize,
          rfa->width,
          rfa->height,
          rfa->depth,
          rfa->spacing);

  for (xn = 0; xn < rfa->width; xn++)
    for (yn = 0; yn < rfa->height; yn++)
      for (zn = 0; zn < rfa->depth; zn++) {
        node = &rfa->nodes[xn][yn][zn];
        if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) DiagBreak();
        fprintf(fp, "rn %d %d %d: %d\n", xn, yn, zn, node->nlabels);
        for (n = 0; n < node->nlabels; n++) fprintf(fp, "%d ", node->labels[n]);
        fprintf(fp, "\n");
        for (n = 0; n < node->nlabels; n++) fprintf(fp, "%2.3f ", node->label_priors[n]);
        fprintf(fp, "\n");
        if (node->nlabels > 1) RFwriteInto(node->rf, fp);
      }

  fclose(fp);
  return (NO_ERROR);
}
RFA *RFAread(char *fname)
{
  RFA *rfa;
  FILE *fp;
  int xn, yn, zn, n, nfeatures, total_training, first = 1;
  RF_NODE *node;
  RFA_PARMS parms;

  memset(&parms, 0, sizeof(parms));
  rfa = (RFA *)calloc(1, sizeof(RFA));
  if (rfa == NULL) ErrorExit(ERROR_NOMEMORY, "RFAalloc: could not allocate rfa");
  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_BADPARM, "RFAwrite(%s): could not open file", fname));
  if (fscanf(fp,
             "%d %d %d %d %d %d %d %f\n",
             &total_training,
             &parms.wsize,
             &nfeatures,
             &parms.wsize,
             &parms.width,
             &parms.height,
             &parms.depth,
             &parms.spacing) != 8) {
    ErrorPrintf(ERROR_BAD_FILE, "RFAwrite: could not read file");
  }

  parms.width *= parms.spacing;
  parms.height *= parms.spacing;
  parms.depth *= parms.spacing;
  rfa = RFAalloc(&parms, 0);
  for (xn = 0; xn < rfa->width; xn++)
    for (yn = 0; yn < rfa->height; yn++)
      for (zn = 0; zn < rfa->depth; zn++) {
        int nlabels, label;
        char line[MAX_LINE_LEN], *cp;

        node = &rfa->nodes[xn][yn][zn];
        if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) DiagBreak();
        cp = fgetl(line, MAX_LINE_LEN, fp);
        sscanf(cp, "rn %*d %*d %*d: %d\n", &nlabels);
        node->nlabels = node->max_labels = nlabels;
        free(node->labels);
        free(node->label_priors);
        node->labels = (unsigned short *)calloc(node->max_labels, sizeof(node->labels[0]));
        if (node->labels == NULL)
          ErrorExit(ERROR_NOMEMORY, "could not reallocate %d node labels @(%d %d %d)", node->max_labels, xn, yn, zn);
        node->label_priors = (float *)calloc(node->max_labels, sizeof(node->label_priors[0]));
        if (node->label_priors == NULL)
          ErrorExit(ERROR_NOMEMORY, "could not reallocate %d node priors @(%d %d %d)", node->max_labels, xn, yn, zn);
        for (n = 0; n < node->nlabels; n++) {
          if (fscanf(fp, "%d ", &label) != 1) {
            ErrorPrintf(ERROR_BAD_FILE, "RFAwrite: could not read file");
          }
          node->labels[n] = label;
        }
        if (fscanf(fp, "\n") != 0) {
          ErrorPrintf(ERROR_BAD_FILE, "RFAwrite: could not read file");
        }
        for (n = 0; n < node->nlabels; n++) {
          if (fscanf(fp, "%f ", &node->label_priors[n]) != 1) {
            ErrorPrintf(ERROR_BAD_FILE, "RFAwrite: could not read file");
          }
        }
        if (fscanf(fp, "\n") != 0) {
          ErrorPrintf(ERROR_BAD_FILE, "RFAwrite: could not read file");
        }
        node->max_labels = node->nlabels;
        if (node->nlabels > 1) {
          if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z) DiagBreak();
          node->rf = RFreadFrom(fp);
          if (first) {
            rfa->nfeatures = node->rf->nfeatures;
            first = 0;
          }
        }
      }

  rfa->ninputs = rfa->nfeatures / (rfa->wsize * rfa->wsize * rfa->wsize);
  fclose(fp);

  return (rfa);
}

#define MAX_FEATURE_LEN 10000
MRI *RFAlabel(MRI *mri_in, RFA *rfa, MRI *mri_labeled, TRANSFORM *transform)
{
  int x, y, z, xn, yn, zn;
  int label;
  RF_NODE *node;
  double feature[MAX_FEATURE_LEN], xatlas, yatlas, zatlas;

  if (mri_labeled == NULL) {
    mri_labeled = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_SHORT);
    MRIcopyHeader(mri_in, mri_labeled);
  }

  // for every input voxel
  for (x = 0; x < mri_in->width; x++) {
    for (y = 0; y < mri_in->height; y++)
      for (z = 0; z < mri_in->depth; z++) {
        RFAsourceVoxelToNode(rfa, mri_in, transform, x, y, z, &xn, &yn, &zn);
        node = &rfa->nodes[xn][yn][zn];
        RFAsourceVoxelToAtlas(rfa, mri_in, transform, x, y, z, &xatlas, &yatlas, &zatlas);
        extract_feature(mri_in, rfa->wsize, x, y, z, feature, xatlas, yatlas, zatlas);
        if (x == Gx && y == Gy && z == Gz) {
          int j;
          printf("voxel (%d, %d, %d) -> node (%d, %d, %d), labels = ", x, y, zn, xn, yn, zn);
          for (j = 0; j < node->nlabels; j++) printf("%d (%.2f) ", node->labels[j], node->label_priors[j]);
          printf("\nfeature = ");
          for (j = 0; j < rfa->nfeatures; j++) printf("%.0f ", feature[j]);
          printf("\n");
          Gdiag |= DIAG_VERBOSE;
          DiagBreak();
        }
        if (node->nlabels == 0)
          label = 0;
        else if (node->nlabels == 1)
          label = node->labels[0];
        else {
          double pval;
          label = RFclassify(node->rf, feature, &pval, -1);
          label = node->labels[label];  // convert to CMA label
        }
        if (x == Gx && y == Gy && z == Gz) Gdiag &= ~DIAG_VERBOSE;
        MRIsetVoxVal(mri_labeled, x, y, z, 0, label);
      }
  }

  return (mri_labeled);
}

int extract_feature(MRI *mri_in, int wsize, int x, int y, int z, double *feature, int xatlas, int yatlas, int zatlas)
{
  int xi, yi, zi, xk, yk, zk, whalf, n, f;

  whalf = (wsize - 1) / 2;
  for (f = 0, xk = -whalf; xk <= whalf; xk++) {
    xi = mri_in->xi[x + xk];
    for (yk = -whalf; yk <= whalf; yk++) {
      yi = mri_in->yi[y + yk];
      for (zk = -whalf; zk <= whalf; zk++) {
        zi = mri_in->zi[z + zk];
        for (n = 0; n < mri_in->nframes; n++) feature[f++] = MRIgetVoxVal(mri_in, xi, yi, zi, n);
      }
    }
  }
  feature[f++] = xatlas;
  feature[f++] = yatlas;
  feature[f++] = zatlas;
  return (NO_ERROR);
}

int extract_long_features(
    MRI *mri_in, MRI *mri_seg, TRANSFORM *transform, GCA *gca, int wsize, int x, int y, int z, double *feature)
{
  int xi, yi, zi, xk, yk, zk, whalf, n, f;

  whalf = (wsize - 1) / 2;
  for (f = 0, xk = -whalf; xk <= whalf; xk++) {
    xi = mri_in->xi[x + xk];
    for (yk = -whalf; yk <= whalf; yk++) {
      yi = mri_in->yi[y + yk];
      for (zk = -whalf; zk <= whalf; zk++) {
        zi = mri_in->zi[z + zk];
        for (n = 0; n < mri_in->nframes; n++) feature[f++] = MRIgetVoxVal(mri_in, xi, yi, zi, n);
      }
    }
  }
  feature[f++] = MRIcountCSFInNbhd(mri_seg, 5, x, y, z);
  feature[f++] = 100 * gm_prior(gca, mri_in, transform, x, y, z);
  feature[f++] = 100 * wm_prior(gca, mri_in, transform, x, y, z);
  feature[f++] = 100 * csf_prior(gca, mri_in, transform, x, y, z);
  feature[f++] = 100 * MRIcountValInNbhd(mri_seg, 3, x, y, z, WM_hypointensities);
  return (NO_ERROR);
}

static int csf_labels[] = {CSF,
                           Unknown,
                           Left_Lateral_Ventricle,
                           Right_Lateral_Ventricle,
                           Left_Inf_Lat_Vent,
                           Right_Inf_Lat_Vent,
                           Third_Ventricle,
                           Fourth_Ventricle};
#define NCSF_LABELS (sizeof(csf_labels) / sizeof(csf_labels[0]))

int MRIcountCSFInNbhd(MRI *mri_seg, int wsize, int x, int y, int z)
{
  unsigned int total, n;

  for (n = total = 0; n < NCSF_LABELS; n++) total += MRIcountValInNbhd(mri_seg, wsize, x, y, z, csf_labels[n]);

  return (total);
}
