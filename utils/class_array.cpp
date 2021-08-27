/**
 * @brief utilities for using SVMs to reclassify aseg voxels
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
#include "class_array.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "diag.h"
#include "error.h"
#include "fio.h"
#include "macros.h"
#include "mri.h"
#include "mrinorm.h"
#include "voxlist.h"

extern const char *Progname;

static MATRIX *compute_ras_basis_vectors(
    MRI *mri_aseg_orig, MRI *mri_aseg_edit, int target_label, int width, int height, int depth, int pad);
static int clUpdateStatistics(CA *ca, CLASSIFIER *cl, float *inputs, int output);
static int clCompleteTraining(CA *ca, CLASSIFIER *cl);
static MRI *caComputeSurfaceNormals(MRI *mri_aseg, MRI *mri_normals, int label);

float **CAbuildTrainingData(VOXEL_LIST *vl_total,
                            int target_label,
                            float *ca_classes,
                            MRI **mri_smooth,
                            MRI **mri_grad,
                            MRI **mri_laplacian,
                            MRI *mri_dtrans,
                            MRI **mri_2nd_deriv_s,
                            int wsize,
                            int nscales,
                            int which_inputs)
{
  float **ca_inputs;
  int i, ninputs;
  MRI *mri_dtrans_grad;

  ca_inputs = (float **)calloc(vl_total->nvox, sizeof(float *));
  if (ca_inputs == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d SVM input vector", Progname, vl_total->nvox);

  ninputs = NINPUTS(which_inputs, wsize, nscales);
  if (which_inputs & CA_INPUT_DTRANS_GRAD) mri_dtrans_grad = MRIsobel(mri_dtrans, NULL, NULL);

  for (i = 0; i < vl_total->nvox; i++) {
    if (ca_classes) ca_classes[i] = ((int)vl_total->vdst[i] == target_label) ? 1 : -1;

    ca_inputs[i] = CAbuildInputsAtVoxel(vl_total,
                                        i,
                                        mri_smooth,
                                        mri_grad,
                                        mri_laplacian,
                                        mri_dtrans,
                                        mri_dtrans_grad,
                                        mri_2nd_deriv_s,
                                        wsize,
                                        nscales,
                                        NULL,
                                        which_inputs);
    if (ca_inputs[i] == NULL)
      ErrorExit(
          ERROR_NOMEMORY, "%s: could not allocate %dth CA input vector size %d", Progname, vl_total->nvox, ninputs);
  }
  if (which_inputs & CA_INPUT_DTRANS_GRAD) MRIfree(&mri_dtrans_grad);
  return (ca_inputs);
}
float *CAbuildInputsAtVoxel(VOXEL_LIST *vl,
                            int i,
                            MRI **mri_smooth,
                            MRI **mri_grad,
                            MRI **mri_laplacian,
                            MRI *mri_dtrans,
                            MRI *mri_dtrans_grad,
                            MRI **mri_2nd_deriv_s,
                            int wsize,
                            int nscales,
                            float *svm_inputs,
                            int which_inputs)
{
  int s, xk, yk, zk, xi, yi, zi, x, y, z, ninputs, whalf, input;

  whalf = (wsize - 1) / 2;
  ninputs = NINPUTS(which_inputs, wsize, nscales);
  if (svm_inputs == NULL) {
    svm_inputs = (float *)calloc(ninputs, sizeof(float));
    if (svm_inputs == NULL) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d SVM input vector", Progname, ninputs);
  }

  x = vl->xi[i];
  y = vl->yi[i];
  z = vl->zi[i];

  input = 0;
  if (which_inputs & CA_INPUT_DTRANS) svm_inputs[input++] = MRIgetVoxVal(mri_dtrans, x, y, z, 0);
  if (!std::isfinite(svm_inputs[input - 1])) DiagBreak();
  if (which_inputs & CA_INPUT_DTRANS_GRAD) {
    svm_inputs[input++] = MRIgetVoxVal(mri_dtrans_grad, x, y, z, 0);
    if (!std::isfinite(svm_inputs[input - 1])) DiagBreak();
    svm_inputs[input++] = MRIgetVoxVal(mri_dtrans_grad, x, y, z, 1);
    if (!std::isfinite(svm_inputs[input - 1])) DiagBreak();
    svm_inputs[input++] = MRIgetVoxVal(mri_dtrans_grad, x, y, z, 2);
    if (!std::isfinite(svm_inputs[input - 1])) DiagBreak();
  }
  for (xk = -whalf; xk <= whalf; xk++)
    for (yk = -whalf; yk <= whalf; yk++)
      for (zk = -whalf; zk <= whalf; zk++) {
        for (s = 0; s < nscales; s++) {
          xi = mri_smooth[s]->xi[x + xk];
          yi = mri_smooth[s]->yi[y + yk];
          zi = mri_smooth[s]->zi[z + zk];
          if (which_inputs & CA_INPUT_INTENSITY) {
            svm_inputs[input++] = MRIgetVoxVal(mri_smooth[s], xi, yi, zi, 0);
            if (!std::isfinite(svm_inputs[input - 1])) DiagBreak();
          }
          if (which_inputs & CA_INPUT_GRADIENT) {
            svm_inputs[input++] = MRIgetVoxVal(mri_grad[s], xi, yi, zi, 0);
            svm_inputs[input++] = MRIgetVoxVal(mri_grad[s], xi, yi, zi, 1);
            svm_inputs[input++] = MRIgetVoxVal(mri_grad[s], xi, yi, zi, 2);
          }

          if (which_inputs & CA_INPUT_LAPLACIAN) {
            svm_inputs[input++] = MRIgetVoxVal(mri_laplacian[s], xi, yi, zi, 0);
            if (!std::isfinite(svm_inputs[input - 1])) DiagBreak();
          }
          if (which_inputs & CA_INPUT_D2I_S) {
            svm_inputs[input++] = MRIgetVoxVal(mri_2nd_deriv_s[s], xi, yi, zi, 0);
            if (!std::isfinite(svm_inputs[input - 1])) DiagBreak();
          }
        }
      }

  for (input = 0; input < ninputs; input++)
    if (!std::isfinite(svm_inputs[input])) DiagBreak();
  return (svm_inputs);
}
#define MAX_SCALES 50
float **CAbuildInputs(VOXEL_LIST *vl_total,
                      MRI *mri_intensity,
                      MRI *mri_labels,
                      int target_label,
                      int which_inputs,
                      int wsize,
                      int nscales,
                      float *sigmas)
{
  float **svm_inputs;
  MRI *mri_grad[MAX_SCALES], *mri_kernel, *mri_smooth[MAX_SCALES], *mri_laplacian[MAX_SCALES], *mri_dtrans,
      *mri_2nd_deriv_s[MAX_SCALES];
  int i;

  if (which_inputs & (CA_INPUT_DTRANS | CA_INPUT_DTRANS_GRAD))
    mri_dtrans = MRIdistanceTransform(mri_labels, NULL, target_label, 10, DTRANS_MODE_SIGNED, NULL);
  for (i = 0; i < nscales; i++) {
    mri_kernel = MRIgaussian1d(sigmas[i], -1);
    mri_smooth[i] = MRIconvolveGaussian(mri_intensity, NULL, mri_kernel);
    if (which_inputs & CA_INPUT_GRADIENT) mri_grad[i] = MRIsobel(mri_smooth[i], NULL, NULL);
    if (which_inputs & CA_INPUT_LAPLACIAN) mri_laplacian[i] = MRIlaplacian(mri_smooth[i], NULL);
    if (which_inputs & CA_INPUT_D2I_S) {
      mri_2nd_deriv_s[i] = MRI2ndDirectionalDerivative(mri_smooth[i], NULL, 0, -1, 0);
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) MRIwrite(mri_2nd_deriv_s[i], "deriv.mgz");
    }
    MRIfree(&mri_kernel);
  }
  svm_inputs = CAbuildTrainingData(vl_total,
                                   0,
                                   NULL,
                                   mri_smooth,
                                   mri_grad,
                                   mri_laplacian,
                                   mri_dtrans,
                                   mri_2nd_deriv_s,
                                   wsize,
                                   nscales,
                                   which_inputs);

  for (i = 0; i < nscales; i++) {
    MRIfree(&mri_smooth[i]);
    if (which_inputs & CA_INPUT_GRADIENT) MRIfree(&mri_grad[i]);
    if (which_inputs & CA_INPUT_LAPLACIAN) MRIfree(&mri_laplacian[i]);
    if (which_inputs & CA_INPUT_D2I_S) MRIfree(&mri_2nd_deriv_s[i]);
  }
  if (which_inputs & (CA_INPUT_DTRANS | CA_INPUT_DTRANS_GRAD)) MRIfree(&mri_dtrans);
  return (svm_inputs);
}

int ca_ninputs(int which_inputs)
{
  int num = 0;

  if (which_inputs & CA_INPUT_INTENSITY) num++;
  if (which_inputs & CA_INPUT_LAPLACIAN) num++;
  if (which_inputs & CA_INPUT_D2I_S) num++;
  if (which_inputs & CA_INPUT_GRADIENT) num += 3;
  if (which_inputs & CA_INPUT_DTRANS) num++;
  if (which_inputs & CA_INPUT_GRADIENT) num += 3;

  return (num);
}

CA *CAalloc(int width,
            int height,
            int depth,
            MATRIX *m_vox2index,
            int type,
            int which_inputs,
            int wsize,
            int nscales,
            char *c1_name,
            char *c2_name,
            float *sigmas)
{
  CLASSIFIER_ATLAS *ca;
  int x, y, z, o, ninputs;

  ninputs = NINPUTS(which_inputs, wsize, nscales);

  ca = (CA *)calloc(1, sizeof(CLASSIFIER_ATLAS));
  if (ca == NULL) ErrorExit(ERROR_NOMEMORY, "CAalloc: could not allocate CA");
  ca->ninputs = ninputs;
  ca->width = width;
  ca->height = height;
  ca->depth = depth;
  ca->type = type;
  ca->wsize = wsize;
  ca->nscales = nscales;
  ca->sigmas = (float *)calloc(nscales, sizeof(float));
  ca->which_inputs = which_inputs;
  memmove(sigmas, ca->sigmas, nscales * sizeof(float));
  ca->m_vox2index = MatrixCopy(m_vox2index, NULL);
  strcpy(ca->c1_name, c1_name);
  strcpy(ca->c2_name, c2_name);

  ca->classifiers = (CLASSIFIER ****)calloc(width, sizeof(CLASSIFIER ***));
  if (ca->classifiers == NULL) ErrorExit(ERROR_NOMEMORY, "CAalloc: could not allocate %d classifers", width);
  for (x = 0; x < width; x++) {
    ca->classifiers[x] = (CLASSIFIER ***)calloc(height, sizeof(CLASSIFIER **));
    if (ca->classifiers[x] == NULL) ErrorExit(ERROR_NOMEMORY, "CAalloc: could not allocate X %d classifers", height);
    for (y = 0; y < height; y++) {
      ca->classifiers[x][y] = (CLASSIFIER **)calloc(depth, sizeof(CLASSIFIER *));
      if (ca->classifiers[x][y] == NULL)
        ErrorExit(ERROR_NOMEMORY, "CAalloc: could not allocate Y %d classifers", depth);
      for (z = 0; z < depth; z++) {
        ca->classifiers[x][y][z] = (CLASSIFIER *)calloc(ICO0_NVERTICES, sizeof(CLASSIFIER));
        if (ca->classifiers[x][y][z] == NULL)
          ErrorExit(ERROR_NOMEMORY, "CAalloc: could not allocate Z %d classifers", ICO0_NVERTICES);
        for (o = 0; o < ICO0_NVERTICES; o++) {
          ca->classifiers[x][y][z][o].type = type;
          switch (type) {
            case CA_SVM:
              ca->classifiers[x][y][z][o].svm = SVMalloc(ninputs, c1_name, c2_name);
              break;
            case CA_GAUSSIAN:
              ca->classifiers[x][y][z][o].c1_means = (double *)calloc(ninputs, sizeof(double));
              ca->classifiers[x][y][z][o].c1_vars = (double *)calloc(ninputs, sizeof(double));
              ca->classifiers[x][y][z][o].c2_means = (double *)calloc(ninputs, sizeof(double));
              ca->classifiers[x][y][z][o].c2_vars = (double *)calloc(ninputs, sizeof(double));
              if (ca->classifiers[x][y][z][o].c1_means == NULL || ca->classifiers[x][y][z][o].c1_vars == NULL ||
                  ca->classifiers[x][y][z][o].c2_means == NULL || ca->classifiers[x][y][z][o].c2_vars == NULL)
                ErrorExit(ERROR_NOMEMORY, "CAalloc: could not alloc (%d, %d, %d, %d)", x, y, z, o);
              break;
          }
        }
      }
    }
  }

  return (ca);
}

int CAtrain(CA *ca, VOXEL_LIST *vl, MRI *mri_norm, MRI *mri_aseg, int source_label, int target_label)
{
  float *classes, **inputs, nx, ny, nz, mag;
  int i, xi, yi, zi, o, vertices[3], v;
  MRI *mri_normals;
  CLASSIFIER *cl;

  classes = (float *)calloc(vl->nvox, sizeof(float));
  if (classes == NULL) ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d class vector", Progname, vl->nvox);
  for (i = 0; i < vl->nvox; i++) classes[i] = ((int)vl->vdst[i] == target_label) ? 1 : -1;
  inputs = CAbuildInputs(vl, mri_norm, mri_aseg, target_label, ca->which_inputs, ca->wsize, ca->nscales, ca->sigmas);

  ca->c1_label = target_label;
  ca->c2_label = source_label;
  VLSTtransformCoords(vl, ca->m_vox2index, 0);
  mri_normals = caComputeSurfaceNormals(mri_aseg, NULL, source_label);

  for (i = 0; i < vl->nvox; i++) {
    xi = nint(vl->xd[i]);
    yi = nint(vl->yd[i]);
    zi = nint(vl->zd[i]);
    if (xi == Gx && yi == Gy && zi == Gz) DiagBreak();
    if (vl->xi[i] == Gx && vl->yi[i] == Gy && vl->zi[i] == Gz) DiagBreak();
    if (i == Gdiag_no) DiagBreak();
    if (xi < 0 || yi < 0 || zi < 0 || xi >= ca->width || yi >= ca->height || zi >= ca->depth) continue;

    nx = MRIgetVoxVal(mri_normals, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    ny = MRIgetVoxVal(mri_normals, vl->xi[i], vl->yi[i], vl->zi[i], 1);
    nz = MRIgetVoxVal(mri_normals, vl->xi[i], vl->yi[i], vl->zi[i], 2);
    mag = sqrt(nx * nx + ny * ny + nz * nz);
    if (FZERO(mag))
      nz = 1;  // arbitrary
    else {
      nx /= mag;
      ny /= mag;
      nz /= mag;
    }
    IcoFindNClosestVertices(ic0_vertices, ICO0_NVERTICES, nx, ny, nz, 3, vertices);
    for (v = 0; v < 3; v++) {
      o = vertices[v];
      if (o < 0) DiagBreak();
      if (vl->xi[i] == Gx && vl->yi[i] == Gy && vl->zi[i] == Gz) {
        printf("CAtrain voxel (%d, %d, %d) at %d, CL = (%d, %d, %d, %d)\n", Gx, Gy, Gz, i, xi, yi, zi, o);
        DiagBreak();
      }
      cl = &ca->classifiers[xi][yi][zi][o];
      switch (ca->type) {
        case CA_SVM:
          SVMtrain(cl->svm, inputs, classes, vl->nvox, ca->svm_C, ca->tol, ca->max_iter);
          break;
        default:
          clUpdateStatistics(ca, cl, inputs[i], classes[i]);
          if (vl->xi[i] == Gx && vl->yi[i] == Gy && vl->zi[i] == Gz) {
            printf("CAtrain: ntraining (%d, %d)\n", cl->c1_ntraining, cl->c2_ntraining);
            DiagBreak();
          }
          break;
      }
    }
  }

  for (i = 0; i < vl->nvox; i++) free(inputs[i]);
  free(classes);
  free(inputs);
  MRIfree(&mri_normals);
  return (NO_ERROR);
}

int CAsetSVMparms(CA *ca, double svm_C, double svm_tol, int svm_max_iter)
{
  ca->svm_C = svm_C;
  ca->tol = svm_tol;
  ca->max_iter = svm_max_iter;
  return (NO_ERROR);
}
int CAvoxelToIndex(CA *ca, double x, double y, double z, double *pxd, double *pyd, double *pzd)
{
  static VECTOR *v1 = NULL, *v2;

  if (v1 == NULL) {
    v1 = VectorAlloc(4, MATRIX_REAL);
    v2 = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0;
  }

  V3_X(v1) = x;
  V3_Y(v1) = y;
  V3_Z(v1) = z;
  MatrixMultiply(ca->m_vox2index, v1, v2);
  *pxd = V3_X(v2);
  *pyd = V3_Y(v2);
  *pzd = V3_Z(v2);
  return (NO_ERROR);
}

static int clUpdateStatistics(CA *ca, CLASSIFIER *cl, float *inputs, int output)
{
  int i;

  if (cl->c1_means == NULL) {
    cl->c1_means = (double *)calloc(ca->ninputs, sizeof(double));
    cl->c1_vars = (double *)calloc(ca->ninputs, sizeof(double));
    if (cl->c1_means == NULL || cl->c1_vars == NULL)
      ErrorExit(ERROR_NOMEMORY, "clUpdateStatistics: couldn't allocate array");
    cl->c2_means = (double *)calloc(ca->ninputs, sizeof(double));
    cl->c2_vars = (double *)calloc(ca->ninputs, sizeof(double));
    if (cl->c2_means == NULL || cl->c2_vars == NULL)
      ErrorExit(ERROR_NOMEMORY, "clUpdateStatistics: couldn't allocate array");
  }
  for (i = 0; i < ca->ninputs; i++) {
    if (output >= 0) {
      cl->c1_means[i] += inputs[i];
      cl->c1_vars[i] += (inputs[i] * inputs[i]);
    }
    else {
      cl->c2_means[i] += inputs[i];
      cl->c2_vars[i] += (inputs[i] * inputs[i]);
    }
  }
  if (output >= 0)
    cl->c1_ntraining++;
  else
    cl->c2_ntraining++;

  return (NO_ERROR);
}

static int clCompleteTraining(CA *ca, CLASSIFIER *cl)
{
  int i;

  if (cl->c1_ntraining > 0) {
    for (i = 0; i < ca->ninputs; i++) {
      cl->c1_means[i] /= cl->c1_ntraining;
      cl->c1_vars[i] = cl->c1_vars[i] / cl->c1_ntraining - cl->c1_means[i] * cl->c1_means[i];
    }
  }

  if (cl->c2_ntraining > 0) {
    for (i = 0; i < ca->ninputs; i++) {
      cl->c2_means[i] /= cl->c2_ntraining;
      cl->c2_vars[i] = cl->c2_vars[i] / cl->c2_ntraining - cl->c2_means[i] * cl->c2_means[i];
    }
  }
  return (NO_ERROR);
}

#define MAX_INPUTS 10000
int CAcompleteTraining(CA *ca)
{
  int x, y, z, o, i, nvars;
  double vars[MAX_INPUTS];
  CLASSIFIER *cl;

  memset(vars, 0, sizeof(vars));
  for (nvars = x = 0; x < ca->width; x++) {
    for (y = 0; y < ca->height; y++) {
      for (z = 0; z < ca->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        for (o = 0; o < ICO0_NVERTICES; o++) {
          cl = &ca->classifiers[x][y][z][o];
          switch (ca->type) {
            case CA_SVM:
              break;
            case CA_GAUSSIAN:
              clCompleteTraining(ca, cl);
              if (cl->c1_ntraining > 0) {
                for (i = 0; i < ca->ninputs; i++)
                  if (!FZERO(cl->c1_vars[i])) {
                    if (i == 0) nvars++;
                    vars[i] += cl->c1_vars[i];
                  }
              }
              break;
          }
        }
      }
    }
  }
  if (ca->type == CA_GAUSSIAN) {
    if (nvars > 0)
      for (i = 0; i < ca->ninputs; i++) vars[i] /= nvars;

    for (x = 0; x < ca->width; x++) {
      for (y = 0; y < ca->height; y++) {
        for (z = 0; z < ca->depth; z++) {
          if (x == Gx && y == Gy && z == Gz) DiagBreak();
          for (o = 0; o < ICO0_NVERTICES; o++) {
            cl = &ca->classifiers[x][y][z][o];
            if (cl->c1_ntraining > 0) {
              for (i = 0; i < ca->ninputs; i++) {
                if (FZERO(cl->c1_vars[i])) cl->c1_vars[i] = vars[i];
              }
            }
            if (cl->c2_ntraining > 0) {
              for (i = 0; i < ca->ninputs; i++) {
                if (FZERO(cl->c2_vars[i])) cl->c2_vars[i] = vars[i];
              }
            }
          }
        }
      }
    }
  }
  return (NO_ERROR);
}
int CAwrite(CA *ca, char *fname)
{
  FILE *fp;
  int i, x, y, z, o;
  CLASSIFIER *cl;

  fp = fopen(fname, "wb");
  if (fp == NULL) ErrorExit(ERROR_NOFILE, "CAwrite(%s): could not open file for writing", fname);

  fwriteInt(ca->width, fp);
  fwriteInt(ca->height, fp);
  fwriteInt(ca->depth, fp);
  fwriteInt(ca->ninputs, fp);
  fwriteInt(ca->nscales, fp);
  fwriteInt(ca->wsize, fp);
  fwriteInt(ca->which_inputs, fp);
  fwriteInt(ca->type, fp);
  for (i = 0; i < ca->nscales; i++) fwriteDouble(ca->sigmas[i], fp);

  fwriteInt(ca->max_iter, fp);
  fwriteInt(ca->c1_label, fp);
  fwriteInt(ca->c2_label, fp);
  fwriteDouble(ca->svm_C, fp);
  fwriteDouble(ca->tol, fp);

  for (x = 0; x < ca->width; x++) {
    for (y = 0; y < ca->height; y++) {
      for (z = 0; z < ca->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) DiagBreak();
        for (o = 0; o < ICO0_NVERTICES; o++) {
          cl = &ca->classifiers[x][y][z][o];
          switch (ca->type) {
            case CA_GAUSSIAN:
              break;
            case CA_SVM:
              fwriteInt(cl->c1_ntraining, fp);
              fwriteInt(cl->c2_ntraining, fp);
              for (i = 0; i < ca->ninputs; i++) {
                fwriteDouble(cl->c1_means[i], fp);
                fwriteDouble(cl->c1_vars[i], fp);
                fwriteDouble(cl->c2_means[i], fp);
                fwriteDouble(cl->c2_vars[i], fp);
              }
              break;
          }
        }
      }
    }
  }

  fclose(fp);
  return (NO_ERROR);
}
CA *CAread(char *fname)
{
  CA *ca = NULL;
#if 0
  float  sigmas[MAX_SCALES] ;
  int    width, height, depth, wsize, which_inputs, type, nscales ;
  char   c1_name[STRLEN], c2_name[STRLEN] ;

  ca = CAalloc(width, height, depth, m_vox2index, type, which_inputs, 
            wsize, nscales, c1_name, c2_name, sigmas) ;
#endif
  return (ca);
}

MRI *CAclassifyBorder(CA *ca, MRI *mri_norm, MRI *mri_aseg, MRI *mri_output, int border, int label)
{
  MRI *mri_border, *mri_tmp = NULL, *mri_normals;
  int i;
  VOXEL_LIST *vl;
  float **ca_inputs, out;
  MATRIX *m_xform;

  m_xform = compute_ras_basis_vectors(mri_aseg, mri_aseg, label, ca->width, ca->height, ca->depth, 4);

  MatrixCopy(m_xform, ca->m_vox2index);
  MatrixFree(&m_xform);

  // compute segmentation surface normal directions
  mri_normals = caComputeSurfaceNormals(mri_aseg, NULL, label);

  mri_border = MRImarkLabelBorderVoxels(mri_aseg, NULL, label, 1, 1);
  while (border-- > 1) {
    mri_tmp = MRIdilate(mri_border, NULL);
    MRIcopy(mri_tmp, mri_border);
    MRIfree(&mri_tmp);
  }
  mri_output = MRIcloneDifferentType(mri_aseg, MRI_FLOAT);

  vl = VLSTcreate(mri_border, 1, 1, NULL, 0, 0);

  ca_inputs = CAbuildInputs(vl, mri_norm, mri_aseg, label, ca->which_inputs, ca->wsize, ca->nscales, ca->sigmas);
  for (i = 0; i < vl->nvox; i++) {
    out = CAclassifyVoxel(ca, mri_normals, vl->xi[i], vl->yi[i], vl->zi[i], ca_inputs[i]);
    MRIsetVoxVal(mri_output, vl->xi[i], vl->yi[i], vl->zi[i], 0, out);
  }

  return (mri_output);
}

float CAclassifyVoxel(CA *ca, MRI *mri_normals, int x, int y, int z, float *inputs)
{
  double xd, yd, zd, nx, ny, nz, mag, output, c1d, c2d;
  int xi, yi, zi, o, i;
  CLASSIFIER *cl;

  CAvoxelToIndex(ca, x, y, z, &xd, &yd, &zd);
  xi = nint(xd);
  yi = nint(yd);
  zi = nint(zd);
  nx = MRIgetVoxVal(mri_normals, x, y, z, 0);
  ny = MRIgetVoxVal(mri_normals, x, y, z, 1);
  nz = MRIgetVoxVal(mri_normals, x, y, z, 2);
  mag = sqrt(nx * nx + ny * ny + nz * nz);
  if (FZERO(mag))
    o = 0;
  else
    o = IcoFindClosestVertex(ic0_vertices, ICO0_NVERTICES, nx, ny, nz);
  if (xi < 0 || xi >= ca->width || yi < 0 || yi >= ca->height || zi < 0 || zi >= ca->depth) return (0.0);

  if (x == Gx && y == Gy && z == Gz && ((Gdiag_no >= 0 && Gdiag_no == o) || Gdiag_no < 0)) DiagBreak();
  cl = &ca->classifiers[xi][yi][zi][o];
  if (xi == Gx && yi == Gy && zi == Gz) {
    printf("classifying voxel (%d, %d, %d), CL = (%d, %d, %d, %d), NT=(%d, %d)\n",
           Gx,
           Gy,
           Gz,
           xi,
           yi,
           zi,
           o,
           cl->c1_ntraining,
           cl->c2_ntraining);
    DiagBreak();
  }
  if (cl->c1_ntraining == 0 || cl->c2_ntraining == 0) return (0.0);

  // compute inputs to classifier
  for (i = 0, output = 0; i < ca->ninputs; i++) {
    c1d = SQR(inputs[i] - cl->c1_means[i]) / cl->c1_vars[i];
    c2d = SQR(inputs[i] - cl->c2_means[i]) / cl->c2_vars[i];
    output += (sqrt(c1d) - sqrt(c2d));
  }
  return (output);
}
static MRI *caComputeSurfaceNormals(MRI *mri_aseg, MRI *mri_normals, int label)
{
  MRI *mri_tmp, *mri_tmp2, *mri_ctrl;
  int i;

  mri_ctrl = MRIcloneDifferentType(mri_aseg, MRI_UCHAR);
  mri_normals = MRIsegmentationSurfaceNormals(mri_aseg, NULL, label, &mri_ctrl);
  for (i = 0; i < 3; i++)  // dx, dy, dz
  {
    mri_tmp = MRIcopyFrame(mri_normals, NULL, i, 0);
    mri_tmp2 = MRIsoapBubble(mri_tmp, mri_ctrl, NULL, 50, -1);
    MRIsoapBubble(mri_tmp2, mri_ctrl, mri_tmp, 50, -1);
#if 0
    MRIsoapBubble(mri_tmp, mri_ctrl, mri_tmp2, 50, -1) ;
    MRImean(mri_tmp2, mri_tmp, 5) ;
#endif
    MRIcopyFrame(mri_tmp, mri_normals, 0, i);
    MRIfree(&mri_tmp);
    MRIfree(&mri_tmp2);
  }
  MRInormalizeFrameVectorLength(mri_normals, mri_normals);
  MRIfree(&mri_ctrl);
  return (mri_normals);
}
static MATRIX *compute_ras_basis_vectors(
    MRI *mri_aseg_orig, MRI *mri_aseg_edit, int label, int width, int height, int depth, int pad)
{
  MRI_REGION box1, box2, box;
  MRI *mri_aligned, *mri_tmp;
  MATRIX *m_xform, *m_trans, *m_tmp, *m_id, *m_inv, *m_targ, *m_src, *m_tmp2, *m_evectors;
  int i, j, max_col;
  double dot, max_dot, means[3];
  float evalues[3];

  m_evectors = MatrixAlloc(3, 3, MATRIX_REAL); /* eigenvectors of label */
  MRIprincipleComponentsRange(mri_aseg_orig, m_evectors, evalues, means, label, label);
  m_xform = MatrixIdentity(4, NULL);
  m_trans = MatrixIdentity(4, NULL);
  *MATRIX_RELT(m_trans, 1, 4) = -means[0];
  *MATRIX_RELT(m_trans, 2, 4) = -means[1];
  *MATRIX_RELT(m_trans, 3, 4) = -means[2];
  m_inv = MatrixInverse(m_evectors, NULL);

  // rearrange vectors to be as close to RAS as possible
  m_id = MatrixIdentity(3, NULL);
  m_tmp = MatrixMultiply(m_id, m_inv, NULL);  // matrix of dot products
  for (i = 1; i <= 3; i++) {
    // find col in the ithrow that is max
    max_col = 0;
    max_dot = 0;
    for (j = 1; j <= 3; j++) {
      dot = *MATRIX_RELT(m_tmp, i, j);
      if (fabs(dot) > fabs(max_dot)) {
        max_dot = dot;
        max_col = j;
      }
    }
    MatrixCopyRegion(m_inv, m_xform, 1, max_col, 3, 1, 1, i);
    if (max_dot < 0)  // reverse eigenvector
      for (j = 1; j <= 3; j++) *MATRIX_RELT(m_xform, j, i) *= -1;
  }
  MatrixFree(&m_tmp);
  MatrixFree(&m_inv);
  m_tmp = MatrixMultiply(m_xform, m_trans, NULL);
  *MATRIX_RELT(m_trans, 1, 4) = means[0];
  *MATRIX_RELT(m_trans, 2, 4) = means[1];
  *MATRIX_RELT(m_trans, 3, 4) = means[2];
  MatrixMultiply(m_trans, m_tmp, m_xform);

  // now compute a transform that takes the bounding box to the desired width/height/depth
  mri_tmp = MRIlinearTransformInterp(mri_aseg_edit, NULL, m_xform, SAMPLE_NEAREST);
  mri_aligned = MRIdilateLabel(mri_tmp, NULL, label, pad);
  MRIlabelBoundingBox(mri_aligned, label, &box1);
  MRIfree(&mri_aligned);
  MRIfree(&mri_tmp);
  mri_tmp = MRIlinearTransformInterp(mri_aseg_orig, NULL, m_xform, SAMPLE_NEAREST);
  mri_aligned = MRIdilateLabel(mri_tmp, NULL, label, pad);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_aligned, "a.mgz");
  }
  MRIlabelBoundingBox(mri_aligned, label, &box2);
  MRIfree(&mri_aligned);
  MRIfree(&mri_tmp);

  box.x = MIN(box1.x, box2.x);
  box.y = MIN(box1.y, box2.y);
  box.z = MIN(box1.z, box2.z);
  box.dx = MAX(box1.dx + box1.x, box2.dx + box2.x) - box.x;
  box.dy = MAX(box1.dy + box1.y, box2.dy + box2.y) - box.y;
  box.dz = MAX(box1.dz + box1.z, box2.dz + box2.z) - box.z;

  // now compute transform that takes corners of bounding box to corners of atlas
  m_targ = MatrixAlloc(4, 5, MATRIX_REAL);
  m_src = MatrixAlloc(4, 5, MATRIX_REAL);
  *MATRIX_RELT(m_targ, 1, 1) = 0;
  *MATRIX_RELT(m_targ, 2, 1) = 0;
  *MATRIX_RELT(m_targ, 3, 1) = 0;
  *MATRIX_RELT(m_src, 1, 1) = box.x;
  *MATRIX_RELT(m_src, 2, 1) = box.y;
  *MATRIX_RELT(m_src, 3, 1) = box.z;

  *MATRIX_RELT(m_targ, 1, 2) = 0;
  *MATRIX_RELT(m_targ, 2, 2) = 0;
  *MATRIX_RELT(m_targ, 3, 2) = depth;
  *MATRIX_RELT(m_src, 1, 2) = box.x;
  *MATRIX_RELT(m_src, 2, 2) = box.y;
  *MATRIX_RELT(m_src, 3, 2) = box.z + box.dz;

  *MATRIX_RELT(m_targ, 1, 3) = 0;
  *MATRIX_RELT(m_targ, 2, 3) = height;
  *MATRIX_RELT(m_targ, 3, 3) = 0;
  *MATRIX_RELT(m_src, 1, 3) = box.x;
  *MATRIX_RELT(m_src, 2, 3) = box.y + box.dy;
  *MATRIX_RELT(m_src, 3, 3) = box.z;

  *MATRIX_RELT(m_targ, 1, 4) = width;
  *MATRIX_RELT(m_targ, 2, 4) = 0;
  *MATRIX_RELT(m_targ, 3, 4) = 0;
  *MATRIX_RELT(m_src, 1, 4) = box.x + box.dx;
  *MATRIX_RELT(m_src, 2, 4) = box.y;
  *MATRIX_RELT(m_src, 3, 4) = box.z;

  *MATRIX_RELT(m_targ, 1, 5) = width;
  *MATRIX_RELT(m_targ, 2, 5) = height;
  *MATRIX_RELT(m_targ, 3, 5) = depth;
  *MATRIX_RELT(m_src, 1, 5) = box.x + box.dx;
  *MATRIX_RELT(m_src, 2, 5) = box.y + box.dy;
  *MATRIX_RELT(m_src, 3, 5) = box.z + box.dz;

  for (i = 1; i <= m_src->cols; i++) {
    *MATRIX_RELT(m_src, 4, i) = *MATRIX_RELT(m_targ, 4, i) = 1.0;
  }

  m_inv = MatrixSVDPseudoInverse(m_src, NULL);
  m_tmp = MatrixMultiply(m_targ, m_inv, NULL);
  m_tmp2 = MatrixMultiply(m_tmp, m_xform, NULL);
  MatrixCopy(m_tmp2, m_xform);

  MatrixFree(&m_targ);
  MatrixFree(&m_src);
  MatrixFree(&m_inv);
  MatrixFree(&m_tmp);
  MatrixFree(&m_tmp2);
  MatrixFree(&m_evectors);
  MatrixFree(&m_trans);
  return (m_xform);
}
