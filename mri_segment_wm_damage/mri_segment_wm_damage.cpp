/**
 * @file  mri_segment_wm_damage.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:24 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


//
// mri_segment_wm_damage.c
//
// written by Bruce Fischl
// Nov. 9th ,2000
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2011/03/02 00:04:24 $
// Revision       : $Revision: 1.3 $
//
////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gcamorph.h"
#include "mri.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "utils.h"
#include "matrix.h"
#include "gca.h"
#include "cma.h"
#include "version.h"
#include "label.h"

const char *Progname ;
static void usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;
static MRI *label_damaged_wm(MRI *mri_T1, MRI *mri_PD, MRI *mri_T2, MRI *mri_aseg,
                             VECTOR **v_means, MATRIX **m_covariances,
                             TRANSFORM *transform, GCA *gca) ;
static int compute_damaged_wm_statistics(MRI *mri_T1, MRI *mri_PD,
    MRI *mri_T2, LABEL *l_damaged_wm,
    VECTOR **pv_mean,
    MATRIX **pm_cov);
static int compute_label_statistics(MRI *mri_T1, MRI *mri_PD, MRI *mri_T2,
                                    MRI *mri_aseg, int label,
                                    VECTOR **pv_mean,
                                    MATRIX **pm_inv_cov, int erode) ;


#define COULD_BE_DAMAGED_WM(l)  ((l == Left_Cerebral_White_Matter) || \
                 (l == Left_Cerebral_Cortex) || \
                 (l == Left_Caudate) || \
                 (l == Left_Lateral_Ventricle) || \
                 (l == WM_hypointensities) || \
                 (l == Right_Caudate) || \
                 (l == Right_Lateral_Ventricle) || \
                 (l == Right_Cerebral_Cortex) || \
                 (l == Right_Cerebral_White_Matter))


#if 0
static int labels[] = {
                        Left_Cerebral_White_Matter,
                        Right_Cerebral_White_Matter,
                        Left_Caudate,
                        Right_Caudate,
                        Left_Lateral_Ventricle,
                        Right_Lateral_Ventricle
                      } ;
#endif
int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs ;
  MRI          *mri_T1, *mri_T2, *mri_PD, *mri_aseg, *mri_damaged_wm ;
  LABEL        *l_damaged_wm ;
  GCA          *gca ;
  TRANSFORM    *transform ;
  Timer start ;
  int          msec, minutes, seconds ;
  char         *PD_fname, *T2_fname, *T1_fname, *aseg_fname, *gca_fname, *xform_fname, *label_fname, *out_fname ;
  MATRIX       *m_inv_covariances[MAX_CMA_LABELS] ;
  VECTOR       *v_means[MAX_CMA_LABELS] ;

  memset(v_means, 0, sizeof(v_means)) ;
  memset(m_inv_covariances, 0, sizeof(m_inv_covariances)) ;
  start.reset() ;
  setRandomSeed(-1L) ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  Progname = argv[0] ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 9)
    usage_exit(1) ;

  PD_fname = argv[1] ;
  T2_fname = argv[2] ;
  T1_fname = argv[3] ;
  aseg_fname = argv[4] ;
  gca_fname = argv[5] ;
  xform_fname = argv[6] ;
  label_fname = argv[7] ;
  out_fname = argv[8] ;
  fprintf(stderr, "reading %s...\n", PD_fname) ;
  mri_PD = MRIread(PD_fname) ;
  if (!mri_PD)
    ErrorExit(ERROR_NOFILE, "%s: could not read PD volume %s",
              Progname, PD_fname) ;


  fprintf(stderr, "reading %s...\n", T2_fname) ;
  mri_T2 = MRIread(T2_fname) ;
  if (!mri_T2)
    ErrorExit(ERROR_NOFILE, "%s: could not read T2 volume %s",
              Progname, T2_fname) ;

  fprintf(stderr, "reading %s...\n", T1_fname) ;
  mri_T1 = MRIread(T1_fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 volume %s",
              Progname, T1_fname) ;

  fprintf(stderr, "reading %s...\n", aseg_fname) ;
  mri_aseg = MRIread(aseg_fname) ;
  if (!mri_aseg)
    ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume %s",
              Progname, aseg_fname) ;

  fprintf(stderr, "reading %s...\n", label_fname) ;
  l_damaged_wm = LabelRead(NULL, label_fname) ;
  if (!l_damaged_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read label from %s",
              Progname, label_fname) ;

  fprintf(stderr, "reading %s...\n", gca_fname) ;
  gca = GCAread(gca_fname) ;
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s: could not read gca from %s",
              Progname, gca_fname) ;

  fprintf(stderr, "reading %s...\n", xform_fname) ;
  transform = TransformRead(xform_fname) ;
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s: could not read xform from %s",
              Progname, xform_fname) ;
  TransformInvert(transform, mri_T1) ;

  compute_damaged_wm_statistics(mri_T1, mri_PD, mri_T2, l_damaged_wm, &v_means[WM_hypointensities],
                                &m_inv_covariances[WM_hypointensities]) ;
  compute_label_statistics(mri_T1, mri_PD, mri_T2, mri_aseg, Left_Cerebral_White_Matter,
                           &v_means[Left_Cerebral_White_Matter], &m_inv_covariances[Left_Cerebral_White_Matter],2) ;
  compute_label_statistics(mri_T1, mri_PD, mri_T2, mri_aseg, Left_Caudate,
                           &v_means[Left_Caudate], &m_inv_covariances[Left_Caudate],1) ;
  compute_label_statistics(mri_T1, mri_PD, mri_T2, mri_aseg, Left_Lateral_Ventricle,
                           &v_means[Left_Lateral_Ventricle], &m_inv_covariances[Left_Lateral_Ventricle],2) ;
  compute_label_statistics(mri_T1, mri_PD, mri_T2, mri_aseg, Right_Cerebral_White_Matter,
                           &v_means[Right_Cerebral_White_Matter], &m_inv_covariances[Right_Cerebral_White_Matter],2) ;
  compute_label_statistics(mri_T1, mri_PD, mri_T2, mri_aseg, Right_Caudate,
                           &v_means[Right_Caudate], &m_inv_covariances[Right_Caudate],1) ;
  compute_label_statistics(mri_T1, mri_PD, mri_T2, mri_aseg, Right_Lateral_Ventricle,
                           &v_means[Right_Lateral_Ventricle], &m_inv_covariances[Right_Lateral_Ventricle],2) ;

  mri_damaged_wm = label_damaged_wm(mri_T1, mri_PD, mri_T2, mri_aseg, v_means, m_inv_covariances, transform, gca) ;
  printf("writing output to %s...\n", out_fname) ;
  MRIwrite(mri_damaged_wm, out_fname) ;

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "wm damage segmentation took %d minutes and %d seconds.\n",
          minutes, seconds) ;
  exit(0) ;
  return(0) ;
}

static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else switch (*option) {
    case '?':
    case 'U':
      usage_exit(1);
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      usage_exit(1) ;
      break ;
    }
  return(nargs) ;
}

static void
usage_exit(int ecode) {
  printf("usage: %s <source> <target> <output xform>\n",Progname) ;
  exit(ecode) ;
}
static int
compute_damaged_wm_statistics(MRI *mri_T1, MRI *mri_PD, MRI *mri_T2, LABEL *l_damaged_wm, VECTOR **pv_mean,
                              MATRIX **pm_cov) {
  MATRIX *m_cov, *m_T1_ras2vox, *m_T2_ras2vox, *m_tmp ;
  VECTOR *v_mean, *v_T1, *v_T2, *v_ras, *v_X, *v_XT ;
  double  xv, yv, zv, val_T1, val_T2, val_PD ;
  int     i, dofs_T1, dofs_T2 ;

  m_cov = MatrixAlloc(3, 3, MATRIX_REAL) ;
  v_mean = VectorAlloc(3, MATRIX_REAL) ;
  v_T1 = VectorAlloc(4, MATRIX_REAL) ;
  v_T2 = VectorAlloc(4, MATRIX_REAL) ;
  v_ras = VectorAlloc(4, MATRIX_REAL) ;
  v_X = VectorAlloc(3, MATRIX_REAL) ;
  v_XT = NULL ;
  m_tmp = NULL ;
  VECTOR_ELT(v_T1, 4) = 1.0 ;
  VECTOR_ELT(v_T2, 4) = 1.0 ;
  VECTOR_ELT(v_ras, 4) = 1.0 ;

  m_T1_ras2vox = MRIgetRasToVoxelXform(mri_T1) ;
  m_T2_ras2vox = MRIgetRasToVoxelXform(mri_T2) ;

  // compute mean for 3 volumes
  for (dofs_T1 = dofs_T2 = i = 0 ; i < l_damaged_wm->n_points ; i++) {
    V3_X(v_ras) = l_damaged_wm->lv[i].x ;
    V3_Y(v_ras) = l_damaged_wm->lv[i].y ;
    V3_Z(v_ras) = l_damaged_wm->lv[i].z ;
    MatrixMultiply(m_T1_ras2vox, v_ras, v_T1) ;
    MatrixMultiply(m_T2_ras2vox, v_ras, v_T2) ;
    xv = V3_X(v_T1) ;
    yv = V3_Y(v_T1) ;
    zv = V3_Z(v_T1) ;
    if (MRIindexNotInVolume(mri_T1, xv, yv, zv) == 0) {
      MRIsampleVolume(mri_T1, xv, yv, zv, &val_T1) ;
      VECTOR_ELT(v_mean, 1) += val_T1 ;
      dofs_T1++ ;
    }
    xv = V3_X(v_T2) ;
    yv = V3_Y(v_T2) ;
    zv = V3_Z(v_T2) ;
    if (MRIindexNotInVolume(mri_T2, xv, yv, zv) == 0) {
      MRIsampleVolume(mri_PD, xv, yv, zv, &val_PD) ;
      VECTOR_ELT(v_mean, 2) += val_PD ;
      MRIsampleVolume(mri_T2, xv, yv, zv, &val_T2) ;
      VECTOR_ELT(v_mean, 3) += val_T2 ;
      dofs_T2++ ;
    }
  }
  VECTOR_ELT(v_mean, 1) /= (float)dofs_T1 ;
  VECTOR_ELT(v_mean, 2) /= (float)dofs_T2 ;
  VECTOR_ELT(v_mean, 3) /= (float)dofs_T2 ;
  printf("mean damaged wm:\n") ;
  MatrixPrint(stdout, v_mean) ;


  // compute covariances for 3 volumes
  for (dofs_T1 = i = 0 ; i < l_damaged_wm->n_points ; i++) {
    V3_X(v_ras) = l_damaged_wm->lv[i].x ;
    V3_Y(v_ras) = l_damaged_wm->lv[i].y ;
    V3_Z(v_ras) = l_damaged_wm->lv[i].z ;
    MatrixMultiply(m_T1_ras2vox, v_ras, v_T1) ;
    MatrixMultiply(m_T2_ras2vox, v_ras, v_T2) ;
    xv = V3_X(v_T1) ;
    yv = V3_Y(v_T1) ;
    zv = V3_Z(v_T1) ;
    if (MRIindexNotInVolume(mri_T1, xv, yv, zv) == 0) {
      MRIsampleVolume(mri_T1, xv, yv, zv, &val_T1) ;
      xv = V3_X(v_T2) ;
      yv = V3_Y(v_T2) ;
      zv = V3_Z(v_T2) ;
      if (MRIindexNotInVolume(mri_T2, xv, yv, zv) == 0)  // only if in both sets of volumes
      {
        MRIsampleVolume(mri_PD, xv, yv, zv, &val_PD) ;
        MRIsampleVolume(mri_T2, xv, yv, zv, &val_T2) ;
        VECTOR_ELT(v_X,1) = val_T1 - VECTOR_ELT(v_mean,1) ;
        VECTOR_ELT(v_X,2) = val_PD - VECTOR_ELT(v_mean,2) ;
        VECTOR_ELT(v_X,3) = val_T2 - VECTOR_ELT(v_mean,3) ;
        v_XT = VectorTranspose(v_X, v_XT) ;
        m_tmp = MatrixMultiply(v_X, v_XT, m_tmp) ;
        MatrixAdd(m_tmp, m_cov, m_cov) ;
        dofs_T1++ ;
      }
    }
  }

  if (dofs_T1 > 0)
    MatrixScalarMul(m_cov, 1.0/(float)dofs_T1, m_cov) ;
  printf("covariance matrix:\n") ;
  MatrixPrint(stdout, m_cov) ;
  VectorFree(&v_T1) ;
  VectorFree(&v_T2) ;
  VectorFree(&v_ras) ;
  VectorFree(&v_X) ;
  VectorFree(&v_XT) ;
  MatrixFree(&m_T1_ras2vox) ;
  MatrixFree(&m_T2_ras2vox) ;
  MatrixFree(&m_tmp) ;
  *pv_mean = v_mean ;
  *pm_cov = MatrixInverse(m_cov, NULL) ;
  MatrixFree(&m_cov) ;
  return(NO_ERROR) ;
}

static int
compute_label_statistics(MRI *mri_T1, MRI *mri_PD, MRI *mri_T2, MRI *mri_aseg, int label, VECTOR **pv_mean,
                         MATRIX **pm_inv_cov, int erode) {
  MATRIX *m_cov, *m_T1_to_T2, *m_tmp ;
  VECTOR *v_mean, *v_T1, *v_T2, *v_X, *v_XT ;
  double  xv, yv, zv, val_T1, val_T2, val_PD ;
  int     dofs_T1, dofs_T2, x, y, z ;
  MRI     *mri_label ;

  mri_label = MRIclone(mri_aseg, NULL) ;
  MRIcopyLabel(mri_aseg, mri_label, label) ;
  while (erode-- > 0)
    MRIerode(mri_label, mri_label) ;

  m_cov = MatrixAlloc(3, 3, MATRIX_REAL) ;
  v_mean = VectorAlloc(3, MATRIX_REAL) ;
  v_T1 = VectorAlloc(4, MATRIX_REAL) ;
  v_T2 = VectorAlloc(4, MATRIX_REAL) ;
  v_X = VectorAlloc(3, MATRIX_REAL) ;
  v_XT = NULL ;
  m_tmp = NULL ;
  VECTOR_ELT(v_T1, 4) = 1.0 ;
  VECTOR_ELT(v_T2, 4) = 1.0 ;

  m_T1_to_T2 = MRIgetVoxelToVoxelXform(mri_T1, mri_T2) ;

  // compute mean for 3 volumes
  for (dofs_T1 = dofs_T2 = x = 0 ; x < mri_T1->width ; x++) {
    V3_X(v_T1) = x ;
    for (y = 0 ; y < mri_T1->height ; y++) {
      V3_Y(v_T1) = y ;
      for (z = 0 ; z < mri_T1->depth ; z++) {
        if (nint(MRIgetVoxVal(mri_label, x, y, z, 0)) != label)
          continue ;
        V3_Z(v_T1) = z ;
        MatrixMultiply(m_T1_to_T2, v_T1, v_T2) ;
        xv = V3_X(v_T2) ;
        yv = V3_Y(v_T2) ;
        zv = V3_Z(v_T2) ;
        val_T1 = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
        VECTOR_ELT(v_mean, 1) += val_T1 ;
        dofs_T1++ ;

        xv = V3_X(v_T2) ;
        yv = V3_Y(v_T2) ;
        zv = V3_Z(v_T2) ;
        if (MRIindexNotInVolume(mri_T2, xv, yv, zv) == 0) {
          MRIsampleVolume(mri_PD, xv, yv, zv, &val_PD) ;
          VECTOR_ELT(v_mean, 2) += val_PD ;
          MRIsampleVolume(mri_T2, xv, yv, zv, &val_T2) ;
          VECTOR_ELT(v_mean, 3) += val_T2 ;
          dofs_T2++ ;
        }
      }
    }
  }
  VECTOR_ELT(v_mean, 1) /= (float)dofs_T1 ;
  VECTOR_ELT(v_mean, 2) /= (float)dofs_T2 ;
  VECTOR_ELT(v_mean, 3) /= (float)dofs_T2 ;
  printf("mean %s (%d):\n", cma_label_to_name(label), label) ;
  MatrixPrint(stdout, v_mean) ;


  // compute covariances for 3 volumes
  for (dofs_T1 = dofs_T2 = x = 0 ; x < mri_T1->width ; x++) {
    V3_X(v_T1) = x ;
    for (y = 0 ; y < mri_T1->height ; y++) {
      V3_Y(v_T1) = y ;
      for (z = 0 ; z < mri_T1->depth ; z++) {
        if (nint(MRIgetVoxVal(mri_label, x, y, z, 0)) != label)
          continue ;
        V3_Z(v_T1) = z ;
        MatrixMultiply(m_T1_to_T2, v_T1, v_T2) ;
        val_T1 = MRIgetVoxVal(mri_T1, x, y, z, 0) ;
        xv = V3_X(v_T2) ;
        yv = V3_Y(v_T2) ;
        zv = V3_Z(v_T2) ;
        if (MRIindexNotInVolume(mri_T2, xv, yv, zv) == 0)  // only if in both sets of volumes
        {
          MRIsampleVolume(mri_PD, xv, yv, zv, &val_PD) ;
          MRIsampleVolume(mri_T2, xv, yv, zv, &val_T2) ;
          VECTOR_ELT(v_X,1) = val_T1 - VECTOR_ELT(v_mean,1) ;
          VECTOR_ELT(v_X,2) = val_PD - VECTOR_ELT(v_mean,2) ;
          VECTOR_ELT(v_X,3) = val_T2 - VECTOR_ELT(v_mean,3) ;
          v_XT = VectorTranspose(v_X, v_XT) ;
          m_tmp = MatrixMultiply(v_X, v_XT, m_tmp) ;
          MatrixAdd(m_tmp, m_cov, m_cov) ;
          dofs_T1++ ;
        }
      }
    }
  }

  if (dofs_T1 > 0)
    MatrixScalarMul(m_cov, 1.0/(float)dofs_T1, m_cov) ;
  printf("covariance matrix:\n") ;
  MatrixPrint(stdout, m_cov) ;
  VectorFree(&v_T1) ;
  VectorFree(&v_T2) ;
  VectorFree(&v_X) ;
  VectorFree(&v_XT) ;
  MatrixFree(&m_T1_to_T2) ;
  MatrixFree(&m_tmp) ;
  *pv_mean = v_mean ;
  *pm_inv_cov = MatrixInverse(m_cov, NULL) ;
  MatrixFree(&m_cov) ;
  MRIfree(&mri_label) ;
  return(NO_ERROR) ;
}

#ifdef MIN_PRIOR
#undef MIN_PRIOR
#endif
#define MIN_PRIOR 0.01
static MRI *
label_damaged_wm(MRI *mri_T1, MRI *mri_PD, MRI *mri_T2, MRI *mri_aseg,
                 VECTOR **v_means, MATRIX **m_inv_covariances, TRANSFORM *transform, GCA *gca) {
  MRI       *mri_damaged_wm, *mri_mask ;
  int       x, y, z, l, xv, yv, zv, min_dist_label, ndamaged ;
  float     dists[MAX_CMA_LABELS], pwm_damage, prior, min_dist, dets[MAX_CMA_LABELS] ;
  MATRIX    *m_T2_to_T1 ;
  VECTOR    *v_T1, *v_T2, *v_vals, *v_vals_minus_means, *v_tmp ;
  GCA_PRIOR *gcap ;
  double    val ;

  for (l = 0 ; l < MAX_CMA_LABELS ; l++) {
    if (m_inv_covariances[l]  != NULL) {
      dets[l] = 1/MatrixDeterminant(m_inv_covariances[l]) ;
      dets[l] = log(sqrt(dets[l])) ; // for log likelihood
    }
  }
  mri_mask = MRIclone(mri_aseg, NULL) ;
  MRIcopyLabel(mri_aseg, mri_mask, Left_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_mask, Right_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_mask, Left_Caudate) ;
  MRIcopyLabel(mri_aseg, mri_mask, Right_Caudate) ;
  MRIcopyLabel(mri_aseg, mri_mask, Left_Lateral_Ventricle) ;
  MRIcopyLabel(mri_aseg, mri_mask, Right_Lateral_Ventricle) ;
  MRIcopyLabel(mri_aseg, mri_mask, WM_hypointensities) ;
  MRIcopyLabel(mri_aseg, mri_mask, Left_Thalamus) ;
  MRIcopyLabel(mri_aseg, mri_mask, Right_Thalamus) ;
  MRIcopyLabel(mri_aseg, mri_mask, Left_Putamen) ;
  MRIcopyLabel(mri_aseg, mri_mask, Right_Putamen) ;
  MRIcopyLabel(mri_aseg, mri_mask, Left_Pallidum) ;
  MRIcopyLabel(mri_aseg, mri_mask, Right_Pallidum) ;
  MRIclose(mri_mask, mri_mask) ;
  MRIerode(mri_mask, mri_mask) ;
  MRIerode(mri_mask, mri_mask) ;  // don't let it occur near cortex
  MRIbinarize(mri_mask, mri_mask, 1, 0, 128) ;
  MRIwrite(mri_mask, "m.mgz") ;
  v_tmp = NULL ;
  v_T1 = VectorAlloc(4, MATRIX_REAL) ;
  v_T2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_T1,4) = VECTOR_ELT(v_T2,4) = 1.0 ;
  m_T2_to_T1 = MRIgetVoxelToVoxelXform(mri_T2, mri_T1) ;
  v_vals = VectorAlloc(3, MATRIX_REAL) ;
  v_vals_minus_means = VectorAlloc(3, MATRIX_REAL) ;

  mri_damaged_wm = MRIclone(mri_PD, NULL) ;

  for (ndamaged = x = 0 ; x < mri_PD->width ; x++) {
    V3_X(v_T2) = x ;
    for (y = 0 ; y < mri_PD->height ; y++) {
      V3_Y(v_T2) = y ;
      for (z = 0 ; z < mri_PD->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_Z(v_T2) = z ;
        MatrixMultiply(m_T2_to_T1, v_T2, v_T1) ;
        xv = nint(V3_X(v_T1)) ;
        yv = nint(V3_Y(v_T1)) ;
        zv = nint(V3_Z(v_T1)) ;
        if (MRIindexNotInVolume(mri_T1, xv, yv, zv) != 0)
          continue ;
        if (nint(MRIgetVoxVal(mri_mask, xv, yv, zv,0)) == 0)
          continue ;

        gcap = getGCAP(gca, mri_T1, transform, xv, yv, zv) ;
        pwm_damage =
          getPrior(gcap, Left_Cerebral_White_Matter) +
          getPrior(gcap, Right_Cerebral_White_Matter) +
          getPrior(gcap, Left_Caudate) +
          getPrior(gcap, Right_Caudate) +
          getPrior(gcap, WM_hypointensities) ;
        if (pwm_damage < MIN_PRIOR)   // only occurs where one of the above classes is possible
          continue ;
        memset(dists, 0, sizeof(dists)) ;
        MRIsampleVolume(mri_T1, xv, yv, zv, &val) ;
        VECTOR_ELT(v_vals, 1) = val ;
        VECTOR_ELT(v_vals, 2) = MRIgetVoxVal(mri_PD, x, y, z, 0) ;
        VECTOR_ELT(v_vals, 3) = MRIgetVoxVal(mri_T2, x, y, z, 0) ;

        min_dist = -1 ;
        min_dist_label = 0 ;
        for (l = 0 ; l < MAX_CMA_LABELS ; l++)
          if (v_means[l]) {
            prior = getPrior(gcap, l) ;
            if (prior < MIN_PRIOR && l != WM_hypointensities)  // don't force prior check for hypos - they occur all over
              continue ;
            VectorSubtract(v_vals, v_means[l], v_vals_minus_means) ;
            v_tmp = MatrixMultiply(m_inv_covariances[l], v_vals_minus_means,v_tmp);
            dists[l] = VectorDot(v_vals_minus_means, v_tmp) + dets[l] ;
            if ((dists[l] < min_dist) || min_dist < 0) {
              min_dist = dists[l] ;
              min_dist_label = l ;
            }
          }
        if (min_dist_label == WM_hypointensities)
          ndamaged++ ;
        MRIsetVoxVal(mri_damaged_wm, x, y, z, 0, min_dist_label) ;
      }
    }
  }

  printf("%d damaged white matter voxels detected\n", ndamaged) ;
  MatrixFree(&m_T2_to_T1) ;
  VectorFree(&v_T1) ;
  VectorFree(&v_T2) ;
  VectorFree(&v_vals) ;
  VectorFree(&v_vals_minus_means) ;
  VectorFree(&v_tmp) ;
  return(mri_damaged_wm) ;
}

